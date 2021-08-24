	!> Element for pos processing variables and store in Param
	!! @param CommonPar = ( ... , a, a*[bi,ci]) where : 
	!! @param iShiftU = CommonPar(1)
	!! @param iFemType = CommonPar(2)
	!! @param iFtype = CommonPar(3)
	!! @param iMaterial = CommonPar(4) 
	!! @param iDamageParType = CommonPar(5)
	!! @param Nelem = CommonPar(6)
	!! @param NFields = CommonPar(7)  ===> a
	!! @param CodeField(i) = CommonPar(see structure) ==> bi
	!! @param PosField(i) = CommonPar(see structure) ==> ci	
	!! @author Rocha, Felipe Figueredo
	
 Subroutine posProcElem(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use finiteStrainLib !! New way
	use damageNewLib
	use multiscaleLib
	use globalVariables, only : getF, NdimE, getDamagePar, getMaterial, maxDamagePar, maxMatPar, &
								addVol
	use ptsGaussLib
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodElt ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
         
	Integer :: nG, i! counters 
	Integer , parameter :: iShiftField = 8,  NFieldPar = 2
	Real*8 :: GradU(NdimE,NdimE) , F(NdimE,NdimE), NvecP(3)
	Real*8 :: PK(NdimE,NdimE), D(NdimE,NdimE,NdimE,NdimE), sigma(NdimE,NdimE) , detF, p_ng
	Real*8 :: MatPar(maxMatPar) ,  damagePar(maxDamagePar), VField(4) , energy, energyMean, dV
	real*8 , allocatable :: SolU(:) ,  Xel(:) , SolP(:)
	integer :: pOrderU , pOrderP, NGP, NodGU, NodGP, iFEMtypeU , iFEMtypeP, iShiftU, iShiftP 
	integer :: constlaw, iFemType,  iSimplex,  iBubble, iFtype, iMaterial, iDamageParType , NFields, CodeField(5) , PosField(5) !! maximum 5 at now
	type(ptGaussClass) :: PtGU, PtGP

	matPar = 0.0d0

	iShiftU= nint(commonPar(1))
	iFemTypeU = nint(commonPar(2))
	iShiftP= nint(commonPar(3))
	iFemTypeP = nint(commonPar(4))
	iFtype = nint(commonPar(5))
	iMaterial = nint(commonPar(6))	 
	iDamageParType = nint(commonPar(7))
	NFields = nint(CommonPar(8))
	do i = 1, NFields
		CodeField(i) = nint(CommonPar( iShiftField + (i-1)*NFieldPar + 1) )
		PosField(i)  = nint(CommonPar( iShiftField + (i-1)*NFieldPar + 2) )
	end do
	
	if(iDamageParType > 0) call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw, matPar, iMaterial)
	call setFEMtype(iFEMtypeU,NodGU,pOrderU,NGP,iSimplex,iBubble)
	call setFEMtype(iFEMtypeP,NodGP,pOrderP,NGP,iSimplex,iBubble)
	
	call getSliceAllocate(SolU,Sol1,1,NodGU,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodGU,1 ,NdimE, Ndim)
	call getSliceAllocate(SolP,Sol1,1,NodGP,iShiftP + 1 ,iShiftP + 1, iDofT)
	
	NGP = 1
	nG = 1
	
	call PtGU%init(Xel,NodGU,NdimE,NGP,pOrderU,iBubble, iSimplex)
	call PtGP%init(Xel,NodGP,NdimE,NGP,pOrderP,iBubble, iSimplex)
	
	call PtGU%calcGradU(GradU,SolU,nG)	
	call getF(F,iFtype)
	F = F + GradU
	call calcPpk(PK,F,NdimE,matPar,constLaw)		
	call calcD(D,F,NdimE,matPar,constLaw)	 
!~ 	if(iDamageParType > 0) call damageModifyTangent(D,PK,damagePar,Param,nG, 4) ! deprecated
	call strainEnergy(energy,F,matPar,constLaw)
	
	call PtGP%getNvecScalar(NvecP,nG)
	P_ng = dot_product(NvecP,SolP)
		
	call MatDet(detF,F)
	
	sigma = (1.0d0/detF)*matmul(PK,transpose(F)) + p_ng*deltaKron(1:NdimE,1:NdimE)
	
	dV=ptGU%dV(nG)
		
	do i = 1 , NFields
		select case(CodeField(i))
			case(1) !! energy 
				VField(i) = energy
			case(2) !! vol
				VField(i) = dV
			case(3) !! sigmaXX
				VField(i) = sigma(1,1)
			case(4) !! sigmaYY
				VField(i) = sigma(2,2)
			case(5) !! sigmaXY
				VField(i) = sigma(1,2)
			case(6) !! pressure
				VField(i) = P_ng
			case(7) !! detF
				VField(i) = detF
			case(8) !! assimetry
				VField(i) = sigma(2,1) - sigma(1,2)
			case(11) !! P11
				VField(i) = PK(1,1)
			case(12) !! P12
				VField(i) = PK(1,2)
			case(13) !! P21
				VField(i) = PK(2,1)
			case(14) !! P22
				VField(i) = PK(2,2)
				
			case default
				write(0,*) "Code Field not found"
				VField(i) = 0
		end select
		
!~ 		write(0,*) "=============> PK =", VField(i)
		Param(LengthParam + PosField(i) ) = VField(i)

	end do
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)	
	
End Subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine posProcElemS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow , i
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)

	Coupling=0 
	
	Coupling(1,1) = 1
	
	return
end Subroutine
