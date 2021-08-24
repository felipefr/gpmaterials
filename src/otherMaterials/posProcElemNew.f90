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
	
 Subroutine posProcElemNew(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
 
	use funcAux
	use globalVariables, only : getF, NdimE, getDamagePar, getMaterial, maxDamagePar, maxMatPar
!~ 	use finiteStrainLib
	use genericConstitutiveLib
	use ptsGaussLib
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodElt ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
         
	Integer :: nG,i,j,k,l,m,n,p,q, ip, A, B, ApRow, BpCol , App,Bqp, ipDim, ipDimE, ipDof,ipDofmod, NodG ! counters 
	Integer , parameter :: iShiftField = 8,  NFieldPar = 2
	Real*8 :: GradU(NdimE,NdimE) , F(NdimE,NdimE)
	Real*8 :: Ppk(NdimE,NdimE), D(NdimE,NdimE,NdimE,NdimE), PKmean(NdimE,NdimE), DMean(NdimE,NdimE,NdimE,NdimE)
	Real*8 :: MatPar(maxMatPar) ,  damagePar(maxDamagePar), VField(4) , energy, energyMean, dVng, dVtot
	real*8 , allocatable :: SolU(:) ,  Xel(:) 
	integer :: NGP,pOrder , IisPeriodic , iShiftU, constLaw , NFields, CodeField(10) , PosField(10) !! maximum 5 at now
	integer :: iFemType,  iSimplex,  iBubble, iFtype, iMaterial, iDamageParType,  damageMethod,isAnalytic 
	type(ptGaussClass) :: PtG

	matPar = 0.0d0

	iShiftU= nint(commonPar(1))
	iFemType = nint(commonPar(2))
	iFtype = nint(commonPar(3))
	iMaterial = nint(commonPar(4))	 
	iDamageParType = nint(commonPar(5))
	damageMethod = nint(commonPar(6))
	isAnalytic = nint(commonPar(7))
	NFields = nint(CommonPar(8))
	do i = 1, NFields
		CodeField(i) = nint(CommonPar( iShiftField + (i-1)*NFieldPar + 1) )
		PosField(i)  = nint(CommonPar( iShiftField + (i-1)*NFieldPar + 2) )
	end do
	
	if(iDamageParType > 0) call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw, matPar, iMaterial)
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)
	
	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
	
	PKmean = 0.0d0
	DMean = 0.0d0
	energyMean = 0.0d0
	dVtot = 0.0d0
	
	do nG = 1 , NGP
		call PtG%calcGradU(GradU,SolU,nG)
		
		call getF(F,iFtype)
		
		F = F + GradU
		
		call SUPD_continuum(Ppk,D,F,damagePar,Param,matPar,constlaw,NDimE,damageMethod,isAnalytic)
		
		call strainEnergy(energy,F,matPar,constLaw)
		
		dVng=ptG%dV(nG)
		
		dVtot = dVtot + dVng
		PKmean = PKmean + dVng*Ppk
		DMean = DMean + dVng*D
		energyMean = energyMean + dVng*energy 	
	end do
	
	PKmean = PKmean/dVtot
	DMean = DMean/dVtot
	energyMean = energyMean/dVtot
	
	do i = 1 , NFields
		select case(CodeField(i))
			case(1) !! energy 
				VField(i) = energy
			case(2) !! vol
				VField(i) = dVtot 
			case(3) !! P11
				VField(i) = PKmean(1,1)
			case(4) !! P12
				VField(i) = PKmean(1,2)
			case(5) !! P21
				VField(i) = PKmean(2,1)
			case(6) !! P22
				VField(i) = PKmean(2,2)
								
			case default
				write(0,*) "Code Field not found"
				VField(i) = 0
		end select
		

		Param(LengthParam + PosField(i) ) = VField(i)

	end do
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)	
	
End Subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine posProcElemNewS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow , i
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)

	Coupling=0 
	
	Coupling(1,1) = 1
	
	return
end Subroutine
