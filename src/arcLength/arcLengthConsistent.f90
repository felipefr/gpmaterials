    !> Generic element for finite strains
    !!
    !!
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iShift_dUf = nint(CommonPar(2))
	!! @param iFemType  = nint(CommonPar(3)) 
	!! @param iFType  = nint(commonPar(4))
	!! @param iMaterial  = nint(commonPar(5))
	!! @param iDamageParType  = nint(commonPar(6))
	!!
    !! @author Rocha, Felipe Figueredo

Subroutine arcLengthConsistent(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : getF, NdimE, getDamagePar, getMaterial, maxDamagePar, maxMatPar
	use ptsGaussLib
	use genericConstitutiveLib , only : SUPD_continuum
	
	use TIME_STEP_ITER
	use DETERMINANT, only : iElem
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
         
	Integer :: nG, i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ipDim, ipDimE, ipDof,ipDofmod ! counters 
	Integer :: constLaw ! if NparamDamage = 0 , no damage
	Real*8 :: dV , Det , MatPar(maxMatPar), damagePar(maxDamagePar)
	Real*8 :: GradU(NdimE,NdimE) , F(NdimE,NdimE)
	Real*8 :: D(NdimE,NdimE,NdimE,NdimE), Ppk(NdimE,NdimE), Dmat(NdimE*NdimE,NdimE*NdimE) , Pmat(NdimE*NdimE)
	Real*8 :: Bmat(NdimE*NdimE,NodElt*NdimE), BmatT(NodElt*NdimE,NdimE*NdimE), Auu(NodElt*NdimE,NodElt*NdimE), Bu(NodElt*NdimE)
	real*8 ::  SolUf(NodElt*NdimE) ,  Sol_dUf_acc(NodElt*NdimE), Xel(NodElt*NdimE) !! actually they should be NodG*NdimE
	integer :: pOrder , NGP, NodG, iFEMtype , iFtype, iMaterial, iSimplex, iBubble, sizeSolU, iDamageParType, & 
				damageMethod,isAnalytic
	integer :: iShift_Uf, iShift_dUf, iShift_dUf_acc, iShift_dLarcLength, iShift_LarcLength, iShift_dLarcLength_acc
	type(ptGaussClass) :: PtG
	 
	real*8 :: SolGy(NodElt*NdimE) ,  Sol_dU_acc(NodElt*NdimE) , Nmat(NdimE,NdimE*NodElt)  &
			, dU_acc(NdimE) , SolGy_gp(NdimE)!! actually they should be NodG*NdimE
	real*8 :: A_ll(1,1) , Alu(1,NodElt*NdimE), Aul(NodElt*NdimE,1), Bl(1), Gbar(NdimE,NdimE), lambda, y(NdimE), &
			  arcLength, dlambda_acc
	
	
	arcLength = CommonPar(1)
	iShift_Uf = nint(CommonPar(2))
	iShift_dUf = nint(CommonPar(3))
	iShift_dUf_acc = nint(CommonPar(4))
	iShift_LarcLength = nint(CommonPar(5))
	iShift_dLarcLength = nint(CommonPar(6))
	iShift_dLarcLength_acc = nint(CommonPar(7))
	iFEMtype =  nint(commonPar(8)) 
	iFtype = nint(commonPar(9)) 
	iMaterial = nint(commonPar(10))
	iDamageParType = nint(commonPar(11))
	damageMethod = nint(commonPar(12))
	isAnalytic = nint(commonPar(13))
	
!~ 	iShift_dLarcLength = 4
 	
	if(iDamageParType>0) call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw, matPar, iMaterial)
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)

	call getSlice(SolUf,Sol1,1,NodG,iShift_Uf + 1 ,iShift_Uf + NdimE, iDofT)
	call getSlice(Sol_dUf_acc,Sol1,1,NodG, iShift_dUf_acc + 1 , iShift_dUf_acc + NdimE, iDofT)
	call getSlice(Xel,XLL,1,NodG,1 ,NdimE, Ndim)

	
	sizeSolU = NodG*NdimE
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
	
	AE =  0.0d0
	BE =  0.0d0
	
	Auu = 0.0d0
	Bu = 0.0d0
	
	Gbar = 0.0d0
	Gbar(2,2) = 1.0d0
	
	SolGy = 0.0d0
	
	lambda = Sol1(NodG*idofT + iShift_LarcLength + 1)
	dlambda_acc = Sol1(NodG*idofT + iShift_dLarcLength_acc + 1)
	
	do i = 1 , NodG
		y = Xel((i-1)*NdimE + 1 : (i-1)*NdimE + 2)
		SolGy((i-1)*NdimE + 1 : (i-1)*NdimE + 2) = matmul(Gbar,y)
	end do
	
	Do nG = 1, NGP ! LoopGauss
		
		call PtG%calcGradU(GradU,SolUf(1:sizeSolU),nG)
				
!~ 		call PtG%getNmat(Nmat(:,1:sizeSolU),nG)

!~ 		dU_acc = matmul(Nmat, dLarcLength_acc*SolGy + Sol_dUf_acc ) !!! watch for dimensions
!~ 		SolGy_gp = matmul(Nmat, SolGy ) !!! watch for dimensions
		
!~ 		call getF(F,iFtype) !!! with identity summed
		
		F = deltaKron2 + lambda*Gbar + GradU 
		
		call SUPD_continuum(Ppk,D,F,damagePar,Param,matPar,constlaw,NDimE,damageMethod,isAnalytic)
		
		call PtG%getBmat(Bmat,nG)
		BmatT = transpose(Bmat)				
		call voigtTen4toTen2(Dmat,D)
		call voigtTen2toVec(Pmat,Ppk)
		
		dV=ptG%dV(nG)
		
		Auu = Auu + matmul(BmatT,matmul(Dmat,Bmat))*dV
		Aul(:,1) = Aul(:,1) + matmul(BmatT,matmul(Dmat,matmul(Bmat,SolGy)))*dV  !! different 
		
		Bu = Bu - matmul(BmatT,Pmat)*dV
		
		Bl  = Bl + arcLength*arcLength*dV 
	EndDo !LoopGauss
	
	Sol_dU_acc = dlambda_acc*SolGy + Sol_dUf_acc
	Alu(1,:) = 2.0d0*Sol_dU_acc
	A_ll(1,1) = 2.0d0*dot_product(Sol_dU_acc, SolGy)
	
	if(Nconverg > 1) then
		Bl = Bl - dot_product(Sol_dU_acc,Sol_dU_acc)
	end if
	
	if(Nstep==1 .and. Nconverg ==1) then
		A_ll(1,1) = 1.0d0
	end if


!~ 		Aul = 0.0d0
!~ 		Alu = 0.0d0
!~ 		A_ll = 1.0d0
!~ 		Bl = 1.0d-4*arcLength
!~ 	end if
	
	
!~ 	Bl = arcLength/(2.0d0**real(Nconverg-1))
	
!~ 	if(Nconverg>5) then
!~ 		Bl = 0.0d0
!~ 	end if
		
	call setBlockToVector(BE,Bu,NodG,iDofT, NdimE, iShift_dUf ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	call setBlockToMatrixSquare(AE,Auu, NodG,iDofT, NdimE, iShift_dUf ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A

	call setBlockToVector(BE,Bl,1,iDofT, 1, NodG*idofT + iShift_dLarcLength) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	call setBlockToMatrix(AE,Aul, NodG, 1, iDofT, iDofT, NdimE, 1, iShift_dUf , NodG*idofT + iShift_dLarcLength)
	call setBlockToMatrix(AE,Alu, 1 , NodG , iDofT, iDofT,  1, NdimE, NodG*idofT + iShift_dLarcLength, iShift_dUf)
	call setBlockToMatrix(AE,A_ll, 1 , 1 , iDofT, iDofT,  1, 1 , NodG*idofT + iShift_dLarcLength, NodG*idofT + iShift_dLarcLength)
	
 end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine arcLengthConsistentS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    use ptsGaussLib
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: iShift_dLarcLength, iShift_dUf, iFemtype, NodG
    integer , parameter :: NdimE = 2

    Coupling = 0
    
	iShift_dUf = nint(CommonPar(3))
	iShift_dLarcLength = nint(CommonPar(6))
	iFEMtype =  nint(commonPar(8)) 

    call setNodG(ifemtype, NodG)
    
    write(0,*) iAdd
    pause
    
	call setSimbolicBlockToMatrixSquare(Coupling,NodG,iDofT, NdimE, iShift_dUf )
	
	
	call setSimbolicBlockToMatrix(Coupling,NodG, 1, iDofT, iDofT, NdimE, 1, iShift_dUf , NodG*idofT + iShift_dLarcLength)
	call setSimbolicBlockToMatrix(Coupling,1 , NodG , iDofT, iDofT,  1, NdimE, NodG*idofT + iShift_dLarcLength, iShift_dUf)
	call setSimbolicBlockToMatrix(Coupling,1 , 1 , iDofT, iDofT,  1, 1 , &
	                              NodG*idofT + iShift_dLarcLength, NodG*idofT + iShift_dLarcLength)

end Subroutine
