    !> Generic element for finite strains
    !!
    !!
	!! @param iShift_Uf = nint(CommonPar(1))
	!! @param iShift_dUf = nint(CommonPar(2))
	!! @param iFemType  = nint(CommonPar(3)) 
	!! @param iFType  = nint(commonPar(4))
	!! @param iMaterial  = nint(commonPar(5))
	!! @param iDamageParType  = nint(commonPar(6))
	!!
    !! @author Rocha, Felipe Figueredo

Subroutine arcLengthInconsistent(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : getF, NdimE, getDamagePar, getMaterial, maxDamagePar, maxMatPar
	use ptsGaussLib
	use genericConstitutiveLib , only : SUPD_continuum
	use TIME_STEP_ITER
		
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
         
	Integer :: nG, i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ipDim, ipDimE, ipDof,ipDofmod ! counters 
	Integer :: constLaw ! if NparamDamage = 0 , no damage
	Real*8 :: dV , Det , MatPar(maxMatPar), damagePar(maxDamagePar)
	Real*8 :: Grad_Uf(NdimE,NdimE) , F(NdimE,NdimE)
	Real*8 :: D(NdimE,NdimE,NdimE,NdimE), Ppk(NdimE,NdimE), Dmat(NdimE*NdimE,NdimE*NdimE) , Pmat(NdimE*NdimE)
	Real*8 :: Bmat(NdimE*NdimE,NodElt*NdimE), BmatT(NodElt*NdimE,NdimE*NdimE), Auu(NodElt*NdimE,NodElt*NdimE), Bu(NodElt*NdimE)
	real*8 , allocatable, dimension(:)::  Sol_Uf , Xel 
	integer :: pOrder , NGP, NodG, iFEMtype , iFtype, iMaterial, iSimplex, iBubble, lU, iDamageParType, & 
				damageMethod,isAnalytic
	integer :: iShift_Uf, iShift_dUf, iShift_LarcLength, computeStar_or_Bar
	type(ptGaussClass) :: PtG
	 
	real*8 :: Gbar(NdimE,NdimE),  Gbar_vec(NdimE*NdimE), lambda, y(NdimE)
	
	iShift_Uf = nint(CommonPar(1))
	iShift_dUf = nint(CommonPar(2))
	iShift_LarcLength = nint(CommonPar(3))
	iFEMtype =  nint(commonPar(4)) 
	iFtype = nint(commonPar(5)) 
	iMaterial = nint(commonPar(6))
	iDamageParType = nint(commonPar(7))
	damageMethod = nint(commonPar(8))
	isAnalytic = nint(commonPar(9))
	computeStar_or_Bar = nint(commonPar(10))

	if(iDamageParType>0) call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw, matPar, iMaterial)
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)

	call getSliceAllocate(Sol_Uf,Sol1,1,NodG,iShift_Uf + 1 ,iShift_Uf + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	lambda = Sol1(NodG*idofT + iShift_LarcLength + 1)
	
	lU = NodG*NdimE
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
	
	AE =  0.0d0
	BE =  0.0d0
	
	Auu = 0.0d0
	Bu = 0.0d0
	
	Gbar = 0.0d0
	Gbar(2,2) = 1.0d0
		
	Do nG = 1, NGP ! LoopGauss
		
		call PtG%calcGradU(Grad_Uf,Sol_Uf(1:lU),nG)
		
		F = deltaKron2 + lambda*Gbar + Grad_Uf 
		
		call SUPD_continuum(Ppk,D,F,damagePar,Param,matPar,constlaw,NDimE,damageMethod,isAnalytic)
		
		call PtG%getBmat(Bmat,nG)
		BmatT = transpose(Bmat)				
		call voigtTen4toTen2(Dmat,D)
		call voigtTen2toVec(Pmat,Ppk)
		call voigtTen2toVec(Gbar_vec,Gbar)
		
		dV=ptG%dV(nG)
		
		if(computeStar_or_Bar == 1) then !! Star
			if(Nconverg == 1) then
				Bu = 0.0d0
			else
				Bu(1:lU) = Bu(1:lU) - matmul(BmatT(1:lU,:),Pmat)*dV
			end if
		else if(computeStar_or_Bar == 2) then !! Bar
			Bu(1:lU) = Bu(1:lU) - matmul(BmatT(1:lU,:),matmul(Dmat,Gbar_vec))*dV  !! different 
		end if
		
		Auu(1:lU,1:lU) = Auu(1:lU,1:lU) + matmul(BmatT(1:lU,:),matmul(Dmat,Bmat(:,1:lU)))*dV
		
	EndDo !LoopGauss
	
	call setBlockToVector(BE,Bu,NodG,iDofT, NdimE, iShift_dUf ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	call setBlockToMatrixSquare(AE,Auu, NodG,iDofT, NdimE, iShift_dUf ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine arcLengthInconsistentS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    use ptsGaussLib , only : setNodG
	use globalVariables, only : NdimE 
	
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol , iShift_dUf, iFemtype, NodG, flagIncrementInside

    Coupling = 0
    
    iShift_dUf = nint(commonPar(2))
    iFemType = nint(commonPar(4)) 
    
    call setNodG(iFemtype, NodG)
    
	call setSimbolicBlockToMatrixSquare(Coupling,NodG,iDofT, NdimE, iShift_dUf )


end Subroutine
