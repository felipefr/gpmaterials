	!> Element that solves the trivial diagonal problem of incrementing variables
	!! @param CommonPar = ( a , a*[bi,ci,di,ei] , dummy ) where : 
	!! @param Nvar = CommonPar(1) ==> a
	!! @param iShiftUi = CommonPar(see structure) ===> bi
	!! @param iShiftDeltaUi = CommonPar(see structure) ===> ci
	!! @param nodei (for each material = see structure) ===> ei 
	!! @author Rocha, Felipe Figueredo


!! last update: lenU given by parameter and not by difference between iShiftU and iShiftDeltaU. NodMax is given, rather than assumed
Subroutine computeIncrementEstimation(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
									Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------

	use funcAux
	use globalVariables, only : NdimE
	use ptsGaussLib
	use TIME_STEP_ITER
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
    integer :: i, ip
	integer :: 	iShift_c1, iShift_c2 , iShift_c3 , iShift_dUf_Star, iShift_dUf_Bar, &
				iShift_DUf_acc, iShift_dLarcLength_acc, iShift_AngleEst1 , iShift_AngleEst2, &
				iFEMtype 
	integer :: NodG,pOrder,NGP,iSimplex,iBubble, nG, sizeSolU
	real*8 :: c1,c2,c3, dLarcLength_acc, dlambdaIncrement1 , dlambdaIncrement2, AngleEst1, AngleEst2 
	real*8, allocatable, dimension(:) :: Sol_DUf_acc, Sol_DUf_acc_old, Sol_dUf_Bar, Sol_dUf_Star, Sol_dUf , Xel
	real*8 :: Gbar(NdimE,NdimE), y(NdimE), Nmat(NdimE, NdimE*NodElt), dV
	
	type(ptGaussClass) :: PtG

	iShift_c1 = nint(CommonPar(1))
	iShift_c2 = nint(CommonPar(2))
	iShift_c3 = nint(CommonPar(3))
	iShift_dUf_Star = nint(CommonPar(4))	
	iShift_dUf_Bar = nint(CommonPar(5))	
	iShift_DUf_acc = nint(CommonPar(6))	
	iShift_dLarcLength_acc = nint(CommonPar(7))	
	iShift_AngleEst1 = nint(CommonPar(8))
	iShift_AngleEst2 = nint(CommonPar(9))
	iFEMtype = nint(CommonPar(10))

	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)
	
	ip = NodG*idofT + 1
	c1 = Sol1(iShift_c1 + ip)
	c2 = Sol1(iShift_c2 + ip)
	c3 = Sol1(iShift_c3 + ip)
	dLarcLength_acc = Sol1(iShift_dLarcLength_acc + ip)

	dlambdaIncrement1 = solveQuadraticEquation(c1,c2,c3,1) !! choose that have to change
	dlambdaIncrement2 = solveQuadraticEquation(c1,c2,c3,-1) !! choose that have to change
	
	call getSliceAllocate(Sol_DUf_acc,Sol1,1,NodG, iShift_DUf_acc + 1 , iShift_DUf_acc + NdimE, iDofT)
	call getSliceAllocate(Sol_DUf_acc_old,Sol0,1,NodG, iShift_DUf_acc + 1 , iShift_DUf_acc + NdimE, iDofT)
	call getSliceAllocate(Sol_dUf_Star,Sol1,1,NodG, iShift_dUf_Star + 1 , iShift_dUf_Star + NdimE, iDofT)
	call getSliceAllocate(Sol_dUf_Bar,Sol1,1,NodG, iShift_dUf_Bar + 1 , iShift_dUf_Bar + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	allocate(Sol_dUf(NodG*NdimE))
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
	
	Gbar = 0.0d0
	Gbar(2,2) = 1.0d0
	
	sizeSolU = NodG*NdimE

	AngleEst1 = 0.0d0
	AngleEst2 = 0.0d0

	if(Nconverg == 1) then ! in this case compute estimative for the signal in AngleEst1
		AngleEst1 = dot_product(Sol_dUf_Bar, Sol_DUf_acc_old)
		AngleEst2 = 0.0d0
	else
		Sol_dUf = Sol_dUf_Star + dlambdaIncrement1*Sol_dUf_Bar
		AngleEst1 = dot_product(Sol_DUf_acc + Sol_dUf, Sol_DUf_acc) ! it is possible to be the old
		
		Sol_dUf = Sol_dUf_Star + dlambdaIncrement2*Sol_dUf_Bar
		AngleEst2 = dot_product(Sol_DUf_acc + Sol_dUf, Sol_DUf_acc)   ! it is possible to be the old
		
	end if
	
	ip = NodG*idofT + 1
	AE(iShift_AngleEst1 + ip, iShift_AngleEst1 + ip) = 1.0d0
	AE(iShift_AngleEst2 + ip, iShift_AngleEst2 + ip) = 1.0d0
	
	BE(iShift_AngleEst1 + ip) = AngleEst1
	BE(iShift_AngleEst2 + ip) = AngleEst2 
	
end subroutine 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine computeIncrementEstimationS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use ptsGaussLib, only : setNodG
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
    
    integer	:: iShift_AngleEst1, iShift_AngleEst2, ip, NodG, iFEMtype
    
	iShift_AngleEst1 = nint(CommonPar(8))
	iShift_AngleEst2 = nint(CommonPar(9))
	iFEMtype = nint(CommonPar(10))

	call setNodG(iFEMtype,NodG)
	
	ip = NodG*idofT + 1
	Coupling(iShift_AngleEst1 + ip, iShift_AngleEst1 + ip) = 1
	Coupling(iShift_AngleEst2 + ip, iShift_AngleEst2 + ip) = 1
		
!~ 	call numprint(Coupling)
!~ 	pause
end Subroutine

