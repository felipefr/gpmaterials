	!> Element that solves the trivial diagonal problem of incrementing variables
	!! @param CommonPar = ( a , a*[bi,ci,di,ei] , dummy ) where : 
	!! @param Nvar = CommonPar(1) ==> a
	!! @param iShiftUi = CommonPar(see structure) ===> bi
	!! @param iShiftDeltaUi = CommonPar(see structure) ===> ci
	!! @param nodei (for each material = see structure) ===> ei 
	!! @author Rocha, Felipe Figueredo


!! last update: lenU given by parameter and not by difference between iShiftU and iShiftDeltaU. NodMax is given, rather than assumed
Subroutine computeIncrementEstimation_simple(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
									Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use TIME_STEP_ITER
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
    integer :: ip, NodG, iShift_c1, iShift_c2, iShift_c3, iShift_da_Star, &
				iShift_da_Bar, iShift_da_acc, iShift_dLarcLength_acc, & 
				iShift_AngleEst1, iShift_AngleEst2
				
	real*8 :: c1,c2,c3, dLarcLength_acc, da_bar(2) , da_Star(2), Da_acc(2), da(2), &
	          Da_acc_old(2), dlambdaIncrement1, dlambdaIncrement2, AngleEst1, AngleEst2

	iShift_c1 = nint(CommonPar(1))
	iShift_c2 = nint(CommonPar(2))
	iShift_c3 = nint(CommonPar(3))
	iShift_da_Star = nint(CommonPar(4))	
	iShift_da_Bar = nint(CommonPar(5))	
	iShift_da_acc = nint(CommonPar(6))	
	iShift_dLarcLength_acc = nint(CommonPar(7))	
	iShift_AngleEst1 = nint(CommonPar(8))
	iShift_AngleEst2 = nint(CommonPar(9))

	NodG = 2
	
	ip = NodG*idofT + 1
	c1 = Sol1(iShift_c1 + ip)
	c2 = Sol1(iShift_c2 + ip)
	c3 = Sol1(iShift_c3 + ip)
	dLarcLength_acc = Sol1(iShift_dLarcLength_acc + ip)

	dlambdaIncrement1 = solveQuadraticEquation(c1,c2,c3,1) !! choose that have to change
	dlambdaIncrement2 = solveQuadraticEquation(c1,c2,c3,-1) !! choose that have to change
	
	da_bar(1) = Sol1(iShift_da_Bar + 1)
	da_bar(2) = Sol1(iDofT + iShift_da_Bar + 1)
	da_star(1) = Sol1(iShift_da_Star + 1)
	da_star(2) = Sol1(iDofT + iShift_da_Star + 1)
	Da_acc(1) = Sol1(iShift_Da_acc + 1)
	Da_acc(2) = Sol1(iDofT + iShift_Da_acc + 1)
	Da_acc_old(1) = Sol0(iShift_Da_acc + 1)
	Da_acc_old(2) = Sol0(iDofT + iShift_Da_acc + 1)
	
	
	AngleEst1 = 0.0d0
	AngleEst2 = 0.0d0
	
	
	if(Nconverg == 1) then ! in this case compute estimative for the signal in AngleEst1
		AngleEst1 = dot_product(da_bar, Da_acc_old)
	else	
		da = da_Star + dlambdaIncrement1*da_Bar
		AngleEst1 = dot_product(Da_acc + da, Da_acc) ! it is possible to be the old
	 
		da = da_Star + dlambdaIncrement2*da_Bar
		AngleEst2 = dot_product(Da_acc + da, Da_acc)  ! it is possible to be the old
	
	end if
	ip = NodG*idofT + 1
	AE(iShift_AngleEst1 + ip, iShift_AngleEst1 + ip) = 1.0d0
	AE(iShift_AngleEst2 + ip, iShift_AngleEst2 + ip) = 1.0d0
		
	BE(iShift_AngleEst1 + ip) = AngleEst1
	BE(iShift_AngleEst2 + ip) = AngleEst2 
	
end subroutine 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine computeIncrementEstimation_simpleS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
    
    integer	:: iShift_AngleEst1, iShift_AngleEst2, ip, NodG
    
	iShift_AngleEst1 = nint(CommonPar(8))
	iShift_AngleEst2 = nint(CommonPar(9))

	NodG = 2
	
	ip = NodG*idofT + 1
	Coupling(iShift_AngleEst1 + ip, iShift_AngleEst1 + ip) = 1
	Coupling(iShift_AngleEst2 + ip, iShift_AngleEst2 + ip) = 1
		
end Subroutine

