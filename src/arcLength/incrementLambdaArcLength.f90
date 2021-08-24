	!> Element that solves the trivial diagonal problem of incrementing variables
	!! @param CommonPar = ( a , a*[bi,ci,di,ei] , dummy ) where : 
	!! @param Nvar = CommonPar(1) ==> a
	!! @param iShiftUi = CommonPar(see structure) ===> bi
	!! @param iShiftDeltaUi = CommonPar(see structure) ===> ci
	!! @param nodei (for each material = see structure) ===> ei 
	!! @author Rocha, Felipe Figueredo


!! last update: lenU given by parameter and not by difference between iShiftU and iShiftDeltaU. NodMax is given, rather than assumed
Subroutine incrementLambdaArcLength(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
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
    
    integer	:: iShift_dLarcLength, iShift_LarcLength, iShift_dLarcLength_acc, iShift_a, iShift_b, iShift_c
    real*8 :: lambdaIncrement, dLarcLength_acc, a, b, c, arcLength
	
	arcLength = CommonPar(1)
	iShift_dLarcLength = nint(CommonPar(2))
	iShift_LarcLength = nint(CommonPar(3))
	iShift_dLarcLength_acc = nint(CommonPar(4))	
	iShift_a = nint(CommonPar(5))
	iShift_b = nint(CommonPar(6))
	iShift_c = nint(CommonPar(7))
	

	a = Sol1(iShift_a + 1)
	b = Sol1(iShift_b + 1)
	c = Sol1(iShift_c + 1)

	if(Nconverg == 1) then
		lambdaIncrement = arcLength/dsqrt(a)
		dLarcLength_acc = lambdaIncrement  
	else 
		lambdaIncrement = solveQuadraticEquation(a,b,c,-1) !! choose that have to change
		dLarcLength_acc = Sol1(iShift_dLarcLength_acc + 1) + lambdaIncrement
	end if
	
	
	AE(iShift_dLarcLength + 1, iShift_dLarcLength + 1) = 1.0d0
	AE(iShift_LarcLength + 1, iShift_LarcLength + 1) = 1.0d0
	AE(iShift_dLarcLength_acc + 1, iShift_dLarcLength_acc + 1) = 1.0d0
	
	
	BE(iShift_dLarcLength + 1) = lambdaIncrement
	BE(iShift_LarcLength + 1) = Sol1(iShift_LarcLength + 1) + lambdaIncrement
	BE(iShift_dLarcLength_acc + 1) = dLarcLength_acc
	
end subroutine 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine incrementLambdaArcLengthS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
    
    integer	:: iShift_dLarcLength, iShift_LarcLength, iShift_LarcLength_acc
    
	iShift_dLarcLength = nint(CommonPar(2))
	iShift_LarcLength = nint(CommonPar(3))
	iShift_LarcLength_acc = nint(CommonPar(4))

	Coupling(iShift_dLarcLength + 1, iShift_dLarcLength + 1) = 1
	Coupling(iShift_LarcLength + 1, iShift_LarcLength + 1) = 1
	Coupling(iShift_LarcLength_acc + 1, iShift_LarcLength_acc + 1) = 1

!~ 	call numprint(Coupling)
!~ 	pause
end Subroutine

