	!> Element that solves the trivial diagonal problem of incrementing variables
	!! @param CommonPar = ( a , a*[bi,ci,di,ei] , dummy ) where : 
	!! @param Nvar = CommonPar(1) ==> a
	!! @param iShiftUi = CommonPar(see structure) ===> bi
	!! @param iShiftDeltaUi = CommonPar(see structure) ===> ci
	!! @param nodei (for each material = see structure) ===> ei 
	!! @author Rocha, Felipe Figueredo


!! last update: lenU given by parameter and not by difference between iShiftU and iShiftDeltaU. NodMax is given, rather than assumed
Subroutine incrementLambda(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
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
    
    integer	:: iShift_dLarcLength, iShift_LarcLength, iShift_LarcLength_acc 
    real*8 :: lambdaIncrement, LarcLength_acc
	
	lambdaIncrement = CommonPar(1)
	iShift_dLarcLength = nint(CommonPar(2))
	iShift_LarcLength = nint(CommonPar(3))
	iShift_LarcLength_acc = nint(CommonPar(4))	

	lambdaIncrement = lambdaIncrement/(2.0d0**real(Nconverg-1))
	
	if(Nconverg>5) then
		lambdaIncrement = 0.0d0
	end if
	
	if(Nconverg == 1) then
		LarcLength_acc = lambdaIncrement
	else 
		LarcLength_acc = Sol1(iShift_LarcLength_acc + 1) + lambdaIncrement
	end if
	
	
	AE(iShift_dLarcLength + 1, iShift_dLarcLength + 1) = 1.0d0
	AE(iShift_LarcLength + 1, iShift_LarcLength + 1) = 1.0d0
	AE(iShift_LarcLength_acc + 1, iShift_LarcLength_acc + 1) = 1.0d0
	
	
	BE(iShift_dLarcLength + 1) = lambdaIncrement
	BE(iShift_LarcLength + 1) = Sol1(iShift_LarcLength + 1) + lambdaIncrement
	BE(iShift_LarcLength_acc + 1) = LarcLength_acc
	
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Subroutine incrementLambdaS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
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

