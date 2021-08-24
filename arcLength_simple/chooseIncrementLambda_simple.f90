	!> Element that solves the trivial diagonal problem of incrementing variables
	!! @param CommonPar = ( a , a*[bi,ci,di,ei] , dummy ) where : 
	!! @param Nvar = CommonPar(1) ==> a
	!! @param iShiftUi = CommonPar(see structure) ===> bi
	!! @param iShiftDeltaUi = CommonPar(see structure) ===> ci
	!! @param nodei (for each material = see structure) ===> ei 
	!! @author Rocha, Felipe Figueredo


!! last update: lenU given by parameter and not by difference between iShiftU and iShiftDeltaU. NodMax is given, rather than assumed
Subroutine chooseIncrementLambda_simple(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
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
    
    integer	:: iShift_dLarcLength, iShift_LarcLength, iShift_dLarcLength_acc, & 
			   iShift_c1, iShift_c2, iShift_c3, iShift_AngleEst1, iShift_AngleEst2
    real*8 :: dLarcLength, dLarcLength_acc, c1, c2, c3, arcLength, AngleEst1, AngleEst2 , dl1, dl2
	
	arcLength = CommonPar(1)
	iShift_dLarcLength = nint(CommonPar(2))
	iShift_LarcLength = nint(CommonPar(3))
	iShift_dLarcLength_acc = nint(CommonPar(4))	
	iShift_c1 = nint(CommonPar(5))
	iShift_c2 = nint(CommonPar(6))
	iShift_c3 = nint(CommonPar(7))
	iShift_AngleEst1 = nint(CommonPar(8))
	iShift_AngleEst2 = nint(CommonPar(9))
	
	c1 = Sol1(iShift_c1 + 1)
	c2 = Sol1(iShift_c2 + 1)
	c3 = Sol1(iShift_c3 + 1)
	AngleEst1 = Sol1(iShift_AngleEst1 + 1)
	AngleEst2 = Sol1(iShift_AngleEst2 + 1)

	dl1 = solveQuadraticEquation(c1,c2,c3,1)
	dl2 = solveQuadraticEquation(c1,c2,c3,-1)
	
	if(Nconverg == 1) then
		if(AngleEst1<0.0) then 
			dLarcLength = - dabs(dl1)
		else
			dLarcLength = dabs(dl1)
		end if
	else if(AngleEst1>AngleEst2) then
		dLarcLength = dl1
	else 
		dLarcLength = dl2 
	end if
	
!~ 	write(0,*) 'DEBUUUUG = ', Nconverg , dl1, dl2, dLarcLength
!~ 	Pause
	
	AE(iShift_dLarcLength + 1, iShift_dLarcLength + 1) = 1.0d0
	AE(iShift_LarcLength + 1, iShift_LarcLength + 1) = 1.0d0
	AE(iShift_dLarcLength_acc + 1, iShift_dLarcLength_acc + 1) = 1.0d0
	
	BE(iShift_dLarcLength + 1) = dLarcLength
	BE(iShift_LarcLength + 1) = Sol1(iShift_LarcLength + 1) + dLarcLength
	BE(iShift_dLarcLength_acc + 1) = Sol1(iShift_dLarcLength_acc + 1) + dLarcLength
	
end subroutine 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine chooseIncrementLambda_simpleS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
    
    integer	:: iShift_dLarcLength, iShift_LarcLength, iShift_dLarcLength_acc
    
	iShift_dLarcLength = nint(CommonPar(2))
	iShift_LarcLength = nint(CommonPar(3))
	iShift_dLarcLength_acc = nint(CommonPar(4))

	Coupling(iShift_dLarcLength + 1, iShift_dLarcLength + 1) = 1
	Coupling(iShift_LarcLength + 1, iShift_LarcLength + 1) = 1
	Coupling(iShift_dLarcLength_acc + 1, iShift_dLarcLength_acc + 1) = 1

!~ 	call numprint(Coupling)
!~ 	pause
end Subroutine

