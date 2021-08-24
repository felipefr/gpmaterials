	!> Element that solves the trivial diagonal problem of incrementing variables
	!! @param CommonPar = ( a , a*[bi,ci,di,ei] , dummy ) where : 
	!! @param Nvar = CommonPar(1) ==> a
	!! @param iShiftUi = CommonPar(see structure) ===> bi
	!! @param iShiftDeltaUi = CommonPar(see structure) ===> ci
	!! @param nodei (for each material = see structure) ===> ei 
	!! @author Rocha, Felipe Figueredo


!! last update: lenU given by parameter and not by difference between iShiftU and iShiftDeltaU. NodMax is given, rather than assumed
Subroutine zeroVariable_arcLength_simple(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
									Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------

	use TIME_STEP_ITER
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
    integer :: iShift_Da_acc, iShift_Dl_acc, NodG
    
	iShift_Da_acc = nint(CommonPar(1))
	iShift_Dl_acc = nint(CommonPar(2))
	
	NodG = 2 
	
	AE(iShift_Da_acc + 1, iShift_Da_acc + 1) = 1.0d0
	AE(iDofT + iShift_Da_acc + 1, iDofT + iShift_Da_acc + 1) = 1.0d0
	AE(NodG*iDofT + iShift_Dl_acc + 1 , NodG*iDofT + iShift_Dl_acc + 1) = 1.0d0

	if(Nconverg == 1) then
		BE(iShift_Da_acc + 1) = 0.0d0
		BE(iDofT + iShift_Da_acc + 1) = 0.0d0
		BE(NodG*iDofT + iShift_Dl_acc + 1) = 0.0d0
	else 
		BE(iShift_Da_acc + 1) = Sol1(iShift_Da_acc + 1)
		BE(iDofT + iShift_Da_acc + 1) = Sol1(iDofT + iShift_Da_acc + 1)
		BE(NodG*iDofT + iShift_Dl_acc + 1) = Sol1(NodG*iDofT + iShift_Dl_acc + 1)
	end if
	
end subroutine 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine zeroVariable_arcLength_simpleS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
    
    integer :: iShift_Da_acc, iShift_Dl_acc, NodG
    
	iShift_Da_acc = nint(CommonPar(1))
	iShift_Dl_acc = nint(CommonPar(2))
	
	NodG = 2
	
	Coupling(iShift_Da_acc + 1, iShift_Da_acc + 1) = 1
	Coupling(iDofT + iShift_Da_acc + 1, iDofT + iShift_Da_acc + 1) = 1
	Coupling(NodG*iDofT + iShift_Dl_acc + 1 , NodG*iDofT + iShift_Dl_acc + 1) = 1
	
end Subroutine

