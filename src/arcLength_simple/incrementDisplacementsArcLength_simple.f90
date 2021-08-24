	!> Element that solves the trivial diagonal problem of incrementing variables
	!! @param CommonPar = ( a , a*[bi,ci,di,ei] , dummy ) where : 
	!! @param Nvar = CommonPar(1) ==> a
	!! @param iShiftUi = CommonPar(see structure) ===> bi
	!! @param iShiftDeltaUi = CommonPar(see structure) ===> ci
	!! @param nodei (for each material = see structure) ===> ei 
	!! @author Rocha, Felipe Figueredo


!! last update: lenU given by parameter and not by difference between iShiftU and iShiftDeltaU. NodMax is given, rather than assumed
Subroutine incrementDisplacementsArcLength_simple(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
									Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
    integer :: NodG
    integer	:: iShift_a, iShift_Da_acc, iShift_da_Bar, iShift_da_Star, iShift_dLarcLength
	
    real*8 :: dl, a(2) , da_bar(2) , da_Star(2), Da_acc(2), da(2)
    
    
	iShift_dLarcLength = nint(CommonPar(1))
	iShift_a = nint(CommonPar(2))
	iShift_Da_acc = nint(CommonPar(3))
	iShift_da_Star = nint(CommonPar(4))	
	iShift_da_Bar = nint(CommonPar(5))
	
	NodG = 2 
	
	dl = Sol1(NodG*iDofT + iShift_dLarcLength + 1)
	
	a(1) = Sol1(iShift_a + 1)
	a(2) = Sol1(iDofT + iShift_a + 1)
	Da_acc(1) = Sol1(iShift_Da_acc + 1)
	Da_acc(2) = Sol1(iDofT + iShift_Da_acc + 1)
	da_bar(1) = Sol1(iShift_da_Bar + 1)
	da_bar(2) = Sol1(iDofT + iShift_da_Bar + 1)
	da_star(1) = Sol1(iShift_da_Star + 1)
	da_star(2) = Sol1(iDofT + iShift_da_Star + 1)
	
	da = da_Star + dl*da_Bar
	
!~ 	write(0,*) 'DEBUG = ' , da(1) , da(2)
!~ 	pause
	
	AE(iShift_a + 1,iShift_a + 1) = 1.0d0
	AE(iDofT + iShift_a + 1,iDofT + iShift_a + 1) = 1.0d0
	AE(iShift_Da_acc + 1,iShift_Da_acc + 1) = 1.0d0
	AE(iDofT + iShift_Da_acc + 1,iDofT + iShift_Da_acc + 1) = 1.0d0
	
	BE(iShift_a + 1) = a(1) + da(1)
	BE(iDofT + iShift_a + 1) = a(2) + da(2)
	BE(iShift_Da_acc + 1) = Da_acc(1) + da(1)
	BE(iDofT + iShift_Da_acc + 1) = Da_acc(2) + da(2)

end subroutine 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine incrementDisplacementsArcLength_simpleS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux
	use globalVariables, only : NdimE
	use ptsGaussLib, only : setNodG
	
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
    
    integer	:: iShift_a, iShift_Da_acc
    
	iShift_a = nint(CommonPar(2))
	iShift_Da_acc = nint(CommonPar(3))
			
	Coupling(iShift_a + 1,iShift_a + 1) = 1
	Coupling(iDofT + iShift_a + 1,iDofT + iShift_a + 1) = 1
	Coupling(iShift_Da_acc + 1,iShift_Da_acc + 1) = 1
	Coupling(iDofT + iShift_Da_acc + 1,iDofT + iShift_Da_acc + 1) = 1

end Subroutine

