
!     ------------------------------------------------------------------
Subroutine NodalForce(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
						Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use loadingLib
	use finiteStrainLib
	use ptsGaussLib
		
	IMPLICIT NONE
		
	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	! Reals Vectors and matrices
	Real*8  :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), &
		Param(*), JParam(*), CommonPar(*)
	
	!   =====   END ARGUMENTS  =======
     		
	Integer :: iShiftU, loadType
	Real*8 :: force, theta
	
	iShiftU = nint(CommonPar(1))
	theta = CommonPar(2)
	loadType = nint(CommonPar(3))
	
	theta = theta*PI
	call chooseLoad(loadType)
	call pressureLoad(force,Time,DelT,CommonPar(4:9))

	BE(iShiftU + 1) = force*dcos(theta) !! in x
	BE(iShiftU + 2) = force*dsin(theta) !! in y

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine NodalForceS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
		
	integer :: iShiftU
	
	iShiftU = nint(CommonPar(1))
	
	Coupling(iShiftU + 1, iShiftU + 1) = 1
	Coupling(iShiftU + 2, iShiftU + 2) = 1
	
end Subroutine



