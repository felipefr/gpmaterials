    !> Incorporates affine boundary for network 
    !! 
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iShiftLag = nint(CommonPar(2))
	!! @param pen  = commonPar(3)

    !! @author Rocha, Felipe Figueredo

!     ------------------------------------------------------------------
Subroutine affineBoundary(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
						Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux , only : numprint, setAEandBE_lagrangeMultipliers, getSliceAllocate
	use globalVariables, only : NdimE, getMaterial, maxMatPar, getAnisoTensorInv
	use fibresLib
	use fibresMod
		
	IMPLICIT NONE
	
	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	! Reals Vectors and matrices
	Real*8  :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), &
		Param(*), JParam(*), CommonPar(*)
	
	!   =====   END ARGUMENTS  =======
     		
	real*8 , allocatable ::  Xel(:) 
	Real*8 :: pen 
	integer :: iShiftU, iShiftDeltaU
	integer :: flag1 , flag2
	integer :: dofD
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	pen = commonPar(3)
	
	flag1 = nint(Param(Ipos_flag1))
	flag2 = nint(Param(Ipos_flag2))
	
	
	if(flag1 > 0) then
		dofD = iShiftDeltaU 
		AE(dofD + 1,dofD + 1) = pen
		AE(dofD + 2,dofD + 2) = pen
		if(iShiftU /= iShiftDeltaU) then 
!~ 			BE(dofD + 1) = -pen*Sol1(dofD + 1)
!~ 			BE(dofD + 2) = -pen*Sol1(dofD + 2)
			BE(dofD + 1) = -pen*Sol1(dofD + iShiftU - iShiftDeltaU + 1)
			BE(dofD + 2) = -pen*Sol1(dofD + iShiftU - iShiftDeltaU + 2)
		end if
	end if
	
	if(flag2 > 0) then
		dofD = iDofT + iShiftDeltaU
		AE(dofD + 1,dofD + 1) = pen
		AE(dofD + 2,dofD + 2) = pen
		if(iShiftU /= iShiftDeltaU) then 
!~ 			BE(dofD + 1) = -pen*Sol1(dofD + 1)
!~ 			BE(dofD + 2) = -pen*Sol1(dofD + 2)
			BE(dofD + 1) = -pen*Sol1(dofD + iShiftU - iShiftDeltaU + 1)
			BE(dofD + 2) = -pen*Sol1(dofD + iShiftU - iShiftDeltaU + 2)
		end if
	end if
		
			
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine affineBoundaryS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux , only : numprint, setCoupling_lagrangeMultipliers
	
	use globalVariables, only : NdimE
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
	integer ::  iShiftU, iShiftDeltaU 
	integer , parameter :: NodG = 2, NodLag = 3

	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))

	Coupling(iShiftDeltaU + 1,iShiftDeltaU + 1) = 1
	Coupling(iShiftDeltaU + 2,iShiftDeltaU + 2) = 1
	Coupling(iShiftDeltaU + iDofT + 1,iShiftDeltaU + iDofT + 1) = 1
	Coupling(iShiftDeltaU + iDofT + 2,iShiftDeltaU + iDofT + 2) = 1
	
end subroutine





