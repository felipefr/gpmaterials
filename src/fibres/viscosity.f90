    !> Incorporates affine boundary for network 
    !! 
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iShiftLag = nint(CommonPar(2))
	!! @param pen  = commonPar(3)

    !! @author Rocha, Felipe Figueredo

!     ------------------------------------------------------------------
Subroutine Viscosity(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
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
     
    integer, parameter ::  NodG = 2
	real*8 , allocatable ::  Xel(:) 
	Real*8 :: gamma1, gamma2
	integer :: iShiftU, iShiftDeltaU
	integer :: flag1 , flag2
	integer :: i, ip, j , i_ini, i_end, opBoundary
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	gamma1 = commonPar(3) !! numeric viscosity 
	gamma2 = commonPar(4) !! real viscosity 
	
	
	opBoundary = 0 !! viscosity in every node
!~ 	opBoundary = 1 !! viscosity just in Boundary
	
	!! default
	i_ini = 1
	i_end = NodG
	
	if(opBoundary == 1) then
		flag1 = nint(Param(Ipos_flag1))
		flag2 = nint(Param(Ipos_flag2))
	
		if(flag1 == 1 .and. flag2 == 0) then
			i_ini = 1
			i_end = 1
		else if(flag2 == 1 .and. flag1 == 0) then
			i_ini = 2
			i_end = 2
		else if(flag1 == 1 .and. flag2 == 1) then
			i_ini = 1
			i_end = 2 
		else
			return
		end if 
		
	end if
	
	
	do i = i_ini , i_end
		ip = (i-1)*iDofT + iShiftDeltaU
		do j = 1 , NdimE
			AE(ip + j, ip + j ) = gamma1 + gamma2 
			if(iShiftU == iShiftDeltaU) then !!! no increment
				BE(ip + j) = gamma1*Sol1(ip + j) + gamma2*Sol0(ip + j)
			else 
				BE(ip + j) = gamma2*(Sol0(ip + iShiftU - iShiftDeltaU + j) - Sol1(ip + iShiftU - iShiftDeltaU + j))
			end if
		end do 
	end do
				
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine ViscosityS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux , only : numprint, setCoupling_lagrangeMultipliers
	
	use globalVariables, only : NdimE
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
	integer ::  iShiftU, iShiftDeltaU 
	integer , parameter :: NodG = 2, NodLag = 3
	integer :: ip, i, j 
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))

	do i = 1 , NodG
		ip = (i-1)*iDofT + iShiftDeltaU
		do j = 1 , NdimE
			Coupling(ip + j, ip + j ) = 1 
		end do 
	end do
	
end subroutine





