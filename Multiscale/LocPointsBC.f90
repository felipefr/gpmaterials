Subroutine LocPointsBC &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)

	use funcAux
	
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
    integer, parameter :: NdimE = 2
    integer :: i , ipLag, ipDeltaU, iShiftLag, iShiftDeltaU, iShiftFlagLoc
    Real*8 ::  pen, eps
    logical :: locActive
   
	iShiftDeltaU = nint(commonPar(1))
	iShiftLag = nint(commonPar(2))
	iShiftFlagLoc = nint(CommonPar(3)) ! in geometrical nodes
	
	pen = 1.0d5
	eps = 1.0d-10
	
	locActive = .false.
!~ 	write(*,*) Sol0( iShiftFlagLoc + 1)
	
	if(Sol0( iShiftFlagLoc + 1) > 0.0d0) locActive = .true.
	
	do i = 1 , NdimE
		ipLag = iShiftLag + i
		ipDeltaU = iShiftDeltaU + i

		if(locActive) then
			AE(ipLag, ipDeltaU) = pen
			AE(ipDeltaU , ipLag ) = 1.0d0
!~ 			AE(ipDeltaU, ipDeltaU) = pen
			AE(ipLag , ipLag ) = eps
			BE(ipLag ) = eps*Sol1(ipLag )
!~ 				BE(ipLag) = 0.0d0

!~ 			write(*,*) "is active"
		else 
			AE(ipLag , ipLag) = eps
!~ 			write(*,*) "is not active"
		end if
	end do

end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine LocPointsBCS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
	
    implicit none
    
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer ::  Coupling(MaxLRows,MaxLRows)
    Real*8 :: CommonPar(*), eps
    integer, parameter :: NdimE = 2
    integer :: iShiftLag, iShiftDeltaU, i,  ipLag, ipDeltaU

	iShiftDeltaU = nint(commonPar(1))
	iShiftLag = nint(commonPar(2))

	do i = 1 , NdimE
		ipLag = iShiftLag + i
		ipDeltaU = iShiftDeltaU + i

		Coupling(ipLag,ipLag) = 1
		Coupling(ipDeltaU,ipLag) = 1
		Coupling(ipLag,ipDeltaU) = 1
!~ 		Coupling(ipDeltaU,ipDeltaU) = 1
		
	end do
	
	
end Subroutine

