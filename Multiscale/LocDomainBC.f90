Subroutine LocDomainBC &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)

	use funcAux
	
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
    integer, parameter :: NdimE = 2
    integer :: i , j , ipLag, ipDeltaU, iShiftLag, iShiftDeltaU, PosBifDomain, PosDetQ
    Real*8 ::  pen,lamb, eps, bValue
    logical :: locActive

	iShiftDeltaU = nint(commonPar(1))
	iShiftLag = nint(commonPar(2))
	PosBifDomain = nint(CommonPar(3)) ! in the Param
!~ 	PosDetQ = nint(CommonPar(4)) ! in the last node
	
	pen = 1.0d0
	eps = 1.0d-8
	
	locActive = .false.
	
	if( Param(PosBifDomain) > 0.0d0) then
		locActive = .true.
	end if
		
	locActive = .false.
	
	do i = 1 , NodElt - 1
		ipLag = (i-1)*iDofT + iShiftLag
		ipDeltaU = (i-1)*iDofT + iShiftDeltaU
		do j = 1 , NdimE
			if(locActive) then
				AE(ipLag + j , ipDeltaU + j) = pen
				AE(ipDeltaU + j , ipLag + j) = pen
!~ 				AE(ipDeltaU + j , ipDeltaU + j) = pen
				AE(ipLag + j , ipLag + j) = eps
				BE(ipLag + j) = eps*Sol1(ipLag + j)
!~ 				BE(ipLag + j) = 0.0d0
			else 
				AE(ipLag + j , ipLag + j) = eps
			end if
		end do
	end do

!~ 	write(*,*) locActive
!~ 	call numprint(Sol1(iShiftDeltaU + 1 : (NodElt-1)*idofT : idofT))
!~ 	call numprint(Sol1(iShiftDeltaU + 2 : (NodElt-1)*idofT : idofT))
!~ 	
!~ 	if(locActive) pause
	
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine LocDomainBCS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
	
    implicit none
    
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer ::  Coupling(MaxLRows,MaxLRows)
    Real*8 :: CommonPar(*), eps
    integer :: iShiftLag, iShiftDeltaU, NdimE, i, j , ipLag, ipDeltaU, NodG
	
	iShiftDeltaU = nint(commonPar(1))
	iShiftLag = nint(commonPar(2))

	NodG = (MaxLRows/iDofT) - 1

	NdimE = 2	
	do i = 1 , NodG
		ipLag = (i-1)*iDofT + iShiftLag
		ipDeltaU = (i-1)*iDofT + iShiftDeltaU
		do j = 1 , NdimE
			Coupling(ipLag + j , ipDeltaU + j) = iAdd
			Coupling(ipDeltaU + j , ipLag + j) = iAdd
!~ 			Coupling(ipDeltaU + j , ipDeltaU + j) = iAdd
			Coupling(ipLag + j , ipLag + j) = iAdd
		end do
	end do
	
	
end Subroutine

