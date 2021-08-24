    !> Implements periodic boundary condition for 2D continuum problems 
    !!
	!! @param iShiftFluc = nint(CommonPar(1))
	!! @param iShiftLag = nint(CommonPar(2))
	!! @param pen  = CommonPar(3)
	!! @param eps  = commonPar(4)
	!!
    !! @author Rocha, Felipe Figueredo
    
Subroutine enforcePeriodic2D &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)

	use funcAux
	
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
    integer, parameter :: NdimE = 2
    integer :: i , j , ip, iShiftLag, iShiftFluc
    Real*8 ::  pen,lamb, eps

!~ 	write(*,*) "enforcePeriodic2D" 
	iShiftFluc = nint(commonPar(1))
	iShiftLag = nint(commonPar(2))
	pen = commonPar(3)
	eps = commonPar(4)

	do i = 1 , NdimE
		if(eps>0.0d0) BE(iShiftLag + i) = eps*Sol1(iShiftLag + i)
		AE(iShiftLag + i , iShiftLag + i) = eps
		AE(iShiftLag + i , iShiftFluc + i) = pen
		AE(iShiftFluc + i , iShiftLag + i) = pen
		AE(iShiftLag + i , iShiftFluc + iDofT + i ) = -pen
		AE(iShiftFluc + iDofT + i, iShiftLag + i) = -pen
	end do

end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine enforcePeriodic2DS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
	
    implicit none
    
    Integer :: MaxLRows,Ndim,iDofT, iAdd, iShiftLag, iShiftFluc, NdimE, i
    Integer ::  Coupling(MaxLRows,MaxLRows)
    Real*8 :: CommonPar(*), eps
	
	iShiftFluc = nint(commonPar(1))
	iShiftLag = nint(commonPar(2))

	eps = commonPar(4)

	NdimE = 2	
		
 	do i = 1 , NdimE
		if(eps>0.0d0) then 
			Coupling(iShiftLag + i , iShiftLag + i) = iAdd
		end if
		Coupling(iShiftLag + i , iShiftFluc + i) = iAdd
		Coupling(iShiftFluc + i , iShiftLag + i) = iAdd
		Coupling(iShiftLag + i , iShiftFluc + iDofT + i ) = iAdd
		Coupling(iShiftFluc + iDofT + i, iShiftLag + i) = iAdd
	end do
	
!~ 	write(*,*) eps, iAdd
!~ 	write(*,*) Coupling(iShiftLag + 1 , iShiftLag + 1), Coupling(iShiftLag + 2 , iShiftLag + 2)
!~ 	call numprint(coupling)
!~ 	PAUSE
	
end Subroutine
