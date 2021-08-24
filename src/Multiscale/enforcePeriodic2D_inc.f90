    !> Implements periodic boundary condition for 2D continuum problems 
    !!
	!! @param iShiftFluc = nint(CommonPar(1))
	!! @param iShiftLag = nint(CommonPar(2))
	!! @param pen  = CommonPar(3)
	!! @param eps  = commonPar(4)
	!!
    !! @author Rocha, Felipe Figueredo
    
Subroutine enforcePeriodic2D_inc &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)

	use funcAux
	
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
    integer, parameter :: NdimE = 2
    integer :: i , j , ip, iShiftLag, iShiftFluc , iShiftDeltaLag, iShiftDeltaFluc,  flagIncrementInside 
    Real*8 ::  pen,lamb, eps

!~ 	write(*,*) "enforcePeriodic2D" 
	iShiftFluc = nint(commonPar(1))
	iShiftDeltaFluc = nint(commonPar(2))
	iShiftLag = nint(commonPar(3))
	iShiftDeltaLag = nint(commonPar(4))
	pen = commonPar(5)
	eps = commonPar(6)
	flagIncrementInside  = nint(commonPar(7)) ! deprecated
	
	! Assemble stiffness matrix
	do i = 1 , NdimE
		AE(iShiftDeltaLag + i , iShiftDeltaFluc + i) = pen
		AE(iShiftDeltaLag + i , iShiftDeltaFluc + iDofT + i ) = -pen
		
		! transposes block
		AE(iShiftDeltaFluc + i , iShiftDeltaLag + i) = AE(iShiftDeltaLag + i , iShiftDeltaFluc + i)
		AE(iShiftDeltaFluc + iDofT + i, iShiftDeltaLag + i) = AE(iShiftDeltaLag + i , iShiftDeltaFluc + iDofT + i ) 

		AE(iShiftDeltaLag + i , iShiftDeltaLag + i) = eps
		
	end do

	if( iShiftLag /= iShiftDeltaLag .and. flagIncrementInside > -1) then
		do i = 1 , NdimE
			AE(iShiftLag + i , iShiftLag + i) = pen
			AE(iShiftLag + i , iShiftDeltaLag + i ) = -pen		
		end do
	end if


	! Assemble the force vector
	! no increments 
	if(iShiftFluc == iShiftDeltaFluc .and. iShiftLag == iShiftDeltaLag) then	
		do i = 1 , NdimE
			BE(iShiftDeltaLag + i) = eps*Sol1(iShiftLag + i)
		end do
	
	! incremental in displacements
	else if(iShiftFluc /= iShiftDeltaFluc .and. iShiftLag == iShiftDeltaLag) then
		do i = 1 , NdimE
			BE(iShiftDeltaLag + i) = eps*Sol1(iShiftLag + i)  - pen*( Sol1(iShiftFluc + i) - Sol1(iShiftFluc + iDofT + i))
		end do
		
	! incremental in multipliers
	else if(iShiftFluc == iShiftDeltaFluc .and. iShiftLag /= iShiftDeltaLag) then
		do i = 1 , NdimE
			BE(iShiftDeltaFluc + i) = - pen*Sol1(iShiftLag + i)
			BE(iShiftDeltaFluc + iDofT + i) =  pen*Sol1(iShiftLag + i)
			
!~ 			!just to increment variables
			if(flagIncrementInside > -1) then
				BE(iShiftLag + i) =  pen*Sol1(iShiftLag + i)
			end if
		end do
	! incremental
	else
		do i = 1 , NdimE
			BE(iShiftDeltaFluc + i) = - pen*Sol1(iShiftLag + i)
			BE(iShiftDeltaFluc + iDofT + i) =  pen*Sol1(iShiftLag + i)
			BE(iShiftDeltaLag + i) = - pen*( Sol1(iShiftFluc + i) - Sol1(iShiftFluc + iDofT + i))
			
			!just to increment variables
			if(flagIncrementInside > -1) then
				BE(iShiftLag + i) =  pen*Sol1(iShiftLag + i)
			end if
		end do	
	end if

end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine enforcePeriodic2D_incS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
	
    implicit none
    
    Integer :: MaxLRows,Ndim,iDofT, iAdd, iShiftLag, iShiftFluc , iShiftDeltaLag, iShiftDeltaFluc, NdimE, i, flagIncrementInside
    Integer ::  Coupling(MaxLRows,MaxLRows)
    Real*8 :: CommonPar(*)
	
	iShiftFluc = nint(commonPar(1))
	iShiftDeltaFluc = nint(commonPar(2))
	iShiftLag = nint(commonPar(3))
	iShiftDeltaLag = nint(commonPar(4))
!~ 	eps = commonPar(6)
	flagIncrementInside = nint(commonPar(7))

	NdimE = 2	
		
 	do i = 1 , NdimE
		Coupling(iShiftDeltaLag + i , iShiftDeltaLag + i) = iAdd
		Coupling(iShiftDeltaLag + i , iShiftDeltaFluc + i) = iAdd
		Coupling(iShiftDeltaFluc + i , iShiftDeltaLag + i) = iAdd
		Coupling(iShiftDeltaLag + i , iShiftDeltaFluc + iDofT + i ) = iAdd
		Coupling(iShiftDeltaFluc + iDofT + i, iShiftDeltaLag + i) = iAdd
	end do
	
	if( iShiftLag /= iShiftDeltaLag .and. flagIncrementInside > -1) then
		do i = 1 , NdimE
			write(0,*) iShiftLag + i, iShiftLag + i
			write(0,*) iShiftLag + i, iShiftDeltaLag + i
			Coupling(iShiftLag + i , iShiftLag + i) = iAdd
			Coupling(iShiftLag + i , iShiftDeltaLag + i ) = iAdd		
		end do
	end if
	
	
	
!~ 	write(*,*) eps, iAdd
!~ 	write(*,*) Coupling(iShiftLag + 1 , iShiftLag + 1), Coupling(iShiftLag + 2 , iShiftLag + 2)
!~ 	call numprint(coupling)
!~ 	PAUSE
	
end Subroutine
