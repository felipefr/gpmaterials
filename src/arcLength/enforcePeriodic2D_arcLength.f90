    !> Implements periodic boundary condition for 2D continuum problems 
    !!
	!! @param iShiftFluc = nint(CommonPar(1))
	!! @param iShiftLag = nint(CommonPar(2))
	!! @param pen  = CommonPar(3)
	!! @param eps  = commonPar(4)
	!!
    !! @author Rocha, Felipe Figueredo
    
Subroutine enforcePeriodic2D_arcLength &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)

	use funcAux
	
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
    integer, parameter :: NdimE = 2
    integer :: i , j , ip, iShift_dLM, iShift_LM, & 
				iShift_Uf, iShift_dUf,  computeStar_or_Bar
    Real*8 ::  pen,lamb, eps

!~ 	write(*,*) "enforcePeriodic2D" 
	iShift_Uf = nint(commonPar(1))
	iShift_LM = nint(commonPar(2))
	iShift_dUf = nint(commonPar(3)) ! it can be star or bar, never total
	iShift_dLM = nint(commonPar(4)) ! it can be star or bar, never total
	computeStar_or_Bar  = nint(commonPar(5))
	pen = commonPar(6)
	eps = commonPar(7)
	
	! Assemble stiffness matrix
	do i = 1 , NdimE
		AE(iShift_dLM + i , iShift_dUf + i) = pen
		AE(iShift_dLM + i , iShift_dUf + iDofT + i ) = -pen
		
		! transposes block
		AE(iShift_dUf + i , iShift_dLM + i) = pen
		AE(iShift_dUf + iDofT + i, iShift_dLM + i) = -pen 

		AE(iShift_dLM + i , iShift_dLM + i) = eps
		
	end do
	
	! Just the case of increments in both variables
	do i = 1 , NdimE
	
		if(computeStar_or_Bar == 1) then !Star
			BE(iShift_dUf + i) = - pen*Sol1(iShift_LM + i)
			BE(iShift_dUf + iDofT + i) =  pen*Sol1(iShift_LM + i)
			BE(iShift_dLM + i) = - pen*( Sol1(iShift_Uf + i) - Sol1(iShift_Uf + iDofT + i))
		else if(computeStar_or_Bar /= 2) then ! computeStar_or_Bar can only assume 1 or 2 as values, if 2 do nothing 
			write(0,*) 'Option not founnd 	computeStar_or_Bar= ', 	computeStar_or_Bar 
			stop
		end if 
	end do	


end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine enforcePeriodic2D_arcLengthS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
	
    implicit none
    
    Integer :: MaxLRows,Ndim,iDofT, iAdd, iShift_dLM, iShift_dUf, NdimE, i
    Integer ::  Coupling(MaxLRows,MaxLRows)
    Real*8 :: CommonPar(*)
	
	iShift_dUf = nint(commonPar(3))
	iShift_dLM = nint(commonPar(4))

	NdimE = 2	
		
 	do i = 1 , NdimE

		Coupling(iShift_dLM + i , iShift_dUf + i) = 1
		Coupling(iShift_dLM + i , iShift_dUf + iDofT + i ) = 1
		
		Coupling(iShift_dUf + i , iShift_dLM + i) = 1
		Coupling(iShift_dUf + iDofT + i, iShift_dLM + i) = 1
		
		Coupling(iShift_dLM + i , iShift_dLM + i) = 1
	end do
	
	
end Subroutine
