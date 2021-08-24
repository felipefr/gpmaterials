	!> Element that solves the trivial diagonal problem of incrementing variables
	!! @param CommonPar = ( a , a*[bi,ci,di,ei] , dummy ) where : 
	!! @param Nvar = CommonPar(1) ==> a
	!! @param iShiftUi = CommonPar(see structure) ===> bi
	!! @param iShiftDeltaUi = CommonPar(see structure) ===> ci
	!! @param nodei (for each material = see structure) ===> ei 
	!! @author Rocha, Felipe Figueredo


!! last update: lenU given by parameter and not by difference between iShiftU and iShiftDeltaU. NodMax is given, rather than assumed
Subroutine zeroVariable_arcLength_generic(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, & 
											Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------

	use funcAux
	use TIME_STEP_ITER
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
    
    integer , parameter :: NauxParam = 4
    integer :: Nvar, iShiftU, lenU, node, NodG , node1, node2
    integer :: i, ip, jp, j, k, enableConditional

	Nvar = nint(commonPar(1))
	NodG = nint(commonPar(2))
	
	ip = NodG*iDofT + 1
	
	do i = 1 , Nvar
		ip = (i-1)*NauxParam + 2
		iShiftU = nint(commonPar(ip + 1))
		node = nint(commonPar(ip + 2))
		lenU = nint(commonPar(ip + 3))
		enableConditional = nint(commonPar(ip + 4))
				
		if(node == 0) then
			node1 = 1
			node2 = NodG
		else
			node1 = node
			node2 = node		
		end if
		
!~ 		write(0,*) iShiftU, node1, node2, lenU, enableConditional
		
		do j = node1,node2
			jp = (j - 1)*iDofT + iShiftU
			
			do k = 1,lenU
				AE(jp + k, jp + k) = 1.0d0
				
				if(enableConditional == 1) then
				
					if(Nconverg == 1) then
						BE(jp + k) = 0.0d0 
					else
						BE(jp + k) = Sol1(jp + k)
					end if
				
				else
					BE(jp + k) = 0.0d0
				end if
				
			end do	
		end do
		
	end do
	
!~ 	pause
	
End subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine zeroVariable_arcLength_genericS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)

    integer , parameter :: NauxParam = 4	
    integer :: Nvar, iShiftU, lenU, node, NodG , node1, node2
    integer :: i, ip, jp, j, k
    
	Nvar = nint(commonPar(1))
	NodG = nint(commonPar(2))
	
	do i = 1 , Nvar
		ip = (i-1)*NauxParam + 2
		iShiftU = nint(commonPar(ip + 1))
		node = nint(commonPar(ip + 2))
		lenU = nint(commonPar(ip + 3))
				
		if(node == 0) then
			node1 = 1
			node2 = NodG
		else
			node1 = node
			node2 = node		
		end if
				
		do j = node1,node2
			jp = (j - 1)*iDofT + iShiftU
		
			do k = 1,lenU
				Coupling(jp + k, jp + k) = 1
			end do	
		end do
		
	end do

!~ 	call numprint(Coupling)
!~ 	pause
end Subroutine

