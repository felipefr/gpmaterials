Subroutine zeroVariable &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)

	use funcAux
	use loadingLib
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
   
    Real*8 :: pen
    integer :: Nnodes, NdofPerNode 
	integer , parameter :: iShift = 3
	integer :: i, j, ip, dof, node

	Nnodes = nint(commonPar(1))
	NdofPerNode = nint(commonPar(2))
	pen = commonPar(3)

	do i = 1, Nnodes
		node = nint(commonPar(iShift + i))
		do j = 1, NdofPerNode
			dof = nint(commonPar(iShift + Nnodes + j))
			
			ip = (node-1)*idofT + dof
						
			AE(ip,ip) = pen
			BE(ip) = 0.0d0
			
		end do
	end do 
	
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine zeroVariableS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    implicit none
    
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    integer :: Nnodes, NdofPerNode 
	integer , parameter :: iShift = 3
	integer :: i, j, ip, dof, node

	Nnodes = nint(commonPar(1))
	NdofPerNode = nint(commonPar(2))
	
	do i = 1, Nnodes
		node = nint(commonPar(iShift + i))
		do j = 1, NdofPerNode
			dof = nint(commonPar(iShift + Nnodes + j))
			
			ip = (node-1)*idofT + dof
						
			Coupling(ip,ip) = 1
			
		end do
	end do 
	
end Subroutine
