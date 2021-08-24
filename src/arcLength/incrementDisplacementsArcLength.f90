	!> Element that solves the trivial diagonal problem of incrementing variables
	!! @param CommonPar = ( a , a*[bi,ci,di,ei] , dummy ) where : 
	!! @param Nvar = CommonPar(1) ==> a
	!! @param iShiftUi = CommonPar(see structure) ===> bi
	!! @param iShiftDeltaUi = CommonPar(see structure) ===> ci
	!! @param nodei (for each material = see structure) ===> ei 
	!! @author Rocha, Felipe Figueredo


!! last update: lenU given by parameter and not by difference between iShiftU and iShiftDeltaU. NodMax is given, rather than assumed
Subroutine incrementDisplacementsArcLength(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
									Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------

	use globalVariables, only : NdimE
	use ptsGaussLib, only : setNodG
	use funcAux
	use TIME_STEP_ITER
	use DETERMINANT
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
    integer :: i, j , ipDim, ipDof, NodG, iFEMtype
    integer	:: iShift_Uf, iShift_dUf_Star, iShift_dUf_Bar, iShift_dUf_acc, iShift_dLarcLength
	
    real*8 :: dl, vAux ,  alpha
    real*8 , allocatable, dimension(:) :: Sol_Uf , Sol_dUf , Sol_DUf_acc, Sol_dUf_Star , Sol_dUf_Bar
	
	iShift_dLarcLength = nint(CommonPar(1))
	iShift_Uf = nint(CommonPar(2))
	iShift_DUf_acc = nint(CommonPar(3))
	iShift_dUf_Star = nint(CommonPar(4))
	iShift_dUf_Bar = nint(CommonPar(5))	
	iFEMtype = nint(CommonPar(6))
	
	call setNodG(iFEMtype, NodG)
		
	dl = Sol1(NodG*iDofT + iShift_dLarcLength + 1)

	call getSliceAllocate(Sol_Uf,Sol1, 1, NodG, iShift_Uf + 1 , iShift_Uf + NdimE, iDofT)
	call getSliceAllocate(Sol_DUf_acc,Sol1, 1, NodG, iShift_DUf_acc + 1 , iShift_dUf_acc + NdimE, iDofT)
	call getSliceAllocate(Sol_dUf_Star,Sol1, 1, NodG, iShift_dUf_Star + 1 , iShift_dUf_Star + NdimE, iDofT)
	call getSliceAllocate(Sol_dUf_Bar,Sol1, 1, NodG, iShift_dUf_Bar + 1 , iShift_dUf_Bar + NdimE, iDofT)
	
	allocate(Sol_dUf(NodG*NdimE))
	
	Sol_dUf = Sol_dUf_Star + dl*Sol_dUf_Bar
	
	do i = 1, NodG
		do j = 1, NdimE
			ipDof = (i-1)*iDofT + j
			ipDim = (i-1)*NdimE + j
			AE(ipDof + iShift_Uf, ipDof + iShift_Uf ) = 1.0d0
			AE(ipDof + iShift_dUf_acc , ipDof + iShift_dUf_acc ) = 1.0d0
			
			BE(ipDof + iShift_Uf ) = Sol_Uf(ipDim) +  Sol_dUf(ipDim)
			BE(ipDof + iShift_dUf_acc ) = Sol_dUf_acc(ipDim ) +  Sol_dUf(ipDim)
			
		end do
	end do
	
!~ 	if(iElem == 619 ) then
!~ 		call numprint(BE(1:NodELT*iDofT))
!~ 		pause
!~ 	end if

	
end subroutine 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine incrementDisplacementsArcLengthS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux
	use globalVariables, only : NdimE
	use ptsGaussLib, only : setNodG
	
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
    
    integer	:: i, j, ipDof, iShift_Uf, iShift_DUf_acc, iFEMtype, NodG
    
    iShift_Uf = nint(CommonPar(2))
	iShift_DUf_acc = nint(CommonPar(3))
	iFEMtype = nint(CommonPar(6))
	
	call setNodG(iFEMtype, NodG)
		
	do i = 1, NodG
		do j = 1, NdimE
			ipDof = (i-1)*iDofT + j
			Coupling(ipDof + iShift_Uf, ipDof + iShift_Uf ) = 1
			Coupling(ipDof + iShift_dUf_acc , ipDof + iShift_dUf_acc ) = 1
		end do
	end do

!~ 	call numprint(Coupling)
!~ 	pause
end Subroutine

