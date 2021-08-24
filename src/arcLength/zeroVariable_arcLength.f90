	!> Element that solves the trivial diagonal problem of incrementing variables
	!! @param CommonPar = ( a , a*[bi,ci,di,ei] , dummy ) where : 
	!! @param Nvar = CommonPar(1) ==> a
	!! @param iShiftUi = CommonPar(see structure) ===> bi
	!! @param iShiftDeltaUi = CommonPar(see structure) ===> ci
	!! @param nodei (for each material = see structure) ===> ei 
	!! @author Rocha, Felipe Figueredo


!! last update: lenU given by parameter and not by difference between iShiftU and iShiftDeltaU. NodMax is given, rather than assumed
Subroutine zeroVariable_arcLength(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
									Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------

	use TIME_STEP_ITER
	use globalVariables, only : NdimE 
    use ptsGaussLib , only : setNodG
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
	integer :: i, j , ip
    integer :: iShift_DUf_acc, iShift_Dl_acc, NodG, iFEMtype

	iShift_DUf_acc = nint(CommonPar(1))
	iShift_Dl_acc = nint(CommonPar(2))
	iFEMtype =  nint(commonPar(3)) 

	call setNodG(iFemtype, NodG)
	
	do i = 1, NodG
		do j = 1, NdimE
			ip = (i-1)*iDofT + j
			AE(iShift_DUf_acc + ip, iShift_DUf_acc + ip) = 1.0d0
		end do
	end do
	
	AE(NodG*iDofT + iShift_Dl_acc + 1 , NodG*iDofT + iShift_Dl_acc + 1) = 1.0d0

	if(Nconverg == 1) then
		do i = 1, NodG
			do j = 1, NdimE
				ip = (i-1)*iDofT + j
				BE(iShift_DUf_acc + ip) = 0.0d0
			end do
		end do
		
		BE(NodG*iDofT + iShift_Dl_acc + 1) = 0.0d0

	else 
		do i = 1, NodG
			do j = 1, NdimE
				ip = (i-1)*iDofT + j
				BE(iShift_DUf_acc + ip) = Sol1(iShift_DUf_acc + ip)
			end do
		end do
		
		BE(NodG*iDofT + iShift_Dl_acc + 1) = Sol1(NodG*iDofT + iShift_Dl_acc + 1)

	end if
	
end subroutine 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine zeroVariable_arcLengthS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	use globalVariables, only : NdimE 
    use ptsGaussLib , only : setNodG
	
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
    
    integer :: i, j , ip
    integer :: iShift_DUf_acc, iShift_Dl_acc, NodG, 	iFEMtype 
    
	iShift_DUf_acc = nint(CommonPar(1))
	iShift_Dl_acc = nint(CommonPar(2))
	iFEMtype =  nint(commonPar(3)) 
	
	call setNodG(iFemtype, NodG)
		
	do i = 1, NodG
		do j = 1, NdimE
			ip = (i-1)*iDofT + j
			Coupling(iShift_DUf_acc + ip, iShift_DUf_acc + ip) = 1
		end do
	end do

	Coupling(NodG*iDofT + iShift_Dl_acc + 1 , NodG*iDofT + iShift_Dl_acc + 1) = 1
	
end Subroutine

