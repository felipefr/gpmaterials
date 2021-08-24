	!> Element that solves the trivial diagonal problem of incrementing variables
	!! @param CommonPar = ( a , a*[bi,ci,di,ei] , dummy ) where : 
	!! @param Nvar = CommonPar(1) ==> a
	!! @param iShiftUi = CommonPar(see structure) ===> bi
	!! @param iShiftDeltaUi = CommonPar(see structure) ===> ci
	!! @param nodei (for each material = see structure) ===> ei 
	!! @author Rocha, Felipe Figueredo


!! last update: lenU given by parameter and not by difference between iShiftU and iShiftDeltaU. NodMax is given, rather than assumed
Subroutine zeroVariable_arcLength2(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
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
	integer :: iShift_dUf_Star, iShift_dUf_Bar, iShift_dLp_Star, iShift_dLp_Bar, &
			   iShift_DUf_acc, iShift_dl, iShift_Dl_acc, iShift_c1, iShift_c2, iShift_c3, &
			   iShift_AngleEst1, iShift_AngleEst2, iFEMtype, NodG
	
	iShift_dUf_Star = nint(CommonPar(1))
	iShift_dUf_Bar = nint(CommonPar(2))
	iShift_dLp_Star = nint(CommonPar(3))
	iShift_dLp_Bar = nint(CommonPar(4))
	iShift_DUf_acc = nint(CommonPar(5))
	iShift_dl = nint(CommonPar(6))
	iShift_Dl_acc = nint(CommonPar(7))
	iShift_c1 = nint(CommonPar(8))
	iShift_c2 = nint(CommonPar(9))
	iShift_c3 = nint(CommonPar(10))
	iShift_AngleEst1 = nint(CommonPar(11))
	iShift_AngleEst2 = nint(CommonPar(12))
	iFEMtype =  nint(commonPar(13)) 

	call setNodG(iFemtype, NodG)
	
	do i = 1, NodG
		do j = 1, NdimE
			ip = (i-1)*iDofT + j
			AE(iShift_dUf_Star + ip, iShift_dUf_Star + ip) = 1.0d0
			AE(iShift_dUf_Bar + ip, iShift_dUf_Bar + ip) = 1.0d0
			AE(iShift_dLp_Star + ip, iShift_dLp_Star + ip) = 1.0d0
			AE(iShift_dLp_Bar + ip, iShift_dLp_Bar + ip) = 1.0d0
			AE(iShift_DUf_acc + ip, iShift_DUf_acc + ip) = 1.0d0
		end do
	end do
	
	ip = NodG*iDofT + 1
	AE(iShift_dl + ip, iShift_dl + ip) = 1.0d0
	AE(iShift_Dl_acc + ip, iShift_Dl_acc + ip) = 1.0d0
	AE(iShift_c1 + ip, iShift_c1 + ip) = 1.0d0
	AE(iShift_c2 + ip, iShift_c2 + ip) = 1.0d0
	AE(iShift_c3 + ip, iShift_c3 + ip) = 1.0d0
	AE(iShift_AngleEst1 + ip, iShift_AngleEst1 + ip) = 1.0d0
	AE(iShift_AngleEst2 + ip, iShift_AngleEst2 + ip) = 1.0d0

	if(Nconverg == 1) then
		do i = 1, NodG
			do j = 1, NdimE
				ip = (i-1)*iDofT + j
				BE(iShift_DUf_acc + ip) = 0.0d0
			end do
		end do
		
		ip = NodG*iDofT + 1
		BE(iShift_Dl_acc + ip) = 0.0d0

	else 
		do i = 1, NodG
			do j = 1, NdimE
				ip = (i-1)*iDofT + j
				BE(iShift_DUf_acc + ip) = Sol1(iShift_DUf_acc + ip)
			end do
		end do
		
		ip = NodG*iDofT + 1
		BE(iShift_Dl_acc + ip) = Sol1(iShift_Dl_acc + ip)

	end if
	
	do i = 1, NodG
		do j = 1, NdimE
			ip = (i-1)*iDofT + j
			BE(iShift_dUf_Star + ip) = 0.0d0
			BE(iShift_dUf_Bar + ip) = 0.0d0
			BE(iShift_dLp_Star + ip) = 0.0d0
			BE(iShift_dLp_Bar + ip) = 0.0d0
		end do
	end do
	
	ip = NodG*iDofT + 1
	BE(iShift_dl + ip) = 0.0d0
	BE(iShift_c1 + ip) = 0.0d0
	BE(iShift_c2 + ip) = 0.0d0
	BE(iShift_c3 + ip) = 0.0d0
	BE(iShift_AngleEst1 + ip) = 0.0d0
	BE(iShift_AngleEst2 + ip) = 0.0d0
	
end subroutine 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine zeroVariable_arcLength2S(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	use globalVariables, only : NdimE 
    use ptsGaussLib , only : setNodG
	
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
    
    integer :: i, j , ip
	integer :: iShift_dUf_Star, iShift_dUf_Bar, iShift_dLp_Star, iShift_dLp_Bar, &
			   iShift_DUf_acc, iShift_dl, iShift_Dl_acc, iShift_c1, iShift_c2, iShift_c3, &
			   iShift_AngleEst1, iShift_AngleEst2, iFEMtype, NodG
    
	
	iShift_dUf_Star = nint(CommonPar(1))
	iShift_dUf_Bar = nint(CommonPar(2))
	iShift_dLp_Star = nint(CommonPar(3))
	iShift_dLp_Bar = nint(CommonPar(4))
	iShift_DUf_acc = nint(CommonPar(5))
	iShift_dl = nint(CommonPar(6))
	iShift_Dl_acc = nint(CommonPar(7))
	iShift_c1 = nint(CommonPar(8))
	iShift_c2 = nint(CommonPar(9))
	iShift_c3 = nint(CommonPar(10))
	iShift_AngleEst1 = nint(CommonPar(11))
	iShift_AngleEst2 = nint(CommonPar(12))
	iFEMtype =  nint(commonPar(13)) 
	
	call setNodG(iFemtype, NodG)
		
	do i = 1, NodG
		do j = 1, NdimE
			ip = (i-1)*iDofT + j
			Coupling(iShift_dUf_Star + ip, iShift_dUf_Star + ip) = 1
			Coupling(iShift_dUf_Bar + ip, iShift_dUf_Bar + ip) = 1
			Coupling(iShift_dLp_Star + ip, iShift_dLp_Star + ip) = 1
			Coupling(iShift_dLp_Bar + ip, iShift_dLp_Bar + ip) = 1
			Coupling(iShift_DUf_acc + ip, iShift_DUf_acc + ip) = 1
		end do
	end do
	
	ip = NodG*iDofT + 1
	Coupling(iShift_dl + ip, iShift_dl + ip) = 1
	Coupling(iShift_Dl_acc + ip, iShift_Dl_acc + ip) = 1
	Coupling(iShift_c1 + ip, iShift_c1 + ip) = 1
	Coupling(iShift_c2 + ip, iShift_c2 + ip) = 1
	Coupling(iShift_c3 + ip, iShift_c3 + ip) = 1
	Coupling(iShift_AngleEst1 + ip, iShift_AngleEst1 + ip) = 1
	Coupling(iShift_AngleEst2 + ip, iShift_AngleEst2 + ip) = 1

	Coupling(NodG*iDofT + iShift_Dl_acc + 1 , NodG*iDofT + iShift_Dl_acc + 1) = 1
	
end Subroutine

