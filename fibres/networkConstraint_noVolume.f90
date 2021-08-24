    !> Incorporates minRestriction for fibres network 
    !! 
	!! @param iShiftFluc = nint(CommonPar(1))
	!! @param iShiftLag = nint(CommonPar(2))
	!! @param pen  = commonPar(3)
	!! @param eps  = commonPar(4)

    !! @author Rocha, Felipe Figueredo

!     ------------------------------------------------------------------
Subroutine networkConstraint_noVolume(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
						Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux , only : numprint, setAEandBE_lagrangeMultipliers, getSliceAllocate
	use globalVariables, only : NdimE, getMaterial, maxMatPar, getAnisoTensorInv
	use fibresLib
	use fibresMod
		
	IMPLICIT NONE
	
	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	! Reals Vectors and matrices
	Real*8  :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), &
		Param(*), JParam(*), CommonPar(*)
	
	!   =====   END ARGUMENTS  =======
     		
	Integer :: i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ApDim, AppDim! counters 
	integer, parameter ::  NodG = 2, NodLag = 3 , NdofLag = 4 ! just for 2D
	Real*8 :: MatPar(maxMatPar)
	real*8 , allocatable ::  Xel(:) 
	Real*8 :: eps, pen 
	Real*8 :: afib(NdimE), Areafib, Lfib, Vf, signal, v(NdimE), MatB(NdofLag,NodG*NdimE), MatC(NdofLag,NdofLag)
	integer :: iShiftFluc, iShiftLag , iMaterial, constLaw
	integer :: flag1 , flag2, opConstraint
	
	
	iShiftFluc = nint(CommonPar(1))
	iShiftLag = nint(CommonPar(2))
	pen = commonPar(3)
	eps = commonPar(4)

	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	call get_afib(afib,Lfib,Xel)	

	Areafib = Param(Ipos_Areaf)
	afib = Param(Ipos_af:Ipos_af+1)
	Lfib = Param(Ipos_Lf)
	Vf = Param(Ipos_Vf)
	
	v = -matmul(Bten_invT,afib)
	
!~ 	v = afib
	flag1 = nint(Param(Ipos_flag1))
	flag2 = nint(Param(Ipos_flag2))
	
	MatB = 0.0d0

	if(flag1 > 0) then
		MatB(1,1) = -v(1)*Areafib
		MatB(2,1) = -v(2)*Areafib	
		MatB(3,2) = -v(1)*Areafib
		MatB(4,2) = -v(2)*Areafib
	end if

	if(flag2 > 0) then	
		MatB(1,3) = v(1)*Areafib
		MatB(2,3) = v(2)*Areafib
		MatB(3,4) = v(1)*Areafib
		MatB(4,4) = v(2)*Areafib
	end if
	
	MatC = 0.0d0
	do i = 1,NdofLag
		MatC(i,i) = eps
	end do
	
	call setAEandBE_lagrangeMultipliers(AE,BE,matB,matC,Sol1,eps,pen,NodG,NodLag,idofT,NdofLag,iShiftFluc,iShiftLag,NdimE)

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine networkConstraint_noVolumeS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux , only : numprint, setCoupling_lagrangeMultipliers
	
	use globalVariables, only : NdimE
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
	integer ::  p , q, App, ApRow, BpCol, Bqp, A,B, iShiftFluc, iShiftLag , NdofLag
	integer , parameter :: NodG = 2, NodLag = 3
	
	iShiftFluc = nint(CommonPar(1))
	iShiftLag = nint(CommonPar(2))

	NdofLag = NdimE*NdimE
 	
	call setCoupling_lagrangeMultipliers(Coupling,NodG,NodLag,idofT,NdofLag,iShiftFluc,iShiftLag,NdimE)
	
end subroutine





