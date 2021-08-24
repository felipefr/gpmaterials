    !> Incorporates (minRestriction + zero average) for fibres network 
    !! 
	!! @param iShiftDeltaU = nint(CommonPar(1))
	!! @param iShiftLag = nint(CommonPar(2))
	!! @param pen  = commonPar(3)
	!! @param eps  = commonPar(4)

    !! @author Rocha, Felipe Figueredo

!     ------------------------------------------------------------------
Subroutine networkConstraint_delta(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
						Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux , only : numprint, setAEandBE_LM, getSliceAllocate
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
	integer, parameter ::  NodG = 2, NodLag = 3 , NdofLag = 6 ! just for 2D
	Real*8 :: MatPar(maxMatPar)
	real*8 , allocatable ::  Xel(:), SolU(:)
	Real*8 :: eps, pen 
	Real*8 :: afib(NdimE), Areafib, Lfib, Vf, signal, v(NdimE) 
	Real*8 :: MatB(NodG*NdimE,NdofLag), MatC(NdofLag,NodG*NdimE), MatD(NdofLag,NdofLag), VecG(NdofLag)
	integer :: iShiftDeltaU, iShiftU, iShiftLag , iMaterial, constLaw
	integer :: flag1 , flag2, opConstraint
	
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iShiftLag = nint(CommonPar(3))
	pen = commonPar(4)
	eps = commonPar(5)
	opConstraint = nint(commonPar(6)) ! 0 on boundary and 1 is in all the nodes
	

	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	call get_afib(afib,Lfib,Xel)	
	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)

	Areafib = Param(Ipos_Areaf)
	afib = Param(Ipos_af:Ipos_af+1)
	Lfib = Param(Ipos_Lf)
	Vf = Param(Ipos_Vf)
	
!~ 	v = -matmul(Bten_invT,afib)/measFibres
	
	v = afib
	
	if(opConstraint == 0) then
		flag1 = nint(Param(Ipos_flag1))
		flag2 = nint(Param(Ipos_flag2))
	else
		flag1 = 1
		flag2 = 1
	end if
	
	MatC = 0.0d0

	!!! This uses C convention [Gamma_11,Gamma_12,Gamma_21,Gamma_22]
	if(flag1 > 0) then
		MatC(1,1) = -v(1)*Areafib
		MatC(2,1) = -v(2)*Areafib	
		MatC(3,2) = -v(1)*Areafib
		MatC(4,2) = -v(2)*Areafib
	end if

	if(flag2 > 0) then	
		MatC(1,3) = v(1)*Areafib
		MatC(2,3) = v(2)*Areafib
		MatC(3,4) = v(1)*Areafib
		MatC(4,4) = v(2)*Areafib
	end if

	MatC(5,1) = 0.5d0*Vf
	MatC(5,3) = 0.5d0*Vf
	MatC(6,2) = 0.5d0*Vf
	MatC(6,4) = 0.5d0*Vf
!~ 
!~ 	MatB(5,1) = 0.5d0
!~ 	MatB(5,3) = 0.5d0
!~ 	MatB(6,2) = 0.5d0
!~ 	MatB(6,4) = 0.5d0

	MatB = pen*transpose(MatC)

	MatD = 0.0d0
	do i = 1,NdofLag
		MatD(i,i) = eps
	end do
	
	App = (NodLag-1)*iDofT + iShiftLag
	
	do i = 1, NdofLag
		VecG(i) = eps*Sol1(App + i)
	end do
	
	if(iShiftU/=iShiftDeltaU) VecG = VecG - matmul(matC,SolU)
	
	
	call setAEandBE_LM(AE,BE,matB,matC,matD,vecG,NodG,NodLag,idofT,NdofLag,iShiftDeltaU,iShiftLag,NdimE)

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine networkConstraint_deltaS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux , only : numprint, setCoupling_lagrangeMultipliers
	
	use globalVariables, only : NdimE
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
	integer ::  p , q, App, ApRow, BpCol, Bqp, A,B, iShiftDeltaU, iShiftLag , NdofLag
	integer , parameter :: NodG = 2, NodLag = 3
	
	iShiftDeltaU = nint(CommonPar(2))
	iShiftLag = nint(CommonPar(3))

	NdofLag = NdimE*NdimE + NdimE
 	
	call setCoupling_lagrangeMultipliers(Coupling,NodG,NodLag,idofT,NdofLag,iShiftDeltaU,iShiftLag,NdimE)
	
end subroutine





