    !> Incorporates (minRestriction + zero average) for fibres network 
    !! 
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iShiftLag = nint(CommonPar(2))
	!! @param pen  = commonPar(3)
	!! @param eps  = commonPar(4)

    !! @author Rocha, Felipe Figueredo

!     ------------------------------------------------------------------
Subroutine networkConstraintLinear3D(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
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
     		
	Integer :: i, ip, ipp, j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ApDim, AppDim! counters 
	integer, parameter ::  NodG = 2 
	Real*8 :: MatPar(maxMatPar)
	real*8 , allocatable ::  Xel(:) 
	Real*8 :: eps, pen 
	Real*8 :: MatB(NodG*NdimE,NdimE), MatC(NdimE,NodG*NdimE), MatD(NdimE,NdimE), vecG(NdimE)
	Real*8 :: nConnect_frac
	integer :: iShiftU, iShiftDeltaU, iShiftLag , iMaterial, constLaw
	integer :: flag, opConstraint, NdofLag, NodLag	
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iShiftLag = nint(CommonPar(3))
	NodLag = nint(CommonPar(4))
	pen = commonPar(5)
	eps = commonPar(6)
	opConstraint = nint(commonPar(7)) ! 0 -> Linear Constraint (the only one at moment)
		
	if(NodLag == 1) then
		flag = nint(Param(Ipos_flag1))
		nConnect_frac = 1.0d0/Param(Ipos_nConnectFibre1)
	else if(NodLag == 2) then
		flag = nint(Param(Ipos_flag2))
		nConnect_frac = 1.0d0/Param(Ipos_nConnectFibre2)
	end if
	
	NdofLag = NdimE ! only option at the moment
	
	MatC = 0.0d0

	if(flag > 0) then
		ip = (NodLag-1)*NdimE
		do i = 1, NdimE
				MatC(i,ip + i) = nConnect_frac	
		end do	
	end if
		
	MatB = transpose(MatC)
 
	MatD = 0.0d0
	
	ipp = (NodLag-1)*iDofT + iShiftLag
	do i = 1,NdofLag
		MatD(i,i) = eps
		VecG(i) = eps*Sol1(ipp + i)
	end do
		
	MatC = pen*MatC
	
	call setAEandBE_LM(AE,BE,matB,matC,matD,vecG,NodG,NodLag,idofT,NdofLag,iShiftDeltaU,iShiftLag,NdimE)

end subroutine


!~ Subroutine networkConstraintLinear3D(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
!~ 						Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!~ !     ------------------------------------------------------------------
!~ 	use funcAux , only : numprint, setAEandBE_lagrangeMultipliers, getSliceAllocate
!~ 	use globalVariables, only : NdimE, getMaterial, maxMatPar, getAnisoTensorInv
!~ 	use fibresLib
!~ 	use fibresMod
		
!~ 	IMPLICIT NONE
	
!~ 	!   ===== SUBROUTINE ARGUMENTS  =======
!~ 	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
!~ 	Real*8 :: DelT, DTm,Time ! all reals
!~ 	! Reals Vectors and matrices
!~ 	Real*8  :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), &
!~ 		Param(*), JParam(*), CommonPar(*)
	
!~ 	!   =====   END ARGUMENTS  =======
     		
!~ 	Integer :: i, ipp1, ipp2, j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ApDim, AppDim! counters 
!~ 	integer, parameter ::  NodG = 2, NodLag1 = 1, NodLag2=2 , NdofLagMax = 3 
!~ 	Real*8 :: MatPar(maxMatPar)
!~ 	real*8 , allocatable ::  Xel(:) , SolU(:), SolLM1(:), SolLM2(:)
!~ 	Real*8 :: eps, pen 
!~ 	Real*8 :: MatB1(NodG*NdimE,NdofLagMax), MatB2(NodG*NdimE,NdofLagMax), &
!~ 			  MatC1(NdofLagMax,NodG*NdimE), MatC2(NdofLagMax,NodG*NdimE), &
!~ 			  MatD1(NdofLagMax,NdofLagMax), MatD2(NdofLagMax,NdofLagMax), &
!~ 			  vecG1(NdofLagMax), vecG2(NdofLagMax), vecF1(NdimE*NodG), vecF2(NdimE*NodG) 
!~ 	Real*8 :: nConnect_frac1, nConnect_frac2
!~ 	integer :: iShiftU, iShiftDeltaU, iShiftLag , iMaterial, constLaw
!~ 	integer :: flag1 , flag2, opConstraint, NdofLag
	
	
!~ 	iShiftU = nint(CommonPar(1))
!~ 	iShiftDeltaU = nint(CommonPar(2))
!~ 	iShiftLag = nint(CommonPar(3))
!~ 	pen = commonPar(4)
!~ 	eps = commonPar(5)
!~ 	opConstraint = nint(commonPar(6)) ! 0 -> Linear Constraint (no increment) , 1 -> Linear Constraint (pure increment)
		
!~ 	flag1 = nint(Param(Ipos_flag1))
!~ 	flag2 = nint(Param(Ipos_flag2))
	
!~ 	NdofLag = NdofLagMax ! only option at the moment
		
!~ 	nConnect_frac1 = 1.0d0/Param(Ipos_nConnectFibre1)
!~ 	nConnect_frac2 = 1.0d0/Param(Ipos_nConnectFibre2)
			
!~ 	MatC1 = 0.0d0
!~ 	MatC2 = 0.0d0
	
!~ 	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
!~ 	call getSliceAllocate(SolLM1,Sol1,NodLag1,NodLag1,iShiftLag + 1 ,iShiftLag + NdofLag, iDofT)
!~ 	call getSliceAllocate(SolLM2,Sol1,NodLag2,NodLag2,iShiftLag + 1 ,iShiftLag + NdofLag, iDofT)
	
	
!~ 	if(flag1 > 0) then
!~ 		do i = 1, NdimE
!~ 			MatC1(i,i) = nConnect_frac1
!~ 		end do 	
!~ 	end if

!~ 	if(flag2 > 0) then	
!~ 		do i = 1, NdimE
!~ 			MatC2(i,NdimE + i) = nConnect_frac2
!~ 		end do
!~ 	end if

	
!~ 	MatB1 = transpose(MatC1)
!~ 	MatB2 = transpose(MatC2)
 
!~ 	MatD1 = 0.0d0
!~ 	MatD2 = 0.0d0
	
!~ 	ipp1 = (NodLag1-1)*iDofT + iShiftLag
!~ 	ipp2 = (NodLag2-1)*iDofT + iShiftLag
!~ 	do i = 1,NdofLag
!~ 		MatD1(i,i) = eps
!~ 		MatD2(i,i) = eps
!~ 	end do
		
!~ 	MatC1 = pen*MatC1
!~ 	MatC2 = pen*MatC2
	

	
!~ 	if(opConstraint == 0) then
!~ 		do i = 1,NdofLag	
!~ 			VecG1(i) = eps*Sol1(ipp1 + i)
!~ 			VecG2(i) = eps*Sol1(ipp2 + i)
!~ 		end do
!~ 	else if(opConstraint == 1) then
!~ 		VecF1 = - matmul(MatB1,SolLM1)
!~ 		VecG1 = - matmul(MatC1,SolU)
!~ 		VecF2 = - matmul(MatB2,SolLM2)
!~ 		VecG2 = - matmul(MatC2,SolU)
!~ 	end if

!	write(0,*) size(SolLM1,1), size(SolLM2,1),size(SolU,1)
!~ !~! 	write(0,*) size(MatB1,1), size(MatB1,2)
!	write(0,*) size(MatB2,1), size(MatB2,2)
!	write(0,*) size(MatC1,1), size(MatC1,2)
!~ !~! 	write(0,*) size(MatC2,1), size(MatC2,2)

!~ 		! [..., B ; C, D] [du ; dlambda] = [vecF ; vecG]  	
!~ 	call setAEandBE_LM_full(AE,BE,matB1,matC1,matD1,vecF1,vecG1,NodG,NodLag1,idofT,NdofLag,iShiftDeltaU,iShiftLag,NdimE)
!~ 	call setAEandBE_LM_full(AE,BE,matB2,matC2,matD2,vecF2,vecG2,NodG,NodLag2,idofT,NdofLag,iShiftDeltaU,iShiftLag,NdimE)

!~ end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine networkConstraintLinear3DS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux , only : numprint, setCoupling_lagrangeMultipliers
	
	use globalVariables, only : NdimE
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
	integer ::  iShiftU, iShiftLag , NdofLag, NodLag
	integer , parameter :: NodG = 2
	
	iShiftU = nint(CommonPar(1))
	iShiftLag = nint(CommonPar(3))
	NodLag = nint(CommonPar(4))

	NdofLag = NdimE
 	
	call setCoupling_lagrangeMultipliers(Coupling,NodG,NodLag,idofT,NdofLag,iShiftU,iShiftLag,NdimE)
	
end subroutine













