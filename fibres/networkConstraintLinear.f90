    !> Incorporates (minRestriction + zero average) for fibres network 
    !! 
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iShiftLag = nint(CommonPar(2))
	!! @param pen  = commonPar(3)
	!! @param eps  = commonPar(4)

    !! @author Rocha, Felipe Figueredo

!     ------------------------------------------------------------------
Subroutine networkConstraintLinear(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
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
     		
	Integer :: i, ipp1, ipp2, j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ApDim, AppDim! counters 
	integer, parameter ::  NodG = 2, NodLag1 = 1, NodLag2=2 , NdofLagMax = 2 ! just for 2D
	Real*8 :: MatPar(maxMatPar)
	real*8 , allocatable ::  Xel(:) 
	Real*8 :: eps, pen 
	Real*8 :: MatB1(NodG*NdimE,NdofLagMax), MatB2(NodG*NdimE,NdofLagMax), &
			  MatC1(NdofLagMax,NodG*NdimE), MatC2(NdofLagMax,NodG*NdimE), &
			  MatD1(NdofLagMax,NdofLagMax), MatD2(NdofLagMax,NdofLagMax), &
			  vecG1(NdofLagMax), vecG2(NdofLagMax)
	Real*8 :: nConnect_frac1, nConnect_frac2
	integer :: iShiftU, iShiftDeltaU, iShiftLag , iMaterial, constLaw
	integer :: flag1 , flag2, opConstraint, NdofLag
	
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iShiftLag = nint(CommonPar(3))
	pen = commonPar(4)
	eps = commonPar(5)
	opConstraint = nint(commonPar(6)) ! 0 -> Linear Constraint (the only one at moment)
		
	flag1 = nint(Param(Ipos_flag1))
	flag2 = nint(Param(Ipos_flag2))
	
	NdofLag = NdofLagMax ! only option at the moment
		
	nConnect_frac1 = 1.0d0/Param(Ipos_nConnectFibre1)
	nConnect_frac2 = 1.0d0/Param(Ipos_nConnectFibre2)
			
	MatC1 = 0.0d0
	MatC2 = 0.0d0
	
	
	if(flag1 > 0) then
		MatC1(1,1) = nConnect_frac1
		MatC1(2,2) = nConnect_frac1	
	end if

	if(flag2 > 0) then	
		MatC2(1,3) = nConnect_frac2
		MatC2(2,4) = nConnect_frac2
	end if

	
	MatB1 = transpose(MatC1)
	MatB2 = transpose(MatC2)
 
	MatD1 = 0.0d0
	MatD2 = 0.0d0
	
	ipp1 = (NodLag1-1)*iDofT + iShiftLag
	ipp2 = (NodLag2-1)*iDofT + iShiftLag
	do i = 1,NdofLag
		MatD1(i,i) = eps
		MatD2(i,i) = eps
		VecG1(i) = eps*Sol1(ipp1 + i)
		VecG2(i) = eps*Sol1(ipp2 + i)
	end do
		
	MatC1 = pen*MatC1
	MatC2 = pen*MatC2
	
	call setAEandBE_LM(AE,BE,matB1,matC1,matD1,vecG1,NodG,NodLag1,idofT,NdofLag,iShiftDeltaU,iShiftLag,NdimE)
	call setAEandBE_LM(AE,BE,matB2,matC2,matD2,vecG2,NodG,NodLag2,idofT,NdofLag,iShiftDeltaU,iShiftLag,NdimE)

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine networkConstraintLinearS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux , only : numprint, setCoupling_lagrangeMultipliers
	
	use globalVariables, only : NdimE
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
	integer ::  p , q, App, ApRow, BpCol, Bqp, A,B, iShiftU, iShiftLag , NdofLag, opConstraint
	integer , parameter :: NodG = 2, NodLag1 = 1, NodLag2 = 2
	
	iShiftU = nint(CommonPar(1))
	iShiftLag = nint(CommonPar(3))
	opConstraint = nint(commonPar(6)) ! 0 -> NBC afib , 1 -> NBC normal , 2 -> just Volume

	NdofLag = NdimE
 	
	call setCoupling_lagrangeMultipliers(Coupling,NodG,NodLag1,idofT,NdofLag,iShiftU,iShiftLag,NdimE)
	call setCoupling_lagrangeMultipliers(Coupling,NodG,NodLag2,idofT,NdofLag,iShiftU,iShiftLag,NdimE)
	
end subroutine













