    !> Incorporates (minRestriction + zero average) for fibres network 
    !! 
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iShiftLag = nint(CommonPar(2))
	!! @param pen  = commonPar(3)
	!! @param eps  = commonPar(4)

    !! @author Rocha, Felipe Figueredo

!     ------------------------------------------------------------------
Subroutine networkConstraint_RVEnormal_new_gen(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
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
     		
	Integer :: i, ipp, j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ApDim, AppDim! counters 
	integer, parameter ::  NodG = 2, NodLag = 3 , NdofLag = 6 ! just for 2D
	Real*8 :: MatPar(maxMatPar) 
	real*8 , allocatable ::  SolU(:) , SolL(:)
	Real*8 :: eps, pen 
	Real*8 :: afib(NdimE), Areafib, Lfib, Vf, signal, v1(NdimE), v2(NdimE), &
			  MatB(NodG*NdimE,NdofLag), MatC(NdofLag,NodG*NdimE), MatD(NdofLag,NdofLag), &
			  vecF(NodG*NdimE) , vecG(NdofLag)
	Real*8 :: A1, A2, normal1(NdimE), normal2(NdimE), nConnect_frac1, nConnect_frac2
	integer :: iShiftU, iShiftDeltaU, iShiftLag , iShiftDeltaLag, iMaterial, constLaw, flagIncrementInside
	integer :: flag1 , flag2, opConstraint
	
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iShiftLag = nint(CommonPar(3))
	iShiftDeltaLag = nint(CommonPar(4))
	pen = commonPar(5)
	eps = commonPar(6)
	opConstraint = nint(commonPar(7)) ! 0 -> NBC afib , 1 -> NBC normal
	flagIncrementInside = nint(commonPar(8))
	
	flag1 = nint(Param(Ipos_flag1))
	flag2 = nint(Param(Ipos_flag2))
	Vf = Param(Ipos_Vf)
	
	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(SolL,Sol1,NodLag,NodLag,iShiftLag + 1 ,iShiftLag + NdofLag, iDofT)


	! ============ Prepare to Assemble stiffness matrix , creates C , B and D (no increment depedent!!) ================

	if(opConstraint == 0) then ! NBC afib
	
		Areafib = Param(Ipos_Areaf)
		afib = Param(Ipos_af:Ipos_af+1)		
		
		nConnect_frac1 = 1.0d0
		nConnect_frac2 = 1.0d0
		
		A1 = AreaFib
		A2 = AreaFib
		
!~ 		v2 = -matmul(Bten_invT,afib)/measFibres
		v2 = afib 
		v1 = -v2
		
	else if(opConstraint == 1) then ! NBC normal
		
		nConnect_frac1 = 1.0d0/Param(Ipos_nConnectFibre1)
		nConnect_frac2 = 1.0d0/Param(Ipos_nConnectFibre2)
		
		A1 = Param(Ipos_Abar1)
		A2 = Param(Ipos_Abar2)
		
		v1 = Param(Ipos_normal1 : Ipos_normal1 + 1)
		v2 = Param(Ipos_normal2 : Ipos_normal2 + 1)
	
	end if
	
!~ 	write(*,*) A1, A2, v1, v2 
	
	MatC = 0.0d0

	if(flag1 > 0) then
		MatC(1,1) = nConnect_frac1*A1*v1(1)
		MatC(2,1) = nConnect_frac1*A1*v1(2)	
		MatC(3,2) = nConnect_frac1*A1*v1(1)
		MatC(4,2) = nConnect_frac1*A1*v1(2)
	end if

	if(flag2 > 0) then	
		MatC(1,3) = nConnect_frac2*A2*v2(1)
		MatC(2,3) = nConnect_frac2*A2*v2(2)
		MatC(3,4) = nConnect_frac2*A2*v2(1)
		MatC(4,4) = nConnect_frac2*A2*v2(2)
	end if

	MatC(5,1) = 0.5d0*Vf
	MatC(5,3) = 0.5d0*Vf
	MatC(6,2) = 0.5d0*Vf
	MatC(6,4) = 0.5d0*Vf
	
	MatB = transpose(MatC)

	MatD = 0.0d0
	do i = 1,NdofLag
		MatD(i,i) = eps
	end do


! =================== end creation matrices ===============================

! ====== increment lagrange multiplier inside if requested, matrix part =============
	if(iShiftLag /= iShiftDeltaLag .and. flagIncrementInside > -1) then
		do i = 1 , NdimE
			AE(iShiftLag + i , iShiftLag + i) = pen
			AE(iShiftLag + i , iShiftDeltaLag + i ) = -pen		
		end do
	end if

!============= Creation of Vectors G and F (RHS) , increment dependant

	VecF = 0.0d0
	VecG = 0.0d0
	
	! Assemble the force vector
	! no increments 
	if(iShiftU == iShiftDeltaU .and. iShiftLag == iShiftDeltaLag) then	
		ipp = (NodLag-1)*iDofT + iShiftLag
		
		VecG = eps*SolL
	
	! incremental in displacements
	else if(iShiftU /= iShiftDeltaU .and. iShiftLag == iShiftDeltaLag) then

		VecG = eps*SolL - matmul(matC,SolU)
		
	! incremental in multipliers
	else if(iShiftU == iShiftDeltaU .and. iShiftLag /= iShiftDeltaLag) then
		
		VecF = -pen*matmul(matB,SolL)
		
		!just to increment variables
		if(flagIncrementInside > -1) then
			BE(iShiftLag + 1: iShiftLag + NdofLag) =  pen*SolL
		end if
	! incremental
	else
		VecF = -pen*matmul(matB,SolL)
		VecG = - matmul(matC,SolU)
		
		!just to increment variables
		if(flagIncrementInside > -1) then
			BE(iShiftLag + 1: iShiftLag + NdofLag) =  pen*SolL
		end if
	end if

	
	call setAEandBE_LM_full(AE,BE,matB,matC,matD,vecF,vecG, &
					NodG,NodLag,idofT,NdofLag,iShiftDeltaU,iShiftDeltaLag,NdimE)

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine networkConstraint_RVEnormal_new_genS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux , only : numprint, setCoupling_lagrangeMultipliers
	
	use globalVariables, only : NdimE
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
	integer ::  i, p , q, App, ApRow, BpCol, Bqp, A,B, iShiftDeltaU, iShiftLag,  iShiftDeltaLag , NdofLag, flagIncrementInside
	integer , parameter :: NodG = 2, NodLag = 3
	
	iShiftDeltaU = nint(CommonPar(2))
	iShiftLag = nint(CommonPar(3))
	iShiftDeltaLag = nint(CommonPar(4))
	flagIncrementInside = nint(commonPar(8))

	NdofLag = NdimE*NdimE + NdimE
 	
	call setCoupling_lagrangeMultipliers(Coupling,NodG,NodLag,idofT,NdofLag,iShiftDeltaU,iShiftDeltaLag,NdimE)
	
	
	
	if( iShiftLag /= iShiftDeltaLag .and. flagIncrementInside > -1) then
		do i = 1 , NdimE
			Coupling(iShiftLag + i , iShiftLag + i) = 1
			Coupling(iShiftLag + i , iShiftDeltaLag + i ) = 1		
		end do
	end if
	
end subroutine













