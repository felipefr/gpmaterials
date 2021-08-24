    !> Incorporates (minRestriction + zero average) for fibres network 
    !! 
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iShiftLag = nint(CommonPar(2))
	!! @param pen  = commonPar(3)
	!! @param eps  = commonPar(4)

    !! @author Rocha, Felipe Figueredo

!     ------------------------------------------------------------------
Subroutine networkConstraintGen_pureInc3D(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
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
	integer, parameter ::  NodG = 2, NodLag = 3 , NdofLagMax = 12 ! just for 2D
	Real*8 :: MatPar(maxMatPar)
	real*8 , allocatable ::  Xel(:), SolUf(:), SolLM(:)
	Real*8 :: eps, pen 
	Real*8 :: afib(NdimE), Areafib, Lfib, Vf, signal, v1(NdimE), v2(NdimE), &
			  MatB(NodG*NdimE,NdofLagMax), MatC(NdofLagMax,NodG*NdimE), MatD(NdofLagMax,NdofLagMax), &
			  vecG(NdofLagMax), vecF(NodG*NdimE)
	Real*8 :: A1, A2, normal1(NdimE), normal2(NdimE), nConnect_frac1, nConnect_frac2
	integer :: iShiftUf, iShiftDeltaUf, iShiftLag , iMaterial, constLaw
	integer :: flag1 , flag2, opConstraint, NdofLag
	
	
	iShiftUf = nint(CommonPar(1))
	iShiftDeltaUf = nint(CommonPar(2))
	iShiftLag = nint(CommonPar(3))
	pen = commonPar(4)
	eps = commonPar(5)
	opConstraint = nint(commonPar(6)) ! 0 -> NBC afib , 1 -> NBC normal , 2 -> just Volume
	
	flag1 = nint(Param(Ipos_flag1))
	flag2 = nint(Param(Ipos_flag2))
	Vf = Param(Ipos_Vf)
	
	if(opConstraint<2) then
		NdofLag = NdimE*NdimE + NdimE
	else
		NdofLag = NdimE
	end if
	
	call getSliceAllocate(SolUf,Sol1,1,NodG,iShiftUf + 1 ,iShiftUf + NdimE, iDofT)
	call getSliceAllocate(SolLM,Sol1,NodLag,NodLag,iShiftLag + 1 ,iShiftLag + NdofLag, iDofT)

	if(opConstraint == 0) then ! NBC afib
	
		Areafib = Param(Ipos_Areaf)
		afib = Param(Ipos_af:Ipos_af + NdimE - 1)		
		
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
		
		v1 = Param(Ipos_normal1 : Ipos_normal1 + NdimE - 1)
		v2 = Param(Ipos_normal2 : Ipos_normal2 + NdimE - 1)
	
	end if
	
	
!~ 	write(*,*) A1, A2, v1, v2 
		
	MatC = 0.0d0
	
	do i = 1, NdimE
		MatC(i,i) = 0.5d0*Vf
		MatC(i,NdimE + i) = 0.5d0*Vf
	end do
	
	if(opConstraint < 2) then

		if(flag1 > 0) then
			do i = 1, NdimE
				ipp = NdimE + (i-1)*NdimE
				do j = 1, NdimE
					MatC(ipp + j,i) = nConnect_frac1*A1*v1(j)
				end do 
			end do
		end if

		if(flag2 > 0) then	
			do i = 1, NdimE
				ipp = NdimE + (i-1)*NdimE
				do j = 1, NdimE
					MatC(ipp + j,NdimE + i) = nConnect_frac2*A2*v2(j)
				end do 
			end do
		end if
	
	end if

	MatB = transpose(MatC)
 
	MatD = 0.0d0
	ipp = (NodLag-1)*iDofT + iShiftLag
	do i = 1,NdofLag
		MatD(i,i) = eps
	end do
	
	MatC = pen*MatC
	
	VecF = - matmul(MatB(:,1:NdofLag),SolLM(1:NdofLag))
	VecG(1:NdofLag) = - matmul(MatC(1:NdofLag,:),SolUf)

	! [..., B ; C, D] [du ; dlambda] = [vecF ; vecG]  
	call setAEandBE_LM_full(AE,BE,matB,matC,matD,vecF,vecG,NodG,NodLag,idofT,NdofLag,iShiftDeltaUf,iShiftLag,NdimE)

!~ 	write(0,*) 'Running network constraint' , NdimE 

!~ 	call numprint(AE)
!~ 	call numprint(BE)
	
!~ 	pause

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine networkConstraintGen_pureInc3DS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux , only : numprint, setCoupling_lagrangeMultipliers
	
	use globalVariables, only : NdimE
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
	integer ::  p , q, App, ApRow, BpCol, Bqp, A,B, iShiftDeltaUf, iShiftLag , NdofLag, opConstraint
	integer , parameter :: NodG = 2, NodLag = 3
	
	iShiftDeltaUf = nint(CommonPar(2))
	iShiftLag = nint(CommonPar(3))
	opConstraint = nint(commonPar(6)) ! 0 -> NBC afib , 1 -> NBC normal , 2 -> just Volume

	if(opConstraint<2) then
		NdofLag = NdimE*NdimE + NdimE
	else
		NdofLag = NdimE
	end if
 	
	call setCoupling_lagrangeMultipliers(Coupling,NodG,NodLag,idofT,NdofLag,iShiftDeltaUf,iShiftLag,NdimE)
	
end subroutine













