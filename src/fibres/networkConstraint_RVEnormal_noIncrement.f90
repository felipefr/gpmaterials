    !> Incorporates the minimal restriction model in terms of normal for fibres network 
    !! 
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iShiftDeltaU = nint(CommonPar(2))
	!! @param iShiftLag = nint(CommonPar(3))
	!! @param pen  = commonPar(4)
	!! @param eps  = commonPar(5)

    !! @author Rocha, Felipe Figueredo

!     ------------------------------------------------------------------
Subroutine networkConstraint_RVEnormal_noIncrement(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
						Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux , only : numprint, setAEandBE_LM, getSliceAllocate, getNormal
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
	real*8 , allocatable ::  Xel(:) 
	Real*8 :: eps, pen 
	Real*8 :: afib(NdimE), Areafib, Lfib, Vf, signal, v1(NdimE), v2(NdimE), & 
				MatB(NdofLag,NodG*NdimE), MatC(NdofLag,NdofLag), normalAux(2,NdimE)
	integer :: iShiftU, iShiftLag , iMaterial, constLaw
	integer :: flag1 , flag2, opConstraint
	
	
	iShiftU = nint(CommonPar(1))
	iShiftLag = nint(CommonPar(2))
	pen = commonPar(3)
	eps = commonPar(4)
	opConstraint = nint(commonPar(5)) ! 0 on boundary and 1 is in all the nodes

	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	call get_afib(afib,Lfib,Xel)	

	Areafib = Param(Ipos_Areaf)
	afib = Param(Ipos_af:Ipos_af+1)
!~ 	Lfib = Param(Ipos_Lf)
	Vf = Param(Ipos_Vf)
	
	flag1 = nint(Param(Ipos_flag1))
	flag2 = nint(Param(Ipos_flag2))
	
	
	! 0 =  interior, 1 = bottom, 2 = right, 3 = top, 4 = left
	! 5 = bottom/left , 6 = bottom/right, 7 = right/top, 8 = top/left  
	
	!! choose of v1
	if(flag1>0) then
		call getNormal(normalAux,flag1)
		if(flag1>4) then
			v1 = Areafib*dabs(dot_product(normalAux(:,1),afib))*normalAux(:,1) + &
			 Areafib*dabs(dot_product(normalAux(:,2),afib))*normalAux(:,2)
		else 
			v1 = Areafib*dabs(dot_product(normalAux(:,1),afib))*normalAux(:,1)
		end if
	end if

	if(flag2>0) then
		call getNormal(normalAux,flag2)
		if(flag2>4) then
			v2 = Areafib*dabs(dot_product(normalAux(:,1),afib))*normalAux(:,1) + &
			   Areafib*dabs(dot_product(normalAux(:,2),afib))*normalAux(:,2)
		else 
			v2 = Areafib*dabs(dot_product(normalAux(:,1),afib))*normalAux(:,1)
		end if
	end if	
	
	MatB = 0.0d0

	!!! This uses C convention [Gamma_11,Gamma_12,Gamma_21,Gamma_22]
	if(flag1 > 0) then
		MatB(1,1) = v1(1)
		MatB(2,1) = v1(2)	
		MatB(3,2) = v1(1)
		MatB(4,2) = v1(2)
	end if

	if(flag2 > 0) then	
		MatB(1,3) = v2(1)
		MatB(2,3) = v2(2)
		MatB(3,4) = v2(1)
		MatB(4,4) = v2(2)
	end if

	MatC(5,1) = 0.5d0*Vf
	MatC(5,3) = 0.5d0*Vf
	MatC(6,2) = 0.5d0*Vf
	MatC(6,4) = 0.5d0*Vf
!~ 
!~ 
!~ 	MatB(5,1) = 0.5d0
!~ 	MatB(5,3) = 0.5d0
!~ 	MatB(6,2) = 0.5d0
!~ 	MatB(6,4) = 0.5d0

	MatC = 0.0d0
	do i = 1,NdofLag
		MatC(i,i) = eps
	end do
	
	call setAEandBE_lagrangeMultipliers(AE,BE,matB,matC,Sol1,eps,pen,NodG,NodLag,idofT,NdofLag,iShiftU,iShiftLag,NdimE)

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine networkConstraint_RVEnormal_noIncrementS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux , only : numprint, setCoupling_lagrangeMultipliers
	
	use globalVariables, only : NdimE
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
	integer ::  p , q, App, ApRow, BpCol, Bqp, A,B, iShiftU, iShiftLag , NdofLag
	integer , parameter :: NodG = 2, NodLag = 3
	
	iShiftU = nint(CommonPar(1))
	iShiftLag = nint(CommonPar(2))

	NdofLag = NdimE*NdimE + NdimE
 	
	call setCoupling_lagrangeMultipliers(Coupling,NodG,NodLag,idofT,NdofLag,iShiftU,iShiftLag,NdimE)
	
	
end subroutine





