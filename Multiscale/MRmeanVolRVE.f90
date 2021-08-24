    !> Incorporates (minRestriction + zero average) for fibres network 
    !! 
	!! @param iShiftFluc = nint(CommonPar(1))
	!! @param iShiftLag = nint(CommonPar(2))
	!! @param pen  = commonPar(3)
	!! @param eps  = commonPar(4)

    !! @author Rocha, Felipe Figueredo

!     ------------------------------------------------------------------
Subroutine MRmeanVolRVE(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
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
	integer, parameter ::  NodG = 2, NodLag = 3 , NdofLag = 6 ! just for 2D
	Real*8 :: MatPar(maxMatPar)
	real*8 , allocatable ::  Xel(:) 
	Real*8 :: eps, pen, lengthLine 
	Real*8 :: afib(NdimE), Areafib, Lfib, Vf, signal, v(NdimE), vecU(NdimE), MatB(NdofLag,NodG*NdimE), MatC(NdofLag,NdofLag)
	integer :: iShiftFluc, iShiftLag , iMaterial, constLaw
	integer :: flag1 , flag2, opConstraint
	
	
	iShiftFluc = nint(CommonPar(1))
	iShiftLag = nint(CommonPar(2))
	pen = commonPar(3)
	eps = commonPar(4)
	opConstraint = nint(commonPar(5)) ! 0 on boundary, 1 just volume
	
	MatB = 0.0d0
	MatC = 0.0d0

	if(opConstraint == 0) then

		do i = 1 , NdimE
			vecU(i) = XLL(Ndim + i) - XLL(i)
		end do
	 
		lengthLine = dsqrt(dot_product(vecU,vecU))
		v(1) = vecU(2)/lengthLine
		v(2) = -vecU(1)/lengthLine 

		MatB(1,1) = 0.5d0*lengthLine*v(1)
		MatB(2,1) = 0.5d0*lengthLine*v(2)	
		MatB(3,2) = 0.5d0*lengthLine*v(1)
		MatB(4,2) = 0.5d0*lengthLine*v(2)

		MatB(1,3) = 0.5d0*lengthLine*v(1)
		MatB(2,3) = 0.5d0*lengthLine*v(2)
		MatB(3,4) = 0.5d0*lengthLine*v(1)
		MatB(4,4) = 0.5d0*lengthLine*v(2)	
		
		do i = 1,4
			MatC(i,i) = eps
		end do

	else

		Vf = Param(Ipos_Vf)
		
		MatB(5,1) = 0.5d0*Vf
		MatB(5,3) = 0.5d0*Vf
		MatB(6,2) = 0.5d0*Vf
		MatB(6,4) = 0.5d0*Vf
		
		do i = 4,NdofLag
			MatC(i,i) = eps
		end do

	end if
	

	do i = 1,NdofLag
		MatC(i,i) = eps
	end do
	
	call setAEandBE_lagrangeMultipliers(AE,BE,matB,matC,Sol1,eps,pen,NodG,NodLag,idofT,NdofLag,iShiftFluc,iShiftLag,NdimE)

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine MRmeanVolRVES(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
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

	NdofLag = NdimE*NdimE + NdimE
 	
	call setCoupling_lagrangeMultipliers(Coupling,NodG,NodLag,idofT,NdofLag,iShiftFluc,iShiftLag,NdimE)
	
end subroutine





