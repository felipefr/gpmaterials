    !> Implements minimum restriction BC for continuum
    !!
    !!
	!! @param iShiftFluc = nint(CommonPar(1))
	!! @param iShiftLag = nint(CommonPar(2))
	!! @param pen  = CommonPar(3)
	!! @param eps  = commonPar(4)
	!!
    !! @author Rocha, Felipe Figueredo
!     ------------------------------------------------------------------
Subroutine minRestrictionRVE(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
						Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
		
	IMPLICIT NONE
	
	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	! Reals Vectors and matrices
	Real*8  :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), &
		Param(*), JParam(*), CommonPar(*)
	
	!   =====   END ARGUMENTS  =======
     		
	Integer :: i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp! counters 
	Integer :: NdofLag
	integer, parameter :: Ndim2D = 2, NodG = 2, NodLag = 3
	Real*8 :: eps, pen, lengthLine 
	Real*8 :: normal0(Ndim2D) , vecU(Ndim2D)
	integer :: iShiftFluc, iShiftLag 
	
	iShiftFluc = nint(CommonPar(1))
	iShiftLag = nint(CommonPar(2))
	pen = commonPar(3)
	eps = commonPar(4)

	NdofLag = Ndim2D*Ndim2D
	
	do i = 1 , Ndim2D
		vecU(i) = XLL(Ndim + i) - XLL(i)
	end do
 
	lengthLine = dsqrt(dot_product(vecU,vecU))
	normal0(1) = vecU(2)/lengthLine
	normal0(2) = -vecU(1)/lengthLine 
	 	
	AE = 0.0d0
	BE = 0.0d0
			
	B = NodLag
	BpCol = (B-1)*iDofT + iShiftLag
	
	A = 1
	ApRow  = (A-1)*iDofT + iShiftFluc
	AE(BpCol + 1,ApRow + 1) = pen*0.5d0*lengthLine*normal0(1) ! l11
	AE(BpCol + 2,ApRow + 2) = pen*0.5d0*lengthLine*normal0(2) ! l12
	AE(BpCol + 3,ApRow + 1) = pen*0.5d0*lengthLine*normal0(1) ! l21
	AE(BpCol + 4,ApRow + 2) = pen*0.5d0*lengthLine*normal0(2) ! l22
	AE(BpCol + 1,ApRow + 1) = AE(ApRow + 1,BpCol + 1)  
	AE(BpCol + 2,ApRow + 2) = AE(ApRow + 2,BpCol + 2)
	AE(BpCol + 3,ApRow + 1) = AE(ApRow + 1,BpCol + 3)
	AE(BpCol + 4,ApRow + 2) = AE(ApRow + 2,BpCol + 4)

	A = 2
	ApRow  = (A-1)*iDofT + iShiftFluc
	AE(BpCol + 1,ApRow + 1) = pen*0.5d0*lengthLine*normal0(1) ! l11
	AE(BpCol + 2,ApRow + 2) = pen*0.5d0*lengthLine*normal0(2) ! l12
	AE(BpCol + 3,ApRow + 1) = pen*0.5d0*lengthLine*normal0(1) ! l21
	AE(BpCol + 4,ApRow + 2) = pen*0.5d0*lengthLine*normal0(2) ! l22
	AE(BpCol + 1,ApRow + 1) = AE(ApRow + 1,BpCol + 1)  
	AE(BpCol + 2,ApRow + 2) = AE(ApRow + 2,BpCol + 2)
	AE(BpCol + 3,ApRow + 1) = AE(ApRow + 1,BpCol + 3)
	AE(BpCol + 4,ApRow + 2) = AE(ApRow + 2,BpCol + 4)

	AE(BpCol + 1,BpCol + 1) = eps
	AE(BpCol + 2,BpCol + 2) = eps
	AE(BpCol + 3,BpCol + 3) = eps
	AE(BpCol + 4,BpCol + 4) = eps
	
	BE(BpCol + 1) = eps*Sol1(BpCol + 1)
	BE(BpCol + 2) = eps*Sol1(BpCol + 2)
	BE(BpCol + 3) = eps*Sol1(BpCol + 3)
	BE(BpCol + 4) = eps*Sol1(BpCol + 4)
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine minRestrictionRVES(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
	integer ::  p , q, App, ApRow, BpCol, Bqp, A,B, iShiftFluc, iShiftLag , NdofLag
	integer , parameter :: NodG = 2, NodLag = 3, Ndim2D = 2
	
	iShiftFluc = nint(CommonPar(1))
	iShiftLag = nint(CommonPar(2))
	
	NdofLag = Ndim2D*Ndim2D
	
	B = NodLag
	BpCol = (B-1)*iDofT + iShiftLag
	
	A = 1
	ApRow  = (A-1)*iDofT + iShiftFluc
	Coupling(BpCol + 1,ApRow + 1) = 1 ! l11
	Coupling(BpCol + 2,ApRow + 2) = 1 ! l12
	Coupling(BpCol + 3,ApRow + 1) = 1 ! l21
	Coupling(BpCol + 4,ApRow + 2) = 1 ! l22
	Coupling(BpCol + 1,ApRow + 1) = 1  
	Coupling(BpCol + 2,ApRow + 2) = 1
	Coupling(BpCol + 3,ApRow + 1) = 1
	Coupling(BpCol + 4,ApRow + 2) = 1

	A = 2
	ApRow  = (A-1)*iDofT + iShiftFluc
	Coupling(BpCol + 1,ApRow + 1) = 1 ! l11
	Coupling(BpCol + 2,ApRow + 2) = 1 ! l12
	Coupling(BpCol + 3,ApRow + 1) = 1 ! l21
	Coupling(BpCol + 4,ApRow + 2) = 1 ! l22
	Coupling(BpCol + 1,ApRow + 1) = 1  
	Coupling(BpCol + 2,ApRow + 2) = 1
	Coupling(BpCol + 3,ApRow + 1) = 1
	Coupling(BpCol + 4,ApRow + 2) = 1
	
	Coupling(BpCol + 1,BpCol + 1) = 1
	Coupling(BpCol + 2,BpCol + 2) = 1
	Coupling(BpCol + 3,BpCol + 3) = 1
	Coupling(BpCol + 4,BpCol + 4) = 1
	
end subroutine





