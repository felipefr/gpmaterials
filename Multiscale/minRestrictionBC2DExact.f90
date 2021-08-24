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
Subroutine minRestrictionBC2DExact(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
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
	
	Do A=1,NodG 
		ApRow  = (A-1)*iDofT + iShiftFluc
		
		Do p = 1 , Ndim2D
			App = ApRow+p
												 
			Do q = 1 , NdofLag
				Bqp = BpCol + q
				
				call getIJfromK(i,j,q)
				
				AE(App,Bqp) = AE(App,Bqp) + pen*0.5*lengthLine*deltaKron(p,i)*normal0(j)
				AE(Bqp,App) = AE(App,Bqp)
				AE(Bqp,Bqp) = eps
!~ 				BE(Bqp) = eps*Sol1(Bqp)
			end do
		end do ! loop Ap
	Enddo !LoopRow
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine minRestrictionBC2DExactS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
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
	
	Do A=1,NodG
		ApRow  = (A-1)*iDofT + iShiftFluc
		
		Do p = 1 , Ndim2D
			App = ApRow+p
												 
			Do q = 1 , NdofLag
				Bqp = BpCol + q
			
				Coupling(App,Bqp) = 1
				Coupling(Bqp,App) = 1
				Coupling(Bqp,Bqp) = 1
			end do
		end do ! loop Ap
	Enddo !LoopRow

!~ 	
end subroutine





