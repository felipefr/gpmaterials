

!     ------------------------------------------------------------------
Subroutine NeumannFiniteStrainSpatial(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
						Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use loadingLib
	use finiteStrainLib
	use ptsGaussLib
		
	IMPLICIT NONE
		
	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	! Reals Vectors and matrices
	Real*8  :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), &
		Param(*), JParam(*), CommonPar(*)
	
	!   =====   END ARGUMENTS  =======
     		
	Integer :: nG, i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp , ip, ipp ! counters 
	Integer :: NodG,  pOrder , NGP, iShiftU, iShiftDeltaU, iBu, iSimplex, iFemType
	integer, parameter :: Ndim2D = 2, Ndim1D = 1
	Real*8 :: PhiA, dv , pm,  dPhiA(Ndim2D)
	real*8 , allocatable ::  SolU(:) , Xel(:) !! all have dimension NodG*NdimE
	Real*8 :: normal0(Ndim2D) , vecU(Ndim2D), lengthLine, lengthLine0
	type(ptGaussClass) :: PtG
	
!~ 	pause
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iFemType = nint(CommonPar(3))
	call chooseLoad(nint(CommonPar(4)))
	call pressureLoad(pm,Time,DelT,CommonPar(5:10))
	
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBu)

!~ 	write(0,*) "ps = " , pm
!~ 	pause
			
	call getSliceAllocate(SolU,Sol1,1,NodElt,iShiftU + 1 ,iShiftU + Ndim2D, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodELT,1 ,Ndim2D, Ndim)
	
	VecU = Xel(Ndim2D + 1:2*Ndim2D) - Xel(1:Ndim2D)
	lengthLine0 = dsqrt(dot_product(vecU,vecU))
	

	
	Xel = Xel + SolU
	
	VecU = Xel(Ndim2D + 1:2*Ndim2D) - Xel(1:Ndim2D)
	
	lengthLine = dsqrt(dot_product(vecU,vecU))

!~ 	normal0(1) = vecU(2)/lengthLine
!~ 	normal0(2) = -vecU(1)/lengthLine 
	normal0(1) = vecU(1)/lengthLine
	normal0(2) = vecU(2)/lengthLine 
	
	
	pm = (lengthLine0/lengthLine)*pm
	
!~ 	call numprint(Xel)
!~ 	write(0,*) NGP
			
	call PtG%init(Xel,NodG,Ndim2D,NGP,pOrder,iBu,iSimplex)
	ptG%dV = (lengthLine/sum(PtG%dV))*ptG%dV
!~ 	
	AE = 0.0d0
	BE = 0.0d0
!~ 	
	Do nG = 1, NGP ! LoopGauss

		DV=lengthline 
				
		Do A=1,NodG !LoopRow
			ApRow  = (A-1)*iDofT + iShiftDeltaU
			PhiA = PtG%Phi(A,nG)
			
			Do p = 1 , Ndim2D
				App = ApRow+p

				BE(App) = BE(App) - pm * normal0(p) * PhiA * DV								 
 			end do ! loop Ap
		end do !LoopRow
!~ 		
	end do !LoopGauss
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine NeumannFiniteStrainSpatialS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux
	use ptsGaussLib
	
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
	integer :: NodG, Ndim2D, iFemType
	Integer :: i,j,ipRow,ipCol,k,l, iShiftDeltaU
	
	iShiftDeltaU=nint(CommonPar(2))
	iFemType = nint(commonPar(3)) 
    call setNodG(iFemtype, NodG)
    	
	Ndim2D = 2
	
	Coupling = 0
	do i = 1, NodG
		ipRow = (i-1)*iDofT + iShiftDeltaU
		do j = 1, NodG
			ipCol = (j-1)*iDofT + iShiftDeltaU
			do k = 1,Ndim2D
				do l = 1,Ndim2D
					Coupling(ipRow + k,ipCol + l) =  1
				end do
			end do
		end do
	end do

	call numPrint(Coupling)
!~ 	
end Subroutine




