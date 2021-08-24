    !> Pressure Load for nonlinear elasticity in a total lagrangian formulation
    !!
    !!
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iShiftDeltaU = nint(CommonPar(2))
	!! @param iFemType  = nint(CommonPar(3)) !!! of the associated volume element (one dimension higher)
	!! @param iLoadProg = nint(CommonPar(4))
	!! @param LoadPar = CommonPar(5:10)
	
!     ------------------------------------------------------------------
Subroutine NeumannFiniteStrain(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
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
	Integer :: NodG,  pOrder , NGP, iShiftU, iShiftDeltaU, iBu, iSimplex, iFemType, iLoadProg
	integer, parameter :: Ndim2D = 2
	Real*8 :: PhiA, dv , pm,  dPhiA(Ndim2D)
	real*8 , allocatable ::  SolU(:) , Xel(:) !! all have dimension NodG*NdimE
	Real*8 :: GradU(Ndim2D,Ndim2D) , F(Ndim2D,Ndim2D), FinvT(Ndim2D,Ndim2D), detF
	Real*8 :: normal0(Ndim2D) , vecU(Ndim2D), FinvTn0(Ndim2D) , lengthLine, LoadPar(6)
	type(ptGaussClass) :: PtG
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iFemType = nint(CommonPar(3))
	iLoadProg = nint(CommonPar(4))
	LoadPar = CommonPar(5:10)
	
	call chooseLoad(iLoadProg)
	call pressureLoad(pm,Time,DelT,LoadPar)
	
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBu)

!~ 	write(0,*) "pm = " , pm
!~ 	pause
			
	call getSliceAllocate(SolU,Sol1,1,NodElt,iShiftU + 1 ,iShiftU + Ndim2D, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodELT,1 ,Ndim2D, Ndim)
	
	VecU = Xel(Ndim2D + 1:2*Ndim2D) - Xel(1:Ndim2D)
	
	lengthLine = dsqrt(dot_product(vecU,vecU))
	normal0(1) = vecU(2)/lengthLine
	normal0(2) = -vecU(1)/lengthLine 
			
	call PtG%initCollapsed(Xel,NodElt,Ndim2D,NGP,pOrder,iBu,iSimplex)
	ptG%dV = (lengthLine/sum(PtG%dV))*ptG%dV
	
	AE = 0.0d0
	BE = 0.0d0
	
	Do nG = 1, NGP ! LoopGauss

		DV=ptG%dV(nG) 
		
		F = 0.0d0
		call PtG%calcGradU(GradU,SolU,nG)
		F = deltaKron(1:Ndim2D,1:Ndim2D) + GradU			
		call matinv(FinvT,detF,transpose(F))
		
		FinvTn0 = matmul(FinvT,normal0)
		
		
		Do A=1,NodG !LoopRow
			ApRow  = (A-1)*iDofT + iShiftDeltaU
			PhiA = PtG%Phi(A,nG)
			dPhiA= PtG%dPhi_G(:,A,nG)
			
			Do p = 1 , Ndim2D
				App = ApRow+p

				BE(App) = BE(App) - detF * pm * FinvTn0(p) * PhiA * DV								 

			end do ! loop Ap
		end do !LoopRow
!~ 		
	end do !LoopGauss
	
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine NeumannFiniteStrainS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
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



