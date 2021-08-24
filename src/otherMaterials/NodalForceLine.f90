!> NodelForceInLine
!!  @param iShiftU = nint(commonPar(1))  
!!	@param iShiftDeltaU  = nint(CommonPar(2))
!!	@param iFemType = iFemType
!!	@param loadType = nint(CommonPar(4))
!!	@param LoadParam = CommonPar(5:10)
!! @author Rocha, Felipe Figueredo

!     ------------------------------------------------------------------
Subroutine NodalForceLine(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
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
	integer, parameter :: Ndim2D = 2
	Real*8 :: PhiA, dv , pm,  dPhiA(Ndim2D)
	real*8 , allocatable ::  SolU(:) , Xel(:) !! all have dimension NodG*NdimE
	Real*8 :: normal0(Ndim2D) , vecU(Ndim2D), lengthLine, lengthLine0
	type(ptGaussClass) :: PtG
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iFemType = nint(CommonPar(3))
	call chooseLoad(nint(CommonPar(4)))
	call pressureLoad(pm,Time,DelT,CommonPar(5:10))
	
	call setNodG(iFEMtype,NodG)
			
	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + Ndim2D, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,Ndim2D, Ndim)
	
	VecU = Xel(Ndim2D + 1:2*Ndim2D) - Xel(1:Ndim2D)
	lengthLine0 = dsqrt(dot_product(vecU,vecU))
	vecU = vecU/lengthline0

	AE = 0.0d0
	BE = 0.0d0
	
	if(NodG == 2) then
!~ 		write(0,*) "===========================>>>>>>>>>> pm = " , pm 
		BE(iShiftDeltaU + 1) = 0.5d0*pm*lengthLine0*VecU(1) !! in x
		BE(iShiftDeltaU + 2) = 0.5d0*pm*lengthLine0*VecU(2) !! in y
		BE(iDofT + iShiftDeltaU + 1) = 0.5d0*pm*lengthLine0*VecU(1) !! in x	
		BE(iDofT + iShiftDeltaU + 2) = 0.5d0*pm*lengthLine0*VecU(2) !! in y	
	end if

	if(NodG == 3) then
!~ 		write(0,*) "===========================>>>>>>>>>> pm = " , pm 
		BE(iShiftDeltaU + 1) = 0.16666666667d0*pm*lengthLine0*VecU(1) !! in x
		BE(iShiftDeltaU + 2) = 0.16666666667d0*pm*lengthLine0*VecU(2) !! in y
		BE(iDofT + iShiftDeltaU + 1) = 0.66666666667d0*pm*lengthLine0*VecU(1) !! in x	
		BE(iDofT + iShiftDeltaU + 2) = 0.66666666667d0*pm*lengthLine0*VecU(2) !! in y	
		BE(2*iDofT + iShiftDeltaU + 1) = 0.16666666667d0*pm*lengthLine0*VecU(1) !! in x	
		BE(2*iDofT + iShiftDeltaU + 2) = 0.16666666667d0*pm*lengthLine0*VecU(2) !! in y	
	end if
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine NodalForceLineS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
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
 	
end Subroutine



