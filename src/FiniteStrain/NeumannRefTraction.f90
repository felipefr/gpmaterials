    !> Pressure Load for nonlinear elasticity in a total lagrangian formulation
    !!
    !!
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iShiftDeltaU = nint(CommonPar(2))
	!! @param iFemType  = nint(CommonPar(3)) !!! of the associated volume element (one dimension higher)
	!! @param iLoadProg = nint(CommonPar(4))
	!! @param LoadPar = CommonPar(5:10)
	
!     ------------------------------------------------------------------
Subroutine NeumannRefTraction(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
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
     		
	Integer :: iLoadProg, iShiftDeltaU, iShiftU
	integer, parameter :: NdimE = 2
	real*8 , allocatable ::  Xel(:) !! all have dimension NodG*NdimE
	Real*8 :: normal0(NdimE) , vecU(NdimE), lengthLine, LoadPar(6), pm
	Real*8 :: Traction(2)
	type(ptGaussClass) :: PtG
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iLoadProg = nint(CommonPar(3))
	LoadPar = CommonPar(4:9)
	
	call chooseLoad(iLoadProg)
	call pressureLoad(pm,Time,DelT,LoadPar)
			
	call getSliceAllocate(Xel,XLL,1,NodELT,1 ,NdimE, Ndim)
	
	VecU = Xel(NdimE + 1 : 2*NdimE) - Xel(1:NdimE)
	
	lengthLine = dsqrt(dot_product(vecU,vecU))
	normal0(1) = vecU(2)/lengthLine
	normal0(2) = -vecU(1)/lengthLine 
	
	write(0,*) "============================= > pm = " , pm 
	Traction = -pm*normal0
	
	AE = 0.0d0
	BE = 0.0d0
	
	BE(iShiftDeltaU + 1 ) = 0.5d0*lengthLine*Traction(1)								 
	BE(iShiftDeltaU + 2 ) = 0.5d0*lengthLine*Traction(2)
	
	BE(iDofT + iShiftDeltaU + 1 ) = 0.5d0*lengthLine*Traction(1)
	BE(iDofT + iShiftDeltaU + 2 ) = 0.5d0*lengthLine*Traction(2)
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine NeumannRefTractionS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux
	use ptsGaussLib
	
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
	Integer :: iShiftDeltaU
	
	iShiftDeltaU=nint(CommonPar(2))

	Coupling(iShiftDeltaU + 1, iShiftDeltaU + 1 ) = 1 								 
	Coupling(iShiftDeltaU + 2, iShiftDeltaU + 2 ) = 1

	Coupling(iDofT + iShiftDeltaU + 1, iDofT + iShiftDeltaU + 1 ) = 1 								 
	Coupling(iDofT + iShiftDeltaU + 2, iDofT + iShiftDeltaU + 2 ) = 1	
	
end Subroutine



