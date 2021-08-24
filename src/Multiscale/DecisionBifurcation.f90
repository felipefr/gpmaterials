!     ------------------------------------------------------------------
Subroutine DecisionBifurcation(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use finiteStrainLib !! New Way	
	use ptsGaussLib
	use DETERMINANT , only : lengthParam
	use multiscaleLib, only : getMinDetQ
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
         
	Integer :: nG , ApRow, i, j ! counters 
	integer, parameter :: NdimE = 2 ! this code only makes sense in 2D 
	Real*8 :: GradU(NdimE,NdimE)
	Real*8 , allocatable ::  SolU(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: iShiftDeltaU , iShiftFlag, Nelem , PosDetQ, PosBifDomain, iShiftC
	integer :: pOrder , NGP, NodG, iFEMtype , iSimplex, iBubble , nStepTheta
	real*8 :: normal(NdimE), beta(NdimE), alpha, theta, bifFac
	real*8 , save :: minDetQ = 0.0d0 ,thetaCrit = 0.0d0, betaCrit = 0.0d0
	type(ptGaussClass) :: PtG
	integer , save :: iElem = 0

	iShiftDeltaU = nint(CommonPar(1))
	iShiftC = nint(CommonPar(2))
	iFEMtype = nint(CommonPar(3))
	NElem = nint(commonPar(4))
	PosBifDomain = nint(CommonPar(5)) ! in the Param
	PosDetQ = nint(CommonPar(6)) ! in the last node
	
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)

	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftDeltaU + 1 ,iShiftDeltaU + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)

	NGP = 1
	nG = 1
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
		
	call PtG%calcGradU(GradU,SolU,nG)

	bifFac = 0.0d0
	iElem = mod(iElem + 1,Nelem)
	
	if(iElem == 0) then
		nStepTheta = 100
		call getMinDetQ(minDetQ,thetaCrit,betaCrit,Sol1,NdimE,idofT,NodElt,iShiftC,nStepTheta)
	end if
	
	if(minDetQ < 0.0d0) then
		alpha = 0.0d0*PI 
		theta = datan(0.3d0/1.0d0)
		normal(1) = -dsin(theta)
		normal(2) = dcos(theta)
		beta(1) = dcos(theta + alpha)
		beta(2) = dsin(theta + alpha)
			
		do i = 1, NdimE
		do j = 1, NdimE
			BifFac = BifFac + GradU(i,j)*beta(i)*normal(j)
		end do
		end do
		
	end if
	
	Param(LengthParam + PosBifDomain ) = BifFac
			
	ApRow  = NodG*iDofT + PosDetQ
	AE(ApRow,ApRow) = 1.0d0
	BE(ApRow) = minDetQ

	deallocate(SolU,Xel)
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine DecisionBifurcationS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux
	use ptsGaussLib
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: ApRow, iShiftFlag, iFemType, NodG, PosDetQ
	integer, parameter :: NdimE = 2

	iFEMtype = nint(CommonPar(3))
	PosDetQ = nint(CommonPar(6)) ! in the last node

	call setNodG(iFemType, NodG)
 
	ApRow  = NodG*iDofT + PosDetQ
			
	Coupling(ApRow,ApRow) =  1
	
end Subroutine

