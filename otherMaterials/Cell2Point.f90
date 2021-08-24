Subroutine Cell2Point(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use ptsGaussLib
	use globalVariables, only : NdimE

	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodElt ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
         
    integer , parameter :: maxNFields = 5, iShiftMapping = 2
	Integer :: NodG, NGP, ipPar, ipSol, i, A, ip ! counters 
	real*8 , allocatable :: Xel(:) 
	integer ::  NFields , iCell2PointMap(2,maxNFields)
    integer :: iFemType, pOrder, iSimplex,  iBubble, iFtype
    real*8 :: dVol 
	type(ptGaussClass) :: PtG
	
	iFemType = nint(commonPar(1))
	NFields = nint(CommonPar(2))
	
	do i = 1 , NFields
		ip = iShiftMapping + (i-1)*2
		iCell2PointMap(1,i) = nint(commonPar(ip+1))
		iCell2PointMap(2,i) = nint(commonPar(ip+2))
	end do
	
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)

	!! supposing NGP == 1
	dVol=dabs(ptG%dV(1))

	do i = 1, NFields
		ipPar = iCell2PointMap(1,i)

		do A = 1 , NodG
			ipSol = (A-1)*iDofT + iCell2PointMap(2,i)
			
			AE(ipSol,ipSol) = dVol 
			BE(ipSol) = dVol*dabs(Param(ipPar))
		end do

	end do

End Subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine Cell2PointS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use ptsGaussLib , only : setNodG
	
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
	
	integer , parameter :: maxNFields = 5, iShiftMapping = 2
	integer :: NodG, ipSol, i, A, ip ! counters 
	integer ::  NFields , iCell2PointMap(2,maxNFields)
    integer :: iFemType
    
	iFemType = nint(commonPar(1))
	NFields = nint(CommonPar(2))
	
	do i = 1 , NFields
		ip = iShiftMapping + (i-1)*2
		iCell2PointMap(1,i) = nint(commonPar(ip+1))
		iCell2PointMap(2,i) = nint(commonPar(ip+2))
	end do
	
	call setNodG(iFemtype, NodG)
	
	do i = 1, NFields
	
		do A = 1 , NodG
			ipSol = (A-1)*iDofT + iCell2PointMap(2,i)
			
			Coupling(ipSol,ipSol) = 1 
	
		end do

	end do	
	
end Subroutine
