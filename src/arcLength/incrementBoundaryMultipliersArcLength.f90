	!> Element that solves the trivial diagonal problem of incrementing variables
	!! @param CommonPar = ( a , a*[bi,ci,di,ei] , dummy ) where : 
	!! @param Nvar = CommonPar(1) ==> a
	!! @param iShiftUi = CommonPar(see structure) ===> bi
	!! @param iShiftDeltaUi = CommonPar(see structure) ===> ci
	!! @param nodei (for each material = see structure) ===> ei 
	!! @author Rocha, Felipe Figueredo


!! last update: lenU given by parameter and not by difference between iShiftU and iShiftDeltaU. NodMax is given, rather than assumed
Subroutine incrementBoundaryMultipliersArcLength(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
									Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------

	use globalVariables, only : NdimE
	use ptsGaussLib, only : setNodG
	use funcAux
	use TIME_STEP_ITER
	use DETERMINANT
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
    integer :: i, j , ipDim, ipDof, NodG, iFEMtype
    integer	:: iShift_LM, iShift_dLM_Star, iShift_dLM_Bar, iShift_dLarcLength
	
    real*8 :: dl
    real*8 , dimension(NdimE) :: dLM, LM , dLM_Star , dLM_Bar
	
	iShift_dLarcLength = nint(CommonPar(1))
	iShift_LM = nint(CommonPar(2))
	iShift_dLM_Star = nint(CommonPar(3))
	iShift_dLM_Bar = nint(CommonPar(4))	
	iFEMtype = nint(CommonPar(5))
	
	call setNodG(iFEMtype, NodG)
		
	dl = Sol1(NodG*iDofT + iShift_dLarcLength + 1)
	
	LM = Sol1(iShift_LM + 1 : iShift_LM + NdimE)
	dLM_Star = Sol1(iShift_dLM_Star + 1 : iShift_dLM_Star + NdimE)
	dLM_Bar = Sol1(iShift_dLM_Bar + 1 : iShift_dLM_Bar + NdimE)
	
	dLM = dLM_Star + dl*dLM_Bar
	
	do j = 1, NdimE
		AE( iShift_LM + j, iShift_LM + j ) = 1.0d0					
		BE( iShift_LM + j) = LM(j) +  dLM(j)
		
!~ 		AE( iShift_dLM_Star + j, iShift_dLM_Star + j ) = 1.0d0
!~ 		BE( iShift_LM + j) = 0.0d0	
		
!~ 		AE( iShift_dLM_Bar + j, iShift_dLM_Bar + j ) = 1.0d0
!~ 		BE( iShift_LM + j) = 0.0d0	
	end do

!~ 	if(iElem == 619 ) then
!~ 		call numprint(BE(1:NodELT*iDofT))
!~ 		pause
!~ 	end if

	
end subroutine 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine incrementBoundaryMultipliersArcLengthS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux
	use globalVariables, only : NdimE
	
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
    
    integer	:: j, iShift_LM, iShift_dLM_Star, iShift_dLM_Bar
    
    iShift_LM = nint(CommonPar(2))
    iShift_dLM_Star = nint(CommonPar(3))
	iShift_dLM_Bar = nint(CommonPar(4))	
		
	do j = 1, NdimE
		Coupling( iShift_LM + j, iShift_LM + j ) = 1					
!~ 		Coupling( iShift_dLM_Star + j, iShift_dLM_Star + j ) = 1					
!~ 		Coupling( iShift_dLM_Bar + j, iShift_dLM_Bar + j ) = 1					
	end do

end Subroutine

