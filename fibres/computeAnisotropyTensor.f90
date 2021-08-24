    !> compute the anisotropy tensor
    !! @param iFemType =  nint(commonPar(1))
	!! @param iMaterial = nint(commonPar(2))
    !! @author Rocha, Felipe Figueredo
Subroutine computeAnisotropyTensor(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, &
									Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : getF, NdimE, getMaterial,maxMatPar, contributeAsisotropyTensor
	use fibresLib
	use ptsGaussLib, only : setNodG
		
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
         
	Real*8 :: MatPar(maxMatPar)
	Real*8 :: afib(NdimE), Lfib, AreaFib, Vfib
	real*8 , allocatable ::  Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG, iMaterial
	integer :: constLaw , iFemType

	iFemType = nint(commonPar(1))
	iMaterial = nint(commonPar(2))
		
    call setNodG(iFemtype, NodG)	
	
	call getMaterial(constLaw, matPar, iMaterial)
	AreaFib = matpar(3)
 	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call get_afib(afib,Lfib,Xel)
	
	Vfib = AreaFib * Lfib
	
	call contributeAsisotropyTensor(afib,Vfib)

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine computeAnisotropyTensorS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Coupling = 0
 
end Subroutine



