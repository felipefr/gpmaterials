    !> Element implementation for hyperlastic fibres
    !!
    !! \f[
    !!   \sum_f \mathbf{s}_f \cdot \hat{\mathbf{q}}_f 
    !! \f]
    !!
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iShiftDeltaU = nint(CommonPar(2)) !!! now it's not being used
	!! @param iFemType  = nint(CommonPar(3)) 
	!! @param iFType  = nint(commonPar(4))
	!! @param iMaterial  = nint(commonPar(5))
	!! @param iDamageParType  = nint(commonPar(6))
	!!
    !! @author Rocha, Felipe Figueredo
   
Subroutine nonlinearFibresQuang(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, &
									Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	
	
	use funcAux , only : setAEandBEfib, getSliceAllocate
	use globalVariables, only : getF, NdimE, getMaterial, getDamagePar, maxMatPar, maxDamagePar, addToPKhomGen
	use fibresLib
	use ptsGaussLib, only : setNodG
	use fibresMod
	use genericConstitutiveLib, only : SUPD_fib, SUPD_fib_1D
		
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)
!~ 	!   =====   END ARGUMENTS  =======
    
	! ======== COMMONPAR ARGUMENTS =====================
	integer iShiftU !> @param  = nint(CommonPar(1))
	integer iShiftDeltaU !> @var  = nint(CommonPar(2)) !!! now it's not being used
	integer iFemType !> @var  = nint(CommonPar(2)) = nint(commonPar(3))
	integer iFType !> @var  = nint(commonPar(4))
	integer iMaterial !> @var  = nint(commonPar(5))
	integer iDamageParType !> @var = nint(commonPar(6))
	
	!====== OTHER VARIABLES ==============================
	Real*8 :: MatPar(maxMatPar),  damagePar(maxDamagePar)
	Real*8 :: DeltaU(NdimE), afib(NdimE), qfib(NdimE), Lfib, Res(NdimE) , K(NdimE,NdimE)
	Real*8 :: DSfib(NdimE,NdimE), Sfib(NdimE), signal, Area, lfa, energy, x0
	real*8 , allocatable ::  SolU(:) , SolU0(:), Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG
	integer :: constLaw , damageMethod, isAnalytic

	iShiftU = nint(CommonPar(1))
	iFemType = nint(commonPar(2))
	iMaterial = nint(commonPar(3))
	isAnalytic = 1
	
    call setNodG(iFemtype, NodG)	
	call getMaterial(constLaw, matPar, iMaterial)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(SolU0,Sol0,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)

	lfa = 1.0d0
	Area = 1.0d0
	
	iFtype = 0
	
	!! 1,2,4,6 is already filled by some other materials parameters initialized in a setMaterial subroutine
	matpar(3) = lfa
	
	call setFibreKinematic(qfib,DeltaU,afib,Lfib,SolU,Xel,iFtype)
	
	Param(Ipos_stretch) = norm2(qfib)
	
!~ 	call SUPD_fib(Sfib,DSfib,qfib,damagePar,Param,matPar,constlaw,NDimE,damageMethod,isAnalytic)
	call SUPD_fib_1D(Sfib,DSfib,qfib,damagePar,Param,matPar,constlaw,NDimE,damageMethod)

	Res = -Area*Sfib	
	K = (Area/Lfib)*DSfib
	
	Res = Res + matmul(K,DeltaU)
	
	call setAEandBEfib(AE,BE,K,Res,NodG,idofT,iShiftU,NdimE)
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine nonlinearFibresQuangS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux, only : setCoupling_pureDisplacement, numprint
	use ptsGaussLib, only : setNodG
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol
    integer , parameter :: NdimE = 2
    integer :: iShiftU, iShiftDeltaU, iFemType, NodG

	iShiftU = nint(CommonPar(1))
	iFemType = nint(commonPar(2))
	     
    call setNodG(iFemtype, NodG)

	call setCoupling_pureDisplacement(Coupling,NodG,idofT,iShiftU,NdimE)

end Subroutine


