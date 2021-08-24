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
   
Subroutine NonLinFibresGen(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, &
									Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	
	
	use funcAux , only : setAEandBEfib, getSliceAllocate
	use globalVariables, only : getF, NdimE, getMaterial, getDamagePar, maxMatPar, maxDamagePar, addToPKhomGen
	use fibresLib
	use ptsGaussLib, only : setNodG
	use fibresMod
	use genericConstitutiveLib
		
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
	Real*8 :: afib(NdimE), qfib(NdimE), Lfib, Res(NdimE) , K(NdimE,NdimE)
	Real*8 :: DSfib(NdimE,NdimE), Sfib(NdimE), F(NdimE,NdimE), Area, lfa, energy, x0, DeltaU(NdimE)
	real*8 , allocatable ::  SolU(:) , SolU0(:), Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG
	integer :: constLaw , damageMethod , isAnalytic

	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2)) !!! now it's not being used
	iFemType = nint(commonPar(3))
	iFType = nint(commonPar(4))
	iMaterial = nint(commonPar(5))
	iDamageParType = nint(commonPar(6))
	damageMethod = nint(commonPar(7))
	isAnalytic = nint(commonPar(8))
	
    call setNodG(iFemtype, NodG)	
	call getMaterial(constLaw, matPar, iMaterial)
	call getDamagePar(damagePar, iDamageParType)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(SolU0,Sol0,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)

	Area = Param(Ipos_Areaf)
	afib = Param(Ipos_af:Ipos_af+1)
	Lfib = Param(Ipos_Lf)
	lfa = Param(Ipos_lfa)
	matpar(3) = lfa
	
	call getF(F,iFtype)
	
	call setFibreKinematicSimple2(qfib,F,SolU,afib,Lfib)
	call SUPD_fib(Sfib,DSfib,qfib,damagePar,Param,matPar,constlaw,NDimE,damageMethod , isAnalytic) 
	
	Res = -Area*Sfib	
	K = (Area/Lfib)*DSfib
	
	DeltaU = (qfib - matmul(F,afib))*Lfib
	
	Res = Res + matmul(K,DeltaU)
	
	call setAEandBEfib(AE,BE,K,Res,NodG,idofT,iShiftU,NdimE)
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine NonLinFibresGenS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
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
	iShiftDeltaU = nint(CommonPar(2)) !! not being used at the moment
	iFemType = nint(commonPar(3))
	     
    call setNodG(iFemtype, NodG)

	call setCoupling_pureDisplacement(Coupling,NodG,idofT,iShiftU,NdimE)

end Subroutine



