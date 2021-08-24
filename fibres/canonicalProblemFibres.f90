    !> Element implementation for hyperlastic fibres
    !!
    !! \f[
    !!   \sum_f \mathbf{s}_f \cdot \hat{\mathbf{q}}_f 
    !! \f]
    !!
	!! @param iShiftUf = nint(CommonPar(1))
	!! @param iShiftDeltaU = nint(CommonPar(2)) !!! now it's not being used
	!! @param iFemType  = nint(CommonPar(3)) 
	!! @param iFType  = nint(commonPar(4))
	!! @param iMaterial  = nint(commonPar(5))
	!! @param iDamageParType  = nint(commonPar(6))
	!!
    !! @author Rocha, Felipe Figueredo
   
Subroutine canonicalProblemFibres(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, &
									Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	
	
	use funcAux , only : setAEandBEfib, getSliceAllocate, setBlockToVector, setBlockDiagonalConstantToMatrix
	use globalVariables, only : getF, NdimE, getMaterial, getDamagePar, maxMatPar, maxDamagePar, addToPKhomGen
	use fibresLib
	use ptsGaussLib, only : setNodG
	use fibresMod
	use genericConstitutiveLib, only : SUPD_fib_1D_visc
		
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)
!~ 	!   =====   END ARGUMENTS  =======
    
	! ======== COMMONPAR ARGUMENTS =====================
	integer :: iShiftUf, iShiftFlucKL, iFemType, iFType, iMaterial, iDamageParType, kl, damageMethod, isAnalytic
	
	!====== OTHER VARIABLES ==============================
	Real*8 :: MatPar(maxMatPar),  damagePar(maxDamagePar)
	Real*8 :: DeltaUf(NdimE), afib(NdimE), qfib(NdimE), Lfib , Kfib(NdimE,NdimE), b(NdimE)
	Real*8 :: DSfib(NdimE,NdimE), Sfib(NdimE), signal, Area, lfa, energy, x0, x0min, x0max
	real*8 , allocatable ::  SolUf(:) , SolUf0(:), Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG, constLaw , k,l
	real*8 :: ek(NdimE), el(NdimE)
	
	iShiftFlucKL = nint(CommonPar(1))
	iShiftUf = nint(CommonPar(2))
	kl = nint(CommonPar(3))
	iFemType = nint(commonPar(4))
	iFType = nint(commonPar(5))
	iMaterial = nint(commonPar(6))
	iDamageParType = nint(commonPar(7))
	damageMethod = nint(commonPar(8))
	isAnalytic = nint(commonPar(9))
	x0min = commonPar(10)
	x0max = commonPar(11)
	
	call getIJfromK(k,l,kl)
	
	ek = 0.0d0
	el = 0.0d0
	ek(k) = 1.0d0
	el(l) = 1.0d0
	
    call setNodG(iFemtype, NodG)	
	call getMaterial(constLaw, matPar, iMaterial)
	call getDamagePar(damagePar, iDamageParType)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	call getSliceAllocate(SolUf,Sol1,1,NodG,iShiftUf + 1 ,iShiftUf + NdimE, iDofT)
	call getSliceAllocate(SolUf0,Sol0,1,NodG,iShiftUf + 1 ,iShiftUf + NdimE, iDofT)

	Area = Param(Ipos_Areaf)
	lfa = Param(Ipos_lfa)
	afib = Param(Ipos_af:Ipos_af+1)
	Lfib = Param(Ipos_Lf)
	
	!! 1,2,4,6 is already filled by some other materials parameters initialized in a setMaterial subroutine
	matpar(3) = lfa
	matpar(5) = getLambdaFibreSimple(SolUf0,afib,Lfib,iFtype,NdimE)
	matpar(7) = Param(Ipos_stretch)
	
!~ 	if (damagePar(3)<0.0d0) then
		damagePar(3) = Param(Ipos_r0)
!~ 	end if
	
	call setFibreKinematicSimple(qfib,DeltaUf,afib,Lfib,SolUf,iFtype)
	
	Param(Ipos_stretch) = norm2(qfib)
	
	call calcSfib(Sfib,qfib,matPar,constLaw,NDimE)
	call calcDSfib(DSfib,qfib,matPar,constLaw,NDimE)

	x0 = 0.5*(Xel(1) + Xel(3))
	
	if( .not.( x0 > x0min .and. x0 < x0max) ) then
		damageMethod = 0
	end if
	
!~ 	damagePar(3) = 3.00001d0
	call SUPD_fib_1D_visc(Sfib,DSfib,qfib,damagePar,Param,matPar,Delt,constlaw,NDimE,damageMethod)

	b = -Area*matmul(DSfib,ek)*dot_product(el,afib)	
	Kfib = (Area/Lfib)*DSfib
	
	call setAEandBEfib(AE,BE,Kfib,b,NodG,idofT,iShiftFlucKL,NdimE)
		
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%d%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine canonicalProblemFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux, only : setCoupling_pureDisplacement, setSimbolicBlockDiagonalToMatrix, numprint
	use ptsGaussLib, only : setNodG
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol
    integer , parameter :: NdimE = 2
    integer :: iShiftFlucKL, iFemType, NodG

	iShiftFlucKL = nint(CommonPar(1))
	iFemType = nint(commonPar(4))
	     
    call setNodG(iFemtype, NodG)

	call setCoupling_pureDisplacement(Coupling,NodG,idofT,iShiftFlucKL,NdimE)

end Subroutine



