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
	!!
    !! @author Rocha, Felipe Figueredo
   
Subroutine nonlinearFibres(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, &
									Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	
	
	use funcAux , only : setAEandBEfib, getSliceAllocate, setBlockDiagonalConstantToMatrix, setBlockToVector
	use globalVariables, only : getF, NdimE, getMaterial, getDamagePar, maxMatPar, maxDamagePar, addToPKhomGen
	use fibresLib
	use ptsGaussLib, only : setNodG
	use fibresMod
	use genericConstitutiveLib, only : SUPD_fib
		
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
	
	!====== OTHER VARIABLES ==============================
	Real*8 :: MatPar(maxMatPar),  damagePar(maxDamagePar)
	Real*8 :: DeltaU(NdimE), afib(NdimE), qfib(NdimE), Lfib, Res(NdimE) , K(NdimE,NdimE)
	Real*8 :: DSfib(NdimE,NdimE), Sfib(NdimE), signal, Area, lfa
	real*8 , allocatable ::  SolU(:) , SolU0(:), Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG
	integer :: constLaw 

	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2)) !!! now it's not being used
	iFemType = nint(commonPar(3))
	iFType = nint(commonPar(4))
	iMaterial = nint(commonPar(5))
	
    call setNodG(iFemtype, NodG)	
	call getMaterial(constLaw, matPar, iMaterial)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(SolU0,Sol0,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)

	Area = Param(Ipos_Areaf)
	lfa = Param(Ipos_lfa)
	afib = Param(Ipos_af:Ipos_af+1)
	Lfib = Param(Ipos_Lf)
	
	!! 1,2,4,6 is already filled by some other materials parameters initialized in a setMaterial subroutine
	matpar(3) = lfa
	matpar(5) = getLambdaFibreSimple(SolU0,afib,Lfib,iFtype,NdimE)
	matpar(7) = Param(Ipos_stretch)
	 
	
	call setFibreKinematicSimple(qfib,DeltaU,afib,Lfib,SolU,iFtype)
	
!~ 	Param(Ipos_stretch) = norm2(qfib)
	
!~ 	write(*,*) 'stretch == ' , norm2(qfib)
	
	call calcSfib(Sfib,qfib,matPar,constLaw,NDimE)
	call calcDSfib(DSfib,qfib,matPar,constLaw,NDimE)

	Res = -Area*Sfib	
	K = (Area/Lfib)*DSfib
	
	if(iShiftU == iShiftDeltaU ) Res = Res + matmul(K,DeltaU)

	call setAEandBEfib(AE,BE,K,Res,NodG,idofT,iShiftDeltaU,NdimE)
	
	 if(iShiftU /= iShiftDeltaU) then
		call setBlockDiagonalConstantToMatrix(AE,1.0d0,NodG, iDofT, NdimE, iShiftU, iShiftU)
		call setBlockDiagonalConstantToMatrix(AE,-1.0d0,NodG, iDofT, NdimE, iShiftU, iShiftDeltaU)
		call setBlockToVector(BE,SolU,NodG,iDofT, NdimE, iShiftU ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	end if
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine nonlinearFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux, only : setCoupling_pureDisplacement, numprint, setSimbolicBlockDiagonalToMatrix
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

	call setCoupling_pureDisplacement(Coupling,NodG,idofT,iShiftDeltaU,NdimE)
	
	if(iShiftU /= iShiftDeltaU) then
		call setSimbolicBlockDiagonalToMatrix(Coupling,NodG, iDofT, NdimE, iShiftU, iShiftU)
		call setSimbolicBlockDiagonalToMatrix(Coupling,NodG, iDofT, NdimE, iShiftU, iShiftDeltaU)
	end if

end Subroutine



