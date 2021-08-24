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
   
Subroutine tangentHomFibres(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, &
									Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	
	
	use funcAux , only : setAEandBEfib, getSliceAllocate, setBlockToVector, setBlockDiagonalConstantToMatrix
	use globalVariables, only : getF, NdimE, getMaterial, getDamagePar, maxMatPar, maxDamagePar, addToPKhomGen
	use fibresLib
	use ptsGaussLib, only : setNodG
	use fibresMod
	use genericConstitutiveLib, only : SUPD_fib_1D_visc
	use TIME_STEP_ITER
	use DETERMINANT, only : iElem
	use multiscaleNewLib, only : Dhom 
		
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
	real*8 , allocatable ::  SolUf(:) , SolUf0(:), Xel(:), SolUkl(:,:,:) !! all have dimension NodG*NdimE
	integer :: NodG, constLaw , i,j, k,l
	real*8 :: vAux(NdimE), DeltaUkl(NdimE)
	
	iShiftFlucKL = nint(CommonPar(1))
	iShiftUf = nint(CommonPar(2))
	iFemType = nint(commonPar(3))
	iFType = nint(commonPar(4))
	iMaterial = nint(commonPar(5))
	iDamageParType = nint(commonPar(6))
	damageMethod = nint(commonPar(7))
	isAnalytic = nint(commonPar(8))
	x0min = commonPar(9)
	x0max = commonPar(10)
	
    call setNodG(iFemtype, NodG)	
	call getMaterial(constLaw, matPar, iMaterial)
	call getDamagePar(damagePar, iDamageParType)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	call getSliceAllocate(SolUf,Sol1,1,NodG,iShiftUf + 1 ,iShiftUf + NdimE, iDofT)
	call getSliceAllocate(SolUf0,Sol0,1,NodG,iShiftUf + 1 ,iShiftUf + NdimE, iDofT)


	allocate(SolUkl(NdimE*NodG,NdimE,NdimE))
	
	do k = 1, NdimE
	do l = 1, NdimE
		call getSlice(SolUkl(:,k,l),Sol1,1,NodG, iShiftFlucKL+ iShiftKL(k,l,NdimE)*NdimE + 1 , &
		             iShiftFlucKL + iShiftKL(k,l,NdimE)*NdimE +  NdimE, iDofT)
	end do
	end do

	Area = Param(Ipos_Areaf)
	lfa = Param(Ipos_lfa)
	afib = Param(Ipos_af:Ipos_af+1)
	Lfib = Param(Ipos_Lf)
	
!~ 	if (damagePar(3)<0.0d0) then
		damagePar(3) = Param(Ipos_r0)
!~ 	end if
	
	!! 1,2,4,6 is already filled by some other materials parameters initialized in a setMaterial subroutine
	matpar(3) = lfa
	matpar(5) = getLambdaFibreSimple(SolUf0,afib,Lfib,iFtype,NdimE)
	matpar(7) = Param(Ipos_stretch)
	
	call setFibreKinematicSimple(qfib,DeltaUf,afib,Lfib,SolUf,iFtype)
	
	Param(Ipos_stretch) = norm2(qfib)
	
	call calcSfib(Sfib,qfib,matPar,constLaw,NDimE)
	call calcDSfib(DSfib,qfib,matPar,constLaw,NDimE)

	x0 = 0.5*(Xel(1) + Xel(3))
	
	if( .not.( x0 > x0min .and. x0 < x0max) ) then
		damageMethod = 0
	end if
	
	call SUPD_fib_1D_visc(Sfib,DSfib,qfib,damagePar,Param,matPar,Delt,constlaw,NDimE,damageMethod)
	
	if(iElem == 1) then
		Dhom = 0.0d0
	end if
	
	Dhom = Dhom + Area*Lfib*VodotW_fun(DSfib,VotimesW_fun(afib,afib))
	
	do k = 1,NdimE
	do l = 1,NdimE	
	
		DeltaUkl(1) = SolUkl(3,k,l) - SolUkl(1,k,l)
		DeltaUkl(2) = SolUkl(4,k,l) - SolUkl(2,k,l) 
		
		vAux = matmul(DSfib,DeltaUkl)
			
		Dhom(:,:,k,l) = Dhom(:,:,k,l) + Area*VotimesW_fun(vAux,afib)
	
	end do
	end do
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%d%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine tangentHomFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux, only : setCoupling_pureDisplacement, setSimbolicBlockDiagonalToMatrix, numprint
	use ptsGaussLib, only : setNodG
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)

	Coupling(1,1) = 1

end Subroutine



