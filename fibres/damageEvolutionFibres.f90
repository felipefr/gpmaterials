!Last Modification : 11/03/2014
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Element for Hyperelasticity 
! Still in development
!     ------------------------------------------------------------------

subroutine damageEvolutionFibres(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	
	use funcAux
	use globalVariables, only : NdimE, getDamagePar, getMaterial, maxDamagePar, maxMatPar, getF
	use fibresMod
	use fibresLib
	use ptsGaussLib, only : setNodG
	use genericConstitutiveLib, only : SUPD_fib

	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
         
	Real*8 :: MatPar(maxMatPar),  damagePar(maxDamagePar)
	real*8 , allocatable ::  SolUf(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG, iShiftUf, constLaw, iFemType, iMaterial, iDamageParType , iFtype, damageMethod, isAnalytic
	real*8 :: lambda, energy, x0, x0min, x0max, lfa, & 
			DSFib(NdimE,NdimE), SFib(NdimE), qfib(NdimE) , F(NdimE,NdimE), afib(NdimE), Lfib

	iShiftUf = nint(commonPar(1))
	iFemType = nint(commonPar(2))
	iFtype = nint(commonPar(3))
	iMaterial = nint(commonPar(4))	
	iDamageParType = nint(commonPar(5))
	damageMethod = nint(CommonPar(6))
	isAnalytic = nint(CommonPar(7))
	x0min = commonPar(8)
	x0max = commonPar(9)
	
    call setNodG(iFemtype, NodG)	
	call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw,matPar,iMaterial)	

	call getSliceAllocate(SolUf,Sol1,1,NodG,iShiftUf + 1 ,iShiftUf + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
!~ 	x0 = 0.5d0*(Xel(1) + Xel(3))
	
	afib = Param(Ipos_af:Ipos_af+1)
	Lfib = Param(Ipos_Lf)
	lfa = Param(Ipos_lfa)
	
	!! 1,2 filled by some other materials parameters initialized in a setMaterial subroutine
	matpar(3) = lfa

	
!~ 	if(x0min < 0.0d0 .or. (x0 > x0min .and. x0 < x0max)) then
!~ 		energy = getEnergyFibre(SolUf,Xel,MatPar,constLaw,iFtype,NdimE)
!~ 		write(0,*) '============ energy = ' , energy
!~ 		call damageUpdateFib(damagePar,Param , energy)	
!~ 		call damageUpdateFibOld(damagePar,Param , energy)		
!~ 	end if

	call getF(F,iFtype)
	call setFibreKinematicSimple2(qfib,F,SolUf,afib,Lfib)
	
	call SUPD_fib(Sfib,DSfib,qfib,damagePar,Param,matPar,constlaw,NDimE,damageMethod,isAnalytic)
		
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine damageEvolutionFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow , i
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)

	Coupling=0 
	
	Coupling(1,1) = 1
	
	return
end subroutine
