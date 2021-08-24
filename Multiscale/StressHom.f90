    !> Generic element for homogenise stress
    !!
	!! @param iShiftUf = nint(CommonPar(1))
	!! @param iFemType  = nint(CommonPar(2)) 
	!! @param iFType  = nint(commonPar(3))
	!! @param iMaterial  = nint(commonPar(4))
	!! @param iDamageParType  = nint(commonPar(5))
	
    !! @author Rocha, Felipe Figueredo
    
Subroutine StressHom(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : getF, NdimE, getMaterial, getDamagePar, &
								maxMatPar, maxDamagePar, contributePKhom, writePKhom, allocatePKhom
	
	use fibresLib
	use finiteStrainLib
	use ptsGaussLib , only : setNodG
!~ 	use damageFibLib
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
	Real*8 :: MatPar(maxMatPar),  damagePar(maxDamagePar)
	Real*8 :: Vmu,  PKmu(NdimE,NdimE)
	real*8 , allocatable ::  SolUf(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG, iMaterial, iDamageParType
	integer :: iShiftUf, constLaw , iFemType, iFtype, iElemBegin, Nelem

	iShiftUf = nint(CommonPar(1))
	iFemType = nint(commonPar(2))
	iFType = nint(commonPar(3))
	iMaterial = nint(commonPar(4))
	iDamageParType = nint(commonPar(5))

	if(iShiftUf <0 ) then
		call writePKhom() !! automatically resets PKhom be written
		AE(1,1) = 1.0d0
		BE(1) = Sol1(1)
		return
	end if
	
	call allocatePKhom() !! allocates if necessary
	
    call setNodG(iFemtype, NodG)	
	call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw, matPar, iMaterial)
 	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	call getSliceAllocate(SolUf,Sol1,1,NodG,iShiftUf + 1 ,iShiftUf + NdimE, iDofT)

	if(NodG == 2) then !! lines
		write(0,*) 'Not implemented yet'
!~ 		call getPKmuFromfib(PKmu,Vmu,SolUf,Xel,matpar,Param,constLaw,iFtype,NdimE)
	else !! triangle, quadrangles, etc
		call getPKmuFromMatrix(PKmu,Vmu,SolUf,Xel,matpar,Param,constLaw,iFtype,NdimE,iFemType)
	end if
	
	call contributePKhom(PKmu,Vmu)
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)

end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine stressHomS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)

	Coupling (1, 1) = 1

end Subroutine
