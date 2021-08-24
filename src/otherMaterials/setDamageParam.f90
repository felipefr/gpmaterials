!> Insert information of damage parameters in the globalvariable Module
!!  @param CommonPar = (a,b,b*[a*ci])  
!!	@param maxDamageParArg  = nint(commonPar(1))  ==> a
!!	@param nbDamagePar = nint(commonPar(2))   ==> b
!!	@param damagePar = CommonPar(see structure) ==> ci (for each damage model)
!!  additional information: @param damagePar = [Integration Scheme, Rtype, Rpar(1:1), Htype, Hpar(1:3)]
!! @author Rocha, Felipe Figueredo

Subroutine setDamageParam &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
	
	use globalVariables
	use funcAux
	use damageNewLib
    
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
!~     integer , save :: ElmNumber = 0
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
	
	integer, parameter :: iShiftDamagePar = 2
    integer ::  i ,j, ip, nbDamagePar, maxDamageParArg, constLaw
    real*8 , allocatable :: damagePar(:) 
   
	maxDamageParArg = nint(commonPar(1))
	nbDamagePar = nint(commonPar(2))

	if(.not. isAllocatedDamagePar() ) then
		call allocateDamagePar(nbDamagePar)
		
		allocate(damagePar(maxDamageParArg))
		do i = 1 , nbDamagePar
			ip = (i-1)*maxDamageParArg + iShiftDamagePar
			damagePar = commonPar( ip + 1 : ip + maxDamageParArg )
			call setDamagePar(damagePar,i)
		end do
		
		
		deallocate(damagePar)
	end if

	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)

end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine setDamageParamS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)

	Coupling(1,1) = 1

end Subroutine
