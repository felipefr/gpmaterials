module damageNewLib
use funcAux
use damageMemoryManagementLib
use HyperModelsLib
use fibresModelsLib
use damageModelsNewLib

implicit none

contains

subroutine modifyStressAndTangent(Ppk,D, damagePar,Param,NDimE,damageMethod)
	real*8, intent(inout) :: Ppk(:,:),D(:,:,:,:) 
	real*8, intent(in) :: damagePar(:) , Param(*)
	integer , intent(in) :: NdimE, damageMethod
	real*8 :: PpktenPpk(NdimE,NdimE,NdimE,NdimE) , damage, derDamage
	
	call decideDamageModel(damage, derDamage , damagePar,Param,damageMethod)
	
	call VotimesW(PpktenPpk,Ppk,Ppk)
	
	D = (1.0 - damage)*D - derDamage*PpktenPpk 	
	Ppk = (1.0 - damage)*Ppk

end subroutine


subroutine modifyStressAndTangent_fib(Sfib,DSfib, damagePar,Param,NDimE,damageMethod)
	real*8, intent(inout) :: Sfib(:),DSfib(:,:) 
	real*8, intent(in) :: damagePar(:) , Param(*)
	integer , intent(in) :: NdimE, damageMethod
	real*8 :: SfibtenSfib(NdimE,NdimE) , damage, derDamage
	
	
	call decideDamageModel(damage, derDamage , damagePar,Param,damageMethod)
	
	call VotimesW(SfibtenSfib,Sfib,Sfib)
	
	DSfib = (1.0 - damage)*DSfib - derDamage*SfibtenSfib 	
	Sfib = (1.0 - damage)*Sfib

end subroutine

subroutine modifyStressAndTangent_fib_1D(sfib,dsfib, sfib_e, damagePar,Param, damageMethod)
	real*8, intent(inout) :: sfib,dsfib 
	real*8, intent(in) :: damagePar(:) , Param(*), sfib_e
	integer , intent(in) :: damageMethod
	real*8 :: damage, derDamage
	
	call decideDamageModel(damage, derDamage , damagePar,Param,damageMethod)
	
	dsfib = (1.0 - damage)*dsfib - derDamage*sfib*sfib_e 	
	sfib = (1.0 - damage)*sfib

end subroutine

subroutine stateVariableUpdate_fib(qfib,damagePar,Param,matPar, constlaw)
	real*8, intent(in) :: qfib(:), damagePar(:), matPar(:)
	real*8, intent(inout) :: Param(*)
	integer , intent(in) :: constLaw
	type(damageState) :: sd, sdNew
	real*8 :: H, Hder, r, r0, rtheta, theta, ftrial,energy0, rAux	

	theta = damagePar(1)
	r0 = damagePar(3)
	
	call strainEnergy_fib(energy0,qfib,matPar,constLaw)
	
	call sd%loadDamage(Param, 1, .false.)	
	
	call fixInitialVariables(sd,sdNew,r0,Param)
	
	ftrial = dsqrt(2.0d0*energy0) - sd%rd
	
	if(ftrial <= 0.0d0) then
		sdNew%rd = sd%rd
		sdNew%rdDot = 0.0d0
		sdNew%qd = sd%qd
		sdNew%damage = sd%damage 
!~ 		sd%derDamage = 0.0d0
		sdNew%derDamage = 0.0d0 
	else 
		sdNew%rd = sd%rd + ftrial
		sdNew%rdDot = ftrial
		
		rtheta = theta*sdNew%rd + (1.0d0 - theta)*sd%rd
		
		call getH(H,Hder,rtheta,damagePar)
		
		sdNew%qd = sd%qd + H*(sdNew%rd - sd%rd)
		
		if(sdNew%qd < 0.0d0) sdNew%qd = 1.0d-15
		
		sdNew%damage = 1.0d0 - sdNew%qd/sdNew%rd !! integrating	
		
		sdNew%derDamage = (sdNew%qd - sdNew%rd*(H + theta*Hder*(sdNew%rd - sd%rd) ))/(sdNew%rd**3.0d0)		
		
!~ 		sdNew%derDamage = (sdNew%qd - sdNew%rd*H)/(sdNew%rd**3.0d0)	! wrong		
	end if
	
	call sdNew%putDamage(Param, 1)

end subroutine


subroutine stateVariableUpdate(F,damagePar,Param,matPar, constlaw)
	use fibresMod , only: Ipos_r0

	real*8, intent(in) :: F(:,:), damagePar(:), matPar(:)
	real*8, intent(inout) :: Param(*)
	integer , intent(in) :: constLaw
	type(damageState) :: sd, sdNew
	real*8 :: H, Hder, r, r0, rtheta, theta, ftrial,energy0, rAux	

	theta = damagePar(1)
	r0 = damagePar(3)

	call strainEnergy(energy0,F,matPar,constLaw)
	
	call sd%loadDamage(Param, 1, .false.)	
	
	call fixInitialVariables(sd,sdNew,r0,Param)
	
	ftrial = dsqrt(2.0d0*energy0) - sd%rd
	
	if(ftrial <= 0.0d0) then
		sdNew%rd = sd%rd
		sdNew%rdDot = 0.0d0
		sdNew%qd = sd%qd
		sdNew%damage = sd%damage 
		sdNew%derDamage = 0.0d0 
	else 
		sdNew%rd = sd%rd + ftrial
		sdNew%rdDot = ftrial
		
		rtheta = theta*sdNew%rd + (1.0d0 - theta)*sd%rd
		
		call getH(H,Hder,rtheta,damagePar)
		
		sdNew%qd = sd%qd + H*(sdNew%rd - sd%rd)
		
		if(sdNew%qd < 0.0d0) sdNew%qd = 1.0d-15
		
		sdNew%damage = 1.0d0 - sdNew%qd/sdNew%rd !! integrating	
		
		sdNew%derDamage = (sdNew%qd - sdNew%rd*(H + theta*Hder*(sdNew%rd - sd%rd) ))/(sdNew%rd**3.0d0)		
		
!~ 		sdNew%derDamage = (sdNew%qd - sdNew%rd*H)/(sdNew%rd**3.0d0)	! wrong		
	end if
	
	call sdNew%putDamage(Param, 1)

end subroutine
	
real*8 function getImplexDamage(rd,rddot,qd,damagePar) result(damageTilde)
	real*8, intent(in) :: rd,rddot,qd,damagePar(:)
	real*8 :: rdtilde, qdtilde, H, Hder
	
	rdtilde = rd + rdDot
	call getH(H,Hder,rd,damagePar)
	qdTilde = qd + H*(rdtilde - rd)
	damageTilde = 1.0d0 - qdTilde/rdTilde
	
end function

subroutine fixInitialVariables(sd,sdNew,r0,param)
	type(damageState) , intent(inout) :: sd, sdNew
	real*8 , intent(out) :: Param(*)
	real*8 , intent(in) :: r0
	
	if(r0>sd%rd) then !!! because r should begin in r0
		sdNew%rd = r0
		sd%rd = r0
		sdNew%qd = r0
		sd%qd = r0
		sdNew%derDamage = 0.0d0
		sd%derDamage = 0.0d0
		sdNew%damage = 0.0d0 
		sd%damage = 0.0d0	
		sdNew%rdDot = 0.0e0 
		sd%rdDot = 0.0e0
		call sd%putDamageIntoOld(Param, 1)
		call sdNew%putDamage(Param, 1)
	end if

end subroutine

subroutine decideDamageModel(damage, derDamage , damagePar,Param,damageMethod)
	real*8, intent(out) :: damage, derDamage
	real*8, intent(in) :: damagePar(:) , Param(*)
	integer , intent(in) :: damageMethod
	
	type(damageState) :: sdNew , sd

	call sdNew%loadDamage(Param, 1, .true.)
	call sd%loadDamage(Param, 1, .false.)
	
	select case(damageMethod)
	
		case(1) !! implicit
			derDamage = sdNew%derDamage
			damage = sdNew%damage
			
		case(2) !! explicit
			derDamage = 0.0d0
			damage = sd%damage
		
		case(3) !! implex
			derDamage = 0.0d0
			damage = getImplexDamage(sd%rd,sd%rddot,sd%qd,damagePar)
			
		case(4) !! wrong case, put the one converging
			derDamage = sd%derDamage
			damage = sdNew%damage
		
		case default  !! no damage at all
			derDamage = 0.0d0
			damage = 0.0d0
	
	end select


end subroutine
	
end module

