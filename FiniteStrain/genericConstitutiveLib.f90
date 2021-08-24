module genericConstitutiveLib
use funcAux
use finiteStrainLib
use fibresLib
use damageNewLib

implicit none

public SUPD_continuum, SUPD_fib, SUPD_fib_1D, SUPD_fib_1D_visc

contains

subroutine SUPD_continuum(Ppk,D,F,damagePar,Param,matPar,constlaw,NDimE,damageMethod,isAnalytic) 
	real*8, intent(in) :: F(:,:), damagePar(:), matPar(:)
	real*8, intent(out) :: Ppk(:,:),D(:,:,:,:) 
	real*8, intent(inout) :: Param(*)
	integer , intent(in) :: NdimE, constLaw, damageMethod, isAnalytic

	if(isAnalytic == 0) then
		call calcPpk(Ppk,F,NdimE,matPar,constLaw)		
		call calcD(D,F,NdimE,matPar,constLaw)
	
	else if(isAnalytic == 1) then
		if(constLaw == 5) then
			call calcForCiarlet(Ppk,D,F,matPar,NdimE)
		else
			call calcPpkD_analytic(Ppk,D,F,matPar,constLaw,NdimE)
		end if
		
	else if(isAnalytic == 2) then ! Analytic in Holzapfel way
		call calcPpkD_analytic_Holzapfel(Ppk,D,F,matPar,constLaw,NdimE)
	
	else
		write(0,*) 'SUPD_continuum, option not found for isAnalytic ' , isAnalytic 
	end if
	
	if(damageMethod>0) then ! it allows damage
	
		call stateVariableUpdate(F,damagePar,Param,matPar, constlaw)
		call modifyStressAndTangent(Ppk,D, damagePar,Param,NDimE,damageMethod)
		
	end if

end subroutine



subroutine SUPD_fib(Sfib,DSfib,qfib,damagePar,Param,matPar,constlaw,NDimE,damageMethod,isAnalytic) 
	real*8, intent(in) :: qfib(:), damagePar(:), matPar(:)
	real*8, intent(out) :: Sfib(:),DSfib(:,:) 
	real*8, intent(inout) :: Param(*)
	integer , intent(in) :: NdimE, constLaw , damageMethod, isAnalytic

!~ 	! this is to compute undamaged responses
	if(isAnalytic == 0) then	
		call calcSfib(Sfib,qfib,matPar,constLaw,NDimE) 
		call calcDSfib(DSfib,qfib,matPar,constLaw,NDimE)
	else if(isAnalytic == 1) then
		call calcSfibDSfib_analytic(Sfib,DSfib,qfib,matPar,constLaw)
	else
		write(0,*) 'SUPD_fib, option not found for isAnalytic ' , isAnalytic 
	end if

	if(damageMethod>0) then ! it allows damage
		call stateVariableUpdate_fib(qfib,damagePar,Param,matPar, constlaw)
		call modifyStressAndTangent_fib(Sfib,DSfib, damagePar,Param,NDimE,damageMethod)
	end if
	
end subroutine


subroutine SUPD_fib_1D(Sfib,DSfib,qfib,damagePar,Param,matPar,constlaw,NDimE,damageMethod) 
	real*8, intent(in) :: qfib(:), damagePar(:), matPar(:)
	real*8, intent(out) :: Sfib(:),DSfib(:,:) 
	real*8, intent(inout) :: Param(*)
	integer , intent(in) :: NdimE, constLaw , damageMethod

	real*8 :: dsfib_1D, sfib_1D, l , dPsi, d2Psi, qfib_ten_qfib(NdimE,NdimE)

	call derivatives_strainEnergy_fib(dPsi,d2Psi,qfib,matPar,constLaw) ! by respect of I4
	
	l = norm2(qfib)

	dsfib_1D = 2.0d0*dPsi + 4.0d0*l*l*d2Psi ! second derivative of undamaged energy
	sfib_1D = 2.0d0*l*dPsi ! first derivative of undamaged energy

	if(damageMethod>0) then ! it allows damage
		call stateVariableUpdate_fib(qfib,damagePar,Param,matPar, constlaw)
		call modifyStressAndTangent_fib_1D(sfib_1D,dsfib_1D, sfib_1D, damagePar,Param, damageMethod)
	end if
	
	call VotimesW(qfib_ten_qfib,qfib,qfib)
	
	Sfib = (sfib_1D/l)*qfib
	DSfib = ((dsfib_1D - sfib_1D/l)/(l*l))*qfib_ten_qfib + (sfib_1D/l)*deltaKron(1:NdimE,1:NdimE)
	
end subroutine

subroutine SUPD_fib_1D_visc(Sfib,DSfib,qfib,damagePar,Param,matPar,Dt,constlaw,NDimE,damageMethod) 
	real*8, intent(in) :: qfib(:), damagePar(:), matPar(:), Dt
	real*8, intent(out) :: Sfib(:),DSfib(:,:) 
	real*8, intent(inout) :: Param(*)
	integer , intent(in) :: NdimE, constLaw , damageMethod

	real*8 :: dsfib_1D, sfib_1D, sfib_1D_e, l , dPsi, d2Psi, qfib_ten_qfib(NdimE,NdimE), eta, lold

	call derivatives_strainEnergy_fib(dPsi,d2Psi,qfib,matPar,constLaw) ! by respect of I4
	
	l = norm2(qfib)

	dsfib_1D = 2.0d0*dPsi + 4.0d0*l*l*d2Psi ! second derivative of undamaged energy
	sfib_1D = 2.0d0*l*dPsi ! first derivative of undamaged energy

	eta = matPar(4)
	lold = matPar(5)
	
! viscosity strategy 3
	sfib_1D_e = sfib_1D
	sfib_1D = sfib_1D_e + (eta/Dt)*(l-lOld)
	dsfib_1D = dsfib_1D + eta/Dt


! viscosity strategy 4
!~ 	sfib_1D_e = sfib_1D
!~ 	sfib_1D = sfib_1D + (eta/Dt)*(l-lOld)/l
!~ 	dsfib_1D = dsfib_1D + (eta/Dt)*lOld/(l*l)
	
!~ 	sfib_1D = sfib_1D_e + (eta/Dt)*(l-lOld)/(lOld)
!~ 	dsfib_1D = dsfib_1D + (eta/Dt)/lOld

	if(damageMethod>0) then ! it allows damage
		call stateVariableUpdate_fib(qfib,damagePar,Param,matPar, constlaw)
		call modifyStressAndTangent_fib_1D(sfib_1D,dsfib_1D, sfib_1D_e, damagePar,Param, damageMethod)
	end if
	
	call VotimesW(qfib_ten_qfib,qfib,qfib)
	
!~ ! viscosity strategy 1	
!~ 	sfib_1D = sfib_1D + eta*(l-lOld)/Dt
!~ 	dsfib_1D = dsfib_1D + eta/Dt

! viscosity strategy 2
!~ 	sfib_1D = sfib_1D + (eta/Dt)*(l-lOld)/l
!~ 	dsfib_1D = dsfib_1D + (eta/Dt)*lOld/l
	
	Sfib = (sfib_1D/l)*qfib
	DSfib = ((dsfib_1D - sfib_1D/l)/(l*l))*qfib_ten_qfib + (sfib_1D/l)*deltaKron(1:NdimE,1:NdimE)
	
	
end subroutine



end module
