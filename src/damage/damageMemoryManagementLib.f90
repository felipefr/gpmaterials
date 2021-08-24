module damageMemoryManagementLib
use funcAux
use DETERMINANT , only : lengthParam 

integer :: iShiftDamageParam = 25 !! it has to set up before (just a initial value)
integer , parameter :: Idamage = 1 , Ird = 2 , Iqd = 3,  IrdDot = 4, IderDamage = 5, NdamageParam = 5

public writeDamage , writeDamageFromParam

type damageState
	real*8 :: damage, rd,qd , rdDot,derDamage
	
	contains
	procedure :: loadDamage, putDamage, putDamageIntoOld
end type

contains

subroutine loadDamage(s,Param, nG,inNew)
	class(damageState) :: s
	Real*8 , intent(in) :: Param(*)
	integer , intent(in) :: nG
	logical , intent(in) :: inNew
	integer :: iShiftDamageParamMod
	
	if(inNew) then
		iShiftDamageParamMod = LengthParam + (ng-1)*NdamageParam + iShiftDamageParam 
	else
		iShiftDamageParamMod = (ng-1)*NdamageParam + iShiftDamageParam
	end if
	
	s%damage = Param(iShiftDamageParamMod + Idamage)
	s%rd = Param(iShiftDamageParamMod + Ird)
	s%qd = Param(iShiftDamageParamMod + Iqd)
	s%rdDot = Param(iShiftDamageParamMod + IrdDot)
	s%derDamage = Param(iShiftDamageParamMod + IderDamage)
	
end subroutine

subroutine putDamage(s, Param, nG)
	class(damageState) :: s
	Real*8 , intent(out) :: Param(*)
	integer , intent(in) :: nG
	integer :: iShiftDamageParamMod
	
	iShiftDamageParamMod = LengthParam + (ng-1)*NdamageParam + iShiftDamageParam 
    
!~     write(0,*) "LengthParam = " , LengthParam 
     
	Param(iShiftDamageParamMod + Idamage) = s%damage
	Param(iShiftDamageParamMod + Ird) = s%rd
	Param(iShiftDamageParamMod + Iqd) = s%qd 
	Param(iShiftDamageParamMod + IrdDot) = s%rdDot 
	Param(iShiftDamageParamMod + IderDamage) = s%derDamage
	
end subroutine

subroutine putDamageIntoOld(s, Param, nG)
	class(damageState) :: s
	Real*8 , intent(out) :: Param(*)
	integer , intent(in) :: nG
	integer :: iShiftDamageParamMod
	
	iShiftDamageParamMod = (ng-1)*NdamageParam + iShiftDamageParam 
    
!~     write(0,*) "LengthParam = " , LengthParam 
     
	Param(iShiftDamageParamMod + Idamage) = s%damage
	Param(iShiftDamageParamMod + Ird) = s%rd
	Param(iShiftDamageParamMod + Iqd) = s%qd 
	Param(iShiftDamageParamMod + IrdDot) = s%rdDot 
	Param(iShiftDamageParamMod + IderDamage) = s%derDamage
	
end subroutine

subroutine writeDamageFromParam(Param,energy,nG)
	real*8 , intent(in) :: Param(*), energy
	integer , intent(in) :: nG
	type(damageState) :: sd,sdNew
	
	call sd%loadDamage(Param, nG, .false.)
	call sdNew%loadDamage(Param, nG, .true.)
	
	call writeDamage(sd%damage, sdNew%damage, energy)
	
end subroutine

subroutine writeDamage(damage,damageNew, energy)
	real*8 , intent(in) :: damage, damageNew, energy
	integer , parameter :: OUnitDamage = 15, OUnitDamageNew = 16, OUnitStrainEnergy = 17
	
	open (OUnitDamage, file='damage.txt', Access = 'append')
	open (OUnitDamageNew, file='damageNew.txt', Access = 'append') 
	open (OUnitStrainEnergy, file='strainEnergy.txt', Access = 'append') 
	
	write(OUnitDamage,*) damage
	write(OUnitDamageNew,*) damageNew
	write(OUnitStrainEnergy,*) energy
	
	close(OUnitDamage)
	close(OUnitDamageNew)
	close(OUnitStrainEnergy)
end subroutine

end module
