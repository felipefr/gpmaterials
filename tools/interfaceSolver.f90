subroutine setNdimE(v)
	use globalVariables, only: NdimE

	integer , intent(in) :: v
	
	NdimE = v

end subroutine

subroutine setLengthParam(v)
	use DETERMINANT, only: LengthParam
	
	integer , intent(in) :: v
	
	LengthParam = v
	
end subroutine

subroutine setiShiftDamageParam(v)
	use damageMemoryManagementLib, only: iShiftDamageParam
	integer , intent(in) :: v

	iShiftDamageParam = v

end subroutine 

subroutine initAllStructures(op,param,Nparam)
	use fibresMod, only : initFibresMod
	use MultiscaleLib, only : initMultiscaleLib
	use MultiscaleNewLib, only : initMultiscaleNewLib
	
	implicit none
	
	integer, intent(in) :: op , Nparam
	real*8 , intent(in) :: param(Nparam)
	
	integer NdimE 
	
	select case(op)
		
	case(1) ! case of fibres	
!~ 		write(0,*) 'Congratulations' , param
!~ 		pause
		NdimE = nint(param(1))
		call setNdimE(NdimE)
		
		if(NdimE == 2) then
			call setiShiftDamageParam(21)
			call setLengthParam(30)
		else if(NdimE == 3) then
			call setiShiftDamageParam(26)
			call setLengthParam(36)
		end if 
	
		call initFibresMod()
		call initMultiscaleLib()
		call initMultiscaleNewLib()

	case default
		
		write(0,*) 'Wrong option in initAllStructures ', op
	
	end select

end subroutine


subroutine initAllStructuresFileVersion()
	use funcAux , only : FindKeyword

	implicit none
	
	integer , parameter :: IUnit = 30
	
	real*8, allocatable :: param(:)
	integer :: op, nParam
	
	open(IUnit, file = 'Basparam.txt')
	
	call FindKeyword(IUnit,'*initAllStructures')
	
	read(IUnit,*) op ,  nParam

	allocate(param(nParam))
	
	read(IUnit,*) param
		
	call initAllStructures(op,param,nparam)

	close(IUnit)

end subroutine

