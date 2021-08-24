!> Insert information of constitutive materials in the globalvariable Module
!!  @param CommonPar = (a,b,b*[c,a*di])  
!!	@param maxMatArg  = nint(commonPar(1))  ==> a
!!	@param nbMaterials = nint(commonPar(2))   ==> b
!!	@param constLaw = nint(commonPar(see structure)) ==> c (for each material) 
!!	@param MatPar = CommonPar(see structure) ==> di (for each material)
!! @author Rocha, Felipe Figueredo

Subroutine setMaterialParam &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
	
	use globalVariables
	use funcAux
    
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
!~     integer , save :: ElmNumber = 0
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
	
	integer, parameter :: iShiftMaterial = 2
    integer ::  i ,j, ip, nbMaterials, maxMatArg, constLaw
    real*8 , allocatable :: matpar(:) 
   
	maxMatArg = nint(commonPar(1))
	nbMaterials = nint(commonPar(2))

	if(.not. isAllocatedMaterial() ) then
		call allocateMaterial(nbMaterials)
		
		allocate(matPar(maxMatArg))
		do i = 1 , nbMaterials
			ip = (i-1)*(maxMatArg + 1) + iShiftMaterial
			constLaw = nint(commonPar(ip+1))
			matpar = commonPar( ip + 2 : ip + maxMatArg + 1)
			call setMaterial(constLaw,matpar,i)
		end do
		
		
		deallocate(matPar)
	end if
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)

end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine setMaterialParamS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)

	Coupling(1,1) = 1

end Subroutine
