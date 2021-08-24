!> Insert deformation in the globalvariable Module (very simple) \n
!! example i) \n
!! 	16 2 2 2 \n
!! 	14 2 0.00005 0.0 0.0 0.0 0.5 0.5 \n 
!! 	14 2 0.018 0.0 0.0 0.0 0.5 0.5 \n
!! 	1 1 \n
!! 	0 2 \n
!! example ii) \n
!! 	0 2 1 1 \n
!! 	14 2 0.00005 0.0 0.0 0.0 0.5 0.5 \n 
!! 	1 \n
!!  @param CommonPar = (a,b,c,d, d*[nLoadParam*ei] , c*[d*fi] ) , actually nLoadParam = 8 
!!	@param iShiftMinDetQ = nint(commonPar(1))  ==> a
!!	@param evolStyle = nint(commonPar(2))   ==> b
!!	@param nbFtype = nint(commonPar(3))   ==> c !! related to the number of regions
!!	@param nbLoadType = nint(commonPar(4))   ==> d !! related to before and after bifurcation
!!	@param LoadType = nint(CommonPar(see structure)) ==> e1
!!	@param LoadTypeProg = nint(CommonPar(see structure)) ==> e2
!!	@param LoadPar = nint(CommonPar(see structure)) ==> e3:8
!!  @param mappingRegionLoadType ==> for d regions and c load cases
!! @author Rocha, Felipe Figueredo

Subroutine insertDeformation &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)

	use globalVariables
	use funcAux
	use DETERMINANT , only : iElem,IdElemLib
    
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
!~     integer , save :: ElmNumber = 0
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
	
	integer, parameter :: nLoadParam = 8
    integer ::  i ,j, ip , nbFtype,  nbLoadType, iShiftMinDetQ, evolStyle, LoadProg,LoadType, iLoadType
    Real*8 :: LoadParGen(nLoadParam) 
    Real* 8 :: LoadPar(nLoadParam-2)
    logical , save:: afterBif = .false.
    integer , parameter :: OUnitStrecht = 56
    real*8 :: Factual(NdimE,NdimE) , lambda
    
	evolStyle = nint(commonPar(1))
	LoadParGen = commonPar(2:nLoadParam + 1)
	
	if(.not. isAllocated_Fsave() ) then
		call allocateFsave(NdimE,1)
	end if
		
	LoadType = nint(LoadParGen(1))
	LoadProg = nint(LoadParGen(2))
	LoadPar = LoadParGen(3:8)

	call evolutionF(1,LoadPar,Time,delT,LoadProg,LoadType,NdimE,evolStyle)
	
	open (OUnitStrecht, file='LambdaStrecht.txt', Access = 'append') 
	
	call getF(Factual,1) 
	lambda = dsqrt(dot_product2(Factual-deltaKron(1:NdimE,1:NdimE),Factual-deltaKron(1:NdimE,1:NdimE)))
	write(OUnitStrecht,*) lambda
	
	close(OUnitStrecht)

	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)
	
	
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine insertDeformationS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)

	Coupling(1,1) = 1

end Subroutine
