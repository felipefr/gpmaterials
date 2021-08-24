!> Insert deformation in the globalvariable Module \n
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
!!	@param iShiftMinDetQ = nint(commonPar(1))  ==> a (<0 if not aplicable)
!!	@param evolStyle = nint(commonPar(2))   ==> b ( 0 - absolute, 1 - incremental by absolute, 2- pure incremental)
!!	@param nbFtype = nint(commonPar(3))   ==> c !! related to the number of regions
!!	@param nbLoadType = nint(commonPar(4))   ==> d !! related to before and after bifurcation
!!	@param LoadType = nint(CommonPar(see structure)) ==> e1 
!!	@param LoadTypeProg = nint(CommonPar(see structure)) ==> e2 (1 - for linear, 2 for constant, for more see code)
!!	@param LoadPar = nint(CommonPar(see structure)) ==> e3:8
!!  @param mappingRegionLoadType ==> for d regions and c load cases
!! @author Rocha, Felipe Figueredo

Subroutine insertDeformationGeneral &
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
	
	integer, parameter :: nLoadParam = 8, iShiftLoadPar = 4
    integer ::  i ,j, ip , nbFtype,  nbLoadType, iShiftMinDetQ, evolStyle, LoadProg,LoadType, iLoadType
    Real*8, allocatable :: LoadParGen(:,:) 
    Real* 8 :: LoadPar(6)
    logical , save:: afterBif = .false.
    integer , parameter :: OUnitStrecht = 56
    real*8 :: Factual(NdimE,NdimE) , lambda
    
	iShiftMinDetQ = nint(commonPar(1))
	evolStyle = nint(commonPar(2))
	nbFtype = nint(commonPar(3)) !! related to the number of regions
	nbLoadType = nint(commonPar(4)) !! related to before and after bifurcation

	allocate(LoadParGen(nLoadParam,nbLoadType))
	
	if((.not. afterBif) .and. iShiftMinDetQ>-1 ) then
		if(Sol1( iShiftMinDetQ + 1) < 0.0d0) then
			afterBif = .true.
		end if				
	end if
	
	
!~ 	if(.not. afterBif) then
!~ 		if(Sol1( iShiftMinDetQ + 1) < 0.0d0) then
!~ 			afterBif = .true.
!~ 			
!~ 			call useFOld()
!~ 			
!~ 			open (OUnitStrecht, file='LambdaStrecht.txt', Access = 'append') 
!~ 			do i = 1, nbFtype
!~ 				call getF(Factual,i) 
!~ 				lambda = dsqrt(dot_product2(Factual-deltaKron(1:NdimE,1:NdimE),Factual-deltaKron(1:NdimE,1:NdimE)))
!~ 				write(OUnitStrecht,*) lambda
!~ 			end do
!~ 			close(OUnitStrecht)
!~ 			
!~ 			 	
!~ 			AE(1,1) = 1.0d0
!~ 			BE(1) = Sol1(1)
!~ 			
!~ 			return
!~ 		end if
!~ 	end if

	if(.not. isAllocated_Fsave() ) then
		call allocateFsave(NdimE,nbFtype)
	end if
	
	do i = 1 , nbLoadType
		ip = (i-1)*nLoadParam + iShiftLoadPar
		LoadParGen(:,i) = commonPar(ip + 1 : ip + nLoadParam )
	end do

	ip = nbLoadType*nLoadParam + iShiftLoadPar
	
	if(afterBif) ip = ip + nbFType 

	do i = 1 , nbFtype
		iLoadType = nint(commonPar(ip + i))
		if(iLoadType == 0) cycle
		LoadType = nint(LoadParGen(1,iLoadType))
		LoadProg = nint(LoadParGen(2,iLoadType))
		LoadPar = LoadParGen(3:8,iLoadType)
		call evolutionF(i,LoadPar,Time,delT,LoadProg,LoadType,NdimE,evolStyle)
	end do
	
	open (OUnitStrecht, file='LambdaStrecht.txt', Access = 'append') 
	do i = 1, nbFtype
		call getF(Factual,i) 
		lambda = dsqrt(dot_product2(Factual-deltaKron(1:NdimE,1:NdimE),Factual-deltaKron(1:NdimE,1:NdimE)))
		write(OUnitStrecht,*) lambda
	end do
	close(OUnitStrecht)

	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)
	
	
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine insertDeformationGeneralS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)

	Coupling(1,1) = 1

end Subroutine
