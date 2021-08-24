Subroutine DirichletNode &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)

	use funcAux
	use loadingLib
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
   
    Real*8 :: Vd , pen, LoadPar(6)
    integer :: dofD, LoadTypeProg 

	dofD = nint(commonPar(1))
	pen = commonPar(2)
	LoadTypeProg = nint(commonPar(3))
	LoadPar(:) = CommonPar(4:9)
	
!~ 	write(*,*) LoadPar
	
	call chooseLoad(LoadTypeProg)
	call pressureLoad(Vd,time,DelT,loadPar)
		
	AE(dofD,dofD) = pen
	BE(dofD) = pen*Vd 
		
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine DirichletNodeS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	use funcAux
    implicit none
    
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    integer dofD
    
    dofD = nint(commonPar(1))
	
	Coupling(dofD,dofD) = iAdd
	
	call numPrint(Coupling)
	
end Subroutine
