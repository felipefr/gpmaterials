    !> Dirichlet conditions via Lagrange Multiplier
    !!
	!! @param dofD = nint(CommonPar(1))
	!! @param dofLag = nint(CommonPar(2))
	!! @param pen = nint(CommonPar(3)) 
	!! @param eps = nint(commonPar(4))
	!! @param LoadTypeProg = nint(commonPar(5))
	!! @param LoadPar(1:6) = nint(commonPar(6:11))
	!!
    !! @author Rocha, Felipe Figueredo

Subroutine DirichletNodeLag &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)

	use funcAux
	use loadingLib
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
   
    Real*8 :: Vd , pen, LoadPar(6), eps
    integer :: dofD, LoadTypeProg, dofLag , dofDold

	dofD = nint(commonPar(1))
	dofLag = nint(commonPar(2))
	dofDold = nint(commonPar(3)) !!! >0 in case of incremental
	pen = commonPar(4)
	eps = commonPar(5)
	LoadTypeProg = nint(commonPar(6))
	LoadPar(:) = CommonPar(7:12)
	
	call chooseLoad(LoadTypeProg)
	call pressureLoad(Vd,time,DelT,loadPar)
	
	if(dofDold>0) then
		BE(dofLag) = pen*(Vd - Sol1(dofDold)) + eps*Sol1(dofLag) !!! this is for incremental approach
	else
		BE(dofLag) = pen*Vd + eps*Sol1(dofLag) !!! this is for non incremental approach
	end if   	
	AE(dofD,dofLag) = pen
	AE(dofLag,dofD) = pen
	AE(dofLag,dofLag) = eps
		
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine DirichletNodeLagS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	use funcAux
    implicit none
    
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    integer dofD, dofLag
    
    dofD = nint(commonPar(1))
    dofLag = nint(commonPar(2))
	
	Coupling(dofLag,dofD) = 1
	Coupling(dofD,dofLag) = 1
	Coupling(dofLag,dofLag) = 1
	
	call numPrint(Coupling)
	
end Subroutine
