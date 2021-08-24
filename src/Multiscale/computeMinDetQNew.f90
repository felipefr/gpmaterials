Subroutine computeMinDetQNew &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)

	use funcAux
	use multiscaleNewLib
	
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
    integer :: nStepTheta, storageDetQ

	nStepTheta = nint(commonPar(1))
	storageDetQ = nint(commonPar(2))

	call updateBifurcationState(nStepTheta,storageDetQ) 
	call writeBifurcationState(storageDetQ)
	
!~ 	DhomOld = Dhom ! Old tensor removed from multiscaleNewLib, theoretically it's not necessary anymore
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)
		
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine computeMinDetQNewS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    implicit none
    
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer ::  Coupling(MaxLRows,MaxLRows)
    Real*8 :: CommonPar(*)
	
	Coupling(1,1) = 1
	
end Subroutine

