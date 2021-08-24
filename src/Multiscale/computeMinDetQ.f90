Subroutine computeMinDetQ &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)

	use funcAux
	use multiscaleLib, only : getMinDetQ
	
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
    integer, parameter :: NdimE = 2, nStepTheta = 999
    integer :: ip1 , ip2, iShiftC, iShiftMinDetQ, iShiftThetaCrit
    Real*8 ::  MinDetQ, ThetaCrit, betaCrit

	iShiftC = nint(commonPar(1))
	iShiftMinDetQ = nint(commonPar(2))
	iShiftThetaCrit = nint(commonPar(3))

	call getMinDetQ(minDetQ,ThetaCrit,betaCrit, Sol1,NdimE,idofT,1,iShiftC,nStepTheta) !! is '1' because is just one node
	
	ip1 = iShiftMinDetQ + 1
	ip2 = iShiftThetaCrit + 1
	
	AE(ip1,ip1) = 1.0d0
	AE(ip2,ip2) = 1.0d0
	AE(ip2 + 1,ip2 + 1) = 1.0d0
	
	BE(ip1) =  MinDetQ
	
	if(MinDetQ<0.0d0) then !! ThetaCrit is the immediately pre bifurcation angle
		BE(ip2) =  Sol1(ip2)
	else
		BE(ip2) =  ThetaCrit
		BE(ip2 + 1) =  betaCrit
	end if

		
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine computeMinDetQS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
	
    implicit none
    
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer ::  Coupling(MaxLRows,MaxLRows)
    Real*8 :: CommonPar(*)
    integer :: iShiftMinDetQ, iShiftThetaCrit, ip1, ip2
	
	iShiftMinDetQ = nint(commonPar(2))
	iShiftThetaCrit = nint(commonPar(3))

	ip1 = iShiftMinDetQ + 1	
	ip2 = iShiftThetaCrit + 1	
	
	Coupling(ip1,ip1) = 1
	Coupling(ip2,ip2) = 1
	
end Subroutine

