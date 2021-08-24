    !> Generic element for finite strains
    !!
    !!
	!! @param iShift_Uf = nint(CommonPar(1))
	!! @param iShift_dUf = nint(CommonPar(2))
	!! @param iFemType  = nint(CommonPar(3)) 
	!! @param iFType  = nint(commonPar(4))
	!! @param iMaterial  = nint(commonPar(5))
	!! @param iDamageParType  = nint(commonPar(6))
	!!
    !! @author Rocha, Felipe Figueredo

Subroutine arcLength_simple(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux

	use TIME_STEP_ITER
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
         
    Real*8 :: Bu(2) , Auu(2,2), a(2), lambda, B, g, dg_x
	real*8 :: theta0, w, fext(2)
	integer :: iShift_a, iShift_da, iShift_LarcLength, computeStar_or_Bar, NodG
	
	computeStar_or_Bar = nint(commonPar(1))
	iShift_a = nint(CommonPar(2))
	iShift_da = nint(CommonPar(3))
	iShift_LarcLength = nint(CommonPar(4))
	theta0 = commonPar(5) 
	w = commonPar(6) 
	fext(1:2) = commonPar(7:8)
	
	NodG = 2
	
	AE =  0.0d0
	BE =  0.0d0
	
	Auu = 0.0d0
	Bu = 0.0d0
		
	a(1) = Sol1(iShift_a + 1)
	a(2) = Sol1(iDofT + iShift_a + 1)	
	lambda = Sol1(NodG*iDofT + iShift_LarcLength + 1)
	
	B = 1.0d0 - 2.0d0*a(1)*dsin(theta0) + a(1)**2.0d0
	g = B**(-0.5d0) - 1.0d0
	dg_x = (dsin(theta0) - a(1))*B**(-1.5d0) 	
	 
	Auu(1,1) = dg_x*(dsin(theta0) - a(1)) - g
	! Auu(1,2) = 0
	Auu(2,1) = -w
	Auu(2,2) = w		
		
	if(computeStar_or_Bar == 1) then !! Star
		if(Nconverg == 1) then
			Bu = 0.0d0
		else !! -R
			Bu(1) = - g*(dsin(theta0) - a(1)) + lambda*fext(1)
			Bu(2) = - w*(a(2) - a(1)) + lambda*fext(2)
		end if
	else if(computeStar_or_Bar == 2) then !! Bar
		Bu = fext
	end if
	
	call setBlockToVector(BE,Bu,NodG,iDofT, 1, iShift_da ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	call setBlockToMatrixSquare(AE,Auu, NodG,iDofT, 1, iShift_da ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine arcLength_simpleS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
	integer :: iShift_da, NodG

    Coupling = 0
    
    iShift_da = nint(commonPar(3))

	NodG = 2
    
	call setSimbolicBlockToMatrixSquare(Coupling,NodG,iDofT, 1, iShift_da )

end Subroutine
