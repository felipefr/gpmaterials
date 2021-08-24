    !> Generic element for finite strains
    !!
    !!
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iShiftDeltaU = nint(CommonPar(2))
	!! @param iFemType  = nint(CommonPar(3)) 
	!! @param iFType  = nint(commonPar(4))
	!! @param iMaterial  = nint(commonPar(5))
	!! @param iDamageParType  = nint(commonPar(6))
	!!
    !! @author Rocha, Felipe Figueredo

Subroutine computeCoeficientsSpherical_arcLength_simple(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, & 
											Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
         
	integer :: 	ip, NodG, iShift_c1, iShift_c2, iShift_c3, iShift_da_Star, & 
				iShift_da_Bar, iShift_Da_acc, iShift_dLarcLength_acc
				
	real*8 :: c1,c2,c3, da_bar(2) , da_Star(2), Da_acc(2), & 
			arcLength, psi, fext(2), fextNormSq , dLarcLength_acc
	
	arcLength = CommonPar(1)
	Psi = CommonPar(2)
	fext = CommonPar(3:4)
	iShift_c1 = nint(CommonPar(5))
	iShift_c2 = nint(CommonPar(6))
	iShift_c3 = nint(CommonPar(7))
	iShift_da_Star = nint(CommonPar(8))
	iShift_da_Bar = nint(CommonPar(9))
	iShift_da_acc = nint(CommonPar(10))
	iShift_dLarcLength_acc = nint(CommonPar(11))
	
	
	NodG = 2
	
	da_bar(1) = Sol1(iShift_da_Bar + 1)
	da_bar(2) = Sol1(iDofT + iShift_da_Bar + 1)
	da_star(1) = Sol1(iShift_da_Star + 1)
	da_star(2) = Sol1(iDofT + iShift_da_Star + 1)
	Da_acc(1) = Sol1(iShift_Da_acc + 1)
	Da_acc(2) = Sol1(iDofT + iShift_Da_acc + 1)
	dLarcLength_acc = Sol1(iDofT*NodG + iShift_dLarcLength_acc + 1)
	
	fextNormSq = dot_product(fext,fext)
	
	c1 = dot_product(da_bar,da_bar) + fextNormSq*psi**2.0 
	c2 = 2.0d0*( dot_product(Da_acc + da_star, da_bar) + dLarcLength_acc*fextNormSq*psi**2.0)
	c3 = dot_product( Da_acc + da_star , Da_acc + da_star) + &
		fextNormSq*(dLarcLength_acc*psi)**2.0 - arcLength**2.0d0
	
	ip = NodG*iDofT + 1
	AE(ip + iShift_c1,ip + iShift_c1) =  1.0d0
	AE(ip + iShift_c2,ip + iShift_c2) =  1.0d0
	AE(ip + iShift_c3,ip + iShift_c3) =  1.0d0
	BE(ip + iShift_c1) =  c1
	BE(ip + iShift_c2) =  c2
	BE(ip + iShift_c3) =  c3
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine computeCoeficientsSpherical_arcLength_simpleS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: iShift_c1, iShift_c2, iShift_c3, NodG,  ip

    Coupling = 0

	iShift_c1 = nint(CommonPar(5))
	iShift_c2 = nint(CommonPar(6))
	iShift_c3 = nint(CommonPar(7))
	
	NodG = 2
	
	ip = NodG * idofT + 1
	
	Coupling(ip + iShift_c1,ip + iShift_c1) =  1
	Coupling(ip + iShift_c2,ip + iShift_c2) =  1
	Coupling(ip + iShift_c3,ip + iShift_c3) =  1

end Subroutine
