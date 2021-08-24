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

Subroutine computeCoeficientsSpherical_arcLength(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, & 
											Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : NdimE
	use ptsGaussLib
	use TIME_STEP_ITER
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
         
	Integer :: nG, i, ip, ipDim, ipDimE, ipDof,ipDofmod ! counters 
	Integer :: constLaw ! if NparamDamage = 0 , no damage
	Real*8 :: dV , Det
	real*8, allocatable, dimension(:) ::  Sol_dUf_Star, Sol_dUf_Bar , Sol_DUf_acc, Xel , SolGy, Sol_DU_acc
	integer :: pOrder , NGP, NodG, iFEMtype , iFtype, iSimplex, iBubble, lU 
	type(ptGaussClass) :: PtG
	 
	real*8 :: Gbar(NdimE,NdimE), arcLength, dLarcLength_acc
	integer :: 	iShift_c1, iShift_c2, iShift_c3, iShift_dUf_Star, iShift_dUf_Bar, & 
				iShift_DUf_acc, iShift_dLarcLength_acc
				
	real*8 :: Nmat(NdimE,NdimE*NodElt) , v(NdimE) , w(NdimE) , &
	          c1, c2, c3, psi, GbarNormSq, yG(NdimE), y(NdimE)
	
	integer :: kindArcLength
	
	integer :: it1, it2, it3
	real :: dampFactor
	
	arcLength = CommonPar(1)
	psi = CommonPar(2)
	iShift_c1 = nint(CommonPar(3))
	iShift_c2 = nint(CommonPar(4))
	iShift_c3 = nint(CommonPar(5))
	iShift_dUf_Star = nint(CommonPar(6))
	iShift_dUf_Bar = nint(CommonPar(7))
	iShift_DUf_acc = nint(CommonPar(8))
	iShift_dLarcLength_acc = nint(CommonPar(9))
	iFEMtype =  nint(commonPar(10)) 
	iFtype = nint(commonPar(11)) 
	kindArcLength = nint(commonPar(12))
	yG = commonPar(13:14)

	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)
	
	call getSliceAllocate(Sol_dUf_Bar,Sol1,1,NodG, iShift_dUf_Bar + 1 , iShift_dUf_Bar + NdimE, iDofT)
	call getSliceAllocate(Sol_dUf_Star,Sol1,1,NodG, iShift_dUf_Star + 1 , iShift_dUf_Star + NdimE, iDofT)
	call getSliceAllocate(Sol_DUf_acc,Sol1,1,NodG, iShift_DUf_acc + 1 , iShift_DUf_acc + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	allocate(SolGy(NodG*NdimE))
	allocate(Sol_DU_acc(NodG*NdimE))
	dLarcLength_acc = Sol1(NodG*idofT + iShift_dLarcLength_acc + 1)
	
	Gbar = 0.0d0
	Gbar(2,2) = 1.0d0
	
	GbarNormSq = dot_product2(Gbar,Gbar)
	
	lU = NodG*NdimE
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
	
	c1 = 0.0d0
	c2 = 0.0d0
	c3 = 0.0d0
	
	it1 = 5
	it2 = 20
	it3 = 30
	dampFactor = 0.5d0

	if(Nconverg<it1) then
		arcLength = arcLength*Nconverg/it1
	else if(Nconverg<it2) then
		arcLength = arcLength
	else if(Nconverg>it3) then
		arcLength = arcLength - (1.0d0 - dampFactor)*arcLength*real(Nconverg - it2)/real(it3 - it2)
	else
		arcLength = dampFactor*arcLength
	end if	
	
	do i = 1 , NodG
		ip = (i-1)*NdimE
		y = Xel(ip + 1 : ip + NdimE)
		SolGy(ip + 1 : ip + NdimE) = matmul(Gbar,y-yG)
	end do
	
	
	
	if(kindArcLength == 1) then
		c1 = dot_product(Sol_dUf_Bar,Sol_dUf_Bar) + GbarNormSq*psi**2.0 
		c2 = 2.0d0*(dot_product(Sol_dUf_acc + Sol_dUf_Star, Sol_dUf_Bar) + dLarcLength_acc*GbarNormSq*psi**2.0)
		c3 = dot_product(Sol_dUf_acc + Sol_dUf_Star,Sol_dUf_acc + Sol_dUf_Star) + &
							GbarNormSq*(dLarcLength_acc*psi)**2.0 - arcLength*arcLength 
	
	else if(kindArcLength == 2) then
		
		Sol_DU_acc = dLarcLength_acc*SolGy + Sol_DUf_acc
	
 		c1 = dot_product(SolGy + Sol_dUf_Bar,SolGy + Sol_dUf_Bar)
		c2 = 2.0d0*dot_product(SolGy + Sol_dUf_Bar, Sol_DU_acc + Sol_dUf_Star)
		c3 = dot_product(Sol_DU_acc + Sol_dUf_Star,Sol_DU_acc + Sol_dUf_Star) - arcLength*arcLength

	
	else if(kindArcLength == 3) then
	
		do nG = 1, NGP ! LoopGauss
			call PtG%getNmat(Nmat(:,1:lU),nG)
			
			Sol_DU_acc = dLarcLength_acc*SolGy + Sol_DUf_acc
			
			v = matmul(Nmat(:,1:lU), SolGy +  Sol_dUf_Bar) 
			w = matmul(Nmat(:,1:lU), Sol_DU_acc + Sol_dUf_Star ) 
			
			dV=ptG%dV(nG)
			
			c1 = c1 + dot_product(v,v)*dV
			c2 = c2 + 2.0d0*dot_product(v,w)*dV
			c3 = c3 + (dot_product(w,w) - arcLength*arcLength)*dV		
		end do !LoopGauss
	
	end if

	ip = NodG*iDofT + 1
	AE(ip + iShift_c1,ip + iShift_c1) =  1.0d0
	AE(ip + iShift_c2,ip + iShift_c2) =  1.0d0
	AE(ip + iShift_c3,ip + iShift_c3) =  1.0d0
	BE(ip + iShift_c1) =  c1
	BE(ip + iShift_c2) =  c2
	BE(ip + iShift_c3) =  c3
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine computeCoeficientsSpherical_arcLengthS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    use ptsGaussLib
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: iShift_c1, iShift_c2, iShift_c3, ip, NodG, iFEMtype

    Coupling = 0

    iShift_c1 = nint(CommonPar(3))
	iShift_c2 = nint(CommonPar(4))
	iShift_c3 = nint(CommonPar(5))
	iFEMtype =  nint(commonPar(10)) 
	
	call setNodG(iFEMtype, NodG)
	ip = NodG * idofT + 1
	
	Coupling(ip + iShift_c1,ip + iShift_c1) =  1
	Coupling(ip + iShift_c2,ip + iShift_c2) =  1
	Coupling(ip + iShift_c3,ip + iShift_c3) =  1

end Subroutine
