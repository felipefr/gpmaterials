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

Subroutine FSgenMixedStab(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : getF, NdimE, getDamagePar, getMaterial, maxDamagePar, maxMatPar
	use finiteStrainLib
	use damageNewLib
	use ptsGaussLib
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
         
	Integer :: p,q,nG! counters 
	Integer :: constLaw ! if NparamDamage = 0 , no damage
	Real*8 :: dV , Det , MatPar(maxMatPar), damagePar(maxDamagePar)
	Real*8 :: GradU(NdimE,NdimE) , F(NdimE,NdimE)
	Real*8 :: D(NdimE,NdimE,NdimE,NdimE), Ppk(NdimE,NdimE), Pmat(NdimE*NdimE), Dmat(NdimE*NdimE,NdimE*NdimE)
	Real*8 :: BmatU(NdimE*NdimE,NodElt*NdimE), BmatU_T(NodElt*NdimE,NdimE*NdimE),  BmatP(NdimE,3), BmatP_T(3,NdimE), NvecP(3)
	Real*8 :: Auu(NodElt*NdimE,NodElt*NdimE), Aup(NodElt*NdimE,3) , Apu(3,NodElt*NdimE), App(3,3) , Bu(NodElt*NdimE), Bp(3)
	Real*8 :: FinvTvec(NdimE*NdimE), J,  Finv(NdimE,NdimE), FinvT(NdimE,NdimE), Cinv(NdimE,NdimE), GradP(NdimE), P_ng
	Real*8 :: d_vecE(NodElt*NdimE), CmatE(NodElt*NdimE,3), DmatE(3,3), Amat(NdimE*NdimE,NdimE*NdimE) , Emat(NdimE,NdimE*NdimE)
	real*8 , allocatable ::  SolU(:) , SolP(:), Xel(:)
	integer :: pOrderU , pOrderP, NGP, NodGU, NodGP, iFEMtypeU , iFEMtypeP, iFtype, iMaterial, iDamageParType, iSimplex, iBubble
	integer :: iShiftU, iShiftDeltaU, iShiftDeltaP, iShiftP
	type(ptGaussClass) :: PtGP , PtGU
	real*8 :: kappa, alpha, mu, he, delta_e
	
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iFEMtypeU =  nint(commonPar(3))
	iShiftP = nint(CommonPar(4))
	iShiftDeltaP = nint(CommonPar(5))
	iFEMtypeP =  nint(commonPar(6)) 
	iFtype = nint(commonPar(7)) 
	iMaterial = nint(commonPar(8))
	iDamageParType = nint(commonPar(9))
	alpha = commonPar(10)

	if(iDamageParType>0) call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw, matPar, iMaterial)
	
	call setFEMtype(iFEMtypeP,NodGP,pOrderP,NGP,iSimplex,iBubble)
	call setFEMtype(iFEMtypeU,NodGU,pOrderU,NGP,iSimplex,iBubble) !! we know that NGP_U is bigger than NGP_P 

	call getSliceAllocate(SolU,Sol1,1,NodGU,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(SolP,Sol1,1,NodGP,iShiftP + 1 ,iShiftP + 1, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodGU,1 ,NdimE, Ndim) ! we know that NGP_U is bigger than NGP_P 
	
	call PtGP%init(Xel,NodGP,NdimE,NGP,pOrderP,iBubble, iSimplex)
	call PtGU%init(Xel,NodGU,NdimE,NGP,pOrderU,iBubble, iSimplex)
	
	AE =  0.0d0
	BE =  0.0d0
	
	Bu = 0.0d0
	Bp = 0.0d0	
	Auu = 0.0d0
	Apu = 0.0d0
	Aup = 0.0d0
	App = 0.0d0
	
	kappa = matpar(3)
	mu = matpar(1)
	call calc_he_tri_max(he,Xel + SolU)

	delta_e = 0.5d0*alpha*he*he/mu	

!~ 	write(0,*) "========================== > " , delta_e
	
	d_vecE = 0.0d0
	NvecP = 0.0d0
	
	Do nG = 1, NGP ! LoopGauss
		
		!! B's and N matrices/vectors
		call PtGU%getBmat(BmatU,nG)
		BmatU_T = transpose(BmatU)
		call PtGP%getBmatScalar(BmatP,nG)
		BmatP_T = transpose(BmatP)
		call PtGP%getNvecScalar(NvecP,nG)
		
		! F, gradP, P ... 
		call PtGU%calcGradU(GradU,SolU,nG)
		call getF(F,iFtype) !!! with identity summed
		F = F + GradU
		GradP = matmul(BmatP,solP)
		P_ng = dot_product(NvecP,SolP)
	
		! Finv, Cinv, J, ... comp
		Call MatInv(Finv,J,F)
		FinvT = transpose(Finv)
		call voigtTen2toVec(FinvTvec,FinvT)
		Cinv = matmul(Finv,FinvT)
		
		! auxiliary matrices, d_vecE, CmatE, DmatE
		d_vecE = J*matmul(BmatU_T, FinvTvec)
!~ 		call VotimesW(CmatE,d_vecE,NvecP)		!!!! I dont know, but is not working
		do p = 1,NodElt*NdimE
		do q = 1,3
			CmatE(p,q) = d_vecE(p)*NvecP(q)
		end do
		end do
		
		call VotimesW(DmatE,NvecP,NvecP)

		
		call calcPpk(Ppk,F,NdimE,matPar,constLaw)		
		call calcD(D,F,NdimE,matPar,constLaw)	
		
		if(iDamageParType>0) then
!~ 			call damageUpdateEquilibrium(F,matPar,damagePar,Param,constLaw,nG) ! deprecated 		
!~ 			call damageModifyTangent(D,Ppk,damagePar,Param, nG, 1) ! deprecated 
		end if						
		
		!! getAmat, getEmat, Pmat
		call getAmat(Amat,D,Finv,FinvT,J,NdimE)
		call getEmat(Emat,Finv,FinvT,Cinv,GradP,J,NdimE)
		call voigtTen4toTen2(Amat,D)
		call voigtTen2toVec(Pmat,Ppk)
		
		dV=ptGU%dV(nG)
		
		Auu = Auu + matmul(BmatU_T,matmul(Amat,BmatU))*dV
		Aup = Aup + CmatE*dV
		Apu = Apu + transpose(CmatE)*dV
		App = App - DmatE*(dV/kappa)
		
		
		Bu = Bu - ( matmul(BmatU_T,Pmat)  + p_ng*d_vecE )*dV
		Bp = Bp - NvecP*((J-1.0d0 - 1.0/kappa ) * dV)
		
		
		!! added stabilization terms
		Apu = Apu - matmul(BmatP_T,matmul(Emat,BmatU))*(delta_e*dV)
		App = App - matmul(BmatP_T,matmul(Cinv,BmatP))*(delta_e*J*dV)
		Bp = Bp + matmul(matmul(BmatP_T,Cinv),GradP)*(delta_e*J*dV)
		

	EndDo !LoopGauss
	
	if(iShiftU == iShiftDeltaU) then 
		Bu = Bu + matmul(Auu,SolU) + matmul(Aup,SolP)
		Bp = Bp + matmul(Apu,SolU) + matmul(App,SolP)
	end if
	
	call setBlockToVector(BE,Bu,NodGU,iDofT, NdimE, iShiftDeltaU ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	call setBlockToVector(BE,Bp,NodGP,iDofT, 1, iShiftDeltaP ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A	

	call setBlockToMatrixSquare(AE,Auu, NodGU,iDofT, NdimE, iShiftDeltaU ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	call setBlockToMatrixSquare(AE,App, NodGP,iDofT, 1, iShiftDeltaP ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
!~ 
	call setBlockToMatrix(AE,Aup, NodGU,NodGP, iDofT,iDofT, NdimE,1, iShiftDeltaU, iShiftDeltaP) ! (A,B, NodGi,NodGj, iDofAi,iDofAj, iDofBi,iDofBj, iShiftAi, iShiftAj)
	call setBlockToMatrix(AE,Apu, NodGP, NodGU, iDofT, iDofT, 1, NdimE, iShiftDeltaP, iShiftDeltaU ) ! (A,B, NodGi,NodGj, iDofAi,iDofAj, iDofBi,iDofBj, iShiftAi, iShiftAj)
	
 	if(iShiftU /= iShiftDeltaU) then
		call setBlockDiagonalConstantToMatrix(AE,1.0d0,NodGU, iDofT, NdimE, iShiftU, iShiftU)
		call setBlockDiagonalConstantToMatrix(AE,-1.0d0,NodGU, iDofT, NdimE, iShiftU, iShiftDeltaU)
		call setBlockToVector(BE,SolU,NodGU,iDofT, NdimE, iShiftU ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
		call setBlockDiagonalConstantToMatrix(AE,1.0d0,NodGP, iDofT, 1, iShiftP, iShiftP)
		call setBlockDiagonalConstantToMatrix(AE,-1.0d0,NodGP, iDofT, 1, iShiftP, iShiftDeltaP)
		call setBlockToVector(BE,SolP,NodGP,iDofT, 1, iShiftP ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	end if

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine FSgenMixedStabS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    use ptsGaussLib
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: iShiftU , iShiftDeltaU, NodGU, NodP, NodGP, iShiftDeltaP, iShiftP, iFEMtypeU , iFEMtypeP
    integer , parameter :: NdimE = 2

    Coupling = 0
    
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iFEMtypeU =  nint(commonPar(3))
	iShiftP = nint(CommonPar(4))
	iShiftDeltaP = nint(CommonPar(5))
	iFEMtypeP =  nint(commonPar(6)) 
    
    call setNodG(iFemTypeU, NodGU)
    call setNodG(iFemtypeP, NodGP)
    
	call setSimbolicBlockToMatrixSquare(Coupling,NodGU,iDofT, NdimE, iShiftDeltaU ) 
	call setSimbolicBlockToMatrixSquare(Coupling,NodGP,iDofT, 1, iShiftDeltaP) 
!~ 
	call setSimbolicBlockToMatrix(Coupling, NodGU,NodGP, iDofT,iDofT, NdimE,1, iShiftDeltaU, iShiftDeltaP) 
	call setSimbolicBlockToMatrix(Coupling, NodGP, NodGU, iDofT, iDofT, 1, NdimE, iShiftDeltaP, iShiftDeltaU) 


	if(iShiftU /= iShiftDeltaU) then
		call setSimbolicBlockDiagonalToMatrix(Coupling,NodGU, iDofT, NdimE, iShiftU, iShiftU)
		call setSimbolicBlockDiagonalToMatrix(Coupling,NodGU, iDofT, NdimE, iShiftU, iShiftDeltaU)
		call setSimbolicBlockDiagonalToMatrix(Coupling, NodGP, iDofT, 1, iShiftP, iShiftP)
		call setSimbolicBlockDiagonalToMatrix(Coupling, NodGP, iDofT, 1, iShiftP, iShiftDeltaP)
	end if
	
end Subroutine


subroutine getAmat(Amat,D,Finv,FinvT,Jac,NdimE)
	use funcAux
	IMPLICIT NONE
	integer , parameter :: n2 = 4
	
	real*8, intent(out) :: Amat(n2,n2)
	real*8, intent(in) :: FinvT(NdimE,NdimE), Finv(NdimE,NdimE), Jac, D(NdimE,NdimE,NdimE,NdimE)
	integer , intent(in) :: NdimE
	
	integer :: voigtMapI(n2), voigtMapJ(n2)
	integer :: i,j,k,l,p,q

	voigtMapI = (/1,1,2,2/)
	voigtMapJ = (/1,2,1,2/)
	
	
	do p = 1 , n2
	do q = 1, n2
		i = voigtMapI(p)
		j = voigtMapJ(p)
		k = voigtMapI(q)
		l = voigtMapJ(q)
		 
		Amat(p,q) = D(i,j,k,l) + Jac * FinvT(i,j)* FinvT(k,l) - Jac * FinvT(i,l) * Finv(j,k)
	
	end do
	end do
	
end subroutine

subroutine getEmat(Emat,Finv,FinvT,Cinv,GradP,Jac,NdimE)
	use funcAux
	IMPLICIT NONE
	integer , parameter :: n2 = 4
	
	real*8, intent(out) :: Emat(NdimE,n2)
	real*8, intent(in) :: FinvT(NdimE,NdimE), Finv(NdimE,NdimE), Cinv(NdimE,NdimE), GradP(NdimE), Jac
	integer , intent(in) :: NdimE
	
	integer :: voigtMapI(n2), voigtMapJ(n2)
	integer :: i,j,k,l,p,q

	voigtMapI = (/1,1,2,2/)
	voigtMapJ = (/1,2,1,2/)
	
	Emat = 0.0d0
	
	do p = 1 , NdimE
	do q = 1, n2
		i = p
		k = voigtMapI(q)
		l = voigtMapJ(q)
		
		do j = 1 , NdimE
		
			Emat(p,q) = Emat(p,q) + Jac * ( Cinv(i,j)*FinvT(k,l) - Finv(i,k)*Cinv(j,l) - Cinv(i,l)*Finv(j,k) )*GradP(j)
		
		end do 
	
	end do
	end do
	
end subroutine
