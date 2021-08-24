module finiteStrainLib

use funcAux
use HyperModelsLib

implicit none

public calcGradU , calcU, DGradPhiBdotGradPhiA, DGradGrad, calcPpk, calcD, modifyToSpatial, calcFC,&
		buildFbarTensors, build_tenQs, getPKmuFromMatrix, calcPpkD_analytic, calcForCiarlet,  calcSDAnalytical

contains

subroutine calcGradU(GradU,Sol,dPhi_G,NdimE,NodElT,iDofSol)
	integer, intent(in) :: NdimE, NodElT, iDofSol 
	Real*8, intent(out) :: GradU(NdimE,NdimE)
	Real*8, intent(in) ::dPhi_G(NdimE,NodElT) , Sol(iDofSol*NodElt)
	integer :: i, j , e ,ep 
	
	GradU = 0.0d0
	do e = 1 , NodElT
		ep = (e-1)*iDofSol
		do i = 1, NdimE
			do j = 1, NdimE
				GradU(i,j) = GradU(i,j) + Sol(ep+i)*dPhi_G(j,e)
			end do
		end do
	end do	
end subroutine

subroutine calcU(U,Sol,Phi,NdimE,NodElT,iDofSol)
	integer, intent(in) :: NdimE,NodElT , iDofSol
	Real*8, intent(out) :: U(NdimE)
	Real*8, intent(in) :: Phi(NodElT) , Sol(iDofSol*NodElt)
	integer :: i, e ,ep 
	
	U= 0.0d0
	do e = 1 , NodElT
		ep = (e-1)*iDofSol
		do i = 1, NdimE
			U(i) = U(i) + Sol(ep+i)*Phi(e)
		end do
	end do	
	
end subroutine

real*8 function DGradGrad(DD,dPhiA,dPhiB,pA,pB,NdimE) result(rAux) !! correct order
	real*8 , intent(in) :: DD(NdimE,NdimE,NdimE,NdimE),dPhiB(NdimE),dPhiA(NdimE)
	integer, intent(in) :: pA,pB, NdimE
	integer :: j,l
	
	rAux = 0.0d0
	
	do j=1,NdimE
	do l=1,NdimE
		rAux= rAux + DD(pB,j,pA,l)*dPhiB(j)*dPhiA(l) 
	end do
	end do
	
end function

real*8 function DGradPhiBdotGradPhiA(DD,dPhiB,dPhiA,p,q,NdimE) result(rAux) !! It would be more easy if it would q,p
	real*8 , intent(in) :: DD(NdimE,NdimE,NdimE,NdimE),dPhiB(NdimE),dPhiA(NdimE)
	integer, intent(in) :: p,q, NdimE
	integer :: j,l
	
	rAux = 0.0d0
	
	do j=1,NdimE
	do l=1,NdimE
		rAux= rAux + DD(p,j,q,l)*dPhiA(j)*dPhiB(l) 
	end do
	end do
	
end function
 
subroutine calcPpk(Ppk,F,NdimE,matPar,constLaw)
	integer , intent(in) :: NDimE,constLaw
	real*8, intent(in) :: MatPar(:), F(NDimE,NDimE)
	real*8, intent(out) :: Ppk(NDimE,NDimE)
	real*8 , parameter :: eps = 1.d-4
	real*8 ::  Fp(NDimE,NDimE), inv2eps, energyPr , energyPl
	integer :: i,j 

	inv2eps = 0.5d0/eps
	
	Do i=1,NdimE  
	Do j=1,NdimE  
		
		Fp=F
		Fp(i,j) = Fp(i,j) + eps
		Call StrainEnergy(energyPr,Fp,MatPar,constLaw)

		Fp=F
		Fp(i,j) = Fp(i,j) - eps
		Call StrainEnergy(energyPl,Fp,MatPar,constLaw)

		Ppk(i,j) = (energyPr - energyPl)*inv2eps
	Enddo
	Enddo
	
end subroutine

subroutine calcD(D,F,NdimE,matPar,constLaw)
	integer , intent(in) :: NDimE,constLaw
	real*8, intent(in) :: MatPar(:), F(NDimE,NDimE)
	real*8, intent(out) :: D(NDimE,NDimE,NDimE,NDimE)
	real*8 , parameter :: eps = 1.d-4
	real*8 ::  Fp(NDimE,NDimE), inveps2, inv4eps2, energy, energyPrr , energyPrl, energyPlr, energyPll
	integer :: i,j, k, l

	inv4eps2 = 0.25d0/(eps*eps)
	inveps2 = 1.0d0/(eps*eps)
	
	call strainEnergy(energy,F,MatPar,constLaw)
	
	Do i=1,NdimE  
	Do j=1,NdimE  
	Do k=1,NdimE  
	Do l=1,NdimE  
		
		if(i==k .and. j==l) then
			Fp=F
			Fp(i,j) = Fp(i,j) + eps
			Call StrainEnergy(energyPrr,Fp,MatPar,constLaw)
				
			Fp=F
			Fp(i,j) = Fp(i,j) - eps
			Call StrainEnergy(energyPll,Fp,MatPar,constLaw)
			
			D(i,j,i,j) = (energyPrr + energyPll - 2.0d0*energy)*inveps2
		else
			Fp=F
			Fp(i,j) = Fp(i,j) + eps
			Fp(k,l) = Fp(k,l) + eps
			Call StrainEnergy(energyPrr,Fp,MatPar,constLaw)
			
			Fp=F
			Fp(i,j) = Fp(i,j) + eps
			Fp(k,l) = Fp(k,l) - eps
			Call StrainEnergy(energyPrl,Fp,MatPar,constLaw)

			Fp=F
			Fp(i,j) = Fp(i,j) - eps
			Fp(k,l) = Fp(k,l) + eps
			Call StrainEnergy(energyPlr,Fp,MatPar,constLaw)
		
			Fp=F
			Fp(i,j) = Fp(i,j) - eps
			Fp(k,l) = Fp(k,l) - eps
			Call StrainEnergy(energyPll,Fp,MatPar,constLaw)
			
			D(i,j,k,l) = (energyPrr + energyPll - energyPrl - energyPlr)*inv4eps2
		end if
	end do
	end do
	end do
	end do
	
end subroutine

subroutine modifyToSpatial(Ds,Sigma,F,detF,NdimE)
	Real*8 , intent(out) :: sigma(NdimE,NdimE), Ds(NdimE,NdimE,NdimE,NdimE)
	Real*8 , intent(in) :: F(NdimE,NdimE),detF
	integer , intent(in) :: NDimE 
	
	Real*8 :: D(NdimE,NdimE,NdimE,NdimE), detFtemp
	integer :: i,j,k,l,m,n,p,q
	
	D = Ds
	
	sigma = (1.0d0/detF)*matmul(sigma,transpose(F))
	
	Ds = 0.0d0
	
	do i = 1,NdimE
	do j = 1,NdimE
	do k = 1,NdimE
	do l = 1,NdimE
	do p = 1,NdimE
	do q = 1,NdimE
		Ds(i,j,k,l) = Ds(i,j,k,l) +  D(i,p,k,q)*F(j,p)*F(l,q)/detF
	end do
	end do
	end do
	end do
	end do
	end do

end subroutine


subroutine calcFC(U,GradU,F,FT,FinvT,C,detF,Sol1,Phi,dPhi_G,NdimE,NodElT,iDofSol)
	integer, intent(in) :: NdimE,NodElT , iDofSol !!! Not exactly iDofT
	Real*8, intent(out) :: U(NdimE), GradU(NdimE,NdimE) , F(NdimE,NdimE), FT(NdimE,NdimE) , & 
								FinvT(NdimE,NdimE) , C(NdimE,NdimE) , detF
	Real*8, intent(in) ::Phi(NodElT), dPhi_G(NdimE,NodElT) , Sol1(iDofSol*NodElt)
	
	integer :: i, j , e ,ep 
	
	! Computes U and GradU
	U = 0.0d0
	GradU = 0.0d0
	do e = 1 , NodElT
		ep = (e-1)*iDofSol
		do i = 1, NdimE
			U(i) = U(i) + Sol1(ep+i)*Phi(e)
			do j = 1, NdimE
!~ 				if(i==3) write(*,*) Sol1(ep+i), dPhi_G(j,e)
				GradU(i,j) = GradU(i,j) + Sol1(ep+i)*dPhi_G(j,e)
			end do
		end do
	end do	
	
!~ 	call printVec(Sol1(1:9))
!~ 	call printMat(dPhi_G)
!~ 	! Computes F = I + GradU + F0 
	call addEye(F) !! F = F0 + I
	F = F + GradU !!! F might be initialized before as zeros or with a given deformation
	FT = transpose(F) !~ 	! Computes  FT = (F)^T
	C = matmul(FT,F) ! Cauchy-Grenn tensor
	call MatInv(FinvT,detF,FT)  ! FinvT= (FT)^-1 , defF = det(FT) = det(F)
	
	return
end subroutine

subroutine buildFbarTensors(Astar,PpkStar,Qstar,Q0star,F,F0,NdimE)
	
	real*8 , intent(inout) :: Astar(NdimE,NdimE,NdimE,NdimE), PpkStar(NdimE,NdimE)
	real*8 , intent(inout) :: Qstar(NdimE,NdimE,NdimE,NdimE), Q0star(NdimE,NdimE,NdimE,NdimE)
	real*8 , intent(in) :: F(NdimE,NdimE), F0(NdimE,NdimE)
	integer , intent(in) :: NdimE
	
	real*8 :: Finv(NdimE,NdimE), F0inv(NdimE,NdimE), Jratio, detF, detF0, alpha , beta, gamma, lamb
	integer :: i,j,k,l,p,q	

	call MatInv(Finv,detF,F) !! maybe detF0 is not correct because is 2D
	call MatInv(F0inv,detF0,F0) !! maybe detF0 is not correct because is 2D
	
	if(NdimE == 2) then
		alpha = -0.5d0
		beta = 0.0d0
		gamma = 0.5d0
		lamb = -0.5d0 
	else if(NdimE == 3) then
		alpha =  -2.0d0/3.0d0
		beta = -1.0d0/3.0d0
		gamma = 1.0d0/3.0d0
		lamb = -2.0d0/3.0d0
	end if
	

	Jratio = detF0/detF
	
	PpkStar = (Jratio**alpha)*PpkStar
	Astar = (Jratio**beta)*Astar
	
	Qstar = 0.0d0
	Q0star = 0.0d0
	
	do i = 1 , NdimE
	do j = 1 , NdimE
	do k = 1 , NdimE
	do l = 1 , NdimE
		do p = 1 , NdimE
		do q = 1 , NdimE
			Qstar(i,j,k,l) = Qstar(i,j,k,l) + gamma*Astar(i,j,p,q)*F(p,q)*Finv(l,k) !! for 2D
			Q0star(i,j,k,l) = Q0star(i,j,k,l) + gamma*Astar(i,j,p,q)*F(p,q)*F0inv(l,k) !! for 2D
		end do
		end do
		Qstar(i,j,k,l) = Qstar(i,j,k,l) + lamb*PpkStar(i,j)*Finv(l,k) 
		Q0star(i,j,k,l) = Q0star(i,j,k,l) + lamb*PpkStar(i,j)*F0inv(l,k) 
	end do
	end do
	end do
	end do
	
end subroutine

subroutine build_tenQs(Ds,sigma,tenQs,NdimE)	
	real*8 , intent(inout) :: Ds(NdimE,NdimE,NdimE,NdimE), sigma(NdimE,NdimE)
	real*8 , intent(inout) :: tenQs(NdimE,NdimE,NdimE,NdimE)
	integer , intent(in) :: NdimE
	
	integer :: i,j,k,l,p,q
	real*8 :: gamma, lamb
	
	tenQs = 0.0d0
	
	if(NdimE == 2) then
		gamma = 0.5d0
		lamb = -0.5d0 
	else if(NdimE == 3) then
		gamma = 1.0d0/3.0d0
		lamb = -2.0d0/3.0d0
	end if
	
	
	do i = 1 , NdimE
	do j = 1 , NdimE
	do k = 1 , NdimE
	do l = 1 , NdimE
		do p = 1 , NdimE
		do q = 1 , NdimE
			tenQs(i,j,k,l) = tenQs(i,j,k,l) + gamma*Ds(i,j,p,q)*deltaKron(p,q)*deltaKron(k,l)
		end do
		end do
		tenQs(i,j,k,l) = tenQs(i,j,k,l) + lamb*sigma(i,j)*deltaKron(k,l) 
	end do
	end do
	end do
	end do
	
end subroutine


subroutine getPKmuFromMatrix(PKmu,Vmu,SolUf,Xel,matpar,Param,constLaw,iFtype,NdimE,iFemType)
	use damageNewLib
	use ptsGaussLib
	use globalVariables , only : getF

	implicit none
		
	real*8 , intent(out) :: PKmu(:,:) , Vmu
	real*8, intent(in) :: matPar(:), Param(*), SolUf(:), Xel(:)
	integer, intent(in) :: NdimE, iFtype, constLaw, iFemType

	real*8 :: Fmu(NdimE,NdimE), dVng, PKaux(NdimE,NdimE) , GradUf(NdimE,NdimE)
	type(ptGaussClass) :: PtG
	integer :: NodG,NGP,pOrder,iBubble, iSimplex, nG
	
	PKmu = 0.0d0
	Vmu = 0.0d0
	
	call setFEMtype(iFemtype,NodG,pOrder,NGP,iSimplex,iBubble)
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
	
	do nG = 1 , NGP
		call PtG%calcGradU(GradUf,SolUf,nG)
		
		call getF(Fmu,iFtype)
		
		Fmu = Fmu + GradUf
		
		PKaux = 0.0d0
		call calcPpk(PKaux,Fmu,NdimE,matPar,constLaw)			
!~ 		call damageModifyStress(PKaux,Param,nG)
		
		dVng=ptG%dV(nG)
		
		Vmu = Vmu + dVng
		PKmu = PKmu + dVng*PKaux
	end do	
	
	PKmu = PKmu/Vmu

end subroutine


subroutine calcSDAnalytical(D,S,C,F,NDimE,MatPar,constLaw,Jp)
	! For Piola-Kirchoff (S) see Holzapfel Nonlinear Solid Mechanics , page 248
	! For Elasticity Tensor (D) see Holzapfel Nonlinear Solid Mechanics , page 261
	integer , parameter :: NF = 15 , ND = 6
	integer :: NDimE
	Real*8 , intent(out) :: S(NdimE,NdimE), D(NdimE,NdimE,NdimE,NdimE)
	Real*8 , intent(in) :: C(NdimE,NdimE) , Jp, MatPar(:) , F(NdimE,NdimE) 
	integer :: i,j,k,l,n
	Real*8 :: DPsi(12) , I1,I2,I3 , gamma(4), alpha(10), Dijkl , Cinv(NdimE,NdimE),  a0(3)
	integer , dimension(NF) :: mapI,mapJ,mapK,mapL
	integer , intent(in) :: constLaw
	
	DATA (mapI(I), I=1,NF) /1,1,1,1,1,1,1,1,1,1,1,1,2,2,3/  
	DATA (mapJ(I), I=1,NF) /1,1,1,1,1,2,2,2,2,3,3,3,2,2,3/  
	DATA (mapK(I), I=1,NF) /1,1,2,2,3,1,2,2,3,2,2,3,2,3,2/  
	DATA (mapL(I), I=1,NF) /2,3,2,3,3,3,2,3,3,2,3,3,3,3,3/ 
	
	DPsi = 0.0d0
	call StrainEnergyDPsi(DPsi,C,MatPar,constLaw)
!~ 	call numprint(DPsi)
	call calcI1(I1,C)
	call MatInv(Cinv,I3,C) ! It needs always compute I3 
	
	
!~ 	!a0 = MatPar(7:9) 

	gamma = 0.0d0
	alpha = 0.0d0
	
	gamma(1) = 2.0d0 * (DPsi(1) + I1*DPsi(2) )
	gamma(2) = -2.0d0 * DPsi(2)
	gamma(3) = 2.0d0 * I3 * DPsi(3)
!~ 	gamma(4) = 2.0d0 * DPsi(10) !! added for transverly isotropic
	
	
	alpha(1) = 4.0d0 * ( DPsi(4) + 2.0d0*I1*DPsi(7) + DPsi(2) + I1*I1*DPsi(5) )
	alpha(2) = -4.0d0 * ( DPsi(7) + I1*DPsi(5) )
	alpha(3) = 4.0d0 * ( I3*DPsi(9) + I1*I3*DPsi(8) )
	alpha(4) = 4.0d0 * DPsi(5)
	alpha(5) = -4.0d0 * I3 * DPsi(8)
	alpha(6) = 4.0d0 * ( I3*DPsi(3) + I3*I3*DPsi(6) )
	alpha(7) = -4.0d0 * I3 * DPsi(3)
	alpha(8) = -4.0d0 * DPsi(2)
!~ 	!alpha(9) = 4.0d0 * I3 * DPsi(12) !! added for transverly isotropic
!~ 	!alpha(10) = 4.0d0 * DPsi(11) !! added for transverly isotropic
	
	
!~ !	write(0,*) "Jp=" , Jp
	! volumetrical modification
	alpha(6) = alpha(6) + Jp
	alpha(7) = alpha(7) - 2.0d0*Jp
	gamma(3) = gamma(3) + Jp 
	
	if(NdimE  == 2) then ! I remember I had problems for deltaKron(1:NdimE,1:NdimE)
		S = gamma(1)*deltaKron2 + gamma(2)*C + gamma(3)*Cinv 
	else ! supposed NdimE = 3
		S = gamma(1)*deltaKron + gamma(2)*C + gamma(3)*Cinv 
	end if
	
	Do i=1,NdimE  
	Do j=i,NdimE  
		k=i
		l=j
		
!~ 	!	S(i,j) = S(i,j) + gamma(4)*a0(i)*a0(j)  ! added for transversely isotropic
!~ 	!	S(j,i) = S(i,j) 
		
		Dijkl=0.0d0
		Dijkl = Dijkl + alpha(1)*deltaKron(i,j)*deltaKron(k,l) 
		Dijkl = Dijkl + alpha(2)*( C(i,j)*deltaKron(k,l) + deltaKron(i,j)*C(k,l) ) 
		Dijkl = Dijkl + alpha(3)*( Cinv(i,j)*deltaKron(k,l) + deltaKron(i,j)*Cinv(k,l) ) 
		Dijkl = Dijkl + alpha(4)*C(i,j)*C(k,l)
		Dijkl = Dijkl + alpha(5)*( Cinv(i,j)*C(k,l) + C(i,j)*Cinv(k,l) )
		Dijkl = Dijkl + alpha(6)*Cinv(i,j)*Cinv(k,l)
		Dijkl = Dijkl + alpha(7)*0.5d0*( Cinv(i,k)*Cinv(j,l) + Cinv(i,l)*Cinv(j,k) ) 
		Dijkl = Dijkl + alpha(8)*0.5d0*( deltaKron(i,k)*deltaKron(j,l) + deltaKron(i,l)*deltaKron(j,k) ) 
		
!~ 	!	Dijkl = Dijkl + alpha(9)*( Cinv(i,j)*a0(k)*a0(l) + a0(i)*a0(j)*Cinv(k,l) ) ! added for transversely isotropic
!~ 	!	Dijkl = Dijkl + alpha(10)*a0(i)*a0(j)*a0(k)*a0(l) ! added for transversely isotropic
		
		D(i,j,k,l) = Dijkl
		D(i,j,l,k) = Dijkl
		D(j,i,k,l) = Dijkl
		D(j,i,l,k) = Dijkl
	end do
	end do
		
	do n = 1, NF
		i = mapI(n)
		j = mapJ(n)
		k = mapK(n)
		l = mapL(n)
					
		if(i>NdimE .or. j>NdimE .or. k>NdimE .or. l>NdimE) cycle ! protect against goes through allocated memory
		
		Dijkl=0.0d0
		Dijkl = Dijkl + alpha(1)*deltaKron(i,j)*deltaKron(k,l) 
		Dijkl = Dijkl + alpha(2)*( C(i,j)*deltaKron(k,l) + deltaKron(i,j)*C(k,l) ) 
		Dijkl = Dijkl + alpha(3)*( Cinv(i,j)*deltaKron(k,l) + deltaKron(i,j)*Cinv(k,l) ) 
		Dijkl = Dijkl + alpha(4)*C(i,j)*C(k,l)
		Dijkl = Dijkl + alpha(5)*( Cinv(i,j)*C(k,l) + C(i,j)*Cinv(k,l) )
		Dijkl = Dijkl + alpha(6)*Cinv(i,j)*Cinv(k,l)
		Dijkl = Dijkl + alpha(7)*0.5d0*( Cinv(i,k)*Cinv(j,l) + Cinv(i,l)*Cinv(j,k) ) 
		Dijkl = Dijkl + alpha(8)*0.5d0*( deltaKron(i,k)*deltaKron(j,l) + deltaKron(i,l)*deltaKron(j,k) ) 

!~ !		Dijkl = Dijkl + alpha(9)*( Cinv(i,j)*a0(k)*a0(l) + a0(i)*a0(j)*Cinv(k,l) ) ! added for transversely isotropic
!~ !		Dijkl = Dijkl + alpha(10)*a0(i)*a0(j)*a0(k)*a0(l) ! added for transversely isotropic		
		
		D(i,j,k,l) = Dijkl ! ok
		D(i,j,l,k) = Dijkl ! ok
		D(j,i,k,l) = Dijkl ! ok 
		D(j,i,l,k) = Dijkl ! ok
		D(k,l,i,j) = Dijkl ! ok
		D(k,l,j,i) = Dijkl ! ok
		D(l,k,i,j) = Dijkl ! ok 
		D(l,k,j,i) = Dijkl ! ok 
	enddo
	
	return
end subroutine


subroutine calcDtPpk(Dt,Ppk,Dc,Spk,F,NdimE)
	Real*8 , intent(out) :: Ppk(NdimE,NdimE), Dt(NdimE,NdimE,NdimE,NdimE) ! FS is temporaly S
	Real*8 , intent(in) :: F(NdimE,NdimE),  Spk(NdimE,NdimE), Dc(NdimE,NdimE,NdimE,NdimE)
	integer , intent(in) :: NDimE
	real*8 :: rAux
	integer :: i,j,k,l,m,n
	
	do i = 1,NdimE
	do j = 1,NdimE
	do k = 1,NdimE
	do l = 1,NdimE
		rAux = 0.0d0
		do m = 1,NdimE
		do n = 1,NdimE
			rAux = rAux + Dc(n,j,m,l)*F(k,m)*F(i,n)
		end do
		end do
		
		Dt(i,j,k,l) = rAux + Spk(l,j)*deltaKron(k,i)
	end do
	end do
	end do
	end do
	
	Ppk = matmul(F,Spk) ! before is just S

end subroutine

subroutine calcPpkD_analytic(Ppk,D,F,matPar,constLaw,NdimE)
	Real*8 , intent(out) :: Ppk(:,:), D(:,:,:,:) 
	Real*8 , intent(in) :: F(:,:), matpar(:)
	integer , intent(in) :: NDimE, constLaw
	real*8 :: rAux , Spk(NdimE,NdimE), Dc(NdimE,NdimE,NdimE,NdimE) , C(NdimE, NdimE), Jp 
	integer :: i,j,k,l,m,n
	
	Jp = 0.0d0 ! not a pressure-like formulation
	C = matmul(transpose(F),F)
	
	call calcSDAnalytical(Dc,Spk,C,F,NDimE,MatPar,constLaw,Jp)
	call calcDtPpk(D,Ppk,Dc,Spk,F,NdimE)
	
end subroutine


subroutine calcPpkD_analytic_Holzapfel(Ppk,D,F,matPar,constLaw,NdimE)
	Real*8 , intent(out) :: Ppk(:,:), D(:,:,:,:) 
	Real*8 , intent(in) :: F(:,:), matpar(:)
	integer , intent(in) :: NDimE, constLaw
	real*8 :: rAux , Spk(NdimE,NdimE), Dc(NdimE,NdimE,NdimE,NdimE) 
	integer :: i,j,k,l,m,n
	
	call calcSDAnalytical_Holzapfel(Dc,Spk,F,NDimE,MatPar,constLaw)
	call calcDtPpk(D,Ppk,Dc,Spk,F,NdimE)
	
end subroutine


subroutine calcForCiarlet(Ppk,D,F,matPar,NDimE) 
	real*8, intent(in) :: F(:,:),  matPar(:)
	real*8, intent(out) :: Ppk(:,:),D(:,:,:,:) 
	integer , intent(in) :: NdimE
	real*8 :: inv(3), E, nu , lambda, mu, Jac, cAux , FinvT(NdimE,NdimE)
	integer :: i,j,k,l
	
	E = MatPar(1) 
	nu = MatPar(2)
  
	lambda = E*nu/((1.0d0 + nu)*(1.0d0 - 2.0d0*nu))
	mu = 0.5d0*E/(1.0d0 + nu)
		
	Call MatInv(FinvT,Jac,F) ! still Finv
	FinvT = transpose(FinvT)
	
	cAux = 0.5d0*(Jac*Jac -1.0d0) - mu
	
	Ppk  = (cAux - mu)*FinvT + mu*F 
	
	do i = 1, NdimE
	do j = 1, NdimE
	do k = 1, NdimE
	do l = 1, NdimE
		D(i,j,k,l) = lambda * Jac * Jac * FinvT(i,j) * FinvT(k,l) - cAux * FinvT(i,l) * FinvT(k,j) + mu*deltaKron(i,k)*deltaKron(j,l)
	end do
	end do
	end do
	end do
	
end subroutine



subroutine calcSDAnalytical_Holzapfel(D,S,F,NDimE,MatPar,constLaw)
	! For Piola-Kirchoff (S) see Holzapfel Nonlinear Solid Mechanics , page 248
	! For Elasticity Tensor (D) see Holzapfel Nonlinear Solid Mechanics , page 261
	integer , parameter :: NS = 6
	integer :: NDimE
	Real*8 , intent(out) :: S(NdimE,NdimE), D(NdimE,NdimE,NdimE,NdimE)
	Real*8 , intent(in) :: MatPar(:) , F(NdimE,NdimE)
	integer :: II , JJ, i,j,k,l,p,q,m,n
	Real*8 :: DPsi_vol(12) , DPsi_iso(12), gamma(2), delta(4)
	integer , dimension(NS) :: mapI,mapJ
	integer , intent(in) :: constLaw
	Real*8  :: Siso(NdimE,NdimE), Diso(NdimE,NdimE,NdimE,NdimE)
	Real*8  :: C(NdimE,NdimE), Ciso(NdimE,NdimE) , Cinv(NdimE,NdimE)
	Real*8  :: tenP(NdimE,NdimE,NdimE,NdimE), tenPtilde(NdimE,NdimE,NdimE,NdimE)
	Real*8 :: UderJ , Uder2J, TrSiso , Jinv23, Jac, I1_iso, Dijkl
	
	DATA (mapI(I), I=1,NS) /1,2,3,2,1,1/  
	DATA (mapJ(I), I=1,NS) /1,2,3,3,2,3/  	
	
	C = matmul(transpose(F),F)
	call MatInv(Cinv,Jac,C) ! It needs always compute I3 
	Jac = dsqrt(Jac)
	Jinv23 = Jac**(-2.0d0/3.0d0)
	Ciso = Jinv23*C
	call calcI1(I1_iso,Ciso)
	
	DPsi_vol = 0.0d0
	call StrainEnergyDPsi(DPsi_iso,C,MatPar,0) ! constlaw == 0 means just volumetric contribution
	DPsi_iso = 0.0d0
	call StrainEnergyDPsi_iso(DPsi_iso,Ciso,MatPar,constLaw)
	
	gamma(1) = 2.0d0 * (DPsi_iso(1) + I1_iso*DPsi_iso(2) )
	gamma(2) = -2.0d0 * DPsi_iso(2)
	
	
	delta(1) = 4.0d0 * ( DPsi_iso(4) + 2.0d0*I1_iso*DPsi_iso(7) + DPsi_iso(2) + I1_iso*I1_iso*DPsi_iso(5) )
	delta(2) = -4.0d0 * ( DPsi_iso(7) + I1_iso*DPsi_iso(5) )
	delta(3) = 4.0d0 * DPsi_iso(5)
	delta(4) = -4.0d0 * DPsi_iso(2)
		
	UderJ = 2.0d0*Jac*DPsi_vol(3)
	Uder2J = 2.0d0*DPsi_vol(3) + 4.0d0*Jac*Jac*DPsi_vol(6)

	if(NdimE == 2) then
		Siso = gamma(1)*deltaKron2 + gamma(2)*Ciso
	else
		Siso = gamma(1)*deltaKron + gamma(2)*Ciso
	end if
	
	TrSiso = dot_product2(Siso,C)

	S = UderJ*Jac*Cinv !! volumetric part
	
	! build projection tensors and S total
	do i = 1,NdimE
	do j = 1,NdimE	
	do k = 1,NdimE
	do l = 1,NdimE
	
		tenP(i,j,k,l) = Jinv23*(deltaKron(i,k)*deltaKron(j,l) - (Cinv(i,j)*C(k,l))/3.0d0)
		tenPtilde(i,j,k,l) = 0.5d0*(Cinv(i,k)*Cinv(j,l) + Cinv(i,l)*Cinv(j,k)) + (Cinv(i,j)*C(k,l))/3.0d0
		
		Diso(i,j,k,l) = delta(1)*deltaKron(i,j)*deltaKron(k,l) + delta(2)*( deltaKron(i,j)*Ciso(k,l) + Ciso(i,j)*deltaKron(k,l) ) +&
						delta(3)*Ciso(i,j)*Ciso(k,l) + delta(4)*deltaKron(i,k)*deltaKron(j,l)
		
		! isochoric part
		S(i,j) = S(i,j) + tenP(i,j,k,l)*Siso(j,k)				
	enddo
	enddo
	enddo
	enddo
		
	do II = 1 , NS
		i = mapI(II)
		j = mapJ(II)
		
		if(i>NdimE .or. j>NdimE) cycle ! protect against goes through allocated memory
		
		do JJ = II , NS 
			k = mapI(JJ)
			l = mapJ(JJ)
					
			if(k>NdimE .or. l>NdimE) cycle ! protect against goes through allocated memory
		
	
			Dijkl= 0.0d0
			
			! volumetric part
			Dijkl = Dijkl + (Uder2J*J + UderJ)*J*Cinv(i,j)*Cinv(k,l)
			Dijkl = Dijkl - UderJ*J*(Cinv(i,k)*Cinv(j,l) + Cinv(i,l)*Cinv(j,k))
	
	
			! isochoric part
				
			Dijkl = Dijkl + (2.0d0*Jinv23/3.0d0)*( TrSiso*tenPtilde(i,j,k,l) - Siso(i,j)*Cinv(k,l) - Cinv(i,j)*Siso(k,l)) 

			do m = 1,NdimE
			do n = 1,NdimE
			do p = 1,NdimE
			do q = 1,NdimE
				Dijkl = Dijkl + tenP(i,j,m,n)*tenP(k,l,p,q)*Diso(m,n,p,q) 
			enddo
			enddo
			enddo
			enddo

			D(i,j,k,l) = Dijkl ! ok
			D(i,j,l,k) = Dijkl ! ok
			D(j,i,k,l) = Dijkl ! ok 
			D(j,i,l,k) = Dijkl ! ok
			D(k,l,i,j) = Dijkl ! ok
			D(k,l,j,i) = Dijkl ! ok
			D(l,k,i,j) = Dijkl ! ok 
			D(l,k,j,i) = Dijkl ! ok
		enddo 
	enddo
	
	
end subroutine

end module

