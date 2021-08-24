module ptsGaussLib
	use funcAux
	implicit none

	public GaussRuleQ, setFEMtype, setNodG
	
	!! A class to encapsulate issues related to shape function and its derivatives evaluated in Gauss Points. Provides as well the N, B matrices for scalar and vectors
	type :: PtGaussClass
	
		integer :: NGP, NodG, NdimE
		real*8 , allocatable :: dPhi_G(:,:,:), Phi(:,:) , dV(:)  ! allocate(s%dPhi_G(NdimE,NodG,NGP),s%Phi(NodG,NGP),s%dV(NGP))
	
		Contains
				
		procedure, public :: init, calcGradU, calcU, initCollapsed, initWithPsi, calcTheta, &
							calcGradTheta, getBmat, getBmatSym, getBmatScalar, getNvecScalar, getNmat
		
	end type
	
	contains 	
	
	!! Given numerical parameters and the element coordinates, stores the shape functions and its derivatives, in a gauss rule. Its works for both simplexes and quadrangles
	subroutine init(s,XLL,NodG,NdimE,NGP, pOrder,iBu, iSimplex)
		class(PtGaussClass) :: s
		integer , intent(in) :: NodG, NGP, NdimE, iBu, iSimplex, pOrder
		real*8 , intent(in) :: XLL(NdimE*NodG)
		integer :: i
		real*8 :: Det, W(NGP), Psi(NdimE,NGP), Jac(NdimE,NdimE), dPhi_L(NdimE,NodG)
			
		s%NodG = NodG
		s%NdimE = NdimE
		s%NGP = NGP 	
		
		allocate(s%dPhi_G(NdimE,NodG,NGP),s%Phi(NodG,NGP),s%dV(NGP))
			
		if(iSimplex == 0) then
			Call GaussRule (Psi,W,NdimE,NGP,iSimplex)
		else
			Call GaussRuleQ (Psi,W,NdimE,NGP,iSimplex)
		end if
		
		do i = 1 , NGP
			if(iSimplex == 0) then
				Call LocalShapeDer(NodG,NdimE,dPhi_L,PSI(1,i),pOrder,iBu)
				Call ShapeF (NodG,NdimE,s%Phi(1,i),PSI(1,i),pOrder,iBu)
			else 
				Call LocalShapeDerQ(NodG,NdimE,dPhi_L,PSI(1,i),pOrder,iBu)
				Call ShapeFQ (NodG,NdimE,s%Phi(1,i),PSI(1,i),pOrder,iBu)			 	
			end if
			
			Call Jacobian(Jac,Det,NdimE,NodG,XLL,dPhi_L)
			Call GlobalShapeDer(NodG,NdimE,dPhi_L,s%dPhi_G(1,1,i),Jac)
			s%dV(i) = W(i)*Det
			
		end do
		
	end subroutine
	
	!! Similar to init but the evaluation points and its weights (which may be gauss points) are given externally
	subroutine initWithPsi(s,Psi,W,XLL,NodG,NdimE,NGP, pOrder,iBu, iSimplex)
		class(PtGaussClass) :: s
		integer , intent(in) :: NodG, NGP, NdimE, iBu, iSimplex, pOrder
		real*8 , intent(in) :: XLL(NdimE*NodG), Psi(NdimE,NodG), W(NGP)
		integer :: i
		real*8 :: Det, Jac(NdimE,NdimE), dPhi_L(NdimE,NodG)
			
		s%NodG = NodG
		s%NdimE = NdimE
		s%NGP = NGP 	
		
		allocate(s%dPhi_G(NdimE,NodG,NGP),s%Phi(NodG,NGP),s%dV(NGP))
			
		do i = 1 , NGP
			if(iSimplex == 0) then
				Call LocalShapeDer(NodG,NdimE,dPhi_L,PSI(1,i),pOrder,iBu)
				Call ShapeF (NodG,NdimE,s%Phi(1,i),PSI(1,i),pOrder,iBu)
			else 
				Call LocalShapeDerQ(NodG,NdimE,dPhi_L,PSI(1,i),pOrder,iBu)
				Call ShapeFQ (NodG,NdimE,s%Phi(1,i),PSI(1,i),pOrder,iBu)			 	
			end if
			
			Call Jacobian(Jac,Det,NdimE,NodG,XLL,dPhi_L)
			Call GlobalShapeDer(NodG,NdimE,dPhi_L,s%dPhi_G(1,1,i),Jac)
			s%dV(i) = W(i)*Det
			
		end do
		
	end subroutine
	
	!! Similar to init but the gauss points are collapsed to one dimension below, used when we want gradients but the region is a surface, e.g. to evaluation the Neumann term due to spatial normal pressure in material configuration. 
	subroutine initCollapsed(s,XLL,NodG,NdimE,NGP, pOrder,iBu, iSimplex)
		class(PtGaussClass) :: s
		integer , intent(in) :: NodG, NGP, NdimE, iBu, iSimplex, pOrder
		real*8 , intent(in) :: XLL(NdimE*NodG)
		integer :: i, j
		real*8 :: Det, W(NGP), Psi(NdimE,NGP), Jac(NdimE,NdimE), dPhi_L(NdimE,NodG)
		real*8 :: DetCollap, PsiCollap(NdimE-1,NGP)
			
		s%NodG = NodG
		s%NdimE = NdimE
		s%NGP = NGP 	
		
		allocate(s%dPhi_G(NdimE,NodG,NGP),s%Phi(NodG,NGP),s%dV(NGP))
			
		if(iSimplex == 0) then
			Call GaussRule (PsiCollap,W,NdimE-1,NGP,iSimplex)			
		else
			Call GaussRuleQ (PsiCollap,W,NdimE-1,NGP,iSimplex)
		end if
		
		Psi(1:NdimE-1,:) = PsiCollap(1:NdimE-1,:)
		Psi(NdimE,:) = 0.0d0
		
		do i = 1 , NGP
			if(iSimplex == 0) then
				Call LocalShapeDer(NodG,NdimE,dPhi_L,PSI(1,i),pOrder,iBu)
				Call ShapeF (NodG,NdimE,s%Phi(1,i),PSI(1,i),pOrder,iBu)
			else 
				Call LocalShapeDerQ(NodG,NdimE,dPhi_L,PSI(1,i),pOrder,iBu)
				Call ShapeFQ (NodG,NdimE,s%Phi(1,i),PSI(1,i),pOrder,iBu)			 	
			end if
			
!~ 			write(0,*) dPhi_L
!~ 			pause
			
			Call Jacobian(Jac,Det,NdimE,NodG,XLL,dPhi_L)
			Call GlobalShapeDer(NodG,NdimE,dPhi_L,s%dPhi_G(1,1,i),Jac)
			s%dV(i) = W(i)*Det
			
		end do
		
	end subroutine
	
	!! encapsulates the GradU (U vector) calc, deprecated, you should use B matrix instead
	subroutine calcGradU(s,GradU,SolU,nG) 
		class(PtGaussClass) :: s
		integer, intent(in) :: nG 
		Real*8, intent(out) :: GradU(:,:) !! have NdimE x NdimE dimension
		Real*8, intent(in) :: SolU(:) !! have NdimE*NodG dimension
		integer :: i, j , e ,ep 
		
		GradU = 0.0d0
		do e = 1 , s%NodG
			ep = (e-1)*s%NdimE
			do i = 1, s%NdimE
				do j = 1, s%NdimE
					GradU(i,j) = GradU(i,j) + SolU(ep+i)*s%dPhi_G(j,e,nG)
				end do
			end do
		end do	
		
	end subroutine

	!! encapsulates the U vector calc, deprecated, you should use N matrix instead
	subroutine calcU(s,U,SolU,nG) 
		class(PtGaussClass) :: s
		integer, intent(in) :: nG 
		Real*8, intent(out) :: U(:) !! have NdimE dimension
		Real*8, intent(in) :: SolU(:) !! have NdimE*NodG dimension 
		integer :: i, e ,ep 
		
		U= 0.0d0
		do e = 1 , s%NodG
			ep = (e-1)*s%NdimE
			do i = 1, s%NdimE
				U(i) = U(i) + SolU(ep+i)*s%Phi(e,nG)
			end do
		end do	
		
	end subroutine
	
	!! encapsulates the GradTheta (Theta scalar) calc, deprecated, you should use B matrix instead
	subroutine calcGradTheta(s,GradTheta,SolTheta,nG) !! scalar field
		class(PtGaussClass) :: s
		integer, intent(in) :: nG 
		Real*8, intent(out) :: GradTheta(:) !! have NdimE dimension
		Real*8, intent(in) :: SolTheta(:) !! have NodG dimension
		integer :: i, j , e 
		
		GradTheta = 0.0d0
		do e = 1 , s%NodG
			do i = 1, s%NdimE
				GradTheta(i) = GradTheta(i) + SolTheta(e)*s%dPhi_G(i,e,nG)
			end do
		end do	
		
	end subroutine

	!! encapsulates the Theta (Theta scalar) calc, deprecated, you should use N matrix instead
	subroutine calcTheta(s,Theta,SolTheta,nG)
		class(PtGaussClass) :: s
		integer, intent(in) :: nG 
		Real*8, intent(out) :: Theta
		Real*8, intent(in) :: SolTheta(:) !! have NodG dimension 
		integer :: i, e 
		
		Theta= 0.0d0
		do e = 1 , s%NodG
			Theta = Theta + SolTheta(e)*s%Phi(e,nG)
		end do
		
	end subroutine
	
	subroutine getBmat(s,B,nG)
		class(PtGaussClass) :: s
		integer , intent(in) :: nG
		real*8 , intent(out) :: B(:,:) !! B is NdimE^2 X NodG*NdimE
		
		integer , parameter :: n2 = 4 !! assuming two dimension
		integer :: voigtMapI(n2), voigtMapJ(n2), p,q, i, j, qp
	
		voigtMapI = (/1,1,2,2/)
		voigtMapJ = (/1,2,1,2/)
		
		B = 0.0d0
	
		do p = 1 , n2
			i = voigtMapI(p)
			j = voigtMapJ(p)
			do q = 1 , s%NodG
				qp = (q-1)*s%NdimE + i
				B(p, qp) = s%dPhi_G(j,q,nG)				
			end do
		end do
		
!~ 		B(1,1) = s%dPhi_G(1,1,nG)
!~ 		B(1,3) = s%dPhi_G(1,2,nG)
!~ 		B(1,5) = s%dPhi_G(1,3,nG)
!~ 		
!~ 		B(2,1) = s%dPhi_G(2,1,nG)
!~ 		B(2,3) = s%dPhi_G(2,2,nG)
!~ 		B(2,5) = s%dPhi_G(2,3,nG)
!~ 
!~ 		B(3,2) = s%dPhi_G(1,1,nG)
!~ 		B(3,4) = s%dPhi_G(1,2,nG)
!~ 		B(3,6) = s%dPhi_G(1,3,nG)
!~ 		
!~ 		B(4,2) = s%dPhi_G(2,1,nG)
!~ 		B(4,4) = s%dPhi_G(2,2,nG)
!~ 		B(4,6) = s%dPhi_G(2,3,nG)
	
	end subroutine 
!~ 

	subroutine getBmatSym(s,B,nG) !! Following Hughes pag. 90
		class(PtGaussClass) :: s
		integer , intent(in) :: nG
		real*8 , intent(out) :: B(:,:) 
	
		B = 0.0d0
		
		B(1,1) = s%dPhi_G(1,1,nG)
		B(1,3) = s%dPhi_G(1,2,nG)
		B(1,5) = s%dPhi_G(1,3,nG)
		
		B(2,2) = s%dPhi_G(2,1,nG)
		B(2,4) = s%dPhi_G(2,2,nG)
		B(2,6) = s%dPhi_G(2,3,nG)

		B(3,1) = s%dPhi_G(2,1,nG)
		B(3,3) = s%dPhi_G(2,2,nG)
		B(3,5) = s%dPhi_G(2,3,nG)
		
		B(3,2) = s%dPhi_G(1,1,nG)
		B(3,4) = s%dPhi_G(1,2,nG)
		B(3,6) = s%dPhi_G(1,3,nG)
		
	end subroutine 
	
	subroutine getBmatScalar(s,B,nG)
		class(PtGaussClass) :: s
		integer , intent(in) :: nG
		real*8 , intent(out) :: B(:,:) !! B is NdimE X NodG
		
		B = s%dPhi_G(:,:,nG)	
	
	end subroutine 
 
 	subroutine getNvecScalar(s,N,nG)
		class(PtGaussClass) :: s
		integer , intent(in) :: nG
		real*8 , intent(out) :: N(:) !! N is NodG
		
		N = s%Phi(:,nG)	
	
	end subroutine 

 	subroutine getNmat(s,N,nG) !! particular case of triangle
		class(PtGaussClass) :: s
		integer , intent(in) :: nG
		real*8 , intent(out) :: N(:,:) !! N is NdimE x (NdimE*NodG)
 		
		N(1,1) = s%Phi(1,nG)
		N(2,2) = s%Phi(1,nG)
		N(1,3) = s%Phi(2,nG)
		N(2,4) = s%Phi(2,nG)
		N(1,5) = s%Phi(3,nG)
		N(2,6) = s%Phi(3,nG)
	
	end subroutine 


	subroutine setFEMtype(FEMtype,NodG,pOrder,NGP,iSimplex,iBubble)
		integer , intent(in) :: FEMtype
		integer , intent(out) :: NodG,pOrder,NGP,iSimplex,iBubble
		
		
		select case(FEMtype)
			case(1) !! linear triangle
				NGP = 1
				NodG = 3
				pOrder = 1
				iSimplex = 0
				iBubble = 0
		
			case(2) !! linear quadrangle
				NGP = 4
				NodG = 4
				pOrder = 1
				iSimplex = 1
				iBubble = 0
		
			case(3) !! quadratic triangle with 3 gauss points
				NGP = 3
				NodG = 6
				pOrder = 2
				iSimplex = 0
				iBubble = 0

			case(4) !! quadratic quadrangle with 9 gauss points
				NGP = 16
				NodG = 9
				pOrder = 2
				iSimplex = 1
				iBubble = 0	

			case(5) !! Linear Line with 1 gauss points
				NGP = 1
				NodG = 2
				pOrder = 1
				iSimplex = 0
				iBubble = 0	
			case(6) !! Linear Tetrahedra
				NGP = 4 !! just 1 or 3
				NodG = 4
				pOrder = 1
				iSimplex = 0
				iBubble = 0	
			case(7) !! Point
				NGP = 1 !! just 1 or 3
				NodG = 1
				pOrder = 0
				iSimplex = 0
				iBubble = 0	
			case(8) !! Quadratic Line with 2 gauss points
				NGP = 2
				NodG = 3
				pOrder = 2
				iSimplex = 0
				iBubble = 0	
			case default
				write(*,*) "FEMtype = " , FEMtype,  "not supported at moment"  
				stop
		end select
		
	end subroutine 

	subroutine setNodG(FEMtype,NodG)
		integer , intent(in) :: FEMtype
		integer , intent(out) :: NodG
		integer :: dummy
		
		call setFEMtype(FEMtype,NodG,dummy,dummy,dummy,dummy)

	end subroutine 
	
	subroutine GaussRuleQ(Psi,W,NdimE,NGP, iSimplex)
		integer , intent(in) :: NdimE, NGP,  iSimplex
		real*8, intent(out) :: Psi(NdimE,NGP), W(NGP)
		integer , parameter :: NGP1Dmax = 3, Ndim1D = 1
		integer :: nGi, nGj, nGk, NGP1D, nGpp
		real*8 :: Psi1D(1,NGP1Dmax), W1D(NGP1Dmax) 
	
		NGP1D = nint(real(NGP)**(1.0d0/real(NdimE)))
		
		Call GaussRule(Psi1D(:,1:NGP1D) ,W1D(1:NGP1D),Ndim1D,NGP1D,iSimplex)
		
		select case(NdimE)
			case(1)
				Psi(1,:) = Psi1D(1,1:NGP1D)
				W = W1D(1:NGP1D)
			case(2)
				do nGi = 1,NGP1D
				do nGj = 1,NGP1D
					nGpp = (nGi - 1)*NGP1D + nGj 
					Psi(1,nGpp) = Psi1D(1,nGi)
					Psi(2,nGpp) = Psi1D(1,nGj)
					W(nGpp) = W1D(nGi) * W1D(nGj)
				end do
				end do
			
			case(3)
				do nGi = 1,NGP1D
				do nGj = 1,NGP1D
				do nGk = 1,NGP1D
					nGpp = (nGi - 1)*NGP1D*NGP1D + (nGj-1)*NGP1D + nGk 
					Psi(1,nGpp) = Psi1D(1,nGi)
					Psi(2,nGpp) = Psi1D(1,nGj)
					Psi(3,nGpp) = Psi1D(1,nGk)
					W(nGpp) = W1D(nGi) * W1D(nGj) * W1D(nGk)
				end do
				end do
				end do
		end select
	
	
	end subroutine

end module
