module fibresLib

use funcAux
use fibresModelsLib

implicit none

public calcSfib, calcDSfib, setFibreKinematic, setFibreKinematicSimple, &
		 getLambdaFibre, getLambdaFibreSimple, getEnergyFibre, calc_Kstiff_and_Res_TorsionalSpring, getTorsionalPsi

!~ 	call calcSfib(Sfib,qfib,matPar,constLaw)		
!~ 	call calcDSfib(DSfib,qfib,matPar,constLaw)	

contains

real*8 function getLambdaFibre(SolU,Xel,iFtype,NdimE) result(lambda) 
	Real*8 , intent(in) :: SolU(:) , Xel(:)
	integer , intent(in) :: NdimE, iFtype
	Real*8  :: DeltaU(NdimE), afib(NdimE), qfib(NdimE), Lfib
	
	call setFibreKinematic(qfib,DeltaU,afib,Lfib,SolU,Xel,iFtype)
	
	lambda = norm2(qfib)
	
end function

real*8 function getEnergyFibre(SolU,Xel,MatPar,constLaw,iFtype,NdimE) result(energy) 
	Real*8 , intent(in) :: SolU(:) , Xel(:), MatPar(:)
	integer , intent(in) :: NdimE, iFtype, constLaw
	Real*8  :: DeltaU(NdimE), afib(NdimE), qfib(NdimE), Lfib
	
	call setFibreKinematic(qfib,DeltaU,afib,Lfib,SolU,Xel,iFtype)
	
!~ 	write(0,*) 'qfib ===== ' , matPar(1), matPar(2) , matPar(3)
!~ 	if (norm2(qfib) > matPar(3) )then
!~ 		write(0,*)	(900.0/8.0d0)*(norm2(qfib) - matPar(3))**2.0d0
!~ 	end if
	call strainEnergy_fib(energy,qfib,MatPar,constLaw)
	
end function


real*8 function getLambdaFibreSimple(SolU,afib,Lfib,iFtype,NdimE) result(lambda) 
	Real*8 , intent(in) :: SolU(:) , afib(:), Lfib
	integer , intent(in) :: NdimE, iFtype
	Real*8  :: DeltaU(NdimE), qfib(NdimE)
	
	call setFibreKinematicSimple(qfib,DeltaU,afib,Lfib,SolU,iFtype)
	
	lambda = norm2(qfib)
	
end function

subroutine get_afib(afib,Lfib,Xel)
	use globalVariables, only : NdimE
	Real*8 , intent(out) :: afib(:), Lfib
	Real*8 , intent(in) :: Xel(:)
	integer i 		
	
	do i = 1, NdimE
		afib(i) = Xel(NdimE + i) - Xel(i)
	end do
	
	Lfib = norm2(afib)
	afib = afib/Lfib

end subroutine

subroutine get_afibbar(afibbar,Xel)
	use globalVariables, only : NdimE
	Real*8 , intent(out) :: afibbar(:)
	Real*8 , intent(in) :: Xel(:)
	integer i 
			
	do i = 1, NdimE
		afibbar(i) = 0.5d0*(Xel(NdimE + i) + Xel(i))
	end do
	
end subroutine


subroutine setFibreKinematic(qfib,DeltaU,afib,Lfib,SolU,Xel,iFtype)
	use globalVariables, only : getF, NdimE
	
	Real*8 , intent(out) :: DeltaU(:), afib(:), qfib(:), Lfib
	Real*8 , intent(in) :: SolU(:) , Xel(:)
	integer , intent(in) :: iFtype
	real*8 :: F(NdimE,NdimE)
	integer i
	
	call getF(F,iFtype)

	call get_afib(afib,Lfib,Xel)
	
	do i = 1, NdimE
		DeltaU(i) = SolU(NdimE+i) - SolU(i)
	end do
	
	qfib = matmul(F,afib) + DeltaU/Lfib	

end subroutine

subroutine setFibreKinematicSimple(qfib,DeltaU,afib,Lfib,SolU,iFtype)
	use globalVariables, only : getF, NdimE
	
	Real*8 , intent(out) :: DeltaU(:), qfib(:)
	Real*8 , intent(in) :: SolU(:) , afib(:), Lfib
	integer , intent(in) :: iFtype
	real*8 :: F(NdimE,NdimE)
	integer i
	
	call getF(F,iFtype)
!~ 	call numprint(F)

	do i = 1, NdimE
		DeltaU(i) = SolU(NdimE+i) - SolU(i)
	end do
	
	qfib = matmul(F,afib) + DeltaU/Lfib	

end subroutine

subroutine setFibreKinematicSimple2(qfib,F,SolUf,afib,Lfib)
	use globalVariables, only : getF, NdimE
	
	Real*8 , intent(out) :: qfib(:)
	Real*8 , intent(in) :: SolUf(:) , afib(:), F(:,:), Lfib
	integer :: i
	
	qfib = matmul(F,afib)
	
	do i = 1,NdimE
		qfib(i) = qfib(i) + (SolUf(NdimE + i) - SolUf(i))/Lfib 
	end do
	
end subroutine


subroutine calcSfib(Sfib,qfib,matPar,constLaw,NdimE)
	integer , intent(in) :: constLaw, NDimE
	real*8, intent(in) :: MatPar(:), qfib(NDimE)
	real*8, intent(out) :: Sfib(NDimE)
	logical , parameter :: isAnalitic = .false.
	
	if(isAnalitic) then
		call SfibAnalytic(Sfib,qfib,matPar,constLaw)
	else
		call calcSfibNumeric(Sfib,qfib,matPar,constLaw,NdimE)
	end if

end subroutine


subroutine calcDSfib(DSfib,qfib,matPar,constLaw,NdimE)
	integer , intent(in) :: constLaw, NDimE
	real*8, intent(in) :: MatPar(:), qfib(NDimE)
	real*8, intent(out) :: DSfib(NDimE,NdimE)
	logical , parameter :: isAnalitic = .false.
	
	if(isAnalitic) then
		call DSfibAnalytic(DSfib,qfib,matPar,constLaw)
	else
		call calcDSfibNumeric(DSfib,qfib,matPar,constLaw,NdimE)
	end if

end subroutine

subroutine calcSfibNumeric(Sfib,qfib,matPar,constLaw,NdimE)
	integer , intent(in) :: constLaw, NDimE
	real*8, intent(in) :: MatPar(:), qfib(NDimE)
	real*8, intent(out) :: Sfib(NDimE)
	real*8 , parameter :: eps = 1.d-4
	real*8 ::  qfibP(NDimE), inv2eps, energyPr , energyPl
	integer :: i 

	inv2eps = 0.5d0/eps
	
	Do i=1,NdimE    
		qfibP=qfib
		qfibP(i) = qfibP(i) + eps
		Call strainEnergy_fib(energyPr,qfibP,MatPar,constLaw)

		qfibP=qfib
		qfibP(i) = qfibP(i) - eps
		Call strainEnergy_fib(energyPl,qfibP,MatPar,constLaw)

		Sfib(i) = (energyPr - energyPl)*inv2eps
	Enddo
	
end subroutine

subroutine calcDSfibNumeric(DSfib,qfib,matPar,constLaw,NdimE)
	integer , intent(in) :: constLaw, NDimE
	real*8, intent(in) :: MatPar(:), qfib(NDimE)
	real*8, intent(out) :: DSfib(NDimE,NdimE)
	real*8 , parameter :: eps = 1.d-4
	real*8 ::  qfibP(NDimE), inveps2, inv4eps2, energy, energyPr , energyPl,  &
				energyPrr , energyPrl , energyPlr, energyPll 
	integer :: i,j 

	Call strainEnergy_fib(energy,qfib,MatPar,constLaw)

	inveps2 = 1.0d0/(eps*eps)
	inv4eps2 = 0.25d0/(eps*eps)	

	Do i=1,NdimE    
		qfibP=qfib
		qfibP(i) = qfibP(i) + eps
		Call strainEnergy_fib(energyPr,qfibP,MatPar,constLaw)

		qfibP=qfib
		qfibP(i) = qfibP(i) - eps
		Call strainEnergy_fib(energyPl,qfibP,MatPar,constLaw)

		DSfib(i,i) = (energyPr + energyPl - 2.0d0*energy)*inveps2
	Enddo

	Do i=1,NdimE    
	Do j=i+1,NdimE        
		qfibP=qfib
		qfibP(i) = qfibP(i) + eps
		qfibP(j) = qfibP(j) + eps
		Call strainEnergy_fib(energyPrr,qfibP,MatPar,constLaw)

		qfibP=qfib
		qfibP(i) = qfibP(i) + eps
		qfibP(j) = qfibP(j) - eps
		Call strainEnergy_fib(energyPrl,qfibP,MatPar,constLaw)

		qfibP=qfib
		qfibP(i) = qfibP(i) - eps
		qfibP(j) = qfibP(j) + eps
		Call strainEnergy_fib(energyPlr,qfibP,MatPar,constLaw)

		qfibP=qfib
		qfibP(i) = qfibP(i) - eps
		qfibP(j) = qfibP(j) - eps
		Call strainEnergy_fib(energyPll,qfibP,MatPar,constLaw)

		DSfib(i,j) = (energyPrr + energyPll - energyPrl - energyPlr)*inv4eps2
		DSfib(j,i) = DSfib(i,j)
	Enddo
	Enddo
 	
end subroutine

!~ subroutine SfibAnalytic(Sfib,qfib,matPar,constLaw)
!~ 	real*8, intent(in) :: qfib(:) , matPar(:)
!~ 	real*8 , intent(out) :: Sfib(:)
!~ 	integer , intent(in) :: constLaw
!~ 	real*8 :: I4, dPsi
	
!~ 	Sfib = 0.0d0
	
!~ 	I4 = dot_product(qfib,qfib)
	
!~ 	select case(constLaw) 
!~ 		case(1) 
!~ 			call dFibreExponential(dPsi,I4,matPar(1:3))
!~ 		case(2) 
!~ 			call dFibreLinear(dPsi,I4,matPar(1:3))
!~ 		case(15) 
!~ 			call dFibreLinearComp(dPsi,I4,matPar(1:3))
!~ 		case default
!~ 			write(0,*) constLaw, matpar(1:3), "Strain Energy case unknown"
!~ 			PAUSE
!~ 	end select
	
!~ 	Sfib = 2.0d0*dPsi*qfib

!~ end subroutine

!~ subroutine DSfibAnalytic(DSfib,qfib,matPar,constLaw)
!~ 	real*8, intent(in) :: qfib(:) , matPar(:)
!~ 	real*8 , intent(out) :: DSfib(:,:)
!~ 	integer , intent(in) :: constLaw
!~ 	real*8 :: I4, dPsi, d2Psi
	
!~ 	DSfib = 0.0d0
	
!~ 	I4 = dot_product(qfib,qfib)
	
!~ 	select case(constLaw) 
!~ 		case(1) 
!~ 			call dFibreExponential(dPsi,I4,matPar(1:3))
!~ 			call d2FibreExponential(d2Psi,I4,matPar(1:3))
!~ 		case(2) 
!~ 			call dFibreLinear(dPsi,I4,matPar(1:3))
!~ 			call d2FibreLinear(d2Psi,I4,matPar(1:3))
!~ 		case(15) 
!~ 			call dFibreLinearComp(dPsi,I4,matPar(1:3))
!~ 			call d2FibreLinearComp(d2Psi,I4,matPar(1:3))
!~ 		case default
!~ 			write(0,*) constLaw, matpar(1:3), "Strain Energy case unknown"
!~ 			PAUSE
!~ 	end select
	
!~ 	DSfib(1,1) = 4.0d0*d2Psi*qfib(1)*qfib(1) + 2.0d0*dPsi
!~ 	DSfib(1,2) = 4.0d0*d2Psi*qfib(1)*qfib(2)
!~ 	DSfib(2,1) = 4.0d0*d2Psi*qfib(2)*qfib(1)
!~ 	DSfib(2,2) = 4.0d0*d2Psi*qfib(2)*qfib(2) + 2.0d0*dPsi

!~ end subroutine


subroutine calcSfibDSfib_analytic(Sfib,DSfib,qfib,matPar,constLaw)
	use globalVariables, only : NdimE
	real*8, intent(in) :: qfib(:) , matPar(:)
	real*8 , intent(out) :: Sfib(:), DSfib(:,:)
	integer , intent(in) :: constLaw
	real*8 :: dPsi, d2Psi
	
	call derivatives_strainEnergy_fib(dPsi,d2Psi,qfib,matPar,constLaw)
	
	Sfib = 2.0d0*dPsi*qfib
	
!~ 	DSfib(1,1) = 4.0d0*d2Psi*qfib(1)*qfib(1) + 2.0d0*dPsi
!~ 	DSfib(1,2) = 4.0d0*d2Psi*qfib(1)*qfib(2)
!~ 	DSfib(2,1) = 4.0d0*d2Psi*qfib(2)*qfib(1)
!~ 	DSfib(2,2) = 4.0d0*d2Psi*qfib(2)*qfib(2) + 2.0d0*dPsi

	DSfib = 0.0d0	
	call VotimesW(DSfib,qfib,qfib)
	DSfib = 4.0d0*d2Psi*DSfib + 2.0d0*dPsi*deltaKron(1:NdimE,1:NdimE)

end subroutine

real*8 function getTorsionalPsi(X,U,eta) result(Psi) ! just for 2D
	real*8 , intent(in) :: eta, X(:), U(:)
	integer, parameter :: NdimE = 2, NodG =3
	real*8 :: M1(NdimE), M2(NdimE), q1(NdimE), q2(NdimE), &
			  L1, L2, a, b, c, Dang 

	M1 = X(3:4) - X(1:2)
	M2 = X(5:6) - X(1:2)
	L1 = norm2(M1)
	L2 = norm2(M2)

	M1 = M1/L1
	M2 = M2/L2
	
	q1 = M1 + (U(3:4) - U(1:2))/L1
	q2 = M2 + (U(5:6) - U(1:2))/L2
	
	q1 = q1/norm2(q1)
	q2 = q2/norm2(q2)

	a = dot_product(q1,q2)
	b = dot_product(M1,M2)

	c = a*b + dsqrt((1.0 - a*a)*(1.0 - b*b))
	Dang = dacos(c)

	Psi = 0.5*eta*Dang**2.0

end function

!!! Torsional Springs
subroutine calc_Kstiff_and_Res_TorsionalSpring(Kstiff, Res, X,U,eta) ! just for 2D

	integer, parameter :: NdimE=2, NodG = 3
	real*8, intent(in) :: X(:), U(:), eta
	real*8, intent(out) :: Kstiff(NdimE*NodG,NdimE*NodG) , Res(NdimE*NodG)
	real*8 , parameter :: eps = 1.d-5
	real*8 ::  Up(NodG*NDimE), inv2eps, inveps2, inv4eps2, &
			energy, energyPr , energyPl, energyPrr , energyPrl , energyPlr, energyPll 
	integer :: i ,j, n 

	inv2eps = 0.5d0/eps
	inveps2 = 1.0d0/(eps*eps)
	inv4eps2 = 0.25d0/(eps*eps)	

	energy = getTorsionalPsi(X,U,eta)

	n = NdimE*NodG
	Do i=1,n    
		Up=U
		Up(i) = Up(i) + eps
		energyPr = getTorsionalPsi(X,Up,eta)
		
		Up=U
		Up(i) = Up(i) - eps
		energyPl = getTorsionalPsi(X,Up,eta)

		Res(i) = (energyPr - energyPl)*inv2eps
		
		Kstiff(i,i) = (energyPr + energyPl - 2.0d0*energy)*inveps2
		
		Do j=i+1,n       
			Up=U
			Up(i) = Up(i) + eps
			Up(j) = Up(j) + eps
			energyPrr = getTorsionalPsi(X,Up,eta)

			Up=U
			Up(i) = Up(i) + eps
			Up(j) = Up(j) - eps
			energyPrl = getTorsionalPsi(X,Up,eta)

			Up=U
			Up(i) = Up(i) - eps
			Up(j) = Up(j) + eps
			energyPlr = getTorsionalPsi(X,Up,eta)

			Up=U
			Up(i) = Up(i) - eps
			Up(j) = Up(j) - eps
			energyPll = getTorsionalPsi(X,Up,eta)

			Kstiff(i,j) = (energyPrr + energyPll - energyPrl - energyPlr)*inv4eps2
			Kstiff(j,i) = Kstiff(i,j)
		Enddo
	Enddo
	
end subroutine


end module

