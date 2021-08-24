module fibresModelsLib
use funcAux

implicit none


real*8 , parameter :: kappaReg3 = 1.0d-5
real*8 , parameter :: kappaReg2 = 1.0d-8
real*8 , parameter :: kappaReg = 1.0d-8
real*8 , parameter :: rho = 1000.0d0

contains

subroutine derivatives_strainEnergy_fib(dPsi,d2Psi,qfib,matPar,constLaw)
	real*8, intent(in) :: qfib(:) , matPar(:)
	real*8 , intent(out) :: dPsi, d2Psi
	integer , intent(in) :: constLaw
	real*8 :: I4
	

	I4 = dot_product(qfib,qfib)

	select case(constLaw) 
		case(1) 
			call dFibreExponential(dPsi,I4,matPar(1:3))
			call d2FibreExponential(d2Psi,I4,matPar(1:3))
		case(2) 
			call dFibreLinear(dPsi,I4,matPar(1:3))
			call d2FibreLinear(d2Psi,I4,matPar(1:3))
		case(15) 
			call dFibreLinearComp(dPsi,I4,matPar(1:3))
			call d2FibreLinearComp(d2Psi,I4,matPar(1:3))
		case(21) 
			call dTrussQuang(dPsi,I4,matPar(1:3))
			call d2TrussQuang(d2Psi,I4,matPar(1:3))
		case(22) 
			call dFibreLinearComp_noRecruitment(dPsi,I4,matPar(1:3))
			call d2FibreLinearComp_noRecruitment(d2Psi,I4,matPar(1:3))
		case default
			write(0,*) constLaw, matpar(1:3), "Strain Energy case unknown"
			PAUSE
	end select

end subroutine

subroutine strainEnergy_fib(energy,qfib,matPar,constLaw)
	real*8, intent(in) :: qfib(:) , matPar(:)
	real*8 , intent(out) :: energy
	integer , intent(in) :: constLaw
	real*8 :: I4
	
	energy = 0.0d0

	I4 = dot_product(qfib,qfib)
	
	selectcase(constLaw) 
		case(1) 
			call FibreExponential(energy,I4,matPar(1:3))
		case(2) 
			call FibreLinear(energy,I4,matPar(1:3))
		case(3) 
			call FibreExponentialReg(energy,I4,matPar(1:3))
		case(4) 
			call FibreLinearReg(energy,I4,matPar(1:3))
		case(5) 
			call FibreExponentialReg2(energy,I4,matPar(1:3))
		case(6) 
			call FibreLinearReg2(energy,I4,matPar(1:3))
		case(7) 
			call FibreExponentialReg3(energy,I4,matPar(1:3))
		case(8) 
			call FibreLinearReg3(energy,I4,matPar(1:3))
		case(9) 
			call FibreExponential2(energy,I4,matPar(1:3))
		case(10) 
			call FibreLinear2(energy,I4,matPar(1:3))
		case(11) 
			call FibreExponentialVisc(energy,I4,matPar)
		case(12) 
			call FibreExponentialCompression(energy,I4,matPar)
		case(13) 
			call FibreExponentialNonlinear(energy,I4,matPar)
		case(14) 
			call FibreLinear3(energy,I4,matPar(1:3)) ! energy is linear
		case(15) 
			call FibreLinearComp(energy,I4,matPar(1:3)) ! energy is linear
		case(16) 
			call FibreLinearComp2(energy,I4,matPar(1:3)) ! energy is linear
		case(17) 
			call FibreLinearComp3(energy,I4,matPar(1:3)) ! energy is linear
		case(18) 
			call FibreLinearComp4(energy,I4,matPar(1:3)) ! energy is linear
		case(19) 
			call FibreZullinger(energy,I4,matPar(1:3)) ! energy is linear
		case(20) 
			call FibreLinearRegNew(energy,I4,matPar) ! energy is linear
		case(21) 
			call TrussQuang(energy,I4,matPar)
		case(22) 
			call FibreLinearComp_noRecruitment(energy,I4,matPar(1:3)) ! energy is linear
		case default
			write(0,*) constLaw, matpar(1:3), "Strain Energy case unknown"
			PAUSE
	end select

end subroutine


real*8 function heav(x) result(h)
	real*8, intent(in) :: x

	h = 0.5d0*(dsign(1.0d0,x) + 1.0d0)
	
!~ 	if(x>0.0d0) then
!~ 		h = 1.0d0
!~ 	else if(x<0.0d0) then
!~ 		h = 0.0d0
!~ 	else
!~ 		h = 0.5d0
!~ 	end if
!~ 	    
end function



real*8 function heavReg(x) result(h)
	real*8, intent(in) :: x
	
	h = 0.5d0*(dsign(1.0d0,x)*(1.0d0 - dexp(-dabs(x)*rho)) + 1.0d0)
!~ 	h = 0.5d0*(dsign(1.0d0,x)*(1.0d0 - dexp(-rho*dabs(x)**0.5d0)) + 1.0d0)
	
end function

subroutine TrussQuang(Psi,I4,MatPar) ! See Javier Bonet
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
	real*8 :: E
	E = MatPar(1) !
			
	Psi = (E/8.0d0)*(I4-1.0d0)**2.0d0 
	
end subroutine 

subroutine dTrussQuang(dPsi,I4,MatPar) ! See Javier Bonet
	real*8 , intent(out) :: dPsi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
	real*8 :: E
	E = MatPar(1) !
			
	dPsi = (E/4.0d0)*(I4-1.0d0) 
	
end subroutine 

subroutine d2TrussQuang(d2Psi,I4,MatPar) ! See Javier Bonet
	real*8 , intent(out) :: d2Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
	real*8 :: E
	E = MatPar(1) !
			
	d2Psi = E/4.0d0
	
end subroutine 

subroutine FibreExponentialReg3(Psi,I4,MatPar)
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1,c2,lfa, lrel, l
	
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k2
	lfa = MatPar(3)
	
	l = dsqrt(I4)
	lrel = l - lfa
	
	Psi = 0.5d0*(c1/c2)*( heav(lrel)*(Dexp(c2*lrel**2.0d0)-1.0d0) + kappaReg3*(Dexp(c2*l**2.0d0)-1.0d0))  
	
end subroutine 


subroutine FibreExponentialReg2(Psi,I4,MatPar)
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1,c2,lfa, lrel
	
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k2
	lfa = MatPar(3)
	
	lrel = dsqrt(I4) - lfa
	
	Psi = (heav(lrel) + kappaReg2)*0.5d0*(c1/c2)*( Dexp(c2*lrel**2.0d0)-1.0d0) 
	
end subroutine 


subroutine FibreExponentialReg(Psi,I4,MatPar)

	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1,c2,lfa, lrel
	
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k2
	lfa = MatPar(3)
	
	lrel = dsqrt(I4) - lfa
	
	Psi = heavReg(lrel)*0.5d0*(c1/c2)*( Dexp(c2*lrel**2.0d0)-1.0d0 + kappaReg) 
	
end subroutine 

subroutine FibreLinearReg3(Psi,I4,MatPar)
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1,lfa, lrel, l
	
	c1 = MatPar(1) ! k1
	lfa = MatPar(3)
	
	l = dsqrt(I4)
	lrel = l - lfa
	
	Psi = (c1/8.0d0)*( heav(lrel)*lrel**2.0d0 + kappaReg3*l**2.0d0)  
	
end subroutine 


subroutine FibreLinearReg2(Psi,I4,MatPar) ! see my Wriggers  Sx = E Ex
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
	real*8 :: c1,lfa, lrel
	c1 = MatPar(1) ! k1
	lfa = MatPar(3)
		
	lrel = dsqrt(I4) - lfa
	
	Psi = (heav(lrel) + kappaReg2)*(c1/8.0d0)*lrel**2.0d0 
	
end subroutine 

subroutine FibreLinearReg(Psi,I4,MatPar) ! see my Wriggers  Sx = E Ex
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
	real*8 :: c1,lfa, lrel
	c1 = MatPar(1) ! k1
	lfa = MatPar(3)
		
	lrel = dsqrt(I4) - lfa
	
	Psi = heavReg(lrel)*(c1/8.0d0)*(lrel**2.0d0 + kappaReg) 
	
end subroutine 

subroutine FibreExponential(Psi,I4,MatPar) ! see my dissertation
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1,c2,lfa, lrel
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k2
	lfa = MatPar(3)
	
	lrel = dsqrt(I4) - lfa
	
	Psi = heav(lrel)*0.5d0*(c1/c2)*( Dexp(c2*(lrel)**2.0d0)-1.0d0 ) 
	
end subroutine 

subroutine FibreExponentialCompression(Psi,I4,MatPar) ! see my dissertation
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1,c2,lfa, lrel
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k2
	lfa = MatPar(3)
	
	lrel = I4 - lfa*lfa
	
	Psi = 0.5d0*(c1/c2)*( Dexp(c2*(lrel)**2.0d0)-1.0d0 ) 
	
end subroutine 


subroutine FibreExponential2(Psi,I4,MatPar) ! full I4
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1,c2,lfa, lrel
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k2
	lfa = MatPar(3)
	
	lrel = I4 - lfa*lfa
	
!~ 	write(0,*) "fibreExpVisc = ", c1, c2, lfa
	
	Psi = heav(lrel)*0.5d0*(c1/c2)*( Dexp(c2*(lrel)**2.0d0)-1.0d0 ) 
	
end subroutine 

subroutine FibreExponentialVisc(Psi,I4,MatPar) ! with viscosity
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1,c2,lfa, lrel, nu_n, nu_k, ln, lk ,l
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k2
	lfa = MatPar(3)
	nu_n = MatPar(4)
	ln = MatPar(5)
	nu_k = MatPar(6)
	lk = MatPar(7)
	
!~ 	write(0,*) "fibreExpVisc = ", nu_n, ln, nu_k, lk, nu_k*(l - lk)**2.0d0, nu_n*(l - ln)**2.0d0
	
	lrel = I4 - lfa*lfa
	l = dsqrt(I4)
	
	Psi = heav(lrel)*0.5d0*(c1/c2)*( dexp(c2*(lrel)**2.0d0)-1.0d0 ) 
	Psi = Psi + nu_n*(l - ln)**2.0d0
	Psi = Psi + nu_k*(l - lk)**2.0d0

end subroutine 

subroutine FibreExponentialNonlinear(Psi,I4,MatPar) ! with viscosity
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1,c2,lfa, lrel, nu_n, nu_k, ln, lk, l 
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k2
	lfa = MatPar(3)
	nu_n = MatPar(4)
	ln = MatPar(5)
	nu_k = MatPar(6)
	lk = MatPar(7)
	
	
	lrel = I4 - lfa*lfa
	l = dsqrt(I4)
	
!~ 	write(0,*) "fibreExpVisc = ", nu_n, ln, nu_k, lk, l-lk, lrel
	
	
	Psi = heav(lrel)*0.5d0*(c1/c2)*( dexp(c2*(lrel)**2.0d0)-1.0d0 ) 
	Psi = Psi + heav(-lrel)*0.5d0*(c1/c2)*( dexp(c2*(lrel)**2.0d0)-1.0d0 )*nu_n*(l - ln)**2.0d0
	Psi = Psi + heav(-lrel)*0.5d0*(c1/c2)*( dexp(c2*(lrel)**2.0d0)-1.0d0 )*nu_k*(l - lk)**2.0d0

end subroutine 

subroutine FibreLinear(Psi,I4,MatPar) ! see my Wriggers  Sx = E Ex
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
	real*8 :: c1,lfa, lrel
	c1 = MatPar(1) ! k1
	lfa = MatPar(3)
		
	lrel = dsqrt(I4) - lfa
	
	Psi = heav(lrel)*(c1/8.0d0)*lrel**2.0d0 
	
end subroutine 


subroutine FibreLinear2(Psi,I4,MatPar) ! full I4
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
	real*8 :: c1,lfa, lrel
	c1 = MatPar(1) ! k1
	lfa = MatPar(3)
		
	lrel = I4 - lfa*lfa
	
	Psi = heav(lrel)*(c1/8.0d0)*lrel**2.0d0 
	
end subroutine 


subroutine dFibreExponential(dPsi,I4,MatPar) ! see my dissertation
	real*8 , intent(out) :: dPsi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1,c2,lfa, lrel
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k2
	lfa = MatPar(3)
	
	lrel = dsqrt(I4) - lfa
	
	dPsi = heav(lrel)*c1*lrel*Dexp(c2*(lrel)**2.0d0) 
	
end subroutine

subroutine d2FibreExponential(d2Psi,I4,MatPar) ! see my dissertation
	real*8 , intent(out) :: d2Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
	real*8 :: c1,c2,lfa, lrel
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k2
	lfa = MatPar(3)
	
	lrel = dsqrt(I4) - lfa
	
	d2Psi = heav(lrel)*c1*Dexp(c2*(lrel)**2.0d0)*(1.0d0 + 2.0d0*lrel**2.0d0) 
	
end subroutine

subroutine dFibreLinear(dPsi,I4,MatPar) ! see my Wriggers  Sx = E Ex
	real*8 , intent(out) :: dPsi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1,lfa, lrel
	c1 = MatPar(1) ! k1
	lfa = MatPar(3)
		
	lrel = dsqrt(I4) - lfa
	
	dPsi = heav(lrel)*(c1/4.0d0)*lrel 
	
end subroutine 

subroutine d2FibreLinear(d2Psi,I4,MatPar) ! see my Wriggers  Sx = E Ex
	real*8 , intent(out) :: d2Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1, lfa, lrel
	c1 = MatPar(1) ! k1
	lfa = MatPar(3)
	
	lrel = dsqrt(I4) - lfa
	
	d2Psi = heav(lrel)*c1/4.0d0 
	
end subroutine 


subroutine FibreLinear3(Psi,I4,MatPar) 

	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
	real*8 :: c1, c2, lfa, lrel
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k2
	lfa = MatPar(3)
		
	lrel = dsqrt(I4) - lfa
	
	Psi = heav(lrel)*(c1/8.0d0)*lrel**2.0d0 
	
	lrel = dsqrt(I4) - 1.0
	
	Psi = Psi + (c2/8.0d0)*lrel**2.0d0
	
end subroutine 

subroutine FibreLinearComp(Psi,I4,MatPar) ! works on compression

	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
	real*8 :: c1, c2, lfa, lrel
	c1 = MatPar(1) ! k1
	c2 = MatPar(2)
	lfa = MatPar(3)
		
	lrel = dsqrt(I4) - lfa
	
	if(lrel > 0.0d0) then 
		Psi = (c1/8.0d0)*lrel**2.0d0 
	else 
		Psi = (c2/8.0d0)*lrel**2.0d0
	end if
	
end subroutine 

subroutine FibreLinearComp_noRecruitment(Psi,I4,MatPar) ! works on compression

	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
	real*8 :: c1, c2, lrel
	c1 = MatPar(1) ! k1
	c2 = MatPar(2)
		
	lrel = dsqrt(I4) - 1.0d0
	
	if(lrel > 0.0d0) then 
		Psi = (c1/8.0d0)*lrel**2.0d0 
	else 
		Psi = (c2/8.0d0)*lrel**2.0d0
	end if
	
end subroutine 

subroutine FibreLinearComp2(Psi,I4,MatPar) ! works on compression

	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
	real*8 :: c1, c2, lfa, lrel
	c1 = MatPar(1) ! k1
	c2 = MatPar(2)
	lfa = MatPar(3)
		
	lrel = dsqrt(I4) - lfa
	
	if(lrel > 0.0d0) then 
		Psi = (c1/8.0d0)*lrel
	else 
		Psi = -(c2/8.0d0)*lrel
	end if
	
end subroutine 

subroutine FibreLinearComp3(Psi,I4,MatPar) ! works on compression

	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
	real*8 :: c1, c2, lfa, lrel
	c1 = MatPar(1) ! k1
	c2 = MatPar(2)
	lfa = MatPar(3)
		
	lrel = I4 - lfa*lfa
	
	if(lrel > 0.0d0) then 
		Psi = (c1/8.0d0)*lrel
	else 
		Psi = -(c2/8.0d0)*lrel
	end if
	
end subroutine

subroutine FibreLinearComp4(Psi,I4,MatPar) ! works on compression

	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
	real*8 :: c1, c2, lfa, lrel
	c1 = MatPar(1) ! k1
	c2 = MatPar(2)
	lfa = MatPar(3)
		
	lrel = I4 - lfa*lfa
	
	if(lrel > 0.0d0) then 
		Psi = (c1/8.0d0)*lrel**2.0d0 
	else 
		Psi = (c2/8.0d0)*lrel**2.0d0
	end if
	
end subroutine  


subroutine FibreZullinger(Psi,I4,MatPar) ! works on compression
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
	real*8 :: c1, c2, lfa, lrel
	c1 = MatPar(1) ! k1
	c2 = MatPar(2)
	lfa = MatPar(3)
		
	lrel = dsqrt(I4) - lfa
	
	if(lrel > 0.0d0) then 
		Psi = c1*(lrel - dlog(lrel + 1.0d0))
	else 
		Psi = -c2*lrel
	end if
	
end subroutine


subroutine FibreLinearRegNew(Psi,I4,MatPar) ! works on compression
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
	real*8 :: k1, k2, lfa, l1, l2, l3, eps, delta, l, dl21, dl32, alpha
	real*8 :: a0, a1, a2, b0, b1, b2, b3, c0, c1, c2
	real*8 , parameter :: tol = 1.0e-15
	
	k1 = MatPar(1) 
	k2 = MatPar(2)
	lfa = MatPar(3) !! reserved entry, given by param, should be 0.0 in commonPar
	l1 = MatPar(4)
	delta = MatPar(5)
	eps = MatPar(6)
	
	if(l1<0.0d0) then
		l1 = lfa
	end if

	l2 = lfa - delta*(lfa-l1)
	l3 = lfa + delta*(lfa-l1)
	
	l = dsqrt(I4)
	dl21 = l2 - l1
	dl32 = l3 - l2
	alpha = (k2 - k1)/dl32
		
	if((k2-k1)<tol .and. dl32 < tol) then
		alpha = 0.0d0
	else
		alpha = (k2 - k1)/dl32
	end if
	
	a0 = 0.0d0
	b0 = eps
	a2 = 0.5d0*k1
	b2 = 0.5d0*k1
	c2 = 0.5d0*k2
	b3 = alpha/6.0d0
	if((b0-a0)<tol .and. dl21 < tol) then
		a1 = 0.0d0
	else
		a1 = (b0 - a0)/dl21 - a2*dl21
	end if
	
	b1 = 2.0d0*a2*dl21 + a1
	c1 = 3.0d0*b3*dl32**2.0d0 + 2.0d0*b2*dl32 + b1
	c0 = b3*dl32**3.0d0 + b2*dl32**2.0d0 + b1*dl32 + b0
	
	Psi = 0.0d0
	
!~ 	if(l<l1) then 
!~ 		Psi = a0
	if(l < l2) then
		Psi = a2*(l-l1)**2.0d0 + a1*(l-l1) + a0
	else if (l < l3) then
		Psi = b3*(l-l2)**3.0d0 + b2*(l-l2)**2.0d0 + b1*(l-l2) + b0
	else
		Psi = c2*(l-l3)**2.0d0 + c1*(l-l3) + c0
	end if
	
end subroutine 



subroutine dFibreLinearComp(dPsi,I4,MatPar) ! see my Wriggers  Sx = E Ex
	real*8 , intent(out) :: dPsi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1, c2,lfa, lrel
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k2
	lfa = MatPar(3)
		
	lrel = dsqrt(I4) - lfa
	
	if(lrel >= 0.0d0) then 
		dPsi = (c1/8.0)*lrel/(lrel + lfa)
	else 
		dPsi = (c2/8.0)*lrel/(lrel + lfa)
	end if
	
	
end subroutine 

subroutine d2FibreLinearComp(d2Psi,I4,MatPar) ! see my Wriggers  Sx = E Ex
	real*8 , intent(out) :: d2Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1, c2, lfa, lrel
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k1
	lfa = MatPar(3)
	
	lrel = dsqrt(I4) - lfa
	
	if(lrel >= 0.0d0) then 
		d2Psi = (c1/16.0d0)*lfa*(lrel+lfa)**(-3.0d0) 
	else 
		d2Psi = (c2/16.0d0)*lfa*(lrel+lfa)**(-3.0d0)
	end if
	 
	
end subroutine 


subroutine dFibreLinearComp_noRecruitment(dPsi,I4,MatPar) ! see my Wriggers  Sx = E Ex
	real*8 , intent(out) :: dPsi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1, c2, lrel
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k2
		
	lrel = dsqrt(I4) - 1.0d0
	
	if(lrel >= 0.0d0) then 
		dPsi = (c1/8.0)*lrel/(lrel + 1.0d0)
	else 
		dPsi = (c2/8.0)*lrel/(lrel + 1.0d0)
	end if
	
	
end subroutine 

subroutine d2FibreLinearComp_noRecruitment(d2Psi,I4,MatPar) ! see my Wriggers  Sx = E Ex
	real*8 , intent(out) :: d2Psi
	real*8 , intent(in) :: I4, MatPar(:) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1, c2, lrel
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k1
	
	lrel = dsqrt(I4) - 1.0d0
	
	if(lrel >= 0.0d0) then 
		d2Psi = (c1/16.0d0)*1.0d0*(lrel+1.0d0)**(-3.0d0) 
	else 
		d2Psi = (c2/16.0d0)*1.0d0*(lrel+1.0d0)**(-3.0d0)
	end if
	 
	
end subroutine 


end module

