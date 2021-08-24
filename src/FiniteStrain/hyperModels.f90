module hyperModelsLib
use funcAux
implicit none
private  Delfino, MooneyRivlin, Fibre, vol, calcInvCfromF, calcInvCfromC, calcInvCMfromF, calcInvCMfromC, Ciarlet , &
		DelfinoDpsi, DelfinoCompDPsi, MooneyRivlinDPsi, MooneyRivlinCompDPsi, MooneyRivlinComp2DPsi
public StrainEnergy, StrainEnergyC, strainEnergyDpsi

contains

subroutine Ciarlet(Psi,InvC,MatPar) !! page 259, Computation Inelasticity
	implicit none
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: InvC(:),MatPar(:)
	real*8 :: E, nu , lambda, mu, J, I1
	
	E = MatPar(1) 
	nu = MatPar(2)
 
	lambda = E*nu/((1.0d0 + nu)*(1.0d0 - 2.0d0*nu))
	mu = 0.5d0*E/(1.0d0 + nu)
	
	J = dsqrt(InvC(3)) ! J
	I1 = InvC(1)
	
	Psi = 0.25d0*(J*J-1.0d0) - (0.5d0*lambda + mu)*dlog(J) + 0.5d0*mu*(I1 - 3.0d0)

end subroutine


subroutine vol(Psi,InvC,MatPar)
	implicit none
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: InvC(:),MatPar
	real*8 :: c1,I3, J
	c1 = MatPar ! kappa

	J = dsqrt(InvC(3)) ! J
	
	Psi = c1*(J - 1.0d0)**2.0d0

end subroutine

!~ subroutine vol(Psi,invC,MatPar)
!~ 	implicit none
!~ 	real*8 , intent(out) :: Psi
!~ 	real*8 , intent(in) :: InvC(:),MatPar
!~ 	real*8 :: c1,J
!~ 	c1 = MatPar ! kappa
!~ 	
!~ 	J = dsqrt(invC(3)) ! J
!~ 	
!~ 	Psi = 0.5d0*c1*( J*J - 1.0d0 - 2.0d0*dlog(J) )
!~ 	return
!~ end subroutine

!~ subroutine vol(Psi,invC,MatPar)
!~ 	implicit none
!~ 	real*8 , intent(out) :: Psi
!~ 	real*8 , intent(in) :: InvC(:),MatPar
!~ 	real*8 :: c1,J
!~ 	c1 = MatPar ! kappa
!~ 	
!~ 	J = dsqrt(invC(3)) ! J
!~ 	
!~ 	Psi = 0.5d0*c1*( 0.5d0*(J*J - 1.0d0) - dlog(J) )
!~ 	return
!~ end subroutine

!~ 
!~ subroutine vol(Psi,invC,MatPar)
!~ 	implicit none
!~ 	real*8 , intent(out) :: Psi
!~ 	real*8 , intent(in) :: InvC(:),MatPar
!~ 	real*8 :: c1, J
!~ 	c1 = MatPar ! kappa
!~ 	
!~ 	J = dsqrt(InvC(3)) ! J
	
!~ 	Psi = c1*dlog(J)**2.0d0
!~ 	return
!~ end subroutine


!~ subroutine vol2(Psi,C,MatPar)
!~ 	implicit none
!~ 	real*8 , intent(out) :: Psi
!~ 	real*8 , intent(in) :: C(:,:),MatPar(2)
!~ 	real*8 :: c1,c2,I3
!~ 	
!~ 	c1 = 2.0d0*MatPar(1) ! kappa
!~ 	c2 = MatPar(2) ! kappa
!~ 	
!~ 	write(0,*) "New volumetric penalty function"
!~ 	
!~ 	call calcI3(I3,C)
!~ 	I3 = dlog(dsqrt(I3)) ! J
!~ 	
!~ 	Psi = -c1*I3 + c2*I3*I3
!~ 	return
!~ end subroutine

subroutine Delfino(Psi,invC,MatPar)
! Delfino et al. [10] proposed an (isotropic) rubber-like potential for carotid arteries
          ! [10] A. Delfino, N. Stergiopulos, J.E. Moore and J.-J. Meister, 
          ! Residual strain effects on the stress field in a thick wall
          ! finite element model of the human carotid bifurcation. 
          ! J. Biomech. 30 (1997) 777786.
	implicit none
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: invC(:),MatPar(2)
	real*8 :: c1,c2,c3,I1bar,J,rI33
	c1 = MatPar(1) ! a
	c2 = MatPar(2)*0.5d0 ! b/2
	
	J = dsqrt(InvC(3)) ! J
	rI33 = J**(-2.0D0/3.0D0)
	I1bar = rI33*InvC(1) ! I1bar
	
	Psi = 0.5d0*(c1/c2)*( Dexp(c2*(I1bar-3.0d0))-1.0d0 )
	
end subroutine

!~ subroutine MooneyRivlin(Psi,invC,MatPar)
!~  !Mooney Rivlin Bathe, K.J., Finite Element Procedures, pp.592; see also 6.4 pp.561
!~  !Crisfield M.A. Vol II, Chap 13, pp 62. See also, Section 10.3 pp 07 
!~ 	implicit none
!~ 	real*8 , intent(out) :: Psi
!~ 	real*8 , intent(in) :: invC(:),MatPar(2)
!~ 	real*8 :: c1,c2,I1bar,I2bar,J,rI33
!~ 	c1 = MatPar(1)
!~ 	c2 = MatPar(2)
!~ 
!~ 	J = dsqrt(InvC(3)) ! J
!~ 	rI33 = J**(-2.0D0/3.0D0)
!~ 	I1bar = rI33*InvC(1) ! I1bar
!~ 	I2bar = rI33*rI33*invC(2) ! I2bar
!~ 	Psi = c1*(I1bar-3.0d0) + c2*(I2bar-3.0d0)
!~ 	return
!~ end subroutine


subroutine MooneyRivlin(Psi,invC,MatPar)
 !Mooney Rivlin Bathe, K.J., Finite Element Procedures, pp.592; see also 6.4 pp.561
 !Crisfield M.A. Vol II, Chap 13, pp 62. See also, Section 10.3 pp 07 
	implicit none
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: invC(:),MatPar(2)
	real*8 :: c1,c2,I1bar,I2bar,J,rI33
	c1 = MatPar(1)
	c2 = MatPar(2)

	J = dsqrt(InvC(3)) ! J
	rI33 = J**(-2.0D0/3.0D0)
	I1bar = rI33*InvC(1) ! I1bar
	I2bar = rI33*rI33*invC(2) ! I2bar
	Psi = c1*(I1bar-3.0d0) + c2*(I2bar-3.0d0)
	return
end subroutine

subroutine MooneyRivlinModified(Psi,invC,MatPar)  !! the constants decrease with the deformation
	implicit none
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: invC(:),MatPar(2)
	real*8 :: c1,c2,I1bar,I2bar,J,rI33
	c1 = MatPar(1)
	c2 = MatPar(2)

	J = dsqrt(InvC(3)) ! J
	rI33 = J**(-2.0D0/3.0D0)
	I1bar = rI33*InvC(1) ! I1bar
	I2bar = rI33*rI33*invC(2) ! I2bar
	
	c1 = c1/J
	c2 = c2/J
	
	Psi = c1*(I1bar-3.0d0) + c2*(I2bar-3.0d0)
	return
end subroutine

subroutine Fibre(Psi,invCM,MatPar) ! see my dissertation
	implicit none
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: invCM(:),MatPar(2) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1,c2,I4
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k2
	
	I4 = invCM(2) !! already incompressible
	
	Psi = 0.5d0*(c1/c2)*( Dexp(c2*(I4-1.0d0)**2.0d0)-1.0d0 ) 
	
end subroutine 


subroutine MatrixRegularisation(Psi,invC,MatPar) 
	implicit none
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: invC(:),MatPar(4)
!~  	real*8 :: c1,c2, c3,c4, I1, J
!~ 	c1 = MatPar(1) 
!~ 	c2 = MatPar(2) 
!~ 	c3 = MatPar(3) 
!~ 	c4 = MatPar(4) 
	
!~ 	I1 = invC(1)
!~ 	J = dsqrt(invC(3))
	
!~ 	if(J>c4) then
!~ 		J = c4
!~ 	end if 

!~ 	if(I1>c3) then
!~ 		I1 = c3
!~ 	end if 
	
!~ 	Psi = c1*(I1 - 3.0d0) + c2*( J*J - 1.0d0 - 2.0d0*dlog(J) )  
	
	
	real*8 :: c1,c2,I1bar,I2bar,J,rI33, c3, c4
	c1 = MatPar(1)
	c2 = MatPar(2)
	c3 = MatPar(3)
	c4 = MatPar(4)

	J = dsqrt(InvC(3)) ! J
	rI33 = J**(-2.0D0/3.0D0)
	I1bar = rI33*InvC(1) ! I1bar
	I2bar = rI33*rI33*invC(2) ! I2bar

!~ 	if(I1bar > c2) then
!~ 		c1 = 0.0d0
!~ 	end if

!~ 	if(J > c4) then
!~ 		c3 = 0.0d0
!~ 	end if

!~ 	c1 = c1/(J**2.0d0)
!~ 	c3 = c3/(J**2.0d0)
	
	c1 = c1/J
	c3 = c3/J
	
	
	Psi = c1*(I1bar-3.0d0) + c3*(J - 1.0d0)**2.0d0
	
end subroutine 

subroutine calcInvCMfromF(invCM,F,beta)
	real*8 , intent(in) :: F(:,:), beta
	real*8 , intent(out) :: InvCM(:)
	real*8 :: v(size(F,1)), J , a0(size(F,1))
		
	a0(1) = dcos(beta)
	a0(2) = dsin(beta)
	
	v = matmul(F,a0)
	InvCM(1) = dot_product(v,v)
	
	call calcI3(J,F)
	
	InvCM(2) = (J**(-2.0D0/3.0D0))*InvCM(1) ! incompressible
	
end subroutine


subroutine calcInvCMfromC(invCM,C,beta)
	real*8 , intent(in) :: C(:,:), beta
	real*8 , intent(out) :: InvCM(:)
	real*8 :: v(size(C,1)), J, a0(size(C,1))
	
	a0(1) = dcos(beta)
	a0(2) = dsin(beta)
	
	v = matmul(C,a0)	
	InvCM(1) = dot_product(v,a0)
	
	call calcI3(J,C)
	J = dsqrt(J)
	
	InvCM(2) = (J**(-2.0D0/3.0D0))*InvCM(1) ! incompressible
	
end subroutine

	
subroutine calcInvCfromC(invC,C)
	real*8 , intent(in) :: C(:,:)
	real*8 , intent(out) :: InvC(:)
		
	call calcI1(InvC(1),C)
	call calcI2(InvC(2),C)
	call calcI3(InvC(3),C)
	
end subroutine

	
subroutine calcInvCfromF(invC,F)
	real*8 , intent(in) :: F(:,:)
	real*8 , intent(out) :: InvC(:)
	real*8 :: C(size(F,1),size(F,2))
		
	C = matmul(transpose(F),F)
	
	call calcI1(InvC(1),C)
	call calcI2(InvC(2),C)
	call calcI3(InvC(3),C)
	
end subroutine

subroutine strainEnergy(energy,F,matPar,constLaw)
	real*8, intent(in) :: F(:,:) , matPar(:)
	real*8 , intent(out) :: energy
	integer , intent(in) :: constLaw
	real*8 :: energyTemp, invC(3)
	
	if(constLaw<10) then 
		call calcInvCfromF(invC,F)
	else
		call calcInvCMfromF(invC,F,matPar(3))
	end if
	
	energy = 0.0d0
	energyTemp = 0.0d0
	
	selectcase(constLaw) 
		case(1) 
			call Delfino(energy,invC,matPar(1:2))
		case(2)
			call MooneyRivlin(energy,invC,matPar(1:2))
		case(3) 
			call Delfino(energy,invC,matPar(1:2))
			call vol(energyTemp,invC,matPar(3))
			energy = energy + energyTemp
		case(4)
			call MooneyRivlin(energy,invC,matPar(1:2))
            call vol(energyTemp,invC,matPar(3))
            energy = energy + energyTemp
		case(5)
			call Ciarlet(energy,invC,matPar(1:2))
		case(8)
			call MatrixRegularisation(energy,invC,matPar(1:4))
		case(9)
			call MooneyRivlinModified(energy,invC,matPar(1:2)) !! the constants decrease with the deformation
		case(13)
			call Fibre(energy,invC,matPar(1:2))
		case default
			write(0,*) constLaw, matpar(1:3), "Strain Energy case unknown"
			PAUSE
	end select

end subroutine

subroutine strainEnergyC(energy,C,matPar,constLaw)
	real*8, intent(in) :: C(:,:) , matPar(:)
	real*8 , intent(out) :: energy
	integer , intent(in) :: constLaw
	real*8 ::  energyTemp, invC(3)
	
	if(constLaw<10) then 
		call calcInvCfromC(invC,C)
	else
		call calcInvCMfromC(invC,C,matPar(3))
	end if
	
	energy = 0.0d0
	energyTemp = 0.0d0

	selectcase(constLaw) 
		case(1) 
			call Delfino(energy,invC,matPar(1:2))
		case(2)
			call MooneyRivlin(energy,invC,matPar(1:2))
		case(3) 
			call Delfino(energy,invC,matPar(1:2))
			call vol(energyTemp,invC,matPar(3))
			energy = energy + energyTemp
		case(4)
			call MooneyRivlin(energy,invC,matPar(1:2))
            call vol(energyTemp,invC,matPar(3))
            energy = energy + energyTemp
		case(5)
			call Ciarlet(energy,invC,matPar(1:2))
		case(9)
			call MooneyRivlinModified(energy,invC,matPar(1:2)) !! the constants decrease with the deformation
		case(13)
			call Fibre(energy,invC,matPar(1:2))
		case default
			write(0,*) constLaw, matpar(1:3), "Strain Energy case unknown"
			PAUSE
	end select

end subroutine



!! ============= quadratic ====================
subroutine volDPsi(DPsi,C,MatPar)
	implicit none
	real*8 , intent(out) :: DPsi(:)
	real*8 , intent(in) :: C(:,:),MatPar
	real*8 :: c1,I3
	c1 = MatPar
	call calcI3(I3,C)
	I3 = dsqrt(I3) ! J
	
	DPsi(3) = c1 - c1/I3
	DPsi(6) = 0.5d0*c1*I3**(-3.0d0)
	return
end subroutine

!! ============= mixed quad-log ======================
!~ subroutine vol(Psi,C,MatPar)
!~ 	implicit none
!~ 	real*8 , intent(out) :: Psi
!~ 	real*8 , intent(in) :: C(:,:),MatPar
!~ 	real*8 :: c1,I3
!~ 	c1 = MatPar ! kappa
!~ 		
!~ 	call calcI3(I3,C)
!~ 	I3 = dsqrt(I3) ! J
!~ 	
!~ 	Psi = c1*( I3*I3 - 1.0d0 - 2.0d0*dlog(I3) )
!~ 	return
!~ end subroutine
!~ 

! logarithm
!~ subroutine volDPsi(DPsi,C,MatPar)
!~ 	implicit none
!~ 	real*8 , intent(out) :: DPsi(:)
!~ 	real*8 , intent(in) :: C(:,:),MatPar
!~ 	real*8 :: c1,I3, logJ
!~ 	c1 = MatPar
!~ 	call calcI3(I3,C)
!~ 	I3 = dsqrt(I3) ! J
!~ 	logJ = dlog(I3)
	
!~ 	DPsi(3) = c1*logJ/(I3*I3)
!~ 	DPsi(6) = 0.5d0*c1*(1.0d0 - 2.0d0*logJ)*I3**(-4.0d0)
!~ 	return
!~ end subroutine

 
subroutine DelfinoDPsi(DPsi,C,MatPar) ! derivatives of Psi with respect to I1 and I2
	implicit none
	real*8 , intent(out) :: DPsi(:)
	real*8 , intent(in) :: C(:,:),MatPar(2)
	real*8 :: c1,c2,I1,rexp
	c1 = MatPar(1) ! a
	c2 = MatPar(2)*0.5d0 ! b/2

	call calcI1(I1,C)

	rexp = dexp(c2*(I1-3.0d0))

	DPsi(1) = 0.5d0*c1*rexp
	DPsi(3) = -I1*c1*rexp/(6.0d0)
	DPsi(4) = 0.5d0*c1*c2*rexp
	DPsi(6) = I1*c1*(I1*c2 + 4.0d0)*rexp/18.0d0
	DPsi(9) = -c1*(I1*c2 + 1.0d0)*rexp/(6.0d0)
	return
end subroutine

subroutine DelfinoCompDPsi(DPsi,C,MatPar)
	implicit none
	real*8 , intent(out) :: DPsi(:)
	real*8 , intent(in) :: C(:,:),MatPar(2)
	real*8 :: c1,c2,I1,I3,rI33,rexp
	c1 = MatPar(1) ! a
	c2 = MatPar(2)*0.5d0 ! b/2
	
	call calcI1(I1,C)
	call calcI3(I3,C)
	I3 = dsqrt(I3) ! J
	rI33 = I3**(-2.0D0/3.0D0)
	I1 = rI33*I1 ! I1bar
	rexp = dexp(c2*(I1-3.0d0))

	DPsi(1) = 0.5d0*c1*rexp*rI33
	DPsi(3) = -I1*c1*rexp/(6.0d0*I3*I3)
	DPsi(4) = 0.5d0*c1*c2*rexp*rI33*rI33
	DPsi(6) = I1*c1*(I1*c2 + 4.0d0)*rexp*I3**(-4.0d0)/18.0d0
	DPsi(9) = -c1*(I1*c2 + 1.0d0)*rexp/(6.0d0*I3**(8.0d0/3.0d0))
	return
end subroutine

subroutine MooneyRivlinDPsi(DPsi,C,MatPar) ! Derivatives of Psi with respect to I1 and I2
 	implicit none
	real*8 , intent(out) :: DPsi(:)
	real*8 , intent(in) :: C(:,:),MatPar(2)
	real*8 :: c1,c2
	c1 = MatPar(1)
	c2 = MatPar(2)
	DPsi(1) = c1
	DPsi(2) = c2
	return
end subroutine


subroutine MooneyRivlinCompDPsi(DPsi,C,MatPar)
 !Mooney Rivlin Bathe, K.J., Finite Element Procedures, pp.592; see also 6.4 pp.561
 !Crisfield M.A. Vol II, Chap 13, pp 62. See also, Section 10.3 pp 07 
	implicit none
	real*8 , intent(out) :: DPsi(:)
	real*8 , intent(in) :: C(:,:),MatPar(2)
	real*8 :: c1,c2,I1,I2,I3,rI33
	c1 = MatPar(1)
	c2 = MatPar(2)
	call calcI1(I1,C)
	call calcI2(I2,C)
	call calcI3(I3,C)
	I3 = dsqrt(I3) ! J
	rI33 = I3**(-2.0D0/3.0D0)
	I1 = rI33*I1 ! I1bar
	I2 = rI33*rI33*I2 ! I2bar
	
	DPsi(1) = c1*rI33
	DPsi(2) = c2*rI33*rI33	
	DPsi(3) = -(I1*c1+2.0d0*I2*c2)/(3.0d0*I3*I3)
	DPsi(6) = (4.0d0*I1*c1 + 10.0d0*I2*c2)/(9.0d0*I3**4.0d0)
	DPsi(8) = -2.0d0*c2*rI33**5.0d0/3.0d0	
	DPsi(9) = -c1*rI33**4.0d0/3.0d0	
	
	return
end subroutine

subroutine MooneyRivlinComp2DPsi(DPsi,C,MatPar)
 !Mooney Rivlin Holzapfel , Nonlinear Solid Mechanics pp.247 eq 6.147 
	implicit none
	real*8 , intent(out) :: DPsi(:)
	real*8 , intent(in) :: C(:,:),MatPar(3)
	real*8 :: c1,c2,c3,I1,I2,I3,rI33,d
	c1 = MatPar(1)
	c2 = MatPar(2)
	c3 = MatPar(3)
	d = 2.0d0*(c1 + 2.0d0*c2)
	call calcI1(I1,C)
	call calcI2(I2,C)
	call calcI3(I3,C)
	I3 = dsqrt(I3) ! J
	rI33 = I3**(-2.0D0/3.0D0)
	I1 = rI33*I1 ! I1bar
	I2 = rI33*rI33*I2 ! I2bar
	
	DPsi(1) = c1*rI33
	DPsi(2) = c2*rI33*rI33	
	DPsi(3) = -(I1*c1+2.0d0*I2*c2)/(3.0d0*I3*I3) + c3 - c3/I3 - 0.5d0*d/(I3**2.0d0)
	DPsi(6) = (4.0d0*I1*c1 + 10.0d0*I2*c2)/(9.0d0*I3**4.0d0) + c3/(2.0d0*I3**3.0d0) + 0.5d0*d/(I3**4.0d0)
	DPsi(8) = -2.0d0*c2*rI33**5.0d0/3.0d0	
	DPsi(9) = -c1*rI33**4.0d0/3.0d0
		
	return
end subroutine




subroutine strainEnergyDpsi(Dpsi,C,matPar,constLaw)
	real*8, intent(in) :: C(:,:) , matPar(:)
	real*8 , intent(out) :: Dpsi(12)
	integer , intent(in) :: constLaw
	real*8 :: DpsiTemp(12)
 
	Dpsi = 0.0d0
	DpsiTemp = 0.0d0
 
	selectcase(constLaw)	
		case(0)
			call volDpsi(Dpsi,C,matPar(3))
		case(1) 
			call DelfinoDpsi(Dpsi,C,matPar(1:2))
		case(2) 
			call MooneyRivlinDpsi(Dpsi,C,matPar(1:2))
		case(3) 
			call DelfinoCompDpsi(Dpsi,C,matPar(1:2))
			call volDpsi(DpsiTemp,C,matPar(3))
			DPsi = Dpsi + DpsiTemp
		case(4) 
			call MooneyRivlinCompDpsi(Dpsi,C,matPar(1:2))
			call volDpsi(DpsiTemp,C,matPar(3))
			DPsi = Dpsi + DpsiTemp
		case(6) 
			call MooneyRivlinComp2Dpsi(Dpsi,C,matPar(1:3))
		case default
			write(0,*) constLaw, "Strain Energy Derivatives case unknown"
			PAUSE
	end select

end subroutine



subroutine strainEnergyDpsi_iso(Dpsi,C,matPar,constLaw) !! improve this, just to solve the problem of taking out the volumetric contribution from DPsi
	real*8, intent(in) :: C(:,:) , matPar(:)
	real*8 , intent(out) :: Dpsi(12)
	integer , intent(in) :: constLaw
	real*8 :: DpsiTemp(12)
 
	Dpsi = 0.0d0
	DpsiTemp = 0.0d0
 
	selectcase(constLaw)	
		case(1) 
			call DelfinoDpsi(Dpsi,C,matPar(1:2))
		case(2) 
			call MooneyRivlinDpsi(Dpsi,C,matPar(1:2))
		case(3) 
			call DelfinoCompDpsi(Dpsi,C,matPar(1:2))
		case(4) 
			call MooneyRivlinCompDpsi(Dpsi,C,matPar(1:2))
		case(6) 
			call MooneyRivlinComp2Dpsi(Dpsi,C,matPar(1:3))
		case default
			write(0,*) constLaw, "Strain Energy Derivatives case unknown"
			PAUSE
	end select

end subroutine






end module

