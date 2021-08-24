module damageModelsNewLib
use funcAux

implicit none

contains


subroutine getH(H,Hder,r,damagePar)
	real*8 , intent(out) :: H, Hder
	real*8, intent(in) :: r , damagePar(:)
	integer  :: damageModel
	
	damageModel = nint(damagePar(2))
	
	select case(damageModel)
	
		case(1)
			call exp1Reg(H,Hder,r,damagePar)
		
		case(2)
			call exp2Reg(H,Hder,r,damagePar)
			
		case(3)
			call LinReg(H,Hder,r,damagePar)	
		
		case(4)
			call exp3dinf(H,Hder,r,damagePar)			
		case default
		
			H = 0.0d0
			
	end select

end subroutine

subroutine exp1Reg(H,Hder,r,damagePar)
	real*8 , intent(out) :: H , Hder
	real*8, intent(in) :: r , damagePar(:)
	real*8 :: r0, Gf, lmu, H0
	
	r0 = damagePar(3)
	Gf = damagePar(4)
	lmu = damagePar(5)
	
!~ 	H0 =  (0.5d0-Gf/(r0*r0*lmu))**(-1.0d0)
!~ 	write(0,*) 'H0 == ' ,  H0
	
	H0 =  -r0*r0*lmu/Gf
	
	H = H0*dexp(H0*(r- r0)/r0)

	Hder = H*(H0/r0)

end subroutine


subroutine exp2Reg(H,Hder,r,damagePar)
	real*8 , intent(out) :: H , Hder
	real*8, intent(in) :: r , damagePar(:)
	real*8 :: r0, Gf, lmu, beta
	
	r0 = damagePar(3)
	Gf = damagePar(4)
	lmu = damagePar(5)
		
	beta = 0.5*(-r0 + dsqrt((4.0d0*Gf/lmu)-r0*r0))
	write(0,*) beta
	
	H = dexp(-(r- r0)/beta)*(1.0d0 - r/beta)
	Hder = - dexp(-(r-r0)/beta)*(2.0d0 - r/beta)/beta
	
		
	write(0,*) r0,Gf,lmu
	
	pause 

end subroutine

subroutine LinReg(H,Hder,r,damagePar)
	real*8 , intent(out) :: H, Hder
	real*8, intent(in) :: r , damagePar(:)
	real*8 :: r0, Gf, lmu, H0

	r0 = damagePar(3)
	Gf = damagePar(4)
	lmu = damagePar(5)
	
	H0 = (1.0 - 2.0*Gf/(r0*r0*lmu))**(-1.0)

	H = H0
	Hder = 0.0d0
	
	write(0,*) r0,Gf,lmu
	
	pause 

end subroutine


subroutine exp3dinf(H,Hder,r,damagePar)
	real*8 , intent(out) :: H , Hder
	real*8, intent(in) :: r , damagePar(:)
	real*8 :: r0, dinf, alpha, E
	
	r0 = damagePar(3)
	alpha = damagePar(4)
	dinf = damagePar(5)
	
	E = dexp(-alpha*(r- r0)/r0)
	
	H = E*dinf*(1.0d0 - alpha*r/r0) + 1.0d0 - dinf
	
	Hder = - (alpha*E*dinf/r0)*(2.0d0 - alpha*r/r0)
	
	write(0,*) r0,alpha, dinf
	
	pause 

end subroutine

	
end module

