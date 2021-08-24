module loadingLib

use funcAux
implicit none

private  LinearPressure, ConstantPressure, SinusoidalPressure ,  SinusoidalPressureZero, adaptative, LinearPeriodic
public pressureLoad, chooseLoad

abstract interface
	subroutine PressureLoadLike(p,time,dt,loadPar) ! flexible load
		implicit none
		real*8 , intent(in) :: time,dt,loadPar(:)
		real*8 , intent(out) :: p
	end subroutine
end interface

procedure (PressureLoadLike), pointer :: pressureLoad => null ()

contains

subroutine setF(F,LoadPar,Time,dt,LoadProg,LoadType)
	Real*8 , intent(out):: F(:,:)
	Real*8 , intent(in):: Time , dt, LoadPar(:)
	integer , intent(in) :: LoadProg, LoadType 
	Real*8 :: beta(2), normal(2), tangent(2) ,  theta, alpha , lamb0, lamb1
	Real*8 :: lamb , eps
	integer :: k,l, dirShear
	
	call chooseLoad(LoadProg)
	call pressureLoad(lamb,time,dt,loadPar)
	
	F = 0.0d0
	
	selectcase(LoadType) 
		case(1) 
			F(1,1) = lamb
!~ 			write(0,*) "lxx=",lamb
		case(2) 
			F(2,2) = lamb
!~ 			write(0,*) "lyy=",lamb 
			
		case(3) 
			F(1,2) = lamb
!~ 			write(0,*) "lxy=",lamb
		case(4) 
			F(2,1) = lamb
!~ 			write(0,*) "lyx=",lamb
		case(5)
			F(1,2) = 0.5d0*lamb
			F(2,1) = 0.5d0*lamb
		case(6) ! biaxial
			F(1,1) = 0.5d0*lamb
			F(2,2) = 0.5d0*lamb
		case(7) !! incompressible !! 
			F(1,1) = lamb
			F(2,2) = -lamb/(1.0d0+lamb)
		case(8) !! incompressible
			F(1,1) = lamb
			F(1,2) = dsqrt(lamb)
			F(2,1) = F(1,2)	
		case(9) !! incompressible
			F(1,1) = lamb
			F(2,1) = 0.5*lamb
		case(10) !! prescribed
			F(1,1) = 1.05d0
			F(1,2) = 0.01d0
			F(2,1) = -0.02d0
			F(2,2) = 0.99d0
			
		case(11) !! beta \otimes normal , where beta = lamb * normal
			theta = datan(0.0d0/1.0d0)
			normal(1) = -dsin(theta)
			normal(2) = dcos(theta)
			
			F(1,1) = lamb * normal(1)* normal(1)
			F(1,2) = lamb * normal(1)* normal(2)
			F(2,1) = lamb * normal(2)* normal(1)
			F(2,2) = lamb * normal(2)* normal(2)

		case(12) !! beta \otimes normal , where beta = lamb * tangent
			theta = datan(0.0d0/1.0d0)
			normal(1) = -dsin(theta)
			normal(2) = dcos(theta)
			tangent(1) = dcos(theta)
			tangent(2) = dsin(theta)
			
			F(1,1) = lamb * tangent(1)* normal(1)
			F(1,2) = lamb * tangent(1)* normal(2)
			F(2,1) = lamb * tangent(2)* normal(1)
			F(2,2) = lamb * tangent(2)* normal(2)
			
		case(13) !! beta \otimes normal , where beta = lamb * ( b_t * tangent + b_n * normal)
!~ 			alpha = datan(0.0d0) + PI!! 
			alpha = 0.5d0*PI 
			theta = datan(0.3d0/1.0d0)
			normal(1) = -dsin(theta)
			normal(2) = dcos(theta)
			beta(1) = dcos(theta + alpha)
			beta(2) = dsin(theta + alpha)
	
			F(1,1) = lamb * beta(1)* normal(1)
			F(1,2) = lamb * beta(1)* normal(2)
			F(2,1) = lamb * beta(2)* normal(1)
			F(2,2) = lamb * beta(2)* normal(2)

		case(14) !! beta \otimes normal , where beta angle and normal angle is given by LoadPar
			theta = LoadPar(5)*PI
			alpha = LoadPar(6)*PI 
			normal(1) = dcos(theta)
			normal(2) = dsin(theta)
			beta(1) = dcos(alpha)
			beta(2) = dsin(alpha)
			
			F(1,1) = lamb * beta(1)* normal(1)
			F(1,2) = lamb * beta(1)* normal(2)
			F(2,1) = lamb * beta(2)* normal(1)
			F(2,2) = lamb * beta(2)* normal(2)

		case(15) 
			lamb0 = LoadPar(5)
			
			F = 0.0d0
			
			if(lamb<lamb0) then
				F(1,1) = lamb
				F(2,1) = 0.0d0			
			else 
				F(1,1) = lamb0
				F(2,1) = lamb - lamb0
			end if

		case(16)
			eps = LoadPar(4)
			k = nint(LoadPar(5))
			l = nint(LoadPar(6))
			
			F(1,1) = lamb
			
			if(k>0 .and. l> 0) then
				F(k,l) = F(k,l) + eps 
			end if

		case(17)
			lamb0 = LoadPar(5)
			lamb1 = LoadPar(6)
			
			F = 0.0d0
			
			theta = (lamb - lamb0)/(lamb1 - lamb0)
			
			if(theta<0.0d0) then
				F(1,1) = lamb
				F(2,1) = 0.0d0
			else if(theta < 1.0d0) then 
				F(1,1) = lamb0 + theta*lamb1
				 
				
			else ! > 1.0d0
				F(1,1) = lamb1
				F(2,1) = lamb - lamb1
			end if

			F(1,1) = (1.0d0 - theta)*lamb
			F(2,1) = theta*(lamb - lamb1)
		
		case default 
			write(0,*) 'Loadtype not recognised ' , LoadType
			stop
		
	end select
	
!~ 	call numprint(F)

end subroutine 

subroutine chooseLoad(LoadProg)
	integer, intent(in) :: LoadProg

	selectcase(LoadProg) 
		case(1) 
			pressureLoad => LinearPressure
		case(2) 
			pressureLoad => ConstantPressure
		case(3) 
			pressureLoad => SinusoidalPressure
		case(4) 
			pressureLoad => SinusoidalPressureZero
		case(5) 
			pressureLoad => adaptative
		case(6) 
			pressureLoad => LinearRampPressure
		case(7) 
			pressureLoad => LinearPeriodic
		case default 
			write(0,*) 'LoadProg not recognised ' , LoadProg
			stop
	end select

end subroutine

subroutine LinearPeriodic(p,time,dt,loadPar) ! periodic load, linear periodic
	implicit none
	real*8 , intent(in) :: time,dt,loadPar(:)
	real*8 , intent(out) :: p
	real*8 :: T,pmax0, dpmax, pmin0, dpmin
	
	real*8 :: di, trel, p0, p1
	integer :: i
	
	T = loadPar(1)
	pmax0 = loadPar(2)
	dpmax = loadPar(3)
	pmin0 = loadPar(4)
	dpmin = loadPar(5)
		
	di = dint(time/T)
	trel = time - di*T
	i = int(di)
	
	if(mod(i,2) == 0) then ! even number
		p0 = pmin0 + dpmin*di
		p1 = pmax0 + dpmax*(di+1.0d0)
	else ! odd number
		p0 = pmax0 + dpmax*di
		p1 = pmin0 + dpmin*(di+1.0d0)
	end if
		
	p = p0 + (p1 - p0)*trel/T
	
	write(0,*) "pressure=", p
	
end subroutine


subroutine LinearPressure(p,time,dt,loadPar) ! flexible load
	! p(t0) = 0 , p(t1) = pmax, linearly
	implicit none
	real*8 , intent(in) :: time,dt,loadPar(:)
	real*8 , intent(out) :: p
	real*8 :: t0,t1, pmax
	
	pmax = loadPar(1)
	t0 = loadPar(2)
	t1 = loadPar(3)
	
	if(time<t1) then
		p = (time-t0)*pmax/(t1-t0)
	else 
		p = pmax
	end if
	
	write(0,*) "pressure=", p
	
end subroutine

subroutine ConstantPressure(p,time,dt,loadPar) ! flexible load
	implicit none
	real*8 , intent(in) :: time,dt,loadPar(:)
	real*8 , intent(out) :: p
	
	p = loadPar(1)

end subroutine

subroutine SinusoidalPressure(p,time,dt,loadPar) ! flexible load
	implicit none
	real*8 , intent(in) :: time,dt,loadPar(:)
	real*8 , intent(out) :: p
	real*8 :: t0,t1, p_amp,p_ave, ncycles,phase

	p_ave = loadPar(1)
	p_amp = loadPar(2)
	t0 = loadPar(3)
	t1 = loadPar(4)
	ncycles = loadPar(5)
 	phase = loadPar(6)*2.0d0*PI 
	
	t1 = (time-t0)/(t1 - t0) ! t1 in [0,1]

 	p = p_ave + p_amp*dsin(ncycles*2.0d0*PI*t1+phase)

end subroutine

subroutine SinusoidalPressureZero(p,time,dt,loadPar) ! flexible load
	implicit none
	real*8 , intent(in) :: time,dt,loadPar(:)
	real*8 , intent(out) :: p
	real*8 :: t0,t1, p_amp,p_ave, ncycles,phase

	p_ave = loadPar(1)
	p_amp = loadPar(2)
	t0 = loadPar(3)
	t1 = loadPar(4)
	ncycles = loadPar(5)
	
 	phase = dasin(-p_ave/p_amp) + t0*ncycles*2.0d0*PI/(t1-t0) 
	
	t1 = (time-t0)/(t1 - t0) ! t1 in [0,1]

 	p = p_ave + p_amp*dsin(ncycles*2.0d0*PI*t1+phase)

end subroutine

subroutine adaptative(p,time,dt,loadPar) 
	implicit none
	real*8 , intent(in) :: time,dt,loadPar(:)
	real*8 , intent(out) :: p
	real*8 :: tf,lf,tm,kmin,kini , den, a1, a2, a3

	tf = loadPar(1)
	lf = loadPar(2)
	tm = loadPar(3)
	kmin = loadPar(4)
	kini = loadPar(5)
	
	den = 3.0d0*(tm**2.0d0)*tf**2.0d0 - 2.0d0*tm*tf**3.0d0
	a2 = (3.0d0*(tm**2.0d0)*(lf-kini*tf) - (tf**3.0d0)*(kmin - kini))/den
	a3 = (-2.0d0*tm*(lf-kini*tf) + (tf**2.0d0)*(kmin - kini))/den
	a1 = kini
	
	p = a1*time + a2*time**2.0 + a3*time**3.0
	
end subroutine

subroutine LinearRampPressure(p,time,dt,loadPar) ! flexible load
	! p(t0) = p0 , p(t1) = p1, linearly
	implicit none
	real*8 , intent(in) :: time,dt,loadPar(:)
	real*8 , intent(out) :: p
	real*8 :: t0,t1, p0, p1
	
	p0 = loadPar(1)
	p1 = loadPar(2)
	t0 = loadPar(3)
	t1 = loadPar(4)

	if(time<t0) then
		p = p0	
	else if(time<t1) then
		p = p0 + (time-t0)*(p1-p0)/(t1-t0)
	else 
		p = p1
	end if
	
	write(0,*) "pressure=", p, p0, p1, t0, t1

end subroutine


end module
