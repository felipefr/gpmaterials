module arcLengthLib

use funcAux

implicit none

private c1, c2, c3, Dlambda_acc, dlambda, arcLength
public lambda

real*8 , allocatable :: Gbar(:,:)
real*8 :: c1, c2, c3, Dlambda_acc, dlambda, arcLength, lambda, psi, GbarNormSq 


contains

subroutine initArcLengthLib(psi, GbarType, NdimE)
	real*8, intent(in) :: psi
	integer, intent(in) :: GbarType
	integer, intent(in) :: NdimE
	
	allocate(Gbar(NdimE,NdimE))
	
	Gbar = 0.0d0
	
	select case(GbarType)
		case(1)
			Gbar(1,1) = 1.0d0		
		case(2)
			Gbar(2,2) = 1.0d0
		case(3)
			Gbar(1,1) = 1.0d0		
	end select
	
	GbarNormSq = dot_product2(Gbar,Gbar)

end subroutine

subroutine contributeCoeficientsSpherical_arcLength(dUf_Bar,dUf_Star,DUf_acc,NodG)
	real*8 , intent(in) :: dUf_Bar(:), dUf_Star(:), DUf_acc(:)
	integer , intent(in) :: NodG
		
	c1 =  c1 + dot_product(dUf_Bar,dUf_Bar) + GbarNormSq*psi**2.0 
	c2 = c2 + 2.0d0*(dot_product(dUf_acc + dUf_Star, dUf_Bar) + Dlambda_acc*GbarNormSq*psi**2.0)
	c3 = c3 + dot_product(dUf_acc + dUf_Star, dUf_acc + dUf_Star) + &
		GbarNormSq*(Dlambda_acc*psi)**2.0 - arcLength*arcLength 

end subroutine

end module
