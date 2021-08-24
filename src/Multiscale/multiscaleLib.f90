module multiscaleLib
use funcAux
use globalVariables , only : NdimE
implicit none

public seeDetQ, calcQ, getCtang, computeQ, getMinDetQ

real*8 , allocatable :: Dhom(:,:,:,:) , PKhom(:,:)


contains 

subroutine initMultiscaleLib()
	allocate(Dhom(NdimE,NdimE,NdimE,NdimE) , PKhom(NdimE,NdimE))
end subroutine

subroutine getCtang(C,Sol1,NdimE,idofT,NodElt,iShiftC)
	integer , intent(in) :: NdimE, idofT, NodElt,iShiftC
	real*8, intent(in) :: Sol1(idofT*NodElt)
	real*8, intent(out) :: C(NdimE,NdimE,NdimE,NdimE)
	integer :: A, ApRow, i, j, k, l,  iAux
	
	A = NodElt
	ApRow  = (A-1)*iDofT + iShiftC
	
	do i = 1 , NdimE
	do j = 1 , NdimE
	do k = 1 , NdimE
	do l = 1 , NdimE
		iAux = ApRow + iShiftIJKL(i,j,k,l,NdimE) + 1
		C(i,j,k,l) = Sol1(iAux) 
	end do
	end do
	end do
	end do

end subroutine 

subroutine computeQ(Q,C,normal,NdimE)
	real*8 , intent(out) :: Q(NdimE,NdimE)
	real*8 , intent(in) :: C(NdimE,NdimE,NdimE,NdimE), normal(NdimE)
	integer , intent(in) :: NdimE
	integer :: i,j,k,l
	
	Q = 0.0d0
	do i = 1 , NdimE
	do j = 1 , NdimE
	do k = 1 , NdimE
	do l = 1 , NdimE
		Q(i,k) = Q(i,k) + C(i,j,k,l)*normal(j)*normal(l) 
	end do
	end do
	end do
	end do

end subroutine

subroutine getMinDetQ(minDetQ,thetaCrit,betaCrit,Sol1,NdimE,idofT,NodElt,iShiftC,nStepTheta) !!! Silent version
	integer , intent(in) :: NdimE, idofT, NodElt,iShiftC, nStepTheta
	real*8, intent(in) :: Sol1(idofT*NodElt)
	real*8, intent(out) :: minDetQ,thetaCrit, betaCrit
	real*8 :: C(NdimE,NdimE,NdimE,NdimE), Q(NdimE,NdimE), normal(NdimE), detQ, theta, stepTheta
	integer :: n
	
	call getCtang(C,Sol1,NdimE,idofT,NodElt,iShiftC)
	
	stepTheta = PI/real(NstepTheta)
	
	theta = 0.0d0
	
	do n = 1 , NstepTheta
		
		normal(1) = dcos(theta); normal(2) = dsin(theta)		
		
		call computeQ(Q,C,normal,NdimE)
		
		detQ = Q(1,1)*Q(2,2) - Q(1,2)*Q(2,1) !!! just in case of NdimE == 2
			
		if(n==1 .or. detQ<minDetQ) then
			minDetQ = detQ
			thetaCrit = theta
			betaCrit = datan2(Q(1,1),-Q(1,2))
		end if
				
		theta = theta + stepTheta
		
	end do

end subroutine


!~ subroutine getMinDetQ(minDetQ,thetaCrit,Sol1,NdimE,idofT,NodElt,iShiftC,nStepTheta)
!~ 	integer , intent(in) :: NdimE, idofT, NodElt,iShiftC, nStepTheta
!~ 	real*8, intent(in) :: Sol1(idofT*NodElt)
!~ 	real*8, intent(out) :: minDetQ,thetaCrit
!~ 	real*8 :: C(NdimE,NdimE,NdimE,NdimE), Q(NdimE,NdimE), normal(NdimE), detQ, theta, stepTheta
!~ 	integer :: n
!~ 	integer , parameter :: OUnit_histDetQ = 35
!~ 	 
!~ 	open (OUnit_histDetQ, file='histDetQ.txt', Access = 'append')  
!~ 	
!~ 	call getCtang(C,Sol1,NdimE,idofT,NodElt,iShiftC)
!~ 	
!~ 	stepTheta = PI/real(NstepTheta)
!~ 	
!~ 	theta = 0.0d0
!~ 	
!~ 	minDetQ = 1e10
!~ 	
!~ 	do n = 1 , NstepTheta
!~ 		
!~ 		normal(1) = dcos(theta); normal(2) = dsin(theta)		
!~ 		
!~ 		call computeQ(Q,C,normal,NdimE)
!~ 		
!~ 		detQ = Q(1,1)*Q(2,2) - Q(1,2)*Q(2,1) !!! just in case of NdimE == 2
!~ 		
!~ 		if(detQ<minDetQ) then
!~ 			minDetQ = detQ
!~ 			thetaCrit = theta
!~ 		end if
!~ 				
!~ 		
!~ 		write(OUnit_histDetQ,*) detQ
!~ 		write(OUnit_histDetQ,*) theta
!~ 
!~ 		theta = theta + stepTheta		
!~ 		
!~ 	end do
!~ 
!~ 	close(OUnit_histDetQ)
!~ 	
!~ end subroutine

subroutine seeDetQ(Sol1,NdimE,idofT,NodElt)
	integer , intent(in) :: NdimE, idofT, NodElt
	real*8, intent(in) :: Sol1(idofT*NodElt)
	real*8 :: C(NdimE,NdimE,NdimE,NdimE)
	integer :: A, ApRow, i, j, k, l,  iShiftC, iAux
	integer , parameter :: OUnitChom = 47
	integer , parameter :: OUnitP11hom = 48, OUnitP12hom = 49, OUnitP21hom = 50, OUnitP22hom = 51
	
	open (OUnitChom, file='Chom.txt', Access = 'append')  
	open (OUnitP11hom, file='P11hom.txt', Access = 'append')  
	open (OUnitP12hom, file='P12hom.txt', Access = 'append')  
	open (OUnitP21hom, file='P21hom.txt', Access = 'append')  
	open (OUnitP22hom, file='P22hom.txt', Access = 'append')  

	iShiftC = 0
	
	A = NodElt
	ApRow  = (A-1)*iDofT + iShiftC
	
	do i = 1 , NdimE
	do j = 1 , NdimE
	do k = 1 , NdimE
	do l = 1 , NdimE
		iAux = ApRow + iShiftIJKL(i,j,k,l,NdimE) + 1
		C(i,j,k,l) = Sol1(iAux) 
	end do
	end do
	end do
	end do
	
!~ 	call printSebaNotation(C)
!~ 	call writeSebaNotation(C,OUnitChom)
!~ 	
!~ 	write(OUnitP11hom,*) Sol1(ApRow + 1 )
!~ 	write(OUnitP12hom,*) Sol1(ApRow + 2 )
!~ 	write(OUnitP21hom,*) Sol1(ApRow + 3 )
!~ 	write(OUnitP22hom,*) Sol1(ApRow + 4 )
 
	call calcQ(C,NdimE,100)
	
	close(OUnitP11hom)
	close(OUnitP12hom)
	close(OUnitP21hom)
	close(OUnitP22hom)	
	close(OUnitChom)

end subroutine


subroutine calcQ (C,NdimE,NstepTheta)
	integer , intent(in) :: NdimE,NstepTheta
	real*8 , intent(in) :: C(NdimE,NdimE,NdimE,NdimE) 
	real*8 :: Q(NdimE,NdimE), normal(NdimE), theta, stepTheta , detQ
	integer :: i, j, k, l, n
	integer , parameter :: OUnitDetQ = 45, OUnitTheta = 46
	
	open (OUnitDetQ, file='DetQ.txt', Access = 'append')  
	open (OUnitTheta, file='Theta.txt', Access = 'append')

!~ 	stepTheta = PI/real(NstepTheta)
	stepTheta = 2.0d0*PI/real(NstepTheta)
	
!~ 	theta = -0.5d0*PI
	theta = 0.0d0
	do n = 1 , NstepTheta
		
		normal(1) = dcos(theta); normal(2) = dsin(theta)
!~ 		normal(1) = -dsin(theta); normal(2) = dcos(theta)
		
		
		Q = 0.0d0
		do i = 1 , NdimE
		do j = 1 , NdimE
		do k = 1 , NdimE
		do l = 1 , NdimE
			Q(i,k) = Q(i,k) + C(i,j,k,l)*normal(j)*normal(l) 
!~ 			Q(j,k) = Q(j,k) + C(i,j,k,l)*normal(i)*normal(l) 
		end do
		end do
		end do
		end do
		
		detQ = Q(1,1)*Q(2,2) - Q(1,2)*Q(2,1)
		
		write(OUnitDetQ,*) DetQ
		write(OUnitTheta,*) Theta
		
		theta = theta + stepTheta
	
	end do
	
	close(OUnitDetQ)	
	close(OUnitTheta)

end subroutine


end module
