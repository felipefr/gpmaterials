module multiscaleNewLib
use funcAux
use globalVariables , only : NdimE
use TIME_STEP_ITER
use DETERMINANT, only : iElem

implicit none

public updateBifurcationState ,  Dhom, PKhom , writeBifurcationState
private getQ, getQbar

real*8 , allocatable :: Dhom(:,:,:,:) , PKhom(:,:)
real*8 :: minDetQ, minDetQbar, thetaCrit, betaCrit, thetaBarCrit, betaBarCrit

real*8, allocatable :: detQ_vec(:), detQbar_vec(:), theta_vec(:), betaBar_vec(:)

contains 

subroutine initMultiscaleNewLib()
	allocate(Dhom(NdimE,NdimE,NdimE,NdimE) , PKhom(NdimE,NdimE))
end subroutine

subroutine getQ(Q,normal)
	real*8 , intent(out) :: Q(:,:)
	real*8 , intent(in) :: normal(:)
	integer :: i,j,k,l
	
	Q = 0.0d0
	do i = 1 , NdimE
	do j = 1 , NdimE
	do k = 1 , NdimE
	do l = 1 , NdimE
		Q(i,k) = Q(i,k) + Dhom(i,j,k,l)*normal(j)*normal(l) 
	end do
	end do
	end do
	end do

end subroutine


subroutine getQbar(Qbar,beta)
	real*8 , intent(out) :: Qbar(:,:)
	real*8 , intent(in) :: beta(:)
	integer :: i,j,k,l
	
	Qbar = 0.0d0
	do i = 1 , NdimE
	do j = 1 , NdimE
	do k = 1 , NdimE
	do l = 1 , NdimE
		Qbar(j,l) = Qbar(j,l) + Dhom(i,j,k,l)*beta(i)*beta(k) 
	end do
	end do
	end do
	end do

end subroutine


subroutine updateBifurcationState(nStepAngle,storeDetQ_vec)
	integer , intent(in) :: nStepAngle, storeDetQ_vec

	real*8 :: Q(NdimE,NdimE), Qbar(NdimE,NdimE), normal(NdimE), betaVec(NdimE), detQ, detQbar, angle, stepAngle
	integer :: n
	
	if((storeDetQ_vec == 1) .and. (.not. allocated(detQ_vec))) then
		allocate(detQ_vec(nStepAngle), detQbar_vec(nStepAngle), theta_vec(nStepAngle), betaBar_vec(nStepAngle))
	end if
	
	stepAngle = PI/real(NstepAngle)
	
	angle = -0.5*PI
	
	do n = 1 , NstepAngle
		
		normal(1) = dcos(angle)
		normal(2) = dsin(angle)		
		
		betaVec(1) = dcos(angle)
		betaVec(2) = dsin(angle)
		
		call getQ(Q,normal)
		call getQbar(Qbar,betaVec)
		
		call MatDet(detQ,Q) 
		call MatDet(detQbar,Qbar)
		
!~ 		call numprint(Q)
!~ 		call numprint(Qbar) 
		
!~ 		pause
			
		if(storeDetQ_vec == 1) then
			theta_vec(n) = angle
			betaBar_vec(n) = angle
			detQ_vec(n) = detQ
			detQbar_vec(n) = detQbar
		end if
		
		if(n==1 .or. detQ<minDetQ) then
			minDetQ = detQ
			thetaCrit = angle
			betaCrit = datan2(Q(1,1),-Q(1,2)) ! only valid for 2D
		end if
		
		if(n==1 .or. detQbar<minDetQbar) then
			minDetQbar = detQbar
			betaBarCrit = angle
			thetaBarCrit = datan2(Qbar(1,1),-Qbar(1,2)) ! only valid for 2D
		end if
				
		angle = angle + stepAngle
		
	end do

end subroutine

subroutine writeBifurcationState(storeDetQ_vec)
	integer , intent(in) :: storeDetQ_vec
	integer, parameter :: OUnit_Q = 15, OUnit_Qbar = 16, OUnit_detQ = 17, OUnit_detQbar=18 , OUnit_Dhom = 19
	
	if(Nstep == 0) then ! Nstep starts in zero
		open (OUnit_Q, file='bifState.txt')  
		open (OUnit_Qbar, file='bifStateBar.txt')  
		open (OUnit_Dhom, file='Dhom.txt')  
		
		if(storeDetQ_vec == 1) then
			open (OUnit_detQ, file='detQ.txt')
			open (OUnit_detQbar, file='detQbar.txt')
			write(OUnit_detQ,*) theta_vec
			write(OUnit_detQbar,*) betaBar_vec
		end if
	end if

	write(OUnit_Q,*) minDetQ , thetaCrit, betaCrit
	write(OUnit_Qbar,*) minDetQbar , thetaBarCrit, betaBarCrit
	write(OUnit_Dhom,*) transposeT4(Dhom)
	
	if(storeDetQ_vec == 1) then 
		write(OUnit_detQ,*) detQ_vec
		write(OUnit_detQbar,*) detQbar_vec
	end if
	
end subroutine

end module
