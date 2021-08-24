module globalVariables

use loadingLib

implicit none

private Fsave, Ftemp, constLaw_save, matPar_save
public getF, allocateFsave, evolutionF, isAllocated_Fsave, &
	   setMaterial, getMaterial, isAllocatedMaterial, allocateMaterial, maxMatPar, addToTcoh, addToPKhom, addVol, addToPKhomGen

real*8, allocatable :: Fsave(:,:,:), Ftemp(:,:), Fold(:,:,:)
integer , allocatable :: constLaw_save(:)
real*8 , allocatable :: matPar_save(:,:)
integer , parameter :: maxMatPar = 9, maxDamagePar = 7
integer :: NdimE 
real*8, allocatable :: damagePar_save(:,:)
real*8 :: Vol, VolL
real*8 , allocatable :: PKhom(:,:), Tcoh(:)
real*8 , allocatable :: AnisoTensor(:,:), AnisoTensorInv(:,:), AnisoBarTensor(:,:), AnisoBarTensorInv(:,:)
real*8 :: VolFibTot
logical :: flagAnisoTensorInv

contains

subroutine contributeAsisotropyTensor(afib,Vfib)
	real*8 , intent(in) :: afib(:), Vfib
	integer :: i,j

	if(.not. allocated(AnisoTensor)) then
		allocate(AnisoTensor(NdimE,NdimE),AnisoTensorInv(NdimE,NdimE))
		AnisoTensor = 0.0d0
		VolFibTot = 0.0d0
		flagAnisoTensorInv = .false.
	end if

	if(flagAnisoTensorInv) return
	
	do i = 1 , NdimE
	do j = 1 , NdimE
		AnisoTensor(i,j) =  AnisoTensor(i,j) + Vfib*afib(i)*afib(j) 
	end do
	end do
	 
	 VolFibTot = VolFibTot + Vfib
	 
end subroutine

subroutine contributeVolFibTot(AreaFib,Lfib)
	real*8 , intent(in) :: AreaFib, Lfib

	if(flagAnisoTensorInv) return
	
	 VolFibTot = VolFibTot + AreaFib*Lfib
	 
end subroutine


subroutine contributeAsisotropyTensor_2(afib,y,AreaFib,Phi_if)
	real*8 , intent(in) :: afib(:), y(:), AreaFib, Phi_if
	integer :: i,j

	if(.not. allocated(AnisoTensor)) then
		allocate(AnisoTensor(NdimE,NdimE),AnisoTensorInv(NdimE,NdimE))
		AnisoTensor = 0.0d0
		VolFibTot = 0.0d0
		flagAnisoTensorInv = .false.
	end if

	if(flagAnisoTensorInv) return
	
	do i = 1 , NdimE
	do j = 1 , NdimE
		AnisoTensor(i,j) =  AnisoTensor(i,j) + Phi_if*AreaFib*y(i)*afib(j) 
	end do
	end do
	 
end subroutine

subroutine contributeAsisoTensors(afib,afibbar,Vfib)
	real*8 , intent(in) :: afib(:), afibbar(:), Vfib
	integer :: i,j

	if(.not. allocated(AnisoTensor)) then
		allocate(AnisoTensor(NdimE,NdimE),AnisoTensorInv(NdimE,NdimE))
		allocate(AnisoBarTensor(NdimE,NdimE),AnisoBarTensorInv(NdimE,NdimE))
		AnisoTensor = 0.0d0
		AnisoBarTensor = 0.0d0
		AnisoTensorInv = 0.0d0
		AnisoBarTensorInv = 0.0d0
		VolFibTot = 0.0d0
		flagAnisoTensorInv = .false.
	end if

	if(flagAnisoTensorInv) return
	
	do i = 1 , NdimE
	do j = 1 , NdimE
		AnisoTensor(i,j) =  AnisoTensor(i,j) + Vfib*afib(i)*afib(j)
		AnisoBarTensor(i,j) =  AnisoTensor(i,j) + Vfib*afib(i)*afibbar(j) 
	end do
	end do
	 
	 VolFibTot = VolFibTot + Vfib
	 
end subroutine


subroutine setAnisoTensorInv()
	real*8 :: dummy
	
	if(flagAnisoTensorInv) return
	
	AnisoTensor = AnisoTensor/VolFibTot
	
	write(0,*) "Aniso Tensor:"
	call numprint(AnisoTensor)
	
	call matinv(AnisoTensorInv,dummy, AnisoTensor)
	
	write(0,*) "Aniso Tensor Inv:"
	call numprint(AnisoTensorInv)
	
	flagAnisoTensorInv = .true.
	
end subroutine

subroutine  setAnisoTensors_inverses()
	real*8 :: dummy
	
	if(flagAnisoTensorInv) return
	
	AnisoTensor = AnisoTensor/VolFibTot
	
	call matinv(AnisoTensorInv,dummy, AnisoTensor)
	
	AnisoBarTensor = AnisoBarTensor/VolFibTot
	
	call matinv(AnisoBarTensorInv,dummy, AnisoBarTensor)
	
	flagAnisoTensorInv = .true.
end subroutine

subroutine getAnisoTensorInv(B)
	real*8 , intent(out):: B(:,:)
	
	if(flagAnisoTensorInv) then
		B = AnisoTensorInv
	else
		write(0,*) "AnisoTensorInv is being demanded without being computed"
		stop
	end if
	
end subroutine

subroutine getAnisoTensors_inverses(B,Bbar)
	real*8 , intent(out):: B(:,:), Bbar(:,:)
	
	if(flagAnisoTensorInv) then
		B = AnisoTensorInv
		Bbar = AnisoBarTensorInv
	else
		write(0,*) "AnisoTensorInv is being demanded without being computed"
		stop
	end if
	
end subroutine


subroutine exportDetQandTheta(detQ,theta,iElem,Nelem)
	real*8 , intent(in) :: detQ,theta 
	integer , intent(in) :: iElem, NElem
	integer , parameter :: OUnitTheta = 86, OUnitDetQ = 85 
	
	if(iElem == Nelem) then
		open (OUnitTheta, file='theta.txt', Access = 'append')
		open (OUnitDetQ, file='detQ.txt', Access = 'append')
	
		write(OUnitTheta,*) theta
		write(OUnitDetQ,*) detQ
		
		close(OUnitTheta)
		close(OUnitDetQ)
	end if	

end subroutine


subroutine addToTcoh(PKmu,normal,dVmu,iElem,IelemBegin,IelemEnd)
	real*8 , intent(in) :: dVmu, PKmu(:,:), normal(:)
	integer , intent(in) :: iElem, IelemBegin,IelemEnd
	integer , parameter :: OUnitTcoh = 85
	integer :: i
	
	if(iElem == IelemBegin) then
		if(.not. allocated(Tcoh)) then
			allocate(Tcoh(NdimE)) ! Temporary 
		end if
		VolL = 0.0
		Tcoh = 0.0
	end if
	
	Tcoh = Tcoh + dVmu*matmul(PKmu,normal)
	VolL = VolL + dVmu
	
	if(iElem == IelemEnd) then
		open (OUnitTcoh, file='Tcoh.txt', Access = 'append')
		
		Tcoh = Tcoh/VolL
	
		do i = 1, NdimE 
			write(OUnitTcoh,*) Tcoh(i)
		end do
				
		close(OUnitTcoh)
	end if	

end subroutine

subroutine addToPKhom(PKmu,dVmu,iElem,Nelem)
	real*8 , intent(in) :: dVmu, PKmu(:,:)
	integer , intent(in) :: iElem,Nelem

	call addToPKhomGen(PKmu,dVmu,iElem,1,Nelem)
	
end subroutine

subroutine addToPKhomGen(PKmu,dVmu,iElem,iElemBegin,Nelem)
	real*8 , intent(in) :: dVmu, PKmu(:,:)
	integer , intent(in) :: iElem,Nelem, iElemBegin

	if(iElem == IelemBegin) then
		call allocatePKhom()
	end if
	
	PKhom = PKhom + dVmu*PKmu
	Vol = Vol + dVmu
	
	if(iElem == Nelem) then
		call writePKhom()
	end if	

end subroutine

subroutine contributePKhom(PKmu,dVmu)
	real*8 , intent(in) :: dVmu, PKmu(:,:)
	
	PKhom = PKhom + dVmu*PKmu
	Vol = Vol + dVmu
	
end subroutine

subroutine allocatePKhom()
	if(.not. allocated(PKhom)) then
		allocate(PKhom(NdimE,NdimE)) ! Temporary 
		Vol = 0.0d0
		PKhom = 0.0d0
	end if
end subroutine 

subroutine writePKhom()
	integer,  parameter :: OUnitPKhom = 85
	integer :: i,j
	open (OUnitPKhom, file='PKhom.txt', Access = 'append')
		
	PKhom = PKhom/Vol
		
	!!! python convention ==> by lines (column are the fast index)
	do i = 1, NdimE
	do j = 1, NdimE
		write(OUnitPKhom,*) PKhom(i,j)
	end do 
	end do
	
	PKhom = 0.0d0
	Vol = 0.0d0
	 
	close(OUnitPKhom)

end subroutine

subroutine addVol(dVol,iElem,Nelem)
	real*8 , intent(in) :: dVol
	integer , intent(in) :: iElem,Nelem
	integer , parameter :: OUnitVol = 86
	
	if(iElem == 1) Vol = 0.0
	Vol = Vol + dVol
	
	if(iElem == Nelem) then
		open (OUnitVol, file='Vol.txt', Access = 'append')
		write(OUnitVol,*) Vol
		close(OUnitVol) 
	end if

end subroutine

subroutine useFold()
	Fsave = Fold
end subroutine

logical function isAllocated_Fsave() result(flag)
	flag = allocated(Fsave) 
end function

subroutine evolutionF(iFtype,LoadPar,Time,dt,LoadTypeProg,LoadType, NdimE, evolStyle)
	Real*8 , intent(in):: Time , dt, LoadPar(:)
	integer , intent(in) :: LoadType, LoadTypeProg, iFType, NdimE, evolStyle
	
	call setF(Ftemp,LoadPar,Time,dt,LoadTypeProg,LoadType)
		
	Fold = Fsave
	select case(evolStyle)
		case(0) !! absolute
			Fsave(:,:,iFtype) = deltaKron(1:NdimE,1:NdimE) + Ftemp
		case(1) !! by increment but using absolute functions
			Fsave(:,:,iFtype) = Fsave(:,:,iFtype) + Ftemp 
			call setF(Ftemp,LoadPar,Time-dt,dt,LoadTypeProg,LoadType)
			Fsave(:,:,iFtype) = Fsave(:,:,iFtype) - Ftemp
		case(2) !! by increment pure
!~ 			write(*,*) 'Fsave(:,:,', iFtype, ') changed'  
			Fsave(:,:,iFtype) = Fsave(:,:,iFtype) + Ftemp(:,:)
	end select
	
end subroutine

subroutine allocateFsave(NdimE,n)
	integer, intent(in) :: NdimE, n
	integer :: i
	
	if(allocated(Fsave)) then
		write(*,*) "trying to allocate Fsave already allocated"
		stop
	else
		allocate(Fsave(NdimE,NdimE,n), Ftemp(NdimE,NdimE), Fold(NdimE,NdimE,n))
		do i = 1,n
			Fsave(1:NdimE,1:NdimE,i) = deltaKron(1:NdimE,1:NdimE)
		end do
	end if
	
	Fold = Fsave
	
end subroutine

subroutine getF(F,iFtype)
	real*8 , intent(out) :: F(:,:)
	integer , intent(in) :: iFtype
	
	if((.not. allocated(Fsave)) .or. iFtype < 1) then
		F = deltaKron(1:NdimE,1:NdimE) 
	else if(iFtype <= size(Fsave,3)  ) then
!~ 		write(*,*) 'Fsave(:,:,', iFtype, ') got'
		F = Fsave(:,:,iFtype)
	else
		write(*,*) iFtype, "No F case is found"
		stop
	end if
	
end subroutine

subroutine getMaterial(constLaw,matPar,iMatType)
	integer , intent(out) :: constLaw
	real*8 , intent(out) :: matPar(:)
	integer , intent(in) :: iMatType
	integer n

	constLaw = constLaw_save(iMatType)
	matPar(:) = matPar_save(:,iMatType)
	
end subroutine

subroutine setMaterial(constLaw,matPar,iMatType)
	integer , intent(in) :: constLaw
	real*8 , intent(in) :: matPar(:)
	integer , intent(in) :: iMatType

	constLaw_save(iMatType) = constLaw
	matPar_save(:,iMatType) = matPar(:)
	
end subroutine

logical function isAllocatedMaterial() result(flag)
	flag = allocated(constLaw_save) 
end function

subroutine allocateMaterial(n)
	integer, intent(in) ::  n
	
	if(allocated(constLaw_save)) then
		write(*,*) "trying to allocate Material already allocated"
		stop
	else
		allocate(constLaw_save(n), matPar_save(maxMatPar,n))
	end if
	
end subroutine

subroutine getDamagePar(damagePar,idamageParType)
	real*8 , intent(out) :: damagePar(:)
	integer , intent(in) :: idamageParType
	integer n

	if(idamageParType>0) then
		damagePar(:) = damagePar_save(:,idamageParType)
	end if

end subroutine

subroutine setDamagePar(damagePar,iDamageParType)
	real*8 , intent(in) :: damagePar(:)
	integer , intent(in) :: iDamageParType

	damagePar_save(:,iDamageParType) = damagePar(:)
	
end subroutine

logical function isAllocatedDamagePar() result(flag)
	flag = allocated(damagePar_save) 
end function

subroutine allocateDamagePar(n)
	integer, intent(in) ::  n
	
	if(allocated(damagePar_save)) then
		write(*,*) "trying to allocate Damage already allocated"
		stop
	else
		allocate( damagePar_save(maxdamagePar,n))
	end if
	
end subroutine


end module
