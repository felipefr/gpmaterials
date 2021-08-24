module fibresMod

use globalVariables , only: NdimE 

implicit none

! New Param
!~ integer , parameter :: 	Ipos_flag1 = 1 , Ipos_flag2 = 2	, Ipos_Lf = 3, &
!~ 						Ipos_Areaf = 4, Ipos_Vf = 5, Ipos_lfa = 6 , Ipos_af = 7, & 
!~ 						Ipos_yrel1 = 9, Ipos_yrel2 = 11, &


public  initFibresMod
						
integer :: Ipos_flag1 , Ipos_flag2, Ipos_Lf, Ipos_Areaf, Ipos_Vf ,Ipos_lfa, Ipos_af, & 
		   Ipos_yrel1, Ipos_yrel2 , Ipos_Abar1 , Ipos_Abar2, &
		   Ipos_normal1, Ipos_normal2, Ipos_nConnectFibre1, Ipos_nConnectFibre2, &
		   Ipos_stress, Ipos_normStress, Ipos_stretch, Ipos_r0

real*8 , allocatable:: yG(:), Bten(:,:), Bten_invT(:,:)
real*8 :: measRVE=0.0d0, measFibres=0.0d0

contains 

subroutine initFibresMod()

	integer , parameter :: OUnitBten = 50
	
	allocate(Bten(NdimE,NdimE),Bten_invT(NdimE,NdimE),yG(NdimE))
	
	open (OUnitBten, file='Bten_yG.txt')
	read(OUnitBten,*) Bten
	read(OUnitBten,*) Bten_invT
	read(OUnitBten,*) measFibres
	read(OUnitBten,*) yG
	close(OUnitBten)  
	
	Ipos_flag1 = 1
	Ipos_flag2 = 2
	Ipos_Lf = 3				
	Ipos_Areaf = 4
	Ipos_Vf = 5
	Ipos_lfa = 6
	Ipos_af = 7 

	if(NdimE == 2) then
		Ipos_yrel1 = 9 
		Ipos_yrel2 = 11				
		Ipos_Abar1 = 13 
		Ipos_Abar2 = 14
		Ipos_normal1 = 15
		Ipos_normal2 = 17
		Ipos_nConnectFibre1 = 19
		Ipos_nConnectFibre2 = 20
		Ipos_r0 = 21
		Ipos_stress = 27
		Ipos_normStress = 29
		Ipos_stretch = 30


	else if(NdimE == 3) then
		Ipos_yrel1 = 10 
		Ipos_yrel2 = 13				
		Ipos_Abar1 = 16 
		Ipos_Abar2 = 17
		Ipos_normal1 = 18
		Ipos_normal2 = 21
		Ipos_nConnectFibre1 = 24
		Ipos_nConnectFibre2 = 25
		Ipos_r0 = 26
		Ipos_stress = 32
		Ipos_normStress = 35
		Ipos_stretch = 36

	else
		write(0,*) 'Wrong number of dimensions initFibresModVariables' , NdimE
	end if

end subroutine

end module

