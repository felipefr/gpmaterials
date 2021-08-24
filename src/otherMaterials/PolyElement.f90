	!> Element (Material) that assemble other elemental (Material) matrices using the same group of elements (Geometrical)
	!! @param CommonPar = ( a , a*[bi,ci] , a*[ci*dj] , dummy ) where : 
	!! @param NelemFamily = CommonPar(1) ==> a
	!! @param id_Element_Family(i) = CommonPar(see structure) ===> bi
	!! @param iSizeCommonPar(i) = CommonPar(see structure) ===> ci
	!! @param CommonPar (for each material = see structure) ===> dj 
	!! @author Rocha, Felipe Figueredo

Subroutine polyElement(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------

	use funcAux
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
    
    integer , parameter :: MaxNelemFamily = 5 , NauxParam = 2  	
    integer :: NelemFamily, id_Element_Family(MaxNelemFamily),  iSizeCommonPar(MaxNelemFamily), iShiftCommonPar  
    integer :: i, ip, ipp
    real*8 :: AEtemp(MaxLRows,MaxLRows), BEtemp(MaxLRows)
         	
	NelemFamily = nint(commonPar(1))
	
!~ 	call numprint(commonPar(1:18))
!~ 	
	do i = 1 , NelemFamily
		ip = (i-1)*NauxParam + 1
		id_Element_Family(i) = nint(commonPar(ip + 1))
		iSizeCommonPar(i) = nint(commonPar(ip + 2))
	end do
	
	iShiftCommonPar = 1 + NelemFamily*NauxParam
	
	do i = 1 , NelemFamily
		ip = iShiftCommonPar + 1
		ipp = ip + iSizeCommonPar(i) - 1
		AEtemp = 0.0d0
		BEtemp = 0.0d0
		call executerElement( id_Element_Family(i),AEtemp, BEtemp, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, &
								CommonPar(ip:ipp), Param, JParam, DelT, DTm, Time)
	
		AE = AE + AEtemp
		BE = BE + BEtemp
		
		iShiftCommonPar = iShiftCommonPar + iSizeCommonPar(i)
	end do
	
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine polyElementS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)

    integer , parameter :: MaxNelemFamily = 5 , NauxParam = 2  	
    integer :: NelemFamily, id_Element_Family(MaxNelemFamily),  &
				iSizeCommonPar(MaxNelemFamily), iShiftCommonPar , CouplingTemp(MaxLRows,MaxLRows)
    integer :: i, ip,ipp , j
    
    write(0,*) "Poly Element "
         	
	NelemFamily = nint(commonPar(1))

	do i = 1 , NelemFamily
		ip = (i-1)*NauxParam + 1
		id_Element_Family(i) = nint(commonPar(ip + 1))
		iSizeCommonPar(i) = nint(commonPar(ip + 2))
	end do
		
	iShiftCommonPar = 1 + NelemFamily*NauxParam
!~ 	
	do i = 1 , NelemFamily
		ip = iShiftCommonPar + 1
		ipp = ip + iSizeCommonPar(i) - 1
		
		CouplingTemp = 0
		call executerSymbolic(id_Element_Family(i),CouplingTemp,CommonPar(ip:ipp),iDofT,Ndim,MaxLRows,iAdd)
		Coupling = Coupling + CouplingTemp
		
		iShiftCommonPar = iShiftCommonPar + iSizeCommonPar(i)
	end do

	do i = 1 , MaxLRows
		do j = 1 , MaxLRows
			if(Coupling(i,j) > 1) Coupling(i,j) = 1 
		end do
	end do

end Subroutine

