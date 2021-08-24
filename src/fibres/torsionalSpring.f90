! ===============================================================================
Subroutine torsionalSpring(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, &
									Param, JParam, DelT, DTm, Time)
! ===============================================================================

	use funcAux , only : getSlice, setBlockToVector, setBlockDiagonalConstantToMatrix
	use globalVariables, only : getF, NdimE
	use fibresLib

	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)
!~ 	!   =====   END ARGUMENTS  =======
    
	! ======== COMMONPAR ARGUMENTS =====================
	integer iShiftUf !> @param  = nint(CommonPar(1))
	integer iShiftDeltaUf !> @var  = nint(CommonPar(2))
	integer iFType !> @var  = nint(commonPar(4))
	real*8 :: eta
	integer :: flagIncrementInside
	
	!====== OTHER VARIABLES ==============================
	integer , parameter :: NodG = 3 
	Real*8 :: SolU(NdimE*NodG), SolUf(NdimE*NodG), Xel(NdimE*NodG)
 	Real*8 :: Kstiff(NdimE*NodG,NdimE*NodG), Res(NdimE*NodG), F(NdimE,NdimE)
 	integer :: i, ip
 	
	iShiftUf = nint(CommonPar(1))
	iShiftDeltaUf = nint(CommonPar(2)) !!! now it's not being used
	iFType = nint(commonPar(3))
	eta = commonPar(4)
	flagIncrementInside = nint(commonPar(5))
	
	call getSlice(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	call getSlice(SolUf,Sol1,1,NodG,iShiftUf + 1 ,iShiftUf + NdimE, iDofT)

	call getF(F,iFtype) !!! with identity summed
	
	F = F - deltaKron2(1:NdimE,1:NdimE)

	do i = 1, NodG
		ip = (i-1)*NdimE
		SolU(ip+1:ip+NdimE) = matmul(F , Xel(ip+1:ip+NdimE) )
	end do
	
	SolU = SolU + SolUf
	
	call calc_Kstiff_and_Res_TorsionalSpring(Kstiff, Res, Xel, SolU, eta)

	if(iShiftUf == iShiftDeltaUf .and. flagIncrementInside == -1) then ! solver python updates outside
		Res = Res + matmul(Kstiff,SolUf)
	end if
	
	call setBlockToVector(BE,Res,NodG,iDofT, NdimE, iShiftDeltaUf ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	call setBlockToMatrixSquare(AE,Kstiff,NodG,iDofT, NdimE, iShiftDeltaUf) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
		
 	if(iShiftUf /= iShiftDeltaUf .and. flagIncrementInside > -1) then
		call setBlockDiagonalConstantToMatrix(AE,1.0d0,NodG, iDofT, NdimE, iShiftUf, iShiftUf)
		call setBlockDiagonalConstantToMatrix(AE,-1.0d0,NodG, iDofT, NdimE, iShiftUf, iShiftDeltaUf)
		call setBlockToVector(BE,SolUf,NodG,iDofT, NdimE, iShiftUf ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	end if
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%d%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine torsionalSpringS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux, only : setCoupling_pureDisplacement, setSimbolicBlockDiagonalToMatrix, numprint
	  
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    integer , parameter :: NdimE = 2, NodG = 3
    integer :: iShiftUf, iShiftDeltaUf, flagIncrementInside

	iShiftUf = nint(CommonPar(1))
	iShiftDeltaUf = nint(CommonPar(2)) !! not being used at the moment
	flagIncrementInside = nint(commonPar(5))
	
	call setCoupling_pureDisplacement(Coupling,NodG,idofT,iShiftDeltaUf,NdimE)

	if(iShiftUf /= iShiftDeltaUf .and. flagIncrementInside > -1) then		
		call setSimbolicBlockDiagonalToMatrix(Coupling,NodG, iDofT, NdimE, iShiftUf, iShiftUf)
		call setSimbolicBlockDiagonalToMatrix(Coupling,NodG, iDofT, NdimE, iShiftUf, iShiftDeltaUf)
	end if

end Subroutine



