!     ------------------------------------------------------------------
Subroutine TangentHom(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use finiteStrainLib !! New Way	
	use globalVariables, only: getF, NdimE, getDamagePar, getMaterial, maxDamagePar, maxMatPar
	use damageNewLib
	use ptsGaussLib
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
         
	Integer :: nG, i,j,k,l, A, ApRow! counters 
	integer ::  ipDim, ipDimE, ipDofU,ipDofUmod, ipDofKL,ipDofKLmod , p, q
	Integer :: NparamElem, constLaw ! if NparamDamage = 0 , no damage
	Real*8 :: PhiA,PhiB, dPhiA(NdimE), dPhiB(NdimE), dv , Det , MatPar(maxMatPar),  damagePar(maxDamagePar)
	Real*8 :: GradUf(NdimE,NdimE) , Fmu(NdimE,NdimE)
	Real*8 :: Dhom(NdimE,NdimE,NdimE,NdimE) , DGradDeltaUkl(NdimE,NdimE), GradDeltaUkl(NdimE,NdimE)
	Real*8 :: D(NdimE,NdimE,NdimE,NdimE), Ppk(NdimE,NdimE)
	real*8 , allocatable ::  SolUf(:) , Xel(:), SolDeltaUkl(:,:,:) !! all have dimension NodG*NdimE
	integer :: iShiftDeltaU , iShiftUf , iShiftC, iAux, ipDofDeltaU, NElem,  IdamageModify
	integer :: pOrder , NGP, NodG, iFEMtype , iFtype, iSimplex, iBubble, iMaterial ,iDamageParType
	type(ptGaussClass) :: PtG
	
	iShiftC = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iShiftUf = nint(CommonPar(3))
	iFEMtype = nint(CommonPar(4))
	iFtype = nint(CommonPar(5))
	iMaterial = nint(CommonPar(6))
	iDamageParType = nint(commonPar(7))
	NElem = nint(commonPar(8))
	
	call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw, matPar, iMaterial)
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)

	call getSliceAllocate(SolUf,Sol1,1,NodG,iShiftUf + 1 ,iShiftUf + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	allocate(SolDeltaUkl(NdimE*NodG,NdimE,NdimE))
	
	do k = 1, NdimE
	do l = 1, NdimE
		call getSlice(SolDeltaUkl(:,k,l),Sol1,1,NodG,iShiftDeltaU + iShiftKL(k,l,NdimE)*NdimE + 1 , &
		             iShiftDeltaU + iShiftKL(k,l,NdimE)*NdimE +  NdimE, iDofT)
	end do
	end do
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
		
	!--------------------------
	Do nG = 1, NGP ! LoopGauss
		call PtG%calcGradU(GradUf,SolUf,nG)
		
		call getF(Fmu,iFtype)
		Fmu =  Fmu + GradUf
		
		call calcPpk(Ppk,Fmu,NdimE,matPar,constLaw)		
		call calcD(D,Fmu,NdimE,matPar,constLaw)	
		
!~ 		call damageModifyTangent(D,Ppk,damagePar,Param, nG, 3) # deprecated
		
	    dV=ptG%dV(nG)
						
		Dhom = D
 		
		do k = 1, NdimE
		do l = 1, NdimE
			call PtG%calcGradU(GradDeltaUkl,SolDeltaUkl(:,k,l),nG)
			call T4xT2(DGradDeltaUKl,D,GradDeltaUkl)
			Dhom(:,:,k,l) = Dhom(:,:,k,l) + DGradDeltaUkl  
		end do
		end do 
		
		A = NodG + 1
		ApRow  = (A-1)*iDofT + iShiftC
			
		do i = 1 , NdimE
		do j = 1 , NdimE
		do k = 1 , NdimE
		do l = 1 , NdimE
			iAux = ApRow + iShiftIJKL(i,j,k,l,NdimE) + 1
			BE(iAux) = Dhom(i,j,k,l)*dV
			AE(iAux,iAux) = 1.0d0/real(NElem) 
		end do
		end do
		end do
		end do
!~ 	
!~ 		BE(ApRow + 1) = Ppk(1,1)*dV
!~ 		BE(ApRow + 2) = Ppk(1,2)*dV
!~ 		BE(ApRow + 3) = Ppk(2,1)*dV
!~ 		BE(ApRow + 4) = Ppk(2,2)*dV
	
	end do !LoopGauss

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine TangentHomS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux
	use ptsGaussLib
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A, ApRow, iShiftC,  i, j, k , l, iAux, femType, NodG
	integer, parameter :: NdimE = 2

    Coupling = 0
    
	iShiftC = nint(commonPar(1))
	femType = nint(commonPar(4)) 
       
	call setNodG(femtype, NodG)
    
	A = NodG + 1
	ApRow  = (A-1)*iDofT + iShiftC
		
	do i = 1 , NdimE
	do j = 1 , NdimE
	do k = 1 , NdimE
	do l = 1 , NdimE
		iAux = ApRow + iShiftIJKL(i,j,k,l,NdimE) + 1
		Coupling(iAux,iAux) = 1
	end do
	end do
	end do
	end do
		
!~ 	write(*,*) "TLflucTangHom Symbolic"
!~ 	call IprintMat(Coupling)
end Subroutine
