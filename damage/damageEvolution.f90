    !> Generic element for damage
    !!
	!! @param iShiftUf = nint(CommonPar(1))
	!! @param iFemType  = nint(CommonPar(2)) 
	!! @param iFType  = nint(commonPar(3))
	!! @param iMaterial  = nint(commonPar(4))
	!! @param iDamageParType  = nint(commonPar(5))
	!!
    !! @author Rocha, Felipe Figueredo
    
Subroutine damageEvolution(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	
	use funcAux
	use globalVariables, only : NdimE, getDamagePar, getMaterial, maxDamagePar, maxMatPar, getF
	use finiteStrainLib
	use damageNewLib
	use ptsGaussLib

	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
         
	Integer :: nG, i,j
	Real*8 :: GradUf(NdimE,NdimE) , Fmu(NdimE,NdimE)
	Real*8 ::energy, MatPar(maxMatPar),  damagePar(maxDamagePar)
	real*8 , allocatable ::  SolUf(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG, iShiftUf, constLaw, iFemType, iFtype, iMaterial, iDamageParType, iBubble, iSimplex, pOrder, NGP
	type(ptGaussClass) :: PtG

	iShiftUf = nint(commonPar(1))
	iFemType = nint(commonPar(2))
	iFtype = nint(commonPar(3))
	iMaterial = nint(commonPar(4))	
	iDamageParType = nint(commonPar(5))
	
	call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw,matPar,iMaterial)	
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)

	call getSliceAllocate(SolUf,Sol1,1,NodG,iShiftUf + 1 ,iShiftUf + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex) 		

	do nG = 1 , NGP
		call PtG%calcGradU(GradUf,SolUf,nG)	
		call getF(Fmu,iFtype)
		Fmu = Fmu + GradUf
		
		call strainEnergy(energy,Fmu,matPar,constLaw)
		
!~ 		call damageUpdate(damagePar,Param , energy , nG)	deprecated	
		
		call writeDamageFromParam(Param,energy,nG)
		
	end do
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)
	
	
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine damageEvolutionS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow , i
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)

	Coupling=0 
	
	Coupling(1,1) = 1
	
	return
end Subroutine
