!Last Modification : 11/03/2014
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Element for Hyperelasticity 
! Still in development
!     ------------------------------------------------------------------

!~ 
Subroutine posProcFibres(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : getF, NdimE, getMaterial, getDamagePar, maxMatPar, maxDamagePar, addToPKhomGen
	use fibresLib
	use ptsGaussLib, only : setNodG
	use fibresMod
	use DETERMINANT , only : lengthParam
	use genericConstitutiveLib
		
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
    
    integer :: i
	Integer , parameter :: iShiftField = 10,  NFieldPar = 2
	Real*8 :: MatPar(maxMatPar),  damagePar(maxDamagePar)
	Real*8 :: F(NdimE,NdimE), afib(NdimE), qfib(NdimE), Lfib
	Real*8 :: Sfib(NdimE), DSfib(NdimE,NdimE), Area, lfa
	real*8 , allocatable ::  SolUf(:) , Xel(:), SolUf0(:) !! all have dimension NodG*NdimE
	integer :: NodG, iMaterial, iDamageParType
	integer :: iShiftUf, constLaw , iFemType, iFtype ,  damageMethod, isAnalytic
	real*8 :: VField(5) , x0min, x0max , x0
	integer ::  NFields, CodeField(5) , PosField(5)
	
	matPar = 0.0d0

	iShiftUf = nint(CommonPar(1))
	iFemType = nint(commonPar(2))
	iFType = nint(commonPar(3))
	iMaterial = nint(commonPar(4))
	iDamageParType = nint(commonPar(5))
	damageMethod = nint(commonPar(6))
	isAnalytic = nint(commonPar(7))
	x0min = commonPar(8)
	x0max = commonPar(9)
	NFields = nint(CommonPar(10))
	
	do i = 1, NFields
		CodeField(i) = nint(CommonPar( iShiftField + (i-1)*NFieldPar + 1) )
		PosField(i)  = nint(CommonPar( iShiftField + (i-1)*NFieldPar + 2) )
	end do	
	
	call setNodG(iFemtype, NodG)	
	call getMaterial(constLaw, matPar, iMaterial) 
	call getDamagePar(damagePar, iDamageParType)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	call getSliceAllocate(SolUf,Sol1,1,NodG,iShiftUf + 1 ,iShiftUf + NdimE, iDofT)
	call getSliceAllocate(SolUf0,Sol0,1,NodG,iShiftUf + 1 ,iShiftUf + NdimE, iDofT)

	Area = Param(Ipos_Areaf)
	lfa = Param(Ipos_lfa)
	afib = Param(Ipos_af:Ipos_af+1)
	Lfib = Param(Ipos_Lf)
	
	matpar(3) = lfa
	matpar(5) = getLambdaFibreSimple(SolUf0,afib,Lfib,iFtype,NdimE)
	matpar(7) = Param(Ipos_stretch)
	
!~ 	if (damagePar(3)<0.0d0) then
		damagePar(3) = Param(Ipos_r0)
!~ 	end if
	
	call getF(F,iFtype)
	call setFibreKinematicSimple2(qfib,F,SolUf,afib,Lfib)

	x0 = 0.5*(Xel(1) + Xel(3))
	
	if( .not.( x0 > x0min .and. x0 < x0max) ) then
		damageMethod = 0
	end if
	
!~ 	damagePar(3) = 3.00001d0
	call SUPD_fib_1D_visc(Sfib,DSfib,qfib,damagePar,Param,matPar,DelT,constlaw,NDimE,damageMethod) 		
		
	if(Nfields < 1 ) then !! default mode
		Param(LengthParam + Ipos_stretch) = norm2(qfib)
		Param(LengthParam + Ipos_normStress) = dsign(1.0d0,dot_product(Sfib,afib))*norm2(Sfib)
		Param(LengthParam + Ipos_stress : LengthParam + Ipos_stress + 1 ) = Sfib 

		Param(Ipos_stretch) = norm2(qfib)
		Param(Ipos_normStress) = dsign(1.0d0,dot_product(Sfib,afib))*norm2(Sfib)
		Param(Ipos_stress : Ipos_stress + 1 ) = Sfib 
	else
		do i = 1 , NFields
			select case(CodeField(i))
				case(1) !! plot fibres strecht
					VField(i) = norm2(qfib)
				case(2) !! plot fibre tension
					VField(i) = dsign(1.0d0,dot_product(Sfib,afib))*norm2(Sfib) 
				case(3) !! Sx
					VField(i) = Sfib(1)
				case(4) !! Sy
					VField(i) = Sfib(2)
				case(5) !! plot fibres activation flag
					if(norm2(qfib) > lfa) then
						VField(i) = 1.0d0
					else
						VField(i) = 0.0d0
					end if
				case(6)
					VField(i) = getEnergyFibre(SolUf,Xel,MatPar,constLaw,iFtype,NdimE)
				case default
					write(0,*) "Code Field not found"
					VField(i) = 0
			end select

			Param(LengthParam + PosField(i) ) = VField(i) 
			
		end do
	end if
!~ 	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)	
	
End Subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine posProcFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow , i
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
!~ 
	Coupling=0 
!~ 	
	Coupling(1,1) = 1
	
	return
end Subroutine
