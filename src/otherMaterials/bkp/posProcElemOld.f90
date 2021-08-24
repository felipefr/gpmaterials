	!> Element for pos processing variables and store in Param
	!! @param CommonPar = ( ... , a, a*[bi,ci]) where : 
	!! @param iShiftU = CommonPar(1)
	!! @param iFemType = CommonPar(2)
	!! @param iFtype = CommonPar(3)
	!! @param iMaterial = CommonPar(4) 
	!! @param iDamageParType = CommonPar(5)
	!! @param Nelem = CommonPar(6)
	!! @param NFields = CommonPar(7)  ===> a
	!! @param CodeField(i) = CommonPar(see structure) ==> bi
	!! @param PosField(i) = CommonPar(see structure) ==> ci	
	!! @author Rocha, Felipe Figueredo
	
 Subroutine posProcElemOld(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use finiteStrainLib !! New way
	use damageLib
	use multiscaleLib
	use globalVariables, only : getF, NdimE, getDamagePar, getMaterial, maxDamagePar, maxMatPar, &
								addVol, exportDetQandTheta, addToPKhom, addToTcoh
	use ptsGaussLib
	use DETERMINANT , only : iElem

	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodElt ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
         
	Integer :: nG,i,j,k,l,m,n,p,q, ip, A, B, ApRow, BpCol , App,Bqp, ipDim, ipDimE, ipDof,ipDofmod, NodG ! counters 
	Integer , parameter :: iShiftField = 7,  NFieldPar = 2
	Real*8 :: GradU(NdimE,NdimE) , F(NdimE,NdimE)
	Real*8 :: Ppk(NdimE,NdimE), D(NdimE,NdimE,NdimE,NdimE), PKmean(NdimE,NdimE), DMean(NdimE,NdimE,NdimE,NdimE)
	Real*8 :: MatPar(maxMatPar) ,  damagePar(maxDamagePar), VField(4) , energy, energyMean, dVng, dVtot
	real*8 , allocatable :: SolU(:) ,  Xel(:) 
	integer :: NGP,pOrder , IisPeriodic , iShiftU, constLaw , NFields, CodeField(5) , PosField(5) !! maximum 5 at now
    integer , parameter :: OUnitP11hom = 48, OUnitP12hom = 49, OUnitP21hom = 50
    integer , parameter :: OUnitP22hom = 51, OUnitDamage = 52, OUnitVonMises = 53
    integer , parameter :: OUnitDisplacement = 54, OUnitRd = 55, OUnitQd = 56
	integer :: iFemType,  iSimplex,  iBubble, iFtype, iMaterial, iDamageParType, Nelem 
	Real*8 :: normal(2), tangent(2) ,  theta , alpha, beta(2)
	type(ptGaussClass) :: PtG

	matPar = 0.0d0

	iShiftU= nint(commonPar(1))
	iFemType = nint(commonPar(2))
	iFtype = nint(commonPar(3))
	iMaterial = nint(commonPar(4))	 
	iDamageParType = nint(commonPar(5))
	Nelem = nint(commonPar(6))
	NFields = nint(CommonPar(7))
	do i = 1, NFields
		CodeField(i) = nint(CommonPar( iShiftField + (i-1)*NFieldPar + 1) )
		PosField(i)  = nint(CommonPar( iShiftField + (i-1)*NFieldPar + 2) )
	end do
	
	if(iDamageParType > 0) call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw, matPar, iMaterial)
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)
	
	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
	
	PKmean = 0.0d0
	DMean = 0.0d0
	energyMean = 0.0d0
	dVtot = 0.0d0
	
	do nG = 1 , NGP
		call PtG%calcGradU(GradU,SolU,nG)
		
		call getF(F,iFtype)
		
		F = F + GradU
		
		call calcPpk(Ppk,F,NdimE,matPar,constLaw)		
		call calcD(D,F,NdimE,matPar,constLaw)	
		if(iDamageParType > 0) call damageModifyTangent(D,Ppk,damagePar,Param,nG, 4)
		
		call strainEnergy(energy,F,matPar,constLaw)
		
		dVng=ptG%dV(nG)
		
		dVtot = dVtot + dVng
		PKmean = PKmean + dVng*Ppk
		DMean = DMean + dVng*D
		energyMean = energyMean + dVng*energy 	
	end do
	
	PKmean = PKmean/dVtot
	DMean = DMean/dVtot
	energyMean = energyMean/dVtot
	
	do i = 1 , NFields
		select case(CodeField(i))
			case(1) !! energy 
				VField(i) = energy
			case(2) !! vol
				VField(i) = dVtot
			case(3) !! Pnn				
				theta = datan(0.3d0/1.0d0)
				normal(1) = -dsin(theta)
				normal(2) = dcos(theta)				
				VField(i) = dot_product(matmul(PKmean,normal),normal)
			case(4) !! Pnt
				theta = datan(0.3d0/1.0d0)
				normal(1) = -dsin(theta)
				normal(2) = dcos(theta)
				tangent(1) = dcos(theta)
				tangent(2) = dsin(theta)
				VField(i) = dot_product(matmul(PKmean,normal),tangent)
			
			case(5) !! Pnb
				alpha = 0.5d0*PI
				theta = datan(0.3d0/1.0d0)
				normal(1) = -dsin(theta)
				normal(2) = dcos(theta)
				beta(1) = dcos(theta + alpha)
				beta(2) = dsin(theta + alpha)
				
				VField(i) = dot_product(matmul(PKmean,normal),beta)
!~ 				VField(i) = dot_product(matmul(PKmean,normal),normal)
				
			case(6) !! Pnb update (in case of fibers for example)
				alpha = 0.5d0*PI
				theta = datan(0.3d0/1.0d0)
				normal(1) = -dsin(theta)
				normal(2) = dcos(theta)
				beta(1) = dcos(theta + alpha)
				beta(2) = dsin(theta + alpha)
				
				VField(i) = Param(LengthParam + PosField(i) ) + dot_product(matmul(PKmean,normal),beta)
				
			case(7) !! Bifurcation Indicator
				alpha = 0.5d0*PI
				theta = datan(0.3d0/1.0d0)
				normal(1) = -dsin(theta)
				normal(2) = dcos(theta)
				beta(1) = dcos(theta + alpha)
				beta(2) = dsin(theta + alpha)
				
				VField(i) = 0.0d0
				do j = 1, NdimE
				do k = 1, NdimE
					VField(i) = VField(i) + GradU(j,k)*beta(j)*normal(k)
				end do
				end do
				write(*,*) VField(i)
!~ 				PAUSE
			
			case(8) !! energy accumulated
				VField(i) = energy + Param(LengthParam + PosField(i))
!~ 				write(0,*) 'LengthParam =' , LengthParam
!~ 				PAUSE 
				 
			case(11) !! P11
				VField(i) = PKmean(1,1)
			case(12) !! P12
				VField(i) = PKmean(1,2)
			case(13) !! P21
				VField(i) = PKmean(2,1)
			case(14) !! P22
				VField(i) = PKmean(2,2)

			case(15) !! T !! temporary
				ip = NodG*idofT + 18
				theta = Sol1(ip) 
				normal(1) = dcos(theta)
				normal(2) = dsin(theta)
				
				VField(i) = dot_product(matmul(PKmean,normal),normal)
				
			case(101) !! test put in global volume
				call addVol(dVtot,iElem,Nelem)
			
			case(102) !! test homogenize PKhom
				call addToPKhom(PKmean,dVtot,iElem,Nelem)
				
			case(103) !! test homogenize Tcoh

				ip = NodG*idofT + 18 !! theta
				theta = Sol1(ip) 
				normal(1) = dcos(theta)
				normal(2) = dsin(theta)
				!! using PosField as IelemBegin and Nelem as IelemEnd
				call addToTcoh(PKmean,normal,dVtot,iElem,PosField(i),Nelem) 

			case(104) !! export detQ and theta
				ip = NodG*idofT + 17 !! detQ
				call exportDetQandTheta(Sol1(ip),Sol1(ip+1),iElem,Nelem) 
				
			case default
				write(0,*) "Code Field not found"
				VField(i) = 0
		end select
		
		if(.not.(CodeField(i) == 101 .or. CodeField(i) == 102 .or. CodeField(i) == 103))  then  
			Param(LengthParam + PosField(i) ) = VField(i)
		end if

	end do
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)	
	
End Subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine posProcElemOldS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow , i
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)

	Coupling=0 
	
	Coupling(1,1) = 1
	
	return
end Subroutine
