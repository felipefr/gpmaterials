!     ------------------------------------------------------------------
Subroutine computeTangent(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use finiteStrainLib
	use damageNewLib
	use ptsGaussLib
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
         
	Integer :: nG, i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ipDim, ipDimE, ipDof,ipDofmod ! counters 
	Integer :: constLaw ! if NparamDamage = 0 , no damage
	Integer , parameter :: NdimE = 2
	Real*8 :: PhiA,PhiB, dPhiA(NdimE), dPhiB(NdimE), dv , Det , MatPar(9)
	Real*8 :: U(NdimE), GradU(NdimE,NdimE) , F(NdimE,NdimE)
	Real*8 :: DGradU(NdimE,NdimE), DFcan(NdimE,NdimE), Fcan(NdimE,NdimE) 
	Real*8 :: Ppk(NdimE,NdimE), D(NdimE,NdimE,NdimE,NdimE)
	real*8 , allocatable ::  SolU(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: pOrder , NGP, NodG, iFEMtype , iSimplex, iBubble
	integer :: iShiftFlucKL , iShiftU, kl, IdamageModify
	!! commonPar for the damage starts after the commonPar of the actual function
	integer , parameter :: iShiftDamageParam = 0 , ishiftDamageCommonPar = 9
	type(ptGaussClass) :: PtG
	
!~ 	!   =====   END ARGUMENTS  =======
	
	iShiftFlucKL = nint(CommonPar(1))
	iShiftU = nint(CommonPar(2))
	kl = nint(CommonPar(3))
	iFEMtype = nint(CommonPar(4))
	constLaw = nint(CommonPar(5))
	matPar(1:3) = commonPar(6:8)
	IdamageModify = nint(commonPar(9))  
	
	call getIJfromK(k,l,kl)
	
	Fcan = 0.0d0
	Fcan(k,l) = 1.0d0
			
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)

	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
		
	!--------------------------
	Do nG = 1, NGP ! LoopGauss
		call PtG%calcGradU(GradU,SolU,nG)

		F =  deltaKron(1:NdimE,1:NdimE) + GradU
		
		call calcPpk(Ppk,F,NdimE,matPar,constLaw)		
		call calcD(D,F,NdimE,matPar,constLaw)	
		
		call damageModifyTangentNew(D,Ppk,commonPar,Param, NdimE, &
								iShiftDamageCommonPar,iShiftDamageParam,nG,IdamageModify)
 		
 		dV=ptG%dV(nG)
 		
		DFcan = 0.0d0		
		call T4xT2(DFcan,D,Fcan)

		do A=1, NodG !LoopRow
			ApRow  = (A-1)*iDofT + iShiftFlucKL
			PhiA = PtG%Phi(A,nG)
			dPhiA= PtG%dPhi_G(:,A,nG)
			
			do p = 1 , NdimE
				App = ApRow+p
				
				BE(App) = BE(App) - dot_product(DFcan(p,:),dphiA)*DV  
				
				do B= 1 , NodG ! LoopCol ! not considering the simmetry
					BpCol  = (B-1)*iDofT + iShiftFlucKL
					PhiB = PtG%Phi(B,nG)
					dPhiB= PtG%dPhi_G(:,B,nG)
 					
					do q=1,NdimE
						Bqp = BpCol+q
						
						AE(App,Bqp) = AE(App,Bqp) + DGradPhiBdotGradPhiA(D,dPhiB,dPhiA,p,q,NdimE)*DV  
				
					end do ! loop Bq
				end do !LoopCol
			end do ! loop Ap
		end do !LoopRow
	end do !LoopGauss
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine computeTangentS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    use ptsGaussLib
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol, iShiftFlucKL, Femtype, NodG    
    integer , parameter :: NdimE = 2

	iShiftFlucKL = nint(commonPar(1))
    femType = nint(commonPar(4)) 
    
    call setNodG(femtype, NodG)
    
    Coupling = 0
    
    do A = 1,NodG
        ApRow = (A-1) * iDofT + iShiftFlucKL
		
		do B = 1,NodG
			BpCol = (B-1) * iDofT + iShiftFlucKL

			do p=1,nDimE
			do q=1,nDimE

				Coupling (ApRow + p, BpCol + q) = 1

			enddo
            enddo
		enddo
	enddo

end Subroutine
