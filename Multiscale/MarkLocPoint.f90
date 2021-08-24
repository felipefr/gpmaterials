Subroutine MarkLocPoint &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)

	use funcAux
	use finiteStrainLib !! New Way	
	use ptsGaussLib
	use damageNewLib !! just to use the lengthParam
	
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
    integer, parameter :: NdimE = 2
    integer :: i, j, ip, nG , ApRow, ipLag, ipDeltaU
    Real*8 ::  pen,lamb, eps, bValue
    logical :: locActive 
	Real*8 :: GradU(NdimE,NdimE)
	Real*8 , allocatable ::  SolU(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: iShiftDeltaU , iShiftMinDetQ, iShiftBifDomain, iShiftFlagLoc
	integer :: pOrder , NGP, NodG, iFEMtype , iSimplex, iBubble
	real*8 :: normal(NdimE), beta(NdimE), alpha, theta, bifFac
	real*8 , save :: minDetQ = 0.0d0 ,thetaCrit = 0.0d0
	type(ptGaussClass) :: PtG
	integer , save :: iElem = 0
	
	iShiftDeltaU = nint(commonPar(1))
	iShiftBifDomain = nint(CommonPar(2)) ! in the Param
	iShiftMinDetQ = nint(CommonPar(3)) ! in the last node
	iShiftFlagLoc = nint(CommonPar(4)) ! in geometrical nodes
	iFEMtype = nint(CommonPar(5))
	
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)

!~ 	if(Param(iShiftBifDomain + 1 )>0.0d0) then
!~ 		do i = 1 , NodG
!~ 			ip = (i-1)*iDofT + iShiftFlagLoc + 1
!~ 		
!~ 			AE(ip,ip) = 1.0d0
!~ 			BE(ip) = 1.0d0
!~ 		end do
!~ 		
!~ 		return
!~ 	end if
	
	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftDeltaU + 1 ,iShiftDeltaU + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	NGP = 1
	nG = 1
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
		
	call PtG%calcGradU(GradU,SolU,nG)

	bifFac = 0.0d0
	
	ip = NodG*idofT + iShiftMinDetQ + 1
	minDetQ = Sol1(ip)
	
	if(minDetQ < 0.0d0) then
		alpha = 0.5d0*PI 
		theta = datan(0.3d0/1.0d0)
		normal(1) = -dsin(theta)
		normal(2) = dcos(theta)
		beta(1) = dcos(theta + alpha)
		beta(2) = dsin(theta + alpha)
			
		do i = 1, NdimE
		do j = 1, NdimE
			BifFac = BifFac + GradU(i,j)*beta(i)*normal(j)
		end do
		end do
		
	end if
	
	Param(LengthParam + iShiftBifDomain + 1 ) = BifFac + Param(iShiftBifDomain + 1 )
	BifFac = Param(LengthParam + iShiftBifDomain + 1 )
!~ 	
!~ 	Param(LengthParam + iShiftBifDomain + 1 ) = BifFac 
			
	do i = 1 , NodG
		ip = (i-1)*iDofT + iShiftFlagLoc + 1
		
		AE(ip,ip) = 1.0d0
		
		if(BifFac > 0.0d0) then
			BE(ip) = 1.0d0
		end if

	end do

	deallocate(SolU,Xel)
	
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine MarkLocPointS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    use ptsGaussLib
	
    implicit none
    
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer ::  Coupling(MaxLRows,MaxLRows)
    Real*8 :: CommonPar(*), eps
    integer :: i, ip, iShiftFlagLoc, iFEMtype, NdimE, NodG
	
	iShiftFlagLoc = nint(CommonPar(4)) ! in geometrical nodes
	iFEMtype = nint(CommonPar(5))
	
	call setNodG(iFemType, NodG)
	
	do i = 1 , NodG
		ip = (i-1)*iDofT + iShiftFlagLoc + 1
		Coupling(ip,ip) = 1
	end do
	
end Subroutine

