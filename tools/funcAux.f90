module funcAux
implicit none

private printVec , printMat, IprintVec, IprintMat, printTen
public AddEye, MatInv3, MatInv2, MatInv, calcI1, calcI2 , calcI3 , dot_product2, &
crossProduct,  crossProduct2, T4xT2, calc_a0, calc_he, calc_he_tri, calc_he_tri_max, & 
init_seed, PI, deltaKron, deltaKron2, VotimesW,  getIJfromK, iShiftKL, iShiftIJKL, numPrint, getSlice, getSliceAllocate, &
getT2fromV, getT4fromV,  getVfromT2, getVfromT4, &
setAEandBEfib, setAEandBE_lagrangeMultipliers, setAEandBE_LM, setAEandBE_LM_full, &
 setCoupling_lagrangeMultipliers, setCoupling_pureDisplacement, & 
isOnBoundary, getNormalBoundary, norm2, setBlockToMatrix, setBlockToVector, setBlockToMatrixSquare, &
setBlockDiagonalConstantToMatrix, setSimbolicBlockToMatrixSquare, setSimbolicBlockToMatrix, setSimbolicBlockDiagonalToMatrix, &
voigtTen4toTen2, voigtTen2toVec, MatDet, getNormal, solveQuadraticEquation, &
VotimesW_fun, VodotW_fun, VodotW, transposeT4, FindKeyword

interface numPrint !! overloading of subroutines
	module procedure :: IprintMat
	module procedure :: IprintVec
	module procedure :: printMat
	module procedure :: printVec
	module procedure :: printTen
end interface			

interface VotimesW !! overloading of subroutines
	module procedure :: VotimesW_vec_vec
	module procedure :: VotimesW_mat_mat
end interface			

interface VotimesW_fun !! overloading of function
	module procedure :: VotimesW_vec_vec_fun
	module procedure :: VotimesW_mat_mat_fun
end interface			


private
integer :: I,J
real*8 , dimension(3,3) :: deltaKron
real*8 , dimension(2,2) :: deltaKron2
DATA((deltaKron(I,J),I=1,3),J=1,3) /1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0/
DATA((deltaKron2(I,J),I=1,2),J=1,2) /1.0d0,0.0d0,0.0d0,1.0d0/


real*8, parameter :: PI=4.D0*DATAN(1.D0)

contains 

subroutine getT2fromV(T,V)
	real*8 , intent(in) :: V(:)
	real*8 , intent(out) :: T(:,:)
	integer :: i,j , ip, n
	 
	n = size(T,1)

	do j = 1 , n
	do i = 1 , n
		ip = iShiftKL(i,j,n)
		T(i,j) = V(ip + 1)  
	end do
	end do

end subroutine


subroutine getT4fromV(T,V)
	real*8 , intent(in) :: V(:)
	real*8 , intent(out) :: T(:,:,:,:)
	integer :: i,j,k,l , ip, n
	 
	n = size(T,1)
	
	do l = 1 , n
	do k = 1 , n
	do j = 1 , n
	do i = 1 , n
		ip = iShiftIJKL(i,j,k,l,n)
		T(i,j,k,l) = V(ip + 1)   
	end do
	end do
	end do
	end do

end subroutine 



subroutine getVfromT2(V,T)
	real*8 , intent(out) :: V(:)
	real*8 , intent(in) :: T(:,:)
	integer :: i,j , ip, n
	 
	n = size(T,1)

	do j = 1 , n
	do i = 1 , n
		ip = iShiftKL(i,j,n)
		V(ip + 1) = T(i,j)  
	end do
	end do

end subroutine


subroutine getVfromT4(V,T)
	real*8 , intent(out) :: V(:)
	real*8 , intent(in) :: T(:,:,:,:)
	integer :: i,j,k,l , ip, n
	 
	n = size(T,1)
	
	do l = 1 , n
	do k = 1 , n
	do j = 1 , n
	do i = 1 , n
		ip = iShiftIJKL(i,j,k,l,n)
		V(ip + 1) = T(i,j,k,l)  
	end do
	end do
	end do
	end do

end subroutine 

subroutine getSlice(Vout,Vin,row0,row1,col0,col1,Ncol)
	real*8 , intent(out) :: Vout(:)
	real*8 , intent(in) :: Vin(*)
	integer , intent(in) :: row0, row1, col0, col1, Ncol
	integer :: i,j,ip, k
	
	k = 0
	do i = row0, row1
		ip = (i-1)*Ncol
		do j = col0, col1
			k = k + 1
			Vout(k) = Vin(ip + j) 
		end do
	end do
	
end subroutine

subroutine getSliceAllocate(Vout,Vin,row0,row1,col0,col1,Ncol)
	real*8 , allocatable, intent(out) :: Vout(:)
	real*8 , intent(in) :: Vin(*)
	integer , intent(in) :: row0, row1, col0, col1, Ncol
	integer :: i,j,ip, k, alocErr

	allocate(Vout((row1-row0 + 1)*(col1-col0 + 1)), stat = alocErr)
	
!~ 	write(*,*) "alocErr = " , alocErr
	
	k = 0
	do i = row0, row1
		ip = (i-1)*Ncol
		do j = col0, col1
			k = k + 1
			Vout(k) = Vin(ip + j) 
		end do
	end do
	
end subroutine

integer function iShiftIJKL(i,j,k,l,NdimE)
	integer, intent(in) :: i,j,k,l,NdimE
	iShiftIJKL = ( ( (i-1)*NdimE + j - 1 )*NdimE + k - 1 )*NdimE + l - 1
end function

integer function iShiftKL(k,l,NdimE) 
	integer, intent(in) :: k,l,NdimE
	iShiftKL = (k-1)*NdimE + l - 1 
end function

subroutine getIJfromK(i,j,k)
	integer , intent(out) :: i,j
	integer , intent(in) :: k
	
	if(k==1) then
		i=1; j=1
	else if(k==2) then
		i=1; j=2
	else if(k==3) then
		i=2; j=1
	else if(k==4) then
		i=2; j=2
	end if
	
end subroutine

subroutine calc_he(XLL,he) !! just for tetrahedra
	integer,parameter :: NodElt = 4 , Ndim = 3
	real*8 , intent(in) :: XLL(NodElt*Ndim)
	real*8 , intent(out) :: he 
	real*8 :: edgeVec(Ndim)
	
	edgeVec = XLL(1:3) - XLL(4:6) !! 1 - 2
	he = dsqrt(dot_product(edgeVec,edgeVec))
	edgeVec = XLL(1:3) - XLL(7:9) !! 1 - 3
	he = he + dsqrt(dot_product(edgeVec,edgeVec))
	edgeVec = XLL(1:3) - XLL(10:12) !! 1 - 4
	he = he + dsqrt(dot_product(edgeVec,edgeVec))
	edgeVec = XLL(4:6) - XLL(7:9) !! 2 - 3
	he = he + dsqrt(dot_product(edgeVec,edgeVec))
	edgeVec = XLL(4:6) - XLL(10:12) !! 2 - 4
	he = he + dsqrt(dot_product(edgeVec,edgeVec))
	edgeVec = XLL(7:9) - XLL(10:12) !! 3 - 4
	he = he + dsqrt(dot_product(edgeVec,edgeVec))
	
	he = he/6.0d0
	
end subroutine

subroutine calc_he_tri(he,XLL) !! just for triangles
	integer,parameter :: NodElt = 3 , Ndim = 2
	real*8 , intent(in) :: XLL(NodElt*Ndim)
	real*8 , intent(out) :: he 
	real*8 :: V(Ndim),U(Ndim), area
	
	U = XLL(3:4) - XLL(1:2) !! 2 - 1
	V = XLL(5:6) - XLL(1:2) !! 3 - 1
	
	call crossProduct2(area,U,V)
	area = 0.5d0*dabs(area)
	
!~ 	he = dsqrt(area)
	he = dsqrt(2.0d0*area)
	
end subroutine

subroutine calc_he_tri_max(he,XLL) !! just for triangles
	integer,parameter :: NodElt = 3 , Ndim = 2
	real*8 , intent(in) :: XLL(NodElt*Ndim)
	real*8 , intent(out) :: he 
	real*8 :: V(Ndim), area, maxHe2, HeV2
	integer :: i , j, ip,jp
	
	maxHe2 = 0.0d0
	
	do i = 1,NodElt
		j = mod(i,NodElt) + 1
		ip = (i-1)*Ndim
		jp = (j-1)*Ndim
		 
		V = XLL(ip + 1 : ip + Ndim) - XLL(jp + 1 : jp + Ndim) 
		
		HeV2 = dot_product(V,V)
		if(HeV2 > maxHe2) maxHe2 = HeV2		
	end do
	
	he = dsqrt(maxHe2)
	
end subroutine

subroutine calc_a0(a0,XLL,betaAngle,Ndim,NodElt)
	integer, intent(in) :: Ndim,NodElt
	real*8 , intent(in) :: XLL(NodElt*Ndim), betaAngle
	real*8, intent(out) :: a0(3)
	real*8 ::  G(Ndim), sinTheta, cosTheta, rAux ! angle from the baricentric
	integer :: i,j,ip
	
	G(:)=0.0d0
	do i = 1,NodElt !!!! Maybe is just until 4
		ip = (i-1)*Ndim
		do j = 1 , Ndim
			G(j) = G(j) + XLL(ip+j) 
		end do
	end do
	
	G=G/real(NodElt)
	
	rAux = dsqrt(G(1)*G(1) + G(2)*G(2))
	sinTheta = G(1)/rAux
	cosTheta = G(2)/rAux
	
	rAux = PI*betaAngle/180.0d0 
	
	a0(1) = - dcos(rAux) * sinTheta
	a0(2) = dcos(rAux) * cosTheta
	a0(3) = dsin(rAux)
	
end subroutine

subroutine calcI1(I1,M)
	implicit none
	real*8 , intent(out) :: I1
	real*8 , intent(in) :: M(:,:)
	integer n
	n = size(M,1)
	
	if(n==3) then 
		I1 = M(1,1) + M(2,2) + M(3,3)
	else if(n==2) then
		I1 = M(1,1) + M(2,2) + 1.0d0
	end if
	
	return
end subroutine

subroutine calcI2(I2,M)
	implicit none
	real*8 , intent(out) :: I2
	real*8 , intent(in) :: M(:,:)
	integer n
	n = size(M,1)
	
	if(n==3) then 
		I2 = M(1,1)*M(2,2) + M(1,1)*M(3,3) + M(2,2)*M(3,3) - M(1,2)*M(2,1) - M(1,3)*M(3,1) - M(2,3)*M(3,2)
	else if(n==2) then
		I2 = M(1,1)*M(2,2) + M(1,1) + M(2,2) - M(1,2)*M(2,1)
	end if
	
	return
end subroutine

subroutine calcI3(I3,M)
	implicit none
	real*8 , intent(out) :: I3
	real*8 , intent(in) :: M(:,:)
	integer n
	n = size(M,1)
	
	if(n==3) then 
		I3 = M(1,1)*M(2,2)*M(3,3) + M(1,2)*M(2,3)*M(3,1) + M(1,3)*M(2,1)*M(3,2) &
			- M(1,1)*M(2,3)*M(3,2) - M(1,2)*M(2,1)*M(3,3) - M(1,3)*M(2,2)*M(3,1)
	else if(n==2) then
		I3 = M(1,1)*M(2,2) - M(1,2)*M(2,1) 
	end if
	
	return
end subroutine

real*8 function dot_product2(A,B) result(dot)
	real*8, intent(in) :: A(:,:) , B(:,:) 
	integer :: n , i, j
	n = size(A,1)
	
	dot = 0.0d0
	do i = 1,n
		do j = 1,n
			dot = dot + A(i,j)*B(i,j)
		end do
	end do
	
end function

subroutine addEye(M)
	implicit none
	real*8, intent(inout) :: M(:,:)
	integer :: i,n
	n = size(M,1)
	
	do i=1,n
		M(i,i) = M(i,i) + 1.0d0
	end do
	
	return
end subroutine

SubRoutine T4xT2(S,C,T) ! S = C : T
!~  Implicit none
 integer NDim
 real*8 ::  C(:,:,:,:), S(:,:) , T(:,:)
 integer :: i,j,k,l
	
 Ndim = size(S,1) !! all the same dimension (2 or 3)
	
	S = 0.0d0 
	do i = 1 , NDim
	do j = 1 , NDim
	do k = 1 , NDim
	do l = 1 , NDim
			S(i,j) = S(i,j) + C(i,j,k,l)*T(k,l)
	end do
	end do
	end do
	end do
	  
Return
End subroutine

subroutine crossProduct(w,normW,u,v) ! w = (u x v)/|u x v|
	real*8 , intent(out) :: w(3), normW
	real*8 , intent(in) :: u(3),v(3)
	w(1) = u(2)*v(3) - v(2)*u(3)
	w(2) = u(3)*v(1) - v(3)*u(1)
	w(3) = u(1)*v(2) - v(1)*u(2)
	
	normW = dsqrt(dot_product(w,w)) 
	w = w/normW
	
end subroutine

subroutine crossProduct2(w,u,v) 
	real*8 , intent(out) :: w 
	real*8 , intent(in) :: u(2),v(2)
	w = u(1)*v(2) - v(1)*u(2)
	
end subroutine


function VotimesW_vec_vec_fun(v,w) result(M)
	real*8, intent(in) :: v(:), w(:)
	real*8, allocatable :: M(:,:)
	integer :: n1, n2
	
	n1 = size(v,1)
	n2 = size(w,1)
	
	allocate(M(n1,n2))

	call VotimesW_vec_vec(M,v,w)
	
end function 

function VotimesW_mat_mat_fun(V,W) result(M)
	real*8, intent(in) :: V(:,:), W(:,:)
	real*8, allocatable :: M(:,:,:,:)
	integer :: n1, n2, n3, n4
	
	n1 = size(V,1)
	n2 = size(V,2)
	n3 = size(W,1)
	n4 = size(W,2)
	
	allocate(M(n1,n2,n3,n4))

	call VotimesW_mat_mat(M,V,W)
	
end function 

subroutine VotimesW_vec_vec(M,v,w)
	implicit none
	real*8, intent(in) :: v(:), w(:)
	real*8, intent(out) :: M(:,:)
	integer :: n1, n2 , j, i
	
	n1 = size(v,1)
	n2 = size(w,1)
	
	M = 0.0d0
	
	do i = 1,n1
		do j = 1,n2
			M(i,j) = v(i)*w(j)
		end do
	end do

end subroutine


subroutine VotimesW_mat_mat(M,V,W)
	implicit none
	real*8, intent(in) :: V(:,:), W(:,:)
	real*8, intent(out) :: M(:,:,:,:)
	integer :: n1, n2 , n3, n4, i, j, k, l
	
	n1 = size(V,1)
	n2 = size(V,2)
	n3 = size(W,1)
	n4 = size(W,2)
		
	M = 0.0d0
	
	do i = 1,n1
	do j = 1,n2
	do k = 1,n3
	do l = 1,n4
		M(i,j,k,l) = V(i,j)*W(k,l)
	end do
	end do
	end do
	end do

end subroutine

function VodotW_fun(V,W) result(M) ! always with matrices
	real*8, intent(in) :: V(:,:), W(:,:)
	real*8, allocatable :: M(:,:,:,:)
	integer :: n1, n2, n3, n4
	
	n1 = size(V,1)
	n2 = size(W,1)
	n3 = size(V,2)
	n4 = size(W,2)
	
	allocate(M(n1,n2,n3,n4))

	call VodotW(M,V,W)
	
end function 

subroutine VodotW(M,V,W) ! always with matrices
	implicit none
	real*8, intent(in) :: V(:,:), W(:,:)
	real*8, intent(out) :: M(:,:,:,:)
	integer :: n1, n2 , n3, n4, i, j, k, l
	
	n1 = size(V,1)
	n2 = size(W,1)
	n3 = size(V,2)
	n4 = size(W,2)
		
	M = 0.0d0
	
	do i = 1,n1
	do j = 1,n2
	do k = 1,n3
	do l = 1,n4
		M(i,j,k,l) = V(i,k)*W(j,l)
	end do
	end do
	end do
	end do

end subroutine

Subroutine MatDet(detM,M)
	implicit none
	real*8, intent(in) :: M(:,:)
	real*8, intent(out) :: detM
	real*8 :: A1, A2, A3, B1, B2, B3, C1, C2, C3 
	integer :: n 
	
	n = size(M,1)
	
	if (n == 2) then
		detM= M(1,1)*M(2,2) - M(1,2)*M(2,1) 
	else if(n==3) then
		A1 = M(3,3)*M(2,2) - M(3,2)*M(2,3)
		B1 = M(3,3)*M(2,1) - M(3,1)*M(2,3)
		C1 = M(3,2)*M(2,1) - M(3,1)*M(2,2)

		A2 = M(3,3)*M(1,2) - M(3,2)*M(1,3)
		B2 = M(3,3)*M(1,1) - M(3,1)*M(1,3)
		C2 = M(3,2)*M(1,1) - M(3,1)*M(1,2)

		A3 = M(2,3)*M(1,2) - M(2,2)*M(1,3)
		B3 = M(2,3)*M(1,1) - M(2,1)*M(1,3)
		C3 = M(2,2)*M(1,1) - M(2,1)*M(1,2)

		detM= M(1,1)*A1 - M(1,2)*B1 + M(1,3)*C1 
	else 
		write(0,*) 'WARNING, compute determinant is not possible (n=2 or n=3)'
	end if
		  
end subroutine

subroutine MatInv(Minv,detM,M)
	implicit none
	real*8, intent(in) :: M(:,:)
	real*8, intent(out) :: Minv(:,:), detM
	integer :: n
	
	n = size(M,1)
	
	if(n==3) then
		call MatInv3(Minv,detM,M)
	else if(n==2) then
		call MatInv2(Minv,detM,M)
	end if
	
end subroutine

Subroutine MatInv2(Minv,detM,M)
	implicit none
	real*8, intent(in) :: M(2,2)
	real*8, intent(out) :: Minv(2,2), detM
	real*8 , parameter :: tol = 10.0e-10

	detM= M(1,1)*M(2,2) - M(1,2)*M(2,1) 

	if (dabs(detM) > tol ) then
		Minv(1,1)    = M(2,2) / detM
		Minv(1,2)    =-M(1,2) / detM
		Minv(2,1)    =-M(2,1) / detM
		Minv(2,2)    = M(1,1) / detM
	else
		write(0,*) 'WARNING, Matrix Determinat = 0 in MatInv2() Felipe'
	end if
		  
end subroutine

 Subroutine MatInv3(Minv,detM,M)
	implicit none
	real*8, intent(in) :: M(3,3)
	real*8, intent(out) :: Minv(3,3), detM
	real*8 :: A1,B1,C1,A2,B2,C2,A3,B3,C3 ! minors
	real*8 , parameter :: tol = 10.0e-10

	A1 = M(3,3)*M(2,2) - M(3,2)*M(2,3)
	B1 = M(3,3)*M(2,1) - M(3,1)*M(2,3)
	C1 = M(3,2)*M(2,1) - M(3,1)*M(2,2)

	A2 = M(3,3)*M(1,2) - M(3,2)*M(1,3)
	B2 = M(3,3)*M(1,1) - M(3,1)*M(1,3)
	C2 = M(3,2)*M(1,1) - M(3,1)*M(1,2)

	A3 = M(2,3)*M(1,2) - M(2,2)*M(1,3)
	B3 = M(2,3)*M(1,1) - M(2,1)*M(1,3)
	C3 = M(2,2)*M(1,1) - M(2,1)*M(1,2)

	detM= M(1,1)*A1 - M(1,2)*B1 + M(1,3)*C1 

	if (dabs(detM) > tol ) then
		Minv(1,1)    = A1 / detM
		Minv(1,2)    =-A2 / detM
		Minv(1,3)    = A3 / detM
		Minv(2,1)    =-B1 / detM
		Minv(2,2)    = B2 / detM
		Minv(2,3)    =-B3 / detM
		Minv(3,1)    = C1 / detM
		Minv(3,2)    =-C2 / detM
		Minv(3,3)    = C3 / detM
	else
		write(0,*) 'WARNING, Matrix Determinat = 0 in MatInv3() Felipe'
!~ 		call stop_error('WARNING, Matrix Determinat = 0 in MatInv3()')
	end if
		 
end subroutine
 
subroutine printVec(v)
	real*8 :: v(:)
	integer n , i
	
	n = size(v)
	write(0,*) "=============="
	do i = 1, n
		write(0,"(e15.9)") v(i)
	end do
	write(0,*) "=============="
	
end subroutine

subroutine IprintVec(v)
	integer :: v(:)
	integer n , i
	
	n = size(v)
	write(0,*) "=============="
	do i = 1, n
		write(0,"(g15.5)") v(i)
	end do
	write(0,*) "=============="
	
end subroutine
!~ 
subroutine printMat(A)
	real*8 :: A(:,:)
	integer n , m, i , j
	
	n = size(A,1)
	m = size(A,2)
	
	write(0,*) "============================="
	do i = 1, n
		write(0,"(100(e15.6,1X))") ( A(i,j) , j = 1 , m)
!~ 		write(0,"(100(e15.5,2X))") ( A(i,j) , j = 1 , m)
	end do
	write(0,*) "============================="
	
end subroutine

subroutine IprintMat(A)
	integer :: A(:,:)
	integer n , m, i , j
	
	n = size(A,1)
	m = size(A,2)
	
	write(0,*) "=============================="
	do i = 1, n
		write(0,"(100I3)") ( A(i,j) , j = 1 , m)
!~ 		write(0,"(100g10.5)") ( A(i,j) , j = 1 , m)
	end do
	write(0,*) "=============================="
	
end subroutine

subroutine printTen(C)
	real*8 :: C(:,:,:,:)
	integer :: n , i , j , k , l 
	
	n = size(C,1)
	write(0,*) "=============================="
	
	do i = 1, n
	do j = 1, n
	do k = 1, n
	do l = 1, n
		write(0,"(4(I1,1X),'=',e15.5)") i,j,k,l, C(i,j,k,l)
	end do
	end do
	end do
	end do
	write(0,*) "=============================="
	
end subroutine

Subroutine init_seed()
	implicit none
	Integer i , size_seed , clock
    Integer, dimension(:), allocatable :: seed

	call random_seed(size = size_seed)
	allocate(seed(size_seed))
	call system_clock(COUNT=clock)
	seed = clock + 37 * (/ (i - 1, i = 1, size_seed) /)
	call random_seed(put = seed)  
	deallocate(seed)
	
End subroutine

subroutine setAEandBEfib(AE,BE,K,Res,NodG,idofT,iShiftU,NdimE)
	integer, intent(in) :: NodG,idofT,iShiftU,NdimE
	real*8 , intent(in) :: K(:,:), Res(:)
	real*8 , intent(out) :: AE(:,:), BE(:)
	Integer :: A, B, p, q, ApRow, BpCol, App, Bqp, signal
	
	Do A=1, NodG !LoopRow
		ApRow  = (A-1)*iDofT + iShiftU
				
		Do p = 1 , NdimE
			App = ApRow+p
			
			signal = (-1.0d0)**A
			
			BE(App) = signal*Res(p) 
			
			do B= 1 , NodG ! LoopCol ! not considering the simmetry
				BpCol  = (B-1)*iDofT + iShiftU
				
				signal = (-1.0d0)**(A + B)
				do q=1,NdimE
					Bqp = BpCol+q
					
					AE(App,Bqp) = signal*K(p,q)  
			
				end do ! loop Bq
			Enddo !LoopCol
		end do ! loop Ap
	Enddo !LoopRow

end subroutine

subroutine setAEandBE_lagrangeMultipliers(AE,BE,matB,matC,Sol1,eps,pen,NodG,NodLag,idofT,NdofLag,iShiftU,iShiftLag,NdimE)
	integer, intent(in) :: NodG,NodLag,idofT,NdofLag,iShiftU,iShiftLag,NdimE
	real*8 , intent(in) :: matB(:,:), matC(:,:), Sol1(*), eps,pen
	real*8 , intent(out) :: AE(:,:), BE(:)
	Integer :: A, B, p, q, ApRow, BpCol, App, Bqp, Apdim, AppDim

	B = NodLag
	BpCol = (B-1)*iDofT + iShiftLag

	Do A=1,NodG 
		ApRow  = (A-1)*iDofT + iShiftU
		ApDim = (A-1)*NdimE 
		
		Do p = 1 , NdimE
			App = ApRow+p
			AppDim = ApDim + p
												 
			Do q = 1 , NdofLag
				Bqp = BpCol + q
				 
				AE(Bqp,App) = pen*MatB(q,AppDim)
				AE(App,Bqp) = MatB(q,AppDim)
				
				AE(Bqp,Bqp) = MatC(q,q)
				
				BE(Bqp) = eps*Sol1(Bqp)

			end do
		end do ! loop Ap
	Enddo !LoopRow
	
end subroutine

! better version of setAEandBE_lagrangeMultipliers AE = [0 , B ; C , D] , BE = [0 , G] 
subroutine setAEandBE_LM(AE,BE,matB,matC,matD,vecG,NodG,NodLag,idofT,NdofLag,iShiftU,iShiftLag,NdimE) !
	integer, intent(in) :: NodG,NodLag,idofT,NdofLag,iShiftU,iShiftLag,NdimE
	real*8 , intent(in) :: matB(:,:), matC(:,:), matD(:,:), vecG(:)
	real*8 , intent(out) :: AE(:,:), BE(:)
	Integer :: A, B, p, q, ApRow, BpCol, App, Bqp, Apdim, AppDim

	B = NodLag
	BpCol = (B-1)*iDofT + iShiftLag

	Do A=1,NodG 
		ApRow  = (A-1)*iDofT + iShiftU
		ApDim = (A-1)*NdimE 
		
		Do p = 1 , NdimE
			App = ApRow+p
			AppDim = ApDim + p
												 
			Do q = 1 , NdofLag
				Bqp = BpCol + q
				 
				AE(Bqp,App) = MatC(q,AppDim)
				AE(App,Bqp) = MatB(AppDim,q)
				
				AE(Bqp,Bqp) = MatD(q,q)
				
				BE(Bqp) = vecG(q)

			end do
		end do ! loop Ap
	Enddo !LoopRow
	
end subroutine


! more complete version of setAEandBE_LM AE = [0 , B ; C , D] , BE = [F , G] 
subroutine setAEandBE_LM_full(AE,BE,matB,matC,matD,vecF,vecG,NodG,NodLag,idofT,NdofLag,iShiftU,iShiftLag,NdimE) !
	integer, intent(in) :: NodG,NodLag,idofT,NdofLag,iShiftU,iShiftLag,NdimE
	real*8 , intent(in) :: matB(:,:), matC(:,:), matD(:,:), vecG(:), vecF(:)
	real*8 , intent(out) :: AE(:,:), BE(:)
	Integer :: A, B, p, q, ApRow, BpCol, App, Bqp, Apdim, AppDim

	B = NodLag
	BpCol = (B-1)*iDofT + iShiftLag

	Do A=1,NodG 
		ApRow  = (A-1)*iDofT + iShiftU
		ApDim = (A-1)*NdimE 
		
		Do p = 1 , NdimE
			App = ApRow+p
			AppDim = ApDim + p
			
			BE(App) = VecF(AppDim)
												 
			Do q = 1 , NdofLag
				Bqp = BpCol + q
				 
				AE(Bqp,App) = MatC(q,AppDim)
				AE(App,Bqp) = MatB(AppDim,q)
				
				AE(Bqp,Bqp) = MatD(q,q)
				
				BE(Bqp) = vecG(q)

			end do
		end do ! loop Ap
	Enddo !LoopRow
	
end subroutine


subroutine setCoupling_lagrangeMultipliers(Coupling,NodG,NodLag,idofT,NdofLag,iShiftU,iShiftLag,NdimE)

	integer, intent(in) :: NodG,NodLag,idofT,NdofLag,iShiftU,iShiftLag,NdimE
	integer, intent(out) :: Coupling(:,:)
	Integer :: A, B, p, q, ApRow, BpCol, App, Bqp
	
	B = NodLag
	BpCol = (B-1)*iDofT + iShiftLag
			
	Do A=1,NodG
		ApRow  = (A-1)*iDofT + iShiftU
		
		Do p = 1 , NdimE
			App = ApRow+p
												 
			Do q = 1 , NdofLag
				Bqp = BpCol + q
!~ 			
				Coupling(App,Bqp) = 1
				Coupling(Bqp,App) = 1
				Coupling(Bqp,Bqp) = 1
			end do
		end do ! loop Ap
	Enddo !LoopRow
	
end subroutine

subroutine setCoupling_pureDisplacement(Coupling,NodG,idofT,iShiftU,NdimE)
	integer, intent(in) :: NodG,idofT,iShiftU,NdimE
	integer, intent(out) :: Coupling(:,:)
	Integer :: A, B, p, q, ApRow, BpCol
			
    do A = 1,NodG
        ApRow = (A-1) * iDofT + iShiftU
		
		do B = 1,NodG
			BpCol = (B-1) * iDofT + iShiftU

			do p=1,nDimE
			do q=1,nDimE

				Coupling (ApRow + p, BpCol + q) = 1

			enddo
            enddo
		enddo
	enddo
	
end subroutine


logical function isOnBoundary(y) result(flag)
	real*8, intent(in) :: y(:) 
	real*8, parameter :: a = 0.0d0, b = 0.0d0 , c = 1.0d0 , d = 1.0d0, eps = 1.0e-10
	
	flag = .false.

	if(dabs(y(1) - a) < eps) then
		flag = .true.
	else if(dabs(y(2) - b) < eps) then
		flag = .true.
	else if(dabs(y(1) - c) < eps) then
		flag = .true.
	else if(dabs(y(2) - d) < eps) then
		flag = .true.	
	end if
		
end function

subroutine getNormalBoundary(normal,flag,y)
	real*8, intent(out) :: normal(:)
	logical , intent(out) :: flag 
	real*8, intent(in) :: y(:) 
	real*8, parameter :: a = 0.0d0, b = 0.0d0 , c = 1.0d0 , d = 1.0d0, eps = 1.0e-10
	
	flag = .false.

	if(dabs(y(1) - a) < eps) then
		flag = .true.
		normal(1) = -1.0d0
		normal(2) = 0.0d0
	else if(dabs(y(2) - b) < eps) then
		flag = .true.
		normal(1) = 0.0d0
		normal(2) = -1.0d0
	else if(dabs(y(1) - c) < eps) then
		flag = .true.
		normal(1) = 1.0d0
		normal(2) = 0.0d0
	else if(dabs(y(2) - d) < eps) then
		flag = .true.
		normal(1) = 0.0d0
		normal(2) = 1.0d0	
	end if
		
end subroutine

real*8 function norm2(v) result(normV)
	real*8, intent(in) :: v(:) 
	
	normV = dsqrt(dot_product(v,v))
	
end function

subroutine setBlockToVector(A,B, NodG, iDofA, iDofB, iShiftA)
	real*8, intent(out) :: A(:)
	real*8, intent(in) :: B(:)
	integer , intent(in) :: NodG, iDofA, iDofB, iShiftA
	integer :: i,p,ipA,ipB 	
	
	do i=1, NodG !LoopRow
		ipA  = (i-1)*iDofA + iShiftA
		ipB  = (i-1)*iDofB 
		
		do p = 1 , iDofB
			
			A(ipA + p) = B(ipB + p)
		
		end do ! loop Ap
	end do !LoopRow
	
end subroutine


subroutine setBlockToMatrixSquare(A,B, NodG, iDofA, iDofB, iShiftA)
	real*8, intent(out) :: A(:,:)
	real*8, intent(in) :: B(:,:)
	integer , intent(in) :: NodG, iDofA, iDofB, iShiftA
	
	call setBlockToMatrix(A,B,NodG,NodG, iDofA, iDofA, iDofB, iDofB, iShiftA, iShiftA)
	
end subroutine

subroutine setBlockToMatrix(A,B, NodGi,NodGj, iDofAi,iDofAj, iDofBi,iDofBj, iShiftAi, iShiftAj)
	real*8, intent(out) :: A(:,:)
	real*8, intent(in) :: B(:,:)
	integer , intent(in) :: NodGi,NodGj, iDofAi,iDofAj, iDofBi,iDofBj, iShiftAi, iShiftAj
	integer :: i,j,p,q,ipA,ipB,jqA,jqB 	
	
	do i=1, NodGi !LoopRow
		ipA  = (i-1)*iDofAi + iShiftAi
		ipB  = (i-1)*iDofBi 
		
		do p = 1 , iDofBi
			
			do j= 1 , NodGj !
				jqA  = (j-1)*iDofAj + iShiftAj
				jqB  = (j-1)*iDofBj
				
				do q=1, iDofBj
					
					A(ipA + p, jqA + q) = B(ipB + p, jqB + q)
						
				end do ! loop Bq
			end do !LoopCol
		end do ! loop Ap
	end do !LoopRow
	
end subroutine

subroutine setBlockDiagonalConstantToMatrix(A,cte,NodG, iDofA, iDofB, iShiftAi, iShiftAj)
	real*8, intent(out) :: A(:,:)
	real*8, intent(in) :: cte
	integer , intent(in) :: NodG, iDofA, iDofB, iShiftAi, iShiftAj
	integer :: i,p,ipAi,ipAj 	
	

	Do i=1, NodG !LoopRow
		ipAi  = (i-1)*iDofA + iShiftAi
		ipAj  = (i-1)*iDofA + iShiftAj
		
		Do p = 1 , iDofB

			A(ipAi + p , ipAj + p) = cte
			
		end do ! loop Ap
	end do !LoopRow
	
end subroutine 

subroutine setSimbolicBlockToMatrixSquare(A, NodG, iDofA, iDofB, iShiftA)
	integer, intent(out) :: A(:,:)
	integer , intent(in) :: NodG, iDofA, iDofB, iShiftA
	
	call setSimbolicBlockToMatrix(A,NodG,NodG, iDofA, iDofA, iDofB, iDofB, iShiftA, iShiftA)
	
end subroutine

subroutine setSimbolicBlockToMatrix(A,NodGi,NodGj, iDofAi,iDofAj, iDofBi,iDofBj, iShiftAi, iShiftAj)
	integer, intent(out) :: A(:,:)
	integer , intent(in) :: NodGi,NodGj, iDofAi,iDofAj, iDofBi,iDofBj, iShiftAi, iShiftAj
	integer :: i,j,p,q,ipA,jqA 	
	
	do i=1, NodGi !LoopRow
		ipA  = (i-1)*iDofAi + iShiftAi
		
		do p = 1 , iDofBi
			
			do j= 1 , NodGj !
				jqA  = (j-1)*iDofAj + iShiftAj
				
				do q=1, iDofBj
					
					A(ipA + p, jqA + q) = 1
						
				end do ! loop Bq
			end do !LoopCol
		end do ! loop Ap
	end do !LoopRow
	
end subroutine

subroutine setSimbolicBlockDiagonalToMatrix(A,NodG, iDofA, iDofB, iShiftAi, iShiftAj)
	integer, intent(out) :: A(:,:)
	integer , intent(in) :: NodG, iDofA, iDofB, iShiftAi, iShiftAj
	integer :: i,p,ipAi,ipAj 	
	

	Do i=1, NodG !LoopRow
		ipAi  = (i-1)*iDofA + iShiftAi
		ipAj  = (i-1)*iDofA + iShiftAj
		
		Do p = 1 , iDofB

			A(ipAi + p , ipAj + p) = 1
			
		end do ! loop Ap
	end do !LoopRow
	
end subroutine 

subroutine voigtTen4toTen2(V,W)
	real*8 , intent(out) :: V(:,:)
	real*8 , intent(in) :: W(:,:,:,:)
	
	integer , parameter :: n = 2, n2 = 4
	integer :: voigtMapI(n2), voigtMapJ(n2)
	
	voigtMapI = (/1,1,2,2/)
	voigtMapJ = (/1,2,1,2/)
	
	do i = 1,n2
	do j = 1,n2
		V(i,j) = W(voigtMapI(i),voigtMapJ(i),voigtMapI(j),voigtMapJ(j))
	end do
	end do

end subroutine

subroutine voigtTen2toVec(V,W)
	real*8 , intent(out) :: V(:)
	real*8 , intent(in) :: W(:,:)
	
	integer , parameter :: n = 2, n2 = 4
	integer :: voigtMapI(n2), voigtMapJ(n2)
	
	voigtMapI = (/1,1,2,2/)
	voigtMapJ = (/1,2,1,2/)
	
	do i = 1,n2
		V(i) = W(voigtMapI(i),voigtMapJ(i))
	end do

end subroutine

subroutine getNormal(n,flag)
	implicit none
	integer, intent(in) :: flag
	real*8 , intent(out) :: n(2,2) ! it may have two normals

	n = 0.0d0
	
	select case(flag)
		
		case(1)
			n(1,1) = 0.0d0 ; n(2,1) = -1.0d0
		case(2)
			n(1,1) = 1.0d0 ; n(2,1) = 0.0d0
		case(3)
			n(1,1) = 0.0d0 ; n(2,1) = 1.0d0
		case(4)
			n(1,1) = -1.0d0 ; n(2,1) = 0.0d0
		case(5)
			n(1,1) = -1.0d0 ; n(2,1) = 0.0d0
			n(1,2) = 0.0d0 ; n(2,2) = -1.0d0
		case(6)
			n(1,1) = 1.0d0 ; n(2,1) = 0.0d0
			n(1,2) = 0.0d0 ; n(2,2) = -1.0d0
		case(7)
			n(1,1) = 1.0d0 ; n(2,1) = 0.0d0
			n(1,2) = 0.0d0 ; n(2,2) = 1.0d0
		case(8)
			n(1,1) = -1.0d0 ; n(2,1) = 0.0d0
			n(1,2) = 0.0d0 ; n(2,2) = 1.0d0
	end select

end subroutine

real*8 function solveQuadraticEquation(a,b,c,sig) result(x)
	real*8 , intent(in) :: a,b,c
	integer , intent(in) :: sig
	real*8 :: Delta
		
	Delta = b*b - 4.0d0*a*c
	x = (-b + sig * dsqrt(Delta))/(2.0d0*a)

end function

function transposeT4(A) result(B)
	real*8, intent(in) :: A(:,:,:,:)
	real*8 , allocatable :: B(:,:,:,:) 
	integer :: n1 , n2, n3, n4, i, j, k, l

	n1 = size(A,1)
	n2 = size(A,2)
	n3 = size(A,3)
	n4 = size(A,4)
	
	allocate(B(n4,n3,n2,n1))
	
	do i = 1, n1
	do j = 1, n2
	do k = 1, n3
	do l = 1, n4
		B(l,k,j,i) = A(i,j,k,l)
	end do
	end do
	end do
	end do
	
end function


subroutine FindKeyword(ioUnit,str1)
! Str1: String to find, ioUnit: Unit assigned to Input File
	Character(*) str1,str2*120
	integer :: n , iError, ioUnit
	
!	rewind(IoUnit)   
	iError=0
	n = Len_Trim(str1)
	!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	do while (.true.)
		read(ioUnit, "(1A120)",iostat=iError) str2
		if (iError.lt.0 ) then
			write(0,*) "Keyword not finded" 
			STOP
		end if
		if ( str2(1:n) == str1(1:n) )  exit 
	end do
	
end subroutine

end module
