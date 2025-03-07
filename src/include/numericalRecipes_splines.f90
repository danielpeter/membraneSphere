! numerical recipes:
! routines from
! http://www.physics.louisville.edu/help/nr/bookf90pdf/chap3f9.pdf
! and
! http://www.physics.louisville.edu/help/nr/bookfpdf/f3-6.pdf
! and
! http://www.phys.uu.nl/DU/num_recipes/Fortran.208/f90/recipes/


function splin2(x1a,x2a,ya,y2a,x1,x2)
!Given x1a, x2a, ya as described in splie2 and y2a as produced by that routine; and given
!a desired interpolating point x1,x2; this routine returns an interpolated function value by
!bicubic spline interpolation.
	USE nrtype; USE nrutil, only: assert_eq
	USE nr, only: spline,splint
	implicit none
	REAL(SP), DIMENSION(:), INTENT(IN) :: x1a,x2a
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya,y2a
	REAL(SP), INTENT(IN) :: x1,x2
	REAL(SP) :: splin2
	INTEGER(I4B) :: j,m,ndum
	REAL(SP), DIMENSION(size(x1a)) :: yytmp,y2tmp2
	m = assert_eq(size(x1a),size(ya,1),size(y2a,1),'splin2: m')
	ndum = assert_eq(size(x2a),size(ya,2),size(y2a,2),'splin2: ndum')
	do j = 1,m
		yytmp(j)=splint(x2a,ya(j,:),y2a(j,:),x2)
	enddo
	call spline(x1a,yytmp,1.0e30_sp,1.0e30_sp,y2tmp2)
	splin2=splint(x1a,yytmp,y2tmp2,x1)
	end function splin2

subroutine splie2(x1a,x2a,ya,y2a)
!Given an MxN tabulated function ya, and N tabulated independent variables x2a, this
!routine constructs one-dimensional natural cubic splines of the rows of ya and returns the
!second derivatives in the MxN array y2a. (The array x1a is included in the argument
!list merely for consistency with routine splin2.)
	USE nrtype; USE nrutil, only: assert_eq
	USE nr, only: spline
	implicit none
	REAL(SP), DIMENSION(:), INTENT(IN) :: x1a,x2a
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: y2a
	INTEGER(I4B) :: j,m,ndum
	m = assert_eq(size(x1a),size(ya,1),size(y2a,1),'splie2: m')
	ndum = assert_eq(size(x2a),size(ya,2),size(y2a,2),'splie2: ndum')
	do j = 1,m
		call spline(x2a,ya(j,:),1.0e30_sp,1.0e30_sp,y2a(j,:))
	enddo
	end subroutine splie2


subroutine spline(x,y,yp1,ypn,y2)
!Given arrays x and y of length N containing a tabulated function, i.e., y
!i = f(xi),with x1 < x2 < ... < xN, and given values yp1 and ypn for the first derivative of the interpolating
!function at points 1 and N, respectively, this routine returns an array y2 of length N
!that contains the second derivatives of the interpolating function at the tabulated points
! xi. If yp1 and/or ypn are equal to 1x10^30 or larger, the routine is signaled to set the
!corresponding boundary condition for a natural spline, with zero second derivative on that
!boundary.
	USE nrtype; USE nrutil, only: assert_eq
	USE nr, only: tridag
	implicit none
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(SP), INTENT(IN) :: yp1,ypn
	REAL(SP), DIMENSION(:), INTENT(OUT) :: y2
	INTEGER(I4B) :: n
	REAL(SP), DIMENSION(size(x)) :: a,b,c,r
	n = assert_eq(size(x),size(y),size(y2),'spline')
	c(1:n-1)=x(2:n)-x(1:n-1)
	r(1:n-1)=6.0_sp*((y(2:n)-y(1:n-1))/c(1:n-1))
	r(2:n-1)=r(2:n-1)-r(1:n-2)
	a(2:n-1)=c(1:n-2)
	b(2:n-1)=2.0_sp*(c(2:n-1)+a(2:n-1))
	b(1)=1.0
	b(n)=1.0
	if (yp1 > 0.99e30_sp) then
		r(1)=0.0
		c(1)=0.0
	else
		r(1)=(3.0_sp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
		c(1)=0.5
	endif
	if (ypn > 0.99e30_sp) then
		r(n)=0.0
		a(n)=0.0
	else
		r(n)=(-3.0_sp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
		a(n)=0.5
	endif
	call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
	end subroutine spline

function splint(xa,ya,y2a,x)
!Given the arrays xa and ya, which tabulate a function (with the xai's in increasing or
!decreasing order), and given the array y2a, which is the output from spline above, and
!given a value of x, this routine returns a cubic-spline interpolated value. The arrays xa, ya
!and y2a are all of the same size.
	USE nrtype; USE nrutil, only: assert_eq,nrerror
	USE nr, only: locate
	implicit none
	REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: splint
	INTEGER(I4B) :: khi,klo,n
	REAL(SP) :: a,b,h
	n = assert_eq(size(xa),size(ya),size(y2a),'splint')
	klo = max(min(locate(xa,x),n-1),1)
  ! We will find the right place in the table by means of locate's bisection algorithm. This is
  ! optimal if sequential calls to this routine are at random values of x. If sequential calls are in
  ! order, and closely spaced, one would do better to store previous values of klo and khi and
  ! test if they remain appropriate on the next call.
	khi = klo+1
	h = xa(khi)-xa(klo)
	if (h == 0.0) call nrerror('bad xa input in splint')
	a=(xa(khi)-x)/h
	b=(x-xa(klo))/h
	splint = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_sp
	end function splint

function locate(xx,x)
!Given an array xx(1:N), and given a value x, returns a value j such that x is between
!xx(j) and xx(j+1). xx must be monotonic, either increasing or decreasing. j =0 or
!j =N is returned to indicate that x is out of range.
	USE nrtype
	implicit none
	REAL(SP), DIMENSION(:), INTENT(IN) :: xx
	REAL(SP), INTENT(IN) :: x
	INTEGER(I4B) :: locate
	INTEGER(I4B) :: n,jl,jm,ju
	LOGICAL :: ascnd
	n=size(xx)
	ascnd = (xx(n) >= xx(1))
	jl = 0
	ju = n+1
	do
		if (ju-jl <= 1) exit
		jm=(ju+jl)/2
		if (ascnd .eqv. (x >= xx(jm))) then
			jl = jm
		else
			ju = jm
		endif
	enddo
	if (x == xx(1)) then
		locate = 1
	else if (x == xx(n)) then
		locate = n-1
	else
		locate = jl
	endif
	end function locate

subroutine tridag_ser(a,b,c,r,u)
	USE nrtype; USE nrutil, only: assert_eq,nrerror
	implicit none
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
	REAL(SP), DIMENSION(:), INTENT(OUT) :: u
	REAL(SP), DIMENSION(size(b)) :: gam
	INTEGER(I4B) :: n,j
	REAL(SP) :: bet
	n = assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_ser')
	bet=b(1)
	if (bet == 0.0) call nrerror('tridag_ser: Error at code stage 1')
	u(1)=r(1)/bet
	do j = 2,n
		gam(j)=c(j-1)/bet
		bet = b(j)-a(j-1)*gam(j)
		if (bet == 0.0) &
			call nrerror('tridag_ser: Error at code stage 2')
		u(j)=(r(j)-a(j-1)*u(j-1))/bet
	enddo
	do j = n-1,1,-1
		u(j)=u(j)-gam(j+1)*u(j+1)
	enddo
	end subroutine tridag_ser

RECURSIVE subroutine tridag_par(a,b,c,r,u)
	USE nrtype; USE nrutil, only: assert_eq,nrerror
	USE nr, only: tridag_ser
	implicit none
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
	REAL(SP), DIMENSION(:), INTENT(OUT) :: u
	INTEGER(I4B), PARAMETER :: NPAR_TRIDAG=4
	INTEGER(I4B) :: n,n2,nm,nx
	REAL(SP), DIMENSION(size(b)/2) :: y,q,piva
	REAL(SP), DIMENSION(size(b)/2-1) :: x,z
	REAL(SP), DIMENSION(size(a)/2) :: pivc
	n = assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_par')
	if (n < NPAR_TRIDAG) then
		call tridag_ser(a,b,c,r,u)
	else
		if (maxval(abs(b(1:n))) == 0.0) &
			call nrerror('tridag_par: possible singular matrix')
		n2=size(y)
		nm=size(pivc)
		nx=size(x)
		piva = a(1:n-1:2)/b(1:n-1:2)
		pivc = c(2:n-1:2)/b(3:n:2)
		y(1:nm) = b(2:n-1:2)-piva(1:nm)*c(1:n-2:2)-pivc*a(2:n-1:2)
		q(1:nm) = r(2:n-1:2)-piva(1:nm)*r(1:n-2:2)-pivc*r(3:n:2)
		if (nm < n2) then
			y(n2) = b(n)-piva(n2)*c(n-1)
			q(n2) = r(n)-piva(n2)*r(n-1)
		endif
		x = -piva(2:n2)*a(2:n-2:2)
		z = -pivc(1:nx)*c(3:n-1:2)
		call tridag_par(x,y,z,q,u(2:n:2))
		u(1) = (r(1)-c(1)*u(2))/b(1)
		u(3:n-1:2) = (r(3:n-1:2)-a(2:n-2:2)*u(2:n-2:2) &
			-c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
		if (nm == n2) u(n)=(r(n)-a(n-1)*u(n-1))/b(n)
	endif
	end subroutine tridag_par


