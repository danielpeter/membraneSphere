!=====================================================================
!
!       m e m b r a n e S p h e r e
!       --------------------------------------------------
!
!      Daniel Peter
!      (c) 2025
!
!      Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!
!=====================================================================

! cubic spline routines for respresentation of a discrete function
!
! handling: first call the...
! subroutine mydrspln: calculates the coefficients of the representation
!
! then find the function values by...
! function mydrsple: calculates the function value at a specific location

!-----------------------------------------------------------------------
subroutine mydrspln(i1, i2, x, y, q, f)
!-----------------------------------------------------------------------

  ! computes cubic spline interpolation coefficients
  ! for y(x) between grid points i1 and i2 saving them in q. The
  ! interpolation is continuous with continuous first and second
  ! derivatives. It agrees exactly with y at grid points and with the
  ! three point first derivatives at both end points (i1 and i2).
  ! x must be monotonic but if two successive values of x are equal
  ! a discontinuity is assumed and separate interpolation is done on
  ! each strictly monotonic segment. The arrays must be dimensioned at
  ! least - x(i2), y(i2), q(3,i2), and f(3,i2). f is working storage.

  implicit none

  ! Arguments
  integer, intent(in) :: i1, i2
  real(kind=8), intent(in) :: x(i2), y(i2)
  real(kind=8), intent(out) :: q(3,i2)
  real(kind=8), intent(out) :: f(3,i2)

  ! Local variables
  integer :: i, j, j1, j2, k
  real(kind=8) :: a0, b0, b1, h, h2, h2a, h2b, h3a, ha, y0
  real(kind=8) :: yy(3) = [0.0d0, 0.0d0, 0.0d0]
  real(kind=8), parameter :: small = 1.0d-10

  j1 = i1 + 1
  y0 = 0.0d0

  ! Bail out if there are less than two points total
  if (i2 < i1) then
    goto 13
  else if (i2 == i1) then
    goto 17
  else
    goto 8
  endif

  ! Position label 8 at the beginning of discontinuity search
8 a0 = x(j1-1)

  ! Search for discontinuities
  do i = j1, i2
    b0 = a0
    a0 = x(i)
    if (abs((a0-b0) / max(a0, b0)) < small) goto 4
  enddo
  ! This is label 3 (end of loop) - flow to label 17
  goto 17

17 j1 = j1 - 1
  j2 = i2 - 2
  goto 5

  4 j1 = j1 - 1
  j2 = i - 3

  ! See if there are enough points to interpolate (at least three)
5 if (j2 + 1 < j1) then
    goto 9
  else if (j2 + 1 == j1) then
    goto 10
  else
    goto 11
  endif

  ! Only two points. Use linear interpolation
10 j2 = j2 + 2

  y0 = (y(j2) - y(j1)) / (x(j2) - x(j1))
  do j = 1, 3
    q(j, j1) = yy(j)
    q(j, j2) = yy(j)
  enddo
  goto 12

  ! More than two points. Do spline interpolation
11 a0 = 0.0d0

  h = x(j1+1) - x(j1)
  h2 = x(j1+2) - x(j1)
  y0 = h * h2 * (h2 - h)
  h = h * h
  h2 = h2 * h2

  ! Calculate derivative at near end
  b0 = (y(j1) * (h - h2) + y(j1+1) * h2 - y(j1+2) * h) / y0
  b1 = b0

  ! Explicitly reduce banded matrix to an upper banded matrix
  do i = j1, j2
    h = x(i+1) - x(i)
    y0 = y(i+1) - y(i)
    h2 = h * h
    ha = h - a0
    h2a = h - 2.0d0 * a0
    h3a = 2.0d0 * h - 3.0d0 * a0
    h2b = h2 * b0
    q(1, i) = h2 / ha
    q(2, i) = -ha / (h2a * h2)
    q(3, i) = -h * h2a / h3a
    f(1, i) = (y0 - h * b0) / (h * ha)
    f(2, i) = (h2b - y0 * (2.0d0 * h - a0)) / (h * h2 * h2a)
    f(3, i) = -(h2b - 3.0d0 * y0 * ha) / (h * h3a)
    a0 = q(3, i)
    b0 = f(3, i)
  enddo

  ! Take care of last two rows
  i = j2 + 1
  h = x(i+1) - x(i)
  y0 = y(i+1) - y(i)
  h2 = h * h
  ha = h - a0
  h2a = h * ha
  h2b = h2 * b0 - y0 * (2.0d0 * h - a0)
  q(1, i) = h2 / ha
  f(1, i) = (y0 - h * b0) / h2a
  ha = x(j2) - x(i+1)
  y0 = -h * ha * (ha + h)
  ha = ha * ha

  ! Calculate derivative at far end
  y0 = (y(i+1) * (h2 - ha) + y(i) * ha - y(j2) * h2) / y0
  q(3, i) = (y0 * h2a + h2b) / (h * h2 * (h - 2.0d0 * a0))
  q(2, i) = f(1, i) - q(1, i) * q(3, i)

  ! Solve upper banded matrix by reverse iteration
  do j = j1, j2
    k = i - 1
    q(1, i) = f(3, k) - q(3, k) * q(2, i)
    q(3, k) = f(2, k) - q(2, k) * q(1, i)
    q(2, k) = f(1, k) - q(1, k) * q(3, k)
    i = k
  enddo
  q(1, i) = b1

  ! Fill in the last point with a linear extrapolation
  9 j2 = j2 + 2
  do j = 1, 3
    q(j, j2) = yy(j)
  enddo

  ! See if this discontinuity is the last
12 if (j2 < i2) then
    goto 6
  else
    goto 13
  endif

  ! No. Go back for more
6 j1 = j2 + 2
  if (j1 < i2) then
    goto 8
  else if (j1 == i2) then
    goto 8
  else
    goto 7
  endif

  ! There is only one point left after the latest discontinuity
7 do j = 1, 3
    q(j, i2) = yy(j)
  enddo

  ! Fini
13 return

  ! No need for these labels as they're now correctly placed in the code flow above
end subroutine mydrspln


!-----------------------------------------------------------------------
double precision function mydrsple(i1, i2, x, y, q, s)
!-----------------------------------------------------------------------

  ! returns the value of the function y(x) evaluated at point s
  ! using the cubic spline coefficients computed by mydrspln and saved in q.
  ! If s is outside the interval (x(i1), x(i2)), it extrapolates
  ! using the first or last interpolation polynomial.
  ! The arrays must be dimensioned at least: x(i2), y(i2), and q(3,i2).

  implicit none

  ! Arguments
  integer, intent(in) :: i1, i2
  double precision, intent(in) :: x(i2), y(i2), q(3,i2), s

  ! Local variables
  integer :: i, ii
  double precision :: h

  ! Initialize i (preserved from original algorithm)
  i = 1
  ii = i2 - 1

  ! Guarantee i within bounds
  i = max(i, i1)
  i = min(i, ii)

  ! See if x is increasing or decreasing
  if (x(i2) - x(i1) < 0.0d0) then
    ! X is decreasing
1   if (s - x(i) <= 0.0d0) then
3     if (s - x(i+1) < 0.0d0) then
        i = i + 1
        if (i - ii < 0) then
          goto 3
        else if (i - ii == 0) then
          continue ! Proceed to calculation
        else
          i = ii   ! Label 7
        endif
      endif
    else
      i = i - 1
      if (i - i1 < 0) then
        i = i1     ! Label 11
      else if (i - i1 == 0) then
        continue   ! Proceed to calculation
      else
        goto 1
      endif
    endif
  else
    ! X is increasing
2   if (s - x(i+1) < 0.0d0) then
8     if (s - x(i) <= 0.0d0) then
        i = i - 1
        if (i - i1 < 0) then
          i = i1   ! Label 11
        else if (i - i1 == 0) then
          continue ! Proceed to calculation
        else
          goto 8
        endif
      endif
    else
      i = i + 1
      if (i - ii < 0) then
        goto 2
      else if (i - ii == 0) then
        continue   ! Proceed to calculation
      else
        i = ii     ! Label 7
      endif
    endif
  endif

  ! Calculate drsple using spline coefficients in y and q (Label 6)
  h = s - x(i)
  mydrsple = y(i) + h*(q(1,i) + h*(q(2,i) + h*q(3,i)))

  return
end function mydrsple


