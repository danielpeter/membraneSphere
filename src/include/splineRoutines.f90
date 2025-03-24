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
! subroutine cubicspline_setup: calculates the coefficients of the representation
!
! then find the function values by...
! function cubicspline_eval: calculates the function value at a specific location

!-----------------------------------------------------------------------
  subroutine cubicspline_setup()
!-----------------------------------------------------------------------
! computes cubic spline interpolation coefficients
!
! for y(x) between grid points i1 and i2 saving them in q. The
! interpolation is continuous with continuous first and second
! derivatives. It agrees exactly with y at grid points and with the
! three point first derivatives at both end points (i1 and i2).
! x must be monotonic but if two successive values of x are equal
! a discontinuity is assumed and separate interpolation is done on
! each strictly monotonic segment. The arrays must be dimensioned at
! least - x(i2), y(i2), q(3,i2), and f(3,i2). f is working storage.
  use splinefunction, only: i1,i2,X,Y,Q,F
  implicit none
  ! local parameters
  integer :: i, j, j1, j2, k
  double precision:: a0, b0, b1, h, h2, h2a, h2b, h3a, ha, y0
  double precision:: yy(3)
  double precision, parameter :: small = 1.0d-10

  ! initializes
  Q(:,:) = 0.0d0
  F(:,:) = 0.0d0
  yy(:) = 0.0d0

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

  end subroutine cubicspline_setup


!-----------------------------------------------------------------------
  double precision function cubicspline_eval(s)
!-----------------------------------------------------------------------
! returns the value of the function y(x) evaluated at point s
! using the cubic spline coefficients computed by cubicspline_setup and saved in q.
! If s is outside the interval (x(i1), x(i2)), it extrapolates
! using the first or last interpolation polynomial.
! The arrays must be dimensioned at least: x(i2), y(i2), and q(3,i2).
  use splinefunction, only: i1,i2,X,Y,Q
  implicit none
  double precision, intent(in) :: s
  ! local parameters
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
  cubicspline_eval = y(i) + h*(q(1,i) + h*(q(2,i) + h*q(3,i)))

  return
  end function


!-----------------------------------------------------------------------
  function myspline_func(trace)
!-----------------------------------------------------------------------
  implicit none
  double precision, dimension(:),intent(in) :: trace
  double precision, dimension(size(trace)) :: myspline_func
  ! local parameters
  integer:: i
  double precision,external:: cubicspline_eval
  !interface
  !  function splineRepresentation(location)
  !  USE nrtype
  !  implicit none
  !  REAL(DP),INTENT(IN) :: location
  !  REAL(DP):: splineRepresentation
  !  end function splineRepresentation
  !end interface

  do i = 1,size(trace)
    !spline_func(i) = splineRepresentation(x(i))
    myspline_func(i) = cubicspline_eval(trace(i))
  enddo

  return
  end function


!-----------------------------------------------------------------------
  double precision function cubicspline_derivative(x, step_size, error_estimate)
!-----------------------------------------------------------------------
! Computes the first derivative of a function stored in the spline module
! at a given location using Richardson extrapolation.
!
! Input:
!   x            - Point at which to compute the derivative.
!   step_size    - Initial step size for numerical differentiation.
!
! Output:
!   spline_derivative - Computed first derivative at `x`.
!   error_estimate    - Estimated error in the computed derivative.
!
! Method:
!   - Uses Richardson extrapolation with decreasing step sizes.
!   - Evaluates cubic spline representation at adjusted points.
!   - Uses an adaptive approach to minimize error.
!
! Dependencies:
!   - Uses `cubicspline_eval` for function evaluation.
!   - `assert` ensures valid inputs.
  implicit none
  double precision, intent(in) :: x, step_size
  double precision, intent(out) :: error_estimate
  ! local parameters
  ! Constants
  integer, parameter :: max_iterations = 10
  real, parameter :: convergence_factor = 1.4
  real, parameter :: convergence_squared = convergence_factor * convergence_factor
  real, parameter :: large_value = 1.e+30
  real, parameter :: safety_factor = 2.0
  ! Internal variables
  integer :: best_index, i, j
  integer :: imin(1)
  real :: current_step
  real, dimension(max_iterations - 1) :: error_table, scaling_factors
  real, dimension(max_iterations, max_iterations) :: richardson_table
  interface
    function geometric_progression(first, factor, n)
      implicit none
      real, intent(in):: first, factor
      integer, intent(in):: n
      real, dimension(n):: geometric_progression
    end function
  end interface
  ! External function for cubic spline evaluation
  double precision, external :: cubicspline_eval

  ! Ensure step_size is not zero
  if (step_size == 0.0) call stopProgram('cubicspline_derivative: step_size must be nonzero    ')

  ! Initialize step size and error estimate
  current_step = step_size
  error_estimate = large_value ! Set large initial error

  ! Compute initial difference quotient
  richardson_table(1,1) = (cubicspline_eval(dble(x + current_step)) - cubicspline_eval(dble(x - current_step))) &
                          / (2.0 * current_step)



  ! Compute scaling factors for Richardson extrapolation
  scaling_factors(1:max_iterations-1) = geometric_progression(convergence_squared, convergence_squared, max_iterations-1)

  ! Iterative refinement using Richardson extrapolation
  do i = 2, max_iterations
    current_step = current_step / convergence_factor

    ! Compute next finite difference approximation
    richardson_table(1,i) = (cubicspline_eval(dble(x + current_step)) - cubicspline_eval(dble(x - current_step))) &
                            / (2.0 * current_step)

    ! Extrapolation loop
    do j = 2, i
      richardson_table(j,i) = (richardson_table(j-1,i) * scaling_factors(j-1) - richardson_table(j-1,i-1)) &
                              / (scaling_factors(j-1) - 1.0)
    enddo

    ! Compute error estimates
    error_table(1:i-1) = max(abs(richardson_table(2:i,i) - richardson_table(1:i-1,i)), &
                              abs(richardson_table(2:i,i) - richardson_table(1:i-1,i-1)))

    ! Find index of minimum error
    imin = minloc(error_table(1:i-1))
    best_index = imin(1)

    ! Update best estimate if error decreases
    if (error_table(best_index) <= error_estimate) then
      error_estimate = error_table(best_index)
      cubicspline_derivative = dble(richardson_table(1 + best_index, i))
    endif

    ! Convergence check
    if (abs(richardson_table(i,i) - richardson_table(i-1,i-1)) >= safety_factor * error_estimate) then
      return
    endif
  enddo

  return
  end function

!-----------------------------------------------------------------------
  function geometric_progression(first, factor, n)
!-----------------------------------------------------------------------
! Computes the first `n` terms of a geometric sequence.
!
! Input:
!   first   - First term of the sequence.
!   factor  - Common ratio between terms.
!   n       - Number of terms to generate.
!
! Output:
!   geometric_progression - Array of `n` terms in the geometric sequence.
!
! Method:
!   - Uses direct multiplication for small `n`
  implicit none
  ! Input parameters
  real, intent(in):: first, factor
  integer, intent(in):: n
  ! Output array
  real, dimension(n):: geometric_progression
  ! local parameters
  integer:: k

  ! Ensure n is positive
  if (n <= 0) return

  ! Initialize first term
  geometric_progression(1) = first

  do k = 2, n
    geometric_progression(k) = geometric_progression(k - 1) * factor
  enddo

  return
  end function
