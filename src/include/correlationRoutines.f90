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


!-----------------------------------------------------------------------
subroutine correlation_traces(trace1,trace2,correl,N)
!-----------------------------------------------------------------------

! Computes the correlation of two real datasets inputArray1 and inputArray2 of length N (including any user-supplied zeropadding).
! N must be an integer power of 2. The answer is returned as the function correl, an array of length N.
! The answer is stored in wrap-around order, i.e., correlations at increasingly negative lags are in correl(N) on down to
! correl(N/2+1), while correlations at increasingly positive lags are in correl(1) (zero lag) on up to correl(N/2).
! Sign convention of this routine: if inputArray1 lags inputArray2, i.e., is shifted to the right of it, then correl will show
! a peak at positive lags.

  use constants, only: WP,IOUT
  use propagationStartup, only: datadirectory
  implicit none
  integer, intent(in):: N
  real(WP), dimension(N), intent(in) :: trace1,trace2
  real(WP), dimension(N), intent(out) :: correl
  ! local parameters
  complex(kind=8), dimension(N/2) :: fourierTransformed1,fourierTransformed2
  double precision :: dnorm
  integer:: nhalf,i
  real(WP), dimension(:), allocatable :: tmp_trace1,tmp_trace2
  logical, parameter:: DEBUG_TEST = .false.

  ! size checks
  if (N <= 0) call stopProgram('correlation_traces needs N > 0 arrays    ')
  if (mod(N,2) /= 0) call stopProgram('correlation_traces N must be multiple of 2   ')
  ! tests if n is a power of 2 (assuming n > 0)
  if (iand(N,N-1) /= 0) call stopProgram('correlation_traces N must be a power of 2    ')

  ! initializes
  correl(:) = 0.0_WP
  fourierTransformed1(:) = cmplx(0.d0,0.d0)
  fourierTransformed2(:) = cmplx(0.d0,0.d0)

  ! fourier transform arrays
  nhalf = N/2
  dnorm = 1.d0 / dble(nhalf) ! normalization factor

  ! debugging
  if (DEBUG_TEST) then
    print *,'  testing correlation traces:'
    print *,'    size N = ',N
    allocate(tmp_trace1(N),tmp_trace2(N))
    ! forward
    call FFT_real(trace1,N,+1,fourierTransformed1)
    call FFT_real(trace2,N,+1,fourierTransformed2)
    ! normalize for inverse FFT
    fourierTransformed1(:) = fourierTransformed1(:) * dnorm
    fourierTransformed2(:) = fourierTransformed2(:) * dnorm
    ! backward
    call FFT_real(tmp_trace1,N,-1,fourierTransformed1)
    call FFT_real(tmp_trace2,N,-1,fourierTransformed2)
    ! file output
    print *,'    printing to file: ',trim(datadirectory)//'CorrelTest_trace1.dat'
    open(IOUT,file=trim(datadirectory)//'CorrelTest_trace1.dat')
    write(IOUT,*) "# format: #i #trace #trace-fft-reconstructed"
    do i = 1,N
      write(IOUT,*) i,trace1(i),tmp_trace1(i)
    enddo
    close(IOUT)
    print *,'    printing to file: ',trim(datadirectory)//'CorrelTest_trace2.dat'
    open(IOUT,file=trim(datadirectory)//'CorrelTest_trace2.dat')
    write(IOUT,*) "# format: #i #trace #trace-fft-reconstructed"
    do i = 1,N
      write(IOUT,*) i,trace2(i),tmp_trace2(i)
    enddo
    close(IOUT)
    deallocate(tmp_trace1,tmp_trace2)
    ! re-initializes
    fourierTransformed1(:) = cmplx(0.d0,0.d0)
    fourierTransformed2(:) = cmplx(0.d0,0.d0)
    print *,'    test done'
  endif

  ! forward FFT
  call FFT_real(trace1,N,+1,fourierTransformed1)
  call FFT_real(trace2,N,+1,fourierTransformed2)

  ! correlation
  fourierTransformed1(1) = cmplx(real(fourierTransformed1(1))*real(fourierTransformed2(1)), &
                                 aimag(fourierTransformed1(1))*aimag(fourierTransformed2(1)), kind=8)
  fourierTransformed1(2:) = fourierTransformed1(2:) * conjg(fourierTransformed2(2:))

  ! normalization for inverse FFT
  fourierTransformed1(:) = fourierTransformed1(:) * dnorm

  ! inverse FFT
  call FFT_real(correl,N,-1,fourierTransformed1)

end subroutine


!-----------------------------------------------------------------------
function simpsons_integral(func, lower_bound, upper_bound)
!-----------------------------------------------------------------------
! Computes the integral of a given function using Simpson’s rule.
!
! This function applies iterative refinement based on the trapezoidal rule
! to approximate the integral of `func` over the interval [lower_bound, upper_bound].
!
! Parameters:
!   func         - The function to be integrated (external real function).
!   lower_bound  - The lower limit of integration.
!   upper_bound  - The upper limit of integration.
!
! returns:
!   simpsons_integral - The computed integral value.
!
! Notes:
!   - `EPS` sets the desired fractional accuracy.
!   - `MAX_ITERATIONS` controls the maximum number of refinement steps.
!   - The function stops execution if convergence is not achieved.
  use precisions
  implicit none
  ! function result
  real(WP) :: simpsons_integral
  ! Input parameters
  double precision, intent(in) :: lower_bound, upper_bound
  real(WP), external :: func  ! function to integrate
  ! local parameters
  ! Constants for accuracy and iteration limits
  integer, parameter :: MAX_ITERATIONS = 40
  double precision, parameter :: EPSILON = 1.d-4

  ! Local variables
  integer :: iteration
  double precision :: previous_approx, previous_trapz, current_trapz, current_approx

  ! Initialize previous values
  previous_trapz = 0.d0
  previous_approx = 0.d0

  ! Iteratively refine the integral using Simpson's rule
  do iteration = 1, MAX_ITERATIONS
    ! Compute the trapezoidal approximation
    call trapezoidal_rule_refine(func, lower_bound, upper_bound, current_trapz, iteration)

    ! Apply Simpson’s rule refinement
    current_approx = (4.d0 * current_trapz - previous_trapz) / 3.d0

    ! Check for convergence after a few iterations
    if (iteration > 5) then
      if (abs(current_approx - previous_approx) < EPSILON * abs(previous_approx) .or. &
          (current_approx == 0.d0 .and. previous_approx == 0.d0)) then
        simpsons_integral = current_approx
        return
      endif
    endif

    ! Update previous values for next iteration
    previous_approx = current_approx
    previous_trapz = current_trapz
  enddo

  ! Print error message if the method fails to converge
  print *, 'Error: simpsons_integral did not converge:', lower_bound, upper_bound, current_approx, current_trapz
  call stopProgram("simpsons_integral: exceeded maximum iterations    ")

  return
end function

!-----------------------------------------------------------------------
subroutine trapezoidal_rule_refine(func, lower_bound, upper_bound, integral, iteration)
!-----------------------------------------------------------------------
! Performs iterative refinement of the integral using the trapezoidal rule.
!
! This subroutine computes the `iteration`-th stage of refinement for an
! extended trapezoidal rule approximation of the integral of `func`
! over the interval [lower_bound, upper_bound].
!
! Parameters:
!   func         - The function to be integrated (external real function).
!   lower_bound  - The lower limit of integration.
!   upper_bound  - The upper limit of integration.
!   integral     - The refined integral estimate (output).
!   iteration    - The iteration step (1 for initial estimate, higher values for refinement).
!
! Notes:
!   - On the first call (`iteration` = 1), the function computes the basic
!     trapezoidal approximation.
!   - For subsequent calls (`iteration` ≥ 2), it refines the approximation by
!     adding `2^(iteration-2)` interior points.
  use precisions, only: WP
  implicit none
  ! Input parameters
  integer, intent(in) :: iteration
  double precision, intent(in) :: lower_bound, upper_bound
  double precision, intent(out) :: integral
  real(WP), external :: func  ! function to integrate
  ! local parameters
  integer :: num_points, index
  double precision :: step_size, sum_values, current_x

  ! Compute integral estimate
  if (iteration == 1) then
    ! Compute the initial trapezoidal estimate with only endpoints
    integral = 0.5 * (upper_bound - lower_bound) * (func(real(lower_bound,kind=WP)) + func(real(upper_bound,kind=WP)))
  else
    ! Compute the number of new interior points to add
    num_points = 2**(iteration - 2)
    step_size = (upper_bound - lower_bound) / dble(num_points)

    ! Initialize midpoint sum
    current_x = lower_bound + 0.5d0 * step_size
    sum_values = 0.d0

    ! Compute function values at new interior points
    do index = 1, num_points
      sum_values = sum_values + func(real(current_x,kind=WP))
      current_x = current_x + step_size
    enddo

    ! Update the integral estimate with new values
    integral = 0.5d0 * (integral + step_size * sum_values)
  endif

  return
end subroutine

!-----------------------------------------------------------------------
subroutine simpsons_integral_degree(func, polynomial_degree, lower_bound, upper_bound, integral, parameter_mu)
!-----------------------------------------------------------------------
! Computes the integral of a function using Simpson’s rule with a specified degree.
!
! This subroutine applies iterative refinement based on the trapezoidal rule
! to approximate the integral of `func` over the interval [lower_bound, upper_bound].
!
! Parameters:
!   func               - The function to be integrated (external real function).
!   polynomial_degree  - Degree of the polynomial approximation.
!   lower_bound        - The lower limit of integration.
!   upper_bound        - The upper limit of integration.
!   integral           - The computed integral value (output).
!   parameter_mu       - Additional parameter used in function evaluation.
!
! Notes:
!   - `EPSILON` sets the desired fractional accuracy.
!   - `MAX_ITERATIONS` controls the maximum number of refinement steps.
!   - The function stops execution if convergence is not achieved.
  use precisions
  implicit none
  ! Input parameters
  integer, intent(in) :: polynomial_degree
  real(WP), intent(in) :: lower_bound, upper_bound, parameter_mu
  real(WP), intent(out) :: integral
  real(WP), external :: func  ! function to integrate
  ! local parameters
  ! Constants for accuracy and iteration limits
  integer, parameter :: MAX_ITERATIONS = 40
  real(WP), parameter :: EPSILON = 1.e-4

  ! Local variables
  integer :: iteration
  double precision :: previous_approx, previous_trapz, current_trapz

  ! Initialize previous values
  previous_trapz = 0.0
  previous_approx = 0.0

  ! Iteratively refine the integral using Simpson's rule
  do iteration = 1, MAX_ITERATIONS
    ! Compute the trapezoidal approximation
    call trapezoidal_rule_refine_degree(func, polynomial_degree, dble(lower_bound), dble(upper_bound), &
                                        current_trapz, iteration, dble(parameter_mu))

    ! Apply Simpson’s rule refinement
    integral = (4.0 * current_trapz - previous_trapz) / 3.0

    ! Check for convergence after a few iterations
    if (iteration > 5) then
      if (abs(integral - previous_approx) < EPSILON * abs(previous_approx) .or. &
          (integral == 0.0 .and. previous_approx == 0.0)) return
    endif

    ! Update previous values for next iteration
    previous_approx = integral
    previous_trapz = current_trapz
  enddo

  ! Print error message if the method fails to converge
  print *, 'Error: simpsons_integral_degree did not converge:', polynomial_degree, lower_bound, upper_bound, integral
  call stopProgram("Exceeded maximum iterations in simpsons_rule_degree    ")

end subroutine

!-----------------------------------------------------------------------
subroutine trapezoidal_rule_refine_degree(func, polynomial_degree, lower_bound, upper_bound, &
                                          integral, iteration, parameter_mu)
!-----------------------------------------------------------------------
! Performs iterative refinement of the trapezoidal rule for numerical integration.
!
! This subroutine refines the integral approximation by progressively adding more points.
!
! Parameters:
!   func               - The function to be integrated (external real function).
!   polynomial_degree  - Degree of the polynomial approximation used in `func`.
!   lower_bound        - The lower limit of integration.
!   upper_bound        - The upper limit of integration.
!   integral           - The computed integral value (output).
!   iteration          - Current refinement iteration.
!   parameter_mu       - Additional parameter used in function evaluation.
!
! Notes:
!   - If called with `refinement_step = 1`, it computes the initial rough integral.
!   - For higher values of `refinement_step`, it refines the estimate.
  use precisions
  implicit none
  ! Input parameters
  integer, intent(in) :: iteration, polynomial_degree
  double precision, intent(in) :: lower_bound, upper_bound
  double precision, intent(in) :: parameter_mu
  ! Output parameter
  double precision, intent(out) :: integral
  ! function to be integrated
  real(WP), external :: func
  ! local parameters
  integer :: num_intervals, interval_index
  double precision :: interval_width, sum, num_points, current_x

  ! Compute the initial estimate when refinement_step = 1
  if (iteration == 1) then
    integral = 0.5 * (upper_bound - lower_bound) * &
              (func(polynomial_degree, real(lower_bound,kind=WP), real(parameter_mu,kind=WP)) + &
               func(polynomial_degree, real(upper_bound,kind=WP), real(parameter_mu,kind=WP)))

  else  ! Refinement step for higher iterations
    num_intervals = 2**(iteration - 2)
    num_points = num_intervals
    interval_width = (upper_bound - lower_bound) / num_points

    ! Initialize sum for midpoint evaluations
    current_x = lower_bound + 0.5 * interval_width
    sum = 0.0

    do interval_index = 1, num_intervals
      sum = sum + func(polynomial_degree, real(current_x,kind=WP), real(parameter_mu,kind=WP))
      current_x = current_x + interval_width
    enddo

    ! Update integral using refined trapezoidal estimate
    integral = 0.5 * (integral + (upper_bound - lower_bound) * sum / num_points)
  endif

  return
end subroutine

!-----------------------------------------------------------------------
function romberg_integral(func, lower_bound, upper_bound)
!-----------------------------------------------------------------------
! Performs numerical integration using Romberg's method.
!
! This function computes the integral of `func` from `lower_bound` to `upper_bound`
! using Romberg integration, which applies Richardson extrapolation on the trapezoidal rule.
!
! Parameters:
!   func           - The function to be integrated (external real function).
!   lower_bound    - The lower limit of integration.
!   upper_bound    - The upper limit of integration.
!
! returns:
!   The computed integral value.
!
! Constants:
!   EPS            - Desired fractional accuracy.
!   JMAX           - Maximum number of refinement steps.
!   K              - Number of points used for polynomial extrapolation.
!   EPS_ZERO       - Small threshold to prevent premature convergence.
!
! Dependencies:
!   - Calls `trapezoidal_rule_refine` for trapezoidal rule integration.
!   - Uses `polynomial_interpolation` for extrapolation.
  use precisions
  implicit none
  ! Output: Computed integral value
  real(WP) :: romberg_integral
  ! Input parameters
  double precision, intent(in) :: lower_bound, upper_bound
  ! function to be integrated
  real(WP), external :: func
  ! local parameters
  ! Constants
  integer, parameter :: JMAX = 40  ! Maximum number of steps
  integer, parameter :: JMAXP = JMAX + 1  ! Array size limit for step storage
  integer, parameter :: K = 5  ! Order of extrapolation
  integer, parameter :: KM = K - 1  ! Offset for array indexing in extrapolation
  real(WP), parameter :: EPS = 1.0e-4  ! Desired accuracy
  real(WP), parameter :: EPS_ZERO = 1.0e-4  ! Small threshold for zero-check
  ! Arrays to store step sizes and integral estimates
  double precision, dimension(JMAXP) :: step_sizes, integral_values
  real(WP) :: error_estimate

  ! Loop index
  integer :: iteration

  ! Initialize the first step size
  step_sizes(1) = 1.0

  ! Iterate to refine integration using the trapezoidal rule
  do iteration = 1, JMAX
    call trapezoidal_rule_refine(func, lower_bound, upper_bound, integral_values(iteration), iteration)

    if (iteration >= K) then
      ! Apply polynomial interpolation to improve accuracy
      call polynomial_interpolation(step_sizes(iteration-KM:iteration),integral_values(iteration-KM:iteration),K, &
                                    0.d0, romberg_integral, error_estimate)

      ! Check for convergence based on extrapolation error
      if (abs(error_estimate) <= EPS * abs(romberg_integral)) return

      ! Additional safeguard against premature convergence
      if (iteration > 5) then
        if (abs(error_estimate) <= EPS * abs(romberg_integral) .or. &
            (abs(romberg_integral) < EPS_ZERO .and. abs(error_estimate) < EPS_ZERO)) return
      endif
    endif

    ! Store the last integral estimate and update step size
    integral_values(iteration + 1) = integral_values(iteration)
    step_sizes(iteration + 1) = 0.25 * step_sizes(iteration)
  enddo

  ! Error handling: Too many iterations without convergence
  print *, 'Error: romberg_integral did not converge:', lower_bound, upper_bound, romberg_integral, error_estimate
  call stopProgram("Romberg integration exceeded maximum steps    ")

  return
end function


!-----------------------------------------------------------------------
subroutine polynomial_interpolation(x_values, y_values, num_points, x_target, y_result, error_estimate)
!-----------------------------------------------------------------------
! Performs polynomial interpolation using Neville's algorithm.
!
! This subroutine computes an interpolated value `y_result` at `x_target`
! based on `num_points` given data points `(x_values, y_values)`.
!
! Parameters:
!   x_values      - Array of x-coordinates of known data points.
!   y_values      - Array of corresponding y-coordinates.
!   num_points    - Number of data points (should not exceed `max_points`).
!   x_target      - The x-value at which interpolation is required.
!   y_result      - Output: Interpolated y-value at `x_target`.
!   error_estimate - Output: Estimated error of interpolation.
!
! Constraints:
!   - `num_points` must be ≤ `max_points` (default: 10).
!   - `x_values` should contain distinct values.
  use precisions
  implicit none
  ! Input parameters
  integer, intent(IN) :: num_points
  double precision, intent(IN) :: x_values(num_points), y_values(num_points), x_target
  ! Output parameters
  real(WP), intent(OUT) :: y_result, error_estimate
  ! local parameters
  ! Constants
  integer, parameter :: max_points = 10  ! Maximum allowed data points
  real(WP) :: coeff_c(max_points), coeff_d(max_points)
  real(WP) :: min_difference, temp_difference, ho, hp, w, denominator
  integer :: i, nearest_index, step

  ! Initialize to find the nearest data point to x_target
  nearest_index = 1
  min_difference = abs(x_target - x_values(1))

  do i = 1, num_points
    temp_difference = abs(x_target - x_values(i))
    if (temp_difference < min_difference) then
      nearest_index = i
      min_difference = temp_difference
    endif
    coeff_c(i) = y_values(i)  ! Copy data into coefficient arrays
    coeff_d(i) = y_values(i)
  enddo

  ! Set initial interpolated value
  y_result = y_values(nearest_index)
  nearest_index = nearest_index - 1

  ! Perform Neville's algorithm for polynomial interpolation
  do step = 1, num_points - 1
    do i = 1, num_points - step
      ho = x_values(i) - x_target
      hp = x_values(i + step) - x_target
      w = coeff_c(i + 1) - coeff_d(i)
      denominator = ho - hp

      if (denominator == 0.0) stop 'Error in polynomial interpolation: Division by zero'

      denominator = w / denominator
      coeff_d(i) = hp * denominator
      coeff_c(i) = ho * denominator
    enddo

    ! Choose the next correction term
    if (2 * nearest_index < num_points - step) then
      error_estimate = coeff_c(nearest_index + 1)
    else
      error_estimate = coeff_d(nearest_index)
      nearest_index = nearest_index - 1
    endif

    ! Update interpolated result
    y_result = y_result + error_estimate
  enddo

end subroutine


!-----------------------------------------------------------------------
function simpsons_integral_vecFunc(vecFunc, lower_bound, upper_bound)
!-----------------------------------------------------------------------
! Adaptive integration using Simpson’s rule.
!
! This function estimates the integral of `func` over the interval [a, b]
! using Simpson’s rule with recursive refinement via the trapezoidal rule.
!
! Input:
!   vecFunc                    - function to integrate, passed as an external function.
!   lower_bound, upper_bound   - Integration limits.
!
! Output:
!   qsimp  - Approximate integral of `func` over [a, b].
!
! Parameters:
!   JMAX   - Maximum number of iterations (refinement steps).
!   EPS    - Desired fractional accuracy.
!
! Method:
!   - Uses repeated trapezoidal rule refinements.
!   - Applies Richardson extrapolation to improve convergence.
!   - Stops early if convergence criteria are met.
  implicit none
  ! Input parameters
  double precision, intent(in) :: lower_bound, upper_bound
  double precision :: simpsons_integral_vecFunc
  ! function interface
  interface
    function vecFunc(x)
      double precision, dimension(:), intent(in) :: x
      double precision, dimension(size(x)) :: vecFunc
    end function
  end interface
  ! local parameters
  ! Constants
  integer, parameter :: MAX_ITERATIONS = 20
  double precision, parameter :: EPSILON = 1.0d-6
  ! Internal variables
  integer :: iteration
  double precision :: previous_approx, previous_trapz, current_trapz, current_approx

  ! Initialization
  previous_trapz = 0.d0
  previous_approx = 0.d0

  ! Iterative refinement
  do iteration = 1, MAX_ITERATIONS
    ! Compute trapezoidal rule approximation
    call trapezoidal_rule_refine_vecFunc(vecFunc, lower_bound, upper_bound, current_trapz, iteration)

    ! Apply Simpson's rule via extrapolation
    current_approx = (4.d0 * current_trapz - previous_trapz) / 3.d0

    ! Convergence check (after 5 iterations to avoid premature stopping)
    if (iteration > 5) then
      if (abs(current_approx - previous_approx) < EPSILON * abs(previous_approx) .or. &
          (current_approx == 0.0 .and. previous_approx == 0.0)) then
        simpsons_integral_vecFunc = current_approx
        return
      endif
    endif

    ! Update previous values
    previous_approx = simpsons_integral_vecFunc
    previous_trapz = current_trapz
  enddo

  ! Error handling: Maximum iterations reached
  print *, 'Error: simpsons_integral_vecFunc did not converge:', lower_bound, upper_bound, current_approx, current_trapz
  call stopProgram("simpsons_integral_vecFunc: exceeded maximum iterations    ")

  return
end function

!-----------------------------------------------------------------------
subroutine trapezoidal_rule_refine_vecFunc(vecFunc, lower_bound, upper_bound, integral, iteration)
!-----------------------------------------------------------------------
! Adaptive trapezoidal rule for numerical integration.
!
! This subroutine computes successive refinements of the integral of
! `func` over [a, b] using the extended trapezoidal rule.
!
! Input:
!   vecFunc                    - function to integrate, passed as an external function.
!   lower_bound, upper_bound   - Integration limits.
!   iteration                  - Refinement level (number of iterations).
!
! Input/Output:
!   integral                   - Integral estimate (refined in each call).
!
! Method:
!   - First call (n=1) computes the simplest trapezoidal rule estimate.
!   - Subsequent calls refine the estimate by adding new points.
!   - Uses the arithmetic sequence generator `arth` for efficient computation.
  implicit none
  ! Input parameters
  double precision, intent(in) :: lower_bound, upper_bound
  integer, intent(in) :: iteration
  ! Input/Output parameter
  double precision, intent(inout) :: integral
  ! function interface
  interface
    function vecFunc(x)
      double precision, dimension(:), intent(in) :: x
      double precision, dimension(size(x)) :: vecFunc
    end function
  end interface
  ! local parameters
  double precision :: step_size, fsum
  integer :: num_points

  ! Compute integral estimate
  if (iteration == 1) then
    ! Compute the initial trapezoidal estimate with only endpoints
    integral = 0.5d0 * (upper_bound - lower_bound) * sum(vecFunc((/ lower_bound, upper_bound /)))
  else
    ! Compute the number of new interior points to add
    num_points = 2**(iteration - 2)
    step_size = (upper_bound - lower_bound) / dble(num_points)

    ! Sum function evaluations at new points
    fsum = sum(vecFunc(arithmetic_sequence(lower_bound + 0.5d0 * step_size, step_size, num_points)))

    ! Refine previous integral estimate
    integral = 0.5d0 * (integral + step_size * fsum)
  endif

contains

  !-----------------------------------------------------------------------
  function arithmetic_sequence(first_term, step_size, num_terms)
  !-----------------------------------------------------------------------
  ! Generates an arithmetic sequence with a given first term, step size,
  ! and number of terms.
  !
  ! Arguments:
  !   first_term - The first term in the sequence.
  !   step_size  - The common difference between consecutive terms.
  !   num_terms  - The number of terms in the sequence.
  !
  ! returns:
  !   arithmetic_sequence - An array containing the arithmetic sequence.
    implicit none
    double precision, intent(in) :: first_term, step_size          ! Input: First term and step size
    integer, intent(in) :: num_terms                               ! Input: Number of terms in sequence
    double precision, dimension(num_terms) :: arithmetic_sequence  ! Output: Generated sequence
    ! local parameters
    integer :: index

    ! Initialize first term if sequence has elements
    if (num_terms > 0) arithmetic_sequence(1) = first_term

    ! use simple loop
    do index = 2, num_terms
      arithmetic_sequence(index) = arithmetic_sequence(index - 1) + step_size
    enddo

  end function

end subroutine
