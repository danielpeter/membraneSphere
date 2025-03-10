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
function mycorrel(inputArray1,inputArray2)
!-----------------------------------------------------------------------

! Computes the correlation of two real datasets inputArray1 and inputArray2 of length N (including any user-supplied zeropadding).
! N must be an integer power of 2. The answer is returned as the function correl, an array of length N.
! The answer is stored in wrap-around order, i.e., correlations at increasingly negative lags are in correl(N) on down to
! correl(N/2+1), while correlations at increasingly positive lags are in correl(1) (zero lag) on up to correl(N/2).
! Sign convention of this routine: if inputArray1 lags inputArray2, i.e., is shifted to the right of it, then correl will show
! a peak at positive lags.

  !USE nrtype; USE nrutil, only: assert,assert_eq
  !USE nr, only: realft
  use precisions, only: WP
  implicit none
  real(WP), dimension(:), intent(inout) :: inputArray1,inputArray2
  real(WP), dimension(size(inputArray1)) :: mycorrel
  complex(kind=WP), dimension(size(inputArray1)/2) :: fourierTransformed1,fourierTransformed2
  integer:: n,halfSize

  !print *,'processing correlation...'
  n = size(inputArray1)
  halfSize = n/2

  !print *,'checking...'
  if (n <= 0) call stopProgram('mycorrel needs n > 0 arrays    ')
  if (n /= size(inputArray2)) call stopProgram('mycorrel needs arrays with same lengths   ')
  ! tests if n is a power of 2 (assuming n > 0)
  if (iand(n,n-1) == 0) call stopProgram('mycorrel n must be a power of 2    ')

  !print *,'fourier transform...'
  call realft(inputArray1,1,fourierTransformed1)
  call realft(inputArray2,1,fourierTransformed2)

  fourierTransformed1(1) = cmplx(real(fourierTransformed1(1),kind=WP)*real(fourierTransformed2(1),kind=WP)/halfSize, &
                                 aimag(fourierTransformed1(1))*aimag(fourierTransformed2(1))/halfSize, kind=WP)
  fourierTransformed1(2:) = fourierTransformed1(2:) * conjg(fourierTransformed2(2:))/halfSize

  !print *,'back fourier transform...'
  call realft(mycorrel,-1,fourierTransformed1)

end function mycorrel

