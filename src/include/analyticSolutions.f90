!=====================================================================
!
!       m e m b r a n e S p h e r e
!       --------------------------------------------------
!      (c) 2025
!
!=====================================================================
!
! Description:
!   calculations for the analytic solutions by
!   Carl Tape (2003, chap. 3)
!
!-----------------------------------------------------------------------
module ray_orbits
!-----------------------------------------------------------------------
  use precisions
  implicit none
  integer:: ray_orbit
  real(WP):: ray_time,ray_cphase,ray_sigma,ray_epiDelta
  real(WP):: ray_sigma2half,ray_cphaseSquare,ray_twoPI,ray_PIby4
  real(WP):: ray_radiusbycSquared,ray_sineDelta
end module

!-----------------------------------------------------------------------
  function u_shape(colatitude,longitude,time,cphase,mu,ACCURACY)
!-----------------------------------------------------------------------
! calculation of the analytic solution for a source with only
! a given initial shape displacement (no time term)
!
! inputs:
!   colatitude, longitude - position of displacement location
!   time                       - given time step
!
! returns: displacement u_shape
  use precisions
  implicit none
  real(WP),intent(in):: colatitude,longitude,time,cphase,mu
  integer,intent(in):: ACCURACY
  real(WP):: u_shape
  ! local parameters
  integer:: l
  real(WP):: eigenfrequency
  real(WP):: imu
  real(WP):: dummy
  double precision,external:: legendrePolynomial
  real(WP),external:: integrand

  !print *,"u_shape: ", colatitude, longitude,time
  !integral boundaries

  u_shape = 0.0
  do l = 0,ACCURACY
    !get eigenfrequency for this degree l
    eigenfrequency = cphase*sqrt(l*(l+1.0))/EARTHRADIUS
    !print *,"eigenfrequency = ",eigenfrequency

    !get integral value
    call simpsons_integral_degree(integrand,l,0.0_WP,PI,imu,mu)

    !calculate displacement
    u_shape = u_shape + &
             (l+0.5) * imu * cos(eigenfrequency*time) &
             * legendrePolynomial(l,dcos(dble(colatitude)) )
  enddo

  ! to avoid compiler warnings
  dummy = longitude

  return
  end function


!-----------------------------------------------------------------------
  function u_f2(epiDelta,time,cphase,sigma,mu,ACCURACY)
!-----------------------------------------------------------------------
! calculation of the analytic solution for a source with only a forcing term 2
! (solution by Tape, 2003: eq. 3.34)
!
! inputs:
!   epiDelta                 - epicentral distance of receiver/position of
!                              displacement location (in radian)
!   time                     - given time step
!   cphase                   - phase velocity
!   sigma                    - timeparametersigma
!   mu                       - widthparametermu
!
! returns: displacement u_f2
  use precisions
  !use filterType, only: bw_waveperiod
  implicit none
  real(WP),intent(in):: epiDelta,time,cphase,sigma,mu
  integer,intent(in):: ACCURACY
  real(WP):: u_f2
  ! local parameters
  integer:: l,l_min,l_max
  real(WP):: eigenfrequency,sigma2,reciEarthradius,cphaseSquare
  real(WP):: freq_min,freq_max,center_freq,result,integralvalue
  real(WP):: termSigma,termIntegral,termSpreading
  real(WP),external:: integrand
  double precision,external:: legendrePolynomial
  ! timing
  !integer:: cpustart,cpuend,cpurate,cpumax
  !real:: cputime

  ! check to have double precision
  sigma2 = sigma*sigma
  reciEarthradius = 1.0d0/EARTHRADIUS
  cphaseSquare = cphase*cphase

  l_min = 0
  l_max = ACCURACY

  freq_min = 0.0
  freq_max = huge(freq_max)
  center_freq = 0.0

!      ! determine angular degree depending on desired frequency range
!      if (ONLY_FREQUENCY_RANGE) then
!        center_freq = 2.0*PI/WAVEPERIOD_L150
!        freq_min = center_freq - 2.0*PI*BW_HALFFREQUENCY
!        freq_max = center_freq + 2.0*PI*BW_HALFFREQUENCY
!
!        ! from w_l = c * sqrt( l(l+1) ) / r_earth it follows the solution for l
!        ! quadratic equation: l**2 + l - (w_l * r_earth / c )**2
!        !     -> l+ = -1/2 + 1/2*sqrt( 1+4( w_l * r_earth/c )**2
!        ! ( +0.5 for rounding up to integer)
!        l_min = int( -0.5 + 0.5*sqrt(1.0 + 4.0*(freq_min*EARTHRADIUS/cphase)**2) + 0.5)
!        l_max = int(-0.5 + 0.5*sqrt(1.0 + 4.0*(freq_max*EARTHRADIUS/cphase)**2) + 0.5)
!      endif

!      ! angular frequency range (output only at begining)
!      if (abs(time) <= 1.e-4) then
!        print *,'solution f2:'
!        print *,'    angular degrees min/max:',l_min,l_max
!        print *,'    source sigma/mu:',sngl(sigma),sngl(mu)
!        print *,'    epicentral distance:',sngl(epiDelta*180.0/PI)
!        print *,'      desired frequency range min/max:',sngl(freq_min),sngl(freq_max)
!        print *,'         period range min/max:',sngl(2*PI/freq_min),sngl(2*PI/freq_max)
!        print *,'      actual frequency range min/max:', &
!                sngl(cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS), &
!                sngl(cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS)
!        print *,'         period range min/max:', &
!                sngl(2*PI/(cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS)), &
!                sngl(2*PI/(cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS))
!        if (l_min > 0) then
!          print *,'    center frequency/period:',sngl(center_freq), &
!                      sngl(2*PI/center_freq)
!        endif
!      endif

  ! summation over angular degrees
  u_f2 = 0.0d0
  do l = l_min,l_max
    ! take time
    !call system_clock(COUNT=cpustart,COUNT_RATE=cpurate,COUNT_MAX=cpumax)

    !get eigenfrequency for this degree l
    eigenfrequency = cphase*sqrt(l*(l+1.0d0))*reciEarthradius

    !get integral value with integral boundaries [0,PI]
    call simpsons_integral_degree(integrand,l,0.0_WP,PI,integralvalue,mu)

    !calculate displacement (sum of Tape, 2003: eq. 3.34)
    termSigma = exp(-sigma2*(eigenfrequency**2)*0.5d0)
    termIntegral = cphaseSquare*(l+0.5d0)*integralvalue
    termSpreading = legendrePolynomial(l,dcos(dble(epiDelta)) )

    result = termIntegral*termSigma*termSpreading*cos(eigenfrequency*time)
    u_f2 = u_f2 + result

    ! skip if accuracy enough
    !print *,'u_f2 =',u_f2,l,result
    !if (mod(l,2)==0 .and. abs( result) < MACHINE_PRECISION) exit
    !call system_clock(cpuend)
    !cputime = real(cpuend-cpustart)/real(cpurate)
    !print *,'u_f2=',u_f2,l,cputime
    !if (abs(cputime) > 0.1) then
    !  print *,'  getting slower...'
    !  print *,'    u_f2 = ',u_f2,l,result
    !endif
    !if (abs(cputime) > 0.3) then
    !  print *,'  getting too slow...',l
    !  !print *,'    u_f2 = ',u_f2,l,result
    !  exit
    !endif
  enddo

  return
  end function


!-----------------------------------------------------------------------
  function integrand(degree,theta,mu)
!-----------------------------------------------------------------------
! integrand function defined by Tape, chapter 3, eq. (3.26), p. 23
!
! input:
!     degree  - degree of legendre polynomial to take
!     theta     - value to calculate at
!
! returns: integrand value
  use precisions
  implicit none
  integer,intent(in):: degree
  real(WP),intent(in):: theta,mu
  real(WP):: integrand
  ! local parameters
  real(WP):: mu2
  double precision,external:: legendrePolynomial

  mu2 = mu*mu
  integrand = legendrePolynomial(degree,dcos(dble(theta)) ) &
              *sin(theta)*exp(-theta*theta/(2.0d0*mu2))/mu2
  return
  end function


!-----------------------------------------------------------------------
  subroutine u_f2_asymptoticSolutions(epiDelta,time,cphase,sigma,mu, &
                                      ACCURACY,u_f2_asymptotic,u_f2_R1,u_f2_R2,u_f2_R3,u_f2_R4)
!-----------------------------------------------------------------------
! calculation of the analytic solution for a source with only a forcing term 2
! (solution by Tape, 2003: eq. 3.34)
! uses an asymptotic formulation for the legendre polynom -> ray theory formulism
!
! inputs:
!   epiDelta                 - epicentral distance of receiver/position of
!                              displacement location (in radian)
!   time                     - given time step
!   cphase                   - phase velocity
!   sigma                    - timeparametersigma
!   mu                       - widthparametermu
!
! returns: displacement u_f2
  use precisions
  use filterType, only: bw_waveperiod
  implicit none
  real(WP),intent(in):: epiDelta,time,cphase,sigma,mu
  integer,intent(in):: ACCURACY
  real(WP),intent(out):: u_f2_asymptotic,u_f2_R1,u_f2_R2,u_f2_R3,u_f2_R4
  ! local parameters
  integer:: l, l_min,l_max
  real(WP):: eigenfrequency,sigma2,reciEarthradius
  real(WP):: A_l,termSigma,termSpreading,termIntegral,integralvalue
  real(WP):: phase_t,phase_k,PI2,PIby4,cphaseSquare,lplus
  real(WP):: freq_min,freq_max,center_freq,result
  real(WP),external:: integrand,integrand_asymptotic
  ! timing
  integer:: cpustart,cpuend,cpurate,cpumax
  real:: cputime
  ! parameters
  logical,parameter:: ONLY_FREQUENCY_RANGE = .false.
  real(WP),parameter:: MACHINE_PRECISION = 1.e-8

  ! check to have double precision
  sigma2 = sigma*sigma
  reciEarthradius = 1.0/EARTHRADIUS
  PI2 = 2.0d0*PI
  PIby4 = PI/4.0d0

  l_min = 0
  l_max = ACCURACY
  freq_min = 0.0
  freq_max = huge(freq_max)
  center_freq = 0.0
  ! determine angular degree depending on desired frequency range
  if (ONLY_FREQUENCY_RANGE) then
    center_freq = 2.0*PI/bw_waveperiod
    freq_min = center_freq - 2.0*PI*BW_HALFFREQUENCY
    freq_max = center_freq + 2.0*PI*BW_HALFFREQUENCY

    ! from w_l = c * sqrt( l(l+1) ) / r_earth it follows the solution for l
    ! quadratic equation: l**2 + l - (w_l * r_earth / c )**2
    !     -> l+ = -1/2 + 1/2*sqrt( 1+4( w_l * r_earth/c )**2
    ! ( +0.5 for rounding up to integer)
    l_min = int( -0.5 + 0.5*sqrt(1.0 + 4.0*(freq_min*EARTHRADIUS/cphase)**2) + 0.5)
    l_max = int(-0.5 + 0.5*sqrt(1.0 + 4.0*(freq_max*EARTHRADIUS/cphase)**2) + 0.5)
  endif

  ! angular frequency range (output only at begining)
  if (abs(time) <= 1.e-4) then
  !  l = 0
  !  freq_min = cphase*sqrt(l*(l+1.0d0))*reciEarthradius
  !  l = ACCURACY
  !  freq_max = cphase*sqrt(l*(l+1.0d0))*reciEarthradius
    print *,'asymptotic solution f2:'
    print *,'    angular degrees min/max        :',l_min,'/',l_max
    print *,'    desired frequency range min/max:',freq_min,'/',freq_max
    print *,'               period range min/max:',sngl(2.d0*PI/freq_min),'/',sngl(2.d0*PI/freq_max)

    ! check
    print *,'    actual frequency range min/max:', &
            sngl(cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS),'/', &
            sngl(cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS)
    print *,'         period range min/max:', &
            sngl(2*PI/(cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS)),'/', &
            sngl(2*PI/(cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS))
    if (l_min > 0) then
      print *,'    center frequency/period:',center_freq,sngl(2.d0*PI/center_freq)
    endif
  endif

  cphaseSquare = cphase*cphase
  u_f2_asymptotic = 0.0d0
  u_f2_R1 = 0.0d0
  u_f2_R2 = 0.0d0
  u_f2_R3 = 0.0d0
  u_f2_R4 = 0.0d0
  do l = l_min,l_max
    ! take time
    call system_clock(COUNT=cpustart,COUNT_RATE=cpurate,COUNT_MAX=cpumax)

    !get eigenfrequency for this degree l
    eigenfrequency = cphase*sqrt(l*(l+1.0d0))*reciEarthradius

    !get integral value with integral boundaries [0,PI]
    call simpsons_integral_degree(integrand,l,0.0_WP,PI,integralvalue,mu)

    ! get asymptotic integral value
    ! does not work, gives NaN
    !call simpsons_integral_degree(integrand_asymptotic,l,0.0_WP,PI,integralvalue,mu)

    ! coefficients A
    lplus = l+0.5d0

    termSigma = exp(-sigma2*(eigenfrequency**2)*0.5d0)
    termSpreading = sqrt(2.0/(PI*lplus*sin(epiDelta)))
    termIntegral = cphaseSquare*lplus*integralvalue

    A_l = termIntegral*termSpreading*termSigma


    ! solution by Tape, 2003: eq. 3.34
    !calculate displacement (sum of Tape, 2003: eq. 3.34)
    phase_t = eigenfrequency*time
    phase_k = lplus*epiDelta

    result = A_l*cos(phase_t)*cos(phase_k - PIby4)
    u_f2_asymptotic = u_f2_asymptotic + result


    ! displacement for different orbits: first orbit R1, ...
    A_l = 0.5*A_l
    u_f2_R1 = u_f2_R1 + &
            A_l*cos(phase_t - phase_k + PIby4)
    u_f2_R2 = u_f2_R2 + &
            A_l*cos(phase_t + phase_k - PIby4)
    u_f2_R3 = u_f2_R3 + &
            A_l*cos(phase_t - phase_k + PIby4 + PI2)
    u_f2_R4 = u_f2_R4 + &
            A_l*cos(phase_t + phase_k - PIby4 + PI2)


    ! skip if accuracy enough
    !if (mod(l,2)==0 .and. abs( result) < MACHINE_PRECISION) exit
    call system_clock(cpuend)
    cputime = real(cpuend-cpustart)/real(cpurate)
    !print *,'u_f2=',u_f2,l,cputime
    !if (abs(cputime) > 0.1) then
    !  print *,'  getting slower...'
    !  print *,'    u_f2 = ',u_f2,l,result
    !endif
    if (abs(cputime) > 0.3) then
      !print *,'  getting too slow...',l
      !print *,'    u_f2_asymptotic = ',u_f2_asymptotic,l,result
      exit
    endif

  enddo

  !debug
  !print *,'u_f2_asymptotic =',u_f2_asymptotic
  !print *,'    integral = ',integralvalue
  !print *,'    phase velocity = ',cphase
  !print *,'    epicentral distance = ',epiDelta
  !print *,'    time = ',time

  end subroutine


!-----------------------------------------------------------------------
  function integrand_asymptotic(degree,theta,mu)
!-----------------------------------------------------------------------
! integrand function defined by Tape, chapter 3, eq. (3.26), p. 23 and
! using the asymptotic approach for the legendre polynomial:
! P_l (cos beta ) = sqrt( 2 / (PI * (l + 1/2) * sin beta ) ) * cos ( (l + 1/2) * beta - PI/4 )

! input:
!     degree  - degree of legendre polynomial to take
!     theta     - value to calculate at
!
! returns: integrand value
!
!
! attention: may be not defined, gives NaN
!
  use precisions
  implicit none
  integer,intent(in):: degree
  real(WP),intent(in):: theta,mu
  real(WP):: integrand_asymptotic
  ! local parameters
  real(WP):: mu2
  integer:: l
  !real(WP), parameter :: PI = 3.1415926535897931d0

  !print *,"asymptotic integrand function for ", degree, theta

  mu2 = mu*mu
  l = degree

  integrand_asymptotic =  sqrt(2.d0/(PI*(l+0.5) * sin(theta))) &
                             * cos( (l+0.5)*theta - PI/4.d0 ) &
                             * sin(theta)*exp(-theta**2/(2.d0*mu2))/mu2

  return
  end function


!-----------------------------------------------------------------------
  function u_f2_orbit(orbit,epiDelta,time,cphase,sigma,mu,ACCURACY)
!-----------------------------------------------------------------------
! calculation of the analytic solution for a source with only a forcing term 2
! (solution by Tape, 2003: eq. 3.34)
! uses an asymptotic formulation for the legendre polynom -> ray theory formulism
! and takes only first arrival solution R1
!
! inputs:
!   epiDelta                 - epicentral distance of receiver/position of
!                              displacement location (in radian)
!   time                     - given time step
!   cphase                   - phase velocity
!   sigma                    - timeparametersigma
!   mu                       - widthparametermu
!
! returns: displacement u_f2
  use precisions
  use filterType, only: bw_waveperiod
  use splinefunction
  use ray_orbits
  implicit none
  integer,intent(in):: orbit
  real(WP),intent(in):: epiDelta,time,cphase,sigma,mu
  integer,intent(in):: ACCURACY
  real(WP):: u_f2_orbit
  ! local parameters
  integer:: ier,length,i,icount
  integer:: l,l_min,l_max
  double precision:: degreel
  real(WP):: integralvalue
  real(WP):: freq_min,freq_max,center_freq,result
  real(WP),external:: u_f2_asymptoticdisplacement
  real(WP),external:: integrand
  double precision,external:: cubicspline_eval

  real(WP):: lambda,omega_k
  real(WP):: dummy
  double precision:: deltal
  ! timing
  !integer:: cpustart,cpuend,cpurate,cpumax
  !real:: cputime

  logical:: do_file_out

  real(WP),external:: romberg_integral
  real(WP),external:: simpsons_integral
  ! filter in frequency range
  logical,parameter:: ONLY_FREQUENCY_RANGE = .false.
  integer,parameter:: DEGREE_SHIFT = 0
  logical,parameter:: DO_BY_SUMMATION = .false.

  ! angular degree
  l_min = 0
  l_max = ACCURACY

  ! set parameters
  ray_orbit = orbit
  ray_time = time

  ! determine angular degree depending on desired frequency range
  if (ONLY_FREQUENCY_RANGE) then
    center_freq = 2.0*PI/bw_waveperiod
    freq_min = center_freq - 2.0*PI*BW_HALFFREQUENCY
    freq_max = center_freq + 2.0*PI*BW_HALFFREQUENCY

    ! from w_l = c * sqrt( l(l+1) ) / r_earth it follows the solution for l
    ! quadratic equation: l**2 + l - (w_l * r_earth / c )**2
    !     -> l+ = -1/2 + 1/2*sqrt( 1+4( w_l * r_earth/c )**2
    ! ( +0.5 for rounding up to integer)
    l_min = int(-0.5 + 0.5*sqrt(1.0 + 4.0*(freq_min*EARTHRADIUS/cphase)**2) + 0.5)
    l_max = int(-0.5 + 0.5*sqrt(1.0 + 4.0*(freq_max*EARTHRADIUS/cphase)**2) + 0.5)
  endif

  ! angular frequency range (output only at begining)
  if (abs(time) <= 1.e-8 .and. orbit == 1) then
    ! determine frequencies
    freq_min = cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS
    freq_max = cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS !huge(freq_max)

    print *,'  asymptotic solution f2 for R1:'
    print *,'    angular degrees min/max:',l_min,'/',l_max
    print *,'    frequency min/max      :',freq_min,'/',freq_max
    print *,'       period min/max      :',sngl(2.d0*PI/freq_min),sngl(2.d0*PI/freq_max)
  endif

  if (.not. allocated(X)) then
    ! for spline representation
    length = l_max - l_min + 1

    ! initialize
    i1 = 1
    i2 = length
    allocate(X(i2),Y(i2),Q(3,i2),F(3,i2),stat=ier)
    if (ier /= 0) call stopProgram('u_f2_orbit - error allocating spline arrays ')

    ! get spline representation
    icount = 0
    open(IOUT,file='OUTPUT/integrals_discretel.dat')
    write(IOUT,*) '#format: degree_l integralvalue'
    do_file_out = .true.

    do i = l_min,l_max
      ! starting degrees from 0 to degreeaccuracy, but array indices from 1 to length
      icount = icount + 1
      X(icount) = dble(i)

      ! get integral value with integral boundaries [0,PI]
      call simpsons_integral_degree(integrand,i,0.0_WP,PI,integralvalue,mu)
      Y(icount) = dble(integralvalue)

      if (do_file_out) write(IOUT,*) i,Y(icount)
    enddo
    close(IOUT)

    ! calculates spline arrays
    call cubicspline_setup()

    ! test
    open(IOUT,file='OUTPUT/integrals_nonintegerl.test.dat')
    write(IOUT,*) '#format: degree_l integralvalue'

    ! old style
    !do degreel = l_min,l_max,0.2d0
    ! new style with loop over integers
    do l = int(l_min/0.2d0),int(l_max/0.2d0)
      degreel = l * 0.2d0
      integralvalue = cubicspline_eval(degreel)
      write(IOUT,*) sngl(degreel),integralvalue
    enddo
    close(IOUT)

    ! debug output
    open(88,file='OUTPUT/integrals_omega.dat')
    write(88,*) '#format: degree_l integralvalue'
  else
    do_file_out = .false.
  endif

  ! take time
  !call system_clock(COUNT=cpustart,COUNT_RATE=cpurate,COUNT_MAX=cpumax)

  ! summation or integration
  if (DO_BY_SUMMATION) then
    ! summation
    deltal = 0.5d0
    u_f2_orbit = 0.0d0
    ! old style
    !do degreel = l_min,l_max,deltal
    ! new style with loop over integers
    do l = int(l_min/deltal),int(l_max/deltal)
      degreel = l * deltal
      ! determines frequency
      lambda = degreel + 0.5d0
      omega_k = cphase * sqrt(lambda**2 - 0.25d0) / EARTHRADIUS

      ! calculates solution
      result = u_f2_asymptoticdisplacement(omega_k)
      u_f2_orbit = u_f2_orbit + result*deltal

      if (do_file_out) write(88,*) omega_k,u_f2_orbit
    enddo
    if (do_file_out) close(88)
  else
    ! integration
    l_min = l_min + DEGREE_SHIFT
    freq_min = cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS
    freq_max = cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS

    ! integration
    u_f2_orbit = simpsons_integral(u_f2_asymptoticdisplacement,dble(freq_min),dble(freq_max))
    ! qromb integration
    !u_f2_orbit = romberg_integral(u_f2_asymptoticdisplacement,dble(freq_min),dble(freq_max))

    u_f2_orbit = u_f2_orbit*EARTHRADIUS/cphase

    !print *,'u_f2_orbit:',l_max,cphase,time
    !print *,'    frequencies:',freq_min,freq_max
    !print *,'    simp:',u_f2_orbit
    if (do_file_out) write(88,*) freq_max,u_f2_orbit
    if (do_file_out) close(88)
  endif

  ! skip if accuracy enough
  !call system_clock(cpuend)
  !cputime = real(cpuend-cpustart)/real(cpurate)
  !if (abs(cputime) > 0.3) then
  !  print *,'  getting too slow...',l
  !  print *,'    u_f2_orbit = ',u_f2_orbit
  !  stop
  !endif

  !debug
  !print *,'u_f2_R1 =',u_f2_R1
  !print *,'    integral = ',integralvalue
  !print *,'    phase velocity = ',cphase
  !print *,'    epicentral distance = ',epiDelta
  !print *,'    time = ',time

  ! to avoid compiler warnings
  dummy = epiDelta
  dummy = sigma

  return
  end function


!-----------------------------------------------------------------------
  function u_f2_asymptoticdisplacement(omega)
!-----------------------------------------------------------------------
  use precisions; use splinefunction
  use ray_orbits
  implicit none
  real(WP),intent(in):: omega
  real(WP):: u_f2_asymptoticdisplacement
  ! local parameters
  integer:: orbit,iorb
  real(WP):: epiDelta,time,omegaSquare,sineDelta,radiusbycSquared
  real(WP):: sigma2half,integralvalue,lambda,termSigma, &
            termSpreading,termIntegral,cphaseSquare,twoPI,PIby4
  double precision:: degreel
  complex(WP) :: A_lplus,A_lminus,A_even,A_odd,A_l
  double precision,external:: cubicspline_eval

  logical:: ODD_ORBIT
  complex(WP),parameter:: II = ( 0.0_WP, 1.0_WP )

  ! use ray parameters
  orbit = ray_orbit
  epiDelta = ray_epiDelta
  time = ray_time
  sigma2half = ray_sigma2half
  cphaseSquare = ray_cphaseSquare
  twoPI = ray_twoPI
  PIby4 = ray_PIby4
  radiusbycSquared = ray_radiusbycSquared
  sineDelta = ray_sineDelta

  ! check
  if (orbit < 1) stop "u_f2_asymptoticdisplacement - wrong orbit"

  ! get corresponding lambda and degree
  omegaSquare = omega**2
  lambda = sqrt(omegaSquare*radiusbycSquared + 0.25)
  degreel = lambda - 0.5d0

  ! get integral value with integral boundaries [0,PI]
  integralvalue = cubicspline_eval(degreel)

  ! coefficients A
  termSigma = exp(-omegaSquare*sigma2half)
  termSpreading = sqrt( 1.0/(twoPI*lambda*sineDelta))
  termIntegral = cphaseSquare*lambda*integralvalue

  !!! until now it is approximated, if we want to use it without approximation l >> 1
  !!! we have to add this term
  !!! termIntegral = termIntegral * sqrt(1.0 - 1.0/(4*lambda*lambda))

  ! orbit adding: epiDelta + 2PI + 2PI + ...
  if (mod(orbit,2) == 1) then
    ! iorbit = 0,1,2,3,4,.. for R1,R3,R5,R7,...
    iorb = (orbit-1)/2
    ODD_ORBIT = .true.
  else
    ! iorbit = 1,2,3,4,.. for R2,R4,R6,...
    iorb = orbit/2
    ODD_ORBIT = .false.
  endif
  if (ODD_ORBIT) then
    ! R1,R3,.. orbit
    A_lminus = termSigma*termSpreading*termIntegral*exp(+II*PIby4)
    A_odd = (-1)**iorb*A_lminus*exp(-II*lambda*(twoPI*iorb+epiDelta) )
    A_l = A_odd
  else
    ! R2,R4,.. orbit
    A_lplus = termSigma*termSpreading*termIntegral*exp(-II*PIby4)
    A_even = (-1)**iorb*A_lplus*exp(-II*lambda*(twoPI*iorb-epiDelta) )
    A_l = A_even
  endif

  u_f2_asymptoticdisplacement = real( A_l*exp(II*omega*time) )

  return
  end function

!-----------------------------------------------------------------------
  function integrand_cont(degree,theta,mu)
!-----------------------------------------------------------------------
! integrand function defined by Tape, chapter 3, eq. (3.26), p. 23
!
! input:
!     degree  - degree of legendre polynomial to take
!     theta     - value to calculate at
!
! returns: integrand value
  use precisions
  implicit none
  real(WP),intent(in):: degree
  real(WP),intent(in):: theta,mu
  real(WP):: integrand_cont
  ! local parameters
  double precision,external:: legendrePolynomial_cont
  real(WP):: mu2

  mu2 = mu*mu
  integrand_cont = legendrePolynomial_cont(degree,cos(theta)) &
              *sin(theta)*exp(-theta*theta/(2.0d0*mu2))/mu2

  return
  end function


!-----------------------------------------------------------------------
  function u_f2_R2(epiDelta,time,cphase,sigma,mu,ACCURACY)
!-----------------------------------------------------------------------
! calculation of the analytic solution for a source with only a forcing term 2
! (solution by Tape, 2003: eq. 3.34)
! uses an asymptotic formulation for the legendre polynom -> ray theory formulism
! and takes only first arrival solution R1
!
! inputs:
!   epiDelta                 - epicentral distance of receiver/position of
!                              displacement location (in radian)
!   time                     - given time step
!   cphase                   - phase velocity
!   sigma                    - timeparametersigma
!   mu                       - widthparametermu
!
! returns: displacement u_f2
  use precisions
  use filterType, only: bw_waveperiod
  implicit none
  real(WP),intent(in):: epiDelta,time,cphase,sigma,mu
  integer,intent(in):: ACCURACY
  real(WP):: u_f2_R2
  ! local parameters
  integer:: l, l_min,l_max
  real(WP):: integralvalue
  real(WP):: eigenfrequency,sigma2,reciEarthradius
  real(WP):: freq_min,freq_max,center_freq, A_l
  real(WP),external:: integrand,integrand_asymptotic
  logical,parameter:: ONLY_FREQUENCY_RANGE = .false.

  ! check to have double precision
  sigma2 = sigma*sigma
  reciEarthradius = 1.0/EARTHRADIUS
  l_min = 0
  l_max = ACCURACY
  freq_min = 0.0
  freq_max = huge(freq_max)
  center_freq = 0.0
  ! determine angular degree depending on desired frequency range
  if (ONLY_FREQUENCY_RANGE) then
    center_freq = 2.0*PI/bw_waveperiod
    freq_min = center_freq - 2.0*PI*BW_HALFFREQUENCY
    freq_max = center_freq + 2.0*PI*BW_HALFFREQUENCY

    ! from w_l = c * sqrt( l(l+1) ) / r_earth it follows the solution for l
    ! quadratic equation: l**2 + l - (w_l * r_earth / c )**2
    !     -> l+ = -1/2 + 1/2*sqrt( 1+4( w_l * r_earth/c )**2
    ! ( +0.5 for rounding up to integer)
    l_min = int( -0.5 + 0.5*sqrt(1.0 + 4.0*(freq_min*EARTHRADIUS/cphase)**2) + 0.5)
    l_max = int(-0.5 + 0.5*sqrt(1.0 + 4.0*(freq_max*EARTHRADIUS/cphase)**2) + 0.5)
  endif

  ! angular frequency range (output only at begining)
  if (abs(time) <= 1.e-4) then
  !  l = 0
  !  freq_min = cphase*sqrt(l*(l+1.0d0))*reciEarthradius
  !  l = DEGREEACCURACY
  !  freq_max = cphase*sqrt(l*(l+1.0d0))*reciEarthradius
    print *,'asymptotic solution f2 for R2:'
    print *,'    angular degrees min/max        :',l_min,'/',l_max
    print *,'    desired frequency range min/max:',freq_min,'/',freq_max
    print *,'               period range min/max:',sngl(2.d0*PI/freq_min),'/',sngl(2.d0*PI/freq_max)

    ! check
    print *,'    actual frequency range min/max:', &
            sngl(cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS),'/', &
            sngl(cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS)
    print *,'         period range min/max:', &
            sngl(2*PI/(cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS)),'/', &
            sngl(2*PI/(cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS))
    if (l_min > 0) then
      print *,'    center frequency/period:',center_freq,sngl(2.d0*PI/center_freq)
    endif
  endif

  u_f2_R2 = 0.0d0
  do l = l_min,l_max
    !get eigenfrequency for this degree l
    eigenfrequency = cphase*sqrt(l*(l+1.0d0))*reciEarthradius

    !get integral value with integral boundaries [0,PI]
    call simpsons_integral_degree(integrand,l,0.0_WP,PI,integralvalue,mu)

    ! asymptotic gives NaN
    !call simpsons_integral_degree(integrand_asymptotic,l,0.0_WP,PI,integralvalue,mu)

    !calculate displacement (sum of Tape, 2003: eq. 3.34, see notes)
    ! coefficients A
    A_l = cphase*cphase*(l+0.5)*integralvalue &
          *sqrt(2.0/(PI*(l+0.5)*sin(epiDelta))) &
          *exp(-sigma2*(eigenfrequency**2)*0.5d0)

    ! displacement
    u_f2_R2 = u_f2_R2 + &
            A_l/2.0*cos(eigenfrequency*time + (l+0.5)*epiDelta - PI/4.0)
  enddo

  !debug
  !print *,'u_f2_R2 =',u_f2_R2
  !print *,'    integral = ',integralvalue
  !print *,'    phase velocity = ',cphase
  !print *,'    epicentral distance = ',epiDelta
  !print *,'    time = ',time

  return
  end function

!-----------------------------------------------------------------------
  function u_f2_R3(epiDelta,time,cphase,sigma,mu,ACCURACY)
!-----------------------------------------------------------------------
! calculation of the analytic solution for a source with only a forcing term 2
! (solution by Tape, 2003: eq. 3.34)
! uses an asymptotic formulation for the legendre polynom -> ray theory formulism
! and takes only first arrival solution R1
!
! inputs:
!   epiDelta                 - epicentral distance of receiver/position of
!                              displacement location (in radian)
!   time                     - given time step
!   cphase                   - phase velocity
!   sigma                    - timeparametersigma
!   mu                       - widthparametermu
!
! returns: displacement u_f2
  use precisions
  use filterType, only: bw_waveperiod
  implicit none
  real(WP),intent(in):: epiDelta,time,cphase,sigma,mu
  integer,intent(in):: ACCURACY
  real(WP):: u_f2_R3
  ! local parameters
  integer:: l, l_min,l_max
  real(WP):: integralvalue
  real(WP):: eigenfrequency,sigma2,reciEarthradius
  real(WP):: freq_min,freq_max,center_freq, A_l
  real(WP),external:: integrand,integrand_asymptotic
  logical,parameter:: ONLY_FREQUENCY_RANGE = .false.

  ! check to have double precision
  sigma2 = sigma*sigma
  reciEarthradius = 1.0/EARTHRADIUS
  l_min = 0
  l_max = ACCURACY
  freq_min = 0.0
  freq_max = huge(freq_max)
  center_freq = 0.0
  ! determine angular degree depending on desired frequency range
  if (ONLY_FREQUENCY_RANGE) then
    center_freq = 2.0*PI/bw_waveperiod
    freq_min = center_freq - 2.0*PI*BW_HALFFREQUENCY
    freq_max = center_freq + 2.0*PI*BW_HALFFREQUENCY

    ! from w_l = c * sqrt( l(l+1) ) / r_earth it follows the solution for l
    ! quadratic equation: l**2 + l - (w_l * r_earth / c )**2
    !     -> l+ = -1/2 + 1/2*sqrt( 1+4( w_l * r_earth/c )**2
    ! ( +0.5 for rounding up to integer)
    l_min = int( -0.5 + 0.5*sqrt(1.0 + 4.0*(freq_min*EARTHRADIUS/cphase)**2) + 0.5)
    l_max = int(-0.5 + 0.5*sqrt(1.0 + 4.0*(freq_max*EARTHRADIUS/cphase)**2) + 0.5)
  endif

  ! angular frequency range (output only at begining)
  if (abs(time) <= 1.e-4) then
  !  l = 0
  !  freq_min = cphase*sqrt(l*(l+1.0d0))*reciEarthradius
  !  l = DEGREEACCURACY
  !  freq_max = cphase*sqrt(l*(l+1.0d0))*reciEarthradius
    print *,'asymptotic solution f2 for R3:'
    print *,'    angular degrees min/max        :',l_min,'/',l_max
    print *,'    desired frequency range min/max:',freq_min,'/',freq_max
    print *,'               period range min/max:',sngl(2.d0*PI/freq_min),'/',sngl(2.d0*PI/freq_max)

    ! check
    print *,'    actual frequency range min/max:', &
            sngl(cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS),'/', &
            sngl(cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS)
    print *,'         period range min/max:', &
            sngl(2*PI/(cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS)),'/', &
            sngl(2*PI/(cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS))
    if (l_min > 0) then
      print *,'    center frequency/period:',center_freq,sngl(2.d0*PI/center_freq)
    endif
  endif

  u_f2_R3 = 0.0d0
  do l = l_min,l_max
    !get eigenfrequency for this degree l
    eigenfrequency = cphase*sqrt(l*(l+1.0d0))*reciEarthradius

    !get integral value with integral boundaries [0,PI]
    call simpsons_integral_degree(integrand,l,0.0_WP,PI,integralvalue,mu)

    ! asymptotic gives NaN
    !call simpsons_integral_degree(integrand_asymptotic,l,0.0_WP,PI,integralvalue,mu)

    !calculate displacement (sum of Tape, 2003: eq. 3.34, see notes)
    ! coefficients A
    A_l = cphase*cphase*(l+0.5)*integralvalue &
          *sqrt(2.0/(PI*(l+0.5)*sin(epiDelta))) &
          *exp(-sigma2*(eigenfrequency**2)*0.5d0)

    ! displacement
    u_f2_R3 = u_f2_R3 + &
        A_l/2.0*cos(eigenfrequency*time - (l+0.5)*epiDelta + PI/4.0 + 2.0*PI)
  enddo

  !debug
  !print *,'u_f2_R3 =',u_f2_R3
  !print *,'    integral = ',integralvalue
  !print *,'    phase velocity = ',cphase
  !print *,'    epicentral distance = ',epiDelta
  !print *,'    time = ',time
  return
  end function


!-----------------------------------------------------------------------
  function u_f2_R4(epiDelta,time,cphase,sigma,mu,ACCURACY)
!-----------------------------------------------------------------------
! calculation of the analytic solution for a source with only a forcing term 2
! (solution by Tape, 2003: eq. 3.34)
! uses an asymptotic formulation for the legendre polynom -> ray theory formulism
! and takes only first arrival solution R1
!
! inputs:
!   epiDelta                 - epicentral distance of receiver/position of
!                              displacement location (in radian)
!   time                     - given time step
!   cphase                   - phase velocity
!   sigma                    - timeparametersigma
!   mu                       - widthparametermu
!
! returns: displacement u_f2
  use precisions
  use filterType, only: bw_waveperiod
  implicit none
  real(WP),intent(in):: epiDelta,time,cphase,sigma,mu
  integer,intent(in):: ACCURACY
  real(WP):: u_f2_R4
  ! local parameters
  integer:: l, l_min,l_max
  real(WP):: integralvalue
  real(WP):: eigenfrequency,sigma2,reciEarthradius
  real(WP):: freq_min,freq_max,center_freq, A_l
  real(WP),external:: integrand,integrand_asymptotic
  logical,parameter:: ONLY_FREQUENCY_RANGE = .false.

  ! check to have double precision
  sigma2 = sigma*sigma
  reciEarthradius = 1.0/EARTHRADIUS
  l_min = 0
  l_max = ACCURACY
  freq_min = 0.0
  freq_max = huge(freq_max)
  center_freq = 0.0

  ! determine angular degree depending on desired frequency range
  if (ONLY_FREQUENCY_RANGE) then
    center_freq = 2.0*PI/bw_waveperiod
    freq_min = center_freq - 2.0*PI*BW_HALFFREQUENCY
    freq_max = center_freq + 2.0*PI*BW_HALFFREQUENCY

    ! from w_l = c * sqrt( l(l+1) ) / r_earth it follows the solution for l
    ! quadratic equation: l**2 + l - (w_l * r_earth / c )**2
    !     -> l+ = -1/2 + 1/2*sqrt( 1+4( w_l * r_earth/c )**2
    ! ( +0.5 for rounding up to integer)
    l_min = int( -0.5 + 0.5*sqrt(1.0 + 4.0*(freq_min*EARTHRADIUS/cphase)**2) + 0.5)
    l_max = int(-0.5 + 0.5*sqrt(1.0 + 4.0*(freq_max*EARTHRADIUS/cphase)**2) + 0.5)
  endif

  ! angular frequency range (output only at begining)
  if (abs(time) <= 1.e-4) then
  !  l = 0
  !  freq_min = cphase*sqrt(l*(l+1.0d0))*reciEarthradius
  !  l = ACCURACY
  !  freq_max = cphase*sqrt(l*(l+1.0d0))*reciEarthradius
    print *,'asymptotic solution f2 for R4:'
    print *,'    angular degrees min/max        :',l_min,'/',l_max
    print *,'    desired frequency range min/max:',freq_min,'/',freq_max
    print *,'               period range min/max:',sngl(2.d0*PI/freq_min),'/',sngl(2.d0*PI/freq_max)

    ! check
    print *,'    actual frequency range min/max:', &
            sngl(cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS),'/', &
            sngl(cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS)
    print *,'         period range min/max:', &
            sngl(2*PI/(cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS)),'/', &
            sngl(2*PI/(cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS))
    if (l_min > 0) then
      print *,'    center frequency/period:',center_freq,sngl(2.d0*PI/center_freq)
    endif
  endif

  u_f2_R4 = 0.0d0
  do l = l_min,l_max
    !get eigenfrequency for this degree l
    eigenfrequency = cphase*sqrt(l*(l+1.0d0))*reciEarthradius

    !get integral value with integral boundaries [0,PI]
    call simpsons_integral_degree(integrand,l,0.0_WP,PI,integralvalue,mu)

    ! asymptotic gives NaN
    !call simpsons_integral_degree(integrand_asymptotic,l,0.0_WP,PI,integralvalue,mu)

    !calculate displacement (sum of Tape, 2003: eq. 3.34, see notes)
    ! coefficients A
    A_l = cphase*cphase*(l+0.5)*integralvalue &
          *sqrt(2.0/(PI*(l+0.5)*sin(epiDelta))) &
          *exp(-sigma2*(eigenfrequency**2)*0.5d0)

    ! displacement
    u_f2_R4 = u_f2_R4 + &
        A_l/2.0*cos(eigenfrequency*time + (l+0.5)*epiDelta - PI/4.0 + 2.0*PI)
  enddo

  !debug
  !print *,'u_f2_R4 =',u_f2_R4
  !print *,'    integral = ',integralvalue
  !print *,'    phase velocity = ',cphase
  !print *,'    epicentral distance = ',epiDelta
  !print *,'    time = ',time
  return
  end function

