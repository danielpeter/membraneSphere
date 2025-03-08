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
!
! Description:
!   calculations for the analytic solutions by
!   Carl Tape (2003, chap. 3)
!
module ray_orbits
  use precisions
  integer:: ray_orbit
  real(WP):: ray_time,ray_cphase,ray_sigma,ray_epiDelta
  real(WP):: ray_sigma2half,ray_cphaseSquare,ray_twoPI,ray_PIby4
  real(WP):: ray_radiusbycSquared,ray_sineDelta
end module

!-----------------------------------------------------------------------
      function u_shape(colatitude,longitude,time,cphase, &
                                       mu,ACCURACY)
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
      real(WP):: u_shape
      real(WP),intent(in):: colatitude,longitude,time,cphase,mu
      integer,intent(in):: ACCURACY
      real(WP):: eigenfrequency
      real(WP):: imu
      double precision,external:: legendrePolynomial
      real(WP):: fromA,toB
      integer:: l
      real(WP),external:: integrand

      !print *,"u_shape: ", colatitude, longitude,time
      !integral boundaries

      u_shape = 0.0
      do l = 0,ACCURACY
        !get eigenfrequency for this degree l
        eigenfrequency= cphase*sqrt(l*(l+1.0))/EARTHRADIUS
        !print *,"eigenfrequency = ",eigenfrequency

        !get integral value
        call qsimp_degree(integrand,l,0.0d0,PI,imu,mu)

        !calculate displacement
        u_shape= u_shape + &
                 (l+0.5)*imu*cos(eigenfrequency*time) &
                 *legendrePolynomial(l,dcos(dble(colatitude)) )
      enddo

      return
      end

!-----------------------------------------------------------------------
      function u_f2(epiDelta,time,cphase,sigma,mu,ACCURACY)
!-----------------------------------------------------------------------
! calculation of the analytic solution for a source with only a forcing term 2
! (solution by Tape, 2003: eq. 3.34)
!
! inputs:
!   epiDelta                  - epicentral distance of receiver/position of displacement
!                                    location (in radian)
!   time                       - given time step
!   cphase                  - phase velocity
!   sigma                     - timeparametersigma
!   mu                          - widthparametermu
!
! returns: displacement u_f2
      use precisions
      use filterType, only: bw_waveperiod
      implicit none
      real(WP):: u_f2
      real(WP),intent(in):: epiDelta,time,cphase,sigma,mu
      integer,intent(in):: ACCURACY
      integer:: l,l_min,l_max
      real(WP):: eigenfrequency,sigma2,reciEarthradius,cphaseSquare
      real(WP):: freq_min,freq_max,center_freq,result,integralvalue
      real(WP):: termSigma,termIntegral,termSpreading
      real(WP),external:: integrand
      double precision,external:: legendrePolynomial
      integer:: cpustart,cpuend,cpurate,cpumax
      real:: cputime

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
!      if ( ONLY_FREQUENCY_RANGE ) then
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
!      if ( abs(time) <= 1.e-4) then
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
!        if ( l_min > 0 ) then
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
        eigenfrequency= cphase*sqrt(l*(l+1.0d0))*reciEarthradius

        !get integral value with integral boundaries [0,PI]
        call qsimp_degree(integrand,l,0.0,PI,integralvalue,mu)

        !calculate displacement (sum of Tape, 2003: eq. 3.34)
        termSigma = exp(-sigma2*(eigenfrequency**2)*0.5d0)
        termIntegral = cphaseSquare*(l+0.5d0)*integralvalue
        termSpreading = legendrePolynomial(l,dcos(dble(epiDelta)) )

        result= termIntegral*termSigma*termSpreading*cos(eigenfrequency*time)
        u_f2 = u_f2 + result

        ! skip if accuracy enough
        !print *,'u_f2 =',u_f2,l,result
        !if ( mod(l,2)==0 .and. abs( result ) < MACHINE_PRECISION) exit
        !call system_clock(cpuend)
        !cputime = real(cpuend-cpustart)/real(cpurate)
        !print *,'u_f2=',u_f2,l,cputime
        !if ( abs(cputime) > 0.1 ) then
        !  print *,'  getting slower...'
        !  print *,'    u_f2 = ',u_f2,l,result
        !endif
        !if ( abs(cputime) > 0.3 ) then
        !  print *,'  getting too slow...',l
        !  !print *,'    u_f2 = ',u_f2,l,result
        !  exit
        !endif
      enddo

      return
      end

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
      real(WP):: integrand
      integer,intent(in):: degree
      real(WP),intent(in):: theta,mu
      double precision,external:: legendrePolynomial
      real(WP):: mu2

      mu2 = mu*mu
      integrand = legendrePolynomial(degree,dcos(dble(theta)) ) &
                  *sin(theta)*exp(-theta*theta/(2.0d0*mu2))/mu2
      return
      end


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
!                             displacement location (in radian)
!   time                       - given time step
!   cphase                  - phase velocity
!   sigma                     - timeparametersigma
!   mu                          - widthparametermu
!
! returns: displacement u_f2
      use precisions
      use filterType, only: bw_waveperiod
      implicit none
      real(WP),intent(in):: epiDelta,time,cphase,sigma,mu
      integer,intent(in):: ACCURACY
      real(WP),intent(out):: u_f2_asymptotic,u_f2_R1,u_f2_R2,u_f2_R3,u_f2_R4
      integer:: l, l_min,l_max
      real(WP):: eigenfrequency,sigma2,reciEarthradius
      real(WP):: A_l,termSigma,termSpreading,termIntegral,integralvalue
      real(WP):: phase_t,phase_k,PI2,PIby4,cphaseSquare,lplus
      real(WP):: freq_min,freq_max,center_freq,result
      real(WP),external:: integrand,integrand_asymptotic
      integer:: cpustart,cpuend,cpurate,cpumax
      real:: cputime
      ! parameters
      logical,parameter:: ONLY_FREQUENCY_RANGE         = .false.
      real(WP), parameter :: MACHINE_PRECISION        = 1.e-8

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
      if ( ONLY_FREQUENCY_RANGE ) then
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
      if ( abs(time) <= 1.e-4) then
      !  l = 0
      !  freq_min = cphase*sqrt(l*(l+1.0d0))*reciEarthradius
      !  l = ACCURACY
      !  freq_max = cphase*sqrt(l*(l+1.0d0))*reciEarthradius
        print *,'asymptotic solution f2:'
        print *,'    angular degrees min/max:',l_min,l_max
        print *,'    desired frequency range min/max:',sngl(freq_min),sngl(freq_max)
        print *,'         period range min/max:',sngl(2*PI/freq_min),sngl(2*PI/freq_max)

        ! check
        print *,'    actual frequency range min/max:', &
                sngl(cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS), &
                sngl(cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS)
        print *,'         period range min/max:', &
                sngl(2*PI/(cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS)), &
                sngl(2*PI/(cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS))
        if ( l_min > 0 ) then
          print *,'    center frequency/period:',sngl(center_freq), &
                      sngl(2*PI/center_freq)
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
        eigenfrequency= cphase*sqrt(l*(l+1.0d0))*reciEarthradius

        !get integral value with integral boundaries [0,PI]
        call qsimp_degree(integrand,l,0.0,PI,integralvalue,mu)

        ! get asymptotic integral value
        ! does not work, gives NaN
        !call qsimp_degree(integrand_asymptotic,l,0.0d0,PI,integralvalue,mu)

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
        !if ( mod(l,2)==0 .and. abs( result ) < MACHINE_PRECISION) exit
        call system_clock(cpuend)
        cputime = real(cpuend-cpustart)/real(cpurate)
        !print *,'u_f2=',u_f2,l,cputime
        !if ( abs(cputime) > 0.1 ) then
        !  print *,'  getting slower...'
        !  print *,'    u_f2 = ',u_f2,l,result
        !endif
        if ( abs(cputime) > 0.3 ) then
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

      return
      end

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
      real(WP):: integrand_asymptotic
      integer,intent(in):: degree
      real(WP),intent(in):: theta,mu
      real(WP):: mu2
      integer:: l
      !real(WP), parameter :: PI = 3.1415926535897931d0

      !print *,"asymptotic integrand function for ", degree, theta

      mu2 = mu*mu
      l = degree

      integrand_asymptotic =  sqrt(2/(PI*(l+0.5)*sin(theta))) &
                  *cos( (l+0.5)*theta - PI/4.0 ) &
                  *sin(theta)*exp(-theta**2/(2.0d0*mu2))/mu2
      return
      end


!-----------------------------------------------------------------------
      function u_f2_orbit(orbit,epiDelta,time,cphase,sigma,mu, &
                                              ACCURACY)
!-----------------------------------------------------------------------
! calculation of the analytic solution for a source with only a forcing term 2
! (solution by Tape, 2003: eq. 3.34)
! uses an asymptotic formulation for the legendre polynom -> ray theory formulism
! and takes only first arrival solution R1
!
! inputs:
!   epiDelta                 - epicentral distance of receiver/position of
!                             displacement location (in radian)
!   time                       - given time step
!   cphase                  - phase velocity
!   sigma                     - timeparametersigma
!   mu                          - widthparametermu
!
! returns: displacement u_f2
      use precisions
      use filterType, only: bw_waveperiod
      use splinefunction
      use ray_orbits
      implicit none
      real(WP):: u_f2_orbit
      integer,intent(in):: orbit
      real(WP),intent(in):: epiDelta,time,cphase,sigma,mu
      integer,intent(in):: ACCURACY
      integer:: ierror,length,i,icount
      real(WP):: degreel,l_min,l_max
      real(WP):: integralvalue
      real(WP):: freq_min,freq_max,center_freq,result
      real(WP),external:: u_f2_asymptoticdisplacement
      real(WP),external:: integrand
      double precision,external:: mydrsple
      real(WP):: lambda,deltal,omega_k
      integer:: cpustart,cpuend,cpurate,cpumax
      real:: cputime
      logical:: do_file_out
      real(WP),external:: qromb77
      real(WP),external:: qsimp_simple
      ! filter in frequency range
      logical,parameter:: ONLY_FREQUENCY_RANGE        = .false.
      real(WP),parameter:: DEGREE_SHIFT              = 0.0
      logical,parameter:: DO_BY_SUMMATION             = .false.

      ! angular degree
      l_min = 0.0
      l_max = ACCURACY

      ! set parameters
      ray_orbit = orbit
      ray_time = time

      ! determine angular degree depending on desired frequency range
      if ( ONLY_FREQUENCY_RANGE ) then
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
      if ( abs(time) <= 1.e-8 .and. orbit == 1) then
        ! determine frequencies
        freq_min = cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS
        freq_max = cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS !huge(freq_max)

        print *,'  asymptotic solution f2 for R1:'
        print *,'    angular degrees min/max:',l_min,l_max
        print *,'    frequency min/max :',sngl(freq_min),sngl(freq_max)
        print *,'    period min/max    :',sngl(2*PI/freq_min),sngl(2*PI/freq_max)
      endif

      if (.not. allocated(X) ) then
        ! for spline representation
        length = l_max-l_min+1
        allocate(X(length),Y(length),Q(3,length),F(3,length),stat=ierror)
        if ( ierror /= 0) call stopProgram('u_f2_orbit - error allocating spline arrays ')

        ! get spline representation
        ! initialize
        ilength = length
        i1 = 1
        i2 = length
        icount = 0
        open(88,file='OUTPUT/integrals_discretel.dat')
        write(88,*) '#format: degree_l integralvalue'
        do_file_out = .true.
        do i = l_min,l_max
          ! starting degrees from 0 to degreeaccuracy, but array indices from 1 to length
          icount = icount + 1
          X(icount)=dble(i)

          ! get integral value with integral boundaries [0,PI]
          call qsimp_degree(integrand,i,0.0,PI,integralvalue,mu)
          Y(icount)=dble(integralvalue)

          if ( do_file_out ) write(88,*) i,Y(icount)
        enddo
        close(88)

        ! calculates spline arrays
        Q(:,:)=0.0d0
        F(:,:)=0.0d0
        call mydrspln(i1,i2,X,Y,Q,F,ilength)

        ! test
        open(88,file='OUTPUT/integrals_nonintegerl.test.dat')
        write(88,*) '#format: degree_l integralvalue'
        do degreel = l_min,l_max,0.2d0
          integralvalue = mydrsple(i1,i2,X,Y,Q,ilength,dble(degreel))
          write(88,*) degreel,integralvalue
        enddo
        close(88)

        ! debug output
        open(88,file='OUTPUT/integrals_omega.dat')
        write(88,*) '#format: degree_l integralvalue'
      else
        do_file_out = .false.
      endif

      ! take time
      !call system_clock(COUNT=cpustart,COUNT_RATE=cpurate,COUNT_MAX=cpumax)

      ! summation or integration
      if ( DO_BY_SUMMATION ) then
        ! summation
        deltal = 0.5d0
        u_f2_orbit = 0.0d0
        do degreel = l_min,l_max,deltal
          ! determines frequency
          lambda = degreel + 0.5d0
          omega_k= cphase * sqrt(lambda**2 - 0.25d0) / EARTHRADIUS

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

        ! qsimp integration
        u_f2_orbit = qsimp_simple(u_f2_asymptoticdisplacement,dble(freq_min),dble(freq_max))
        ! qromb integration
        !u_f2_orbit = qromb77(u_f2_asymptoticdisplacement,freq_min,freq_max)


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
      !if ( abs(cputime) > 0.3 ) then
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
      return
      end

!-----------------------------------------------------------------------
      function u_f2_asymptoticdisplacement(omega)
!-----------------------------------------------------------------------
      use precisions
      use splinefunction
      use ray_orbits
      implicit none
      real(WP):: u_f2_asymptoticdisplacement
      real(WP),intent(in):: omega
      integer:: orbit,iorb,length
      real(WP):: cphase,epiDelta,time,omegaSquare,sineDelta,radiusbycSquared
      real(WP):: sigma2half,integralvalue,lambda,degreel,termSigma, &
                termSpreading,termIntegral,cphaseSquare,twoPI,PIby4
      complex(WP) :: A_lplus,A_lminus,A_even,A_odd,A_l
      double precision,external:: mydrsple
      logical:: ODD_ORBIT
      complex(WP),parameter:: II = ( 0.0, 1.0 )

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
      if ( orbit < 1 ) stop "u_f2_asymptoticdisplacement - wrong orbit"

      ! get corresponding lambda and degree
      omegaSquare = omega**2
      lambda = sqrt(omegaSquare*radiusbycSquared + 0.25)
      degreel = lambda - 0.5

      !get integral value with integral boundaries [0,PI]
      integralvalue = mydrsple(i1,i2,X,Y,Q,ilength,dble(degreel))
      ! coefficients A
      termSigma = exp(-omegaSquare*sigma2half)
      termSpreading = sqrt( 1.0/(twoPI*lambda*sineDelta))
      termIntegral = cphaseSquare*lambda*integralvalue

      !!! until now it is approximated, if we want to use it without approximation l >> 1
      !!! we have to add this term
      !!! termIntegral = termIntegral * sqrt(1.0 - 1.0/(4*lambda*lambda))

      ! orbit adding: epiDelta + 2PI + 2PI + ...
      if ( mod(orbit,2) == 1 ) then
        ! iorbit = 0,1,2,3,4,.. for R1,R3,R5,R7,...
        iorb = (orbit-1)/2
        ODD_ORBIT = .true.
      else
        ! iorbit = 1,2,3,4,.. for R2,R4,R6,...
        iorb = orbit/2
        ODD_ORBIT = .false.
      endif
      if ( ODD_ORBIT ) then
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
      real(WP):: integrand_cont
      real(WP),intent(in):: degree
      real(WP),intent(in):: theta,mu
      double precision,external:: legendrePolynomial_cont
      real(WP):: mu2

      mu2 = mu*mu
      integrand_cont = legendrePolynomial_cont(degree,cos(theta)) &
                  *sin(theta)*exp(-theta*theta/(2.0d0*mu2))/mu2
      return
      end



!-----------------------------------------------------------------------
      function u_f2_R2(epiDelta,time,cphase,sigma,mu, &
                                              ACCURACY)
!-----------------------------------------------------------------------
! calculation of the analytic solution for a source with only a forcing term 2
! (solution by Tape, 2003: eq. 3.34)
! uses an asymptotic formulation for the legendre polynom -> ray theory formulism
! and takes only first arrival solution R1
!
! inputs:
!   epiDelta                 - epicentral distance of receiver/position of
!                               displacement location (in radian)
!   time                       - given time step
!   cphase                  - phase velocity
!   sigma                     - timeparametersigma
!   mu                          - widthparametermu
!
! returns: displacement u_f2
      use precisions
      use filterType, only: bw_waveperiod
      implicit none
      real(WP):: u_f2_R2
      real(WP),intent(in):: epiDelta,time,cphase,sigma,mu
      integer,intent(in):: ACCURACY
      integer:: l, l_min,l_max
      real(WP):: integralvalue
      real(WP):: eigenfrequency,sigma2,reciEarthradius
      real(WP):: freq_min,freq_max,center_freq, A_l
      logical,parameter:: ONLY_FREQUENCY_RANGE = .false.
      real(WP),external:: integrand,integrand_asymptotic

      ! check to have double precision
      sigma2 = sigma*sigma
      reciEarthradius = 1.0/EARTHRADIUS
      l_min = 0
      l_max = ACCURACY
      freq_min = 0.0
      freq_max = huge(freq_max)
      center_freq = 0.0
      ! determine angular degree depending on desired frequency range
      if ( ONLY_FREQUENCY_RANGE ) then
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
      if ( abs(time) <= 1.e-4) then
      !  l = 0
      !  freq_min = cphase*sqrt(l*(l+1.0d0))*reciEarthradius
      !  l = DEGREEACCURACY
      !  freq_max = cphase*sqrt(l*(l+1.0d0))*reciEarthradius
        print *,'asymptotic solution f2 for R2:'
        print *,'    angular degrees min/max:',l_min,l_max
        print *,'    desired frequency range min/max:',sngl(freq_min),sngl(freq_max)
        print *,'         period range min/max:',sngl(2*PI/freq_min),sngl(2*PI/freq_max)

        ! check
        print *,'    actual frequency range min/max:', &
                sngl(cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS), &
                sngl(cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS)
        print *,'         period range min/max:', &
                sngl(2*PI/(cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS)), &
                sngl(2*PI/(cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS))
        if ( l_min > 0 ) then
          print *,'    center frequency/period:',sngl(center_freq), &
                      sngl(2*PI/center_freq)
        endif
      endif


      u_f2_R2 = 0.0d0
      do l = l_min,l_max
        !get eigenfrequency for this degree l
        eigenfrequency= cphase*sqrt(l*(l+1.0d0))*reciEarthradius

        !get integral value with integral boundaries [0,PI]
        call qsimp_degree(integrand,l,0.0,PI,integralvalue,mu)

        ! asymptotic gives NaN
        !call qsimp_degree(integrand_asymptotic,l,0.0d0,PI,integralvalue,mu)

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
      end

!-----------------------------------------------------------------------
      function u_f2_R3(epiDelta,time,cphase,sigma,mu, &
                                              ACCURACY)
!-----------------------------------------------------------------------
! calculation of the analytic solution for a source with only a forcing term 2
! (solution by Tape, 2003: eq. 3.34)
! uses an asymptotic formulation for the legendre polynom -> ray theory formulism
! and takes only first arrival solution R1
!
! inputs:
!   epiDelta                 - epicentral distance of receiver/position of
!                               displacement location (in radian)
!   time                       - given time step
!   cphase                  - phase velocity
!   sigma                     - timeparametersigma
!   mu                          - widthparametermu
!
! returns: displacement u_f2
      use precisions
      use filterType, only: bw_waveperiod
      implicit none
      real(WP):: u_f2_R3
      real(WP),intent(in):: epiDelta,time,cphase,sigma,mu
      integer,intent(in):: ACCURACY
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
      if ( ONLY_FREQUENCY_RANGE ) then
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
      if ( abs(time) <= 1.e-4) then
      !  l = 0
      !  freq_min = cphase*sqrt(l*(l+1.0d0))*reciEarthradius
      !  l = DEGREEACCURACY
      !  freq_max = cphase*sqrt(l*(l+1.0d0))*reciEarthradius
        print *,'asymptotic solution f2 for R3:'
        print *,'    angular degrees min/max:',l_min,l_max
        print *,'    desired frequency range min/max:',sngl(freq_min),sngl(freq_max)
        print *,'         period range min/max:',sngl(2*PI/freq_min),sngl(2*PI/freq_max)

        ! check
        print *,'    actual frequency range min/max:', &
                sngl(cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS), &
                sngl(cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS)
        print *,'         period range min/max:', &
                sngl(2*PI/(cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS)), &
                sngl(2*PI/(cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS))
        if ( l_min > 0 ) then
          print *,'    center frequency/period:',sngl(center_freq), &
                      sngl(2*PI/center_freq)
        endif
      endif


      u_f2_R3 = 0.0d0
      do l = l_min,l_max
        !get eigenfrequency for this degree l
        eigenfrequency= cphase*sqrt(l*(l+1.0d0))*reciEarthradius

        !get integral value with integral boundaries [0,PI]
        call qsimp_degree(integrand,l,0.0,PI,integralvalue,mu)

        ! asymptotic gives NaN
        !call qsimp_degree(integrand_asymptotic,l,0.0d0,PI,integralvalue,mu)

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
      end

!-----------------------------------------------------------------------
      function u_f2_R4(epiDelta,time,cphase,sigma,mu, &
                                              ACCURACY)
!-----------------------------------------------------------------------
! calculation of the analytic solution for a source with only a forcing term 2
! (solution by Tape, 2003: eq. 3.34)
! uses an asymptotic formulation for the legendre polynom -> ray theory formulism
! and takes only first arrival solution R1
!
! inputs:
!   epiDelta                 - epicentral distance of receiver/position of
!                               displacement location (in radian)
!   time                       - given time step
!   cphase                  - phase velocity
!   sigma                     - timeparametersigma
!   mu                          - widthparametermu
!
! returns: displacement u_f2
      use precisions
      use filterType, only: bw_waveperiod
      implicit none
      real(WP):: u_f2_R4
      real(WP),intent(in):: epiDelta,time,cphase,sigma,mu
      integer,intent(in):: ACCURACY
      integer:: l, l_min,l_max
      real(WP):: integralvalue
      real(WP):: eigenfrequency,sigma2,reciEarthradius
      real(WP):: freq_min,freq_max,center_freq, A_l
      logical,parameter:: ONLY_FREQUENCY_RANGE = .false.
      real(WP),external:: integrand,integrand_asymptotic

      ! check to have double precision
      sigma2 = sigma*sigma
      reciEarthradius = 1.0/EARTHRADIUS
      l_min = 0
      l_max = ACCURACY
      freq_min = 0.0
      freq_max = huge(freq_max)
      center_freq = 0.0

      ! determine angular degree depending on desired frequency range
      if ( ONLY_FREQUENCY_RANGE ) then
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
      if ( abs(time) <= 1.e-4) then
      !  l = 0
      !  freq_min = cphase*sqrt(l*(l+1.0d0))*reciEarthradius
      !  l = ACCURACY
      !  freq_max = cphase*sqrt(l*(l+1.0d0))*reciEarthradius
        print *,'asymptotic solution f2 for R4:'
        print *,'    angular degrees min/max:',l_min,l_max
        print *,'    desired frequency range min/max:',sngl(freq_min),sngl(freq_max)
        print *,'         period range min/max:',sngl(2*PI/freq_min),sngl(2*PI/freq_max)

        ! check
        print *,'    actual frequency range min/max:', &
                sngl(cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS), &
                sngl(cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS)
        print *,'         period range min/max:', &
                sngl(2*PI/(cphase*sqrt(l_min*(l_min+1.0d0))/EARTHRADIUS)), &
                sngl(2*PI/(cphase*sqrt(l_max*(l_max+1.0d0))/EARTHRADIUS))
        if ( l_min > 0 ) then
          print *,'    center frequency/period:',sngl(center_freq), &
                      sngl(2*PI/center_freq)
        endif
      endif


      u_f2_R4 = 0.0d0
      do l = l_min,l_max
        !get eigenfrequency for this degree l
        eigenfrequency= cphase*sqrt(l*(l+1.0d0))*reciEarthradius

        !get integral value with integral boundaries [0,PI]
        call qsimp_degree(integrand,l,0.0,PI,integralvalue,mu)

        ! asymptotic gives NaN
        !call qsimp_degree(integrand_asymptotic,l,0.0d0,PI,integralvalue,mu)

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
      end



!-----------------------------------------------------------------------
!     originally from: Sample page from NUMERICAL RECIPES IN Fortran
!     Copyright (C) 1986-1992 by Cambridge University
!     chapter 4
!     http://www.nr.com
!-----------------------------------------------------------------------
      function qsimp_simple(func,a,b)
!-----------------------------------------------------------------------
      ! USES trapzd
      ! returns as s the integral of the function func from a to b.
      ! The parameters EPS can be set to the desired fractional accuracy
      ! and JMAX so that 2tothepowerJMAX-1 is the
      ! maximum allowed number of steps.
      ! Integration is performed by Simpson's rule.
      use precisions
      implicit none
      real(WP):: qsimp_simple
      double precision,intent(in):: a,b
      real(WP),external:: func
      integer,parameter:: JMAX = 40
      double precision:: EPS = 1.e-4
      integer:: j
      double precision:: os,ost,st,s

      ost= 0.0d0  ! originally: -1.e30
      os= 0.0d0   ! originally: -1.e30
      do j = 1,JMAX
        ! trapezoid rule
        call trapzdf77_simple(func,a,b,st,j)
        s=(4.0*st-ost)/3.0
        if (j > 5) then
          !Avoidspurious earlyconvergence.
          if (abs(s-os) < EPS*abs(os).or.(s == 0.0 .and. os == 0.0)) then
            qsimp_simple = s
            return
          endif
        endif
        os = s
        ost = st
      enddo
      print *,'qsimp_simple not converging:',a,b,s,st
      stop "too many steps in qsimp_simple"
      END

!-----------------------------------------------------------------------
      subroutine trapzdf77_simple(func,a,b,s,n)
!-----------------------------------------------------------------------
      use precisions
      implicit none
      integer,intent(in):: n
      double precision,intent(in):: a,b
      double precision,intent(out):: s
      real(WP),external:: func
      !integer degree
      !This routine computes then th stage of refinementof
      !an extended trapezoidal rule.
      ! func is input as the name of the function to be integrated
      ! between limits a and b, also input.
      !When called with n=1, the routine returns as s the
      ! crudest estimate of b a f(x)dx.
      ! Subsequent calls with n=2,3,... (in that sequential order)will
      ! improve the accuracy of s by adding 2n-2 additional interior points.
      ! s should not be modified between sequential calls.
      integer:: it,j
      double precision:: del,sum,tnm,x,f_a,f_b

      if (n == 1) then
        if ( WP == 4 ) then
          s = 0.5*(b-a)*dble(func(sngl(a))+func(sngl(b)))
        else
          s = 0.5*(b-a)*(func(a)+func(b))
        endif
      else
        it=2**(n-2)
        tnm = it
        !Thisisthespacingof thepointstobeadded.
        del=(b-a)/tnm
        x = a+0.5*del
        sum = 0.
        do j = 1,it
          if ( WP == 4 ) then
            sum = sum+dble(func(sngl(x)))
          else
            sum=sum+func(x)
          endif
          x = x+del
        enddo
        !This replaces s byitsrefinedvalue.
        s = 0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END


!-----------------------------------------------------------------------
!     originally from: Sample page from NUMERICAL RECIPES IN Fortran
!     Copyright (C) 1986-1992 by Cambridge University
!     chapter 4
!     http://www.nr.com
!-----------------------------------------------------------------------
      subroutine qsimp_degree(func,degree,a,b,s,mu)
!-----------------------------------------------------------------------
      use precisions
      implicit none
      integer,intent(in):: degree
      real(WP),intent(in):: a,b,mu
      real(WP),intent(out):: s
      real(WP),external:: func
      integer,parameter:: JMAX = 40
      real(WP),parameter:: EPS = 1.e-4
      ! USES trapzd
      ! returns as s the integral of the function func from a to b.
      ! The parameters EPS can be set to the desired fractional accuracy
      ! and JMAX so that 2tothepowerJMAX-1 is the
      ! maximum allowed number of steps.
      ! Integration is performed by Simpson's rule.
      integer:: j
      double precision:: os,ost,st

      ost= 0.0  ! originally: -1.e30
      os= 0.0   ! originally: -1.e30
      do j = 1,JMAX
        ! trapezoid rule
        call trapzdf77(func,degree,dble(a),dble(b),st,j,dble(mu))
        s=(4.0*st-ost)/3.0
        if (j > 5) then
          !Avoidspurious earlyconvergence.
          if (abs(s-os) < EPS*abs(os).or.(s == 0.0 .and. os == 0.0)) return
        endif
        os = s
        ost = st
      enddo
      print *,'qsimp not converging:',degree,a,b,s
      stop "too many steps in qsimp"
      END

!-----------------------------------------------------------------------
      subroutine trapzdf77(func,degree,a,b,s,n,mu)
!-----------------------------------------------------------------------
      use precisions
      implicit none
      integer,intent(in):: n,degree
      double precision,intent(in):: a,b
      double precision,intent(in):: mu
      double precision,intent(out):: s
      real(WP),external:: func
      !integer degree
      !This routine computes then th stage of refinementof
      !an extended trapezoidal rule.
      ! func is input as the name of the function to be integrated
      ! between limits a and b, also input.
      !When called with n=1, the routine returns as s the
      ! crudest estimate of b a f(x)dx.
      ! Subsequent calls with n=2,3,... (in that sequential order)will
      ! improve the accuracy of s by adding 2n-2 additional interior points.
      ! s should not be modified between sequential calls.
      integer:: it,j
      double precision:: del,sum,tnm,x,f_a,f_b

      if (n == 1) then
        if ( WP == 4 ) then
          s = 0.5*(b-a)*dble(func(degree,sngl(a),sngl(mu))+func(degree,sngl(b),sngl(mu)))
        else
          s = 0.5*(b-a)*(func(degree,a,mu)+func(degree,b,mu))
        endif
      else
        it=2**(n-2)
        tnm = it
        !Thisisthespacingof thepointstobeadded.
        del=(b-a)/tnm
        x = a+0.5*del
        sum = 0.
        do j = 1,it
          if ( WP == 4 ) then
            sum = sum+dble(func(degree,sngl(x),sngl(mu)))
          else
            sum=sum+func(degree,x,mu)
          endif
          x = x+del
        enddo
        !This replaces s byitsrefinedvalue.
        s = 0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END


!-----------------------------------------------------------------------
      function qromb77(func,a,b)
!-----------------------------------------------------------------------
      !USE nrtype; USE nrutil, only: nrerror
      !USE nr, only: polint,trapzd
      use precisions
      implicit none
      REAL(WP), INTENT(IN) :: a,b
      REAL(WP) ::qromb77
      !INTERFACE
      !function func(x)
      !USE nrtype
      !REAL(SP), DIMENSION(:), INTENT(IN) :: x
      !REAL(SP), DIMENSION(size(x)) :: func
      !end function func
      !END INTERFACE
      real(WP),external:: func
      INTEGER, PARAMETER :: JMAX = 40,JMAXP = JMAX+1,K = 5,KM = K-1
      REAL(WP), PARAMETER :: EPS=1.0e-4
      real(WP), parameter:: EPS_ZERO = 1.e-4
      ! returns the integral of the function func from a to b.
      ! Integration is performed by Romberg's
      ! method of order 2K, where, e.g., K=2 is Simpson's rule.
      ! Parameters: EPS is the fractional accuracy desired, as determined by the extrapolation er-
      ! ror estimate; JMAX limits the total number of steps; K is the number of points used in the
      ! extrapolation.
      REAL(WP), DIMENSION(JMAXP) :: h,s ! These store the successive trapezoidal ap-
                                           ! proximations and their relative step sizes.
      REAL(WP) ::dqromb77
      INTEGER:: j

      h(1)=1.0
      do j = 1,JMAX
        call trapzdf77_simple(func,a,b,s(j),j)
        if (j >= K) then
          call polint77(h(j-KM:j),s(j-KM:j),0.0,qromb77,dqromb77)

          ! org:
          if (abs(dqromb77) <= EPS*abs(qromb77)) return

          ! from qsimp:
          if (j > 5) then
            !Avoid spurious early convergence.
            if (abs(dqromb77) <= EPS*abs(qromb77) .or. (abs(qromb77) < EPS_ZERO .and. abs(dqromb77) < EPS_ZERO)) return
          endif
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)    ! This is a key step: The factor is 0.25 even
                            ! though the step size is decreased by only
                            ! 0.5. This makes the extrapolation a poly-
                            ! nomial in h2 as allowed by equation(4.2.1),
                            ! not just a polynomial in h.
      enddo
      !call nrerror('qromb: too many steps')
      print *,'qromb77:',a,b,qromb77,dqromb77
      stop "qromb77 too many steps"
      end function qromb77



!-----------------------------------------------------------------------
      subroutine polint77(xa,ya,n,x,y,dy)
!-----------------------------------------------------------------------
      use precisions
      implicit none
      integer,parameter:: nmax = 10
      integer,intent(in):: n
      real(WP),intent(in):: xa(n),ya(n),x
      real(WP),intent(out):: y,dy
      real(WP):: c(nmax),d(nmax),dif,dift,ho,hp,w,den
      integer:: i,ns,m

      ns = 1
      dif = abs(x-xa(1))
      do i = 1,n
        dift = abs(x-xa(i))
        if (dift < dif) then
          ns = i
          dif = dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
      enddo

      y=ya(ns)
      ns = ns-1
      do m = 1,n-1
        do i = 1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w = c(i+1)-d(i)
          den = ho-hp
          if (den == 0.) pause "Problem in Polint"
          den = w/den
          d(i)=hp*den
          c(i)=ho*den
        enddo
        if (2*ns < n-m) then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns = ns-1
        endif
        y = y+dy
      enddo
      return
      end
