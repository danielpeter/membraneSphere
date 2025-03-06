!=====================================================================
!
!       m e m b r a n e S p h e r e  1 . 0
!       --------------------------------------------------
!
!      Daniel Peter
!      ETH Zurich - Institute of Geophysics
!      (c) ETH July 2006
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

!-----------------------------------------------------------------------
      double precision function u_shape(colatitude,longitude,time,phasevelocity)
!-----------------------------------------------------------------------
! calculation of the analytic solution for a source with only 
! a given initial shape displacement (no time term)
!
! inputs:
!   colatitude, longitude - position of displacement location
!   time                       - given time step
!
! returns: displacement u_shape
      use verbosity; use precision
      implicit none
      double precision:: colatitude,longitude,time,phasevelocity
      double precision:: integrand, eigenfrequency
      double precision:: imu, legendrePolynomial
      double precision:: cphase,fromA,toB
      integer:: l
      external:: integrand
      
      cphase = phasevelocity
       
      !print*,"u_shape: ", colatitude, longitude,time
      !integral boundaries
      fromA = 0.0
      toB = PI
       
      u_shape=0.0
      do l=0,DEGREEACCURACY
        !get eigenfrequency for this degree l
        eigenfrequency= cphase*sqrt(l*(l+1.0))/EARTHRADIUS
        !print*,"eigenfrequency = ",eigenfrequency
        
        !get integral value
        call qsimp_degree(integrand, l, fromA, toB, imu)
        
        !debug
        if(DEBUG) print*,"integral value = ",imu
                
        !calculate displacement
        u_shape=u_shape+(l+0.5)*imu*cos(eigenfrequency*time)*legendrePolynomial(l,dcos(colatitude))        
      enddo
    
      return
      end
      
!-----------------------------------------------------------------------
      double precision function u_f2(epiDelta,time,phasevelocity)
!-----------------------------------------------------------------------
! calculation of the analytic solution for a source with only a forcing term 2
!
! inputs:
!   epiDelta                  - epicentral distance of receiver/position of displacement location (in radian)
!   time                       - given time step
!
! returns: displacement u_f2
      use verbosity; use precision
      implicit none
      real(WP):: epiDelta,time,phasevelocity
      integer:: l
      double precision:: integrand,integralvalue,cphase,legendrePolynomial
      double precision:: depiDelta,dtime,eigenfrequency,dpi,dsigma2,dreciEarthradius
      external:: integrand
            
      ! check to have double precision
      cphase = dble(phasevelocity)
      depiDelta = dble(epiDelta)
      dtime = dble(time) 
      dpi = dble(PI)
      dsigma2 = dble(TimeParameterSigma)*dble(TimeParameterSigma) 
      dreciEarthradius = 1.0d0/dble( EARTHRADIUS )
      
      !print*,"u_shape: ", colatitude, longitude,time
       
      u_f2=0.0d0

      do l=0,DEGREEACCURACY
        !get eigenfrequency for this degree l
        eigenfrequency= cphase*dsqrt(l*(l+1.0d0))*dreciEarthradius
        
        !get integral value with integral boundaries [0,PI]
        call qsimp_degree(integrand, l, 0.0d0, dpi, integralvalue)                
        
        !calculate displacement
        u_f2=u_f2+(l+0.5d0)*integralvalue*dcos(eigenfrequency*dtime)   &
              *dexp(-dsigma2*(eigenfrequency**2)*0.5d0)*legendrePolynomial(l,dcos(depiDelta))        
      enddo
      
      u_f2=cphase*cphase*u_f2

      !debug
      if(DEBUG) print*,'u_f2 =',u_f2,integralvalue,cphase,depiDelta,dtime

      return
      end
      
!-----------------------------------------------------------------------      
      double precision function integrand(degree,theta)
!-----------------------------------------------------------------------
! integrand function defined by Tape, chapter 3, eq. (3.26), p. 23 
!
! input:
!     degree  - degree of legendre polynomial to take
!     theta     - value to calculate at
!
! returns: integrand value
      use precision
      implicit none
      double precision:: legendrePolynomial, theta  
      double precision:: mu2
      integer:: degree
      
      !print*,"integrand function for ", degree, theta
      
      mu2 = WidthParameterMu**2
      integrand=legendrePolynomial(degree,dcos(theta))*dsin(theta)*dexp(-theta**2/(2.0d0*mu2))/mu2      
      return
      end

!-----------------------------------------------------------------------      
!     originally from: Sample page from NUMERICAL RECIPES IN FORTRAN
!     Copyright (C) 1986-1992 by Cambridge University 
!     chapter 4
!     http://www.nr.com 
!-----------------------------------------------------------------------      
      SUBROUTINE qsimp_degree(func,degree,a,b,s) 
!-----------------------------------------------------------------------
      INTEGER JMAX 
      double precision:: EPS 
      double precision:: func, a, b, s
      EXTERNAL:: func 
      PARAMETER (EPS=1.e-4, JMAX=40)
      integer:: degree 
      ! USES trapzd 
      ! Returns as s the integral of the function func from a to b. 
      ! The parameters EPS can be set to the desired fractional accuracy 
      ! and JMAX so that 2tothepowerJMAX-1is the 
      ! maximum allowed number of steps. 
      ! Integration is performed by Simpson’s rule. 
      INTEGER j 
      double precision os,ost,st 
            
      ost= 0.0d0  ! originally: -1.e30
      os= 0.0d0   ! originally: -1.e30
      do j=1,JMAX 
        !print*,"integration step ",j        
        ! trapezoid rule
        call trapzdf77(func,degree,a,b,st,j) 
        s=(4.0d0*st-ost)/3.0d0 
        if (j.gt.5) then 
          !Avoidspurious earlyconvergence. 
          if(abs(s-os).lt.EPS*abs(os).or.(s.eq.0.0 .and. os.eq.0.0)) return 
        endif 
        os=s 
        ost=st 
      enddo  
      print*,'not converging:',degree,a,b,s
      call stopProgram( "too many steps in qsimp" )
      END 
      
!-----------------------------------------------------------------------
      SUBROUTINE trapzdf77(func,degree,a,b,s,n) 
!-----------------------------------------------------------------------
      integer n,degree 
      double precision a,b,s,func 
      EXTERNAL func 
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
      integer it,j 
      double precision del,sum,tnm,x, f_a, f_b 
      
      !common/exactsolution/ degree
      !integer degree
      
      !degree = deg
      
      !print*,"trapzd with", degree, a, b, s, n
      
      if (n.eq.1) then 
        s=0.5*(b-a)*(func(degree,a)+func(degree,b)) 
      else 
        it=2**(n-2) 
        tnm=it 
        !Thisisthespacingof thepointstobeadded. 
        del=(b-a)/tnm 
        
        x=a+0.5*del 
        sum=0. 
        do j=1,it 
          sum=sum+func(degree,x) 
          x=x+del 
        enddo 
        !This replaces s byitsrefinedvalue. 
        s=0.5*(s+(b-a)*sum/tnm) 
        
      endif 
      
      !print*,"trapzd ", s
      return 
      END
      
