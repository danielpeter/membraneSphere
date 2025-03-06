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

!----------------------------------------------------------------------
      double precision function sphericalHarmonics(vertex)
!-----------------------------------------------------------------------
! calculates the real spherical harmonic value for a given vertex
! with fixed degrees l and m from header file
! (formula used: Tape, 2003, App. B, B.14)
!
! input:
!    vertex  - index to voronoi cell center
!
! returns: sphericalHarmonics value
      use verbosity; use precision
      implicit none
      integer::vertex, l, m
      double precision:: colatitude,longitude,harmonical,generalizedLegendre
      
      !debug
      if(DEBUG) print*,'vertex',vertex
      
      !translate cartesian to spherical coordinates for a
      ! given triangle corner on the sphere
      call getSphericalCoord(vertex, colatitude, longitude)    
      
      !calculate the real part of spherical harmonic function
      l= DEGREE_L
      m= DEGREE_M
      harmonical = generalizedLegendre(l,m,colatitude)      
      if( m .ne. 0) then
        harmonical = harmonical*dsin( m*longitude)
      endif
      !debug
      if(DEBUG) print*,l,m,'-harmonical:',harmonical
      
      ! return value
      sphericalHarmonics = harmonical      
      return
      end
      
            
!----------------------------------------------------------------------
      double precision function generalizedLegendre(l,m,x) 
!----------------------------------------------------------------------
! Computes the generalized Legendre function X_m l (x)
! Here m and l are integers satisfying -l <= m <= l
!
! and: http://mathworld.wolfram.com/LegendrePolynomial.html
!
! input:
!    l,m - spherical degrees
!    x    - value to calculate the legendre function with
!
! returns: generalizedLegendre function value
      use verbosity; use precision
      implicit none
      integer:: l,m
      double precision:: x,fact1,fact2,associatedLegendre
      integer:: factorial
      
      ! range of m
      if( m .lt. 0) then
        m = abs(m)
      endif

      !multiplication factors
      fact1=dsqrt((2.0*l+1.0)/(4.0*dble(PI)))
      fact2=1.0*factorial(l-m)/(1.0*factorial(l+m))
      fact2= dsqrt( fact2)      
      !debug
      if(DEBUG) print*,'factors ',fact1,fact2
      
      !calculates the generalized legendre function value
      generalizedLegendre = fact1*fact2*associatedLegendre(l,m,dcos(x))
      
      !debug
      if(DEBUG)print*,x,'generalizedLegendre',generalizedLegendre
      return
      end
      
      
!----------------------------------------------------------------------
      double precision function legendrePolynomial(l,x) 
!----------------------------------------------------------------------
! Computes the Legendre polynomial P_l (x)
! Here l is an integer, while x lies in the
! range  -1 <= x <= 1. 
!     
! input:
!    l     - degree
!    x    - value to calculate the legendre function with
!
! returns: legendre polynomial function value
      implicit none
      integer:: l
      double precision:: x, associatedLegendre
      
      legendrePolynomial = associatedLegendre(l,0,x)
      return
      end

      
!----------------------------------------------------------------------
      double precision function associatedLegendre(l,m,x) 
!----------------------------------------------------------------------
! Computes the associated Legendre polynomial P_m l (x)
! Here m and l are integers satisfying 0 <= m <= l, while x lies in the
! range  -1 <= x <= 1. 
!
! routine from: numerical recipes
!        http://www.library.cornell.edu/nr/cbookfpdf.html
! and: http://mathworld.wolfram.com/LegendrePolynomial.html
!     
! input:
!    l,m - spherical degrees
!    x    - value to calculate the legendre function with
!
! returns: associatedLegendre function value
      use verbosity
      implicit none
      integer:: l,m,i,ll 
      double precision:: pLegendre,x,fact,pll,pmm,pmmp1,somx2 
      
      if(m.lt.0 .or. m.gt.l .or. abs(x).gt.1.0) call stopProgram('abort -associatedLegendre   ')

      !Compute P_m m 
      pmm=1.0d0 
      if(m .gt. 0) then 
        somx2=dsqrt((1.0d0-x)*(1.0d0+x))
        fact=1.0d0 
        do i=1,m 
          pmm=-pmm*fact*somx2 
          fact=fact+2.0d0 
        enddo
      endif 
      
      if(m .eq. l) then 
        pLegendre=pmm 
      else 
        pmmp1=x*(m+m+1.0d0)*pmm 
        !Compute P_m m+1 . 
        if(m+1 .eq. l) then 
          pLegendre=pmmp1 
        else 
          !Compute P_m l , l >m+1. 
          do ll=m+2,l
            pll=(x*(ll+ll-1.0d0)*pmmp1-(ll+m-1.0d0)*pmm)/(ll-m) 
            pmm=pmmp1
            pmmp1=pll 
          enddo 
        pLegendre=pll 
        endif 
      endif 
      
      associatedLegendre = pLegendre
      !debug
      !if(DEBUG)print*,x,'associatedLegendre', associatedLegendre
      return 
      end
      
      
      
!----------------------------------------------------------------------
      integer function factorial(N) 
!----------------------------------------------------------------------
! calculates the factorial product 1*2*...*N
!
! input:
!    N - integer 
!
! return: factorial Product
      use verbosity
      implicit none
      integer N,i
      
      if( N .lt. 0) stop 'abort - factorial'
           
      factorial = 1
      do i=2,N
        factorial = factorial*i
      enddo
      !debug
      if(DEBUG)print*,'factorial', factorial
      return      
      end      
      
      
