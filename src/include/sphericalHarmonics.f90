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
      use precisions
      implicit none
      integer,intent(in)::vertex
      integer:: l, m
      double precision:: colatitude,longitude,harmonical,generalizedLegendre
      real(WP):: colat,lon
      external:: generalizedLegendre

      !translate Cartesian to spherical coordinates for a
      ! given triangle corner on the sphere
      call getSphericalCoord(vertex, colat, lon)
      colatitude = colat
      longitude = lon

      !calculate the real part of spherical harmonic function
      l= DEGREE_L
      m= DEGREE_M
      harmonical = generalizedLegendre(l,m,colatitude)
      if ( m /= 0) then
        harmonical = harmonical*dsin( m*longitude)
      endif

      !debug
      !print *,l,m,'-harmonical:',harmonical

      ! return value
      sphericalHarmonics = harmonical
      return
      end


!----------------------------------------------------------------------
      double precision function generalizedLegendre(l,m,x)
!----------------------------------------------------------------------
! Computes the generalized Legendre function X_m l (x)
! Here m and l are integers satisfying -l <= m <= l
! (see Tape, thesis, 2003, B4)
! and: http://mathworld.wolfram.com/LegendrePolynomial.html
!
! input:
!    l,m - spherical degrees
!    x    - value to calculate the legendre function with
!
! returns: generalizedLegendre function value
      implicit none
      integer,intent(in):: l,m
      double precision,intent(in):: x
      double precision:: fact1,fact2,associatedLegendre
      integer:: factorial,m_positive
      double precision,parameter:: PI = 3.1415926535897931d0
      external:: associatedLegendre

      ! range of m
      if ( m_positive < 0) then
        m_positive = abs(m)
      else
        m_positive = m
      endif

      !multiplication factors
      fact1 = dsqrt( (2*l+1.0)/(4*PI) )
      fact2 = 1.0*factorial(l-m_positive)/(1.0*factorial(l+m_positive))
      fact2= dsqrt( fact2 )

      !calculates the generalized legendre function value
      generalizedLegendre = fact1*fact2*associatedLegendre(l,m_positive,dcos(x))
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
      integer,intent(in):: l
      double precision,intent(in):: x
      double precision:: associatedLegendre
      external:: associatedLegendre

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
! http://gershwin.ens.fr/vdaniel/Doc-Locale/Langages-Program-Scientific/Numerical_Recipies/bookf90pdf/chap6f9.pdf
! or
! http://www.physics.louisville.edu/help/nr/bookf90pdf/chap6f9.pdf
!
! and theory: http://mathworld.wolfram.com/LegendrePolynomial.html
!
! input:
!    l,m - spherical degrees
!    x    - value to calculate the legendre function with
!
! returns: associatedLegendre function value
      implicit none
      integer,intent(in):: l,m
      double precision,intent(in):: x
      integer:: i,ll
      double precision:: pLegendre,fact,pll,pmm,pmmp1,somx2

      !print *,'associatedlegendre:',l,m,x

      ! check
      if (m < 0 .or. m > l .or. abs(x) > 1.0) then
        print *,'Error: associatedLegendre error:'
        print *,'       l / m :',l,m
        print *,'       input x:',x
        call stopProgram('abort -associatedLegendre   ')
      endif

      !Compute P_m m
      pmm = 1.0d0
      if (m > 0) then
        somx2 = dsqrt((1.0d0-x)*(1.0d0+x))
        fact = 1.0d0
        do i = 1,m
          pmm=-pmm*fact*somx2
          fact = fact+2.0d0
        enddo
      endif

      if (m == l) then
        pLegendre = pmm
      else
        pmmp1=x*(m+m+1.0d0)*pmm
        !Compute P_m m+1 .
        if (m+1 == l) then
          pLegendre = pmmp1
        else
          !Compute P_m l , l >m+1.
          do ll = m+2,l
            pll=(x*(ll+ll-1.0d0)*pmmp1-(ll+m-1.0d0)*pmm)/(ll-m)
            pmm = pmmp1
            pmmp1 = pll
            ! debug
            !if (abs(pmmp1-pmm) < 1.e-8) exit
          enddo
        pLegendre = pll
        endif
      endif

      associatedLegendre = pLegendre
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
      implicit none
      integer,intent(in):: N
      integer:: i

      !check
      if ( N < 0) stop 'abort - factorial'

      factorial = 1
      do i = 2,N
        factorial = factorial*i
      enddo
      return
      end


!----------------------------------------------------------------------
      double precision function legendrePolynomial_cont(l,x)
!----------------------------------------------------------------------
! Computes the Legendre polynomial P_l (x)
! Here l can be a continuous value, while x lies in the
! range  -1 <= x <= 1.
!
! input:
!    l     - degree
!    x    - value to calculate the legendre function with
!
! returns: legendre polynomial function value
      implicit none
      double precision,intent(in):: l
      double precision,intent(in):: x
      double precision:: associatedLegendre_cont
      external:: associatedLegendre_cont

      legendrePolynomial_cont = associatedLegendre_cont(l,0.0d0,x)
      return
      end


!----------------------------------------------------------------------
      double precision function associatedLegendre_cont(l,m,x)
!----------------------------------------------------------------------
! Computes the associated Legendre polynomial P_m l (x)
! Here m and l are double precisions, satisfying 0 <= m <= l, while x lies in the
! range  -1 <= x <= 1.
!
! routine from: numerical recipes
! http://gershwin.ens.fr/vdaniel/Doc-Locale/Langages-Program-Scientific/Numerical_Recipies/bookf90pdf/chap6f9.pdf
! or
! http://www.physics.louisville.edu/help/nr/bookf90pdf/chap6f9.pdf
!
! and theory: http://mathworld.wolfram.com/LegendrePolynomial.html
!
! input:
!    l,m - spherical degrees
!    x    - value to calculate the legendre function with
!
! returns: associatedLegendre function value
      implicit none
      double precision,intent(in):: l,m
      double precision,intent(in):: x
      double precision:: ll
      integer:: i
      double precision:: pLegendre,fact,pll,pmm,pmmp1,somx2

      !print *,'associatedlegendre:',l,m,x

      ! check
      if (m < 0.0d0 .or. m > l .or. abs(x) > 1.0) &
            call stopProgram('abort -associatedLegendre_cont   ')

      !Compute P_m m
      pmm = 1.0d0
      if (m > 0.0d0) then
        somx2 = dsqrt((1.0d0-x)*(1.0d0+x))
        fact = 1.0d0
        do i = 1,int(m)
          pmm=-pmm*fact*somx2
          fact = fact+2.0d0
        enddo
      endif

      if (abs(m-l) < 1.e-5) then
        pLegendre = pmm
      else
        pmmp1=x*(m+m+1.0d0)*pmm
        !Compute P_m m+1 .
        if ( abs(m+1.0d0-l) < 1.e-5 ) then
          pLegendre = pmmp1
        else
          !Compute P_m l , l >m+1.

          ! example: m == 2, l == 5
          ! old Fortran:
          ! -> loop from ll==4.0, 4.2, .., 5.0
          !do ll = m+2.0d0,l,0.2d0
          ! new Fortran: uses integer loop variables, otherwise gives a warning
          ! -> loop from i = (m+2)/0.2==20, 21, .., (l/0.2)==25
          !              and ll = 20*0.2==4, 21*0.2==4.2,.., 25*0.2==5.0
          do i = int((m+2.d0)/0.2d0), int(l/0.2d0)
            ll = i * 0.2d0
            pll=(x*(ll+ll-1.0d0)*pmmp1-(ll+m-1.0d0)*pmm)/(ll-m)
            pmm = pmmp1
            pmmp1 = pll
            ! debug
            !if (abs(pmmp1-pmm) < 1.e-8) exit
          enddo
        pLegendre = pll
        endif
      endif

      associatedLegendre_cont = pLegendre
      return
      end
