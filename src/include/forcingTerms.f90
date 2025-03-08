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
!   calculations for the forcing terms defined by
!   Carl Tape (2003, chap. 3)
!

!-----------------------------------------------------------------------
      double precision function forceTerm1(colatitude,longitude,time)
!-----------------------------------------------------------------------
! calculation of the forcing term one
!
! inputs:
!   colatitude, longitude - position of displacement location
!   time                       - given time step
!
! returns: forceTerm1
      implicit none
      double precision:: colatitude,longitude,time
      double precision:: timeTerm1, widthTerm

      forceTerm1 = timeTerm1(time)*widthTerm(colatitude)
      return
      end function


!-----------------------------------------------------------------------
      double precision function forceTerm2(colatitude,longitude,time)
!-----------------------------------------------------------------------
! calculation of the forcing term two, assumes source on north pole
!
! inputs:
!   colatitude, longitude - position of displacement location
!   time                       - given time step
!
! returns: forceTerm2
      use precisions
      use verbosity
      implicit none
      double precision,intent(in):: colatitude,longitude,time
      double precision:: widthTerm
      real(WP) :: time_wp
      real(WP),external:: timeTerm2

      ! approximative this term is too small to have an influence beyond this time
      ! if ( time > 7000.0d0) then
      !   forceTerm2 = 0.0d0
      !   return
      ! endif
      time_wp = real(time,kind=WP)

      forceTerm2 = timeTerm2(time_wp) * widthTerm(colatitude)
      return
      end function


!-----------------------------------------------------------------------
      double precision function initialShapeTerm(colatitude,longitude,time)
!-----------------------------------------------------------------------
! calculation of the forcing term two
!
! inputs:
!   colatitude, longitude - position of displacement location
!   time                       - given time step
!
! returns: forceTerm2
      implicit none
      double precision:: colatitude,longitude,time
      double precision:: widthTerm

      if ( abs(time) < 1.e-7) then
        initialShapeTerm= widthTerm(colatitude)
      else
        initialShapeTerm = 0.0
      endif
      return
      end function


!-----------------------------------------------------------------------
      double precision function timeTerm1(time)
!-----------------------------------------------------------------------
! calculation of the time term one
!
! inputs:
!   time                       - given time step
!
! returns: timeTerm1
      use propagationStartup;use verbosity
      implicit none
      double precision:: colatitude,longitude,time,sigma2
      double precision,parameter:: sqrt2pi = 2.506628274631d0


      sigma2 = timeParameterSigma*timeParameterSigma
      timeTerm1 = dexp(-time*time/(sigma2+sigma2))/(sqrt2pi*timeParameterSigma)
      return
      end function

!-----------------------------------------------------------------------
      function timeTerm2(time)
!-----------------------------------------------------------------------
! calculation of the time term two
!
! inputs:
!   time                       - given time step
!
! returns: timeTerm2
      use precisions
      use propagationStartup, only: timeParameterSigma
      implicit none
      real(WP),intent(in):: time
      real(WP):: timeTerm2,sigma2
      real(WP),parameter:: sqrt2pi = 2.506628274631_WP

      sigma2 = timeParameterSigma*timeParameterSigma
      timeTerm2 = (-time/sigma2)*exp(-time*time/(sigma2+sigma2))/(sqrt2pi*timeParameterSigma)
      return
      end function

!-----------------------------------------------------------------------
      double precision function widthTerm(theta)
!-----------------------------------------------------------------------
! calculation of the width term
!
! inputs:
!   theta               - epicentral distance of position of displacement location to source
!                            ( given in radian)
! returns: widthTerm
      use propagationStartup, only: muSquare,muTwo
      implicit none
      double precision,intent(in):: theta
      double precision:: col2

      col2 = theta*theta
      widthTerm = dexp(-col2/(muTwo))/muSquare
      return
      end function

!-----------------------------------------------------------------------
      function forceTerm2Source(refVertex,time,sourceVertex)
!-----------------------------------------------------------------------
! calculation of the forcing term two
!
! inputs:
!   refVertex         - position of displacement location
!   time                       - given time step
!   sourceVertex  - position of source
!
! returns: forceTerm2Source
      use precisions
      implicit none
      integer,intent(in):: sourceVertex, refVertex
      real(WP):: forceTerm2Source
      real(WP):: time
      real(WP),external:: timeTerm2,widthTermSource

      ! approximative: this term is too small to have an influence beyond this time
      ! so for speeding up calculations return here
      ! if ( time > 7000.0d0) then
      !   forceTerm2Source = 0.0d0
      !   return
      ! endif

      forceTerm2Source = timeTerm2(time)*widthTermSource(refVertex,sourceVertex)
      return
      end function

!-----------------------------------------------------------------------
      function widthTermSource(refVertex,sourceLocationVertex)
!-----------------------------------------------------------------------
! calculation of the width term with respect to a source located at a source vertex
!
! inputs:
!   refVertex                     - position of displacement location
!   sourceLocationVertex               -  vertex index of source
!
! returns: widthTermSource
      use propagationStartup; use cells
      implicit none
      integer,intent(in):: sourceLocationVertex, refVertex
      real(WP):: widthTermSource,colatitude,longitude
      real(WP):: theta,col2
      real(WP):: vectorS(3),vectorRef(3)

      ! parameter
      !mu2=widthParameterMu*widthParameterMu

      ! distance to source in radian
      vectorS(:) = vertices(sourceLocationVertex,:)
      vectorRef(:) = vertices(refVertex,:)
      call greatCircleDistance(vectorS,vectorRef,theta)

      col2 = theta*theta
      widthTermSource = exp(-col2/(muTwo))/muSquare

      ! single cell source
      !if ( refVertex == sourceLocationVertex ) then
      !  widthTermSource=1.0d0
      !else
      !  widthTermSource=0.0d0
      !endif

      return
      end function

!-----------------------------------------------------------------------
      function forceAdjointSource(refVertex,step)
!-----------------------------------------------------------------------
! calculation of the forcing term for the adjoint source
!
! inputs:
!   refVertex         - position of displacement location
!   step                 - given time index
!
! returns: forceAdjointSource
      use propagationStartup; use adjointVariables; use cells; use verbosity
      implicit none
      integer,intent(in):: refVertex,step
      integer:: i
      real(WP):: forceAdjointSource
      real(WP):: distance,factor,dx
      logical:: isNeighbor
      logical,parameter:: extendedSource = .false.
      logical,parameter:: Gaussian       = .false.
      logical,parameter:: plateau        = .false.
      logical,parameter:: onlyNeighbors  = .false.

      ! calculate adjoint source term: based upon a single cell displacement
      if ( refVertex == adjointSourceVertex) then
        forceAdjointSource = adjointSource(2,step)
      else
        forceAdjointSource = 0.0_WP
      endif

      ! extend adjoint source to multiple cells
      if ( extendedSource ) then
        ! calculate adjoint source term: based upon a Gaussian cell displacement
        if ( Gaussian .or. plateau ) then
          call greatCircleDistance(vertices(refVertex,:),vertices(adjointSourceVertex,:),distance)
          distance = distance*EARTHRADIUS
          dx = averageCellDistance

          if ( Gaussian ) then
            ! factor determined empirically such that it's a nice curve to distances till 280 km (4*dx for grid level 6)
            factor= exp(-distance*distance/(300*dx))
          else
            ! plateau means we take the same value as at the adjointSource vertex location for all the others as well
            factor = 1.0_WP
          endif

          ! extend to a multiple of the cell distances
          if ( distance <= 4.0*dx ) then
            forceAdjointSource=factor*adjointSource(2,step)
          endif
        endif

        ! enlarge the source to its closest neighbor cells
        if ( onlyNeighbors ) then
          isNeighbor = .false.
          do i = 1,cellNeighbors(adjointSourceVertex,0)
            if ( refVertex == cellNeighbors(adjointSourceVertex,i) ) then
              isNeighbor = .true.
              exit
            endif
          enddo
          if ( isNeighbor ) then
            forceAdjointSource=adjointSource(2,step)
          endif
        endif
      endif

      return
      end function

!-----------------------------------------------------------------------
      function forceTermExact(refVertex,time,sourceLat,sourceLon)
!-----------------------------------------------------------------------
! calculation of the forcing term two
!
! inputs:
!   refVertex         - position of displacement location
!   time                - given time step
!   sourcelat/lon  - position of source (in degree)
!
! returns: forceTerm2Source
      use verbosity; use cells
      use propagationStartup, only: muSquare,muTwo
      implicit none
      integer:: sourceVertex, refVertex
      real(WP):: forceTermExact,time,width,sourceLat,sourceLon,theta
      real(WP):: vectorS(3),vectorRef(3),col2
      real(WP),external:: timeTerm2

      ! distance to source in radian
      call getVector(sourceLat,sourceLon,vectorS(1),vectorS(2),vectorS(3))
      vectorRef(:) = vertices(refVertex,:)
      call greatCircleDistance(vectorS,vectorRef,theta)
      if (theta > PI .or. theta < 0.0) then
        print *,"strange forcetermExact:",theta,time,forceTermExact,refVertex,width
      endif

      col2 = theta*theta
      width = exp(-col2/muTwo)/muSquare

      ! force
      forceTermExact = timeTerm2(time)*width
      return
      end function

