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
      double precision colatitude,longitude,time
      double precision timeTerm1, widthTerm
      
      forceTerm1= timeTerm1(time)*widthTerm(colatitude)
      return
      end

      
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
      use verbosity
      implicit none
      double precision colatitude,longitude,time
      double precision timeTerm2, widthTerm
      
      ! approximative this term is too small to have an influence beyond this time
      ! if( time .gt. 7000.0d0) then
      !   forceTerm2 = 0.0d0
      !   return
      ! endif
      
      forceTerm2= timeTerm2(time)*widthTerm(colatitude)
      
      !debug
      if(DEBUG)print*,"force:",colatitude,longitude,time,forceTerm2
      
      return
      end      
      
      
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
      double precision colatitude,longitude,time
      double precision timeTerm2, widthTerm
      
      if( abs( time - 0.0) .lt. 0.0000001) then
        initialShapeTerm= widthTerm(colatitude)
      else
        initialShapeTerm=0.0
      endif
      return
      end      
      
      
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
      double precision colatitude,longitude,time,sigma2
      double precision sqrt2pi/2.506628274631d0/
      
      
      sigma2=timeParameterSigma*timeParameterSigma
      timeTerm1 = dexp(-time*time/(sigma2+sigma2))/(sqrt2pi*timeParameterSigma)
      
      !debug
      if(DEBUG) print*,"time term1:",time,timeTerm1
      return
      end

!-----------------------------------------------------------------------
      function timeTerm2(time)
!-----------------------------------------------------------------------
! calculation of the time term two
!
! inputs:
!   time                       - given time step
!
! returns: timeTerm2
      use propagationStartup; use verbosity
      implicit none
      real(WP):: timeTerm2,time,sigma2
      real(WP),parameter:: sqrt2pi = 2.506628274631_WP
      
      sigma2 = timeParameterSigma*timeParameterSigma
      timeTerm2 =(-time/sigma2)*exp(-time*time/(sigma2+sigma2))/(sqrt2pi*timeParameterSigma)
      !debug
      !if(DEBUG) then
      !  if( time .eq. 100.0) print*,"time term2:", time, timeTerm2
      !endif
      return
      end      
        
!-----------------------------------------------------------------------
      double precision function widthTerm(theta)
!-----------------------------------------------------------------------
! calculation of the width term
!
! inputs:
!   theta               - epicentral distance of position of displacement location to source
!                            ( given in radian)
! returns: widthTerm
      use propagationStartup; use verbosity
      implicit none
      double precision:: theta, mu2, col2
      
      mu2=widthParameterMu*widthParameterMu
      col2 = theta*theta
      widthTerm = dexp(-col2/(mu2+mu2))/mu2

      !debug
      if(DEBUG) print*,"width term:",theta, widthTerm
      
      return
      end  
      
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
      use precision; use verbosity
      implicit none
      integer:: sourceVertex, refVertex
      real(WP):: forceTerm2Source,time,timeTerm2,widthTermSource
      external::timeTerm2,widthTermSource
      
      ! approximative: this term is too small to have an influence beyond this time
      ! so for speeding up calculations return here
      ! if( time .gt. 7000.0d0) then
      !   forceTerm2Source = 0.0d0
      !   return
      ! endif
      
      forceTerm2Source = timeTerm2(time)*widthTermSource(refVertex,sourceVertex)
           
      !debug
      if(DEBUG) then
        if(refVertex .eq. sourceVertex .and. forceTerm2Source .ne. 0.0) print*,"forceterm2source:",time,forceTerm2Source,refVertex
        ! check if Nan
        if( forceTerm2Source .ne. forceTerm2Source)then
          print*,'forceTerm:',forceTerm2Source
          print*,time,refVertex,sourceVertex,timeTerm2(time),widthTermSource(refVertex,sourceVertex)
          call stopProgram( 'forceTerm2Source is NaN      ')
        endif
      endif
      
      return
      end     
       
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
      use propagationStartup; use cells; use verbosity
      implicit none
      real(WP):: widthTermSource,colatitude, longitude
      integer:: sourceLocationVertex, refVertex
      real(WP):: theta,mu2, col2
      real(WP):: vectorS(3),vectorRef(3)
      
      ! parameter
      mu2=widthParameterMu*widthParameterMu
      
      ! distance to source in radian
      vectorS(:) = vertices(sourceLocationVertex,:)
      vectorRef(:) = vertices(refVertex,:)
      call greatCircleDistance(vectorS,vectorRef,theta)
      
      !debug
      if(DEBUG) then
        if( theta .gt. PI) print*,'widthTermSource has impossible theta',theta
      endif
      
      col2 = theta*theta
      widthTermSource = exp(-col2/(mu2+mu2))/mu2
      
      !debug
      if(DEBUG) then
        if(widthTermSource .ne. widthTermSource) then
          print*,'width: not a number:',sourceLocationVertex,refVertex
          print*,'wdith:',vectorS(:),vectorRef(:),theta
          call stopProgram( 'abort widthTermSource - NaN    ')
        endif
        !if(refVertex .eq. sourceLocationVertex .and. widthTermSource .ne. 0.0) print*,"width term source:",refVertex,widthTermSource        
      endif

      ! single cell source
      !if( refVertex .eq. sourceLocationVertex ) then
      !  widthTermSource=1.0d0
      !else
      !  widthTermSource=0.0d0
      !endif

      return
      end  

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
      integer:: refVertex,step,i
      real(WP):: forceAdjointSource,distance,factor,dx
      logical:: isNeighbor
      logical,parameter:: extendedSource = .false.
      logical,parameter:: gaussian       = .false.
      logical,parameter:: plateau        = .false.
      logical,parameter:: onlyNeighbors  = .false.
      
      ! calculate adjoint source term: based upon a single cell displacement
      if( refVertex .eq. adjointSourceVertex) then
        forceAdjointSource=adjointSource(2,step)
      else
        forceAdjointSource=0.0_WP
      endif

      ! extend adjoint source to multiple cells
      if( extendedSource ) then
        ! calculate adjoint source term: based upon a gaussian cell displacement      
        if( gaussian .or. plateau ) then
          call greatCircleDistance(vertices(refVertex,:),vertices(adjointSourceVertex,:),distance)
          distance=distance*EARTHRADIUS
          dx=averageCellDistance
          
          if( gaussian ) then
            ! factor determined empirically such that it's a nice curve to distances till 280 km (4*dx for grid level 6)
            factor= exp(-distance*distance/(300*dx))
          else
            ! plateau means we take the same value as at the adjointSource vertex location for all the others as well
            factor=1.0_WP
          endif

          ! extend to a multiple of the cell distances
          if( distance .le. 4.0*dx ) then
            forceAdjointSource=factor*adjointSource(2,step)  
            ! debug
            if(DEBUG) then
              if( step .eq. 1) print*,'force adjoint: ',refVertex    
            endif
          endif
        endif
                
        ! enlarge the source to its closest neighbor cells
        if( onlyNeighbors ) then
          isNeighbor=.false.
          do i=1,cellNeighbors(adjointSourceVertex,0)
            if( refVertex .eq. cellNeighbors(adjointSourceVertex,i) ) then
              isNeighbor=.true.
              exit
            endif
          enddo
          if( isNeighbor ) then
            forceAdjointSource=adjointSource(2,step)
          endif
        endif
      endif
      return
      end

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
      use precision; use verbosity; use cells
      implicit none
      integer:: sourceVertex, refVertex
      real(WP):: forceTermExact,time,timeTerm2,width,sourceLat,sourceLon,theta
      external::timeTerm2   
      real(WP):: vectorS(3),vectorRef(3),mu2,col2
            
      ! distance to source in radian
      call getVector(sourceLat,sourceLon,vectorS(1),vectorS(2),vectorS(3))
      vectorRef(:) = vertices(refVertex,:)
      call greatCircleDistance(vectorS,vectorRef,theta)      
            
      mu2=WidthParameterMu*WidthParameterMu
      col2 = theta*theta
      width = exp(-col2/(mu2+mu2))/mu2

      ! force
      forceTermExact = timeTerm2(time)*width
           
      !debug
      if(DEBUG) then
        if(theta .gt. PI) then
          print*,"forcetermExact:",time,forceTermExact,refVertex,theta,width
        endif
        ! check if Nan
        if( forceTermExact .ne. forceTermExact)then
          print*,'forceTerm:',forceTermExact
          print*,time,refVertex,sourceLat,sourceLon
          call stopProgram( 'forceTermExact is NaN      ')
        endif
      endif      
      return
      end     
             