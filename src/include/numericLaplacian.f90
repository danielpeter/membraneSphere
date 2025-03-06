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

!-----------------------------------------------------------------------
      function discreteLaplacian(vertex)
!-----------------------------------------------------------------------
! calculates the laplacian of a grid point given by index
!
! (formula used: Tape, 2003, Chap. 4, 4.13 )
      use displacements; use cells; use verbosity; use propagationStartup
      implicit none
      integer:: vertex,i,count
      real(WP):: discreteLaplacian,centerDistances(0:6),edgesLength(0:6),laplace,area,sumNeighbors, sumFactor
      
      !get spherical area of vertex's cell
      call getCellArea(vertex,area)
      !get arc distances from center to neighbors
      call getCellCenterDistances(vertex,centerDistances)
      !get cell edges length
      call getCellEdgesLength(vertex,edgesLength)
      
      !calculate laplacian approximation
      sumFactor = 0.0
      sumNeighbors = 0.0
      count = 0
      do i=1, cellNeighbors(vertex,0)
        !factor for own displacement
        if(abs(centerDistances(i)*area) .lt. 1.0e-9) then
          print*,'invalid size in calculation for laplacian'
          print*,'  vertex:',vertex
          print*,'  neighbor:',i
          print*,'  centerDistance:',centerDistances(i)
          print*,'  area:',area
          call stopProgram( 'discreteLaplacian - centerDistance too small   ')
        endif
        
        sumFactor= sumFactor+                                           &
     &     edgesLength(i)/(centerDistances(i)*area)
    
        !influence by neighbors displacement
        sumNeighbors = sumNeighbors +                                   &
     &      (edgesLength(i)/(centerDistances(i)*area))*                 &
     &        displacement(cellNeighbors(vertex,i)) 
          
        !check
        count = count +1
        if( count .gt. 6 ) call stopProgram( 'abort-discreteLaplacian neighbors   ')
        
        !debug
        if(DEBUG) then    
          if( vertex .eq. receiverVertex ) then
            print*,'vertex ',vertex,' neighbor:',i, sumFactor, sumNeighbors
            print*,'edge ',edgesLength(i),' centers ',centerDistances(i)
            print*,'displacement ',cellNeighbors(vertex,i),displacement(cellNeighbors(vertex,i))
          endif
        endif
          
      enddo
        
      laplace = sumNeighbors - sumFactor*displacement(vertex)
      discreteLaplacian = laplace
      
      !debug
      if(DEBUG) print*,'discreteLaplacian',discreteLaplacian

      return        
      end
      
!-----------------------------------------------------------------------
      function precalc_discreteLaplacian(vertex)
!-----------------------------------------------------------------------
! calculates the laplacian of a grid point given by index, 
! uses precalculated values for cell areas, cell center distances and cell edge lengths
!
! (formula used: Tape, 2003, Chap. 4, 4.13 )
! needs arrays cellAreas(), cellEdgesLength(),cellCenterDistances(),
! cellNeighbors() and displacement()
      use cells; use displacements; use verbosity; use propagationStartup
      implicit none
      integer:: vertex,i,count
      real(WP):: precalc_discreteLaplacian,area,areareci,sumNeighbors,sumFactor
      
      !get spherical area of vertex's cell
      area = cellAreas(vertex)
      areareci = 1.0_WP/area
            
      !calculate laplacian approximation
      sumFactor = 0.0_WP
      sumNeighbors = 0.0_WP
      count = 0
      do i=1, cellNeighbors(vertex,0)
        !factor for own displacement        
        sumFactor= sumFactor+cellEdgesLength(vertex,i)/cellCenterDistances(vertex,i)
    
        !influence by neighbors displacement
        sumNeighbors = sumNeighbors+cellEdgesLength(vertex,i)/cellCenterDistances(vertex,i)*displacement(cellNeighbors(vertex,i)) 
          
        !check
        count = count +1
        
        !debug
        if(DEBUG) then
          if( count .gt. 6 ) call stopProgram( 'abort-discreteLaplacian neighbors     ')                    
          !if( vertex .eq. receiverVertex) then
          !  print*,'vertex ',vertex,' neighbor:',i, sumFactor, sumNeighbors
          !  print*,'edge ',cellEdgesLength(vertex,i),' centers ',cellCenterDistances(vertex,i)
          !  print*,'displacement ',cellNeighbors(vertex,i),displacement(cellNeighbors(vertex,i))        
          !endif
        endif
      enddo
      
      sumFactor = sumFactor*displacement(vertex)*areareci
      sumNeighbors = sumNeighbors*areareci
      precalc_discreteLaplacian = sumNeighbors - sumFactor      
    
      !debug
      if(DEBUG) then
        ! check if NaN
        if( precalc_discreteLaplacian .ne. precalc_discreteLaplacian) then
          print*,'precalculated laplacian error:',precalc_discreteLaplacian,vertex,sumFactor,sumNeighbors,area,areareci,count
          call stopProgram('precalc_discreteLaplacian error     ')
        endif
      endif
      
      return        
      end

!-----------------------------------------------------------------------
      function precalc_backdiscreteLaplacian(vertex)
!-----------------------------------------------------------------------
! calculates the laplacian of a grid point given by index for backward calculation of previous forward displacement fields
!
! (formula used: Tape, 2003, Chap. 4, 4.13 )
! needs arrays cellAreas(), cellEdgesLength(),cellCenterDistances(),
! cellNeighbors() and backwarddisplacement()
      use cells; use propagationStartup; use adjointVariables; use verbosity
      implicit none
      integer:: vertex,i,count
      real(WP):: precalc_backdiscreteLaplacian,area,areareci,sumNeighbors,sumFactor
      
      !get spherical area of vertex's cell
      area = cellAreas(vertex)
      areareci = 1.0_WP/area
            
      !calculate laplacian approximation
      sumFactor = 0.0_WP
      sumNeighbors = 0.0_WP
      count = 0
      do i=1, cellNeighbors(vertex,0)
        !factor for own displacement
        sumFactor= sumFactor+cellEdgesLength(vertex,i)/cellCenterDistances(vertex,i)
    
        !influence by neighbors displacement
        sumNeighbors = sumNeighbors+cellEdgesLength(vertex,i)/cellCenterDistances(vertex,i)*backwardDisplacement(cellNeighbors(vertex,i)) 
          
        !check
        count = count +1        
      enddo
      
      sumFactor = sumFactor*backwardDisplacement(vertex)*areareci
      sumNeighbors = sumNeighbors*areareci
      precalc_backdiscreteLaplacian = sumNeighbors - sumFactor
      
      !debug
      if(DEBUG) then
        ! check if Nan
        if( precalc_backdiscreteLaplacian .ne. precalc_backdiscreteLaplacian ) then
          print*,'precalculated backwardlaplacian error:',precalc_backdiscreteLaplacian,vertex,sumFactor,sumNeighbors,area,areareci,count
          call stopProgram('precalc_backdiscreteLaplacian error     ')
        endif
      endif
      
      
      return        
      end
            
