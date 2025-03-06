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
      subroutine getCellArea(vertex,area)
!-----------------------------------------------------------------------
! calculates for vertex his cell area on unit sphere by summating
! the correspoinding triangles build by the reference
! vertex and his cell corners 
!
! uses haversine formula for excess calculation of area
!
! input:
!    vertex     - index of reference vertex
!    area      - cell area
!
! needs:
!    vertices - voronoi cell center vertices array
!    cellFace - index array of face vertices 
!    cellCorners - index array of cell corner vertices
!
! return: cell area in km^2
! (formula:http://mathforum.org/library/drmath/view/65316.html)
      use cells; use verbosity
      implicit none
      integer:: vertex
      real(WP):: area,vectorA(3),vectorB(3),vectorC(3),angleA,angleB,angleC,normalAB, normalAC
      integer:: n,k, corner1,corner2,corner3
      real(WP):: a,b,c,s,excess,tanexcess4,haversine
      
      !check index
      if( cellFace(vertex,0) .gt. 6 .or. cellFace(vertex,0) .lt. 5) then
        print*,'vertex',vertex,'faces',cellFace(vertex,0)
        call stopProgram( 'abort-cellArea cellFace    ')
      endif
        
      corner1=vertex
      vectorA(:) = vertices(corner1,:)
            
      area = 0.0
      do n=1, cellFace(vertex,0)
        !determine triangles corner indices
        corner2=cellFace(vertex,n)
        corner3=cellFace(vertex,(mod(n,cellFace(vertex,0))+1) )
        !debug
        if(DEBUG) print*,corner1,corner2,corner3

        !get corner vectors
        vectorB(:) = cellCorners(corner2,:)
        vectorC(:) = cellCorners(corner3,:)
        
        !using spherical sides of triangle
        a = haversine(vectorA, vectorB)
        b = haversine(vectorB, vectorC)
        c = haversine(vectorC, vectorA)         
        s = (a + b + c)*0.5_WP
        
        tanexcess4 = sqrt( tan(s*0.5_WP)*tan((s-a)*0.5_WP)*tan((s-b)*0.5_WP)*tan((s-c)*0.5_WP) )       
        excess = 4.0_WP*atan( tanexcess4)
        
        
        area = area+(EARTHRADIUS*EARTHRADIUS)*excess
      enddo
      
      if( area .lt. 1.0) then
        print*,'cellArea calculation has invalid area'
        print*,'  vertex:',vertex
        print*,'  infos:',n,a,b,c,s,excess,tanexcess4
      endif
      
      
      return
      end


!-----------------------------------------------------------------------
      subroutine getCellCenterDistances(vertex,distances)
!-----------------------------------------------------------------------
! calculates for vertex his arc distance to all his
! voronoi cell center neighbors on the unit sphere
!
! input:
!    vertex     - index of reference vertex
!    distances      - distances to cell center neighbors 
!
! needs:
!    vertices - voronoi cell center vertices array
!    cellNeighbors - index array of face neighbor vertices 
!
! return: distances array [in km] 
      use cells; use verbosity
      implicit none
      integer:: vertex,n,k
      real(WP):: distances(0:6), distance
      real(WP):: vectorA(3),vectorB(3)
      
      !cell center reference
      vectorA(:)=vertices(vertex, : )

      !check index
      if( cellNeighbors(vertex,0) .gt. 6 .or. cellNeighbors(vertex,0) .lt. 5) then  
        print*,'vertex',vertex
        call stopProgram( 'abort-cellDistances cellNeighbors    ')
      endif
      
      !get distances
      distance = 0.0
      distances(0) = cellNeighbors(vertex,0)
      do n=1, cellNeighbors(vertex,0)
        !determine distance to neighbor
        vectorB(:) = vertices(cellNeighbors(vertex,n), : )
        
        !calculate arc distance between these two
        call greatCircleDistance(vectorA,vectorB,distance)
        distances(n)=EARTHRADIUS*distance
        !debug
        if(DEBUG) print*,'vertex(',vertex,'),',n,': ',distances(n)        
      enddo

      end
            
      
!-----------------------------------------------------------------------
      subroutine getCellEdgesLength(vertex,lengths)
!-----------------------------------------------------------------------
! calculates for vertex his edge lengths corresponding
! to its neighbor cells
!
! input:
!    vertex     - index of reference vertex
!    lengths      - edges length
!
! needs:
!    vertices - voronoi cell center vertices array
!    cellFace - index array of face vertices 
!    cellCorners - index array of cell corner vertices
!
! return: lengths [in km]
      use cells; use verbosity
      implicit none
      integer:: vertex,n,k,corner1,corner2
      real(WP):: vectorA(3),vectorB(3),lengths(0:6),distance

      !check index
      if( cellFace(vertex,0) .gt. 6 .or. cellFace(vertex,0) .lt. 5) then
        print*,'vertex',vertex,'faces',cellFace(vertex,0)
        call stopProgram( 'abort-cellEdgesLength cellFace   ')  
      endif
                      
      !get cell edge lengths
      !write(1,*) '#',vertex
      lengths(0) = cellFace(vertex,0)
      do n=1, cellFace(vertex,0)        
        ! determine which corners (they are ordered clockwise)
        corner1=cellFace(vertex,n)
        corner2=cellFace(vertex, mod(n,cellFace(vertex,0))+1)
             
        !get corresponding corner vectors
        !do k=1,3
        !  vectorA(k) = cellCorners(corner1,k)
        !  vectorB(k) = cellCorners(corner2,k)
        !enddo
        
        vectorA(:) = cellCorners(corner1,:)
        vectorB(:) = cellCorners(corner2,:)
        
        !calculate arc distance between these two
        call greatCircleDistance(vectorA,vectorB,distance)
        lengths(n)=EARTHRADIUS*distance
        
        !print an output file for gnuplot       
        !write(1,'(3f12.8,f12)') vectorA(1),vectorA(2),vectorA(3),
         !&                             lengths(n)
        !debug
        if(DEBUG) print*,'vertex(',vertex,'),',n,': ',lengths(n),corner1,corner2
      enddo

      !end output with first vertex again
      !write(1,'(3f12.8,f12)') (vectorB(k),k=1,3), lengths(1)
      !write(1,*) '\n'
      
      end

!-----------------------------------------------------------------------
      subroutine getCellArea_alternate(vertex,area)
!-----------------------------------------------------------------------
! calculates for vertex his cell area on unit sphere by summating
! the correspoinding triangles build by the reference
! vertex and his cell corners 
!
! uses spherical angles to calculate the area (supposed to be less exact)
!
! input:
!    vertex     - index of reference vertex
!    area      - cell area
!
! needs:
!    vertices - voronoi cell center vertices array
!    cellFace - index array of face vertices 
!    cellCorners - index array of cell corner vertices
!
! return: cell area [in km^2]
! (formula:http://mathworld.wolfram.com/SphericalTrigonometry.html)
      use cells; use verbosity
      implicit none
      integer:: vertex,n,k, corner1,corner2,corner3
      real(WP):: area,vectorA(3),vectorB(3),vectorC(3),angleA,angleB,angleC,normalAB, normalAC

      !check index
      if( cellFace(vertex,0) .gt. 6 .or. cellFace(vertex,0) .lt. 5) then
        print*,'vertex',vertex,'faces',cellFace(vertex,0)
        call stopProgram( 'abort-cellArea cellFace    ')
      endif
        
      corner1=vertex
      vectorA(:) = vertices(corner1,:)
            
      area = 0.0
      do n=1, cellFace(vertex,0)
        !determine triangles corner indices
        corner2=cellFace(vertex,n)
        corner3=cellFace(vertex,(mod(n,cellFace(vertex,0))+1) )

        !debug
        if(DEBUG) print*,corner1,corner2,corner3

        !get corner vectors
        vectorB(:) = cellCorners(corner2,:)
        vectorC(:) = cellCorners(corner3,:)
        
        !calculate corner angles of spherical triangle
        call sphericalAngle(vectorA,vectorB,vectorC,angleA)        
        call sphericalAngle(vectorB,vectorA,vectorC,angleB)
        call sphericalAngle(vectorC,vectorA,vectorB,angleC)
        !debug
        if(DEBUG) then
          print*,'angleA',angleA/PI*180
          print*,'angleB',angleB/PI*180
          print*,'angleC',angleC/PI*180
        endif
        
        area = area+(EARTHRADIUS*EARTHRADIUS)*(angleA + angleB + angleC - PI) 
        
      enddo

      end

