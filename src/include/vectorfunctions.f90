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
      subroutine midpoint(vector1, vector2, middle)
!-----------------------------------------------------------------------
! calculates the point in the middle of the two endpoints marked by the
! given vectors 1 & 2
! return: vector to middle point
      use precisions
      implicit none
      real(WP):: vector1(3),vector2(3),middle(3)
       
      middle(1) = 0.5*(vector1(1) + vector2(1))
      middle(2) = 0.5*(vector1(2) + vector2(2))
      middle(3) = 0.5*(vector1(3) + vector2(3))
      end
      
!-----------------------------------------------------------------------
      subroutine getCentralpoint(vector1, vector2, vector3, central)
!-----------------------------------------------------------------------
! calculates the point in the center of the three endpoints marked by
! the given vectors 1, 2 & 3
!
! (see also: Carl's thesis, eq. (4.3) on page 38 )
! return: vector to central point
      use precisions
      implicit none
      real(WP):: vector1(3),vector2(3),vector3(3),central(3),vTmpA(3),vTmpB(3),vlength
      
      ! using the formula with one vectorproduct
      call vectordiff(vector3,vector1,vTmpA)
      call vectordiff(vector2,vector1,vTmpB)
      call vectorproduct(vTmpA,vTmpB,central)
      
      ! normalize vector
      !call scalarproduct(central,central,vlength)
      !vlength = sqrt(vlength)      
      !call vectorscale(central,vlength,central)
      
      vlength = sqrt(central(1)**2 + central(2)**2 + central(3)**2 )
      central(:) = central(:) / vlength
      
      ! compare with this calculation
      !central(1) = (vector1(1) + vector2(1)+vector3(1))/3.0
      !central(2) = (vector1(2) + vector2(2)+vector3(2))/3.0
      !central(3) = (vector1(3) + vector2(3)+vector3(3))/3.0
      
      end
      
!-----------------------------------------------------------------------
      subroutine vectorproduct(vector1, vector2, product)
!-----------------------------------------------------------------------
! calculates the vector product of the
! given vectors 1 & 2
! return: vector product
      use precisions
      implicit none
      real(WP):: vector1(3),vector2(3),product(3)
      
      product(1) = vector1(2)*vector2(3) - vector1(3)*vector2(2)
      product(2) = - vector1(1)*vector2(3) + vector1(3)*vector2(1)
      product(3) = vector1(1)*vector2(2) - vector1(2)*vector2(1)
      end

!-----------------------------------------------------------------------
      subroutine scalarproduct(vector1, vector2, scalar)
!-----------------------------------------------------------------------
! calculates the scalar product of the
! given vectors 1 & 2
! return: scalar product
      use precisions
      implicit none
      real(WP):: vector1(3),vector2(3),scalar
      
      scalar=vector1(1)*vector2(1)+vector1(2)*vector2(2)+vector1(3)*vector2(3)
      end
      
!-----------------------------------------------------------------------
      subroutine vectordiff(vector1, vector2, difference)
!-----------------------------------------------------------------------
! calculates the difference vector of the
! given vectors: vector1 - vector2
!
! input:
!    vector1,vector2  - vectors to subtract
!    difference          - vector difference
!
! return: difference vector
      use precisions
      implicit none
      real(WP):: vector1(3),vector2(3),difference(3)
      
      difference(1) = vector1(1)-vector2(1)
      difference(2) = vector1(2)-vector2(2)
      difference(3) = vector1(3)-vector2(3)
      end

!-----------------------------------------------------------------------
      subroutine vectoradd(vector1, vector2, addition)
!-----------------------------------------------------------------------
! calculates the addition vector of the
! given vectors 1 & 2
!
! input:
!    vector1, vector2 - vectors to add
!    addition             - resulting vector
!
! return: addition vector
      use precisions
      implicit none
      real(WP):: vector1(3),vector2(3),addition(3)
      
      addition(1) = vector1(1)+vector2(1)
      addition(2) = vector1(2)+vector2(2)
      addition(3) = vector1(3)+vector2(3)
      end
      
!-----------------------------------------------------------------------
      subroutine vectorscale(vector1, scale, scaledvector)
!-----------------------------------------------------------------------
! multiplies the vector with a given scale factor
! 
! return: scaled vector
      use precisions
      implicit none
      real(WP):: vector1(3),scaledvector(3),scale
      
      scaledvector(1) = scale*vector1(1)
      scaledvector(2) = scale*vector1(2)
      scaledvector(3) = scale*vector1(3)
      end
      
!-----------------------------------------------------------------------
      subroutine normalize(vector)
!-----------------------------------------------------------------------
! normalizes vector to unit length 1
! 
! return: normalized vector
      use precisions
      implicit none
      real(WP):: vector(3), vectorlength
      
      call vectorscale(vector,1.0_WP/vectorlength(vector),vector)
      end      
      
!-----------------------------------------------------------------------
      function distanceToSurface(vector1,vector2,vector3,vCheck)
!-----------------------------------------------------------------------
! calculates distance of point vCheck to surface spanned by first 
! 3 vectors
!
! returns: distance
! check formula: http://mathworld.wolfram.com/Point-PlaneDistance.html
      use precisions
      implicit none
      real(WP):: vector1(3),vector2(3),vector3(3), vCheck(3)
      real(WP):: normale(3), tmpA(3), tmpB(3),tmpC(3)
      real(WP):: distanceToSurface,length, dist
      
      ! get the unit normal of the surface spanned by the first 3 vector
      call vectordiff(vector2,vector1,tmpA)
      call vectordiff(vector3, vector1,tmpB)
      call vectorproduct(tmpA,tmpB, normale)
      call scalarproduct(normale,normale,length)
      length = sqrt(length)
      call vectorscale( normale, 1.0_WP/length, normale)
      
      ! get the distance of the point from vCheck
      call vectordiff(vector1, vCheck, tmpC)
      call scalarproduct(normale,tmpC,dist)
      
      distanceToSurface = dist
      return
      end
      
!-----------------------------------------------------------------------
      function distance(vector1, vector2)
!-----------------------------------------------------------------------
! calculates the euclidian distance between points
! given by the two vectors
!
! return: distance
      use precisions
      implicit none
      real(WP):: vector1(3),vector2(3),vTmp(3),distance
      
      ! get vector length between these two	      
      call vectordiff(vector1, vector2, vTmp)
      call scalarproduct(vTmp,vTmp,distance)
      distance = sqrt(distance)      
      return
      end 
	
!-----------------------------------------------------------------------
      function vectorlength(vector)
!-----------------------------------------------------------------------
! calculates the length of the given vector
! in euclidian norm
! input: 
!   vector - vector on which the length is calculated
!
! return: vectorlength
      use precisions
      implicit none
      real(WP):: vectorlength,vector(3),scalar
      
      ! get scalarproduct
      scalar=vector(1)*vector(1)+vector(2)*vector(2)+vector(3)*vector(3)
      vectorlength = sqrt(scalar)      
      return
      end 

            
!-----------------------------------------------------------------------
      subroutine pentagonCentral(ordered,vectors, central)
!-----------------------------------------------------------------------
! calculates the central point of a given regular pentagon.
! the vectors constituting the pentagon can be 
! ordered or not.
!
! for calculating the central point a middle-line is created
! from a reference point to the middle-point of the distant
! pentagon side.
! the central point lies on this middle-line where the distance
! uses formula: d = a/2 * 1/sin(36°)
!                [R=1/10 * a * sqrt( 50+10*sqrt(5) )]
!			where a is the pentagon side length
!
! input parameters:
!   ordered - zero indicates unordered, one ordered (either
!                    clockwise or counter-clockwise) edge vectors
!   vectors   - holds the edge vectors  
!   central    - vector which takes the central point vector     
!
! return: central vector
       use precisions; use verbosity
       implicit none
       integer:: ordered        
       real(WP):: vectors(5,3),central(3),vRef(3),vTmpA(3), vTmpB(3), vTmpC(3)
       real(WP):: sideMiddle(3),side, minside, maxside,distance, factor
       integer:: farNeighbors(2), count,i,k

       !factor = 1.0_WP/(1.0_WP+0.80901699437495_WP)
       factor = 0.5527864045000412_WP
       
       ! set first point as reference point       
       do k=1,3      
         vRef(k) = vectors(1,k)
       enddo      

      ! determine the distant pentagon side edges
       if( ordered .eq. 1) then
          farNeighbors(1) = 3
          farNeighbors(2) = 4
       else  
          ! get side length of pentagon and far neighbors distance       
          minside = 100000.0
          maxside= -100000.0
          do i=2,5
            do k=1,3
              vTmpA(k) = vectors(i,k)
            enddo
            
            ! get vector length between these two	      
            side = distance(vRef,vTmpA)
    
            !determine min and max of lengths
            if( side .gt. maxside) then
                    maxside = side
            endif
            
            if( side .lt. minside) then
                    minside = side
            endif
          enddo  
          ! get the points which are far neighbors
          count = -1
          do i=2,5
            do k=1,3
              vTmpA(k) = vectors(i,k)
            enddo
          
            if( abs(distance(vRef,vTmpA) -maxside ) .lt. 1.0e-4) then
              count = count+1
              farNeighbors(mod(count,2)+1) = i
            endif
          enddo
          ! check if we have the two points
          if( count .lt. 1) then
            print *,'abort condition met pentagonCentral neighbors!'
            stop 'abort - pentagonCentral'
          endif
       endif !ordered
              
       !debug
       !print *,'far neighbors ', farNeighbors(1), farNeighbors(2)
        
       ! get middle point of far side
       do i=1,3
       	 vTmpB(i) = vectors( farNeighbors(1), i)
         vTmpC(i) = vectors( farNeighbors(2), i)
       enddo      
       call midpoint(vTmpB,vTmpC,sideMiddle)

       !debug 
       ! write(*,'(3f12.8)') sideMiddle(1),sideMiddle(2),sideMiddle(3)
       
       ! calculate the central point
       call vectordiff(sideMiddle,vRef,central) !central points R->middleside
       call vectorscale(central, factor, central)    !central points R->center 
       call vectoradd(vRef, central, central)     !central point Orig->center

       end

!-----------------------------------------------------------------------
      subroutine projectToSphere(vertices,nCount, vector)
!-----------------------------------------------------------------------
! projects the given vector on to the unit sphere , 
! normal to the determined plane given by vertices
!
!  input:
!    vertices - vertex vectors of plane
!    nCount - number of vertices
!    vector   - vector to project
!
! return: projected vector
      use precisions; use verbosity
      implicit none
      integer:: k,nCount
      real(WP):: vertices(nCount,3), vector(3), vectorlength, vNormale(3)
      real(WP):: vTmpA(3),vTmpB(3),vTmpC(3), vProj(3), factor1,factor2, A,B,C,det

      ! be sure vector is not already on sphere
      if( abs(vectorlength( vector) - 1.0) .lt. 1.0e-5) then
        print *,'already projected'
        return
      endif      

      ! get plane normale
      do k=1,3
        vTmpA(k) = vertices(1,k)
        vTmpB(k) = vertices(2,k)
        vTmpC(k) = vertices(3,k)
      enddo
      call vectordiff(vTmpB,vTmpA,vTmpB)
      call vectordiff(vTmpC,vTmpA,vTmpC)
      call vectorproduct(vTmpB,vTmpC, vNormale)
      call vectorscale( vNormale, 1.0_WP/vectorlength(vNormale),         &
     &                         vNormale)

      ! find mulitplication factor of normale to project on unit sphere
      ! (v+f*n)² = 1 = (v(1)+f*n(1))²+(v(2)+f*n(2))²+(v(3)+f*n(3))²
      !                    =  (v)² +2*f*(v(1)n(1)+v(2)n(2)+v(3)n(3)) +f² * (n)² 
      !  (n)²*f² + 2(v*n)*f +(v)²-1 = 0
      !  f =  (-2(v*n) +- det) / (2*(n)²) 
      !       with det= sqrt((2v*n)²-4*(n)²((v)²-1))
      call scalarproduct( vNormale,vNormale,A)
      call scalarproduct( vector,vNormale,B)
      B= 2*B
      call scalarproduct( vector, vector,C)
      C = C - 1
      det = B**2 - 4*A*C
      !check
      if( det .lt. 0.0) then
        print *,'complex calculation found, abort...'
        stop 'abort - projectToSphere'
      endif
      det = sqrt(det)

      factor1= (-B+det)/(2*A)
      factor2= (-B-det)/(2*A)      
      
      ! project on closer side of the sphere
      if( abs(factor1) .lt. abs(factor2) ) then
        call vectorscale( vNormale, factor1,vNormale)
      else
        call vectorscale( vNormale, factor2,vNormale)
      endif

      ! get resulting vector
      call vectoradd(vector, vNormale, vector)

      end

!-----------------------------------------------------------------------
      integer function getClosestPoint(vertices,vDestination)
!-----------------------------------------------------------------------
! finds the closest vertex to the given vector
!
! input:
!     vertices - vertices array
!     vDestination - vector of desired location
!
! return: index of vertex in vertices array which is closest 
!            to the desired location
      use precisions; use verbosity
      implicit none
      real(WP):: vertices(32,3), vDestination(3),vTmp(3),distance, minimum, tmpDist
      integer:: n, k, vertex
      
      vertex = -1
      minimum = 10000.0
      do n=1,32
        ! get vertex vector
        do k=1,3
          vTmp(k)=vertices(n,k)
        enddo
      
        ! get minimal distance between two points
        tmpDist = distance(vTmp, vDestination)
        if( tmpDist .lt. minimum) then
          minimum = tmpDist 
          vertex = n
        endif      
      enddo
      
      getClosestPoint = vertex      
      return
      end
      
                       
!----------------------------------------------------------------------
      subroutine getSphericalCoord(vertex,colatitude,longitude) 
!---------------------------------------------------------------------
! translates the cartesian coordinates of triangle vertex into
! spherical colatitude [0,PI] and longitude [0,2PI[
!
! input:
!    vertex - index of cell center
!    colatitude, longitude
!
! returns: colatitude between [0,PI], longitude between [0,2PI[ in radian
      use cells
      implicit none
      integer:: vertex,k
      real(WP):: colatitude, longitude,vectorlength,vector(3)
      
      !get corresponding vertex vector
      do k=1,3
        vector(k) = vertices(vertex,k)
      enddo
            
      !check if point is on sphere
      if( abs(vectorlength(vector) - 1.0) .gt. 1.0e-4 ) then
         print*,vertex,'vector not normalized ', vectorlength(vector)
         stop 'abort-getSphericalCoord length'
      endif
     
      !calculate spherical coordinates
      colatitude = acos( vector(3) )
      !note: atan2 returns between [-PI,+PI]
      longitude = atan2(vector(2),vector(1))
      if(longitude .lt. 0.0) longitude = longitude+2*PI
      
      !check
      if( longitude .lt. 0.0) stop 'abort-getSphericalCoord longitude'
      
      !debug
      !print*,vertex,' vector: ',vector(1),vector(2),vector(3)
      !print*,'colatitude',colatitude
      !print*,'longitude',longitude
      
      end
      
!----------------------------------------------------------------------
      subroutine getSphericalCoord_Lat(vertex,latitude,longitude) 
!---------------------------------------------------------------------
! translates the cartesian coordinates of triangle vertex into
! spherical latitude [90¼,-90¼] and longitude [0¼,360¼[
!
! (uses vector array vertices)
! (90¼ latitude corresponds to north pole, -90¼ to sourth pole)
!
! input:
!    vertex - index of cell center
!    latitude, longitude
!
! returns: latitude/longitude in degrees
      use cells
      implicit none
      integer::  vertex,k
      real(WP):: latitude,longitude,vectorlength,vector(3)
      
      !get corresponding vertex vector
      vector(:) = vertices(vertex,:)
            
      !check if point is on sphere
      !if( abs(vectorlength(vector) - 1.0) .gt. 1.0e-4 ) then
      !   print*,vertex,'vector not normalized ', vectorlength(vector)
      !   stop 'abort-getSphericalCoord_Lat length'
      !endif
     
      !calculate spherical coordinates
      call getSphericalCoordinates(vector,latitude,longitude)      
      end
                                  
                                                        
!-----------------------------------------------------------------------
      subroutine getSphericalCoordinates(vector,latitude,longitude)
!-----------------------------------------------------------------------
! converts point with vector( in cartesian coordinates (x,y,z) ) to
! spherical coordinates (latitude,longitude)
!
! input:
!    vector    - vector in cartesian coordinates of point
!    latitude, longitude - spherical coordinates to determine
! return: latitude [-90,90], longitude [0, 360[ in degree
      use precisions
      implicit none
      real(WP):: vector(3),x,y,z,latitude,longitude
      
      ! initialize
      x=vector(1)
      y=vector(2)
      z=vector(3)
      
      ! determine longitude and latitude      
      longitude = atan2(y,x)
      latitude = asin(z)
          
      ! convert to degrees    
      longitude = RAD2DEGREE*longitude
      latitude  = RAD2DEGREE*latitude
      
      ! return only positive longitudes
      if( longitude < 0.0) then
        longitude = longitude + 360.0_WP        
      endif            
      end subroutine
      
!-----------------------------------------------------------------------
      subroutine getVector(latitude,longitude,x,y,z)
!-----------------------------------------------------------------------
! converts point with spherical coordinates (latitude,longitude) 
! to cartesian coordinates at (x,y,z)
! http://mathworld.wolfram.com/SphericalCoordinates.html
!
! input:
!    x,y,z    - cartesian coordinates of point
!    latitude, longitude - spherical coordinates (lat[-90/90] lon[0/360] in degree)
! return: x,y,z
      use precisions
      implicit none
      real(WP):: x,y,z,latitude,longitude,colat,lon
      
      lon = DEGREE2RAD*longitude
      colat = PI/2.0_WP - DEGREE2RAD*latitude
      
      x=cos(lon)*sin(colat)
      y=sin(lon)*sin(colat)
      z=cos(colat)      
      end      
      
      
!-----------------------------------------------------------------------
      subroutine sphericalAngle(vectorA,vectorB,vectorC,                &
     &                                          angle)
!-----------------------------------------------------------------------
! calculates the angle in the corner of a spherical triangle
!
! input:
!    vectorA,vectorB,vectorC - triangle corner vectors
!    angle  - corner angle at A
!
! return: angles in radian
! (formula:http://mathworld.wolfram.com/SphericalTrigonometry.html)  
      use precisions
      implicit none
      real(WP):: vectorA(3),vectorB(3),vectorC(3), angle,normalAB(3),normalAC(3),normalBC(3)
      real(WP):: normalResult(3), vectorlength,scalar
      
      !get cross products
      call vectorproduct(vectorA,vectorB, normalAB)
      call vectorproduct(vectorA,vectorC, normalAC)
      call vectorproduct(vectorB,vectorC, normalBC)
      call vectorproduct(normalAB,normalAC,normalResult)
      !calculate angle at A
      scalar = vectorlength(normalResult)/(vectorlength(normalAB)*vectorlength(normalAC))
      angle = asin( scalar )
     
      end       
           
!-----------------------------------------------------------------------
      function haversine(vector1,vector2)
!-----------------------------------------------------------------------
! calculates the exact great circle distance of two points given 
! by the vectors on a unit sphere (radius r=1.0)
!
! input:
!    vector1,vector2 - vectors
!
! return: haversine - great circle distance on unit sphere
! (formula:http://mathforum.org/library/drmath/view/51879.html)    
      use precisions
      implicit none
      real(WP):: vector1(3),vector2(3),dlon,dlonhalf,dlat,dlathalf,lon1,lon2,lat1,lat2
      real(WP):: haversine,a
      
      ! get longitude and latitudes
      lon1 = atan2(vector1(2),vector1(1))
      lat1 = asin(vector1(3))
      lon2 = atan2(vector2(2),vector2(1))
      lat2 = asin(vector2(3))
      !print*,'vector1:',lat1,lon1
      !print*,'vector2:',lat2,lon2
      
      !haversine formula
      dlon = lon2 - lon1
      dlonhalf= dlon*0.5_WP
      dlat = lat2 - lat1
      dlathalf = dlat*0.5_WP
      a = sin(dlathalf)*sin(dlathalf)+cos(lat1)*cos(lat2)*sin(dlonhalf)*sin(dlonhalf)
      haversine = 2.0_WP*atan2(sqrt(a), sqrt(1.0_WP-a)) 
      return
      end
      
      
!-----------------------------------------------------------------------
      subroutine findVertex(desiredLat,desiredLon,vertex)
!-----------------------------------------------------------------------
! finds nearest vertex for desired location (lat/lon in degrees) in
! the vertices array   
!
! input:
!     desiredLat/desiredLon - location to look for (in degree)
!     vertex                          - index of nearest vertex in vertices array
!
! return: vertex which is nearest
      use cells; use propagationStartup
      implicit none
      real(WP):: lat,lon,desiredLat,desiredLon,vectorS(3),vectorV(3),distanceMin,distance
      integer:: vertex,vertexMin,n
      
      !pick vector on unit sphere      
      lat = DEGREE2RAD*desiredLat
      lon = DEGREE2RAD*desiredLon
      vectorS(1) = cos(lat)*cos(lon)
      vectorS(2) = cos(lat)*sin(lon)
      vectorS(3) = sin(lat)
      
      ! determine minimal distance
      distanceMin = 1e+6      
      vertexMin = 0      
      do n=1,numVertices
        vectorV(:) = vertices(n,:)
        call greatCircleDistance(vectorS,vectorV,distance)
        if( distance .lt. distanceMin) then
          distanceMin = distance
          vertexMin = n
        endif
        
        ! for a distance very close to the vertex we quit
        if( distance .lt. distanceQuit ) then
          vertex = n
          return
        endif
      enddo
      vertex = vertexMin      
      end
      
      
!-----------------------------------------------------------------------
      subroutine rotateVector(rot,vector,rotVector)
!-----------------------------------------------------------------------
! calculates the rotation of given vector by rotation matrix: u = R * v
! 
! input:
!     rot       -   rotation matrix
!     vector -    vector to rotate
!     rotVector - rotated vector
!
! returns: rotated vector
      use precisions
      implicit none
      real(WP):: vector(3),rotVector(3),rot(3,3),tmp(3)
      integer:: i,j
      
      ! matrix - vector multiplication
      do i=1,3
        tmp(i)=0.0_WP
        do j=1,3
          tmp(i)= tmp(i)+rot(i,j)*vector(j)
        enddo
      enddo
      
      rotVector=tmp      
      end
      
!-----------------------------------------------------------------------
      subroutine matrixMultip(mat1,mat2, mat)
!-----------------------------------------------------------------------
! calculates the matrix multiplication
! returns: rotated vector
      use precisions
      implicit none
      real(WP):: mat1(3,3),mat2(3,3),mat(3,3),tmp(3,3)
      integer:: i,j,k
      
      ! matrix - matrix multiplication
      do i=1,3
        do j=1,3
          tmp(i,j)=0.0_WP
          do k=1,3
             tmp(i,j)=tmp(i,j)+mat1(i,k)*mat2(k,j)  
          enddo
        enddo
      enddo
      
      mat=tmp
      
      end
      
!-----------------------------------------------------------------------
      subroutine getRotationMatrix(vsource,vreceiver,rotationMatrix)
!-----------------------------------------------------------------------
! calculates the matrix for rotating the frame (1,0,0)/(0,1,0)/(0,0,1) to source/receiver
! (determines euler angles for rotation source to equator and receiver into equatorial plane)
! input:
!       vsource,vreceiver  - source/receiver vectors
!       rotationMatrix        - matrix to return
!
! returns: rotation matrix
      use precisions
      implicit none
      real(WP):: vsource(3),vreceiver(3),rotationMatrix(3,3)
      real(WP):: Vnormal1(3),Vnormal2(3),tauX,tauY,tauZ,vectorlength
      external:: vectorlength
      
      ! determine euler angles for rotation source to equator and receiver into equatorial plane
      ! get normale on plane of source and receiver
      call vectorproduct(vsource,vreceiver,Vnormal1)
      call vectorscale(Vnormal1,1.0_WP/vectorlength(Vnormal1),Vnormal1)
      
      ! get normale on plane of source and normal 1
      call vectorproduct(vsource,Vnormal1,Vnormal2)
      call vectorscale(Vnormal2,-1.0_WP/vectorlength(Vnormal2),Vnormal2)
            
      tauX = asin(-Vnormal2(3))
      tauY = atan2(vsource(3),Vnormal1(3))
      tauZ = atan2(Vnormal2(1),Vnormal2(2))
      
      !print*,'tau X:',tauX, tauX/PI*180.0
      !print*,'tau Y:',tauY, tauY/PI*180.0
      !print*,'tau Z:',tauZ, tauZ/PI*180.0
      
      ! rotation of frame (1,0,0)/(0,1,0)/(0,0,1) to new frame
      rotationMatrix(1,1)= cos(tauZ)*cos(tauY)+sin(tauZ)*sin(tauX)*sin(tauY)
      rotationMatrix(1,2)= sin(tauZ)*cos(tauX)
      rotationMatrix(1,3)= -cos(tauZ)*sin(tauY)+sin(tauZ)*sin(tauX)*cos(tauY)
      rotationMatrix(2,1)=-sin(tauZ)*cos(tauY)+cos(tauZ)*sin(tauX)*sin(tauY)
      rotationMatrix(2,2)=cos(tauZ)*cos(tauX)
      rotationMatrix(2,3)=sin(tauZ)*sin(tauY)+cos(tauZ)*sin(tauX)*cos(tauY)
      rotationMatrix(3,1)=cos(tauX)*sin(tauY)
      rotationMatrix(3,2)=-sin(tauX)
      rotationMatrix(3,3)=cos(tauX)*cos(tauY)

      
      end      
      
      
!-----------------------------------------------------------------------
      subroutine getInverseRotationMatrix(vsource,vreceiver,invRotationMatrix)
!-----------------------------------------------------------------------
! calculates the matrix for rotating the frame source/receiver to (1,0,0)/(0,1,0)/(0,0,1)
!
! input:
!       vsource,vreceiver  - source/receiver vectors
!       invRotationMatrix        - matrix to return
!
! returns: inverse rotation matrix
      use precisions
      implicit none
      real(WP):: vsource(3),vreceiver(3),invRotationMatrix(3,3)
      real(WP):: Vnormal1(3),Vnormal2(3),tauX,tauY,tauZ,vectorlength
      external:: vectorlength
      
      ! determine euler angles for rotation source to equator and receiver into equatorial plane
      ! get normale on plane of source and receiver
      call vectorproduct(vsource,vreceiver,Vnormal1)
      call vectorscale(Vnormal1,1.0_WP/vectorlength(Vnormal1),Vnormal1)
      
      ! get normale on plane of source and normal 1
      call vectorproduct(vsource,Vnormal1,Vnormal2)
      call vectorscale(Vnormal2,-1.0_WP/vectorlength(Vnormal2),Vnormal2)
      
      tauX = asin(-Vnormal2(3))
      tauY = atan2(vsource(3),Vnormal1(3))
      tauZ = atan2(Vnormal2(1),Vnormal2(2))
      
      !print*,'tau X:',tauX, tauX/PI*180.0
      !print*,'tau Y:',tauY, tauY/PI*180.0
      !print*,'tau Z:',tauZ, tauZ/PI*180.0
      
      ! inverse rotation matrix
      invRotationMatrix(1,1)= cos(tauZ)*cos(tauY)+sin(tauZ)*sin(tauX)*sin(tauY)
      invRotationMatrix(2,1)= sin(tauZ)*cos(tauX)
      invRotationMatrix(3,1)= -cos(tauZ)*sin(tauY)+sin(tauZ)*sin(tauX)*cos(tauY)
      invRotationMatrix(1,2)=-sin(tauZ)*cos(tauY)+cos(tauZ)*sin(tauX)*sin(tauY)
      invRotationMatrix(2,2)=cos(tauZ)*cos(tauX)
      invRotationMatrix(3,2)=sin(tauZ)*sin(tauY)+cos(tauZ)*sin(tauX)*cos(tauY)
      invRotationMatrix(1,3)=cos(tauX)*sin(tauY)
      invRotationMatrix(2,3)=-sin(tauX)
      invRotationMatrix(3,3)=cos(tauX)*cos(tauY)    
      
      end            
      
!-----------------------------------------------------------------------
      subroutine newDeltaLocation(latDelta,lonDelta,sourceVertex,receiverVertex,dlat,dlon)
!-----------------------------------------------------------------------
! calculates the delta location with reference to the midpoint of source and receiver
! http://astronomy.swin.edu.au/~pbourke/projection/eulerangle/

! input:
!   latDelta,lonDelta                     - delta location increments with reference to midpoint of source&receiver
!   sourceVertex,receiverVertex - location indices of source and receiver
!   dlat,dlon                                 - new delta location
! return: new delta location dlat,dlon
      use cells
      implicit none
      integer:: sourceVertex,receiverVertex
      real(WP):: latDelta,lonDelta,dlat,dlon,Vmid(3),Vsource(3),Vreceiver(3),Vnormal1(3),Vnormal2(3),vtmp(3)
      real(WP):: vectorlength,tauX,tauY,tauZ,rot(3,3),rotInv(3,3),lat,lon,tmp,rotX(3,3),rotY(3,3),rotZ(3,3)
      
      Vsource(:)=vertices(sourceVertex,:)
      Vreceiver(:)=vertices(receiverVertex,:)
      
      ! get midpoint between source and receiver
      !call greatCircleMidpoint(Vsource,Vreceiver,Vmid)      
      call getRotationMatrix(Vsource,Vreceiver,rot)
              
      ! inverse rotation
      call getInverseRotationMatrix(Vsource,Vreceiver,rotInv)
            
      ! determine delta vector
      call getVector(latDelta,lonDelta,vtmp(1),vtmp(2),vtmp(3))
      
      call getSphericalCoordinates(vtmp,lat,lon)
      print*,'  delta location:',lat,lon
      
      ! rotate delta to new frame (source/../..)
      call rotateVector(rot,vtmp,vtmp)
      
      call getSphericalCoordinates(vtmp,lat,lon)
      print*,'  rotated delta location',lat,lon
      
      ! set return values
      dlat=lat
      dlon=lon
      
      end

!-----------------------------------------------------------------------
      subroutine newDeltaLocation_alt(latDelta,lonDelta,dlat,dlon)
!-----------------------------------------------------------------------
! calculates the delta location with reference to the midpoint of source and receiver
! http://astronomy.swin.edu.au/~pbourke/projection/eulerangle/

! input:
!   latDelta,lonDelta - delta location increments with reference to midpoint of source&receiver
!   dlat,dlon             - new delta location
! return: new delta location dlat,dlon
      use cells
      implicit none
      integer:: midpoint,sourceVertex,receiverVertex
      real(WP):: latDelta,lonDelta,dlat,dlon
      real(WP):: Vmid(3),Vsource(3),Vreceiver(3),Vnormal1(3),Vnormal2(3),vtmp(3)
      real(WP):: vectorlength,tauX,tauY,tauZ,rot(3,3),rotInv(3,3),lat,lon,tmp,rotX(3,3),rotY(3,3),rotZ(3,3)
      real(WP):: sourceLat,sourceLon,receiverLat,receiverLon
      
      sourceLat=0.0
      sourceLon=60.0
      receiverLat=60.0
      receiverLon=150.0
      call findVertex(sourceLat,sourceLon,sourceVertex)      
      call findVertex(receiverLat,receiverLon,receiverVertex)
      
      Vsource(:)=vertices(sourceVertex,:)
      Vreceiver(:)=vertices(receiverVertex,:)
      
      ! get midpoint between source and receiver
      call greatCircleMidpoint(Vsource,Vreceiver,Vmid)
      
      ! determine euler angles for rotation source to equator and receiver into equatorial plane
      ! get normale on plane of source and receiver
      call vectorproduct(Vsource,Vreceiver,Vnormal1)
      call vectorscale(Vnormal1,1.0/vectorlength(Vnormal1),Vnormal1)
      
      ! get normale on plane of source and normal 1
      call vectorproduct(Vsource,Vnormal1,Vnormal2)
      call vectorscale(Vnormal2,-1.0/vectorlength(Vnormal2),Vnormal2)
      
      call getSphericalCoordinates(Vsource,lat,lon)
      !print*,'source:',lat,lon
      call getSphericalCoordinates(Vnormal1,lat,lon)
      !print*,'normal 1:',lat,lon
      call getSphericalCoordinates(Vnormal2,lat,lon)
      !print*,'normal 2:',lat,lon
      
      tauX = asin(-Vnormal2(3))
      tauY = atan2(Vsource(3),Vnormal1(3))
      tauZ = atan2(Vnormal2(1),Vnormal2(2))
      
      ! rotation of frame (1,0,0)/(0,1,0)/(0,0,1) to new frame
      rot(1,1)= cos(tauZ)*cos(tauY)+sin(tauZ)*sin(tauX)*sin(tauY)
      rot(1,2)= sin(tauZ)*cos(tauX)
      rot(1,3)= -cos(tauZ)*sin(tauY)+sin(tauZ)*sin(tauX)*cos(tauY)
      rot(2,1)=-sin(tauZ)*cos(tauY)+cos(tauZ)*sin(tauX)*sin(tauY)
      rot(2,2)=cos(tauZ)*cos(tauX)
      rot(2,3)=sin(tauZ)*sin(tauY)+cos(tauZ)*sin(tauX)*cos(tauY)
      rot(3,1)=cos(tauX)*sin(tauY)
      rot(3,2)=-sin(tauX)
      rot(3,3)=cos(tauX)*cos(tauY)
      
      call getVector(latDelta,lonDelta,vtmp(1),vtmp(2),vtmp(3))
      call getSphericalCoordinates(vtmp,lat,lon)
      !print*,'  delta location:',lat,lon
      call rotateVector(rot,vtmp,vtmp)
      call getSphericalCoordinates(vtmp,lat,lon)
      !print*,'  rotated delta location',lat,lon
      dlat=lat
      dlon=lon
      
      end
      
!-----------------------------------------------------------------------
      subroutine placeSecondDelta(lat,lon,secondlat,secondlon,secondvertex,distance)
!-----------------------------------------------------------------------
! places a second delta location next to the available one
!
! input:
!   lat,lon                        - location to reference
!   secondlat,secondlon - second location
!   secondvertex            - vertex index of second location
!
! return: new SecondLat,SecondLon,SecondVertex for second delta location
      use precisions
      implicit none
      real(WP):: lat,lon,secondlat,secondlon,distance
      integer:: secondvertex
      
      ! new location 
      secondlon = lon+distance
      secondlat  = lat
            
      ! get nearest vertex for new location
      call findVertex(secondlat,secondlon,secondvertex)
      end
      
!-----------------------------------------------------------------------
      subroutine normalizeArray(array,arraylength)
!-----------------------------------------------------------------------
! normalizes a vector array by its maximum value
! 
! input:
!   array           -  1xn array
!   arraylength - n length of array
! 
! returns: normalized array
      use precisions
      implicit none
      integer:: arraylength,i
      real(WP):: array(arraylength),maximum
      
      ! determine maximum value in array
      maximum=0.0_WP
      do i=1,arraylength
        if( abs(array(i)) .gt. maximum) maximum=abs(array(i))
      enddo
      
      ! divide by maximum
      if( maximum .lt. 1.0e-9) then
        print*,'array maximum too small, could not normalize!',maximum
      else
        array(:)=array(:)/maximum
      endif

      end
      
      
!-----------------------------------------------------------------------
      subroutine epicentralDistance(lat1,lon1,lat2,lon2,distance)
!-----------------------------------------------------------------------
! calculates the great circle distance between the two given
! points (assuming they lay on the unit sphere)
!
! uses haversine formula (supposed to be exact)
!
! input:
!   lon1,lat1 - longitude/latitude point 1 (in degree)
!   lon2,lat2 - for point 2
!   distance - great circle distance between point 1 and point 2
!
! return: great circle distance (in degree)
      use precisions
      implicit none
      real(WP):: lon,lonhalf,lat,lathalf,lon1,lon2,lat1,lat2,a,distance
      real(WP):: radlat1,radlat2,radlon1,radlon2
      ! convert to radian
      radlat1=DEGREE2RAD*lat1
      radlat2=DEGREE2RAD*lat2
      radlon1=DEGREE2RAD*lon1
      radlon2=DEGREE2RAD*lon2
      
      !haversine formula
      lon = radlon2 - radlon1
      lonhalf= lon*0.5
      lat = radlat2 - radlat1
      lathalf = lat*0.5
      a = sin(lathalf)*sin(lathalf)+cos(radlat1)*cos(radlat2)*sin(lonhalf)*sin(lonhalf)
     
      !rare case that a becomes too big (numerical possible)
      if( a .gt. 1.0 ) then
        print*,'greatCircleDistance: a=',a,lon1,lat1,lon2,lat2
        a=1.0
      endif      
      distance = 2.0*atan2(sqrt(a),sqrt(1.0-a)) 
      
      ! convert to degrees
      distance = RAD2DEGREE*distance
      
      ! check if Nan
      if( distance .ne. distance) then
        print*,'greatCircleDistance:',distance
        print*,'greatCircleDistance:',lon,lat,a
      endif
      return
      end

     
!-----------------------------------------------------------------------
      subroutine greatCircleDistance(vector1,vector2,distance)
!-----------------------------------------------------------------------
! calculates the great circle distance between the two given
! points (assuming they lay on the unit sphere)
!
! uses haversine formula (supposed to be exact)
!
! input:
!     vector1,vector2 - point vectors in cartesian coordinates
!
! return: great circle distance (in radian)
      use precisions
      implicit none
      real(WP):: vector1(3), vector2(3)
      real(WP):: lon,lonhalf,lat,lathalf,lon1,lon2,lat1,lat2,a
      real(WP):: distance
              
      ! get longitude and latitudes
      lon1 = atan2(vector1(2),vector1(1))
      lat1 = asin(vector1(3))
      lon2 = atan2(vector2(2),vector2(1))
      lat2 = asin(vector2(3))
      !print*,'vector1:',lat1,lon1
      !print*,'vector2:',lat2,lon2
      
      !haversine formula
      lon = lon2 - lon1
      lonhalf= lon*0.5_WP
      lat = lat2 - lat1
      lathalf = lat*0.5_WP
      a = sin(lathalf)*sin(lathalf)+cos(lat1)*cos(lat2)*sin(lonhalf)*sin(lonhalf)
     
      !rare case that a becomes too big (numerical possible)
      if( a .gt. 1.0 ) then
        !if(VERBOSE)print*,'greatCircleDistance: a=',a,dlon,dlat
        a=1.0_WP
      endif
      
      distance = 2.0_WP*atan2(sqrt(a), sqrt(1.0_WP-a)) 
      if( distance .ne. distance) then
        print*,'greatCircleDistance:',distance
        print*,'greatCircleDistance:',lon,lat,a
      endif
      return
      end
      
!-----------------------------------------------------------------------
      double precision function greatCircleDistance_alternate(vector1,vector2)
!-----------------------------------------------------------------------
! calculates the great circle distance between the two given
! points (assuming they lay on the unit sphere)
!
! input:
!     vector1,vector2 - point vectors in cartesian coordinates
!
! return: great circle distance in radian
      use precisions
      implicit none
      real(WP):: vector1(3), vector2(3)
      real(WP):: scalar, vectorlength
      integer               k
      
      if( abs(vector1(1)*vector1(1)+vector1(2)*vector1(2)+              &
     &   vector1(3)*vector1(3) - 1.0d0) .gt. 0.00001d0 .or.             &
     &   abs(vector2(1)*vector2(1)+vector2(2)*vector2(2)+               &
     &   vector2(3)*vector2(3) - 1.0d0) .gt. 0.00001d0 ) then
        print*,'vector1',(vector1(k),k=1,3),vectorlength(vector1)
        print*,'vector2',(vector2(k),k=1,3),vectorlength(vector2)
        stop 'abort-greatCircleDistance vectorlength'     
      endif
      
      call scalarproduct( vector1,vector2,scalar)
      
      greatCircleDistance_alternate = dble( acos( scalar) )
      return
      end
      
!-----------------------------------------------------------------------
      subroutine greatCircleMidpoint(vector1, vector2, midpoint)
!-----------------------------------------------------------------------
! calculates the midpoint of the great circle through vector1
! and vector2
! 
! return: midpoint vector (on the sphere)
      use precisions
      implicit none
      real(WP):: vector1(3),vector2(3),midpoint(3),normale(3)
      real(WP):: distance, vectorlength
      real(WP):: midAngle, vX,vY,vZ, cosAlpha,sinAlpha,factor
      integer:: k
      
      ! get angle for the midpoint between the two vectors
      ! (the radian of this angle for the unit sphere is equal to half
      !   of its great circle distance)
      if( abs(vectorlength(vector1) - 1.0_WP) .gt. 0.00001 .or.          &
     &   abs(vectorlength(vector2) - 1.0_WP) .gt. 0.00001 ) then
        print*,'vector1',(vector1(k),k=1,3),vectorlength(vector1)
        print*,'vector2',(vector2(k),k=1,3),vectorlength(vector2)
        stop 'abort-greatCircleMidpoint vectorlength'     
      endif

      !get the half of the distance and the multiplication factor    
      ! ( http://williams.best.vwh.net/avform.htm#Intermediate )  
      call greatCircleDistance(vector1,vector2,distance)
      midAngle = distance*0.5_WP  
      call greatCircleDistance(vector1,vector2,distance)    
      factor=sin(midAngle)/sin(distance)
      
      midpoint(1) =factor*(vector1(1)+vector2(1))
      midpoint(2) =factor*(vector1(2)+vector2(2))
      midpoint(3) =factor*(vector1(3)+vector2(3))
      
!      !debug
!      !print *,'midpoint2',factor,(midpoint(k),k=1,3)

!      !debug
!      ! geometrical calculation is to rotate a vector by half the angle
!      cosAlpha = dcos( midAngle )
!      sinAlpha = dsqrt( 1.0 - cosAlpha**2 )
!      call vectorproduct(vector1,vector2,normale)
!      call normalize(normale)
!      
!      ! rotate first vector
!      vX = (factor*(normale(1)**2)+cosAlpha)*vector1(1) +
!     & (factor*normale(1)*normale(2)-sinAlpha*normale(3))*vector1(2)+
!     & (factor*normale(1)*normale(3)+sinAlpha*normale(2))*vector1(3) 
!     
!      vY=
!     & (factor*normale(1)*normale(2)+sinAlpha*normale(3))*vector1(1)+
!     & (factor*(normale(2)**2)+cosAlpha)*vector1(2)+
!     & (factor*normale(2)*normale(3)-sinAlpha*normale(1))*vector1(3) 
!
!        vZ =
!     & (factor*normale(1)*normale(3)-sinAlpha*normale(2))*vector1(1)+
!     & (factor*normale(2)*normale(3)+sinAlpha*normale(1))*vector1(2)+
!     & (factor*(normale(3)**2)+cosAlpha)*vector1(3) 
!
!      midpoint(1)=vX
!      midpoint(2)=vY
!      midpoint(3)=vZ
!      !debug
!      print*,'mid angle',midAngle*RAD2DEGREE
!
!      !debug
!      !another possibility
!      call midpoint(vector1,vector2,midpoint)
!      call normalize(midpoint)
!
      end
      
!-----------------------------------------------------------------------
      subroutine determineTimeStep()
!-----------------------------------------------------------------------     
! determines the time step size according to approximative formula
! done by C. Tape (paragraph 5.5)
!
! input:
!     cphaseRef          - phase velocity
!     subdivisions  - grid level
!     timestepsize  - model time step
!     averageCellDistance - grid spacing
!
! returns: dt in seconds, averageCellDistance in km
      use propagationStartup; use cells
      implicit none
      real(WP),parameter:: sqrt2=1.414213562373095_WP

      ! initialize number of time steps
      numofTimeSteps = 0            
            
      ! timestep width relation (Tape, chap. 5, (5.9), p. 66) with space steps (table 4.4) [averageCellDistance in km]
      if( subdivisions .eq. 8) averageCellDistance = 17.3860147304427_WP
      if( subdivisions .eq. 7) averageCellDistance = 34.7719816287154444_WP
      if( subdivisions .eq. 6) averageCellDistance = 69.5435806022331491_WP
      if( subdivisions .eq. 5) averageCellDistance = 139.084100036946921_WP
      if( subdivisions .eq. 4) averageCellDistance = 278.143713101840262_WP
      if( subdivisions .eq. 3) averageCellDistance = 556.091606162381709_WP
      if( subdivisions .eq. 2) averageCellDistance = 1110.61907020550302_WP
      if( subdivisions .eq. 1) averageCellDistance = 2208.80170795688991_WP
      if( subdivisions .eq. 0) averageCellDistance = 4320.48077165147788_WP

      ! for routine findVertex
      distanceQuit = averageCellDistance*0.5/EARTHRADIUS
            
      ! calculate mean cphase velocity for delta phase map
      !if(DELTA ) then
      !  !reset cphase
      !  cphase = 0.0_WP
      !  !delta location source
      !  !lat=deltaLat*PI/180.0_WP
      !  !lon=deltaLon*PI/180.0_WP
      !  !vectorS(1) = cos(lat)*cos(lon)
      !  !vectorS(2) = cos(lat)*sin(lon)
      !  !vectorS(3) = sin(lat)
      !  vectorS(:)= vertices(deltaVertex,:)
      !  
      !  ! range of delta locations
      !  ! calculate mean phase velocity
      !  do n=1,numVertices
      !    vectorV(:) =  vertices(n,:)
      !    distance=EARTHRADIUS*greatCircleDistance(vectorS,vectorV)
      !    if( distance .gt. DELTARADIUS) then ! changed since DELTARADIUS was defined only after this check
      !      cphase = cphase+cphaseRef
      !    else
      !      cphase = cphase+ (cphaseRef + deltaPerturbation)
      !    endif
      !  enddo
      !  !set cphase to mean value
      !  cphase = cphase/numVertices
      !endif
            
      ! time step - trial'n error formula (5.8)              
      dt = averageCellDistance/(cphaseRef*sqrt2)        
      
      
      !if( cphaseref > 4.78 ) dt = dt/2.0
      
      end
      
!-----------------------------------------------------------------------
      function rootmeansquare(x,length)
!-----------------------------------------------------------------------     
! rms function calculates the root-mean-square value
! for the values given by vector x 
!
! returns: rms of x
        use precisions
        implicit none
        real(WP):: rootmeansquare
        integer,intent(in):: length
        real(WP),intent(in):: x(length)
        real(WP):: sum
        integer:: i
        ! see: http://en.wikipedia.org/wiki/Root_mean_square
        sum = 0.0
        do i=1,length
          sum = sum + x(i)*x(i)
        enddo        
        rootmeansquare = sqrt(sum/length)        
      end function

