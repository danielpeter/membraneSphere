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
  integer,intent(in):: vertex
  real(WP),intent(out):: area
  ! local parameters
  real(WP):: vectorA(3),vectorB(3),vectorC(3)
  integer:: n,corner1,corner2,corner3
  real(WP):: a,b,c,s,excess,tanexcess4,haversine

  !check index
  if (cellFace(vertex,0) > 6 .or. cellFace(vertex,0) < 5) then
    print *,'Error: vertex',vertex,'faces',cellFace(vertex,0)
    call stopProgram( 'abort-cellArea cellFace    ')
  endif

  corner1 = vertex
  vectorA(:) = vertices(corner1,:)

  area = 0.0
  do n = 1, cellFace(vertex,0)
    !determine triangles corner indices
    corner2=cellFace(vertex,n)
    corner3 = cellFace(vertex,(mod(n,cellFace(vertex,0))+1) )

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


    area = area+EARTHRADIUS_SQUARED*excess
  enddo

  if (area < 1.0) then
    print *,'cellArea calculation has invalid area'
    print *,'  vertex:',vertex
    print *,'  infos:',n,a,b,c,s,excess,tanexcess4
  endif

  end subroutine


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
  integer,intent(in):: vertex
  real(WP),intent(out):: distances(0:6)
  ! local parameters
  integer:: n
  real(WP):: distance
  real(WP):: vectorA(3),vectorB(3)

  ! cell center reference
  vectorA(:) = vertices(vertex, : )

  ! check index
  if (cellNeighbors(vertex,0) > 6 .or. cellNeighbors(vertex,0) < 5) then
    print *,'Error: vertex',vertex
    call stopProgram( 'abort-cellDistances cellNeighbors    ')
  endif

  distances(:) = 0.0_WP

  ! get distances
  distance = 0.0
  distances(0) = cellNeighbors(vertex,0)
  do n = 1, cellNeighbors(vertex,0)
    ! determine distance to neighbor
    vectorB(:) = vertices(cellNeighbors(vertex,n), : )

    ! calculate arc distance between these two
    call greatCircleDistance(vectorA,vectorB,distance)
    distances(n) = EARTHRADIUS*distance
  enddo

  end subroutine


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
  integer,intent(in):: vertex
  real(WP),intent(out):: lengths(0:6)
  ! local parameters
  integer:: n,corner1,corner2
  real(WP):: vectorA(3),vectorB(3),distance

  ! check index
  if (cellFace(vertex,0) > 6 .or. cellFace(vertex,0) < 5) then
    print *,'Error: vertex',vertex,'faces',cellFace(vertex,0)
    call stopProgram( 'abort-cellEdgesLength cellFace   ')
  endif

  lengths(:) = 0.0_WP

  ! get cell edge lengths
  !write(11,*) '#',vertex
  lengths(0) = cellFace(vertex,0)
  do n = 1, cellFace(vertex,0)
    ! determine which corners (they are ordered clockwise)
    corner1 = cellFace(vertex,n)
    corner2 = cellFace(vertex, mod(n,cellFace(vertex,0))+1)

    !get corresponding corner vectors
    !do k=1,3
    !  vectorA(k) = cellCorners(corner1,k)
    !  vectorB(k) = cellCorners(corner2,k)
    !enddo

    vectorA(:) = cellCorners(corner1,:)
    vectorB(:) = cellCorners(corner2,:)

    !calculate arc distance between these two
    call greatCircleDistance(vectorA,vectorB,distance)
    lengths(n) = EARTHRADIUS*distance

    !print an output file for gnuplot
    !write(11,'(3f12.8,f12)') vectorA(1),vectorA(2),vectorA(3),
     !&                             lengths(n)
  enddo

  !end output with first vertex again
  !write(11,'(3f12.8,f12)') (vectorB(k),k=1,3), lengths(1)
  !write(11,*) '\n'

  end subroutine


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
  integer,intent(in):: vertex
  real(WP),intent(out):: area
  ! local parameters
  integer:: n,corner1,corner2,corner3
  real(WP):: vectorA(3),vectorB(3),vectorC(3),angleA,angleB,angleC

  !check index
  if (cellFace(vertex,0) > 6 .or. cellFace(vertex,0) < 5) then
    print *,'Error: vertex',vertex,'faces',cellFace(vertex,0)
    call stopProgram( 'abort-cellArea cellFace    ')
  endif

  corner1 = vertex
  vectorA(:) = vertices(corner1,:)

  area = 0.0_WP
  do n = 1, cellFace(vertex,0)
    !determine triangles corner indices
    corner2=cellFace(vertex,n)
    corner3 = cellFace(vertex,(mod(n,cellFace(vertex,0))+1) )

    !get corner vectors
    vectorB(:) = cellCorners(corner2,:)
    vectorC(:) = cellCorners(corner3,:)

    !calculate corner angles of spherical triangle
    call sphericalAngle(vectorA,vectorB,vectorC,angleA)
    call sphericalAngle(vectorB,vectorA,vectorC,angleB)
    call sphericalAngle(vectorC,vectorA,vectorB,angleC)
    !debug
    !print *,'angleA',angleA/PI*180
    !print *,'angleB',angleB/PI*180
    !print *,'angleC',angleC/PI*180

    area = area+EARTHRADIUS_SQUARED*(angleA + angleB + angleC - PI)
  enddo

  end subroutine


!-----------------------------------------------------------------------
  subroutine get_dataVariation(timestep,statfile_output)
!-----------------------------------------------------------------------
! calculates heikes&randall ratio for midpoints
! and writes it to files
  implicit none
  integer,intent(in):: timestep
  logical,intent(in):: statfile_output
  ! local parameters
  character(len=24):: charstep

  ! writes precalculated values to files
  if (statfile_output) then
    if (timestep < 0) then
      charstep=""
    else
      write(charstep,*) timestep
    endif

    !for print out of cells for gnuplot
    open(11,file='OUTPUT/cellAreas.'//trim(adjustl(charstep))//'.dat')
    open(12,file='OUTPUT/cellFractionAverage.'//trim(adjustl(charstep))//'.dat')
    open(13,file='OUTPUT/cellEdgesLength.'//trim(adjustl(charstep))//'.dat')
    open(14,file='OUTPUT/cellCenterDistances.'//trim(adjustl(charstep))//'.dat')
    open(15,file='OUTPUT/cellFractions.'//trim(adjustl(charstep))//'.dat')
  endif

  !calculations
  call dataVariation(statfile_output,.false.)

  if (statfile_output) then
    !end print out
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    print *,'  files written:'
    print *,'    OUTPUT/cellAreas.'//trim(adjustl(charstep))//'.dat'
    print *,'    OUTPUT/cellFractionAverage.'//trim(adjustl(charstep))//'.dat'
    print *,'    OUTPUT/cellEdgesLength.'//trim(adjustl(charstep))//'.dat'
    print *,'    OUTPUT/cellCenterDistances.'//trim(adjustl(charstep))//'.dat'
    print *,'    OUTPUT/cellFractions.'//trim(adjustl(charstep))//'.dat'
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine dataVariation(statfile_output,plotfile_output)
!-----------------------------------------------------------------------
! calculates variation of cell areas, cell center distances,
! and cell edge lengths
!
! return: prints output
  use cells
  implicit none
  logical,intent(in):: statfile_output,plotfile_output
  ! local parameters
  real(WP):: centerDistances(0:6),edgesLength(0:6)
  real(WP):: area
  integer::  n,k
  real(WP):: averageArea, averageLength
  real(WP):: averageDist
  integer::  areaCount, lengthCount, distCount
  real(WP):: areaMin, areaMax
  real(WP):: distMin, distMax
  real(WP):: edgeMin, edgeMax
  real(WP):: fractionMin,fractionMax,fractions(0:6),fractionavg
  integer::  areaMaxIndex,areaMinIndex
  integer::  edgeMaxIndex,edgeMinIndex
  integer::  distMaxIndex,distMinIndex
  integer::  fractionMaxIndex,fractionMinIndex
  real(WP):: fractionR

  !initialize
  averageArea   = 0.0_WP
  averageLength = 0.0_WP
  averageDist   = 0.0_WP
  areaCount     = 0
  lengthCount   = 0
  distCount     = 0
  areaMin       = 1000000000.0_WP
  areaMax       = 0.0_WP
  distMin       = 1000000000.0_WP
  distMax       = 0.0_WP
  edgeMin       = 1000000000.0_WP
  edgeMax       = 0.0_WP
  fractionMin   = 1000000000.0_WP
  fractionMax   = 0.0_WP
  fractionR     = 0.0_WP

  !for all vertices (cell centers)
  print *,'calculating averages for ',numVertices,'vertices...'
  do n = 1, numVertices
    !get spherical area of vertex's cell
    call getCellArea(n,area)

    ! output to cellAreas
    if (statfile_output) write(11,'(f18.12)') area

    !print *,'area',area
    if (area > areaMax) then
       areaMax = area
       areaMaxIndex = n
    endif
    if (area < areaMin) then
      areaMin = area
      areaMinIndex = n
    endif

    averageArea = averageArea + area
    areaCount = areaCount + 1

    !get arc distances from center to neighbors
    call getcellCenterDistances(n,centerDistances)

    ! output to cellCenterDistances
    if (statfile_output) write(14,'(7f18.12)') (centerDistances(k),k=0,6)

    do k = 1, int(centerDistances(0))
      !debug
      !print *,'distance one',centerDistances(k)
      if (centerDistances(k) > distMax) then
        distMax = centerDistances(k)
        distMaxIndex = n
      endif
      if (centerDistances(k) < distMin) then
         distMin = centerDistances(k)
         distMinIndex = n
      endif

      averageDist = averageDist + centerDistances(k)
      distCount = distCount + 1
    enddo

    !get cell edges length
    call getCellEdgesLength(n,edgesLength)

    ! output to cellEdgesLength
    if (statfile_output) write(13,'(7f18.12)') (edgesLength(k),k=0,6)

    do k = 1, int(edgesLength(0))
      if (edgesLength(k) > edgeMax) then
        edgeMax = edgesLength(k)
        edgeMaxIndex = n
      endif
      if (edgesLength(k) < edgeMin) then
        edgeMin =edgesLength(k)
        edgeMinIndex = n
      endif

      averageLength = averageLength+edgesLength(k)
      lengthCount = lengthCount +1
    enddo

    !get cell derivative fractions
    call cellDerivativeFractions(n,fractions,fractionavg,plotfile_output)

    ! output to cellFractionAverage.dat
    if (statfile_output) write(12,*) fractionavg
    ! output to cellFractions.dat
    if (statfile_output) write(15,'(f18.12,6e18.8)') (fractions(k),k=0,6)

    do k = 1, int(fractions(0))
      ! global function (see Miura, 2005: eq (7), p. 2822)
      fractionR = fractionR + fractions(k)**4

      ! min / max
      if (fractions(k) > fractionMax) then
        fractionMax =fractions(k)
        fractionMaxIndex = n
      endif
      if (fractions(k) < fractionMin) then
        fractionMin =fractions(k)
        fractionMinIndex = n
      endif

    enddo
  enddo !n

  !info
  print *
  print *,'AREAS'
  print *,'number of face areas: ',areaCount
  print *,'average area        : ',averageArea/areaCount, ' km2'
  print *,'total areas         : ',averageArea, ' km2'
  print *,'            expected: ',4*PI*EARTHRADIUS_SQUARED, ' km2'
  print *,'a_min/a_max         : ', areaMin/areaMax
  print *,'minimum area        : ',areaMin, ' km2 ',areaMinIndex,'(min index)'
  print *,'maximum area        : ',areaMax, ' km2 ',areaMaxIndex,'(max index)'
  print *
  print *,'CELL CENTER DISTANCES'
  print *,'number of distances : ',distCount
  print *,'average distances   : ', averageDist/distCount, ' km'
  print *,'d_min/d_max         : ', distMin/distMax
  print *,'minimum distance    : ',distMin, ' km ',distMinIndex,'(min index)'
  print *,'maximum distance    : ',distMax, ' km ',distMaxIndex,'(max index)'
  print *
  print *,'CELL EDGES'
  print *,'number of edges     : ', lengthCount
  print *,'average lengths     : ',averageLength/lengthCount, ' km'
  print *,'e_min/e_max         : ', edgeMin/edgeMax
  print *,'minimum edge length : ',edgeMin, ' km ',edgeMinIndex,'(min index)', &
                ' (radius fraction',edgeMin/EarthRadius,')'
  print *,'maximum edge length : ',edgeMax, ' km ',edgeMaxIndex,'(max index)', &
                ' (radius fraction',edgeMax/EarthRadius,')'
  print *
  print *,'CELL DERIVATIVE FRACTIONS'
  print *,'min/max fraction    : ',fractionMin,fractionMax,' - ',fractionMinIndex,'(min index)',fractionMaxIndex,'(max index)'
  print *,'global function R   : ',fractionR
  print *

  end subroutine


!-----------------------------------------------------------------------
  subroutine cellDerivativeFractions(vertex,fractions,fractionavg,plotfile_output)
!-----------------------------------------------------------------------
! calculates for vertex the fraction of the distance
! between an edge midpoint and the corresponding
! cell center line midpoint for all cell neighbors
! (defined after Heikes&Randall
!
! input:
!    vertex     - index of reference vertex
!    fractions      - fractions for all neighbors
!
! return: fractions
  use cells
  implicit none
  integer,intent(in):: vertex
  real(WP),intent(out):: fractions(0:6),fractionavg
  logical,intent(in):: plotfile_output
  ! local parameters
  real(WP):: fraction, distance
  real(WP):: vectorA(3),vectorB(3),vectorC(3),vectorD(3)
  integer:: n,k,corner1,corner2
  !real(WP) greatCircleDistance
  real(WP):: midpointCenterLine(3)
  real(WP):: midpointEdge(3)
  real(WP):: edgeLength

  !check index
  if (cellFace(vertex,0) > 6.or.cellFace(vertex,0) < 5) then
    print *,'vertex',vertex
    stop 'Abort - cellDerivativeFraction cellFace'
  endif
  if (cellNeighbors(vertex,0) > 6.or.cellNeighbors(vertex,0) < 5) then
    print *,'vertex',vertex
    stop 'Abort - cellDerivativeFraction cellNeighbors'
  endif

  !set how many neighbors
  fractions(:) = 0.0
  fractions(0) = cellNeighbors(vertex,0)

  !write(21,*) '\n\n#',vertex
  !write(22,*) '\n\n#',vertex
  !write(23,*) '\n\n#',vertex

  !for each neighbor cell
  do n = 1, cellNeighbors(vertex,0)
    !get center vectors
    vectorA(:) = vertices(vertex,:)
    vectorB(:) = vertices(cellNeighbors(vertex,n),:)
    !debug
    !print *,(vectorA(k),k=1,3)
    !print *,(vectorB(k),k=1,3)

    !get cell center line midpoint
    call greatCircleMidpoint(vectorA,vectorB,midpointCenterLine)

    !get midpoint of corresponding cell edge
    !first determine which corners (they are ordered clockwise)
    corner1=cellFace(vertex,n)
    corner2 = cellFace(vertex, mod(n,int(cellFace(vertex,0)))+1)

    !get corresponding corner vectors
    vectorC(:) = cellCorners(corner1,:)
    vectorD(:) = cellCorners(corner2,:)
    !debug
    !print *,(vectorC(k),k=1,3)
    !print *,(vectorD(k),k=1,3)

    !get cell edge midpoint and length
    call greatCircleMidpoint(vectorC,vectorD,midpointEdge)
    call greatCircleDistance(vectorC,vectorD,edgeLength)

    !get distance (in rad) between two midpoints
    call greatCircleDistance(midpointCenterLine,midpointEdge,distance)

    !fraction
    fraction = distance/edgeLength
    fractions(n) = fraction

    !debug
    !print *,n,' fraction:',fraction
    !print *,'midpoint vectors:'
    !print *,(midpointCenterLine(k),k=1,3)
    !print *,(midpointEdge(k), k=1,3)
    !print *,'distance:',distance
    !print *,'edge ',n,' Length:',edgeLength

    if (plotfile_output .eqv. .true.) then
      ! output to cellCenterLineMidpoints.dat
      write(21,'(3f12.8)')(midpointCenterLine(k),k=1,3)
      ! output to cellEdgeMidpoints.dat
      write(22,'(3f12.8)')(midpointEdge(k),k=1,3)
      ! output to cellCenterLineMidpoints.mat
      write(31,'(3f12.8)')(midpointCenterLine(k),k=1,3)
      ! output to cellEdgeMidpoints.mat
      write(32,'(3f12.8)')(midpointEdge(k),k=1,3)
      ! output to cellCenterLines.dat
      write(23,'(4f12.8)') (vectorA(k),k=1,3), fraction
      write(23,'(4f12.8)') (vectorB(k),k=1,3), fraction
      ! output to cellCenterLines.mat
      write(33,'(4f12.8)') (vectorA(k),k=1,3), fraction
      write(33,'(4f12.8)') (vectorB(k),k=1,3), fraction
    endif
  enddo !n

  !for face average
  fractionavg = 0.0
  do n = 1, cellNeighbors(vertex,0)
    fractionavg = fractionavg + fractions(n)
  enddo
  fractionavg = fractionavg/cellNeighbors(vertex,0)

  end subroutine

