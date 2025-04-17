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

! spherical grid construction
!
! dodecahedron - icosahedron grid construction algorithm.
! used in generating the triangular and hexagonal grids.
!
! - Based on original dod5.f and dodecv.h from John Woodhouse and Carl Tape
!
! - Begun in Nov 2001 (John's program); 14-Feb-2003.
!
! - Carl Tape:
!   17-Feb-2003, 18-Dec-2001
!   John's program, first modified on 28-Nov-2001.
!   Original file is dodecvct.f; this file is dod5.f.
!
!   This file generate the grid-points for the geodesic grids,
!   i.e. the triangular grids and the hexagonal grids.
!   It also determines the nearest neighbors, which is needed
!   when using a FDM approximation for the Laplacian on the sphere.
!   Data are written to five files.
!
!   Note: include file is dodecv.h
!   This program is now in double precision (17-Feb-2003).
!   This program is self-contained and does not require a Makefile.
!
! - Daniel Peter:
!   updated for newer Fortran syntax, removing dependency on .h include file, using dynamic array allocations
!   for scaling to larger meshes

!-----------------------------------------------------------------------
module dodeco
!-----------------------------------------------------------------------

  implicit none

  ! dodecahedron-icosahedron solid

  ! spherical grid parameters
  integer :: MAXLEVELS,MAXTRIANGLES,MAXVERTICES

  double precision, dimension(:,:), allocatable :: tvertices
  double precision, dimension(:,:), allocatable :: voronoiVertices

  integer, dimension(:,:), allocatable :: indexTriangleVert
  integer, dimension(:,:), allocatable :: indexTriContainVertex

  integer, dimension(:), allocatable :: ntcv,nvnear
  integer, dimension(:,:), allocatable :: indexvnear

  integer, dimension(:), allocatable :: numTriangles,numVertices
  integer, dimension(:), allocatable :: isons,nTrianglesContainVertex

  ! number of initial solid vertices/triangles
  integer :: nv, nt

  ! hash table search in subdivision for adding midpoint vertices
  ! (much, much faster than brute-force linear search)
  logical, parameter :: USE_HASH_TABLE_SEARCH = .true.

end module dodeco


!-----------------------------------------------------------------------
  program gridConstruction
!-----------------------------------------------------------------------
  use dodeco
  implicit none
  integer :: levels, level, nVertices, nTriangles

  print *,'spherical grid construction'
  print *,'---------------------------'

  ! calculates up to this order grids (max is 8?)
  print *,'number of grid subdivisions?'
  read(*,*) levels

  print *
  print *,'levels up to: ',levels
  print *

  ! setup arrays
  call allocateArrays(levels)

  ! create dodecahedron-icosahedron solid
  call dodec()

  ! subdivide each triangle by the given number of levels
  call subdivide(levels)

  ! info
  print *,'grid data:'
  ! print out numbers of vertices and triangles for each level
  do level = 0, levels
    call getDodProperties(level, nVertices, nTriangles)
    print *,'  properties - level ',level,':',nVertices,' vertices',nTriangles,' triangles'
  enddo

  end program

!-----------------------------------------------------------------------
  subroutine allocateArrays(levels)
!-----------------------------------------------------------------------
  use dodeco
  implicit none
  integer, intent(in) :: levels
  ! local parameters
  integer :: ier

  ! maximum values
  MAXLEVELS    = levels
  MAXTRIANGLES = 20 * (4**(MAXLEVELS + 1) - 1)
  MAXVERTICES  = 30 * (4**MAXLEVELS) + 2

  print *,'allocating arrays'
  print *,'  maximum triangles: ',MAXTRIANGLES
  print *,'  maximum vertices : ',MAXVERTICES

  ! array allocations
  allocate(tvertices(3, MAXVERTICES), &
           voronoiVertices(3, MAXTRIANGLES),stat=ier)
  tvertices(:,:) = 0.d0
  voronoiVertices(:,:) = 0.d0

  allocate(indexTriangleVert(3, MAXTRIANGLES), &
           indexTriContainVertex(6, MAXVERTICES), &
           ntcv(MAXVERTICES), &
           nvnear(MAXVERTICES), &
           indexvnear(6, MAXVERTICES),stat=ier)
  indexTriangleVert(:,:) = 0
  indexTriContainVertex(:,:) = 0
  ntcv(:) = 0
  nvnear(:) = 0
  indexvnear(:,:) = 0

  allocate(numTriangles(0:MAXLEVELS+1), &
           numVertices(0:MAXLEVELS), &
           isons(MAXTRIANGLES), &
           nTrianglesContainVertex(MAXVERTICES),stat=ier)
  numTriangles(:) = 0
  numVertices(:) = 0
  isons(:) = 0
  nTrianglesContainVertex(:) = 0

  print *,'  allocation done'

  end subroutine

!-----------------------------------------------------------------------
  subroutine dodec()
!-----------------------------------------------------------------------
! generates the x,y,z coordinates of the 20 vertices of a dodecohedron.
!
! In the calling this subroutine, x,y,z should be dimensioned 20.
! indexpoly(i,j) gives the index to the i-th vertex of the j-th
! pentagonal face, traced in a clockwise sense, looking from the
! exterior.
!
! indexpoly(6,j) is equal to indexpoly(1,j), completing the pentagon
! by repeating the first vertex.
!
! xcNormal(j), ycNormal(j), zcNormal(j) are the components of the
! faces' unit normals.

  use dodeco
  implicit none
  ! local parameters
  double precision :: a72,a36,c72,s72,c36,s36
  double precision :: xsid2,xsid,xv,zv
  double precision :: a,vx,vy,tx,ty,tz,ta,ts,s,ss
  double precision :: alph,ph,th
  integer :: i,j,i2,i3,i4,ip,iv

  ! alternative views (was equivalence in f77)
  ! equivalence (x,x1), (y,y1), (z,z1)
  ! will use reshape for x1/y1/z1
  double precision, dimension(5,4) :: x, y, z
  double precision, dimension(20) :: x1, y1, z1

  ! polygon
  integer :: indexpoly(6, 12)
  double precision :: xcNormal(12), ycNormal(12), zcNormal(12)
  double precision :: vec1(3, 12)
  double precision :: vec2(2, 12)

  double precision, parameter :: pi  = 4.0*datan(1.d0)

  print *,'creating dodecahedron-icosahedron solid'

  a72 = 72.d0 * pi/180.d0    !arc 72 deg
  a36 = a72 * 0.5d0          !arc 36 deg

  c72 = dcos(a72)
  s72 = dsin(a72)
  c36 = dcos(a36)
  s36 = dsin(a36)

  xsid2 = (1.0 - c72)**2 + s72**2
  xsid = dsqrt(xsid2)
  xv = c72 * xsid2 / (1.0 - c72)
  zv = dsqrt(xsid2 - xv**2)

  ! get vertices of dodecahedron
  do i = 1,5
    a = (i - 1.d0)*a72
    x(i,1) = dcos(a)
    y(i,1) = dsin(a)
    z(i,1) = 0.d0
  enddo

  do i = 1,5
    x(i,2) = x(i,1)*(1.d0 + xv)
    y(i,2) = y(i,1)*(1.d0 + xv)
    z(i,2) = zv
  enddo

  do i = 1,5
    i2 = 1 + mod(i,5)
    i4 = 1 + mod(i + 2, 5)
    vx = 0.5d0 * (x(i2,1) + x(i,1))
    vy = 0.5d0 * (y(i2,1) + y(i,1))
    tx = 0.5d0 * (x(i2,2) + x(i,2)) - vx
    ty = 0.5d0 * (y(i2,2) + y(i,2)) - vy
    tz = 0.5d0 * (z(i2,2) + z(i,2))
    ta = dsqrt(tx**2 + ty**2 + tz**2)
    tx = tx/ta
    ty = ty/ta
    tz = tz/ta
    ts = dsqrt((x(i4,1)-vx)**2 + (y(i4,1)-vy)**2)
    x(i,3) = vx + ts*tx
    y(i,3) = vy + ts*ty
    z(i,3) = ts*tz
  enddo

  do i = 1,5
    i2 = 1 + mod(i,5)
    i3 = 1 + mod(i+1,5)
    i4 = 1 + mod(i+2,5)
    alph = c72 * xsid2 /((x(i,1) + x(i2,1))*(x(i4,1) - x(i3,1)) &
                        +(y(i,1) + y(i2,1))*(y(i4,1) - y(i3,1)))
    tx = alph * (x(i,1) + x(i2,1))
    ty = alph * (y(i,1) + y(i2,1))
    tz = sqrt(xsid2 - tx**2 - ty**2)
    x(i,4) = x(i,3) + tx
    y(i,4) = y(i,3) + ty
    z(i,4) = z(i,3) + tz
  enddo

  tx = 0.d0
  ty = 0.d0
  tz = 0.d0
  do i = 1,5
    do j = 1,4
      tx = tx + x(i,j)
      ty = ty + y(i,j)
      tz = tz + z(i,j)
    enddo
  enddo

  tx = tx/20.d0
  ty = ty/20.d0
  tz = tz/20.d0

  do i = 1,5
    do j = 1,4
      x(i,j) = x(i,j) - tx
      y(i,j) = y(i,j) - ty
      z(i,j) = z(i,j) - tz
      ! normalize vector x,y,z
      s = dsqrt(x(i,j)**2 + y(i,j)**2 + z(i,j)**2)
      x(i,j) = x(i,j)/s
      y(i,j) = y(i,j)/s
      z(i,j) = z(i,j)/s
    enddo
  enddo

  ! update module arrays x1, y1, z1 with shape [20]
  x1 = reshape(x, [20])
  y1 = reshape(y, [20])
  z1 = reshape(z, [20])

  print *,'  vertex coordinates done'

  ! set indices of each dodecahedron face
  indexpoly(1,1) = 1
  indexpoly(2,1) = 2
  indexpoly(3,1) = 3
  indexpoly(4,1) = 4
  indexpoly(5,1) = 5
  indexpoly(6,1) = 1

  indexpoly(1,2) = 11
  indexpoly(2,2) = 7
  indexpoly(3,2) = 2
  indexpoly(4,2) = 1
  indexpoly(5,2) = 6
  indexpoly(6,2) = 11

  indexpoly(1,3) = 12
  indexpoly(2,3) = 8
  indexpoly(3,3) = 3
  indexpoly(4,3) = 2
  indexpoly(5,3) = 7
  indexpoly(6,3) = 12

  indexpoly(1,4) = 13
  indexpoly(2,4) = 9
  indexpoly(3,4) = 4
  indexpoly(4,4) = 3
  indexpoly(5,4) = 8
  indexpoly(6,4) = 13

  indexpoly(1,5) = 14
  indexpoly(2,5) = 10
  indexpoly(3,5) = 5
  indexpoly(4,5) = 4
  indexpoly(5,5) = 9
  indexpoly(6,5) = 14

  indexpoly(1,6) = 15
  indexpoly(2,6) = 6
  indexpoly(3,6) = 1
  indexpoly(4,6) = 5
  indexpoly(5,6) = 10
  indexpoly(6,6) = 15

  indexpoly(1,7) = 6
  indexpoly(2,7) = 15
  indexpoly(3,7) = 20
  indexpoly(4,7) = 16
  indexpoly(5,7) = 11
  indexpoly(6,7) = 6

  indexpoly(1,8) = 7
  indexpoly(2,8) = 11
  indexpoly(3,8) = 16
  indexpoly(4,8) = 17
  indexpoly(5,8) = 12
  indexpoly(6,8) = 7

  indexpoly(1,9) = 8
  indexpoly(2,9) = 12
  indexpoly(3,9) = 17
  indexpoly(4,9) = 18
  indexpoly(5,9) = 13
  indexpoly(6,9) = 8

  indexpoly(1,10) = 9
  indexpoly(2,10) = 13
  indexpoly(3,10) = 18
  indexpoly(4,10) = 19
  indexpoly(5,10) = 14
  indexpoly(6,10) = 9

  indexpoly(1,11) = 10
  indexpoly(2,11) = 14
  indexpoly(3,11) = 19
  indexpoly(4,11) = 20
  indexpoly(5,11) = 15
  indexpoly(6,11) = 10

  indexpoly(1,12) = 18
  indexpoly(2,12) = 17
  indexpoly(3,12) = 16
  indexpoly(4,12) = 20
  indexpoly(5,12) = 19
  indexpoly(6,12) = 18


  ! determine for each dodecahedron face its normal which
  ! projected to the sphere is the new midpoint of this face
  do ip = 1,12
    xcNormal(ip) = 0.d0
    ycNormal(ip) = 0.d0
    zcNormal(ip) = 0.d0
    do iv = 1,5
      xcNormal(ip) = xcNormal(ip) + x1(indexpoly(iv,ip))
      ycNormal(ip) = ycNormal(ip) + y1(indexpoly(iv,ip))
      zcNormal(ip) = zcNormal(ip) + z1(indexpoly(iv,ip))
    enddo

    ! normalize vector xcNormal,ycNormal,zcNormal
    ss = dsqrt(xcNormal(ip)**2 + ycNormal(ip)**2 + zcNormal(ip)**2)
    xcNormal(ip) = xcNormal(ip)/ss
    ycNormal(ip) = ycNormal(ip)/ss
    zcNormal(ip) = zcNormal(ip)/ss

    ! determine Cartesian coordinates of each midpont of a face
    if (ip == 1) then
      ! the first must lie on this position
      vec1(1,ip) = 1.d0
      vec1(2,ip) = 0.d0
      vec1(3,ip) = 0.d0
      vec2(1,ip) = 0.d0
      vec2(2,ip) = 1.d0
    else if (ip >= 2 .and. ip <= 6) then
      ph = datan2(ycNormal(ip), xcNormal(ip))
      th = dacos(zcNormal(ip))
      vec1(1,ip) = - dcos(th)*dcos(ph)
      vec1(2,ip) = - dcos(th)*dsin(ph)
      vec1(3,ip) = dsin(th)
      vec2(1,ip) = - dsin(ph)
      vec2(2,ip) = dcos(ph)
    else if (ip >= 7.and. ip <= 11) then
      ph = datan2(ycNormal(ip), xcNormal(ip))
      th = dacos(zcNormal(ip))
      vec1(1,ip) = dcos(th)*dcos(ph)
      vec1(2,ip) = dcos(th)*dsin(ph)
      vec1(3,ip) = -dsin(th)
      vec2(1,ip) = dsin(ph)
      vec2(2,ip) = -dcos(ph)
    else if (ip == 12) then
      !the last must lie on opposite position to first one
      vec1(1,ip) = -1.d0
      vec1(2,ip) = 0.d0
      vec1(3,ip) = 0.d0
      vec2(1,ip) = 0.d0
      vec2(2,ip) = 1.d0
    else
      stop 'dodec: poly error'
    endif
  enddo

  print *,'  normals done'

  ! fill common vertices array with all vertex vectors
  ! and the new midpoints of each polygon face
  nv = 0
  ! fill in vertex vectors
  do i = 1,20
    nv = nv + 1
    tvertices(1,nv) = x1(i)
    tvertices(2,nv) = y1(i)
    tvertices(3,nv) = z1(i)
  enddo
  ! fill in normal vectors
  do i = 1,12
    nv = nv + 1
    tvertices(1,nv) = xcNormal(i)
    tvertices(2,nv) = ycNormal(i)
    tvertices(3,nv) = zcNormal(i)
  enddo

  print *,'  common vertices done'

  ! fill common index array for all triangles and count
  ! how many there are
  nt = 0
  do i = 1,12
    do j = 1,5
      nt = nt + 1
      indexTriangleVert(1,nt) = 20+i
      indexTriangleVert(2,nt) = indexpoly(j,i)
      indexTriangleVert(3,nt) = indexpoly(j+1,i)
    enddo
  enddo

  print *,'  solid done'

  end subroutine

!-----------------------------------------------------------------------
  subroutine subdivide(levels)
!-----------------------------------------------------------------------
  use dodeco
  implicit none
  integer, intent(in) :: levels

  ! local parameters
  double precision :: vvMidpoint(3,3),vvec(3),azms(6)
  double precision :: ss,svvec
  double precision, external :: dazimv,ddot3ScalarProduct
  integer :: indv(3),idx(6),itmp(6)
  integer :: i,k,is,is1,ic,iv,it,il,iv1,iv2,iv3
  integer :: ii,jj,kk,ii1,it0,it1,ivgot,itc,nlevel
  integer :: index1,index2,indexTmp
  logical :: is_done
  ! hash table
  type :: edge_entry
    integer :: v1, v2       ! The two vertex indices
    integer :: mid_id       ! ID of the midpoint vertex
    type(edge_entry), pointer :: next => null()  ! For collision handling
  end type edge_entry
  type(edge_entry), pointer :: hash_table(:)
  integer :: hash_size,hash_key
  ! hash table search
  double precision :: vvMid(3)
  ! timing
  double precision :: start_time,end_time
  ! tolerance (machine precision for double precision)
  double precision, parameter :: TOL_EPS = 1.0d-12

  print *,'subdividing grid'

  nlevel = 0
  numTriangles(nlevel) = 0

  do while (nlevel <= levels)

    numVertices(nlevel) = nv
    numTriangles(nlevel + 1) = nt

    if (nlevel >= levels) exit

    print *,'  proceeding level ', nlevel, ' - triangles: ',numTriangles(nlevel+1) - numTriangles(nlevel)

    ! Start timing for this level
    call cpu_time(start_time)

    ! way 1 - linear search
    ! becomes very slow for levels > 6 ..
    if (.not. USE_HASH_TABLE_SEARCH) then
      do it = numTriangles(nlevel) + 1, numTriangles(nlevel + 1)
        ! get midpoint of each triangle edge
        do is = 1,3
          is1 = mod(is,3) + 1   ! is == 1 - > is1 == 2 / is == 2 - > is1 == 3 / is == 3 - > is1 == 1
          ! edge indices
          index1 = indexTriangleVert(is,it)
          index2 = indexTriangleVert(is1,it)
          do k = 1,3
            vvMidpoint(k,is) = 0.5d0 * (tvertices(k,index1) + tvertices(k,index2))
          enddo
          ! normalize vector vvMidpoint
          ss = vvMidpoint(1,is)**2 + vvMidpoint(2,is)**2 + vvMidpoint(3,is)**2
          ss = 1.d0/dsqrt(ss)
          do k = 1,3
            vvMidpoint(k,is) = vvMidpoint(k,is)*ss
          enddo
        enddo

        ! add new midpoints to common vertices array
        do is = 1,3

          is_done = .false.

          ! this search goes through all previous vertices and scales like O(n^2).
          ! TODO: to speed up this search, one could try to implement a faster search instead... try hash table search instead
          ! brute-force search
          do iv = 1,nv
            ! check if we don't have it yet
            if (abs(vvMidpoint(1,is) - tvertices(1,iv)) < TOL_EPS) then
              if (abs(vvMidpoint(2,is) - tvertices(2,iv)) < TOL_EPS) then
                if (abs(vvMidpoint(3,is) - tvertices(3,iv)) < TOL_EPS) then
                  indv(is) = iv
                  is_done = .true.
                  exit  ! Early exit from loop once found
                endif
              endif
            endif
          enddo

          ! add vertex
          if (.not. is_done) then
            nv = nv + 1
            do k = 1,3
              tvertices(k,nv) = vvMidpoint(k,is)
            enddo
            indv(is) = nv
          endif
        enddo

        ! add four new sub-triangles to common index array of triangles
        isons(it) = nt

        ! Add first triangle
        nt = nt + 1
        indexTriangleVert(1,nt) = indexTriangleVert(1,it)
        indexTriangleVert(2,nt) = indv(1)
        indexTriangleVert(3,nt) = indv(3)

        ! Add second triangle
        nt = nt + 1
        indexTriangleVert(1,nt) = indexTriangleVert(2,it)
        indexTriangleVert(2,nt) = indv(2)
        indexTriangleVert(3,nt) = indv(1)

        ! Add third triangle
        nt = nt + 1
        indexTriangleVert(1,nt) = indexTriangleVert(3,it)
        indexTriangleVert(2,nt) = indv(3)
        indexTriangleVert(3,nt) = indv(2)

        ! Add fourth triangle
        nt = nt + 1
        indexTriangleVert(1,nt) = indv(2)
        indexTriangleVert(2,nt) = indv(3)
        indexTriangleVert(3,nt) = indv(1)
      enddo

    else
      ! way 2 - hash table search
      ! hash table search scales like ~O(n) instead of ~O(n**2) for linear search

      ! create hash table
      call create_hash_table()

      ! loop over triangles
      do it = numTriangles(nlevel) + 1, numTriangles(nlevel + 1)
        ! get midpoint of each triangle edge
        do is = 1,3
          is1 = mod(is,3) + 1   ! is == 1 - > is1 == 2 / is == 2 - > is1 == 3 / is == 3 - > is1 == 1
          ! triangle edge
          index1 = indexTriangleVert(is,it)
          index2 = indexTriangleVert(is1,it)

          ! makes sure index1 < index2
          if (index1 > index2) then
            indexTmp = index1
            index1 = index2
            index2 = indexTmp
          endif

          ! create unique hash key
          hash_key = create_hash_key(index1,index2)

          ! check hash table entry
          iv = check_hash_table(index1,index2,hash_key)

          if (iv == 0) then
            ! creates midpoint
            vvMid(:) = 0.5d0 * (tvertices(:,index1) + tvertices(:,index2))
            ! normalize vector vvMidpoint
            ss = vvMid(1)*vvMid(1) + vvMid(2)*vvMid(2) + vvMid(3)*vvMid(3)
            ss = 1.d0/dsqrt(ss)
            vvMid(:) = vvMid(:)*ss

            ! add new vertex
            nv = nv + 1
            indv(is) = nv
            tvertices(:,nv) = vvMid(:)

            ! add hash table entry
            call add_hash_table_entry(index1,index2,hash_key,nv)
          else
            ! previously stored
            indv(is) = iv
          endif
        enddo

        ! add four new sub-triangles to common index array of triangles
        isons(it) = nt

        ! Add first triangle
        nt = nt + 1
        indexTriangleVert(1,nt) = indexTriangleVert(1,it)
        indexTriangleVert(2,nt) = indv(1)
        indexTriangleVert(3,nt) = indv(3)

        ! Add second triangle
        nt = nt + 1
        indexTriangleVert(1,nt) = indexTriangleVert(2,it)
        indexTriangleVert(2,nt) = indv(2)
        indexTriangleVert(3,nt) = indv(1)

        ! Add third triangle
        nt = nt + 1
        indexTriangleVert(1,nt) = indexTriangleVert(3,it)
        indexTriangleVert(2,nt) = indv(3)
        indexTriangleVert(3,nt) = indv(2)

        ! Add fourth triangle
        nt = nt + 1
        indexTriangleVert(1,nt) = indv(2)
        indexTriangleVert(2,nt) = indv(3)
        indexTriangleVert(3,nt) = indv(1)
      enddo

      ! clear hash table
      call free_hash_table()

    endif

    ! End timing for this level
    call cpu_time(end_time)
    print *, '    elapsed time',sngl(end_time-start_time), 's'

    ! update level
    nlevel = nlevel + 1

  enddo

  ! finds the Voronoi vertex corresponding to each triangle
  print *,'finding voronoi vertex to each triangle'

  do il = 0,levels
    ! do for each triangle at this level
    do it = numTriangles(il) + 1, numTriangles(il+1)
      ! indices of triangle corners
      iv1 = indexTriangleVert(1,it)
      iv2 = indexTriangleVert(2,it)
      iv3 = indexTriangleVert(3,it)

      ! finds the voronoi center of the three vertices by the equation
      !vvMidpoint = v1 x v2 + v2 x v1 + v3 x v2

      vvec(1) = + tvertices(2,iv1) * tvertices(3,iv2) &
                - tvertices(3,iv1) * tvertices(2,iv2) &
                + tvertices(2,iv2) * tvertices(3,iv3) &
                - tvertices(3,iv2) * tvertices(2,iv3) &
                + tvertices(2,iv3) * tvertices(3,iv1) &
                - tvertices(3,iv3) * tvertices(2,iv1)

      vvec(2) = + tvertices(3,iv1) * tvertices(1,iv2) &
                - tvertices(1,iv1) * tvertices(3,iv2) &
                + tvertices(3,iv2) * tvertices(1,iv3) &
                - tvertices(1,iv2) * tvertices(3,iv3) &
                + tvertices(3,iv3) * tvertices(1,iv1) &
                - tvertices(1,iv3) * tvertices(3,iv1)

      vvec(3) = + tvertices(1,iv1) * tvertices(2,iv2) &
                - tvertices(2,iv1) * tvertices(1,iv2) &
                + tvertices(1,iv2) * tvertices(2,iv3) &
                - tvertices(2,iv2) * tvertices(1,iv3) &
                + tvertices(1,iv3) * tvertices(2,iv1) &
                - tvertices(2,iv3) * tvertices(1,iv1)

      ! normalize
      svvec = 1.d0 / dsqrt(ddot3ScalarProduct(vvec,vvec))
      do i = 1,3
        vvec(i) = vvec(i)*svvec  ! normalize vectors
      enddo

      ! if the vectors are arranged in the wrong direction, then reverse the order
      if (ddot3ScalarProduct(vvec, tvertices(1,iv1)) < 0.d0) then
        do i = 1,3
          vvec(i) = -vvec(i)
        enddo
      endif
      ! save with the index it corresponding to the triangle
      do i = 1,3
        voronoiVertices(i,it) = vvec(i)
      enddo
    enddo   ! it
  enddo  ! il

  ! >> start of il do loop
  do il = 0,levels
    ! vertices from 1 to the max for level il
    do iv = 1, numVertices(il)
      ! Number of Triangles Containing this
      !Vertex
      nTrianglesContainVertex(iv) = 0
    enddo  ! iv

    do it = numTriangles(il) + 1, numTriangles(il+1)
      do ic = 1,3         ! 3 vertex indices
        iv = indexTriangleVert(ic,it)
        nTrianglesContainVertex(iv) = nTrianglesContainVertex(iv) + 1

        if (nTrianglesContainVertex(iv) > 6) stop 'more than 6 triangles'

        indexTriContainVertex(nTrianglesContainVertex(iv),iv) = it
      enddo  ! ic
    enddo  ! it

    !  arrange in increasing order around each vertex
    print *,'  arrange vertex orders'

    do iv = 1, numVertices(il)
      if (nTrianglesContainVertex(iv) /= 5 .and. nTrianglesContainVertex(iv) /= 6) &
        stop 'unexpected vertex order'

      do itc = 1, nTrianglesContainVertex(iv)
        it = indexTriContainVertex(itc,iv)
        azms(itc) = dazimv(tvertices(1,iv), voronoiVertices(1,it))
      enddo  ! itc

      ! sorting routine
      call drsoinc(azms, nTrianglesContainVertex(iv), idx)

      do itc = 1, nTrianglesContainVertex(iv)
        itmp(itc) = indexTriContainVertex(itc,iv)
      enddo  ! itc
      do itc = 1, nTrianglesContainVertex(iv)
        indexTriContainVertex(itc,iv) = itmp(idx(itc))
      enddo  ! itc
      do itc = nTrianglesContainVertex(iv) + 1, 6
        indexTriContainVertex(itc,iv) = 0
      enddo  !itc
    enddo  ! iv

    !  We want to use info in ictv to figure out vertices
    !  corresponding to each voronoi face.
    print *,'  vertices to each voronoi face'

    do iv = 1, numVertices(il)
      ! voronoi vertices around this vertex are indexed by
      !   indexTriContainVertex(,iv)
      ! and these can be used to access the triangles involved
      do ii = 1,nTrianglesContainVertex(iv)
        ii1 = 1 + mod(ii,nTrianglesContainVertex(iv))
        it0 = indexTriContainVertex(ii,iv)
        it1 = indexTriContainVertex(ii1,iv)

        ! find which point these triangles have in common, other than iv
        ivgot = 0
        do jj = 1,3
          do kk = 1,3
            if (indexTriangleVert(jj,it0) == indexTriangleVert(kk,it1) .and. &
                indexTriangleVert(jj,it0) /= iv) then
              if (ivgot /= 0) stop 'There should only be one'
              ivgot = indexTriangleVert(jj,it0)
            endif
          enddo
        enddo
        if (ivgot == 0) stop 'No common vertex found'
        indexvnear(ii,iv) = ivgot
      enddo

      nvnear(iv) = nTrianglesContainVertex(iv)
    enddo  ! iv

    ! file output for this level
    call writeLevelFiles(il)

  enddo  ! il
  ! >> finish of il do loop

  ! file output for triangle grid data
  call writeTriangleFiles(levels)

contains

    !----------------------------------------
    subroutine create_hash_table()
    !----------------------------------------
    implicit none
    ! local parameters
    integer :: i,imax,ier,num_tri,num_vert
    double precision :: memory_size

    ! info
    num_tri = numTriangles(nlevel+1) - numTriangles(nlevel)
    num_vert = nv
    print *,'    hash table: '
    print *,'      number of triangles     = ',num_tri
    print *,'      number of vertices      = ',num_vert

    ! estimated maximum size required
    imax = num_vert + 3 * num_tri
    ! Choose a prime number: 10007,100003,1000003,1000003,..
    if (imax < 10007) then
      hash_size = 10007
    else if (imax < 100003) then
      hash_size = 100003
    else if (imax < 1000003) then
      hash_size = 1000003
    else if (imax < 10000003) then
      hash_size = 10000003
    else
      hash_size = 100000003
    endif

    ! memory required (edge_type minimum size 2 * int + pointer
    memory_size = dble(hash_size) * dble(8) / 1024.d0 / 1024.d0

    print *,'      hash size               = ',hash_size
    print *,'      minimum memory required = ',sngl(memory_size),'Mb'

    ! allocate table
    allocate(hash_table(hash_size),stat=ier)
    if (ier /= 0) stop 'Error allocating hash table'

    ! Initialize hash table
    do i = 1, hash_size
      hash_table(i)%v1 = 0
      hash_table(i)%v2 = 0
      hash_table(i)%mid_id = 0  ! Indicates empty
      nullify(hash_table(i)%next)
    enddo

    end subroutine

    !----------------------------------------
    function create_hash_key(index1,index2) result(hash_key)
    !----------------------------------------
    implicit none
    integer, intent(in) :: index1,index2
    integer :: hash_key
    ! local parameters
    integer(kind=8) :: iv1,iv2,sum_v,min_v,key

    ! create unique hash key
    iv1 = int(index1 - 1,kind=8)     ! indexing starts at 0: {0,1,2,..}
    iv2 = int(index2 - 1,kind=8)

    ! Cantor pairing function
    ! key = ((v1 + v2) * (v1 + v2 + 1)) / 2 + min(v1, v2) with v1,v2 = {0,..,nv}
    sum_v = iv1 + iv2
    min_v = iv1
    key = (sum_v * (sum_v + 1)) / 2 + min_v      ! using integer kind=8

    ! checks integer overflow
    ! 32-bit integer < 2,147,483,647
    !if (dble(sum_v)/2.d0 > (2147483646.d0 - dble(min_v)) / ((dble(sum_v) + 1.d0)/2.d0)) &
    !  stop 'Error hash key might exceed integer limit'

    ! 64-bit integer < 9,223,372,036,854,775,807
    if (dble(sum_v)/2.d0 > (9223372036854775806.d0 - dble(min_v)) / ((dble(sum_v) + 1.d0)/2.d0)) &
      stop 'Error hash key might exceed integer limit'

    if (key < 0) then
      print *,'Error: invalid key ',key,' index1/index2',index1,index2,'iv1/iv2',iv1,iv2,'sum/min',sum_v,min_v
      stop 'Error invalid key'
    endif

    ! Modulo to fit in hash table
    key = mod(key, int(hash_size,kind=8)) + 1     ! indexing starts at 1

    ! converts to 32-bit integer
    hash_key = int(key,kind=4)

    ! check
    if (hash_key < 0) then
      print *,'Error: invalid hash key ',hash_key,' index1/index2',index1,index2,'key',key
      stop 'Error invalid hash key'
    endif

    return
    end function

    !----------------------------------------
    function check_hash_table(index1, index2, key) result(midpoint_id)
    !----------------------------------------
    integer, intent(in) :: index1, index2, key
    integer :: midpoint_id
    ! local parameters
    type(edge_entry), pointer :: current

    ! check
    if (key < 1 .or. key > hash_size) then
      print *,'Error: check_hash_table() has invalid hash key ',hash_key,' - should be > 0 and < ',hash_size
      stop 'Error invalid hash key'
    endif

    ! empty bucket
    if (hash_table(key)%mid_id == 0) then
      midpoint_id = 0
      return
    endif

    ! Check head of list
    if (hash_table(key)%v1 == index1 .and. hash_table(key)%v2 == index2) then
      midpoint_id = hash_table(key)%mid_id
      return
    endif

    ! Check collision chain
    current => hash_table(key)
    do while (associated(current%next))
      current => current%next
      if (current%v1 == index1 .and. current%v2 == index2) then
        midpoint_id = current%mid_id
        return
      endif
    enddo

    ! Not found
    midpoint_id = 0
    return
    end function

    !----------------------------------------
    subroutine add_hash_table_entry(index1,index2,key,id)
    !----------------------------------------
    implicit none
    integer, intent(in) :: index1,index2,key,id
    ! local parameters
    integer :: ier
    type(edge_entry), pointer :: current, new_entry

    ! check
    if (key < 1 .or. key > hash_size) then
      print *,'Error: add_hash_table_entry() has invalid hash key ',hash_key,' - should be > 0 and < ',hash_size
      stop 'Error invalid hash key'
    endif

    ! empty bucket
    if (hash_table(key)%mid_id == 0) then
      ! fill with new entry
      hash_table(key)%v1 = index1
      hash_table(key)%v2 = index2
      hash_table(key)%mid_id = id
      return
    else
      ! goto end of collision chain
      current => hash_table(key)
      do while (associated(current%next))
        current => current%next
      enddo

      ! create new entry
      allocate(new_entry,stat=ier)
      if (ier /= 0) stop 'Error allocating new hash table entry'
      new_entry%v1 = index1
      new_entry%v2 = index2
      new_entry%mid_id = id
      ! no next yet
      nullify(new_entry%next)
      ! add as next to current entry
      current%next => new_entry
    endif

    end subroutine

    !----------------------------------------
    subroutine free_hash_table()
    !----------------------------------------
    ! local parameters
    integer :: i
    type(edge_entry), pointer :: current, next_entry

    ! clean up collision chain entries
    do i = 1, hash_size
      current => hash_table(i)%next
      do while (associated(current))
        next_entry => current%next
        deallocate(current)
        current => next_entry
      enddo
    enddo
    ! clean up table
    deallocate(hash_table)
    end subroutine

  end subroutine

!-----------------------------------------------------------------------
  subroutine sort_vertices(vertices, n, permutation)
!-----------------------------------------------------------------------
  implicit none
  integer, intent(in) :: n
  double precision, intent(inout) :: vertices(3, n)
  integer, intent(inout) :: permutation(n)
  ! local parameters
  integer :: i, j
  double precision :: temp_v(3)
  integer :: temp_p
  ! tolerance (machine precision for double precision)
  double precision, parameter :: TOL_EPS = 1.0d-12

  ! simple bubble sort routine
  do i = 1, n - 1
    do j = i + 1, n
      if (vertices(1, j) < vertices(1, i) .or. &
          (abs(vertices(1, j) - vertices(1, i)) < TOL_EPS .and. &
          (vertices(2, j) < vertices(2, i) .or. &
          (abs(vertices(2, j) - vertices(2, i)) < TOL_EPS .and. &
          vertices(3, j) < vertices(3, i))))) then
        ! Swap vertices
        temp_v = vertices(:, i)
        vertices(:, i) = vertices(:, j)
        vertices(:, j) = temp_v
        ! Swap permutation indices
        temp_p = permutation(i)
        permutation(i) = permutation(j)
        permutation(j) = temp_p
      endif
    enddo
  enddo

  end subroutine

!-----------------------------------------------------------------------
  subroutine binary_search_vertex(midpoint, vertices, n, index, found)
!-----------------------------------------------------------------------
! for binary search in the sorted tvertices
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: midpoint(3)
  double precision, intent(in) :: vertices(3, n)
  integer, intent(out) :: index
  logical, intent(out) :: found
  ! local parameters
  integer :: low, high, mid
  ! tolerance (machine precision for double precision)
  double precision, parameter :: TOL_EPS = 1.0d-12

  low = 1
  high = n
  found = .false.
  index = 0

  do while (low <= high .and. .not. found)
    mid = (low + high) / 2
    if (abs(midpoint(1) - vertices(1, mid)) < TOL_EPS) then
      if (abs(midpoint(2) - vertices(2, mid)) < TOL_EPS) then
        if (abs(midpoint(3) - vertices(3, mid)) < TOL_EPS) then
          found = .true.
          index = mid
        else if (midpoint(3) < vertices(3, mid)) then
          high = mid - 1
        else
          low = mid + 1
        endif
      else if (midpoint(2) < vertices(2, mid)) then
        high = mid - 1
      else
        low = mid + 1
      endif
    else if (midpoint(1) < vertices(1, mid)) then
      high = mid - 1
    else
      low = mid + 1
    endif
  enddo

  end subroutine


!-----------------------------------------------------------------------
  subroutine writeLevelFiles(il)
!-----------------------------------------------------------------------

  use dodeco
  implicit none
  integer, intent(in) :: il
  ! local parameters
  integer :: iv,i,it,kk
  character(len=2) :: strLev

  !  Opens the data files and writes the "near vertices" indices to them
  !  Order q=6 has xxx,xxx triangular vertices;
  !  in the format statements for integers, we leave 8 slots.
  print *,'writing data files for level ',il

  write(strLev,'(i0)') il
  ! left adjust string number
  strLev = adjustl(strLev)

  open(93,file = 'OUTPUT/' // 'Dnear'//trim(strLev)//'.dat')
  do iv = 1, numVertices(il)
    if (nvnear(iv) /= 5 .and. nvnear(iv) /= 6) stop 'something wrong with nvnear'
    if (nvnear(iv) == 5) indexvnear(6,iv) = indexvnear(1,iv)
  enddo
  do iv = 1, numVertices(il)
    write(93, '(6i8)') (indexvnear(kk,iv), kk = 1,6)
  enddo
  close(93)

  ! Opens the data files and writes the voronoi grid data to them.
  open(93,file = 'OUTPUT/' // 'Dvvert'//trim(strLev)//'.dat')

  write(93,'(3f18.14)') ((voronoiVertices(i,it),i=1,3),it=numTriangles(il)+1,numTriangles(il+1))
  close(93)

  open(93, file = 'OUTPUT/' // 'Dvface'//trim(strLev)//'.dat')
  do iv = 1, numVertices(il)
    ! if not a pentagon
    if (indexTriContainVertex(6,iv) /= 0) then
      write(93,'(6i8)') (indexTriContainVertex(i,iv) - numTriangles(il), i=1,6)
    else
      write(93,'(6i8)') (indexTriContainVertex(i,iv) - numTriangles(il), i=1,5), &
                         indexTriContainVertex(1,iv) - numTriangles(il)
    endif
  enddo
  close(93)

  end subroutine

!-----------------------------------------------------------------------
  subroutine writeTriangleFiles(levels)
!-----------------------------------------------------------------------

  use dodeco
  implicit none
  integer, intent(in) :: levels
  ! local parameters
  integer :: il,i,iv,it
  character(len=2) :: strLev

  ! Opens the data files and writes the triangle grid data to them.
  print *,'writing triangle grid datas Dtvert and Dtface'

  do il = 0,levels

    write(strLev,'(i0)') il
    ! left adjust string
    strLev = adjustl(strLev)

    open(93,file = 'OUTPUT/' // 'Dtvert'//trim(strLev)//'.dat')
    write(93,'(3f18.14)') ((tvertices(i,iv),i=1,3), iv=1,numVertices(il))
    close(93)

    open(93,file = 'OUTPUT/' // 'Dtface'//trim(strLev)//'.dat')
    write(93,'(3i8)') ((indexTriangleVert(i,it),i=1,3), it=numTriangles(il)+1,numTriangles(il+1))
    close(93)
  enddo

  end subroutine

!-----------------------------------------------------------------------
  subroutine drsoinc(a,n,idx)
!-----------------------------------------------------------------------
! performs an indirect heapsort in decreasing order such that: a(idx(1)) >= a(idx(2)) >= ... >= a(idx(n))
  implicit none
  integer, intent(in) :: n
  integer, intent(inout) :: idx(n)
  double precision, intent(inout) :: a(n)
  ! local parameters
  integer :: i, ict, n1, n2, n21, np2, nn, ik, jk, ic, it
  double precision :: c, t

  ! Early return condition
  if (n == 1) then
    idx(1) = 1
    return
  endif

  ! check
  if (n <= 0) then
    print *, 'error return from rsoinc - n less than or equal to 1'
    stop
  endif

  ! Initialize index array
  do i = 1, n
    idx(i) = i
  enddo

  ! Initialize heap sort parameters
  n2 = n / 2
  n21 = n2 + 2
  ict = 1
  i = 2

  ! Phase 1: Building the heap
  build_heap: do
    n1 = n21 - i
    nn = n
    ik = n1

    ! Adjust heap for current node
    heapify: block
      c = a(ik)
      ic = idx(ik)

      heap_adjust: do
        jk = 2 * ik
        if (jk > nn) exit heap_adjust

        ! Choose larger child
        if (jk < nn) then
          if (a(jk+1) > a(jk)) jk = jk + 1
        endif

        ! If parent >= child, heap property satisfied
        if (a(jk) <= c) exit heap_adjust

        ! Move child up and continue down the heap
        a(ik) = a(jk)
        idx(ik) = idx(jk)
        ik = jk
      enddo heap_adjust

      ! Place value at final position
      a(ik) = c
      idx(ik) = ic
    end block heapify

    ! Phase transition check
    if (i >= n2) then
      ! Transition to Phase 2: Extract sorted elements
      ict = 2
      np2 = n + 2
      i = 2
      exit build_heap
    endif

    i = i + 1
  enddo build_heap

  ! Phase 2: Extract elements from heap in sorted order
  extract_sorted: do
    n1 = np2 - i
    nn = n1
    ik = 1

    ! Adjust heap after root extraction
    heapify_extract: block
      c = a(ik)
      ic = idx(ik)

      extract_adjust: do
        jk = 2 * ik
        if (jk > nn) exit extract_adjust

        ! Choose larger child
        if (jk < nn) then
            if (a(jk+1) > a(jk)) jk = jk + 1
        endif

        ! If parent >= child, heap property satisfied
        if (a(jk) <= c) exit extract_adjust

        ! Move child up and continue down the heap
        a(ik) = a(jk)
        idx(ik) = idx(jk)
        ik = jk
      enddo extract_adjust

      ! Place value at final position
      a(ik) = c
      idx(ik) = ic
    end block heapify_extract

    ! Swap current max (root) with end of active array
    t = a(1)
    a(1) = a(n1)
    a(n1) = t
    it = idx(1)
    idx(1) = idx(n1)
    idx(n1) = it

    ! Exit condition
    if (i >= n) exit extract_sorted

    i = i + 1
  enddo extract_sorted

  end subroutine

!-----------------------------------------------------------------------
  double precision function dazimv(v1,v2)
!-----------------------------------------------------------------------
  implicit none
  double precision, intent(in) :: v1(3),v2(3)
  double precision, parameter :: pi = 3.1415926535897931d0

  dazimv = pi + datan2(v1(2)*v2(1) - v1(1)*v2(2), &
                       v1(1)*v1(3)*v2(1) + v1(2)*v1(3)*v2(2) - (v1(1)*v1(1) + v1(2)*v1(2))*v2(3) )

  return
  end function

!-----------------------------------------------------------------------
  double precision function ddot3ScalarProduct(x1,x2)
!-----------------------------------------------------------------------
  implicit none
  double precision, intent(in) :: x1(3),x2(3)

  ddot3ScalarProduct = x1(1)*x2(1) + x1(2)*x2(2) + x1(3)*x2(3)

  return
  end function

!-----------------------------------------------------------------------
  subroutine getDodProperties(level,nVertices,nTriangles)
!-----------------------------------------------------------------------

! gets the number of vertices and triangles of a given level
!
! return: nvertices and ntriangles

  use dodeco
  implicit none
  integer, intent(in) :: level
  integer, intent(out) :: nTriangles,nVertices

  ! number of triangles
  nTriangles = numTriangles(level + 1) - numTriangles(level)

  ! number of vertices
  nVertices = numVertices(level)

  end subroutine

