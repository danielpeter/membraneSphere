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
  subroutine constructParallelDomains()
!-----------------------------------------------------------------------
! makes grid domains for parallelization
  use parallel;use griddomain;use propagationStartup; use cells; use verbosity
  implicit none
  ! local parameters
  integer:: range,index,n,domain,numDomainVert(0:nprocesses-1),ier
  integer, external:: getDomain

  ! domains only make sense for more then 1 processors
  if (nprocesses < 2) then
    return
  endif

  ! construct domains
  ! count vertices in each domain
  numDomainVert(:) = 0
  do n = 1,numVertices
    domain = getDomain(n)
    numDomainVert(domain) = numDomainVert(domain)+1
  enddo
  range = numDomainVert(myrank)
  if (range <= 0) then
    print *,'Error: rank ',myrank,': range has no elements',range
    call stopProgram('abort - constructParallelDomains    ')
  endif

  ! allocate region information for this process
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'  allocating parallel domain arrays:'
    print *,'    vertexDomain size   :  ',numVertices * 4/1024./1024.,'Mb' ! assume integer == 4 bytes
    print *,'    domainVertices size :  ',range * 4/1024./1024.,'Mb'
  endif
  allocate(vertexDomain(numVertices),stat=ier)
  if (ier /= 0) call stopProgram('abort - constructParallelDomains:error vertexDomain array   ')
  vertexDomain(:) = -1

  !allocates array which holds all vertices in the boundary to another parallel domain
  allocate(domainVertices(range),stat=ier)
  if (ier /= 0) call stopProgram('abort - constructParallelDomains: error domainVertices array   ')
  domainVertices(:) = -1

  ! divide vertices into domains
  index = 0
  do n = 1,numVertices
    ! fill vertexRegion array (for all vertices)
    vertexDomain(n) = getDomain(n)

    ! fill specific domain array (contains only vertices for this domain/process)
    if (vertexDomain(n) == myrank) then
      index = index+1
      if (index > range) then
       print *,'Error: rank ',myrank,': exceeding range limite in domain', range
       call stopProgram( 'abort - constructParallelDomains    ')
      endif
      domainVertices(index)=n
    endif
  enddo

  ! check
  if (index /= range) then
    print *,'Error: rank ',myrank,': range in domain wrong',range,index
    call stopProgram('abort - constructParallelDomains    ')
  endif

  ! find boundary vertices
  call findBoundaries()

  ! find domain neighbors
  call findDomainNeighbors()

  end subroutine

!-----------------------------------------------------------------------
  integer function getDomain(vertex)
!-----------------------------------------------------------------------
! finds domain in which the vertex lays
!
! input:
!       vertex - index of vertex in vertices array
!
! returns: domain of the grid where vertex is in
  use parallel;use griddomain
  implicit none
  integer,intent(in):: vertex
  ! local parameters
  real(WP):: colat,long,longitudeIncrement,latitudeIncrement
  integer:: n,m


  ! determine domains on the sphere
  select case(nprocesses)
  case(1)
    ! always whole sphere
    getDomain = 0
    return
  case(2) ! for 2 processes
    ! find lat/lon of vertex
    call getSphericalCoord(vertex,colat,long)

    ! domains simply divide into two half-spheres
    if (long >= 0.0 .and. long < PI) then
      getDomain = 0
      return
    else
      getDomain = 1
      return
    endif
  case(4) ! for 4 processes
    ! find lat/lon of vertex
    call getSphericalCoord(vertex,colat,long)

    ! domains are a the half-sphere split
    if (long >= 0.0 .and. long < PI) then
      if (colat < PI/2.0) then
        getDomain = 0
        return
      else
        getDomain = 1
        return
      endif
    else
      if (colat < PI/2.0) then
        getDomain = 2
        return
      else
        getDomain = 3
        return
      endif
    endif
  case(6)
    ! find lat/lon of vertex
    call getSphericalCoord(vertex,colat,long)

    ! domains are split by colatitude PI/4,3PI/4 and longitude 0,PI/2,PI,3PI/2
    if (colat < (PI/4.0)) then
      getDomain = 0
      return
    else if (colat > (3.0*PI/4.0)) then
      getDomain = 1
      return
    else if (long >= 0.0 .and. long < (PI/2.0)) then
      getDomain = 2
      return
    else if (long >= (PI/2.0) .and. long < PI) then
      getDomain = 3
      return
    else if (long >= PI.and.long < (3.0*PI/2.0)) then
      getDomain = 4
      return
    else
      getDomain = 5
      return
    endif
  case(8)
    ! find lat/lon of vertex
    call getSphericalCoord(vertex,colat,long)

    ! domains are the quarter-sphere split (equal to spherical octahedron)
    ! still, this is optimal distribution
    if (long >= 0.0 .and. long < PI/2.0) then
      if (colat < PI/2.0) then
        getDomain = 0
        return
      else
        getDomain = 1
        return
      endif
    else if (long >= PI/2.0 .and. long < PI) then
      if (colat < PI/2.0) then
        getDomain = 2
        return
      else
        getDomain = 3
        return
      endif
    else if (long >= PI .and. long < (3.0*PI/2.0)) then
      if (colat < PI/2.0) then
        getDomain = 4
        return
      else
        getDomain = 5
        return
      endif
    else
      if (colat < PI/2.0) then
        getDomain = 6
        return
      else
        getDomain = 7
        return
      endif
    endif
  case(16)
    ! find lat/lon of vertex
    call getSphericalCoord(vertex,colat,long)

    ! domains are the quarter-sphere split and divide into 2 slices
    ! though, no more optimal
    if (long >= 0.0 .and. long < 0.25*PI) then
      if (colat < 0.5*PI) then
        getDomain = 0
        return
      else
        getDomain = 1
        return
      endif
    else if (long >= 0.25*PI .and. long < 0.5*PI) then
      if (colat < 0.5*PI) then
        getDomain = 2
        return
      else
        getDomain = 3
        return
      endif
    else if (long >= 0.5*PI .and. long < 0.75*PI) then
      if (colat < 0.5*PI) then
        getDomain = 4
        return
      else
        getDomain = 5
        return
      endif
    else if (long >= 0.75*PI .and. long < PI) then
      if (colat < 0.5*PI) then
        getDomain = 6
        return
      else
        getDomain = 7
        return
      endif
    else if (long >= PI .and. long < 1.25*PI) then
      if (colat < 0.5*PI) then
        getDomain = 8
        return
      else
        getDomain = 9
        return
      endif
    else if (long >= 1.25*PI .and. long < 1.5*PI) then
      if (colat < 0.5*PI) then
        getDomain = 10
        return
      else
        getDomain = 11
        return
      endif
    else if (long >= 1.5*PI .and. long < 1.75*PI) then
      if (colat < 0.5*PI) then
        getDomain = 12
        return
      else
        getDomain = 13
        return
      endif
    else
      if (colat < 0.5*PI) then
        getDomain = 14
        return
      else
        getDomain = 15
        return
      endif
    endif
  case(32)
    ! find lat/lon of vertex
    call getSphericalCoord(vertex,colat,long)

    ! splits sphere in eight first, then halves the latitudes and longitudes of each domain
    longitudeIncrement = PI/4.0
    latitudeIncrement = PI/4.0
    if (colat >= 0.0 .and. colat < latitudeIncrement) then
      n = 1
      do m = 1,7
        if (long >= (m-1)*longitudeIncrement .and. long < m*longitudeIncrement) then
          ! starts with index 0 for main, then up to nprocesses - 1
          getDomain = (n-1)*8 + m - 1
          return
        endif
      enddo
      if (long >= 7*longitudeIncrement) then
        m = 8
        ! starts with index 0 for main, then up to nprocesses - 1
        getDomain = (n-1)*8 + m - 1
        return
      endif
    else if (colat >= latitudeIncrement .and. colat < 2*latitudeIncrement) then
      n = 2
      do m = 1,7
        if (long >= (m-1)*longitudeIncrement .and. long < m*longitudeIncrement) then
          ! starts with index 0 for main, then up to nprocesses - 1
          getDomain = (n-1)*8 + m - 1
          return
        endif
      enddo
      if (long >= 7*longitudeIncrement) then
        m = 8
        ! starts with index 0 for main, then up to nprocesses - 1
        getDomain = (n-1)*8 + m - 1
        return
      endif
    else if (colat >= 2*latitudeIncrement .and. colat < 3*latitudeIncrement) then
      n = 3
      do m = 1,7
        if (long >= (m-1)*longitudeIncrement .and. long < m*longitudeIncrement) then
          ! starts with index 0 for main, then up to nprocesses - 1
          getDomain = (n-1)*8 + m - 1
          return
        endif
      enddo
      if (long >= 7*longitudeIncrement) then
        m = 8
        ! starts with index 0 for main, then up to nprocesses - 1
        getDomain = (n-1)*8 + m - 1
        return
      endif
    else if (colat >= 3*latitudeIncrement) then
      n = 4
      do m = 1,7
        if (long >= (m-1)*longitudeIncrement .and. long < m*longitudeIncrement) then
          ! starts with index 0 for main, then up to nprocesses - 1
          getDomain = (n-1)*8 + m - 1
          return
        endif
      enddo
      if (long >= 7*longitudeIncrement) then
        m = 8
        ! starts with index 0 for main, then up to nprocesses - 1
        getDomain = (n-1)*8 + m - 1
        return
      endif
    endif
  case default
    ! find lat/lon of vertex
    call getSphericalCoord(vertex,colat,long)

    ! differs for even/odd number of processes
    ! even number of processes
    if (mod(nprocesses,2) == 0) then
      ! even number of processes distrbuted on north- and south-hemnisphere
      ! grid split into upper/lower part and by longitudes
      longitudeIncrement=2.0*PI/(nprocesses/2)
      if (colat <= PI/2.0) then
        do n = 0,(nprocesses/2)-1
          if ( (long >= longitudeIncrement*n) .and. ( long < longitudeIncrement*(n+1))) then
            getDomain = n
            return
          endif
        enddo
      else
        do n=(nprocesses/2),nprocesses-1
          if ( (long >= longitudeIncrement*(n-nprocesses/2)) .and. ( long < longitudeIncrement*(n+1-nprocesses/2))) then
            getDomain = n
            return
          endif
        enddo
      endif
    endif

    ! uneven number of processes
    if (mod(nprocesses,2) /= 0) then
      ! uneven number of processes
      ! grid split only by longitudes
      getDomain = 0
      longitudeIncrement = 2.0*PI/nprocesses
      do n = 0,nprocesses-1
        if ( (long >= longitudeIncrement*n) .and. ( long < longitudeIncrement*(n+1))) then
          getDomain = n
          return
        endif
      enddo
    endif
  end select

  ! reaching this point should not occur
  print *,'Error: could not split into domains',nprocesses,vertex,colat,long
  call stopProgram( 'abort - getDomain    ')

  return
  end function


!-----------------------------------------------------------------------
  subroutine findBoundaries()
!-----------------------------------------------------------------------
! finds vertices on boundaries for each domain
!
! returns: boundaries is filled such that for each domain
!               boundaries(domain,neighbor,:) is the array of all vertices of that
!               domain which lay on the boundary to the specific neighbor
  use parallel;use griddomain; use propagationStartup; use verbosity; use cells
  implicit none
  ! local parameters
  integer:: numDomainVert(0:nprocesses,0:nprocesses),index,maxRange,n,k,domain,ier
  integer:: neighbors(nprocesses-1)
  integer, external:: getDomain
  logical, external:: isBoundary

  !count boundary vertices in each domain
  numDomainVert(:,:) = 0
  do n = 1,numVertices
    domain = getDomain(n)
    !check if on bounday
    if (isBoundary(n)) then
      !get neighbors
      call getVertexNeighbors(n,neighbors)
      !add corresponding counts
      do k = 1,size(neighbors(:))
        if (neighbors(k) /= (-1)) then
          numDomainVert(domain,neighbors(k)) = numDomainVert(domain,neighbors(k))+1
        endif
      enddo
    endif
  enddo

  !determine maximal range
  maxRange = 0
  do n = 0,nprocesses-1
    do k = 0,nprocesses-1
      if (numDomainVert(n,k) > maxRange) then
        maxRange=numDomainVert(n,k)
      endif
    enddo
  enddo
  boundariesMaxRange = maxRange

  !debug
  !if (VERBOSE .and. MAIN_PROCESS) then
  !do k=0,nprocesses-1
  !  print *,'boundary vertices in(',k,'/.):',numDomainVert(k,:)
  !enddo
  !print *,'maximal number of boundary vertices',maxRange
  !endif

  !allocate array
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'  allocating boundary array:'
    print *,'    size :  ',nprocesses*nprocesses*maxRange * 4/1024./1024.,'Mb'  ! assume integer == 4 bytes
  endif
  allocate(boundaries(0:nprocesses-1,0:nprocesses-1,maxRange), &
            sendDisp(maxRange),receiveDisp(maxRange),stat=ier)
  if (ier /= 0) call stopProgram( 'abort - findBoundaries     ')
  boundaries(:,:,:) = 0
  sendDisp(:) = 0
  receiveDisp(:) = 0

  ! fill boundaries array
  numDomainVert(:,:) = 0
  do n = 1,numVertices
    domain = getDomain(n)
    ! check if on bounday
    if (isBoundary(n)) then
      ! get neighbors
      call getVertexNeighbors(n,neighbors)
      ! add corresponding counts
      do k = 1,size(neighbors(:))
        if (neighbors(k) /= (-1)) then
          numDomainVert(domain,neighbors(k)) = numDomainVert(domain,neighbors(k))+1
          index = numDomainVert(domain,neighbors(k))
          ! check index
          if (index > maxRange) then
            print *,'Error: rank ',myrank,': exceeding range in boundaries',domain,maxRange
            call stopProgram( 'abort - findBoundaries     ')
          endif
          ! fill boundaries array
          boundaries(domain,neighbors(k),index)=n
        endif
      enddo
    endif
  enddo

  end subroutine


!-----------------------------------------------------------------------
  logical function isBoundary(vertex)
!-----------------------------------------------------------------------
! determines if vertex is on a domain boundary
!
! input:
!     vertex - index of vertex in vertices array
!
! returns: true if this vertex has neighbors in a different domain
  use parallel;use griddomain; use cells
  implicit none
  integer,intent(in):: vertex
  ! local parameters
  integer:: k,domain

  ! domain this reference vertex is in
  domain = vertexDomain(vertex)

  !check neighbors
  isBoundary = .false.
  do k = 1, cellNeighbors(vertex,0)
    if (vertexDomain(cellNeighbors(vertex,k)) /= domain) then
      isBoundary = .true.
      return
    endif
  enddo

  return
  end function


!-----------------------------------------------------------------------
  subroutine getVertexNeighbors(vertex,neighbors)
!-----------------------------------------------------------------------
! find domains which are next to the referenced vertex one
!
! input:
!       vertex - reference vertice
!       neighbors - array which holds neighbors (length=nprocesses-1)
!
! returns: neighbors of this domain (each element contains the domain number, -1 for no neighbor)
  use parallel;use griddomain; use cells
  implicit none
  integer,intent(in):: vertex
  integer,intent(out):: neighbors(nprocesses-1)
  ! local parameters
  integer:: domain,domainRef,n,k,index
  logical:: isnew
  integer,external:: getDomain

  !init
  neighbors(:) = -1
  domainRef = getDomain(vertex)

  !find neighbor domain
  index = 0
  do n = 1,cellNeighbors(vertex,0)
    domain = getDomain(cellNeighbors(vertex,n))
    if (domain /= domainRef) then
      ! fill into neighbor if new
      isnew = .true.
      do k = 1,size(neighbors(:))
        if (neighbors(k) == domain) then
          isnew = .false.
          exit
        endif
      enddo
      if (isnew) then
        index = index+1
        neighbors(index)=domain
      endif
    endif
  enddo

  end subroutine


!-----------------------------------------------------------------------
  subroutine findDomainNeighbors()
!-----------------------------------------------------------------------
! find domains which are next to the referenced one
!
! returns: domainNeighbors array (each element contains the domain number, -1 for no neighbor)
  use parallel;use griddomain; use verbosity
  implicit none
  ! local parameters
  integer:: domain,index,k,ier

  !check
  if (nprocesses == 1) return

  !init
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'  allocating domainNeighbors array:'
    print *,'    size :  ',nprocesses*nprocesses * 4/1024./1024.,'Mb'   ! assume integer == 4 bytes
  endif
  allocate(domainNeighbors(0:nprocesses-1,nprocesses-1),stat=ier)
  if (ier /= 0) call stopProgram( 'abort - findDomainNeighbors    ')
  domainNeighbors(:,:)=-1

  !find neighbor domains
  do domain = 0,nprocesses-1
    index = 0
    do k = 0,nprocesses-1
      !check if domain has entries of boundary vertices in this neighbor domain
      if (boundaries(domain,k,1) /= 0) then
        ! first element in boundaries array is a valid vertex index, so add this neighbor
        index = index+1
        if (index > nprocesses-1) then
          print *,'Error: rank ',myrank,':exceeding number of neighbor domains',domain
          call stopProgram( 'abort - findDomainNeighbors    ')
        endif
        ! add this as neighbor
        domainNeighbors(domain,index)=k
      endif
    enddo
  enddo

  end subroutine


