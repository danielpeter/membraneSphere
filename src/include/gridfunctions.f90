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
      subroutine constructParallelDomains()
!-----------------------------------------------------------------------
! makes grid domains for parallelization
      use parallel;use griddomain;use propagationStartup
      implicit none
      integer:: range,index,n,domain,numDomainVert(0:nprocesses-1),getDomain,ierror
      real(WP):: colat,long
      external:: getDomain
      
      ! domains only make sense for more then 1 processors
      if( nprocesses .lt. 2 ) then
        return
      endif
      
      ! construct domains      
      !count vertices in each domain
      !if(MASTER) open(11,file='vertexdomains.dat')
      numDomainVert(:)=0
      do n=1,numVertices
        domain = getDomain(n)
        numDomainVert(domain)=numDomainVert(domain)+1
        !if(MASTER) write(11,*) n, domain
      enddo      
      !if(MASTER) close(11)
      
      ! allocate region information for this process
      allocate(vertexDomain(numVertices),stat=ierror)
      if( ierror .gt. 0 ) call stopProgram('abort - constructParallelDomains:error in allocating vertexRegion array   ')
            
      !allocate domain
      range = numDomainVert(rank)
      if( range .le. 0 ) then
        print*,rank,': range has no elements',range
        call stopProgram('abort - constructParallelDomains    ')
      endif
      allocate(domainVertices(range),stat=ierror)
      if( ierror .gt. 0 ) call stopProgram('abort - constructParallelDomains: error in allocating domain boundary array   ')
      
      ! divide vertices into domains
      index=0
      do n=1,numVertices
        ! fill vertexRegion array (for all vertices)
        vertexDomain(n)=getDomain(n)
        
        ! fill specific domain array (contains only vertices for this domain/process)
        if( vertexDomain(n) .eq. rank) then
          index=index+1
          if( index .gt. range) then
           print*,rank,': exceeding range limite in domain', range
           call stopProgram( 'abort - constructParallelDomains    ')
          endif
          domainVertices(index)=n
        endif
      enddo
      !check
      if( index .ne. range ) then
        print*,rank,': range in domain wrong',range,index
        call stopProgram('abort - constructParallelDomains    ')
      endif
      
      ! find boundary vertices
      call findBoundaries()
            
      ! find domain neighbors
      call findDomainNeighbors()      
            
      end
  
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
      real(WP):: colat,long,longitudeIncrement
      integer:: vertex,n
      logical:: placed

      ! find lat/lon of vertex
      call getSphericalCoord(vertex,colat,long)

      ! determine domains on the sphere
      select case(nprocesses)
      case(1)
        ! always whole sphere
        getDomain=0  
        return
      case(2) ! for 2 processes
        ! domains simply divide into two half-spheres
        if( long .ge. 0.0 .and. long .lt. PI ) then
          getDomain = 0
        else
          getDomain = 1
        endif
        return
      case(4) ! for 4 processes
        ! domains are a the half-sphere split 
        if( long .ge. 0.0 .and. long .lt. PI ) then
          if( colat .lt. PI/2.0) then
            getDomain=0
          else
            getDomain=1
          endif
        else
          if( colat .lt. PI/2.0) then
            getDomain=2
          else
            getDomain=3
          endif
        endif
        return
      case(6)
        ! domains are split by colatitude PI/4,3PI/4 and longitude 0,PI/2,PI,3PI/2
        if( colat.lt.(PI/4.0)) then
          getDomain=0
        else if( colat.gt.(3.0*PI/4.0)) then
          getDomain=1
        else if( long.ge.0.0 .and. long .lt. (PI/2.0)) then
          getDomain=2
        else if( long.ge.(PI/2.0) .and. long.lt.PI) then
          getDomain=3
        else if( long .ge. PI.and.long.lt.(3.0*PI/2.0)) then
          getDomain=4
        else
          getDomain=5
        endif
        return
      case(8)
        ! domains are the quarter-sphere split (equal to spherical octahedron)
        if( long .ge. 0.0 .and. long .lt. PI/2.0 ) then
          if( colat .lt. PI/2.0) then
            getDomain=0
          else
            getDomain=1
          endif
        else if( long .ge. PI/2.0 .and. long .lt. PI) then        
          if( colat .lt. PI/2.0) then
            getDomain=2
          else
            getDomain=3
          endif
        else if( long .ge. PI .and. long .lt. (3.0*PI/2.0)) then
          if( colat .lt. PI/2.0) then
            getDomain=4
          else
            getDomain=5
          endif
        else
          if( colat .lt. PI/2.0) then
            getDomain=6
          else
            getDomain=7
          endif        
        endif
        return        
      case(16)
        ! domains are the quarter-sphere split (equal to spherical octahedron) and divide into 2 slices
        if( long .ge. 0.0 .and. long .lt. 0.25*PI ) then
          if( colat .lt. 0.5*PI) then
            getDomain=0
          else
            getDomain=1
          endif
        else if( long .ge. 0.25*PI .and. long .lt. 0.5*PI) then        
          if( colat .lt. 0.5*PI) then
            getDomain=2
          else
            getDomain=3
          endif
        else if( long .ge. 0.5*PI .and. long.lt.0.75*PI)then
          if( colat .lt. 0.5*PI) then
            getDomain=4
          else
            getDomain=5
          endif
        else if(long.ge.0.75*PI.and.long.lt.PI)then
          if( colat .lt. 0.5*PI) then
            getDomain=6
          else
            getDomain=7
          endif        
        else if(long.ge.PI.and.long.lt.1.25*PI)then
          if( colat .lt. 0.5*PI) then
            getDomain=8
          else
            getDomain=9
          endif        
        else if(long.ge.1.25*PI.and.long.lt.1.5*PI)then
          if( colat .lt. 0.5*PI) then
            getDomain=10
          else
            getDomain=11
          endif        
        else if(long.ge.1.5*PI.and.long.lt.1.75*PI)then
          if( colat .lt. 0.5*PI) then
            getDomain=12
          else
            getDomain=13
          endif        
        else 
          if( colat .lt. 0.5*PI) then
            getDomain=14
          else
            getDomain=15
          endif        
        endif
        return
      case default
        ! split just longitude
        if( MASTER ) print*,'grid is split into a non-optimized grid...'
        
        longitudeIncrement=360.0_WP/nprocesses        
        placed=.false.
        do n=0,nprocesses-1
          if( (long .ge. longitudeIncrement*n ) .and. ( long .lt. longitudeIncrement*(n+1) )) then
            getDomain=n
            placed=.true.
            return
          endif
        enddo
        if( .not. placed ) then
          print*,'could not split into domains',nprocesses,vertex,colat,long
          call stopProgram( 'abort - getDomain    ')
        endif
      end select
      
      return
      end
      
!-----------------------------------------------------------------------
      subroutine findBoundaries()
!-----------------------------------------------------------------------
! finds vertices on boundaries for each domain
! 
! returns: boundaries is filled such that for each domain
!               boundaries(domain,neighbor,:) is the array of all vertices of that
!               domain which lay on the boundary to the specific neighbor
      use parallel;use griddomain; use propagationStartup; use verbosity
      implicit none
      integer:: numDomainVert(0:nprocesses,0:nprocesses),index,maxRange,range,n,k,domain,ierror
      logical:: isBoundary
      integer:: getDomain,neighbors(nprocesses-1)
      external:: isBoundary,getDomain
          
      !count boundary vertices in each domain
      numDomainVert(:,:)=0
      do n=1,numVertices
        domain = getDomain(n)
        !check if on bounday
        if( isBoundary(n) ) then
          !get neighbors 
          call getVertexNeighbors(n,neighbors)
          !add corresponding counts
          do k=1,size(neighbors(:))
            if( neighbors(k) .ne. (-1) ) then
              numDomainVert(domain,neighbors(k)) = numDomainVert(domain,neighbors(k))+1
            endif
          enddo
        endif
      enddo
      
      !determine maximal range
      maxRange=0
      do n=0,nprocesses-1
        do k=0,nprocesses-1
          if( numDomainVert(n,k) .gt. maxRange) then
            maxRange=numDomainVert(n,k)
          endif
        enddo
      enddo
      boundariesMaxRange = maxRange
      
      ! console output
      if(VERBOSE .and. MASTER ) then 
        !debug
        if(DEBUG) then
          do k=0,nprocesses-1
            print*,'boundary vertices in(',k,'/.):',numDomainVert(k,:)
          enddo        
        endif
        print*,'maximal number of boundary vertices',maxRange
      endif
      
      !allocate array
      allocate(boundaries(0:nprocesses-1,0:nprocesses-1,maxRange),      &
                sendDisp(maxRange),receiveDisp(maxRange),stat=ierror)
      if( ierror .gt. 0 ) then
        print*,'error in allocating boundaries arrays'
        call stopProgram( 'abort - findBoundaries     ')
      endif
      boundaries(:,:,:)=0
      
      ! fill boundaries array
      numDomainVert(:,:)=0
      do n=1,numVertices
        domain=getDomain(n)
        !check if on bounday
        if( isBoundary(n) ) then
          !get neighbors 
          call getVertexNeighbors(n,neighbors)
          !add corresponding counts
          do k=1,size(neighbors(:))
            if( neighbors(k) .ne. (-1) ) then
              numDomainVert(domain,neighbors(k)) = numDomainVert(domain,neighbors(k))+1
              index = numDomainVert(domain,neighbors(k))
              !check index
              if( index .gt. maxRange) then
                print*,rank,': exceeding range in boundaries',domain,maxRange
                call stopProgram( 'abort - findBoundaries     ')
              endif
              !fill boundaries array
              boundaries(domain,neighbors(k),index)=n
            endif
          enddo
        endif

      enddo
      
      !debug
      if(DEBUG) then
        if( nprocesses.gt.1.and. MASTER) then
          do n=0,nprocesses-1
            print*,'boundaries(0/',n,'):',boundaries(0,n,:)
          enddo
          do n=0,nprocesses-1
            print*,'boundaries(1/',n,'):',boundaries(1,n,:)
          enddo
        endif
      endif
      
      end
      
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
      real(WP):: colat,long
      integer:: vertex,k, domain
      
      ! domain this reference vertex is in
      domain = vertexDomain(vertex)
      
      !check neighbors
      isBoundary = .false.
      do k=1, cellNeighbors(vertex,0)
        if(vertexDomain(cellNeighbors(vertex,k)).ne.domain)then
          isBoundary = .true.
          return
        endif              
      enddo
      
      return
      end
      
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
      integer:: vertex,domain,domainRef,n,k,index,neighbors(nprocesses-1),getDomain
      logical:: isnew
      external:: getDomain
      
      !init
      neighbors(:)=-1
      domainRef=getDomain(vertex)
      
      !find neighbor domain
      index=0
      do n=1,cellNeighbors(vertex,0)
        domain=getDomain(cellNeighbors(vertex,n))
        if( domain .ne. domainRef) then
          ! fill into neighbor if new
          isnew=.true.
          do k=1,size(neighbors(:))
            if( neighbors(k) .eq. domain) then
              isnew=.false.
              exit
            endif
          enddo
          if( isnew ) then
            index=index+1
            neighbors(index)=domain
          endif
        endif 
      enddo      
      
      end    
      
      
!-----------------------------------------------------------------------
      subroutine findDomainNeighbors()
!-----------------------------------------------------------------------
! find domains which are next to the referenced one
!
! returns: domainNeighbors array (each element contains the domain number, -1 for no neighbor)
      use parallel;use griddomain; use verbosity
      implicit none
      integer:: domain,neighbors(nprocesses-1),index,k,ierror
      
      !check
      if( nprocesses .eq. 1) return
      
      !init
      allocate(domainNeighbors(0:nprocesses-1,nprocesses-1),stat=ierror)
      if( ierror .gt. 0 ) then
        print*,'error in allocating domainNeighbors array'
        call stopProgram( 'abort - findDomainNeighbors    ')
      endif
      domainNeighbors(:,:)=-1
      
      !find neighbor domains
      do domain=0,nprocesses-1
        index=0
        do k=0,nprocesses-1
          !check if domain has entries of boundary vertices in this neighbor domain
          if( boundaries(domain,k,1) .ne. 0 ) then 
            ! first element in boundaries array is a valid vertex index, so add this neighbor
            index=index+1
            if(index .gt. nprocesses-1) then
              print*,rank,':exceeding number of neighbor domains',domain
              call stopProgram( 'abort - findDomainNeighbors    ')
            endif
            ! add this as neighbor
            domainNeighbors(domain,index)=k
          endif
        enddo
      enddo
      
      !debug
      if(DEBUG) then
        if( MASTER)then
          do k=0,nprocesses-1
            print*,'domain neighbors (',k,'/.):',domainNeighbors(k,:)
          enddo
        endif      
      endif      
      end
      

