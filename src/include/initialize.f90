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
      subroutine initialize()
!-----------------------------------------------------------------------     
! initializes program parameters
      use propagationStartup;use parallel;use cells;use phaseVelocityMap;use displacements 
      use traveltime; use griddomain;use phaseBlockData;use loop;use deltaSecondLocation
      use filterType;use verbosity
      implicit none
      
      ! process parallelization by MPI
      call initializeMPI()

      ! parameters
      call initializeParameters()
      
      ! mesh
      call initializeMesh()
      
      ! time
      call initializeWorld()
          
      ! prescribe source
      call initializeSource()              

      ! console output         
      if( MASTER .and. VERBOSE) then
        print*
        print*,'initialization o.k.'
        print*
        
        ! benchmark
        benchAllEnd = MPI_WTIME()      
        print*
        print*,'running time: ',int((benchAllEnd-benchAllStart)/60.0),'min ',mod((benchAllEnd-benchAllStart),60.0),'sec'
        print*
      endif      
                  
      end

!-----------------------------------------------------------------------
      subroutine initializeMPI()
!-----------------------------------------------------------------------     
! initializes program parallelization
      use propagationStartup;use parallel
      implicit none
      integer:: ierror
      
      ! parallelization
      call MPI_Init(ierror)
      if( ierror .ne. 0) then
        print*,'error starting mpi program.'
        call MPI_ABORT(MPI_COMM_WORLD,ierror)
      endif
      
      ! get number of processors
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocesses,ierror)
      if( ierror .ne. 0) stop "MPI_COMM_SIZE failed"      
      
      ! determine idproc, the processor's id
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      if( ierror .ne. 0) stop "MPI_COMM_RANK failed"      

      ! determine master process
      if( rank .eq. 0 ) then 
        MASTER=.true.
      else
        MASTER=.false.
      endif
      
      ! determine representation
      if( WP .eq. 4 ) then
        MPI_CUSTOM = MPI_REAL
      else
        MPI_CUSTOM = MPI_DOUBLE_PRECISION
      endif
      
      ! benchmark
      if( MASTER) benchAllStart = MPI_WTIME()      

      end

!-----------------------------------------------------------------------     
      subroutine initializeParameters()
!-----------------------------------------------------------------------     
! initializes program parallelization
      use propagationStartup;use parallel;use filterType;use verbosity
      use loop; use phaseBlockData; use deltaSecondLocation; use adjointVariables;use cells
      implicit none

      ! initialize grid level
      subdivisions = -1

      ! console output
      if( MASTER .and. VERBOSE) then
        if( Adjoint_InversionProgram ) print*,'heterogeneousInversion - computes inversion matrices by'
        if( Adjoint_Program ) then
          print*,'adjointMethod - computes kernel values'
        else
          if( Phaseshift_Program ) then
            print*,'phaseshift - computes kernel values'
          else
            print*,'propagation - membrane wave simulation'
          endif
        endif
      endif
  
      ! console output
      if( MASTER .and. VERBOSE) print*,'-----------------------------------------------------------------------'  
      
      ! reads model input file 'Parameter_Input' and synchronizes with all processes
      ! read in parameters from file 
      if( MASTER ) then 
        call readParameters()
      endif

      ! adjoint method initialization
      if( Adjoint_Program ) then
        DELTA           = .false. ! no delta scatterers needed
        MOVEDELTA       = .false. ! no looping over delta locations
        SECONDDELTA     = .false. ! no second delta scatterers        
        manyReceivers   = .false. ! no many receivers supported        
        kernelIteration = .false. ! no many kernels
        rotate_frame    = .false. ! no rotation of heterogeneous phase maps
        if( manyKernels ) then
          manyKernels=.false.   ! iteration of kernels must be done differently, no manykernels routine supported 
          kernelIteration=.true. 
        endif        
        
        ! console output
        if( MASTER .and. VERBOSE ) then
          print*
          print*,'adjoint calculation (done without delta scatterers)'
        endif
        
!        ! fix original source location (to 1,0,0/...) for inversion (due to rotation of frame)
!        if( Adjoint_InversionProgram ) then
!          sourceLat = 0.0_WP
!          sourceLon = 0.0_WP
!        endif
      endif
      
      ! synchronize with processes
      call syncParameters(rank,nprocesses)
      
      ! parallel simulation needs processes
      if( nprocesses .lt. 2 ) PARALLELSEISMO = .false.
      
      !set subdivisions & maximum numbers of triangles and vertices
      MaxTriangles = 20*(4**(subdivisions+1) - 1)
      MaxVertices  = 30*(4**subdivisions) + 2

      ! initialize for filtering
      call determineFilterParameters()

      ! determine time step size and averageCellDistance
      call determineTimeStep()
             
      ! console output         
      if( MASTER ) then
        print*
        print*,'number of processes:',nprocesses
        if( .not. FIXED_SOURCEPARAMETER ) print*,'wave parameters: factors ',FACTOR_TIME,FACTOR_WIDTH
        if(DELTA) print*,'delta phase velocity increment:', deltaPerturbation
        print*,'number of subdivisions used:', subdivisions
        if( HETEROGENEOUS ) then
          print*,'heterogeneous phase map used'
        else
          print*,'homogeneous phase map', cphaseRef
        endif
        print*,'-----------------------------------------------------------------------'  
      endif      
      
      ! set additional filter verbosity
      if( MASTER .and. VERBOSE ) then
        beVerbose=.true.
      else
        beVerbose=.false.    
      endif
                  
      end


!-----------------------------------------------------------------------     
      subroutine initializeMesh()
!-----------------------------------------------------------------------     
! initializes the mesh, sets up the heterogeneous phase shift map and parallelizes the mesh into domains
      use propagationStartup;use parallel;use cells;use phaseVelocityMap;use displacements 
      use traveltime; use griddomain;use phaseBlockData;use loop;use deltaSecondLocation
      use filterType;use verbosity; use adjointVariables
      implicit none
      real(WP):: lat,lon,distance
      integer:: ierror,i
      
      !allocate arrays   
      allocate(vertices(MaxVertices,3),cellNeighbors(MaxVertices,0:6),cellFace(MaxVertices,0:6),cellCorners(MaxTriangles,3),stat=ierror )
      if( ierror .ne. 0 ) call stopProgram('error in allocating arrays for cell face,..   ')

      if( Station_Correction ) then
        allocate( cellTriangleFace(MaxTriangles,3),stat=ierror)
        if( ierror .ne. 0 ) call stopProgram('error allocating cellTriangleFace    ')
      endif
            
      !read and share initial grid points
      call syncInitialData()      

      ! wait until all processes reached this point
      call MPI_Barrier( MPI_COMM_WORLD, ierror )
      if( ierror .ne. 0) call stopProgram('initializeMesh() - MPI_Barrier failed    ')      
      
      ! find source & receiver vertex
      desiredSourceLat=sourceLat
      desiredSourceLon=sourceLon
      desiredReceiverLat=receiverLat
      desiredReceiverLon=receiverLon
      call findVertex(sourceLat,sourceLon,sourceVertex)      
      call findVertex(receiverLat,receiverLon,receiverVertex)

      ! find triangle where receiver lies
      if( Station_Correction ) then
        call getClosestTriangle(desiredReceiverlat,desiredReceiverlon,interpolation_triangleIndex,&
                              interpolation_distances,interpolation_corners,interpolation_triangleLengths)      
        ! triangle array no more needed
        !deallocate(cellTriangleFace)
        
        ! synchronize with all the other processes
        !call syncInterpolationData()
      endif


      ! store as origin
      originSourceVertex=sourceVertex
      originReceiverVertex=receiverVertex                      
      call getSphericalCoord_Lat(sourceVertex,sourceLat,sourceLon)
      call getSphericalCoord_Lat(receiverVertex,receiverLat,receiverLon)
      
      
      if(DELTA) then
        !pick vector on unit sphere      
        call findVertex(deltaLat,deltaLon,deltaVertex)        
        
        !second delta
        if( SECONDDELTA ) call placeSecondDelta(deltaLat,deltaLon,deltaSecondLat,&
                            deltaSecondLon,deltaSecondVertex,deltaSecondDistance)        
      endif
      
      ! console output
      if( MASTER .and. VERBOSE) then
        print*
        print*,'source(lat/lon) desired:',desiredSourceLat,desiredSourceLon
        print*,'                got:',sourceLat,sourceLon
        print*,'                index:',sourceVertex
        print*,'receiver(lat/lon) desired:',desiredReceiverLat,desiredReceiverLon
        print*,'                  got:',receiverLat,receiverLon
        print*,'                  index:',receiverVertex   
        call greatCircleDistance(vertices(sourceVertex,:),vertices(receiverVertex,:),distance)     
        print*,'distance source-receiver:',distance*180.0/PI
        print*
        if( DELTA) then
          print*,'delta phase map ', DELTA, (cphaseRef + deltaPerturbation)
          print*,'delta(lat/lon) desired:',deltaLat,deltaLon
          call getSphericalCoord_Lat(deltaVertex,lat,lon)
          print*,'               got:',lat,lon
          print*,'               index:',deltaVertex     
          if( SECONDDELTA ) then
            print*,'second delta(lat/lon):',deltaSecondLat,deltaSecondLon
            call getSphericalCoord_Lat(deltaSecondVertex,lat,lon)
            print*,'               got:',lat,lon            
            print*,'               index:',deltaSecondVertex             
          endif             
        endif
        if( Station_Correction ) then
          print*,'interpolation receiver station:'
          print*,'        triangle:',interpolation_triangleIndex
          print*,'        corners :',(interpolation_corners(i),i=1,3)
          print*,'        side lengths (degrees):',(interpolation_triangleLengths(i)*180.0/PI,i=1,3)          
          print*,'        receiverdistance (degr.):',(interpolation_distances(i)*180.0/PI,i=1,3)
          print*
        endif
        print*        
      endif

      ! allocate arrays for displacement datas  
      allocate(displacement(numVertices),displacement_old(numVertices),newdisplacement(numVertices), stat=ierror )
      if( ierror .ne. 0 ) call stopProgram('error in allocating displacement arrays   ')

      !allocate new arrays for precalculated cell attributes  
      allocate( cellAreas(numVertices),cellEdgesLength(numVertices,0:6),&
              cellCenterDistances(numVertices,0:6),phaseVelocitySquare(numVertices), &
              stat=ierror )
      if( ierror .ne. 0 ) call stopProgram('error in allocating arrays for cell area,..   ')

      ! initialize displacements and arrays
      displacement_old(:) = 0.0_WP
      displacement(:)     = 0.0_WP
      newdisplacement(:)  = 0.0_WP
      cellAreas(:)        = 0.0_WP
      cellEdgesLength(:,:)     = 0.0_WP
      cellCenterDistances(:,:) = 0.0_WP
      phaseVelocitySquare(:)   = 0.0_WP

      ! read and share precalculated arrays and sync them
      if( PRECALCULATED_CELLS ) call syncPrecalculated()

      ! initialize geometry of phase map and parallel domains
      call setupFrame()

      ! free arrays if possible
      if( PRECALCULATED_CELLS ) then 
        deallocate(cellFace,cellCorners)
      else
        deallocate(cellAreas,cellEdgesLength,cellCenterDistances)
      endif

!      if( Adjoint_Program .and. ADJOINT_ONTHEFLY ) then
!        ! allocate arrays for displacement datas  
!        allocate(backwardDisplacement(numVertices),backwardDisplacement_old(numVertices),backwardNewdisplacement(numVertices), stat=ierror )
!        if( ierror .ne. 0 ) call stopProgram('error in allocating displacement arrays   ')
!        ! initialize
!        backwardDisplacement_old(:)=0.0_WP
!        backwardDisplacement(:)=0.0_WP
!        backwardNewdisplacement(:)=0.0_WP        
!      endif            
      end


!-----------------------------------------------------------------------     
      subroutine initializeWorld()
!-----------------------------------------------------------------------     
! initializes phase map, time, source
      use propagationStartup;use parallel;use cells;use phaseVelocityMap;use displacements 
      use traveltime; use griddomain;use phaseBlockData;use loop;use deltaSecondLocation
      use filterType;use verbosity; use adjointVariables
      implicit none
      real(WP):: long,colat
      integer:: ierror
      
      ! console output
      if( MASTER .and. VERBOSE) then
        print*
        print*,'world parameters:'
      endif
            
      ! mostly only squared time step is used for iteration
      dt2 = dt*dt
            
      ! time window
      call determineSimulationsteps()

      ! check for sampling rate and wave period (nyquist frequency limitation )
      if( bw_waveperiod .lt. dt*2 .or. bw_waveperiod*cphaseRef .lt. averageCellDistance ) then
        print*,'sampling rate is too low compared to wave period length!'
        print*,'    wave period  : ',bw_waveperiod
        print*,'    wave length  : ',bw_waveperiod*cphaseRef
        print*,'    sampling rate: ',dt
        print*,'    cell length  : ',averageCellDistance
        call stopProgram('abort initializeWorld() - sampling too low   ')      
      endif
      
      !console output
      if( MASTER .and. VERBOSE) then
        print*,'    time step dt:' ,dt
        print*,'    time range:', firsttimestep*dt, lasttimestep*dt
        print*,'    time iteration steps:',numofTimeSteps
        print*,'    phase velocity:',cphaseRef        
        if(DELTA)print*,'    delta heterogeneity: radius=',DELTARADIUS
      endif


      ! adjoint method has one forward and one backward simulation, so start allocating those wavefields    
      if( Adjoint_Program ) then            
        ! allocate arrays
        call allocateAdjointArrays()        
      endif
      end      

    
!-----------------------------------------------------------------------     
      subroutine initializeSource()
!-----------------------------------------------------------------------     
! prescribe initial source
      use propagationStartup;use parallel;use cells;use phaseVelocityMap;use displacements 
      use traveltime; use griddomain;use phaseBlockData;use loop;use deltaSecondLocation
      use filterType;use verbosity; use adjointVariables
      implicit none
      integer:: n,vertex,timestep,index,ierror     
      real(WP):: forceTerm2Source,forceTermExact,time

      ! source parameters
      if( .not. FIXED_SOURCEPARAMETER ) then
        ! calculate time parameter
        TimeParameterSigma=FACTOR_TIME*dt
      
        ! calculate width parameter (see Carl Tape Thesis , eq. 3.21 )
        WidthParameterMu=FACTOR_WIDTH*averageCellDistance/(EARTHRADIUS*2.0*sqrt(-2.0*log(0.05)))
      endif

      ! console output
      if( MASTER .and. VERBOSE) then
        print*
        print*,'source parameters:'
        print*,'    time parameter sigma:',TimeParameterSigma
        print*,'    width parameter mu:',WidthParameterMu        
      endif
      
      !debug
      if(DEBUG) print*,'    allocating forceTermPrescribed',rank,numDomainVertices,numofTimeSteps
      
      ! allocate force array
      if( allocated(forceTermPrescribed) ) deallocate(forceTermPrescribed)
      if( .not. allocated(forceTermPrescribed) ) then
        allocate(forceTermPrescribed(numDomainVertices,numofTimeSteps),stat=ierror)
        if( ierror .gt. 0 ) then
          print*,'    cannot allocate prescribedForce array. ',rank
          print*,'    using calculation of source on the fly...',rank
          PRESCRIBEDSOURCE=.false.
          if( FILTERINITIALSOURCE ) then
            print*,'    using file to store source array.',rank
            sourceOnFile = .true.
            open(sourceFileID,file=datadirectory(1:len_trim(datadirectory))//'forceTermPrescribed.bin',&
                 access='direct',form='unformatted',recl=4)
            call stopProgram('not yet implemented ...   ')
            !print*,'    cannot filter source. ',rank
            !call stopProgram('unable filtering source     ')
          endif
          !call stopProgram('error in allocating prescribed force array    ')
          return        
        endif
        
        !init (critical: processes may stop here)
        forceTermPrescribed(:,:)=0.0_WP
    
        if(DEBUG) then
          ! wait until all processes reached this point    
          call MPI_Barrier( MPI_COMM_WORLD, ierror )
          if( ierror .ne. 0) call stopProgram('abort - MPI_Barrier kernels failed    ')      
        endif        
      endif
      
      ! prescribe the source
      if( MASTER .and. VERBOSE) print*,'    prescribing source ...'
      index=0
      do timestep=firsttimestep,lasttimestep
        ! model time
        time = timestep*dt      
        index=index+1        

        !debug
        !if(DEBUG) print*,'    force time:',rank,time
        
        ! prescribe source for each vertex
        do n=1,numDomainVertices
          ! get cell vertex
          if( PARALLELSEISMO ) then
            vertex=domainVertices(n)
          else
            vertex=n
          endif          
          
          if( Station_Correction ) then
            forceTermPrescribed(n,index)=forceTermExact(vertex,time,desiredSourceLat,desiredSourceLon)          
          else          
            forceTermPrescribed(n,index)=forceTerm2Source(vertex,time,sourceVertex)          
          endif
        enddo   
      enddo
          
      ! filter source
      if( FILTERINITIALSOURCE) call filterSource()
                                          
      end


!-----------------------------------------------------------------------     
      subroutine filterSource()
!-----------------------------------------------------------------------     
! filter initial source
      use propagationStartup;use parallel;use cells;use phaseVelocityMap;use displacements 
      use traveltime; use griddomain;use phaseBlockData;use loop;use deltaSecondLocation
      use filterType;use verbosity; use adjointVariables
      implicit none
      integer:: vertex,timestep,index,i,n,ierror     
      character*8::vertexstr
      real(WP)::seismo(2,numofTimeSteps)
      real(WP)::min,max
      
      !console output
      if( MASTER .and. VERBOSE) then
        print*
        print*,'    filtering source...'
      endif
            
      ! initialize forces time
      index=0
      do timestep=firsttimestep,lasttimestep
        index=index+1
        seismo(1,index)=timestep*dt
      enddo
      
      ! filter initial source 
      do n=1,numDomainVertices
        ! get cell vertex
        if( PARALLELSEISMO ) then
          vertex=domainVertices(n)
        else
          vertex=n
        endif          

        ! get corresponding source 
        seismo(2,:)=forceTermPrescribed(n,:)

        ! file output for sourceVertex
        if( vertex .eq. sourceVertex .and. beVerbose ) then
          write(vertexstr,'(i8.8)') vertex
          open(22,file=datadirectory(1:len_trim(datadirectory))//'source2Prescribed.'//vertexstr//'.dat',iostat=ierror)
          if( ierror .ne. 0 ) call stopProgram('error opening file: '//datadirectory(1:len_trim(datadirectory))//'source2Prescribed.'//vertexstr//'.dat   ')
          do i=1,numofTimeSteps
            write(22,*) seismo(1,i),forceTermPrescribed(n,i)
          enddo
          close(22)       
        endif
        
        ! filter only when there is displacement
        min=0.0
        max=0.0
        do i=1,numofTimeSteps
          if( seismo(2,i) .lt. min ) min = seismo(2,i)
          if( seismo(2,i) .gt. max ) max = seismo(2,i)          
        enddo        
        
        ! filter
        if( (abs(max)+abs(min)) .gt. 1e-4 ) then
          !print*,'    filter',rank,vertex,max,min
          call dofilterSeismogram(seismo,numofTimeSteps)
        endif
          
        ! retain as new source
        forceTermPrescribed(n,:)=seismo(2,:)          
          
        ! start source at zero time for adjoint methods            
        if( Adjoint_Program .and. ADJOINT_STARTATZERO ) then
          if( MASTER .and. VERBOSE .and. n .eq. 1) print*,'starting source for adjoint method at time: zero'
          do i=1,numofTimeSteps
            if( seismo(1,i) .lt. 0.0 ) then            
              forceTermPrescribed(n,i)=0.0_WP            
            endif
          enddo
        endif        
        
        ! file output
        if( vertex .eq. sourceVertex .and. beVerbose ) then
          write(vertexstr,'(i8.8)') vertex
          open(22,file=datadirectory(1:len_trim(datadirectory))//'source2Prescribed.filtered.'//vertexstr//'.dat',iostat=ierror)
          if( ierror .ne. 0 ) call stopProgram('error opening file: '//datadirectory(1:len_trim(datadirectory))//'source2Prescribed.filtered'//vertexstr//'.dat   ')
          do i=1,numofTimeSteps
            write(22,*) seismo(1,i),forceTermPrescribed(n,i)
          enddo
          close(22)       
        endif
        
        if( DEBUG ) then
          if( mod(n,10000) .eq. 0 ) print*,'   filtered',rank,n
        endif
      enddo            
      end

!-----------------------------------------------------------------------
      subroutine setupFrame()
!-----------------------------------------------------------------------
! setup of the phase map and source/receiver/delta locations
      use propagationStartup; use phaseVelocityMap; use parallel; use phaseBlockData
      use cells; use adjointVariables; use verbosity; use griddomain
      implicit none
      real(WP):: distance,vtmp(3),lat,lon
      integer:: i,ierror

      if( HETEROGENEOUS ) then
        !allocate phaseMap array
        allocate(phaseMap(numVertices), stat=ierror)
        if( ierror .ne. 0) then
          print*,'error allocating phaseMap array'
          call stopProgram( 'abort - setupFrame   ')
        endif
           
        ! read phase map (e.g. Love 150 s) (note: is rotated for non-adjoint simulations)
        if( MASTER ) then
          call constructPhasedata()
        endif  
          
        ! synchronize between master and slave processes  
        !print*,rank,'synchronizing phase map...'
        call syncPhaseMap(rank,nprocesses)        
        
        if( rotate_frame ) then
          ! rotated frame has the source sitting at equator (1,0,0)
          ! rotate source & receiver (only for 1 receiver station setup)
          call getRotatedIndex(sourceVertex)
          call getRotatedIndex(receiverVertex)
            
          ! output
          if( MASTER) then
            call getSphericalCoord_Lat(sourceVertex,lat,lon)
            print*, '  rotated source vertex:',sourceVertex,lat,lon
            call getSphericalCoord_Lat(receiverVertex,lat,lon)
            print*, '  rotated receiver vertex:',receiverVertex,lat,lon 
            call greatCircleDistance(vertices(sourceVertex,:),vertices(receiverVertex,:),distance)    
            print*,'  distance source-receiver:',distance*180.0/PI
            print*                      
          endif
        endif
      endif ! HETEROGENEOUS

        
      ! parallelize grid or delta phase maps (when using different processors for a single simulation) 
      ! or determine where the deltaVertex lies (parallel different simulations)
      if( PARALLELSEISMO ) then
        ! construct grid domains for parallelization of single simulation
        call constructParallelDomains()
        
        ! consider number of vertices for finite-difference iteration scheme
        numDomainVertices=size(domainVertices(:)) ! take only vertices in this domain
      else            
        ! set for each process a new delta location
        call parallelizeDeltaLocations()
        
        ! consider number of vertices for finite-difference iteration scheme
        numDomainVertices=numVertices   ! propagate through all vertices
      endif

      ! precalculate the phase velocities for all grid points
      call constructPhaseVelocitySquare()      

      ! phase map no more needed
      if( HETEROGENEOUS ) deallocate(phaseMap)      

!      ! output phase speed model to file
!      if( MASTER ) then
!        ! output
!        open(10,file=datadirectory(1:len_trim(datadirectory))//'tmpPhaseMap.dat')
!        do i=1, numVertices
!          call getSphericalCoord_Lat(i,lat,lon)
!          write(10,*) real(lon),real(lat),sqrt(phaseVelocitySquare(i))
!        enddo
!        close(10)
!      endif



      ! console output
      if( MASTER .and. VERBOSE .and. PARALLELSEISMO) then
        print*,'master'
        print*,'    infos: vertex from',1,'to',numDomainVertices        
      endif      
      ! debug
      if(DEBUG .and. .not. MASTER .and. PARALLELSEISMO) then
        print*,'slave', rank, ' of ',nprocesses
        print*,'    infos:',numVertices, numNeighbors,numCorners,numFaces
        print*,'    infos:vertex from',1,'to',numDomainVertices
        print*      
      endif            
      end
      
!-----------------------------------------------------------------------
      subroutine parallelizeDeltaLocations()
!-----------------------------------------------------------------------
! places delta location according to the number of available processes
!
! return: new common deltaLat,deltaLon,deltaVertex for delta location
      use loop;use propagationStartup;use parallel; use deltaSecondLocation
      implicit none
      
      ! check
      if( .not. DELTA ) return
      
      ! increment latitude corresponding for this process
      if( (deltaLat - (latitudeEnd - (nprocesses-1)*deltaMoveIncrement)) .gt. deltaMoveIncrement/2 ) then
        print*,'delta location too close to end'
        call stopProgram( 'abort - parallelizeDeltaLocations    ')
      endif

      ! next process gets higher delta latitude
      deltaLat=deltaLat+rank*deltaMoveIncrement
      
      ! get nearest vertex for new delta location
      call findVertex(deltaLat,deltaLon,deltaVertex)
      
      ! for 2 delta scatterers
      if( SECONDDELTA) call placeSecondDelta(deltaLat,deltaLon,deltaSecondLat,deltaSecondLon,deltaSecondVertex,deltaSecondDistance)
      end      
          
      
!-----------------------------------------------------------------------
      subroutine getRotatedIndex(n)
!-----------------------------------------------------------------------
! returns the rotated cell index of the rotated frame source/receiver to (1,0,0)/(0,1,0)/(0,0,1)
! (needs original source and receiver vertex)
!
! input:
!       n       - cell index in original frame
!
! returns: rotated cell index
      use propagationStartup; use cells
      implicit none
      integer:: n
      real(WP):: Vsource(3),Vreceiver(3),vtmp(3),rot(3,3),lat,lon
      
      ! determine rotation matrix from source/receiver to (1,0,0)/... frame
      Vsource(:)= vertices(originSourceVertex,:)
      Vreceiver(:)= vertices(originReceiverVertex,:)
      call getInverseRotationMatrix(Vsource,Vreceiver,rot)
      
      ! determine new cell index
      vtmp(:)= vertices(n,:)
      call rotateVector(rot,vtmp,vtmp)
      call getSphericalCoordinates(vtmp,lat,lon)
      call findVertex(lat,lon,n)
      !print*,'  rotated location',n,lat,lon
      end
      
      
!-----------------------------------------------------------------------
      subroutine constructPhasedata()
!-----------------------------------------------------------------------
! calculates for each vertex of the spherical grid a corresponding phase velocity
! which is taken from a pixel map 
! notice: the heterogeneous phase map is rotated for the non-adjoint simulations! 
!           (with source/receiver vertices lying on the equator) 
!
! returns: phaseMap() array
      use propagationStartup;use phaseVelocityMap;use phaseBlockData
      use cells; use adjointVariables
      implicit none
      integer:: index,i
      real(WP):: lat,lon,phasevelocity,vtmp(3),rot(3,3),Vsource(3),Vreceiver(3)

      ! check if phase reference is same as the one used for the propagation
      if(phaseBlockVelocityReference .ne. cphaseRef) then
        print*,'using a false phase velocity data file for this referenced velocity:',cphaseRef,phaseBlockVelocityReference,phaseBlockFile
        call stopProgram( 'abort - constructPhasedata   ')
      endif
      
      ! check if we read a gsh-file
      phaseBlockFile=trim(phaseBlockFile)
      print*,'reading phase map file-type: '//phaseBlockFile(len_trim(phaseBlockFile)-2:len_trim(phaseBlockFile))
      if( phaseBlockFile(len_trim(phaseBlockFile)-2:len_trim(phaseBlockFile)) .eq. "gsh") then
        ! read in phaseMap() array data
        ! notice: using a heterogeneous phase map always rotates the map for a non-adjoint simulation!
        !             such that source/receiver will lie on the equator
        call readGSHPhasemap(phaseBlockFile)
        return
      endif
            
      !read in phase anomalies in percent
      call readPixelPhasemap(phaseBlockFile)
      
      ! recalculate the absolute phase velocity (uses a reference velocity )
      do i=1, numBlocks
        phasevelocity = phaseBlockVelocityReference+phaseBlock(i)*phaseBlockVelocityReference/100.0
        phaseBlock(i) = phasevelocity
      enddo

      if( rotate_frame ) then
        ! determine rotation matrix from (1,0,0)/... to source/receiver frame
        Vsource(:)= vertices(originSourceVertex,:)
        Vreceiver(:)= vertices(originReceiverVertex,:)
        call getRotationMatrix(Vsource,Vreceiver,rot)
      endif
  
      !print*
      do i=1, numVertices
        ! vector to vertex
        vtmp(:)= vertices(i,:)
      
        ! get rotated vector
        if( rotate_frame ) then
          call rotateVector(rot,vtmp,vtmp)
        endif
        
        ! determine lat/lon
        call getSphericalCoordinates(vtmp,lat,lon)
        
        !get vertex block index
        call determineBlock(real(lat),real(lon),index)
        if( index .lt. 0 .or. index .gt. numBlocks) then
          print*,'   block error:',lat,lon
          print*,'   index violates boundaries:',index, i,numBlocks
          call stopProgram( 'abort - constructPhasedata   ')
        endif
        
        ! set corresponding phase velocity
        phaseMap(i)=phaseBlock(index)
      enddo
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
      use precision; use propagationStartup
      implicit none
      real(WP),parameter:: sqrt2=1.414213562373095_WP
            
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
      
      ! time step modification
      !if( timestepsize .gt. 1.0_WP) timestepsize = 53.0_WP
      !!!!!!!!!!!!!!!!!!!! fixed time step: 
      !!!!!!!!!!!!!!!!!!!!    11.0 for phase velocity L35
      !!!!!!!!!!!!!!!!!!!!    9.4 for phase velocity L300
      !timestepsize = 2.2_WP
      end
            
!-----------------------------------------------------------------------
!      subroutine rotatePhaseMap()
!-----------------------------------------------------------------------
! rotates phase map from frame source/receiver to (1,0,0)/(0,1,0)/(0,0,1)
!
! returns: rotated cell index
!      use propagationStartup;use phaseVelocityMap
!      implicit none
!      double precision,allocatable,dimension(:):: rotPhaseMap
!      integer i,index,ierror
!      double precision Vsource(3),Vreceiver(3),vtmp(3)
!      double precision rot(3,3),lat,lon
!        
!      allocate(rotPhaseMap(numVertices),stat=ierror)
!      if( ierror .gt. 0 ) then
!        print*,'error in allocating rotPhaseMap array'
!        stop 'abort - rotatePhaseMap'
!      endif
! 
!      ! determine rotation matrix from source/receiver to (1,0,0)/... frame
!      Vsource(:)=vertices(originSourceVertex,:)
!      Vreceiver(:)=vertices(originReceiverVertex,:)      
!      call getInverseRotationMatrix(Vsource,Vreceiver,rot)
!                         
!      ! create rotated phase map
!      do i=1,numVertices
!        index=i
!        vtmp(:)=vertices(i,:)
!        call rotateVector(rot,vtmp,vtmp)
!        call getSphericalCoordinates(vtmp,lat,lon)
!        call findVertex(lat,lon,index)
!
!        !call getRotatedIndex(index)
!        !print*,'   index:',i,index
!        !print*,'   phase value:',phaseMap(i)
!        !print*,'   rotated value:',rotPhaseMap(index)
!        if( index .lt. 0 .or. index .gt. numVertices) then
!          print*,'   rotated index violates boundaries:',index, i
!          stop 'abort - rotatePhaseMap'
!        endif
!        
!        if( rotPhaseMap(index) .lt. 0.000000001) then
!          rotPhaseMap(index)=phaseMap(i)
!        else
!          print*,'phase map rotation: ',i,' already mapped to', index, rotPhaseMap(index)
!          rotPhaseMap(index)=phaseMap(i)
!          !stop
!        endif
!      enddo
!      
!      ! copy it over
!      phaseMap=rotPhaseMap
!      
!      ! output
!      open(10,file='tmpRotatedPhaseMap.dat')
!      do i=1, numVertices
!        write(10,*) i, phaseMap(i)
!      enddo
!      close(10)
!      end
      
      
!-----------------------------------------------------------------------
      subroutine initializeKernels()
!-----------------------------------------------------------------------     
! initializes kernels and receivers
      use propagationStartup;use parallel;use cells;use phaseVelocityMap;use displacements; 
      use traveltime; use griddomain; use phaseBlockData
      use loop;use deltaSecondLocation;use filterType;use verbosity
      implicit none
      integer:: m,rounded,ierror
      real(WP):: exact,distance
       
      ! place all receiver stations
      if( manyReceivers ) then
        if( manyKernels ) then
          ! determine how many kernels
          numofKernels=int(kernelEndDistance-kernelStartDistance+1)          
          if( numofKernels .le. 0) then
            print*,'kernels cannot be found correctly:',kernelStartDistance,kernelEndDistance
            call stopProgram( 'abort - initializeKernels()   ')            
          endif
          ! allocate kernel arrays
          allocate( kernelsReceivers(numofKernels,numofReceivers), &
            kernelsReceiversSeismogram(numofKernels,numofReceivers+1,numofTimeSteps), &
            kernelsReceiversSeismogramRef(numofKernels,numofReceivers+1,numofTimeSteps), &
            stat=ierror)
          if( ierror .ne. 0 ) call stopProgram('cannot allocate receivers array   ')

          ! console output
          if( VERBOSE .and. MASTER) then
            print*,'arrays for receivers allocated successfully',size(receiversSeismogram),size(kernelsReceiversSeismogramRef)
            print*,'    epicentral distances between:',kernelStartDistance,kernelEndDistance
            print*,'    number of kernels:',numofKernels
          endif
          
          ! setup the stations
          do m=1,numofKernels          
            currentKernel=m
            ! determine kernels station locations (with respect to a source at north pole)
            receiverLat=90.0 - (kernelStartDistance + m-1)
            if( receiverLat .lt. -90.0) then
              print*,'kernels to far away from source (north pole):',receiverLat
              call stopProgram( 'abort - initializeKernels()   ')              
            endif
            receiverLon=0.0 ! initial first station
            call setupReceivers()            
          enddo
          ! console output
          if( VERBOSE .and. MASTER ) then
            print*,'    set up kernels for latitude: ',90-kernelStartDistance,90-kernelEndDistance            
          endif          
        else
          call setupReceivers()
        endif
      endif
                     
      ! prepare kernel data files
      if( manyKernels) then
        do m=1,numofKernels
          currentKernel=m
          call createKernelFiles()
        enddo      
      else
        call createKernelFiles()
      endif
      
      ! output epicentral distances
      if( MASTER ) then
        open(10,file=datadirectory(1:len_trim(datadirectory))//'tmpEpiDistances.dat')
        write(10,*)'# epicentral distances'
        write(10,*)'# rounded  exactDistance (in degrees)'
        if( VERBOSE ) print*,'epicentral distances:'
        if( manyKernels ) then  
          do m=1,numofKernels
            rounded=int(kernelStartDistance+m-1)
            call greatCircleDistance(vertices(sourceVertex,:),vertices(kernelsReceivers(m,1),:),distance)
            exact=distance*180.0/PI                  
            write(10,*) rounded, exact
            if( VERBOSE ) print*,'    ',rounded,exact            
          enddo
        else
          call greatCircleDistance(vertices(sourceVertex,:),vertices(receiverVertex,:),distance)
          rounded=nint(distance)*180.0/PI
          call greatCircleDistance(vertices(sourceVertex,:),vertices(receiverVertex,:),distance)
          exact=distance*180.0/PI
          write(10,*) rounded, exact          
          if( VERBOSE ) print*,'    ',rounded,exact                      
        endif
        close(10)
        if( VERBOSE ) print*
      endif
      
      ! output     
      if( MASTER .and. VERBOSE ) then
        print*
        print*,'kernel initialization o.k.'
        print*
      endif               
      
      end

!-----------------------------------------------------------------------
      subroutine setupReceivers()
!-----------------------------------------------------------------------   
! places all receiver stations around the same latitude
!
! returns: new receivers/... arrays
      use propagationStartup; use parallel; use verbosity
      implicit none
      integer:: i,ierror
      real(WP):: lat,lon
      character*128:: datafile
      character*3:: kernelstr
      
      ! useless for heterogeneous case
      if( HETEROGENEOUS ) call stopProgram('multiple receiver stations and heterogeneous phase map not applicable!    ')
      ! check number of receivers
      if( numofReceivers .le. 0 ) then
        print*,'number of receiver stations invalid!',numofReceivers
        call stopProgram('setupReceivers() - abort    ')
      endif
    
      ! store for all receivers the displacements and add one seismogram to store the corresponding time   
      if( .not. manyKernels ) then
        ! allocate receiver arrays
        allocate( receivers(numofReceivers), &
                receiversSeismogram(numofReceivers+1,numofTimeSteps),&
                receiversSeismogramRef(numofReceivers+1,numofTimeSteps),&
                stat=ierror)      
        if( ierror .ne. 0 ) call stopProgram('cannot allocate receivers array   ')
        !console output
        if( MASTER .and. VERBOSE ) print*,'arrays for receivers allocated successfully',size(receiversSeismogram),size(receiversSeismogramRef)
      endif
      
      ! assumes the source to be on north pole, the receiver station placed along the same latitude
      if( importKernelsReceivers ) then      
        if( manyReceivers ) then
          ! read in values
          if( manyKernels ) then
            write(kernelstr,'(i3.3)')int(kernelStartDistance+currentKernel-1)
            datafile=datadirectory(1:len_trim(datadirectory))//'tmpReceiverStations'//kernelstr//'.dat'
            open(10,file=datafile)
            do i=1, numofReceivers
              read(10,*) lat,lon,kernelsReceivers(currentKernel,i)
            enddo
            close(10)
          else
            if( VERBOSE .and. MASTER) print*,'importing kernels only for manyKernels supported'
            importKernelsReceivers=.false.
          endif
          
          ! output
          if( VERBOSE .and. MASTER .and. currentKernel .eq. 1) then
            print*,'read in receiver stations from file:', datafile
          endif
        else
          if( VERBOSE .and. MASTER) print*,'importing kernels only for manyReceivers supported'
          importKernelsReceivers=.false.
        endif
      endif
      
      ! determine receiver locations
      if( .not. importKernelsReceivers ) then
        ! check if needed
        if( .not. manyReceivers ) return
        
        ! place stations
        do i=1, numofReceivers
          ! vary longitude
          receiverLon=(i-1)*360.0_WP/numofReceivers
          
          ! vary latitude
          !receiverLat=(i-1)*360.0_WP/numofReceivers
          
          ! get corresponding vertex for this receiver station on our grid
          if( manyKernels) then
            call findVertex(receiverLat,receiverLon,kernelsReceivers(currentKernel,i))
          else
            call findVertex(receiverLat,receiverLon,receivers(i))
          endif
        enddo
        
        ! debug file output
        if( MASTER) then
          if( manyKernels) then
            write(kernelstr,'(i3.3)')int(kernelStartDistance+currentKernel-1)
            datafile=datadirectory(1:len_trim(datadirectory))//'tmpReceiverStations'//kernelstr//'.dat'
            open(10,file=datafile)
            do i=1, numofReceivers
              call getSphericalCoord_Lat(kernelsReceivers(currentKernel,i),lat,lon)
              write(10,*) lat,lon,kernelsReceivers(currentKernel,i)
            enddo
            close(10)
          else
            datafile=datadirectory(1:len_trim(datadirectory))//'tmpReceiverStations.dat'
            open(10,file=datafile)
            do i=1, numofReceivers
              call getSphericalCoord_Lat(receivers(i),lat,lon)
              write(10,*) lat,lon,receivers(i)
            enddo
            close(10)        
          endif
        endif              
      endif
      
      end
      
      
!-----------------------------------------------------------------------   
      subroutine createKernelFiles()
!-----------------------------------------------------------------------   
! creates new data files in datadirectory to write kernel values in
!
! returns: new ttkernel.rank***.dat and ttkernel.rot.rank***.dat files
      use propagationStartup; use parallel; use verbosity
      implicit none
      integer:: i,n,ioerror
      character*3:: rankstr,kernelstr
      character*128:: datafile
      
      ! only master process creates the files
      ! (but all processes must have access to these files)    
      if( MASTER) then
        if( VERBOSE .and. currentKernel .le. 1) then
          print*,'creating new output files:'
          print*,'   ',datadirectory(1:len_trim(datadirectory))//'ttkernel.rank***.dat'
          print*,'   ',datadirectory(1:len_trim(datadirectory))//'ttkernel.rot.rank***.dat'
        endif
        
        do i=0,nprocesses-1
          write(rankstr,'(i3.3)') i
          
          ! create new 'ttkernel.dat' and 'ttkernel.rot.dat' files
          do n=1,2
            ! determine filename            
            if( n .eq. 1) then
              if( manyKernels) then            
                write(kernelstr,'(i3.3)') int(kernelStartDistance+currentKernel-1)              
                datafile=datadirectory(1:len_trim(datadirectory))//'ttkernel'//kernelstr//'.rank'//rankstr//'.dat'
              else
                datafile=datadirectory(1:len_trim(datadirectory))//'ttkernel.rank'//rankstr//'.dat'              
              endif
            else
              if( manyKernels) then
                write(kernelstr,'(i3.3)') int(kernelStartDistance+currentKernel-1)              
                datafile=datadirectory(1:len_trim(datadirectory))//'ttkernel'//kernelstr//'.rot.rank'//rankstr//'.dat'          
              else
                datafile=datadirectory(1:len_trim(datadirectory))//'ttkernel.rot.rank'//rankstr//'.dat'                        
              endif
            endif
            datafile=trim(datafile)
            ! delete if existing
            !open(10,file=datafile,iostat=ioerror)
            !if( ioerror .eq. 0 ) close(10,status='DELETE')
            
            ! put two comment lines at beginning
            open(10,file=datafile,iostat=ioerror)
            if( ioerror .ne. 0) then
              ! check if really we can not create the file
              print*,'could not open file: '//datafile
              print*,'  try again...'
              open(10,file=datafile,status='unknown',iostat=ioerror)
              if( ioerror .ne. 0 ) call stopProgram('createKernelFiles() - still not possible. shutting down    ')
              print*,'  successed opening file'
            endif
            ! determine comment line text            
            if( n .eq. 1 ) then
              write(10,*) '# phase shift - sensitivity kernel'
            else
              write(10,*) '# rotated phase shift - sensitivity kernel'
            endif
            
            ! add a values legend
            write(10,'(1x,''# lon lat kernel kernelAnalytic receiverVertex'',''timelag timelagAnayltic vperturbation receiverLat receiverLon'')')
            !write(100,*) '# lon  lat  kernel  kernelAnalytic  receiverVertex  timelag  timelagAnalytic  vperturbation receiverLat  receiverLon'
            close(10)
          enddo
        enddo
      endif
      end

!-----------------------------------------------------------------------
      subroutine prepareCouple(epla,eplo,stla,stlo)
!-----------------------------------------------------------------------
! prepare model for a new simulation with the setup of a new source and receiver location
!
! input:
!   epla,eplo    - event/source latitude and longitude
!   stla,stlo       - station/receiver latitude and longitude
!
! returns: sets new desired locations
      use propagationStartup; use cells
      implicit none
      real:: epla,eplo,stla,stlo
      real(WP):: distance
      
      ! source will be placed
      desiredSourceLat = epla
      desiredSourceLon = eplo
      call setupNewSource( desiredSourceLat, desiredSourceLon )

      ! receivers will be placed
      desiredReceiverLat = stla
      desiredReceiverLon = stlo
      call setupNewReceiver( desiredReceiverLat, desiredReceiverLon )

      ! determine epicentral distance
      call greatCircleDistance(vertices(sourceVertex,:),vertices(receiverVertex,:),distance)

      ! change simulation time according to epicentral distance or antipode 
      call setupNewSimulationtime(distance)      

      ! recalculate the force array
      call initializeSource()             

            
      end

!-----------------------------------------------------------------------
      subroutine setupNewSource(lat,lon)
!-----------------------------------------------------------------------
! setup new receiver location
      use propagationStartup; use parallel; use verbosity
      implicit none
      real(WP):: distance,lat,lon
      
      ! find receiver location on the grid
      call findVertex(lat,lon,sourceVertex)
      call getSphericalCoord_Lat(sourceVertex,sourceLat,sourceLon)
      
      if( MASTER .and. VERBOSE) then
        print*,'source(lat/lon) desired:',lat,lon
        print*,'                  got:',sourceLat,sourceLon
        print*,'                  index:',sourceVertex
      endif 
      
      end      


!-----------------------------------------------------------------------
      subroutine setupNewReceiver(lat,lon)
!-----------------------------------------------------------------------
! setup new receiver location
      use propagationStartup; use parallel; use verbosity; use cells
      implicit none
      real(WP):: lat,lon
      real(WP):: distance
      integer:: i
      
      ! find receiver location on the grid
      call findVertex(lat,lon,receiverVertex)
      call getSphericalCoord_Lat(receiverVertex,receiverLat,receiverLon)

      ! store as origin
      originReceiverVertex=receiverVertex

      ! determine epicentral distance
      call greatCircleDistance(vertices(sourceVertex,:),vertices(receiverVertex,:),distance)
      
      ! find triangle where receiver lies
      if( Station_Correction ) then
        call getClosestTriangle(desiredReceiverlat,desiredReceiverlon,interpolation_triangleIndex,&
                              interpolation_distances,interpolation_corners,interpolation_triangleLengths)      
      endif

      
      if( MASTER .and. VERBOSE) then
        print*,'receiver(lat/lon) desired:',desiredReceiverLat,desiredReceiverLon
        print*,'                  got:',receiverLat,receiverLon
        print*,'                  index:',receiverVertex
        print*,'distance source-receiver:',distance*180.0/PI
        print*
        if( Station_Correction ) then
          print*,'interpolation receiver station:'
          print*,'        triangle:',interpolation_triangleIndex
          print*,'        corners :',(interpolation_corners(i),i=1,3)
          print*,'        side lengths (degrees):',(interpolation_triangleLengths(i)*180.0/PI,i=1,3)          
          print*,'        receiverdistance (degr.):',(interpolation_distances(i)*180.0/PI,i=1,3)
          print*
        endif        
      endif        
      end            

!-----------------------------------------------------------------------
      subroutine setupStation(lat,lon)
!-----------------------------------------------------------------------
! setup new receiver location and adapt simulation time
      use adjointVariables; use propagationStartup; use parallel; use cells; use verbosity
      implicit none
      real(WP):: distance,lat,lon
      character*64,parameter::distanceFile                  = 'tmpEpiDistances.dat'
      integer:: ierror

      ! place new receiver
      call setupNewReceiver(lat,lon)
            
      ! determine epicentral distance
      call greatCircleDistance(vertices(sourceVertex,:),vertices(receiverVertex,:),distance)
      
      if( MASTER .and. VERBOSE) then
        if( Adjoint_Program ) then
          if( kernelIteration ) then
            ! print epicentral distances into file
            open(10,file=datadirectory(1:len_trim(datadirectory))//distanceFile(1:len_trim(distanceFile)),status='old',position='append',iostat=ierror)    
            if( ierror .ne. 0 ) then
              open(10,file=datadirectory(1:len_trim(datadirectory))//distanceFile(1:len_trim(distanceFile)),status='new',iostat=ierror)        
              if( ierror .eq. 0 ) then
                write(10,*) '# epicentral distances'
                write(10,*) '# rounded exactDistance (in degrees)'
              endif
            endif
            if( ierror .eq. 0 ) then
              write(10,*) int(desiredReceiverLon),distance*180.0/PI
              close(10)
            endif        
          endif
        endif
      endif
      
      ! change simulation time according to epicentral distance or antipode 
      call setupNewSimulationtime(distance)      

      end            


!-----------------------------------------------------------------------
      subroutine setupNewSimulationtime(distance)
!-----------------------------------------------------------------------
! change times according to epicentral distance or antipode 
      use propagationStartup; use parallel; use verbosity
      implicit none
      real(WP)::distance,arrivalTime,antipodeTime
      
      ! set last time
      antipodeTime = PI*EARTHRADIUS/cphaseRef
      arrivalTime = distance*EARTHRADIUS/cphaseRef
      if( USEOVERTIME ) then
        LASTTIME = arrivalTime + SIMULATIONOVERTIMEPERCENT*arrivalTime
        if( LASTTIME > antipodeTime ) LASTTIME = antipodetime
      else
        LASTTIME = antipodeTime
      endif      
      if( MASTER .and. VERBOSE) then
        print*,'simulation time will end at:',LASTTIME
      endif            

      ! set number of time steps
      call determineSimulationsteps()      
      
      end
      
      
!-----------------------------------------------------------------------
      subroutine determineSimulationsteps()
!-----------------------------------------------------------------------
! determine the time steping for the simulation
      use propagationStartup
      implicit none
      integer:: ierror
      
      ! stepings    
      firsttimestep = int( FIRSTTIME/dt) !default: -6500/ almost same results for -800
      lasttimestep = int( LASTTIME/dt) !default: 22000
      numofTimeSteps=lasttimestep-firsttimestep+1
      
      
      ! newly allocate seismogram for receiver
      if( allocated(seismogramReceiver) ) deallocate(seismogramReceiver)
      allocate(seismogramReceiver(2,numofTimeSteps),stat=ierror)
      if( ierror .gt. 0 ) call stopProgram('error in allocating seismogram array    ')
      
      ! fft size
      call determineFFTWindowsize()
      
      end

      