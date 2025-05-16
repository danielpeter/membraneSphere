!=====================================================================
!
!       m e m b r a n e S p h e r e
!       --------------------------------------------------
!      (c) 2025
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
  real(WP), external:: syncWtime

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

  ! wait until all processes reached this point
  call syncProcesses()

  ! console output
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'initialization ok'

    ! benchmark
    benchAllEnd = syncWtime()
    print *,'initialization running time: ',int((benchAllEnd-benchAllStart)/60.0),'min ', &
                      mod((benchAllEnd-benchAllStart),60.0),'sec'
    print *
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine initializeMPI()
!-----------------------------------------------------------------------
! initializes program parallelization
  use propagationStartup; use parallel
  implicit none
  ! local parameters
  real(WP), external:: syncWtime

  ! parallelization
  call syncInitMPI()

  ! get number of processors
  call syncGetNumberProcesses(nprocesses)

  ! determine idproc, the processor's id
  call syncGetRank(myrank)

  ! determine main process
  if (myrank == 0) then
    MAIN_PROCESS = .true.
  else
    MAIN_PROCESS = .false.
  endif

  ! determine representation
  call syncDetermineMPICustom(WP,MPI_CUSTOM)

  ! benchmark
  if (MAIN_PROCESS) benchAllStart = syncWtime()

  end subroutine


!-----------------------------------------------------------------------
  subroutine initializeParameters()
!-----------------------------------------------------------------------
! initializes program parallelization
  use propagationStartup; use parallel; use filterType; use verbosity
  use loop; use phaseBlockData; use deltaSecondLocation; use adjointVariables;
  use cells; use precisions
  implicit none

  ! initialize grid level
  subdivisions = -1

  ! console output
  if (MAIN_PROCESS .and. VERBOSE) then
    if (Adjoint_InversionProgram) print *,'heterogeneousInversion - computes inversion matrices by'
    if (Adjoint_Program) then
      print *,'adjointMethod - computes kernel values'
    else if (Phaseshift_Program) then
      print *,'phaseshift - computes kernel values'
    else if (HetPhaseshift_Program) then
      print *,'phaseshift measurements - computed on the membrane'
    else
      print *,'propagation - membrane wave simulation'
    endif
    !endif
  endif

  ! console output
  if (MAIN_PROCESS .and. VERBOSE) print *,'-----------------------------------------------------------------------'

  ! reads model input file 'Parameter_Input' and synchronizes with all processes
  ! read in parameters from file
  if (MAIN_PROCESS) then
    call readParameters()
  endif
  ! synchronize with processes
  call syncParameters(myrank,nprocesses)

  ! adjoint method initialization
  if (Adjoint_Program) then
    DELTA           = .false. ! no delta scatterers needed
    MOVEDELTA       = .false. ! no looping over delta locations
    SECONDDELTA     = .false. ! no second delta scatterers
    manyReceivers   = .false. ! no many receivers supported
    kernelIteration = .false. ! no many kernels
    if (manyKernels) then
      manyKernels = .false.   ! iteration of kernels must be done differently, no manykernels routine supported
      kernelIteration = .true.
    endif

    ! check that no rotation of heterogeneous phase maps
    if (ROTATE_FRAME) call stopProgram('adjoint method needs ROTATE_FRAME set to .false.    ')

    ! console output
    if (MAIN_PROCESS .and. VERBOSE) then
      print *
      print *,'adjoint calculation (done without delta scatterers)'
    endif

!        ! fix original source location (to 1,0,0/...) for inversion (due to rotation of frame)
!        if (Adjoint_InversionProgram) then
!          sourceLat = 0.0_WP
!          sourceLon = 0.0_WP
!        endif
  endif

  ! parallel simulation needs processes
  if (nprocesses < 2) PARALLELSEISMO = .false.

  ! initialize for filtering
  call determineFilterParameters()

  ! determine time step size and averageCellDistance
  call determineTimeStep()

  ! console output
  if (MAIN_PROCESS) then
    print *
    print *,'number of processes            : ',nprocesses
    if (.not. FIXED_SOURCEPARAMETER) &
      print *,'wave parameters                : ',THETA_WIDTH,' theta_width factor'
    if (DELTA) &
      print *,'delta phase velocity increment : ', deltaPerturbation
    print *,'number of subdivisions used    : ', subdivisions
    if (HETEROGENEOUS) then
      print *,'phase map                      : heterogeneous map used'
    else
      print *,'phase map                      : ', cphaseRef,' (homogeneous map)'
    endif
    print *,'-----------------------------------------------------------------------'
    print *
  endif

  ! set additional filter verbosity
  if (MAIN_PROCESS .and. VERBOSE) then
    beVerbose = .true.
  else
    beVerbose = .false.
  endif

  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'parameters ok'
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine initializeMesh()
!-----------------------------------------------------------------------
! initializes the mesh, sets up the heterogeneous phase shift map and parallelizes the mesh into domains
  use propagationStartup;use parallel;use cells;use phaseVelocityMap;use displacements
  use traveltime; use griddomain;use phaseBlockData;use loop;use deltaSecondLocation
  use filterType;use verbosity; use adjointVariables
  implicit none
  ! local parameters
  real(WP):: lat,lon,distance
  real(WP):: vectorA(3),vectorB(3)
  integer:: i

  ! allocate grid arrays
  call allocateMesh()

  ! read and share initial grid points
  if (MAIN_PROCESS) call readData(VERBOSE)

  call syncInitialData()

  ! wait until all processes reached this point
  call syncProcesses()

  ! find source & receiver vertex
  desiredSourceLat = sourceLat
  desiredSourceLon = sourceLon
  desiredReceiverLat = receiverLat
  desiredReceiverLon = receiverLon
  call findVertex(sourceLat,sourceLon,sourceVertex)
  call findVertex(receiverLat,receiverLon,receiverVertex)

  ! find triangle where receiver lies
  if (STATION_CORRECTION) then
    call getClosestTriangle(desiredReceiverlat,desiredReceiverlon,interpolation_triangleIndex, &
                            interpolation_distances,interpolation_corners,interpolation_triangleLengths)
    ! triangle array no more needed
    !deallocate(cellTriangleFace)

    ! synchronize with all the other processes
    !call syncInterpolationData()
  endif

  ! store as origin
  originSourceVertex = sourceVertex
  originReceiverVertex = receiverVertex
  call getSphericalCoord_Lat(sourceVertex,sourceLat,sourceLon)
  call getSphericalCoord_Lat(receiverVertex,receiverLat,receiverLon)

  if (DELTA) then
    !pick vector on unit sphere
    call findVertex(deltaLat,deltaLon,deltaVertex)

    !second delta
    if (SECONDDELTA) call placeSecondDelta(deltaLat,deltaLon,deltaSecondLat, &
                                           deltaSecondLon,deltaSecondVertex,deltaSecondDistance)
  endif

  ! console output
  if (MAIN_PROCESS .and. VERBOSE) then
    print *
    print *,'source(lat/lon) desired   : ',desiredSourceLat,desiredSourceLon
    print *,'                got       : ',sourceLat,sourceLon
    print *,'                index     : ',sourceVertex
    print *,'receiver(lat/lon) desired : ',desiredReceiverLat,desiredReceiverLon
    print *,'                  got     : ',receiverLat,receiverLon
    print *,'                  index   : ',receiverVertex
    ! distance
    vectorA(:) = vertices(sourceVertex,:)
    vectorB(:) = vertices(receiverVertex,:)
    call greatCircleDistance(vectorA,vectorB,distance)
    print *,'distance source-receiver  : ',distance*180.0/PI
    print *
    if (DELTA) then
      print *,'delta phase map       : ', DELTA,'/',(cphaseRef + deltaPerturbation)
      print *,'delta(lat/lon) desired: ',deltaLat,deltaLon
      call getSphericalCoord_Lat(deltaVertex,lat,lon)
      print *,'                   got: ',lat,lon
      print *,'                 index: ',deltaVertex
      if (SECONDDELTA) then
        print *,'second delta(lat/lon): ',deltaSecondLat,deltaSecondLon
        call getSphericalCoord_Lat(deltaSecondVertex,lat,lon)
        print *,'                  got: ',lat,lon
        print *,'                index: ',deltaSecondVertex
      endif
    endif
    if (STATION_CORRECTION) then
      print *,'interpolation receiver station:'
      print *,'  triangle              : ',interpolation_triangleIndex
      print *,'  corners               : ',(interpolation_corners(i),i=1,3)
      print *,'  side lengths (deg)    : ',(interpolation_triangleLengths(i)*180.0/PI,i=1,3)
      print *,'  receiverdistance (deg): ',(interpolation_distances(i)*180.0/PI,i=1,3)
      print *
    endif
    print *
  endif

  ! allocate arrays for displacement datas
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'initializing data arrays...'
  endif
  call allocateData()

  ! read and share precalculated arrays and sync them
  if (PRECALCULATED_CELLS) call syncPrecalculated()

  ! initialize geometry of phase map and parallel domains
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'initializing phase map...'
  endif
  call setupFrame()

  ! free arrays if possible
  if (PRECALCULATED_CELLS) then
    deallocate(cellFace,cellCorners)
  else
    deallocate(cellAreas,cellEdgesLength,cellCenterDistances)
  endif

!      if (Adjoint_Program .and. ADJOINT_ONTHEFLY) then
!        ! allocate arrays for displacement datas
!        allocate(backwardDisplacement(numVertices),backwardDisplacement_old(numVertices), &
!                      backwardNewdisplacement(numVertices), stat=ier )
!        if (ier /= 0) call stopProgram('error in allocating displacement arrays   ')
!        ! initialize
!        backwardDisplacement_old(:) = 0.0_WP
!        backwardDisplacement(:) = 0.0_WP
!        backwardNewdisplacement(:) = 0.0_WP
!      endif
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'  meshing ok'
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine initializeWorld()
!-----------------------------------------------------------------------
! initializes phase map, time, source
  use propagationStartup;use parallel;use cells;use phaseVelocityMap;use displacements
  use traveltime; use griddomain;use phaseBlockData;use loop;use deltaSecondLocation
  use filterType;use verbosity; use adjointVariables
  implicit none
  ! local parameters
  real(WP):: distance,wavelength,gridpoints_per_wavelength
  integer:: iorbit

  ! console output
  if (MAIN_PROCESS .and. VERBOSE) then
    print *
    print *,'world parameters:'
  endif

  ! mostly only squared time step is used for iteration
  dt2 = dt*dt

  ! determines simulation time
  if (Adjoint_Program .and. ADJOINT_ANTIPODE_TIME) then
    if (MAIN_PROCESS .and. VERBOSE) &
      print *,'  using antipode time as last time'
    call determineOrbitDistance(distance,iorbit)
    call setupNewSimulationtime(distance,iorbit)
  endif

  ! sets number of time steps
  call determineSimulationsteps()

  ! reference wavelength (in km)
  wavelength = bw_waveperiod * cphaseRef
  ! grid sampling
  gridpoints_per_wavelength = wavelength / averageCellDistance

  ! console output
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'  simulation:'
    print *,'    reference phase velocity          : ',cphaseRef,'km/s'
    print *,'              wave period             : ',bw_waveperiod,'s'
    print *,'              wave length             : ',wavelength,'km'
    print *,'    average cell center distance      : ',averageCellDistance,'km'
    print *,'            gridpoints per wavelength : ',gridpoints_per_wavelength
    print *,'    time step dt                      : ',dt
    print *,'    time range                        : ',firsttimestep*dt,lasttimestep*dt
    print *,'    time iteration steps              : ',numofTimeSteps
    if (DELTA) &
      print *,'    delta heterogeneity - radius      : ',DELTARADIUS
  endif

  ! check for sampling rate and wave period (nyquist frequency limitation )
  if (bw_waveperiod < dt*2 .or. bw_waveperiod*cphaseRef < averageCellDistance) then
    print *
    print *,'Warning: sampling is very low compared to wave period or wave length!'
    print *,'         time sampling rate is ',dt
    print *,'         average cell distance is ',averageCellDistance
    print *
    !call stopProgram('abort initializeWorld() - sampling too low   ')
  endif

  ! adjoint method has one forward and one backward simulation,
  ! allocates and initializes those wavefields
  if (Adjoint_Program) then
    call initializeAdjointArrays()
  endif

  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'  world ok'
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine initializeSource()
!-----------------------------------------------------------------------
! prescribe initial source
  use propagationStartup; use parallel; use griddomain
  use cells, only: numDomainVertices
  use verbosity
  implicit none
  ! local parameters
  integer:: n,vertex,timestep,index,ier,jrec
  real(WP),external:: forceTerm2Source,forceTermExact
  real(WP):: time,force,arraysize
  real(WP):: forcetrace(numofTimeSteps)
  character(len=3):: rankstr
  real(WP), parameter:: EPS = 1.0e-5
  real(WP),parameter:: RAM_LIMIT = 2000   ! in MB ; e.g.  2 GB equal 2000 MB
  character(len=128),parameter:: temporaryDir   = '/scratch/tmp/'  ! choose something appropriate for your system

  ! console output
  if (MAIN_PROCESS .and. VERBOSE) then
    print *
    print *,'initializing source...'
  endif

  ! source parameters
  call determineSourceParameters()

  ! pre-calculates the force terms
  if (PRESCRIBED_SOURCE) then
    ! prescribe the source
    if (MAIN_PROCESS .and. VERBOSE) print *,'  prescribing source...'

    ! check if to reallocate force array
    if (allocated(forceTermPrescribed)) then
      if (size(forceTermPrescribed(:,1)) /= numDomainVertices .or. &
         size(forceTermPrescribed(1,:)) /= numofTimeSteps) then
        deallocate(forceTermPrescribed)
      endif
    endif

    ! reallocates array
    if (.not. allocated(forceTermPrescribed)) then
      arraysize = numDomainVertices*numofTimeSteps*WP/1024./1024.
      ! console output
      if (MAIN_PROCESS .and. VERBOSE) then
        print *,'  allocating prescribedForce array: '
        print *,'      vertices            : ',numDomainVertices
        print *,'      steps               : ',numofTimeSteps
        print *,'      size                : ',arraysize,'(MB)'
      endif

      ! allocate array
      if (arraysize > RAM_LIMIT .and. FILTER_INITIALSOURCE) then
        if (MAIN_PROCESS .and. VERBOSE) print *,'    using file to store source array.'
        sourceOnFile = .true.
      else
        ! allocates memory
        allocate(forceTermPrescribed(numDomainVertices,numofTimeSteps),stat=ier)
        if (ier /= 0) then
          print *,'    cannot allocate prescribedForce array: ',myrank
          ! change prescribed flag if necessary
          if (FILTER_INITIALSOURCE) then
            if (MAIN_PROCESS .and. VERBOSE) print *,'    using file to store source array.'
            sourceOnFile = .true.
          else
            ! abort
            print *,'Error: exceeding memory limit, invalid prescribe source parameter'
            print *,'       Please set PRESCRIBED_SOURCE = .false. to use calculation of source on the fly...'
            call stopProgram('Exceeding memory limit for PRESCRIBED_SOURCE    ')
          endif
        endif
      endif
      call syncFlag(myrank,nprocesses,sourceOnFile)

      ! storage on file
      if (sourceOnFile) then
        if (allocated(forceTermPrescribed)) deallocate(forceTermPrescribed)

        ! names file
        write(rankstr,'(i3.3)') myrank
        print *,'      prescribed source file: ',trim(datadirectory)//'forceTermPrescribed'//rankstr//'.bin'
        ! opens file to store data for each process
        open(sourceFileID,file=trim(datadirectory)//'forceTermPrescribed'//rankstr//'.bin', &
             access='direct',form='unformatted',recl=WP)

        ! storage in scratch dir
        !open(sourceFileID,file=trim(temporaryDir)//'forceTermPrescribed'//rankstr//'.bin', &
        !     access='direct',form='unformatted',recl=WP)
      endif
    endif ! not allocated

    !init (critical: processes may stop here)
    if (.not. sourceOnFile) then
      !print *,'    intializing prescribedForce array: rank',myrank
      forceTermPrescribed(:,:) = 0.0_WP
    else
      forcetrace(:) = 0.0_WP
    endif

    ! prescribe source for each vertex
    do n = 1,numDomainVertices
      ! get cell vertex
      if (PARALLELSEISMO) then
        vertex = domainVertices(n)
      else
        vertex = n
      endif

      index = 0
      do timestep = firsttimestep,lasttimestep
        ! model time
        time = timestep*dt
        index = index+1
        if (STATION_CORRECTION) then
          force = forceTermExact(vertex,time,desiredSourceLat,desiredSourceLon)
        else
          force = forceTerm2Source(vertex,time,sourceVertex)
        endif
        ! store values
        if (.not. sourceOnFile) then
          forceTermPrescribed(n,index) = force
        else
          forcetrace(index) = force
        endif

        ! stops if source vanishes
        if (time > 500.0 .and. abs(force) < EPS) exit
      enddo

      ! store on file if needed
      if (sourceOnFile) then
        do index = 1,numofTimeSteps
          jrec = (n-1)*numofTimeSteps + index
          write(sourceFileID,rec=jrec) forcetrace(index)
        enddo
      endif
    enddo

    ! filter source
    if (FILTER_INITIALSOURCE) call filterSource()
  endif

  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'  source ok'
    print *
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine filterSource()
!-----------------------------------------------------------------------
! filter initial source
  use precisions
  use propagationStartup; use parallel; use cells; use phaseVelocityMap; use displacements
  use traveltime; use griddomain; use phaseBlockData; use loop; use deltaSecondLocation
  use filterType; use verbosity; use adjointVariables
  implicit none
  ! local parameters
  integer:: vertex,timestep,index,i,n,ier,jrec
  character(len=9):: vertexstr
  real(WP):: seismo(2,numofTimeSteps)
  real(WP):: originalT,originalf
  real(WP)::amp_min,amp_max,amp_min_all,amp_max_all,amp_min_all_filtered,amp_max_all_filtered
  ! amplitude threshold to determine whether source is zero or gets filtered
  real(WP),parameter :: THRESHOLD_SOURCE_AMPLITUDE = 1e-4

  ! checks if anything to do
  if (.not. FILTER_INITIALSOURCE) return

  ! filter around a specified (different) period
  if (FILTER_AT_NEWPERIOD) then
    originalf = bw_frequency
    originalT = bw_waveperiod

    bw_waveperiod = BW_NEWPERIOD

    ! get (upper) half-bandwidth frequency for filtering
    if (BW_FIXFREQUENCY) then
      bw_frequency = BW_HALFFREQUENCY
    else
      bw_frequency = BW_PERCENT/((BW_PERCENT+1.0)*bw_waveperiod)
    endif
  endif

  ! console output
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'  filtering source...'
    if (FILTER_AT_NEWPERIOD) then
      print *,'    using filter at new period'
      print *,'    period              : ',bw_waveperiod
    endif
    if (BW_FIXFREQUENCY) then
      print *,'    using fixed halfband-frequency'
    else
      print *,'    filter percent      : ',BW_PERCENT
      print *,'    filter center period: ',bw_waveperiod
    endif
    print *,'    filter frequency    : ',bw_frequency
    print *,'           as period    : ',1.0/bw_frequency
  endif

  ! initialize forces time
  index = 0
  do timestep = firsttimestep,lasttimestep
    index = index+1
    seismo(1,index) = timestep*dt
  enddo

  ! stats
  amp_max_all = 0.0_WP
  amp_min_all = 0.0_WP
  amp_max_all_filtered = 0.0_WP
  amp_min_all_filtered = 0.0_WP

  ! filter initial source
  do n = 1,numDomainVertices
    ! get cell vertex
    if (PARALLELSEISMO) then
      vertex = domainVertices(n)
    else
      vertex = n
    endif

    ! get corresponding source
    if (.not. sourceOnFile) then
      seismo(2,:) = forceTermPrescribed(n,:)
    else
      do i = 1,numofTimeSteps
        jrec = (n-1)*numofTimeSteps + i
        read(sourceFileID,rec=jrec) seismo(2,i)
      enddo
    endif

    ! file output for sourceVertex
    if (vertex == sourceVertex .and. fileOutput .and. .not. sourceOnFile) then
      write(vertexstr,'(i9.9)') vertex
      print *,'  writing to file: ',trim(datadirectory)//'source2Prescribed.'//vertexstr//'.dat'
      open(IOUT,file=trim(datadirectory)//'source2Prescribed.'//vertexstr//'.dat',iostat=ier)
      if (ier /= 0) &
        call stopProgram('error opening file: '//trim(datadirectory)//'source2Prescribed.'//vertexstr//'.dat   ')
      do i = 1,numofTimeSteps
        write(IOUT,*) seismo(1,i),forceTermPrescribed(n,i)
      enddo
      close(IOUT)
    endif

    ! stats
    ! get +/- amplitudes
    amp_min = minval(seismo(2,:))
    amp_max = maxval(seismo(2,:))
    if (amp_min < amp_min_all) amp_min_all = amp_min
    if (amp_max > amp_max_all) amp_max_all = amp_max

    ! filter only when there is some displacement in the trace
    if ( (abs(amp_max)+abs(amp_min)) > THRESHOLD_SOURCE_AMPLITUDE) then
      !print *,'    filter',myrank,vertex,max,min
      call dofilterSeismogram(seismo,numofTimeSteps)
      ! apply hanning window to smooth ends
      call taperSeismogram(seismo,numofTimeSteps,numofTimeSteps,beVerbose)
    endif

    ! start source at zero time for adjoint methods
    if (Adjoint_Program .and. ADJOINT_STARTATZERO) then
      if (MAIN_PROCESS .and. VERBOSE .and. n == 1) print *,'starting source for adjoint method at time: zero'
      do i = 1,numofTimeSteps
        if (seismo(1,i) < 0.0) then
          seismo(2,i) = 0.0_WP
        endif
      enddo
    endif

    ! retain as new source
    if (.not. sourceOnFile) then
      forceTermPrescribed(n,:) = seismo(2,:)
    else
      do i = 1,numofTimeSteps
        jrec = (n-1)*numofTimeSteps + i
        write(sourceFileID,rec=jrec) seismo(2,i)
      enddo
    endif

    ! stats for filtered source
    ! get +/- amplitudes
    amp_min = minval(seismo(2,:))
    amp_max = maxval(seismo(2,:))
    if (amp_min < amp_min_all_filtered) amp_min_all_filtered = amp_min
    if (amp_max > amp_max_all_filtered) amp_max_all_filtered = amp_max

    ! file output
    if (vertex == sourceVertex .and. fileOutput .and. .not. sourceOnFile) then
      write(vertexstr,'(i9.9)') vertex
      print *,'  writing to file: ',trim(datadirectory)//'source2Prescribed.filtered.'//vertexstr//'.dat'
      open(IOUT,file=trim(datadirectory)//'source2Prescribed.filtered.'//vertexstr//'.dat',iostat=ier)
      if (ier /= 0) &
        call stopProgram('error opening file: '//trim(datadirectory)//'source2Prescribed.filtered'//vertexstr//'.dat   ')
      do i = 1,numofTimeSteps
        write(IOUT,*) seismo(1,i),forceTermPrescribed(n,i)
      enddo
      close(IOUT)
    endif
  enddo

  ! get minimum/maximum over all processes
  call syncMin(amp_min_all)
  call syncMax(amp_max_all)
  call syncMin(amp_min_all_filtered)
  call syncMax(amp_max_all_filtered)

  ! console output
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'    original source amplitude min/max : ',amp_min_all,'/',amp_max_all
    print *,'    filtered source amplitude min/max : ',amp_min_all_filtered,'/',amp_max_all_filtered
  endif

  ! filter around a specified (different) period
  if (FILTER_AT_NEWPERIOD) then
    bw_waveperiod = originalT
    bw_frequency = originalf
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine setupFrame()
!-----------------------------------------------------------------------
! sets up of the phase map and source/receiver/delta locations
  use propagationStartup; use phaseVelocityMap; use parallel; use phaseBlockData
  use cells; use adjointVariables; use verbosity; use griddomain
  implicit none
  ! local parameters
  real(WP):: distance,lat,lon
  real(WP):: vectorA(3),vectorB(3)
  integer:: ier

  ! reads in a heterogeneous phase velocity map
  if (HETEROGENEOUS) then
    if (MAIN_PROCESS .and. VERBOSE) then
      print *,'  initializing heterogeneous phase map...'
    endif

    !allocate phaseMap array
    allocate(phaseMap(numVertices), stat=ier)
    if (ier /= 0) then
      print *,'Error: allocating phaseMap array'
      call stopProgram( 'abort - setupFrame   ')
    endif

    ! read phase map (e.g. Love 150 s) (note: is rotated for non-adjoint simulations)
    if (MAIN_PROCESS) then
      call constructPhasedata()
    endif

    ! synchronize between main and secondary processes
    !print *,'rank ',myrank,'synchronizing phase map...'
    call syncPhaseMap()

    if (ROTATE_FRAME) then
      ! rotated frame has the source sitting at equator (1,0,0)
      ! rotate source & receiver (only for 1 receiver station setup)
      call getRotatedIndex(sourceVertex)
      call getRotatedIndex(receiverVertex)

      ! output
      if (MAIN_PROCESS) then
        call getSphericalCoord_Lat(sourceVertex,lat,lon)
        print *, '  rotated source vertex:',sourceVertex,lat,lon
        call getSphericalCoord_Lat(receiverVertex,lat,lon)
        print *, '  rotated receiver vertex:',receiverVertex,lat,lon
        ! distance
        vectorA(:) = vertices(sourceVertex,:)
        vectorB(:) = vertices(receiverVertex,:)
        call greatCircleDistance(vectorA,vectorB,distance)
        print *,'  distance source-receiver:',distance*180.0/PI
        print *
      endif
    endif
  endif ! HETEROGENEOUS


  ! parallelize grid or delta phase maps (when using different processors for a single simulation)
  ! or determine where the deltaVertex lies (parallel different simulations)
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'  initializing parallelized maps...'
  endif
  if (PARALLELSEISMO) then
    ! construct grid domains for parallelization of single simulation
    call constructParallelDomains()

    ! consider number of vertices for finite-difference iteration scheme
    numDomainVertices = size(domainVertices(:)) ! take only vertices in this domain
  else
    ! set for each process a new delta location
    if (DELTA) call parallelizeDeltaLocations()

    ! consider number of vertices for finite-difference iteration scheme
    numDomainVertices = numVertices   ! propagate through all vertices
  endif

  ! precalculate the phase velocities for all grid points
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'  initializing phase velocities...'
  endif
  call constructPhaseVelocitySquare()

  ! phase map no more needed
  if (HETEROGENEOUS) deallocate(phaseMap)

  ! console output
  if (MAIN_PROCESS .and. VERBOSE .and. PARALLELSEISMO) then
    print *,'  main process info: vertex from',1,'to',numDomainVertices
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine parallelizeDeltaLocations()
!-----------------------------------------------------------------------
! places delta location according to the number of available processes
!
! return: new common deltaLat,deltaLon,deltaVertex for delta location
  use loop;use propagationStartup;use parallel; use deltaSecondLocation
  implicit none

  ! increment latitude corresponding for this process
  if ( (deltaLat - (latitudeEnd - (nprocesses-1)*deltaMoveIncrement)) > deltaMoveIncrement/2) then
    print *,'Error: delta location too close to end'
    call stopProgram( 'abort - parallelizeDeltaLocations    ')
  endif

  ! next process gets higher delta latitude
  deltaLat = deltaLat + myrank * deltaMoveIncrement

  ! get nearest vertex for new delta location
  call findVertex(deltaLat,deltaLon,deltaVertex)

  ! for 2 delta scatterers
  if (SECONDDELTA) call placeSecondDelta(deltaLat,deltaLon,deltaSecondLat, &
                                         deltaSecondLon,deltaSecondVertex,deltaSecondDistance)

  end subroutine


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
  integer,intent(inout):: n
  ! local parameters
  real(WP):: Vsource(3),Vreceiver(3),vtmp(3),rot(3,3),lat,lon

  ! determine rotation matrix from source/receiver to (1,0,0)/... frame
  Vsource(:) = vertices(originSourceVertex,:)
  Vreceiver(:) = vertices(originReceiverVertex,:)
  call getInverseRotationMatrix(Vsource,Vreceiver,rot)

  ! determine new cell index
  vtmp(:) = vertices(n,:)
  call rotateVector(rot,vtmp,vtmp)
  call getSphericalCoordinates(vtmp,lat,lon)
  call findVertex(lat,lon,n)
  !print *,'  rotated location',n,lat,lon

  end subroutine


!-----------------------------------------------------------------------
  subroutine constructPhasedata()
!-----------------------------------------------------------------------
! calculates for each vertex of the spherical grid a corresponding phase velocity
! which is taken from a pixel map
!
! notice: ROTATE_FRAME - the heterogeneous phase map is rotated for the non-adjoint simulations!
!                        (with source/receiver vertices lying on the equator)
!
! returns: phaseMap() array
  use propagationStartup;use phaseVelocityMap;use phaseBlockData
  use cells; use adjointVariables; use verbosity; use parallel
  implicit none
  ! local parameters
  integer:: index,i
  real(WP):: lat,lon,phasevelocity,vtmp(3),rot(3,3),Vsource(3),Vreceiver(3)

  ! only main process constructs phase map (will be broadcast afterwards)
  if (.not. MAIN_PROCESS) return

  ! phase velocity map
  if (DO_CHECKERBOARD) then
    ! checkerboard
    call makeCheckerboard()
  else
    ! heterogeneous phase map
    ! check if phase reference is same as the one used for the propagation
    if (phaseBlockVelocityReference /= cphaseRef) then
      print *,'Error: using a false phase velocity data file for this referenced velocity:', &
                    cphaseRef,phaseBlockVelocityReference,trim(phaseBlockFile)
      call stopProgram('Abort - invalid phase reference in constructPhasedata()   ')
    endif

    ! check if we read a gsh-file
    phaseBlockFile = trim(phaseBlockFile)
    print *,'  reading phase map file-type: '//&
          phaseBlockFile(len_trim(phaseBlockFile)-2:len_trim(phaseBlockFile))

    if (phaseBlockFile(len_trim(phaseBlockFile)-2:len_trim(phaseBlockFile)) == "gsh") then
      ! GSH phase map
      ! read in phaseMap() array data
      ! notice: using a heterogeneous phase map always rotates the map for a non-adjoint simulation!
      !             such that source/receiver will lie on the equator
      if (phaseBlockFile(len_trim(phaseBlockFile)-6:len_trim(phaseBlockFile)) == "abs.gsh") then
        call readGSHPhasemap(phaseBlockFile,.false.)
      else
        call readGSHPhasemap(phaseBlockFile,.true.)
      endif
    else
      ! Pixel phase map
      ! read in phase anomalies in percent
      call readPixelPhasemap(phaseBlockFile)

      ! recalculate the absolute phase velocity (uses a reference velocity )
      do i = 1, numBlocks
        phasevelocity = phaseBlockVelocityReference+phaseBlock(i)*phaseBlockVelocityReference/100.0
        phaseBlock(i) = phasevelocity
      enddo

      if (ROTATE_FRAME) then
        ! determine rotation matrix from (1,0,0)/... to source/receiver frame
        Vsource(:) = vertices(originSourceVertex,:)
        Vreceiver(:) = vertices(originReceiverVertex,:)
        call getRotationMatrix(Vsource,Vreceiver,rot)
      endif

      !print *
      do i = 1, numVertices
        ! vector to vertex
        vtmp(:) = vertices(i,:)

        ! get rotated vector
        if (ROTATE_FRAME) then
          call rotateVector(rot,vtmp,vtmp)
        endif

        ! determine lat/lon
        call getSphericalCoordinates(vtmp,lat,lon)

        ! get vertex block index
        call determineBlock(real(lat),real(lon),index)
        if (index < 0 .or. index > numBlocks) then
          print *,'Error: block error:',lat,lon
          print *,'       index violates boundaries:',index, i,numBlocks
          call stopProgram( 'abort - constructPhasedata   ')
        endif

        ! set corresponding phase velocity
        phaseMap(i) = phaseBlock(index)
      enddo
    endif
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine initializeKernels()
!-----------------------------------------------------------------------
! initializes kernels and receivers
  use propagationStartup;use parallel;use cells;use phaseVelocityMap;use displacements;
  use traveltime; use griddomain; use phaseBlockData
  use loop;use deltaSecondLocation;use filterType;use verbosity
  implicit none
  ! local parameters
  integer:: m,rounded,ier
  real(WP):: exact,distance
  real(WP):: vectorA(3),vectorB(3)

  ! place all receiver stations
  if (manyReceivers) then
    if (manyKernels) then
      ! determine how many kernels
      numofKernels=int(kernelEndDistance-kernelStartDistance+1)
      if (numofKernels <= 0) then
        print *,'Error: kernels cannot be found correctly:',kernelStartDistance,kernelEndDistance
        call stopProgram( 'abort - initializeKernels()   ')
      endif
      ! allocate kernel arrays
      allocate( kernelsReceivers(numofKernels,numofReceivers), &
        kernelsReceiversSeismogram(numofKernels,numofReceivers+1,numofTimeSteps), &
        kernelsReceiversSeismogramRef(numofKernels,numofReceivers+1,numofTimeSteps), &
        stat=ier)
      if (ier /= 0) call stopProgram('cannot allocate receivers array   ')

      ! console output
      if (VERBOSE .and. MAIN_PROCESS) then
        print *,'  arrays for receivers allocated successfully',size(receiversSeismogram), &
                          size(kernelsReceiversSeismogramRef)
        print *,'  epicentral distances between:',kernelStartDistance,kernelEndDistance
        print *,'  number of kernels:',numofKernels
      endif

      ! setup the stations
      do m = 1,numofKernels
        currentKernel = m
        ! determine kernels station locations (with respect to a source at north pole)
        receiverLat=90.0 - (kernelStartDistance + m-1)
        if (receiverLat < -90.0) then
          print *,'Error: kernels to far away from source (north pole):',receiverLat
          call stopProgram( 'abort - initializeKernels()   ')
        endif
        receiverLon = 0.0 ! initial first station
        call setupReceivers()
      enddo
      ! console output
      if (VERBOSE .and. MAIN_PROCESS) then
        print *,'  set up kernels for latitude: ',90-kernelStartDistance,90-kernelEndDistance
      endif
    else
      call setupReceivers()
    endif
  endif

  ! prepare kernel data files
  if (manyKernels) then
    do m = 1,numofKernels
      currentKernel = m
      call createKernelFiles()
    enddo
  else
    call createKernelFiles()
  endif

  ! output epicentral distances
  if (MAIN_PROCESS) then
    open(IOUT,file=trim(datadirectory)//'tmpEpiDistances.dat')
    write(IOUT,*)'# epicentral distances'
    write(IOUT,*)'# rounded  exactDistance (in degrees)'
    if (VERBOSE) print *,'epicentral distances:'
    if (manyKernels) then
      do m = 1,numofKernels
        rounded = int(kernelStartDistance+m-1)
        ! distance
        vectorA(:) = vertices(sourceVertex,:)
        vectorB(:) = vertices(kernelsReceivers(m,1),:)
        call greatCircleDistance(vectorA,vectorB,distance)
        exact = distance*180.0/PI
        write(IOUT,*) rounded, exact
        if (VERBOSE) print *,'    ',rounded,exact
      enddo
    else
      ! distance
      vectorA(:) = vertices(sourceVertex,:)
      vectorB(:) = vertices(receiverVertex,:)
      call greatCircleDistance(vectorA,vectorB,distance)
      rounded = nint(distance)*180.0/PI
      exact = distance*180.0/PI
      write(IOUT,*) rounded, exact
      if (VERBOSE) print *,'    ',rounded,exact
    endif
    close(IOUT)
    if (VERBOSE) print *
  endif

  ! output
  if (MAIN_PROCESS .and. VERBOSE) then
    print *
    print *,'kernel initialization ok'
    print *
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine setupReceivers()
!-----------------------------------------------------------------------
! places all receiver stations around the same latitude
!
! returns: new receivers/... arrays
  use propagationStartup; use parallel; use verbosity
  implicit none
  ! local parameters
  integer:: i,ier
  real(WP):: lat,lon
  character(len=128):: datafile
  character(len=3):: kernelstr

  ! useless for heterogeneous case
  if (HETEROGENEOUS) &
    call stopProgram('multiple receiver stations and heterogeneous phase map not applicable    ')
  ! check number of receivers
  if (numofReceivers <= 0) then
    print *,'Error: number of receiver stations invalid!',numofReceivers
    call stopProgram('setupReceivers() - abort    ')
  endif

  ! store for all receivers the displacements and add one seismogram to store the corresponding time
  if (.not. manyKernels) then
    ! allocate receiver arrays
    allocate( receivers(numofReceivers), &
            receiversSeismogram(numofReceivers+1,numofTimeSteps), &
            receiversSeismogramRef(numofReceivers+1,numofTimeSteps), &
            stat=ier)
    if (ier /= 0) call stopProgram('cannot allocate receivers array   ')
    !console output
    if (MAIN_PROCESS .and. VERBOSE) &
      print *,'arrays for receivers allocated successfully',size(receiversSeismogram),size(receiversSeismogramRef)
  endif

  ! assumes the source to be on north pole, the receiver station placed along the same latitude
  if (importKernelsReceivers) then
    if (manyReceivers) then
      ! read in values
      if (manyKernels) then
        write(kernelstr,'(i3.3)')int(kernelStartDistance+currentKernel-1)
        datafile = trim(datadirectory)//'tmpReceiverStations'//kernelstr//'.dat'
        open(IIN,file=trim(datafile))
        do i = 1, numofReceivers
          read(IIN,*) lat,lon,kernelsReceivers(currentKernel,i)
        enddo
        close(IIN)
      else
        if (VERBOSE .and. MAIN_PROCESS) print *,'importing kernels only for manyKernels supported'
        importKernelsReceivers = .false.
      endif

      ! output
      if (VERBOSE .and. MAIN_PROCESS .and. currentKernel == 1) then
        print *,'read in receiver stations from file: ',trim(datafile)
      endif
    else
      if (VERBOSE .and. MAIN_PROCESS) print *,'importing kernels only for manyReceivers supported'
      importKernelsReceivers = .false.
    endif
  endif

  ! determine receiver locations
  if (.not. importKernelsReceivers) then
    ! check if needed
    if (.not. manyReceivers) return

    ! place stations
    do i = 1, numofReceivers
      ! vary longitude
      receiverLon=(i-1)*360.0_WP/numofReceivers

      ! vary latitude
      !receiverLat=(i-1)*360.0_WP/numofReceivers

      ! get corresponding vertex for this receiver station on our grid
      if (manyKernels) then
        call findVertex(receiverLat,receiverLon,kernelsReceivers(currentKernel,i))
      else
        call findVertex(receiverLat,receiverLon,receivers(i))
      endif
    enddo

    ! debug file output
    if (MAIN_PROCESS) then
      if (manyKernels) then
        write(kernelstr,'(i3.3)')int(kernelStartDistance+currentKernel-1)
        datafile = trim(datadirectory)//'tmpReceiverStations'//kernelstr//'.dat'
        open(IOUT,file=trim(datafile))
        do i = 1, numofReceivers
          call getSphericalCoord_Lat(kernelsReceivers(currentKernel,i),lat,lon)
          write(IOUT,*) lat,lon,kernelsReceivers(currentKernel,i)
        enddo
        close(IOUT)
      else
        datafile = trim(datadirectory)//'tmpReceiverStations.dat'
        open(IOUT,file=trim(datafile))
        do i = 1, numofReceivers
          call getSphericalCoord_Lat(receivers(i),lat,lon)
          write(IOUT,*) lat,lon,receivers(i)
        enddo
        close(IOUT)
      endif
    endif
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine createKernelFiles()
!-----------------------------------------------------------------------
! creates new data files in datadirectory to write kernel values in
!
! returns: new ttkernel.rank***.dat and ttkernel.rot.rank***.dat files
  use propagationStartup; use parallel; use verbosity
  implicit none
  ! local parameters
  integer:: i,n,ier
  character(len=3):: rankstr,kernelstr
  character(len=128):: datafile

  ! only main process creates the files
  ! (but all processes must have access to these files)
  if (.not. MAIN_PROCESS) return

  if (VERBOSE .and. currentKernel <= 1) then
    print *,'creating new output files:'
    print *,'   ',trim(datadirectory)//'ttkernel.rank***.dat'
    print *,'   ',trim(datadirectory)//'ttkernel.rot.rank***.dat'
  endif

  do i = 0,nprocesses-1
    write(rankstr,'(i3.3)') i

    ! create new 'ttkernel.dat' and 'ttkernel.rot.dat' files
    do n = 1,2
      ! determine filename
      if (n == 1) then
        if (manyKernels) then
          write(kernelstr,'(i3.3)') int(kernelStartDistance+currentKernel-1)
          datafile = trim(datadirectory)//'ttkernel'//kernelstr//'.rank'//rankstr//'.dat'
        else
          datafile = trim(datadirectory)//'ttkernel.rank'//rankstr//'.dat'
        endif
      else
        if (manyKernels) then
          write(kernelstr,'(i3.3)') int(kernelStartDistance+currentKernel-1)
          datafile = trim(datadirectory)//'ttkernel'//kernelstr//'.rot.rank'//rankstr//'.dat'
        else
          datafile = trim(datadirectory)//'ttkernel.rot.rank'//rankstr//'.dat'
        endif
      endif
      datafile = trim(datafile)

      ! delete if existing
      !open(IOUT,file=datafile,iostat=ier)
      !if (ier == 0) close(IOUT,status='DELETE')

      ! put two comment lines at beginning
      open(IOUT,file=trim(datafile),iostat=ier)
      if (ier /= 0) then
        ! check if really we can not create the file
        print *,'could not open file: ',trim(datafile)
        print *,'  try again...'
        open(IOUT,file=trim(datafile),status='unknown',iostat=ier)
        if (ier /= 0) call stopProgram('createKernelFiles() - still not possible. shutting down    ')
        print *,'  successed opening file'
      endif
      ! determine comment line text
      if (n == 1) then
        write(IOUT,*) '# phase shift - sensitivity kernel'
      else
        write(IOUT,*) '# rotated phase shift - sensitivity kernel'
      endif

      ! add a values legend
      write(IOUT,*) '# lon lat kernel kernelAnalytic receiverVertex timelag ' // &
                    'timelagAnayltic vperturbation receiverLat receiverLon'
      close(IOUT)
    enddo
  enddo

  end subroutine


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
  real,intent(in):: epla,eplo,stla,stlo
  ! local parameters
  real(WP):: distance
  real(WP):: vectorA(3),vectorB(3)
  integer:: iorbit

  ! source will be placed
  desiredSourceLat = epla
  desiredSourceLon = eplo
  call setupNewSource( desiredSourceLat, desiredSourceLon )

  ! receivers will be placed
  desiredReceiverLat = stla
  desiredReceiverLon = stlo
  call setupNewReceiver( desiredReceiverLat, desiredReceiverLon )

  ! determine epicentral distance
  vectorA(:) = vertices(sourceVertex,:)
  vectorB(:) = vertices(receiverVertex,:)
  call greatCircleDistance(vectorA,vectorB,distance)

  ! changes simulation time according to epicentral distance or antipode
  iorbit = 1
  call setupNewSimulationtime(distance,iorbit)

  ! sets number of time steps
  call determineSimulationsteps()

  ! recalculate the force array
  call initializeSource()

  end subroutine


!-----------------------------------------------------------------------
  subroutine setupNewSource(lat,lon)
!-----------------------------------------------------------------------
! setup new receiver location
  use propagationStartup; use parallel; use verbosity
  implicit none
  real(WP),intent(in):: lat,lon

  ! find receiver location on the grid
  call findVertex(lat,lon,sourceVertex)
  call getSphericalCoord_Lat(sourceVertex,sourceLat,sourceLon)

  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'    source(lat/lon) desired  : ',lat,lon
    print *,'                        got  : ',sourceLat,sourceLon
    print *,'                      index  : ',sourceVertex
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine setupNewReceiver(lat,lon)
!-----------------------------------------------------------------------
! setup new receiver location
  use propagationStartup; use parallel; use verbosity; use cells
  implicit none
  real(WP),intent(in):: lat,lon
  ! local parameters
  real(WP):: distance
  real(WP):: vectorA(3),vectorB(3)
  integer:: i

  ! find receiver location on the grid
  call findVertex(lat,lon,receiverVertex)
  call getSphericalCoord_Lat(receiverVertex,receiverLat,receiverLon)

  ! store as origin
  originReceiverVertex = receiverVertex

  ! determine epicentral distance
  vectorA(:) = vertices(sourceVertex,:)
  vectorB(:) = vertices(receiverVertex,:)
  call greatCircleDistance(vectorA,vectorB,distance)

  ! find triangle where receiver lies
  if (STATION_CORRECTION) then
    call getClosestTriangle(desiredReceiverlat,desiredReceiverlon,interpolation_triangleIndex, &
                          interpolation_distances,interpolation_corners,interpolation_triangleLengths)
  endif

  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'    receiver(lat/lon) desired: ',desiredReceiverLat,desiredReceiverLon
    print *,'                          got: ',receiverLat,receiverLon
    print *,'                        index: ',receiverVertex
    print *,'    distance source-receiver : ',distance*180.0/PI
    if (STATION_CORRECTION) then
      print *,'    interpolation receiver station:'
      print *,'      triangle               : ',interpolation_triangleIndex
      print *,'      corners                : ',(interpolation_corners(i),i=1,3)
      print *,'      side lengths (degr)    : ',(interpolation_triangleLengths(i)*180.0/PI,i=1,3)
      print *,'      receiverdistance (degr): ',(interpolation_distances(i)*180.0/PI,i=1,3)
    endif
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine setupStation(lat,lon)
!-----------------------------------------------------------------------
! sets up new receiver location and kernel output file
  use adjointVariables; use propagationStartup; use parallel; use cells; use verbosity
  implicit none
  real(WP),intent(in):: lat,lon
  ! local parameters
  real(WP):: distance
  character(len=64),parameter::distanceFile     = 'Kernel_EpiDistances.dat'
  integer:: ier,iorbit

  ! place new receiver
  call setupNewReceiver(lat,lon)

  ! writes epicentral distance to file
  if (MAIN_PROCESS) then
    if (Adjoint_Program) then
      if (kernelIteration) then
        ! determine epicentral distance
        call determineOrbitDistance(distance,iorbit)

        ! print epicentral distances into file
        open(IOUT,file=trim(datadirectory)//trim(distanceFile), &
             status='old',position='append',iostat=ier)
        if (ier /= 0) then
          ! console output
          if (VERBOSE) then
            print *,'  epicentral distances of kernels in file:'
            print *,'    ',trim(datadirectory)//trim(distanceFile)
            print *
          endif
          ! open as new file
          open(IOUT,file=trim(datadirectory)//trim(distanceFile), &
               status='new',iostat=ier)
          if (ier == 0) then
            write(IOUT,*) '# epicentral distances'
            write(IOUT,*) '# rounded exactDistance (in degrees)'
          endif
        endif
        if (ier == 0) then
          write(IOUT,*) int(desiredReceiverLon),distance*180.0/PI
          close(IOUT)
        endif
      endif
    endif
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine setupNewSimulationtime(distance,iorbit)
!-----------------------------------------------------------------------
! changes simulation time according to epicentral distance or antipode given and
! depending on parameter USE_OVERTIME set by commonModules
!
! input:
!   distance - epicentral distance between source & receiver (can be higher orbit)
!   iorbit   - orbit number
!
! returns: LASTTIME is set accordingly
  use propagationStartup; use parallel; use verbosity
  implicit none
  real(WP),intent(in)::distance
  integer,intent(in):: iorbit
  ! local parameters
  real(WP):: arrivalTime,antipodeTime

  ! set last time
  antipodeTime = (PI + (iorbit-1)*PI)*EARTHRADIUS/cphaseRef
  arrivalTime = distance*EARTHRADIUS/cphaseRef

  if (USE_OVERTIME) then
    ! uses overtime percent instead of antipodal time
    LASTTIME = arrivalTime + OVERTIME_PERCENT*arrivalTime
    if (LASTTIME > antipodeTime) LASTTIME = antipodetime
  else if (WINDOWED_INTEGRATION) then
    if (MAIN_PROCESS .and. VERBOSE) then
      print *,'    orbit: ',iorbit
    endif
    LASTTIME = antipodeTime
  else
    LASTTIME = antipodeTime
  endif

  ! console output
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'    new simulation time will end at: ',LASTTIME
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine determineSimulationsteps()
!-----------------------------------------------------------------------
! determines the time steping for the simulation
  use propagationStartup, only: FIRSTTIME,LASTTIME,dt,firsttimestep,lasttimestep, &
                        numofTimeSteps,seismogramReceiver
  use parallel, only: MAIN_PROCESS
  use verbosity, only: VERBOSE
  use filterType, only: WindowSIZE
  implicit none
  ! local parameters
  integer:: ier,isteps

  ! stepings
  firsttimestep = int(FIRSTTIME/dt) !default: -6500/ almost same results for -800
  lasttimestep  = int(LASTTIME/dt)  !default: 22000

  isteps = lasttimestep - firsttimestep + 1

  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'    number of time steps  : ',isteps
  endif

  if (numofTimeSteps /= isteps) then
    numofTimeSteps = isteps
    ! newly allocate seismogram for receiver
    if (allocated(seismogramReceiver)) then
     if (size(seismogramReceiver(1,:)) /= numofTimeSteps) deallocate(seismogramReceiver)
    endif
    if (.not. allocated(seismogramReceiver)) allocate(seismogramReceiver(2,numofTimeSteps),stat=ier)
    if (ier /= 0) call stopProgram('error in allocating seismogram array    ')

    ! fft size
    call determineFFTWindowsize(numofTimeSteps,WindowSIZE)
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine determineOrbitDistance(distance,iorbit)
!-----------------------------------------------------------------------
! determines the travelled distance between source and station, taking account of
! multiple orbits, when the parameter WINDOW_START, set by commonModules,
! indicates so.
!
! input:
!     distance, iorbit      -   epicentral distance between source/receiver, orbit number
!
! returns: distance (in rad) and iorbit are set accordingly
  use precisions
  use cells, only: vertices
  use propagationStartup
  implicit none
  real(WP),intent(out):: distance
  integer,intent(out):: iorbit
  ! local parameters
  real(WP):: minor_distance,arrivalTime
  real(WP):: vectorA(3),vectorB(3)

  ! gets minor arc distance (in rad)
  vectorA(:) = vertices(sourceVertex,:)
  vectorB(:) = vertices(receiverVertex,:)
  call greatCircleDistance(vectorA,vectorB,minor_distance)
  distance = minor_distance

  ! higher orbits
  iorbit = 1
  arrivalTime = distance*EARTHRADIUS/cphaseRef
  do while( WINDOW_START > arrivalTime )
    iorbit = iorbit + 1
    if (mod(iorbit,2) == 0) then
      ! even orbits
      distance = ((iorbit)/2.0)*2.0*PI - minor_distance
      arrivalTime = distance*EARTHRADIUS/cphaseRef
    else
      ! odd orbits
      distance = ((iorbit-1)/2.0)*2.0*PI + minor_distance
      arrivalTime = distance*EARTHRADIUS/cphaseRef
    endif
  enddo

  end subroutine


!-----------------------------------------------------------------------
  subroutine determineSourceParameters()
!-----------------------------------------------------------------------
! sets parameters sigma and mu which govern the Gaussian source characteristics
!
! returns: muSquare,muTwo,TimeParameterSigma,WidthParameterMu
  use precisions
  use propagationStartup; use parallel
  use cells, only: subdivisions
  use filterType; use verbosity
  implicit none
  ! local parameters
  real(WP):: theta_radian,tmp

  ! calculates width parameter (see Carl Tape Thesis , eq. 3.21 )
  if (.not. FIXED_SOURCEPARAMETER) then
    !WidthParameterMu=FACTOR_WIDTH*averageCellDistance/(EARTHRADIUS*2.0*sqrt(-2.0*log(0.05)))
    theta_radian = THETA_WIDTH * PI/180.0
    WidthParameterMu = theta_radian / (2.0*sqrt(-2.0*log(0.05)))
  endif
  muSquare = WidthParameterMu * WidthParameterMu
  muTwo    = muSquare + muSquare

  ! empirical formula to obtain a maximum at the reference period
  if (ADAPT_SOURCE) then
    if ( (20.+WidthParameterMu*7000.) > bw_waveperiod) then
      print *,'Error: estimated minimum period possible: ',20.+WidthParameterMu*7000.
      print *,'                     width parameter mu : ',WidthParameterMu
      print *,'       minimum period should be smaller than bandwidth wave period ',bw_waveperiod
      call stopProgram("source not possible    ")
    endif
    tmp = 1.0 - 8.0/3.0 * (20.+WidthParameterMu*7000.-bw_waveperiod)
    if (tmp < 0.0) then
      call stopProgram("source parameter sigma not defined    ")
    endif
    TimeParameterSigma = 0.75 + 0.75 * sqrt( tmp )
  endif

  tmp = 2./3. * TimeParameterSigma**2 - TimeParameterSigma + 20. + 7000. * WidthParameterMu
  ! console output
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'  source parameters:'
    print *,'    time parameter sigma  :',TimeParameterSigma
    print *,'    width parameter mu    :',WidthParameterMu
    print *,'    empirical spectral maximum around period:',tmp
    print *
  endif

  ! checks width parameter
  if (subdivisions <= 6 .and. WidthParameterMu < 0.008) then
    print *,'    estimated minimum period possible: ',20. + WidthParameterMu*7000.
    print *
    !call stopProgram("source not optimal ")
  else if (subdivisions <= 7 .and. WidthParameterMu < 0.004) then
    print *,'    estimated minimum period possible: ',20. + WidthParameterMu*7000.
    print *
    !call stopProgram("source not optimal ")
  endif

  end subroutine

