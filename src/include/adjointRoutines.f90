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
  subroutine prepareKernel(kernelIncrement)
!-----------------------------------------------------------------------
! prepare model for a new simulation with the setup of a new receiver location
!
! input:
!     kernelIncrement     -   integer of degrees to increment the epicentraldistance
!                                        between source and station
!
! returns: newly initialize station setup
  use adjointVariables; use propagationStartup; use parallel
  use cells; use verbosity
  implicit none
  integer, intent(in):: kernelIncrement
  ! local parameters
  character(len=3):: kernelstr

  ! append kernel number to name of adjoint kernel file
  !(originally something like 'adjointKernel.dat' should become
  ! e.g. 'adjointKernel.023.dat' for epicentral distance 23 degree )
  write(kernelstr,'(i3.3)') int(kernelStartDistance+kernelIncrement-1)
  if (kernelIncrement == 1) then
    adjointKernelName(len_trim(adjointKernelName)-2:len_trim(adjointKernelName)+4)&
            = kernelstr//'.dat'
  else
    adjointKernelName(len_trim(adjointKernelName)-6:len_trim(adjointKernelName))&
            = kernelstr//'.dat'
  endif

  if (MAIN_PROCESS .and. VERBOSE) then
    print *
    print *,'kernel name: ',adjointKernelName
  endif

  ! receivers will be placed on the equator
  desiredReceiverLat = 0.0_WP
  desiredReceiverLon = kernelStartDistance + kernelIncrement - 1
  call setupStation( desiredReceiverLat, desiredReceiverLon )

  ! initialize new model with new simulation time
  if (USE_OVERTIME) then
    call initializeWorld()
    ! build new source (only when number of time steps changed)
    if (size(forceTermPrescribed(1,:)) /= numofTimeSteps) call initializeSource()
  else
    call initializeAdjointArrays()
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine initializeAdjointArrays()
!-----------------------------------------------------------------------
! allocates dynamically the arrays:   adjointKernel, adjointSource, wavefieldForward
! and initializes them
  use propagationStartup;use parallel;use adjointVariables
  use verbosity; use cells
  implicit none
  ! local parameters
  real(WP)::memory
  integer:: ier

  ! allocates array for storing the adjoint kernel values
  if (.not. allocated(adjointKernel)) then
    if (MAIN_PROCESS .and. VERBOSE) then
      print *,'  allocating adjointKernel array...'
      print *,'    size: ',numVertices*WP/1024./1024.,'Mb'
    endif
    allocate(adjointKernel(numVertices), stat=ier)
    if (ier /= 0) call stopProgram('error allocating adjoint kernel arrays    ')
  endif
  ! initializes array
  adjointKernel(:) = 0.0_WP

  ! reallocates velocity seismogram
  ! (maybe the number of time steps changed..)
  if (allocated(adjointSource)) then
    if (size(adjointSource(1,:)) /= numofTimeSteps) deallocate(adjointSource)
  endif

  if (.not. allocated(adjointSource)) then
    allocate(adjointSource(2,numofTimeSteps),stat=ier)
    if (ier /= 0) call stopProgram('error allocating adjoint source arrays   ')
  endif

  ! initializes the adjoint source trace
  adjointSource(:,:)= 0.0_WP

  ! machine memory holds for 2 GB RAM:
  !   level 6: numVertices=122'882, numofTimeSteps~500, double precision 8 byte
  !               -> needs ~ 500 MB per wavefield, still o.k.
  !   level 7: numVertices=491'522, numofTimeSteps~100, dp 8 byte
  !               -> needs ~ 3.8 GB ! per wavefield, too big
  !
  ! for faster computation try to store values in arrays than files:
  ! trys to reallocate wavefield arrays
  ! ( maybe time steps changed for new kernel calculation...)
  if (allocated(wavefieldForward)) then
    if (size(wavefieldForward(1,:)) /= numofTimeSteps) deallocate(wavefieldForward)
  endif

  storeAsFile = .false.
  if (.not. allocated(wavefieldForward)) then
    memory = numDomainVertices*numofTimeSteps*WP
    if (MAIN_PROCESS .and. VERBOSE) then
      print *,'  allocating wavefieldForward array...'
      print *,'    vertices   : ',numDomainVertices
      print *,'    time steps : ',numofTimeSteps
      print *,'    size       : ',memory/1024./1024.,'Mb'
    endif
    allocate(wavefieldForward(numDomainVertices,numofTimeSteps),stat=ier)
    if (ier /= 0) then
      ! sets some reasonable limit to computation:
      ! required harddisk memory shouldn't be bigger than 10 GB,
      ! else the kernels will waste too much space for storage
      if (memory > 10.0e9) call stopProgram('too much memory use to store data   ')
      storeAsFile = .true.
    endif
    ! ensures that flag is equal for all processes
    call syncFlag(myrank, nprocesses,storeAsFile)
  endif

  ! storage of adjoint wavefields
  if (storeAsFile) then
    ! free up space if they were allocated successfully
    if (allocated(wavefieldForward)) deallocate(wavefieldForward)

    ! console output
    if (MAIN_PROCESS .and. VERBOSE) then
      print *
      print *,'Error: adjoint arrays - not available RAM memory size:', memory/1024./1024.,'Mb'
      print *
      print *,'adjoint wavefields would be stored as data-files:'
      print *,'    - slower performance expected'
      print *,'    - large storage files expected'
      print *
      call stopProgram("error adjoint arrays  ")
    endif
  endif

  ! initializes
  wavefieldForward(:,:) = 0.0_WP

  ! stores adjoint, backpropagating wavefield
  if (.not. ADJOINT_ONTHEFLY) then
    ! checks if already an instance exists
    if (allocated(wavefieldAdjoint)) then
      if (size(wavefieldAdjoint(1,:)) /= numofTimeSteps) deallocate(wavefieldAdjoint)
    endif

    ! allocates memory
    if (.not. allocated(wavefieldAdjoint)) then
      if (MAIN_PROCESS .and. VERBOSE) then
        print *,'  allocating wavefieldAdjoint array...'
        print *,'    size       : ',numDomainVertices*numofTimeSteps*WP/1024./1024.,'Mb'
      endif
      allocate(wavefieldAdjoint(numDomainVertices,numofTimeSteps),stat=ier)
      if (ier /= 0) call stopProgram('error allocating adjoint wavefield    ')
    endif

    ! initializes
    wavefieldAdjoint(:,:) = 0.0_WP
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine storeForwardDisplacements(timestep,index)
!-----------------------------------------------------------------------
! saves displacements of a forward simulation for future adjoint calculation
!
! input:  timestep     -  time step in iterative scheme
!            index          - iteration number
! returns: stores newdisplacements either in file or in array wavefieldForward
  use adjointVariables;use displacements;use parallel;use propagationStartup
  use verbosity;use griddomain; use cells
  implicit none
  integer,intent(in):: timestep,index
  ! local parameters
  integer::vertex,n,ier
  character(len=8):: timestepstr
  character(len=3):: rankstr

  if (storeAsFile) then
    ! collect displacements
    !if (PARALLELSEISMO) then
    !  call collectFullNewdisplacement()
    !endif
    !if (MAIN_PROCESS) then
    !endif

    ! output to file for each process
    write(timestepstr,'(i8.7)')timestep
    write(rankstr,'(i3.3)') myrank
    open(IOUT,file=trim(datadirectory)//'wavefield_'//timestepstr//'.rank'//rankstr//'.dat', &
         access='direct',recl=WP,iostat=ier)
    if (ier /= 0) call stopProgram('could not open file: wavefield_'//timestepstr//'.rank'//rankstr//'.dat ....     ')
    ! fill in displacements
    !do i=1,numVertices
    !  call getSphericalCoord_Lat(i,lat,lon)
    !  write(10,*) lon,lat,newdisplacement(i)
    !enddo
    do n = 1, numDomainVertices
      ! choose vertex
      if (PARALLELSEISMO) then
        vertex = domainVertices(n)
      else
        vertex = n
      endif
      !call getSphericalCoord_Lat(vertex,lat,lon)
      !write(10,*) lon,lat,newdisplacement(vertex)
      write(IOUT,rec=n) newdisplacement(vertex)
    enddo
    close(IOUT)
  else
    ! store in wavefield array
    do n = 1, numDomainVertices
      ! choose vertex
      if (PARALLELSEISMO) then
        vertex = domainVertices(n)
      else
        vertex = n
      endif
      wavefieldForward(n,index) = newdisplacement(vertex)
    enddo
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine storeBackwardDisplacements(timestep,index)
!-----------------------------------------------------------------------
! saves displacements for future adjoint calculation
!
! input:  timestep     -  time step in iterative scheme
!            index          - iteration number
! returns: stores newdisplacements either in file or in array wavefieldAdjoint
  use adjointVariables;use displacements;use parallel;use propagationStartup
  use verbosity;use griddomain; use cells
  implicit none
  integer,intent(in):: timestep,index
  ! local parameters
  integer::vertex,n,ier
  character(len=8):: timestepstr
  character(len=3):: rankstr

  if (storeAsFile) then
    ! output to file for each process
    write(timestepstr,'(i8.7)')timestep
    write(rankstr,'(i3.3)') myrank
    open(IOUT,file=trim(datadirectory)//'wavefieldAdj_'//timestepstr//'.rank'//rankstr//'.dat', &
          access='direct',recl=WP,iostat=ier)
    if (ier /= 0) call stopProgram('could not open file wavefieldAdj_'//timestepstr//'.rank'//rankstr//'.dat ....     ')
    ! fill in displacements
    do n = 1, numDomainVertices
      ! choose vertex
      if (PARALLELSEISMO) then
        vertex = domainVertices(n)
      else
        vertex = n
      endif
      write(IOUT,rec=n) newdisplacement(vertex)
    enddo
    close(IOUT)
  else
    ! store in adjoint wavefield array
    do n = 1,numDomainVertices
      ! choose vertex
      if (PARALLELSEISMO) then
        vertex = domainVertices(n)
      else
        vertex = n
      endif
      wavefieldAdjoint(n,index) = newdisplacement(vertex)
    enddo
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine getAdjointSource()
!-----------------------------------------------------------------------
! determine adjoint source
! (Tromp et al., 2005, sec 4.1 eq (45) )
  use propagationStartup; use parallel
  use verbosity; use filterType; use adjointVariables
  implicit none
  ! local parameters
  real(WP)::seismo(numofTimeSteps),seismoDisplacement(numofTimeSteps)
  real(WP),dimension(:,:),allocatable:: window_signal
  integer::n,i,ier
  integer:: window_startindex,window_endindex,window_range
  real(WP):: normFactor

  ! console output
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'creating adjoint source...'
    if (ADJOINT_TAPERSIGNAL) &
      print *,'  using tapering'
    if (WINDOWEDINTEGRATION) &
      print *,'  using window  : start time =',WINDOW_START,' end time = ',WINDOW_END
  endif

  ! synchronize seismogram at receiver location
  if (PARALLELSEISMO) call syncReceiverSeismogram()

  ! filter around corner frequency
  !if (FILTER_SEISMOGRAMS) then
  !  if (MAIN_PROCESS .and. VERBOSE) print *,'    using filtered receiver seismogram...'
  !  if (.not. MAIN_PROCESS) beVerbose=.false.
  !  call dofilterSeismogram(seismogramReceiver,numofTimeSteps)
  !endif

  ! get receiver's velocity seismogram (first time derivative of recorded seismogram at receiver)
  seismo(:) = 0.0_WP
  seismoDisplacement(:) = seismogramReceiver(2,:)

  if (MAIN_PROCESS .and. fileOutput) then
    print *,'  printing to file: ',trim(datadirectory)//'seismo.displacement.dat'
    open(IOUT,file=trim(datadirectory)//'seismo.displacement.dat')
    do n = 1,numofTimeSteps
      write(IOUT,*) seismoDisplacement(n)
    enddo
    close(IOUT)
  endif

  call getFirstDerivative(seismoDisplacement,seismo,dt,numofTimeSteps)

  ! file output
  if (MAIN_PROCESS .and. fileOutput) then
    print *,'  printing to file: ',trim(datadirectory)//'seismo.velocity.dat'
    open(IOUT,file=trim(datadirectory)//'seismo.velocity.dat')
    do n = 1,numofTimeSteps
      write(IOUT,*) seismogramReceiver(1,n),seismo(n)
    enddo
    close(IOUT)
  endif

  ! adjoint source is the velocity seismogram reversed in time
  do i = 1,numofTimeSteps
    adjointSource(1,i) = seismogramReceiver(1,numofTimeSteps-i+1)
    adjointSource(2,i) = seismo(numofTimeSteps-i+1)
  enddo

  ! determines the window for the cross-correlation ( Tromp et al, 2005: eq. 41, parameter w_r(t) )
  if (WINDOWEDINTEGRATION) call determineWindowsize(window_startindex,window_endindex,window_range)

  ! apply hanning window to smooth adjoint source ends
  if (ADJOINT_TAPERSIGNAL) then
    if (WINDOWEDINTEGRATION) then
      ! cuts out a window of the adjoint signal
      !if (MAIN_PROCESS .and. VERBOSE) print *,'    allocating window...',window_startindex,window_endindex
      allocate(window_signal(2,window_startindex-window_endindex-1),stat=ier)
      if (ier /= 0) call stopProgram('error - window_signal ')

      ! sets window which contains signal
      window_signal(1,:) = adjointSource(1,window_endindex+1:window_startindex-1)
      window_signal(2,:) = adjointSource(2,window_endindex+1:window_startindex-1)

      ! tapers window
      call taperSeismogram(window_signal,window_startindex-window_endindex-1, &
                           window_startindex-window_endindex-1,.false.)
      adjointSource(2,window_endindex+1:window_startindex-1) = window_signal(2,:)
    else
      call taperSeismogram(adjointSource,numofTimeSteps,numofTimeSteps,.false.)
    endif
  endif

  ! normalization factor (Tromp et al., 2005, sec. 4.1 eq. (42) )
  ! get second time derivative
  seismo(:) = 0.0
  call getSecondDerivative(seismoDisplacement,seismo,dt,numofTimeSteps)

  ! file output
  if (MAIN_PROCESS .and. fileOutput) then
    print *,'  printing to file: ',trim(datadirectory)//'seismo.acceleration.dat'
    open(IOUT,file=trim(datadirectory)//'seismo.acceleration.dat')
    do n = 1,numofTimeSteps
      write(IOUT,*) seismogramReceiver(1,n),seismo(n)
    enddo
    close(IOUT)
    if (WINDOWEDINTEGRATION) then
      print *,'  printing to file: ',trim(datadirectory)//'seismo.window.dat'
      open(IOUT,file=trim(datadirectory)//'seismo.window.dat')
      do n = 1,numofTimeSteps
        write(IOUT,*) adjointSource(1,n),adjointSource(2,n)
      enddo
      close(IOUT)
    endif
  endif

  ! normalization factor is the time integral
  normFactor = 0.0
  if (WINDOWEDINTEGRATION) then
    ! windows signal ( index = numofTimeSteps-i+1 )
    seismoDisplacement(1:numofTimeSteps-window_startindex+1) = 0.0
    seismoDisplacement(numofTimeSteps-window_endindex+1:numofTimeSteps) = 0.0
  endif
  call getIntegral(seismo,seismoDisplacement,normFactor,dt,numofTimeSteps)

  ! synchronizes
  call syncProcesses()

  ! check for division by zero
  if (abs(normFactor) < 1e-6) then
    if (MAIN_PROCESS) print *,'Error: norm:',normFactor,' index:',window_startindex,window_endindex
    call stopProgram('normalization factor zero   ')
  endif

  ! normalize source
  adjointSource(2,:) = adjointSource(2,:)/normFactor

  ! adjoint source localization is at the receiver
  adjointSourceVertex = receiverVertex

  ! console & file output
  if (MAIN_PROCESS .and. VERBOSE) then
    ! debug normalization output
    print *,'  normalization factor: ',normFactor

    ! store adjoint source seismogram
    if (fileOutput) then
      print *,'  printing to file: ',trim(datadirectory)//'seismo.adjointSource.dat'
      open(IOUT,file=trim(datadirectory)//'seismo.adjointSource.dat')
      do n = 1,numofTimeSteps
        write(IOUT,*) adjointSource(1,n),adjointSource(2,n)
      enddo
      close(IOUT)
    endif
  endif

  ! precalculate the second time derivatives of all vertices belonging to this process
  if (ADJOINT_ONTHEFLY .and. PRECALCULATE_DERIVATIVES) call precalculateSecondDerivatives()

  end subroutine


!-----------------------------------------------------------------------
  subroutine precalculateSecondDerivatives()
!-----------------------------------------------------------------------
! calculates the second derivatives of the forward wavefield seismogram at each vertex
!
! returns: stores the second derivatives in the wavefieldForward array
  use propagationStartup; use parallel; use griddomain; use filterType
  use adjointVariables;use verbosity; use cells
  implicit none
  ! local parameters
  integer:: n,vertex !,index
  real(WP):: seismo(2,numofTimeSteps),seismoTmp(numofTimeSteps)

  ! console output
  if (MAIN_PROCESS .and. VERBOSE) print *,'precalculate time derivatives...'

  ! every vertex has its own seismogram to derive
  seismo(1,:) = seismogramReceiver(1,:) ! time is the same as from the receiver station seismogram
  if (MAIN_PROCESS .and. VERBOSE .and. FILTER_SEISMOGRAMS) print *,'    using filtered forward seismograms...'

  do n = 1,numDomainVertices
    ! get cell vertex
    if (PARALLELSEISMO) then
      vertex = domainVertices(n)
    else
      vertex = n
    endif

    ! get corresponding seismogram
    if (storeAsFile) then
      call getForwardWave(n,seismoTmp)
    else
      seismoTmp(:) = wavefieldForward(n,:)
    endif
    seismo(2,:) = seismoTmp(:)

    ! filter seismogram
    if (FILTER_SEISMOGRAMS) then
      call dofilterSeismogram(seismo,numofTimeSteps)
      seismoTmp(:) = seismo(2,:)
    endif

    ! compute second time derivative (by spline)
    call getSecondDerivative(seismo(2,:),seismoTmp,dt,numofTimeSteps)

    ! compute second time derivative (by central-difference time scheme)
    !do index=2,numofTimeSteps-1
    !  seismoSecondDerivative(vertex,index)=(seismoTmp(2,index+1) &
    !                - 2.0_WP*seismoTmp(2,index) &
    !                + seismoTmp(2,index-1))/dt2
    !enddo

    ! store as new seismogram
    if (storeAsFile) then
      call setForwardWave(n,seismoTmp)
    else
      wavefieldForward(n,:)=seismoTmp(:)
    endif
  enddo

  end subroutine


!-----------------------------------------------------------------------
  subroutine backwardIteration()
!-----------------------------------------------------------------------
! (back-propagation) simulation with the adjoint source
  use propagationStartup; use cells
  use parallel; use displacements; use verbosity
  use griddomain, only: domainVertices
  use phaseVelocityMap, only: phaseVelocitySquare
  use adjointVariables
  implicit none
  ! local parameters
  integer:: n,timestep,vertex,index,ier
  real(WP):: u_t,u_tplus1,u_tminus1,forcing,D2,time
  real(WP), external:: forceAdjointSource,discreteLaplacian,precalc_discreteLaplacian
  real(WP), external:: precalc_backdiscreteLaplacian
  real(WP), external:: syncWtime

  ! benchmark
  if (MAIN_PROCESS .and. VERBOSE) benchstart = syncWtime()

  if (ADJOINT_ONTHEFLY) then
    ! simultaneous backward displacement (same calculation as forward one but reversed) initial start
    !backwardDisplacement_old(:) = newdisplacement(:)
    !backwardDisplacement(:) = newdisplacement(:)
    !backwardNewDisplacement(:) = displacement(:)

    ! look for vertex in middle of source/receiver
    call findVertex(sourceLat+nint((receiverLat-sourceLat)/2.0), &
                    sourceLon+nint((receiverLon-sourceLon)/2.0),midpointVertex)

    if (MAIN_PROCESS .and. VERBOSE) &
      print *,'   storing file: ',trim(datadirectory)//'seismo.integral_source.dat'
    open(adjSourceFileID,file=trim(datadirectory)//'seismo.integral_source.dat')
    if (MAIN_PROCESS .and. VERBOSE) &
      print *,'   storing file: ',trim(datadirectory)//'seismo.integral_rec.dat'
    open(adjRecFileID,file=trim(datadirectory)//'seismo.integral_rec.dat')
    if (MAIN_PROCESS .and. VERBOSE) &
      print *,'   storing file: ',trim(datadirectory)//'seismo.integral_midpoint.dat'
    open(adjMidpointFileID,file=trim(datadirectory)//'seismo.integral_midpoint.dat')
  else
    if (.not. storeAsFile) then
      ! allocate adjoint wavefield
      if (.not. allocated(wavefieldAdjoint)) then
        if (MAIN_PROCESS .and. VERBOSE) then
          print *,'    allocating adjoint wavefield...'
          print *,'      vertices    : ',numDomainVertices
          print *,'      time steps  : ',numofTimeSteps
          print *,'      array size  : ',numDomainVertices*numofTimeSteps*WP/1024./1024.,'Mb'
        endif

        allocate(wavefieldAdjoint(numDomainVertices,numofTimeSteps),stat=ier)
        if (ier /= 0) call stopProgram('error allocating adjoint wavefield    ')
        !initialize
        wavefieldAdjoint(:,:) = 0.0_WP
      endif
    endif
  endif

  ! reset displacements for adjoint propagation
  displacement_old(:) = 0.0_WP
  displacement(:) = 0.0_WP
  newdisplacement(:) = 0.0_WP

  ! time iteration of displacements
  if (MAIN_PROCESS .and. VERBOSE) print *,'  backward iteration...'
  index = 0
  do timestep = lasttimestep,firsttimestep,-1
    ! model time
    time  = timestep*dt
    index = index+1

    ! swap displacement arrays
    displacement_old(:) = displacement(:)
    displacement(:)     = newdisplacement(:)

    ! propagate only corresponding vertices
    do n = 1, numDomainVertices
      ! choose vertex
      if (PARALLELSEISMO) then
        vertex = domainVertices(n)
      else
        vertex = n
      endif

      ! spherical Laplacian
      if (PRECALCULATED_CELLS) then
        D2 = precalc_discreteLaplacian(vertex) !uses displacement(..) field
      else
        D2 = discreteLaplacian(vertex)
      endif

      ! calculate new displacement in time
      u_t       = displacement(vertex)
      u_tminus1 = displacement_old(vertex)

      ! determines adjoint source
      forcing = forceAdjointSource(vertex,index)

      ! gets phase velocity square
      cphase2 = phaseVelocitySquare(vertex)

      !propagation step
      u_tplus1 = u_t + u_t - u_tminus1 + cphase2*dt2*(D2+forcing)

      ! iterated displacements
      newdisplacement(vertex) = u_tplus1

      ! computes adjoint kernel value step by step
      if (ADJOINT_ONTHEFLY) then
        call getAdjointKernel_fly(vertex,timestep)
      endif
    enddo

    ! synchronize new displacement arrays
    if (PARALLELSEISMO) then
      call syncNewdisplacement()
    endif

    ! save displacements at each time step for future adjoint calculation
    if (.not. ADJOINT_ONTHEFLY) then
      call storeBackwardDisplacements(timestep,index)
    endif

    !file output for simulation snapshots
    if (SIMULATIONOUTPUT) then
      call printAdjointWavefield(timestep,time)
    endif
  enddo !timestep

  ! on the fly calculation needs to be scaled finally
  if (ADJOINT_ONTHEFLY) then
    call scaleAdjointKernel_fly()
  endif

  ! free memory to be able to have some to collect data
  !if (PARALLELSEISMO .and. (.not. kernelIteration)) &
  !      deallocate(boundaries,sendDisp,receiveDisp,domainNeighbors)

  ! benchmark output
  if (MAIN_PROCESS .and. VERBOSE) then
    benchend = syncWtime()
    print *,'  benchmark seconds ',benchend-benchstart
    print *
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine getAdjointKernel_fly(vertex,timestep)
!-----------------------------------------------------------------------
! determine sensitivity kernel value when calculating "on the fly"
!
! return: stores values in the adjointKernel array
  use propagationStartup; use displacements; use phaseVelocityMap; use traveltime
  use adjointVariables; use parallel; use griddomain; use cells; use verbosity
  implicit none
  integer,intent(in):: vertex,timestep
  ! local parameters
  integer:: n,index,vertexIndex
  real(WP):: seismo(3),seismoAdjoint(3),seismoTmp(numofTimeSteps) !,timewindow
  real(WP):: kernelVal,val1,val2,for1,for2,for3 !,for4
  !real(WP):: derivativeActual,derivativeNext,
  real(WP):: kernelfactor

  ! current time step index
  index = timestep - firsttimestep  + 1

  ! time window for integration
  if (WINDOWEDINTEGRATION) then
    if (timestep*dt < WINDOW_START .or. timestep*dt > WINDOW_END) then
      return
    endif
  endif

  ! index for vertex location in the wavefield files
  if (PARALLELSEISMO) then
    do n = 1,numDomainVertices
      if (vertex == domainVertices(n)) then
        vertexIndex = n
        exit
      endif
    enddo
  else
    vertexIndex = vertex
  endif

  ! calculate kernel value with actual available displacement values
  if (ADJOINT_INTEGRALBYSPLINE) then
    ! time integral calculated as a sum of discrete rectangles with size dt
    !if ( (index > numofTimeSteps-2) .or. (index < 1)) then
    !  return
    !endif

    ! get seismogram
    !if (storeAsFile) then
    !  call getForwardWave(vertexIndex,seismoTmp)
    !else
    !  seismoTmp(:)=wavefieldForward(vertexIndex,:)
    !endif
    !seismo(1)=seismoTmp(index)
    !seismo(2)=seismoTmp(index+1)
    !seismo(3)=seismoTmp(index+2)

    seismo(1) = backwardNewdisplacement(vertex)
    seismo(2) = backwardDisplacement(vertex)
    seismo(3) = backwardDisplacement_old(vertex)

    seismoAdjoint(1) = newdisplacement(vertex)
    seismoAdjoint(2) = displacement(vertex)
    seismoAdjoint(3) = displacement_old(vertex)

    ! time derivative of seismo
    if (.not. PRECALCULATE_DERIVATIVES) then
      if (seismo(1) == 0.0 .and. seismo(2) == 0.0 .and. seismo(3) == 0.0) then
        continue
      else
        call getSecondDerivative(seismo,seismo,dt,3)
      endif
    endif

    ! kernel value is the time integral
    kernelVal = 0.0
    val1 = seismo(1)*seismoAdjoint(1)
    val2 = seismo(2)*seismoAdjoint(2)
    kernelVal = 0.5_WP*(val1+val2)*dt
  else
    ! time integral calculated as a sum of discrete rectangles with size dt
    !if ( (index > numofTimeSteps-1).or.(index < 2)) then
    !  return
    !endif

    ! get seismogram
    !seismo(1) = backwardNewdisplacement(vertex)
    !seismo(2) = backwardDisplacement(vertex)
    !seismo(3) = backwardDisplacement_old(vertex)
    ! get adjoint seismogram
    !seismoAdjoint(1) = newdisplacement(vertex)
    seismoAdjoint(2) = displacement(vertex)
    !seismoAdjoint(3) = displacement_old(vertex)

    ! time derivative of seismo (by central-differences)
    !val1 = (seismo(1)+seismo(3)-2.0_WP*seismo(2))/dt2
    !kernelval = val1*seismoAdjoint(2)*dt

    ! get seismogram
    if (storeAsFile) then
      call getForwardWave(vertexIndex,seismoTmp)
    else
      seismoTmp(:) = wavefieldForward(vertexIndex,:)
    endif

    ! compute second time derivative (by central-difference time scheme)
    if (PRECALCULATE_DERIVATIVES) then
      val1 = seismoTmp(index)*newdisplacement(vertex)
      val2 = seismoTmp(index+1)*displacement(vertex)
    else
      for1 = seismoTmp(index-1)
      for2 = seismoTmp(index)
      for3 = seismoTmp(index+1)
    !  for4 = seismoTmp(index+2)

    !  !seismoSecondDerivative(vertex,index) = (for1 + for3 - 2.0_WP*for2)/dt2
    !  !seismoSecondDerivative(vertex,index+1) = (for2 + for4 - 2.0_WP*for3)/dt2
    !  derivativeActual = (for1 + for3 - 2.0_WP*for2)/dt2
    !  derivativeNext = (for2 + for4 - 2.0_WP*for3)/dt2

      ! kernel value increment
    !  !val1 = seismoSecondDerivative(vertex,index)*newdisplacement(vertex)
    !  !val2 = seismoSecondDerivative(vertex,index+1)*displacement(vertex)
    !  val1 = derivativeActual*newdisplacement(vertex)
    !  val2 = derivativeNext*displacement(vertex)
      val1 = (for1 + for3 - 2.0_WP*for2)/dt2
    endif
    !kernelVal = 0.5_WP*(val1+val2)*dt
    kernelVal = val1*seismoAdjoint(2)*dt

  endif

  ! calculate kernel factor
  kernelfactor = 2.0_WP/(phaseVelocitySquare(vertex))
  kernelVal = kernelfactor*kernelVal

  ! sum up to build adjointKernel value
  adjointKernel(vertex) = adjointKernel(vertex)+kernelVal

  ! output to files (timestep: correct +1 by back-propagation and to account
  ! for values stored in displacement() and not newdisplacement() )
  if (vertex == sourceVertex) then
    write(adjSourceFileID,'(6e16.4e3)') (timestep+1)*dt,seismo(2), &
        (timestep+1)*dt,seismoAdjoint(2),adjointKernel(vertex),kernelVal/kernelFactor
  endif
  ! debug output
  if (vertex == receiverVertex) then
    write(adjRecFileID,'(6e16.4e3)') (timestep+1)*dt,seismo(2), &
        (timestep+1)*dt,seismoAdjoint(2),adjointKernel(vertex),kernelVal/kernelFactor
  endif
  ! debug output
  if (vertex == midpointVertex) then
    write(adjMidpointFileID,'(6e16.4e3)') (timestep+1)*dt,seismo(2), &
        (timestep+1)*dt,seismoAdjoint(2),adjointKernel(vertex),kernelVal/kernelFactor
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine scaleAdjointKernel_fly()
!-----------------------------------------------------------------------
  use adjointVariables
  use cells; use propagationStartup
  implicit none
  ! local parameters
  real(WP):: scalefactor,vertexCellArea,distance,arrivalTime
  integer:: iorbit

  ! close output files
  close(adjSourceFileID)
  close(adjRecFileID)
  close(adjMidpointFileID)

  ! cell area in [rad^2]
  vertexCellArea = cellAreas(receiverVertex)/EARTHRADIUS_SQUARED

  ! calculate the reference travel time
  call determineOrbitDistance(distance,iorbit)
  arrivalTime = distance*EARTHRADIUS/cphaseRef

  ! scale each entry of the adjoint kernel
  scalefactor      = 1.0/(arrivalTime*vertexCellArea)
  adjointKernel(:) = scalefactor*adjointKernel(:)

  ! free some memory
  deallocate(backwardDisplacement,backwardDisplacement_old,backwardNewdisplacement)

  end subroutine


!-----------------------------------------------------------------------
  subroutine frechetKernel()
!-----------------------------------------------------------------------
! determine adjoint kernel value when calculating NOT "on the fly"
!
! return: stores values in the adjointKernel() array
  use propagationStartup; use parallel; use griddomain
  use phaseVelocityMap; use filterType; use verbosity; use cells; use adjointVariables
  implicit none
  ! local parameters
  integer:: index,indexadjoint,vertex,n
  real(WP):: kernelVal,val1,val2,derivativeActual,derivativeNext,kernelfactor
  real(WP):: arrivalTime,distance,vertexCellArea,vertexCellAreaRad,timewindow,for1,for2,for3
  real(WP)::seismo(2,numofTimeSteps),seismoAdjoint(2,numofTimeSteps),seismoTmp(numofTimeSteps)

  ! console output
  if (MAIN_PROCESS .and. VERBOSE) print *,'  calculating kernel values...'

  ! initialize with time (dealing with newdisplacements means at time steps t+dt)
  do index = 1,numofTimeSteps
    seismo(1,index)=(firsttimestep+index-1)*dt+dt
    seismoAdjoint(1,index) = (lasttimestep-index+1)*dt+dt
  enddo

  ! calculate the reference travel time
  ! attention: for a heterogeneous background it takes here the PREM value as well.
  !                this should be considered when inverting and using these kernels.
  !call determineOrbitDistance(distance,iorbit)
  !arrivalTime=distance*EARTHRADIUS/cphaseRef

  ! gets minor arc distance (in rad)
  call greatCircleDistance(vertices(sourceVertex,:),vertices(receiverVertex,:),distance)
  arrivalTime = distance*EARTHRADIUS/cphaseRef

  ! cell area
  vertexCellArea = cellAreas(receiverVertex)
  ! convert cell area into [rad^2]
  vertexCellAreaRad = vertexCellArea/EARTHRADIUS_SQUARED

  ! console output
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'    reference travel time [s]      : ',arrivalTime
    !for [rad2]: vertexCellArea/(EARTHRADIUS*EARTHRADIUS)
    print *,'    receiver cell area [km2]       : ',vertexCellArea
    print *,'    time integration:'
    print *,'    starting seismogram at         : ',seismo(1,1)
    print *,'    ending seismogram at           : ',seismo(1,numofTimeSteps)
    if (WINDOWEDINTEGRATION) &
    print *,'    windowed adjoint source between: ',WINDOW_START,WINDOW_END
    print *
  endif

  ! look for vertex in middle of source/receiver
  !call findVertex(sourceLat+nint((receiverLat-sourceLat)/2.0), &
  !          sourceLon+nint((receiverLon-sourceLon)/2.0),midpointVertex)

  ! determine kernel value at each grid point for the corresponding processor domain
  do n = 1, numDomainVertices
    ! choose vertex
    if (PARALLELSEISMO) then
      vertex = domainVertices(n)
    else
      vertex = n
    endif

    ! get forward seismogram
    if (storeAsFile) then
      call getForwardWave(n,seismoTmp)
      seismo(2,:) = seismoTmp(:)
    else
      seismo(2,:) = wavefieldForward(n,:)
    endif

    ! get adjoint seismogram
    if (storeAsFile) then
      call getAdjointWave(n,seismoTmp)
      seismoAdjoint(2,:) = seismoTmp(:)
    else
      seismoAdjoint(2,:) = wavefieldAdjoint(n,:)
    endif

    ! filter (and taper) the seismograms
    if (FILTER_SEISMOGRAMS) then
      ! filter
      if (beverbose) print *,'  filtering seismograms...'
      call dofilterSeismogram(seismo,numofTimeSteps)
      call dofilterSeismogram(seismoAdjoint,numofTimeSteps)
      ! be verbose only once
      beverbose = .false.
    endif

    ! kernel value
    kernelVal = 0.0_WP
    seismoTmp(:) = 0.0_WP
    if (ADJOINT_INTEGRALBYSPLINE) then
      ! compute second derivative of forward seismogram
      call getSecondDerivative(seismo(2,:),seismo(2,:),dt,numofTimeSteps)

      ! kernel value is the time integral
      do index = 1,numofTimeSteps-1
        ! time window for integration is taken account of in building the adjoint source

        derivativeActual = seismo(2,index)
        derivativeNext = seismo(2,index+1)

        ! kernel value is the time integral (by a sum of discrete rectangles with size dt)
        val1 = derivativeActual*seismoAdjoint(2,numofTimeSteps-index+1)
        val2 = derivativeNext*seismoAdjoint(2,numofTimeSteps-index+1-1)

        kernelVal = kernelVal+0.5_WP*(val1+val2)*dt*timewindow

        ! check & store temporary
        if (index > 1) then
          if (val1 /= seismoTmp(index)) print *,'    kernel values:',val1,seismoTmp(index)
        endif
        seismoTmp(index) = val1
        seismoTmp(index+1) = val2
      enddo

      ! integral by spline representation: little bit better accuracy (influence on ~ 4. digit; neglectable), little bit slower
      !reverse adjoint
      do index = 1,numofTimeSteps/2
        val1 = seismoAdjoint(2,index)
        val2 = seismoAdjoint(2,numofTimeSteps-index+1)
        seismoAdjoint(2,index) = val2
        seismoAdjoint(2,numofTimeSteps-index+1) = val1
      enddo
      call getIntegral(seismo(2,:),seismoAdjoint(2,:),kernelVal,dt,numofTimeSteps)
    else
      ! central differences - scheme
      do index = 2,numofTimeSteps-1
        ! time window for integration is taken account of in building the adjoint source
        ! compute second derivative of forward seismogram
        for1 = seismo(2,index-1)
        for2 = seismo(2,index)
        for3 = seismo(2,index+1)
        ! second time derivative of forward signal
        val1 = (for1 + for3 - 2.0_WP*for2)/dt2

        ! adjoint wavefield index
        ! careful: the index of the adjoint wavefield corresponds to T-t,
        ! therefore (numofTimeSteps - index) for this second derivative
        indexadjoint = numofTimeSteps - index


        ! kernel value is the time integral (by a sum of discrete trapezoids with size dt)
        kernelVal = kernelVal + val1*seismoAdjoint(2,indexadjoint)*dt
      enddo
    endif

    ! calculate kernel factor for relative phase shift kernel (units in radians;
    ! sign convection depending on time lag definition)
    kernelfactor = 2.0_WP/(phaseVelocitySquare(vertex)*arrivalTime*vertexCellAreaRad)
    kernelVal = kernelfactor*kernelVal

    ! store in array
    adjointKernel(vertex) = kernelVal
  enddo

  end subroutine


!-----------------------------------------------------------------------
  subroutine storeAdjointKernel()
!-----------------------------------------------------------------------
! print the adjointKernel()-array values to a file 'adjointkernel.dat'
  use propagationStartup; use parallel; use adjointVariables;use verbosity; use cells
  implicit none
  ! local parameters
  real(WP):: lat,lon,sum_kern
  integer:: i,ier
  character(len=128):: kernelfile

  ! console output
  if (MAIN_PROCESS .and. VERBOSE) print *,'  writing values to kernel file...'

  ! get complete adjoint kernel array in case we run a single simulation on parallel processors
  if (PARALLELSEISMO) call collectAdjointKernel()

  ! remove obsolete files
  if (storeAsFile) call cleanupWaveFiles()

  ! slaves are done
  if (MAIN_PROCESS) then

    ! open kernel file
    kernelfile = trim(datadirectory)//trim(adjointKernelName)
    if (VERBOSE) then
      print *,'    storing kernel values in file  : ',trim(kernelfile)
    endif
    open(IOUT,file=trim(kernelfile),iostat=ier)
    if (ier /= 0) then
      print *,'Error: could not open '//trim(kernelfile)
      call stopProgram('Abort - storeAdjointKernel    ')
    endif

    ! file header
    write(IOUT,*) '# adjoint method - sensitivity kernel'
    write(IOUT,*) '# lon lat kernel vertexID'

    ! write values to file
    sum_kern = 0.0_WP
    do i = 1,numVertices
      call getSphericalCoord_Lat(i,lat,lon)
      write(IOUT,'(2f8.2,e18.6e3,i12)') lon,lat,adjointKernel(i),i

      ! summate values
      sum_kern = sum_kern+adjointKernel(i)*cellAreas(i)/EARTHRADIUS_SQUARED
    enddo
    close(IOUT)

    if (MAIN_PROCESS .and. VERBOSE) then
      print *,'    integrated over sphere         : ',sum_kern
      print *,'    kernel min/max                 : ',minval(adjointKernel(:)),'/',maxval(adjointKernel(:))
    endif
  endif

  ! wait until all processes reached this point
  call syncProcesses()

  end subroutine


!-----------------------------------------------------------------------
  subroutine getForwardWave(index,seismo)
!-----------------------------------------------------------------------
! get forward seismogram for at given vertex location from the wavefield storage-files
  use propagationStartup; use parallel
  implicit none
  integer,intent(in):: index
  real(WP),intent(out):: seismo(numofTimeSteps)
  ! local parameters
  character(len=8):: timestepstr
  character(len=3):: rankstr
  integer:: i,timestep,ier

  ! input
  write(rankstr,'(i3.3)') myrank
  i = 0
  do timestep = firsttimestep, lasttimestep
    i = i+1
    write(timestepstr,'(i8.7)')timestep
    !if (rank == 1) print *,'    open timestep:',i,timestep
    open(IIN,file=trim(datadirectory)//'wavefield_'//timestepstr//'.rank'//rankstr//'.dat', &
        access='direct',recl=WP,iostat=ier)
    if (ier /= 0) call stopProgram('could not open for input, file wavefield_'//&
        timestepstr//'.rank'//rankstr//'.dat ...    ')
    ! fill in displacements
    !if (rank == 1) print *,'    read timestep:',i,timestep
    read(IIN,rec=index,iostat=ier) seismo(i)
    if (ier /= 0) then
      print *,'Error: read error:',timestep,index,ier
      call stopProgram('could not read input from file wavefield_'//&
              timestepstr//'.rank'//rankstr//'.dat ...    ')
    endif
    !if (rank == 1) print *,'    close timestep:',i,timestep
    close(IIN)
  enddo

  end subroutine


!-----------------------------------------------------------------------
  subroutine setForwardWave(index,seismo)
!-----------------------------------------------------------------------
! set forward seismogram for at given vertex location to the wavefield storage-files
  use propagationStartup; use parallel
  implicit none
  integer,intent(in):: index
  real(WP),intent(in):: seismo(numofTimeSteps)
  ! local parameters
  character(len=8):: timestepstr
  character(len=3):: rankstr
  integer:: i,timestep,ier

  ! output
  write(rankstr,'(i3.3)') myrank
  i = 0
  do timestep = firsttimestep, lasttimestep
    i = i+1
    write(timestepstr,'(i8.7)')timestep
    open(IOUT,file=trim(datadirectory)//'wavefield_'//timestepstr//'.rank'//rankstr//'.dat', &
        access='direct',recl=WP,iostat=ier)
    if (ier /= 0) call stopProgram('could not open for output, file wavefield_'//&
        timestepstr//'.rank'//rankstr//'.dat ...    ')
    ! fill in displacements
    write(IOUT,rec=index,iostat=ier) seismo(i)
    if (ier /= 0) then
      print *,'Error: write error:',timestep,index,ier
      call stopProgram('could not write input from file wavefield_'//&
          timestepstr//'.rank'//rankstr//'.dat ...    ')
    endif
    close(IOUT)
  enddo

  end subroutine


!-----------------------------------------------------------------------
  subroutine getAdjointWave(index,seismo)
!-----------------------------------------------------------------------
! get adjoint seismogram for at given vertex location from the wavefield storage-files
  use propagationStartup; use parallel
  implicit none
  integer,intent(in):: index
  real(WP),intent(out)::seismo(numofTimeSteps)
  ! local parameters
  character(len=8):: timestepstr
  character(len=3):: rankstr
  integer:: i,timestep,ier

  ! input
  write(rankstr,'(i3.3)') myrank
  i = numofTimeSteps
  do timestep = firsttimestep, lasttimestep
    write(timestepstr,'(i8.7)')timestep
    open(IIN,file=trim(datadirectory)//'wavefieldAdj_'//timestepstr//'.rank'//rankstr//'.dat', &
        access='direct',recl=WP,iostat=ier)
    if (ier /= 0) call stopProgram('could not open for input, file wavefieldAdj_'//&
        timestepstr//'.rank'//rankstr//'.dat ...    ')
    ! fill in displacements
    read(IIN,rec=index,iostat=ier) seismo(i)
    if (ier /= 0) then
      print *,'Error: read error:',timestep,index,ier
      call stopProgram('could not read input from file wavefieldAdj_'//&
          timestepstr//'.rank'//rankstr//'.dat ...    ')
    endif
    close(IIN)
    ! adjoint wavefield: stores T-t meaning that firsttimestep value (t=0)
    ! is at seismo(T=numofTimeSteps) and lasttimestep at seismo(1)
    i = i-1
  enddo

  end subroutine


!-----------------------------------------------------------------------
  subroutine cleanupWaveFiles()
!-----------------------------------------------------------------------
! delete wavefield files
  use propagationStartup; use adjointVariables; use parallel
  implicit none
  ! local parameters
  character(len=8):: timestepstr
  character(len=3):: rankstr
  integer:: timestep,ier

  ! get rid of files
  write(rankstr,'(i3.3)') myrank
  do timestep = firsttimestep, lasttimestep
    write(timestepstr,'(i8.7)')timestep
    ! forward wavefield files
    open(IOUT,file=trim(datadirectory)//'wavefield_'//timestepstr//'.rank'//rankstr//'.dat', &
        access='direct',recl=WP,iostat=ier)
    if (ier == 0) close(IOUT,status='DELETE')
    ! adjoint wavefield files
    open(IOUT,file=trim(datadirectory)//'wavefieldAdj_'//timestepstr//'.rank'//rankstr//'.dat', &
        access='direct',recl=WP,iostat=ier)
    if (ier == 0) close(IOUT,status='DELETE')
  enddo

  end subroutine

!-----------------------------------------------------------------------
  subroutine checkKernel()
!-----------------------------------------------------------------------
! calculates the integral over the sphere of the kernel
! this value should have a modulus close to 1
  use propagationStartup; use parallel; use adjointVariables;use verbosity; use cells
  implicit none
  ! local parameters
  integer:: i
  real(WP):: sum_kern
  real(WP),parameter:: EPS = 1.5

  ! get complete adjoint kernel array in case we run a single simulation on parallel processors
  if (PARALLELSEISMO) call collectAdjointKernel()

  ! slaves are done
  if (.not. MAIN_PROCESS) return

  ! summate values
  sum_kern = 0.0_WP
  do i = 1,numVertices
    sum_kern = sum_kern+adjointKernel(i)*cellAreas(i)/EARTHRADIUS_SQUARED
  enddo

  ! check
  if (abs(sum_kern) > EPS) then
    print *
    print *,'Error: kernel integral over sphere:',sum_kern
    call stopProgram('kernel integral error ')
  endif

  ! console output
  print *,'      kernel integral:',sum_kern

  end subroutine

!-----------------------------------------------------------------------
  subroutine determineWindowsize(startindex,endindex,range)
!-----------------------------------------------------------------------
! if halfwidth_estrom is chosen, then it is
! based upon the window size mentioned in Ekstrom et al., 1997: eq. 16 and
! denoted parameter J_i which is one half of the window:
! the total window includes three cycles of the center period
  use precisions; use parallel; use verbosity
  use propagationStartup
  use filtertype, only: bw_waveperiod
  use adjointVariables
  implicit none
  integer,intent(out):: startindex,endindex,range
  ! local parameters
  real(WP):: halfwidth,max
  real(WP):: arrivalTime,distance
  integer,dimension(1):: loc
  integer:: center,i,default_start,default_end,iorbit
  logical,parameter:: HALFWIDTH_EKSTROM = .false.

  ! determines which is the default window's start index
  loc = maxloc( seismogramReceiver(1,:), MASK = seismogramReceiver(1,:) < WINDOW_START )
  default_start = loc(1)
  if (default_start < 1 .or. default_start > numofTimeSteps) default_start = 1
  ! determines which is the default window's end index
  loc = minloc( seismogramReceiver(1,:), MASK = seismogramReceiver(1,:) > WINDOW_END )
  default_end = loc(1)
  if (default_end < 1 .or. default_end > numofTimeSteps) default_end = numofTimeSteps
  if (default_start > default_end) default_end = default_start

  ! this takes the reference travel time as center of the window
  if (HALFWIDTH_EKSTROM) then
    call determineOrbitDistance(distance,iorbit)
    arrivalTime = distance*EARTHRADIUS/cphaseRef
    loc = maxloc( seismogramReceiver(1,default_start:default_end), &
              MASK = seismogramReceiver(1,default_start:default_end) < arrivalTime )
    center = default_start + loc(1) ! take next higher step index
    max = seismogramReceiver(2,center)
  else
    ! alternative:
    ! window centers around the time of maximum of the default seismogram
    !loc = maxloc( abs(seismogramReceiver(2,default_start:default_end)) )
    !center = default_start + loc(1) - 1

    ! center of window_start and window_end
    loc = maxloc( seismogramReceiver(1,default_start:default_end), &
              MASK = seismogramReceiver(1,default_start:default_end) < (WINDOW_START+WINDOW_END)/2.0 )
    center = default_start + loc(1) ! take next higher step index
  endif

  ! half of the window (total window includes three cycles)
  if (HALFWIDTH_EKSTROM) then
    halfwidth = 3.0 * bw_waveperiod / 2.0
  else
    halfwidth = (WINDOW_END - WINDOW_START)/2.0
  endif

  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'  determining cross-correlation window...'
    print *,'    center time (s) / value            : ',seismogramReceiver(1,center),max
    print *,'    reference arrival time (s) / orbit : ',arrivalTime,iorbit
    print *,'    half width (s)                     : ',halfwidth
    print *,'    default window start/end           : ',WINDOW_START,WINDOW_END
    print *,'    default index start/end            : ',default_start,default_end
    if (WINDOW_END < arrivalTime) print *,'  default WINDOW_END may be too small !'
  endif

  ! centers window around maximum
  WINDOW_START = seismogramReceiver(1,center) - halfwidth
  WINDOW_END = seismogramReceiver(1,center) + halfwidth

  ! initializes
  startindex = 0
  endindex = 1

  ! time window sets adjoint signal to zero outside of window
  ! (note: adjointSource is reversed in time)
  do i = 1,numofTimeSteps
    if (adjointSource(1,i) > WINDOW_END) then
      adjointSource(2,i) = 0.0
      endindex = i  ! points to where end begins
    else if (adjointSource(1,i) < WINDOW_START) then
      adjointSource(2,i) = 0.0
      if (startindex == 0) startindex = i ! points to where start ends :)
    endif
  enddo
  WINDOW_START = adjointSource(1,startindex)
  WINDOW_END   = adjointSource(1,endindex)

  ! console output
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'  adjoint source window (in s):'
    print *,'    start              : ',WINDOW_START
    print *,'    end                : ',WINDOW_END
    print *,'    center/halfwidth    : ',seismogramReceiver(1,center),halfwidth
  endif

  ! checks window
  range = startindex - endindex + 1
  if (range < 1) then
      print *,'Error: no window for windowed integration : ',myrank
      print *,'       window start/end                   : ',WINDOW_START,WINDOW_END
      print *,'       adjoint indices start/end             : ',startindex,endindex
      print *,'       number of timesteps                : ',numofTimeSteps
      call stopProgram('error - crosscorrelation window')
  endif

  ! inform for different windows, in cases
  ! for large epicentral distances and simulation time set to antipode time
  if (abs(seismogramReceiver(1,center)-WINDOW_START) < halfwidth &
     .or. abs(seismogramReceiver(1,center)-WINDOW_END) < halfwidth) then
    if (MAIN_PROCESS) then
      print *,'  adjoint window for cross-correlation shortened:'
      print *,'    window start/end       : ',WINDOW_START,WINDOW_END
      print *,'    adjoint indices start/end : ',startindex,endindex
      print *,'    simulation start/end     : ',adjointSource(1,numofTimeSteps),adjointSource(1,1)
      !call stopProgram("window too small    ")
    endif
  endif

  end subroutine

!-----------------------------------------------------------------------
  subroutine printAdjointWavefield(timestep,time)
!-----------------------------------------------------------------------
! prints complete newdisplacement wavefield to file
!
! returns: wavefiled in file #output-dir#/simulation.#time#.dat
  use precisions
  use cells, only: vertices,numVertices
  use parallel, only: MAIN_PROCESS
  use displacements, only: newdisplacement
  use propagationStartup, only: datadirectory,IOUT
  implicit none
  integer,intent(in):: timestep
  real(WP),intent(in):: time
  ! local parameters
  character(len=4):: timestr
  integer:: n,k

  ! only main process writes to files
  if (mod(timestep,10) == 0 .and. time >= 0.0 .and. MAIN_PROCESS) then
    write(timestr,'(i4.4)') timestep
    open(IOUT,file=trim(datadirectory)//'simulationAdjoint.'//timestr//'.dat')
    do n = 1, numVertices
      write(IOUT,'(4f18.6)')(vertices(n,k),k=1,3),newdisplacement(n)
    enddo
    close(IOUT)
    print *,'    file written: ',trim(datadirectory)//'simulationAdjoint.'//timestr//'.dat'
  endif

  end subroutine

