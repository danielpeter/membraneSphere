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
!
! Proper acknowledgement shall be made to the authors of the software in publications and
! presentations resulting from the use of this software:
!
! Peter, D., C. Tape, L. Boschi and J. H. Woodhouse, 2007. Surface wave tomography:
! global membrane waves and adjoint methods, Geophys. J. Int., , 171: p. 1098-1117.


! heterogeneousPhaseshift.f90
!
! based upon heterogeneousInversion.f90 & phaseshift.f90

!-----------------------------------------------------------------------
  module heterogeneousArrays
!-----------------------------------------------------------------------
    use precisions
    real(WP), allocatable, dimension(:):: phaseVelocitySquarePREM,phaseVelocitySquareHET
    real, allocatable, dimension(:,:):: dataStore
    real:: shift_min,shift_max
  end module

!-----------------------------------------------------------------------
  program heterogeneousPhaseshift
!-----------------------------------------------------------------------
! calculates the new synthetic data which can be used for a benchmark exercise
! the synthetic phase shifts for a heterogeneous background earth are plotted as
! a new data file similar as e.g. wei_sum.2.L0150.1.txt
!
! performs a simulation, with a reference earth model (PREM), and calculates the
! phase shift to a second simulation within a heterogeneous background earth
!
! adjoint method: information in tromp et al. (2005)
! finite-difference iteration: information in carl tape thesis (2003), chap 5, (5.7)
    use propagationStartup; use parallel
    use minimize_trace, only: trace_length
    implicit none
    integer,parameter:: dataUnit = 81, newdataUnit = 82
    integer:: datacount

    !-----------------------------------------------------------------------
    ! parameters
    ! most parameters concerning wave propagation set in file
    ! Parameter_Input (&default values in commonModules.f90/initialize.f90)
    HetPhaseshift_Program  = .true.
    rotate_frame           = .false.  ! no rotation of phase map
    !-----------------------------------------------------------------------

    ! initialization of parameters and arrays
    call initialize()

    ! initialize for simplex routine
    trace_length = 0

    ! open file with data
    call openDatafile(dataUnit,newdataUnit,datacount)

    ! store data in array
    call syncDataArray(dataUnit,datacount)

    ! handle events
    call processDatafile(newdataUnit,datacount)

    ! end smoothly
    call finish()
  end


  !-----------------------------------------------------------------------
  subroutine openDatafile(dataUnit,newdataUnit,datacount)
  !-----------------------------------------------------------------------
  ! opens file with given real data (only master process)
  ! (e.g., ETL97, phase anomaly measurements, wei_sum.1.L150.txt )
    use parallel; use phaseBlockData; use propagationStartup
    use filtertype, only: bw_waveperiod
    implicit none
    integer,intent(in):: dataUnit,newdataUnit
    integer,intent(out):: datacount
    integer:: stationscount,eventscount,i,ierror,slashindex
    real:: stla_old,stlo_old,epla_old,eplo_old,epidelta
    real:: stations(10000,2), events(20000,2)
    logical:: stationfound,eventfound
    character:: commentline*47,wave_type*1
    real:: period,v0,omega
    real:: eplo,epla,stla,stlo,datum,error
    real:: epicentraldistance_min,epicentraldistance_max,datummax,datummin
    integer:: epiminIndex,epimaxIndex
    character(len=256):: nameout

    ! checks that only master is going to do it
    if (.not. MASTER ) return

    ! open data file
    print *,'  open datafile: '
    print *,'      ',heterogeneousDataFile(1:len_trim(heterogeneousDataFile))
    open(dataUnit,file=heterogeneousDataFile(1:len_trim(heterogeneousDataFile)), &
            status="old",iostat=ierror)
    if ( ierror /= 0 ) then
      print *,'Error: process ',rank,' could not open data file: ', &
                heterogeneousDataFile(1:len_trim(heterogeneousDataFile))
      call stopProgram( "data file not found     ")
    endif

    !------read its header
    !------header MUST contain period (in s) at which measurements were made,
    !------and reference velocity of that type of sw at that period
    !read(2,1003) wave_type,period,v0
    read(dataUnit,*) wave_type,period,v0
    print *,'    wave type, period, v0:',wave_type,period,v0
    read(dataUnit,'(a47)') commentline

    omega = 1./period

 1003 format(8x,a1,7x,f8.4,10x,f8.5)
 1002 format(4(f11.4,1x),1x,f11.6,1x,f11.7)

    ! count data entries
    stla_old = 0.
    stlo_old = 0.
    epla_old = 0.
    eplo_old = 0.
    stations(:,:)=0.
    events(:,:)=0.

    epicentraldistance_min=1000.0 ! modified d.peter (8.2.2006)
    epicentraldistance_max = 0.0
    datummax = 0.0
    datummin = 0.0


    ierror = 0
    datacount = 0
    stationscount = 0
    eventscount = 0
    do while( ierror == 0 )
      !read(2,1002,iostat=ierror) epla,eplo,stla,stlo,datum,error
      read(dataUnit,*,iostat=ierror) epla,eplo,stla,stlo,datum,error
      if ( ierror == 0 ) then
        datacount = datacount + 1

        ! receivers: file is (normally) ordered by receiver location
        stationfound = .false.
        do i = 1,stationscount
          if ( abs(stla -stations(i,1)) < 1.e-4 .and. &
             abs(stlo - stations(i,2)) < 1.e-4 ) then
            stationfound = .true.
            exit
          endif
        enddo
        if (.not. stationfound ) then
          stationscount = stationscount+1
          stations(stationscount,1)=stla
          stations(stationscount,2)=stlo
        endif


        ! events/sources: file is (normally) ordered by receiver location,
        !             the stations can have different number of events
        eventfound = .false.
        do i = 1,eventscount
          if ( abs(epla-events(i,1)) < 1.e-4 .and. &
             abs(eplo-events(i,2)) < 1.e-4) then
            eventfound = .true.
            exit
          endif
        enddo
        if (.not. eventfound ) then
          eventscount = eventscount+1
          events(eventscount,1)=epla
          events(eventscount,2)=eplo
        endif

        ! some statistics
        if ( datum > datummax ) datummax = datum
        if ( datum < datummin ) datummin = datum

        ! find epicentral distances between source and station
        call epicentralDistance(epla,eplo,stla,stlo,epidelta)

        ! store min/max of great circle distance between epicenter and
        ! station location (modified d.peter)
        if (epidelta < epicentraldistance_min) then
          epicentraldistance_min = epidelta
          epiminIndex = datacount
        endif
        if (epidelta > epicentraldistance_max) then
          epicentraldistance_max = epidelta
          epimaxIndex = datacount
        endif
      endif
    enddo

    ! reset to first data entry
    !rewind(dataUnit)
    !read(dataUnit,*)
    !read(dataUnit,*)

    ! console output
    print *,'  data:'
    print *,'    number of stations:',stationscount
    print *,'    number of events:',eventscount
    print *,'    great circle distances:'
    print *,'      minimum = ',epicentraldistance_min,epiminIndex
    print *,'      maximum = ',epicentraldistance_max,epimaxIndex
    print *,'    datum data:'
    print *,'      min = ',datummin
    print *,'      max = ',datummax

    ! open new output data file
    slashindex = 1
    do i = 1,len_trim(heterogeneousDataFile)
      if ( heterogeneousDataFile(i:i) == '/') slashindex = i
    enddo
    nameout = datadirectory(1:len_trim(datadirectory))// &
        heterogeneousDataFile(slashindex+1:len_trim(heterogeneousDataFile))&
        //'.new'
    open(newdataUnit,file=trim(nameout),iostat=ierror)
    if ( ierror /= 0 ) &
      call stopProgram("error opening new data file .new for output   ")

    print *,'    new datum file will be:'
    print *,'      ',trim(nameout)

    ! write first two line
    write(newdataUnit,1003) wave_type,bw_waveperiod,cphaseRef
    write(newdataUnit,'(a47)') commentline
    call myflush(newdataUnit)

    ! info
    print *,'  data read done'
    call myflush(6)

  end subroutine


  !-----------------------------------------------------------------------
  subroutine syncDataArray(dataUnit,datacount)
  !-----------------------------------------------------------------------
  ! stores the data in an array
    use parallel; use phaseBlockData; use heterogeneousArrays
    implicit none
    integer,intent(in):: dataUnit,datacount
    integer:: stationscount,eventscount,i,ierror,length
    real:: stla_old,stlo_old,epla_old,eplo_old,epidelta
    real:: stations(10000,2), events(20000,2)
    logical:: stationfound,eventfound
    character:: commentline*47,wave_type*1
    real:: period,v0,omega
    real:: eplo,epla,stla,stlo,datum,error

    ! info
    if ( MASTER ) then
      print *,'  data count:',datacount
    endif

    ! wait until all processes reached this point
    call MPI_Barrier( MPI_COMM_WORLD, ierror )
    if ( ierror /= 0) call stopProgram('abort syncDataArray - MPI_Barrier failed    ')

    ! synchronize datacount with processes
    if ( MASTER ) then
      do i = 1,nprocesses-1
        call MPI_Send(datacount,1,MPI_INTEGER,i,i,MPI_COMM_WORLD,ierror)
        if ( ierror /= 0 ) then
          print *,'Error: datacount',rank,ierror
          call stopProgram("syncDataArray()")
        endif
      enddo
    else
      call MPI_Recv(datacount,1,MPI_INTEGER,0,rank,MPI_COMM_WORLD,status,ierror)
      if ( ierror /= 0 ) then
        print *,'Error: datacount',rank,ierror
        call stopProgram("syncDataArray()")
      endif
    endif

    ! allocate array
    allocate(dataStore(6,datacount),stat=ierror)
    if ( ierror /= 0 ) call stopProgram('error allocating data - syncDataArray()')

    ! process data file
    if ( MASTER ) then
      ! reset to first data entry
      rewind(dataUnit)
      read(dataUnit,*) wave_type,period,v0
      read(dataUnit,'(a47)') commentline

      ! read data
      print *,'reading in data array...'
      print *,'    number of data entries:',datacount
      do i = 1,datacount
        read(dataUnit,*,iostat=ierror) epla,eplo,stla,stlo,datum,error
        if ( ierror /= 0 ) call stopProgram('error reading data - syncDataArray()')

        ! store in data array
        dataStore(1,i) = epla
        dataStore(2,i) = eplo
        dataStore(3,i) = stla
        dataStore(4,i) = stlo
        dataStore(5,i) = datum
        dataStore(6,i) = error
      enddo

      ! close file
      close(dataUnit)
    endif


    ! synchronize datacount with processes
    length = 6*datacount
    if ( MASTER ) then
      print *,'    synchronize data ',8*length/1024./1024.,' MB'
      call myflush(6)
      do i = 1,nprocesses-1
        call MPI_Send(dataStore(:,:),length,MPI_REAL,i,i,MPI_COMM_WORLD,ierror)
        if ( ierror /= 0 ) then
          print *,'Error: datastore',rank,ierror
          call stopProgram("syncDataArray()")
        endif
      enddo
    else
      call MPI_Recv(dataStore(:,:),length,MPI_REAL,0,rank,MPI_COMM_WORLD,status,ierror)
      if ( ierror /= 0 ) then
        print *,'Error: datastore',rank,ierror
        call stopProgram("syncDataArray()")
      endif
    endif

    ! TEST console output
    if ( MASTER) then
      print *,'    first datum extract:',(dataStore(i,1),i=1,6)
      call myflush(6)
    endif

    ! wait until all processes reached this point
    call MPI_Barrier( MPI_COMM_WORLD, ierror )
    if ( ierror /= 0) call stopProgram('abort syncDataArray - MPI_Barrier failed    ')

  end subroutine syncDataArray

  !-----------------------------------------------------------------------
  subroutine processDatafile(newdataUnit,datacount)
  !-----------------------------------------------------------------------
  ! handle events
    use parallel; use propagationStartup; use griddomain; use heterogeneousArrays
    use verbosity
    implicit none
    integer,intent(in):: newdataUnit,datacount
    integer:: i,j,n,vertex,dataLoopIndex
    logical:: doCalc
    real:: benchmarkLoopStart,benchmarkLoopEnd,epidelta
    character(len=5)::kernelstr
    integer:: isqre,ierror
    real:: eplo,epla,stla,stlo,datum,error
    real:: time_shift
    external:: isqre,doCalc

    ! prepare phasevelocities for time iterations
    call constructPhaseVelocities()

    ! loop over data
    if ( MASTER) then
      print *,'begin loop over data:',datacount
      call myflush(6)
    endif

    shift_min = 0.0
    shift_max = 0.0
    dataLoopIndex = 0
    do j = 1,datacount
      ! console output
      if ( MASTER .and. ( mod(j,500) == 0 .or. j == 1 )  ) then
        VERBOSE = .true.
        print *
        print *,j," data read"
        call myflush(6)
      else
        VERBOSE = .false.
      endif

      ! check if this process has to do something
      if ( doCalc(j) ) then
        ! take data from array
        epla = dataStore(1,j)
        eplo = dataStore(2,j)
        stla = dataStore(3,j)
        stlo = dataStore(4,j)
        datum = dataStore(5,j)
        error = dataStore(6,j)

        ! benchmark
        if ( MASTER .and. mod(j,500) < nprocesses ) then
          benchmarkLoopStart = MPI_WTIME()
        endif

        ! prepare source/receiver locations for simulation
        call prepareCouple(epla,eplo,stla,stlo)

        ! do the time iteration and phase shift calculation
        call doTimeIteration(j,time_shift)

        ! slave processes transfers data to MASTER which prints them to file
        call outputDataLine(j,datacount,epla,eplo,stla,stlo,time_shift, &
                            error,newdataUnit)

        ! output to file
        if ( MASTER .and. mod(j,500) < nprocesses ) then
          ! find epicentral distances between source and station
          call epicentralDistance(epla,eplo,stla,stlo,epidelta)

          ! console output
          print *,'    data:',j
          print *,'      source:',sourceLat,sourceLon
          print *,'        (original):',epla,eplo
          print *,'      receiver:',receiverLat,receiverLon
          print *,'        (original):',stla,stlo
          print *,'      distance source-receiver:',epidelta
          print *,'      time shift:',time_shift
          ! benchmark for this run
          benchmarkLoopEnd = MPI_WTIME()
          print *,'      benchmark seconds:',benchmarkLoopEnd-benchmarkLoopStart
          ! flush to console output (i/o unit nr. 6 or on SGI system nr. 101)
          call myflush(6)
        endif

        ! for possible exits in a job queue
        if ( j > dataLoopIndex) dataLoopIndex = j
      endif
    enddo

    ! wait until all processes reached this point
    call MPI_Barrier( MPI_COMM_WORLD, ierror )
    if ( ierror /= 0) call stopProgram('abort - final MPI_Barrier failed    ')

    ! close new data file
    if ( MASTER ) then
      close(newdataUnit)
    endif

    ! console output
    if ( MASTER ) then
      VERBOSE = .true.
      print *,'data written:',datacount
      call myflush(6)
    endif
  end subroutine processDatafile


  !-----------------------------------------------------------------------
  subroutine constructPhaseVelocities()
  !-----------------------------------------------------------------------
  ! determine phase velocity array for PREM and heterogeneous run
    use propagationStartup; use phaseVelocityMap; use parallel; use verbosity
    use heterogeneousArrays; use cells
    implicit none
    integer:: ierror

    ! allocate memory
    allocate(phaseVelocitySquarePREM(numVertices), &
            phaseVelocitySquareHET(numVertices),stat=ierror)
    if ( ierror /= 0 ) stop "allocate phasevelocities error"

    ! PREM
    if ( MASTER .and. VERBOSE ) print *,' constructing reference phase map...'
    HETEROGENEOUS = .false. ! homogeneous
    DELTA         = .false. ! no delta scatterer
    call constructPhaseVelocitySquare()
    phaseVelocitySquarePREM(:) = phaseVelocitySquare(:)

    ! heterogeneous phase speed map
    if ( MASTER .and. VERBOSE ) print *,' constructing heterogeneous phase map...'
    HETEROGENEOUS = .true.
    rotate_frame  = .false.  ! no rotation of phase map
    DELTA         = .false.  ! no delta scatterer

    !allocate phaseMap array
    if ( allocated(phaseMap) ) deallocate(phaseMap)
    allocate(phaseMap(numVertices), stat=ierror)
    if ( ierror /= 0) stop 'abort - phase map allocation'

    ! read phase map (e.g. Love 150 s) (note: is rotated for non-adjoint simulations)
    if ( MASTER ) then
      call constructPhasedata()
    endif

    ! synchronize between master and slave processes
    call syncPhaseMap(rank,nprocesses)

    ! get phase velocities
    call constructPhaseVelocitySquare()
    phaseVelocitySquareHET(:) = phaseVelocitySquare(:)

    ! phase map no more needed
    if ( HETEROGENEOUS ) deallocate(phaseMap)

    if ( MASTER .and. VERBOSE) print *

  end subroutine constructPhaseVelocities

  !-----------------------------------------------------------------------
  logical function doCalc(index)
  !-----------------------------------------------------------------------
  ! determines wheter the data entry (index) should be proceeded
  ! and calculations should be done
    use parallel
    implicit none
    integer:: index

    ! default
    doCalc = .true.

    ! parallel processes
    if ( nprocesses < 2 ) return

    ! distribute following calculation over all processes
    if (.not. PARALLELSEISMO ) then
      ! master: j=1,9,17,...; rank 1: j=2,10,...; rank 2: j=3,11,...
      if ( mod((index-rank+2*nprocesses-1),nprocesses) == 0 ) then
        doCalc = .true.
      else
        doCalc = .false.
      endif
    else
      doCalc = .true.
    endif

    return
  end function doCalc


  !-----------------------------------------------------------------------
  subroutine doTimeIteration(jrecord,time_shift)
  !-----------------------------------------------------------------------
    use propagationStartup; use phaseVelocityMap
    use traveltime; use cells; use verbosity; use parallel
    use heterogeneousArrays
    implicit none
    integer,intent(in):: jrecord
    real,intent(out):: time_shift
    real(WP),dimension(2,numofTimeSteps):: seismogramPREM,seismogramHET
    real(WP):: startTime,distance,amplification
    integer:: ierror,i
    character:: jstr*6

    ! be sure only one receiver at a time
    manyReceivers = .false.

    ! PREM reference simulation
    phaseVelocitySquare(:) = phaseVelocitySquarePREM(:)
    call forwardIteration()
    seismogramPREM(:,:) = seismogramReceiver(:,:)
    ! file output
    if ( MASTER .and. mod(jrecord,500) < nprocesses ) then
      write(jstr,'(i6.6)')jrecord
      open(20,file=datadirectory(1:len_trim(datadirectory))//&
                              'seismoPREM.'//jstr//'.txt')
      do i = 1,numofTimeSteps
        write(20,*)seismogramPREM(1,i), seismogramPREM(2,i)
      enddo
      close(20)
    endif

    ! heterogeneous simulation
    phaseVelocitySquare(:) = phaseVelocitySquareHET(:)
    call forwardIteration()
    seismogramHET(:,:) = seismogramReceiver(:,:)
    ! file output
    if ( MASTER .and. mod(jrecord,500) < nprocesses ) then
      open(20,file=datadirectory(1:len_trim(datadirectory))//&
                              'seismoHET.'//jstr//'.txt')
      do i = 1,numofTimeSteps
        write(20,*)seismogramHET(1,i), seismogramHET(2,i)
      enddo
      close(20)
    endif

    ! determine phase shift:
    ! determine startTime
    call getStartTimeSeismogram(seismogramPREM,numofTimeSteps, &
                          seismogramPREM(1,numofTimeSteps),dt,startTime)

    ! calculate the analytical arrival time
    call greatCircleDistance(vertices(sourceVertex,:), &
                             vertices(receiverVertex,:),distance)
    arrivalTime = distance*EARTHRADIUS/cphaseRef

    ! console output
    if ( MASTER .and. VERBOSE ) then
      print *,'    arrival time:',arrivalTime
      print *,'    start time for fourier transformation:',startTime
    endif

    ! calculate time lag (in seconds, when dt in seconds)
    call getTimelagSeismos(t_lag,seismogramHET,seismogramPREM,numofTimeSteps,dt)

    ! find Euler angles and rotation matrix (code by John W.)
    !call pole(epla,eplo,stla,stlo,xlatp,xlonp,azmp,delta)
    !call epicentralDistance(epla,eplo,stla,stlo,epidelta)

    ! t0: traveltime (s) with factor PI/180 * 6371 km ~ 111.2
    !t0=epidelta*111.19492664/v0
    ! relative phase shift = relative traveltime anomaly
    !dphi_phi=datum/t0
    !print *,'    relative traveltime = ',dphi_phi,datum,t0
    !      dphi_phi=datum            !to invert for absolute delta_p
    if ( MASTER .and. VERBOSE ) then
      print *,'    cross-correlation lag:',t_lag
    endif

    ! downhill simplex
    call getTimelagSimplex(t_lag,amplification,seismogramHET,seismogramPREM,numofTimeSteps)

    ! set phase shift
    time_shift = t_lag

    ! console output
    if ( MASTER .and. VERBOSE ) then
      print *,'    time shift:',time_shift
      call myflush(6)
    endif
  end subroutine doTimeIteration

  !-----------------------------------------------------------------------
  subroutine outputDataLine(jrecord,datacount,epla,eplo,stla,stlo,shift, &
                            error,newdataUnit)
  !-----------------------------------------------------------------------
  ! sync with master process and print to file
    use parallel
    implicit none
    integer,intent(in):: jrecord,datacount,newdataUnit
    real,intent(in):: shift,epla,eplo,stla,stlo,error
    integer:: restprocesses,i,ierror
    real:: epla_proc,eplo_proc,stla_proc,stlo_proc,shift_proc,error_proc

    ! single processor/single simulation job
    if ( nprocesses < 2 .or. PARALLELSEISMO) then
      call printDataLine(newdataUnit,epla,eplo,stla,stlo,shift,error)
      return
    endif

    ! more processes: get phaseshift data and print to file
    if ( MASTER ) then
      ! print master info first
      call printDataLine(newdataUnit,epla,eplo,stla,stlo,shift,error)

      ! how many process are there to collect
      restprocesses = datacount - jrecord
      if ( restprocesses >= nprocesses ) then
        restprocesses = nprocesses - 1
      endif

      ! check if something to receive
      if ( restprocesses == 0 ) return

      ! MASTER receives
      do i = 1,restprocesses
        ! get data
        !print *,'getting from',i
        call MPI_Recv(shift_proc,1,MPI_REAL,i,i,MPI_COMM_WORLD,status,ierror)
        if ( ierror /= 0 ) then
          print *,'Error: shift',rank,ierror
          call stopProgram("outputDataline    ")
        endif
        call MPI_Recv(epla_proc,1,MPI_REAL,i,i,MPI_COMM_WORLD,status,ierror)
        if ( ierror /= 0 ) then
          print *,'Error: shift',rank,ierror
          call stopProgram("outputDataline    ")
        endif
        call MPI_Recv(eplo_proc,1,MPI_REAL,i,i,MPI_COMM_WORLD,status,ierror)
        if ( ierror /= 0 ) then
          print *,'Error: shift',rank,ierror
          call stopProgram("outputDataline    ")
        endif
        call MPI_Recv(stla_proc,1,MPI_REAL,i,i,MPI_COMM_WORLD,status,ierror)
        if ( ierror /= 0 ) then
          print *,'Error: shift',rank,ierror
          call stopProgram("outputDataline    ")
        endif
        call MPI_Recv(stlo_proc,1,MPI_REAL,i,i,MPI_COMM_WORLD,status,ierror)
        if ( ierror /= 0 ) then
          print *,'Error: shift',rank,ierror
          call stopProgram("outputDataline    ")
        endif
        call MPI_Recv(error_proc,1,MPI_REAL,i,i,MPI_COMM_WORLD,status,ierror)
        if ( ierror /= 0 ) then
          print *,'Error: shift',rank,ierror
          call stopProgram("outputDataline    ")
        endif

        ! print to file
        call printDataLine(newdataUnit,epla_proc,eplo_proc,stla_proc,stlo_proc, &
                          shift_proc,error_proc)
      enddo
    else
      ! slave process sends information
      call MPI_Send(shift,1,MPI_REAL,0,rank,MPI_COMM_WORLD,ierror)
      if ( ierror /= 0 ) then
        print *,'Error: shift',rank,ierror
        call stopProgram("outputDataline    ")
      endif
      call MPI_Send(epla,1,MPI_REAL,0,rank,MPI_COMM_WORLD,ierror)
      if ( ierror /= 0 ) then
        print *,'Error: shift',rank,ierror
        call stopProgram("outputDataline    ")
      endif
      call MPI_Send(eplo,1,MPI_REAL,0,rank,MPI_COMM_WORLD,ierror)
      if ( ierror /= 0 ) then
        print *,'Error: shift',rank,ierror
        call stopProgram("outputDataline    ")
      endif
      call MPI_Send(stla,1,MPI_REAL,0,rank,MPI_COMM_WORLD,ierror)
      if ( ierror /= 0 ) then
        print *,'Error: shift',rank,ierror
        call stopProgram("outputDataline    ")
      endif
      call MPI_Send(stlo,1,MPI_REAL,0,rank,MPI_COMM_WORLD,ierror)
      if ( ierror /= 0 ) then
        print *,'Error: shift',rank,ierror
        call stopProgram("outputDataline    ")
      endif
      call MPI_Send(error,1,MPI_REAL,0,rank,MPI_COMM_WORLD,ierror)
      if ( ierror /= 0 ) then
        print *,'Error: shift',rank,ierror
        call stopProgram("outputDataline    ")
      endif
    endif
  end subroutine outputDataLine

  !-----------------------------------------------------------------------
  subroutine printDataLine(newdataUnit,epla,eplo,stla,stlo,shift,error)
  !-----------------------------------------------------------------------
  ! writes given data to file
    use heterogeneousArrays
    implicit none
    integer:: newdataUnit
    real:: epla,eplo,stla,stlo,shift,error

    ! new data file output
    write(newdataUnit,1002) epla,eplo,stla,stlo,shift,error
 1002 format(4(f11.4,1x),1x,f11.6,1x,f11.7)

    call myflush(newdataUnit)

    ! statistics
    if ( shift < shift_min ) shift_min = shift
    if ( shift > shift_max ) shift_max = shift

  end subroutine printDataLine

  !-----------------------------------------------------------------------
  subroutine myflush(i)
  !-----------------------------------------------------------------------
  ! should write out buffered data
    implicit none
    integer, intent(in) :: i

    ! flush to console output (i/o unit nr. 6 or on SGI system nr. 101)
    ! if ( i == 6 ) i = 101

    !  #ifdef _IBM_
      !call flush_(i)
    !#else
      call FLUSH(i)
    !#endif
  end subroutine myflush

  !-----------------------------------------------------------------------
  subroutine finish()
  !-----------------------------------------------------------------------
  ! finalize MPI communication
    use parallel; use propagationStartup; use heterogeneousArrays
    implicit none
    integer:: ierror

    ! benchmark
    if ( MASTER ) then
      benchAllEnd = MPI_WTIME()
      print *,'statistics:'
      print *,'  datum min/max:',shift_min,shift_max
      print *
      print *,'number of processors:',nprocesses
      print *,'running time: ',int((benchAllEnd-benchAllStart)/60.0),'min ', &
                                mod((benchAllEnd-benchAllStart),60.0),'sec'
      print *
    endif

    ! wait until all processes reached this point
    call MPI_Barrier( MPI_COMM_WORLD, ierror )
    if ( ierror /= 0) call stopProgram('abort - final MPI_Barrier failed    ')

    ! end parallelization
    call MPI_FINALIZE(ierror)
    if ( ierror /= 0) call stopProgram('abort - finalize failed    ')

  end subroutine finish
