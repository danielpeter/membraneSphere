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
  implicit none
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
! adjoint method: information in Tromp et al. (2005)
! finite-difference iteration: information in Carl Tape thesis (2003), chap 5, (5.7)
  use propagationStartup; use parallel
  use minimize_trace, only: trace_length
  implicit none
  ! local parameters
  integer,parameter:: dataUnit = 81, newdataUnit = 82
  integer:: datacount

  !-----------------------------------------------------------------------
  ! parameters
  ! most parameters concerning wave propagation set in file
  ! Parameter_Input (&default values in commonModules.f90/initialize.f90)
  HetPhaseshift_Program  = .true.
  !-----------------------------------------------------------------------

  ! initialization of parameters and arrays
  call initialize()

  ! initialize for simplex routine
  trace_length = 0

  ! open file with data
  call openDatafile(dataUnit,newdataUnit,datacount)

  ! store data in array
  call storeHeterogeneousDataArray(dataUnit,datacount)

  ! handle events
  call processDatafile(newdataUnit,datacount)

  ! end smoothly
  call finish()

  end program


!-----------------------------------------------------------------------
  subroutine openDatafile(dataUnit,newdataUnit,datacount)
!-----------------------------------------------------------------------
! opens file with given real data (only main process)
! (e.g., ETL97, phase anomaly measurements, wei_sum.1.L150.txt )
  use parallel; use phaseBlockData; use propagationStartup
  use filtertype, only: bw_waveperiod
  implicit none
  integer,intent(in):: dataUnit,newdataUnit
  integer,intent(out):: datacount
  ! local parameters
  integer:: stationscount,eventscount,i,ier,slashindex
  real:: stla_old,stlo_old,epla_old,eplo_old,epidelta
  real,dimension(:,:),allocatable:: stations, events
  logical:: stationfound,eventfound
  character(len=47):: commentline
  character(len=1):: wave_type
  real:: period,v0,omega
  real:: eplo,epla,stla,stlo,datum,error
  real:: epicentraldistance_min,epicentraldistance_max,datummax,datummin
  integer:: epiminIndex,epimaxIndex
  character(len=256):: nameout

  ! checks that only main process is going to do it
  if (.not. MAIN_PROCESS) return

  print *,'getting phaseshift data...'

  ! open data file
  print *,'  open datafile: '
  print *,'    ',trim(heterogeneousDataFile)
  open(dataUnit,file=trim(heterogeneousDataFile),status="old",iostat=ier)
  if (ier /= 0) then
    print *,'Error: process ',myrank,' could not open data file: ',trim(heterogeneousDataFile)
    call stopProgram( "data file not found     ")
  endif

  !------read its header
  !------header MUST contain period (in s) at which measurements were made,
  !------and reference velocity of that type of sw at that period
  !read(2,1003) wave_type,period,v0
  read(dataUnit,*) wave_type,period,v0
  print *,'    wave type/period/v0: ',wave_type,'/',period,'/',v0
  read(dataUnit,'(a47)') commentline

  omega = 1.0/period

1003 format(8x,a1,7x,f8.4,10x,f8.5)
!1002 format(4(f11.4,1x),1x,f11.6,1x,f11.7)

  ! count data entries
  stla_old = 0.0
  stlo_old = 0.0
  epla_old = 0.0
  eplo_old = 0.0

  allocate(stations(10000,2), events(20000,2),stat=ier)
  if (ier /= 0) call stopProgram('Error allocating stations/events    ')
  stations(:,:) = 0.0
  events(:,:) = 0.0

  epicentraldistance_min = 1000.0 ! modified d.peter (8.2.2006)
  epicentraldistance_max = 0.0
  datummax = 0.0
  datummin = 0.0

  datacount = 0
  stationscount = 0
  eventscount = 0
  do while( ier == 0 )
    !read(2,1002,iostat=ier) epla,eplo,stla,stlo,datum,error
    read(dataUnit,*,iostat=ier) epla,eplo,stla,stlo,datum,error
    if (ier == 0) then
      datacount = datacount + 1

      ! receivers: file is (normally) ordered by receiver location
      stationfound = .false.
      do i = 1,stationscount
        if (abs(stla -stations(i,1)) < 1.e-4 .and. &
           abs(stlo - stations(i,2)) < 1.e-4) then
          stationfound = .true.
          exit
        endif
      enddo
      if (.not. stationfound) then
        stationscount = stationscount+1
        stations(stationscount,1) = stla
        stations(stationscount,2) = stlo
      endif

      ! events/sources: file is (normally) ordered by receiver location,
      !             the stations can have different number of events
      eventfound = .false.
      do i = 1,eventscount
        if (abs(epla-events(i,1)) < 1.e-4 .and. &
           abs(eplo-events(i,2)) < 1.e-4) then
          eventfound = .true.
          exit
        endif
      enddo
      if (.not. eventfound) then
        eventscount = eventscount+1
        events(eventscount,1) = epla
        events(eventscount,2) = eplo
      endif

      ! some statistics
      if (datum > datummax) datummax = datum
      if (datum < datummin) datummin = datum

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
  print *,'    number of stations : ',stationscount
  print *,'    number of events   : ',eventscount
  print *,'    great circle distances:'
  print *,'      minimum = ',epicentraldistance_min,epiminIndex
  print *,'      maximum = ',epicentraldistance_max,epimaxIndex
  print *,'    datum data:'
  print *,'      minimum = ',datummin
  print *,'      maximum = ',datummax

  ! open new output data file
  slashindex = 0
  do i = 1,len_trim(heterogeneousDataFile)
    if (heterogeneousDataFile(i:i) == '/') slashindex = i
  enddo
  nameout = trim(datadirectory)//heterogeneousDataFile(slashindex+1:len_trim(heterogeneousDataFile))//'.new'
  open(newdataUnit,file=trim(nameout),iostat=ier)
  if (ier /= 0) call stopProgram("Error opening new data file .new for output   ")

  print *,'    new datum file will be:'
  print *,'      ',trim(nameout)

  ! write first two line
  write(newdataUnit,1003) wave_type,bw_waveperiod,cphaseRef
  write(newdataUnit,'(a47)') commentline
  call myflush(newdataUnit)

  ! info
  print *,'  data read done'
  print *
  call myflush(6)

  end subroutine


!-----------------------------------------------------------------------
  subroutine storeHeterogeneousDataArray(dataUnit,datacount)
!-----------------------------------------------------------------------
! stores the data in an array
  use parallel; use phaseBlockData; use heterogeneousArrays
  implicit none
  integer,intent(in):: dataUnit
  integer,intent(inout):: datacount
  ! local parameters
  integer:: i,ier,length
  character:: commentline*47,wave_type*1
  real:: period,v0
  real:: eplo,epla,stla,stlo,datum,error

  ! info
  if (MAIN_PROCESS) then
    print *,'preparing data array...'
    print *,'  data count: ',datacount
  endif

  ! wait until all processes reached this point
  call syncProcesses()

  ! synchronize datacount from main process to all other processes
  call syncBroadcastSinglei(datacount)
  !or
  !call syncSingleInteger(datacount)

  ! allocate array
  allocate(dataStore(6,datacount),stat=ier)
  if (ier /= 0) call stopProgram('error allocating data - storeHeterogeneousDataArray()    ')

  ! process data file
  if (MAIN_PROCESS) then
    ! reset to first data entry
    rewind(dataUnit)
    read(dataUnit,*) wave_type,period,v0
    read(dataUnit,'(a47)') commentline

    ! read data
    print *,'  reading in data array...'
    print *,'    number of data entries: ',datacount
    do i = 1,datacount
      read(dataUnit,*,iostat=ier) epla,eplo,stla,stlo,datum,error
      if (ier /= 0) call stopProgram('error reading data - storeHeterogeneousDataArray()    ')

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

  ! synchronize dataStore array with processes
  length = 6*datacount
  call syncBroadcastRealArray(dataStore,length)
  !or
  !call syncRealArray(dataStore,length)

  ! TEST console output
  if (MAIN_PROCESS) then
    print *,'    first datum extract: ',(dataStore(i,1),i=1,6)
    print *
    call myflush(6)
  endif

  ! wait until all processes reached this point
  call syncProcesses()

  end subroutine


!-----------------------------------------------------------------------
  subroutine processDatafile(newdataUnit,datacount)
!-----------------------------------------------------------------------
! handle events
  use parallel; use propagationStartup; use griddomain; use heterogeneousArrays
  use verbosity
  implicit none
  integer,intent(in):: newdataUnit,datacount
  ! local parameters
  integer:: j,dataLoopIndex
  real(WP):: benchmarkLoopStart,benchmarkLoopEnd
  real:: epidelta
  real:: eplo,epla,stla,stlo,datum,error
  real:: time_shift
  ! console output with information about current data index
  integer :: info_interval
  logical:: output_info,do_simulation
  logical, external:: doCalc
  real(WP), external:: syncWtime

  if (MAIN_PROCESS) then
    print *,'processing data...'
  endif

  ! prepare phasevelocities for time iterations
  call constructPhaseVelocities()

  ! loop over data
  if (MAIN_PROCESS) then
    print *,'  begin loop over data: ',datacount
    print *
    call myflush(6)
  endif

  ! determine interval for output infos
  info_interval = int(datacount / 20)
  if (info_interval > 500) info_interval = 500
  if (info_interval < 5) info_interval = 5

  ! stats
  shift_min = 0.0
  shift_max = 0.0

  dataLoopIndex = 0
  output_info = .false.

  do j = 1,datacount
    ! flag if process needs to do simulation for this data index
    do_simulation = doCalc(j)

    ! console output
    if (MAIN_PROCESS) then
      if (PARALLELSEISMO) then
        ! parallelized (forward) simulations
        if (mod(j,info_interval) == 0 .or. j == 1) then
          VERBOSE = .true.      ! turns on verbosity of source initialization
          output_info = .true.  ! console info output & file trace output
        else
          ! turn off verbosity
          VERBOSE = .false.
          output_info = .false.  ! console info output
        endif
      else
        ! distributed simulations (each process simulates different source/data entry)
        if (mod(j,info_interval) < nprocesses .and. do_simulation) then
          VERBOSE = .true.      ! turns on verbosity of source initialization
          output_info = .true.  ! console info output & file trace output
        else
          ! turn off verbosity
          VERBOSE = .false.
          output_info = .false.  ! console info output
        endif
      endif
    else
      VERBOSE = .false.
    endif

    ! info header
    if (output_info) then
      print *,'*****************************************************************'
      print *,'*** processing data index : ',j,' out of ',datacount,'***'
      print *,'*****************************************************************'
      call myflush(6)
    endif

    ! check if this process has to do something
    if (do_simulation) then
      ! take data from array
      epla = dataStore(1,j)
      eplo = dataStore(2,j)
      stla = dataStore(3,j)
      stlo = dataStore(4,j)
      datum = dataStore(5,j)
      error = dataStore(6,j)

      ! benchmark
      if (output_info) benchmarkLoopStart = syncWtime()

      ! prepare source/receiver locations for simulation
      call prepareCouple(epla,eplo,stla,stlo)

      ! do the time iteration and phase shift calculation
      call doTimeIteration(j,time_shift,output_info)

      ! secondary processes transfers data to main process which prints them to file
      call outputDataLine(j,datacount,epla,eplo,stla,stlo,time_shift,error,newdataUnit)

      ! output to console
      if (output_info) then
        ! find epicentral distances between source and station
        call epicentralDistance(epla,eplo,stla,stlo,epidelta)

        ! console output
        print *,'    data: ',j
        print *,'      source                  : ',sourceLat,sourceLon
        print *,'        (original)            : ',epla,eplo
        print *,'      receiver                : ',receiverLat,receiverLon
        print *,'        (original)            : ',stla,stlo
        print *,'      distance source-receiver: ',epidelta
        print *,'      time shift              : ',time_shift

        ! benchmark for this run
        benchmarkLoopEnd = syncWtime()
        print *,'    benchmark seconds ',benchmarkLoopEnd-benchmarkLoopStart
        print *
        ! flush to console output (i/o unit nr. 6 or on SGI system nr. 101)
        call myflush(6)
      endif

      ! for possible exits in a job queue
      if (j > dataLoopIndex) dataLoopIndex = j
    endif
  enddo

  ! wait until all processes reached this point
  call syncProcesses()

  ! close new data file
  if (MAIN_PROCESS) then
    close(newdataUnit)
  endif

  ! console output
  if (MAIN_PROCESS) then
    VERBOSE = .true.
    print *,'data written:',datacount
    print *
    call myflush(6)
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine constructPhaseVelocities()
!-----------------------------------------------------------------------
! determine phase velocity array for PREM and heterogeneous run
  use propagationStartup; use phaseVelocityMap; use parallel; use verbosity
  use heterogeneousArrays; use cells
  implicit none
  ! local parameters
  integer:: ier

  ! allocate memory
  allocate(phaseVelocitySquarePREM(numVertices), &
           phaseVelocitySquareHET(numVertices),stat=ier)
  if (ier /= 0) stop "allocate phasevelocities error"
  phaseVelocitySquarePREM(:) = 0.0_WP
  phaseVelocitySquareHET(:) = 0.0_WP

  ! PREM
  if (MAIN_PROCESS .and. VERBOSE) then
    print *
    print *,'  constructing reference phase map...'
  endif
  HETEROGENEOUS = .false. ! homogeneous
  DELTA         = .false. ! no delta scatterer
  call constructPhaseVelocitySquare()
  phaseVelocitySquarePREM(:) = phaseVelocitySquare(:)

  ! heterogeneous phase speed map
  if (MAIN_PROCESS .and. VERBOSE) then
    print *
    print *,'  constructing heterogeneous phase map...'
  endif
  HETEROGENEOUS = .true.
  DELTA         = .false.  ! no delta scatterer

  ! check that no rotation of phase map
  if (ROTATE_FRAME) call stopProgram('heterogeneouse phase map needs rotate_frame set to .false.    ')

  !allocate phaseMap array
  if (allocated(phaseMap)) deallocate(phaseMap)
  allocate(phaseMap(numVertices), stat=ier)
  if (ier /= 0) stop 'Abort - phase map allocation'
  phaseMap(:) = 0.0_WP

  ! read phase map (e.g. Love 150 s) (note: is rotated for non-adjoint simulations)
  if (MAIN_PROCESS) then
    call constructPhasedata()
  endif

  ! synchronize between main and secondary processes
  call syncPhaseMap()

  ! get phase velocities
  call constructPhaseVelocitySquare()
  phaseVelocitySquareHET(:) = phaseVelocitySquare(:)

  ! phase map no more needed
  if (HETEROGENEOUS) deallocate(phaseMap)

  if (MAIN_PROCESS .and. VERBOSE) print *

  end subroutine


!-----------------------------------------------------------------------
  logical function doCalc(index)
!-----------------------------------------------------------------------
! determines wheter the data entry (index) should be proceeded
! and calculations should be done
  use parallel
  implicit none
  integer,intent(in):: index

  ! default
  doCalc = .true.

  ! parallel processing
  if (nprocesses < 2) return
  if (PARALLELSEISMO) return

  ! distribute following calculation over all processes
  ! example: > mpirun -np 8 ./bin/heterogeneousPhaseshift
  !            main process: j=1,9,17,...
  !            rank 1      : j=2,10,...
  !            rank 2      : j=3,11,...
  if (mod((index-myrank+2*nprocesses-1),nprocesses) == 0) then
    doCalc = .true.
  else
    doCalc = .false.
  endif

  return
  end function


!-----------------------------------------------------------------------
  subroutine doTimeIteration(jrecord,time_shift,output_info)
!-----------------------------------------------------------------------
  use propagationStartup; use phaseVelocityMap
  use traveltime; use cells; use verbosity; use parallel
  use heterogeneousArrays
  implicit none
  integer,intent(in):: jrecord
  real,intent(out):: time_shift
  logical,intent(in):: output_info
  ! local parameters
  real(WP),dimension(2,numofTimeSteps):: seismogramPREM,seismogramHET
  real(WP):: startTime,distance,amplification
  integer:: ier,i
  character(len=6):: jstr

  ! be sure only one receiver at a time
  manyReceivers = .false.

  ! PREM reference simulation
  phaseVelocitySquare(:) = phaseVelocitySquarePREM(:)

  call forwardIteration()

  seismogramPREM(:,:) = seismogramReceiver(:,:)

  ! file output
  if (output_info) then
    write(jstr,'(i6.6)') jrecord
    open(IOUT,file=trim(datadirectory)//'seismo.PREM.'//jstr//'.txt',iostat=ier)
    if (ier /= 0) call stopProgram('Error opening seismo.PREM output file    ')
    do i = 1,numofTimeSteps
      write(IOUT,*) seismogramPREM(1,i), seismogramPREM(2,i)
    enddo
    close(IOUT)
  endif

  ! heterogeneous simulation
  phaseVelocitySquare(:) = phaseVelocitySquareHET(:)

  call forwardIteration()

  seismogramHET(:,:) = seismogramReceiver(:,:)

  ! file output
  if (output_info) then
    open(IOUT,file=trim(datadirectory)//'seismo.HET.'//jstr//'.txt',iostat=ier)
    if (ier /= 0) call stopProgram('Error opening seismo.HET output file    ')
    do i = 1,numofTimeSteps
      write(IOUT,*)seismogramHET(1,i), seismogramHET(2,i)
    enddo
    close(IOUT)
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
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'  determining phase shift...'
    print *,'    arrival time                         : ',arrivalTime
    print *,'    start time for fourier transformation: ',startTime
  endif

  ! calculate time lag (in seconds, when dt in seconds)
  call getTimelagSeismos(t_lag,seismogramHET,seismogramPREM,numofTimeSteps,dt)

  ! find Euler angles and rotation matrix (code by John W.)
  !call pole(epla,eplo,stla,stlo,xlatp,xlonp,azmp,delta)
  !call epicentralDistance(epla,eplo,stla,stlo,epidelta)

  ! t0: traveltime (s) with factor PI/180 * 6371 km ~ 111.2
  !t0 = epidelta*111.19492664/v0
  ! relative phase shift = relative traveltime anomaly
  !dphi_phi = datum/t0
  !print *,'    relative traveltime = ',dphi_phi,datum,t0
  !dphi_phi = datum            !to invert for absolute delta_p

  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'    cross-correlation lag : ',t_lag
  endif

  ! downhill simplex
  call getTimelagSimplex(t_lag,amplification,seismogramHET,seismogramPREM,numofTimeSteps)

  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'    downhill simplex  lag : ',t_lag
  endif

  ! set phase shift
  time_shift = t_lag

  ! console output
  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'    using phase time shift: ',time_shift
    print *
    call myflush(6)
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine outputDataLine(jrecord,datacount,epla,eplo,stla,stlo,shift,error,newdataUnit)
!-----------------------------------------------------------------------
! sync with main process and print to file
  use parallel
  implicit none
  integer,intent(in):: jrecord,datacount,newdataUnit
  real,intent(in):: shift,epla,eplo,stla,stlo,error
  ! local parameters
  integer:: restprocesses,i
  real:: epla_proc,eplo_proc,stla_proc,stlo_proc,shift_proc,error_proc

  ! single processor/single simulation job
  if (nprocesses < 2 .or. PARALLELSEISMO) then
    ! only main process writes out data
    if (MAIN_PROCESS) then
      call printDataLine(newdataUnit,epla,eplo,stla,stlo,shift,error)
    endif
    ! all done
    return
  endif

  ! more processes: get phaseshift data and print to file
  if (MAIN_PROCESS) then
    ! print main info first
    call printDataLine(newdataUnit,epla,eplo,stla,stlo,shift,error)

    ! how many process are there to collect
    restprocesses = datacount - jrecord
    if (restprocesses >= nprocesses) then
      restprocesses = nprocesses - 1
    endif

    ! check if something to receive
    if (restprocesses == 0) return

    ! main receives
    do i = 1,restprocesses
      !print *,'getting from',i
      ! get data
      call syncRecvDataFromProc(i,shift_proc,epla_proc,eplo_proc,stla_proc,stlo_proc,error_proc)

      ! print to file
      call printDataLine(newdataUnit,epla_proc,eplo_proc,stla_proc,stlo_proc,shift_proc,error_proc)
    enddo
  else
    ! secondary process sends information
    call syncSendDataFromProc(myrank,shift,epla,eplo,stla,stlo,error)
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine printDataLine(newdataUnit,epla,eplo,stla,stlo,shift,error)
!-----------------------------------------------------------------------
! writes given data to file
  use heterogeneousArrays
  implicit none
  integer,intent(in):: newdataUnit
  real,intent(in):: epla,eplo,stla,stlo,shift,error

  ! new data file output
  write(newdataUnit,1002) epla,eplo,stla,stlo,shift,error

 1002 format(4(f11.4,1x),1x,f11.6,1x,f11.7)

  call myflush(newdataUnit)

  ! statistics
  if (shift < shift_min) shift_min = shift
  if (shift > shift_max) shift_max = shift

  end subroutine


!-----------------------------------------------------------------------
  subroutine myflush(i)
!-----------------------------------------------------------------------
! should write out buffered data
  implicit none
  integer, intent(in) :: i

  ! flush to console output (i/o unit nr. 6 or on SGI system nr. 101)
  ! if (i == 6) i = 101

  !  #ifdef _IBM_
    !call flush_(i)
  !#else
    FLUSH(i)
  !#endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine finish()
!-----------------------------------------------------------------------
! finalize MPI communication
  use parallel; use propagationStartup; use heterogeneousArrays
  implicit none
  ! local parameters
  real(WP), external:: syncWtime

  ! benchmark
  if (MAIN_PROCESS) then
    benchAllEnd = syncWtime()
    print *,'statistics:'
    print *,'  datum min/max:',shift_min,shift_max
    print *
    print *,'number of processors:',nprocesses
    print *,'running time: ',int((benchAllEnd-benchAllStart)/60.0),'min ', &
                              mod((benchAllEnd-benchAllStart),60.0),'sec'
    print *
  endif

  ! end parallelization
  call syncFinalizeMPI()

  end subroutine
