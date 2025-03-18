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
      subroutine runSimulations()
!-----------------------------------------------------------------------
! runs all numerical simulations for membrane waves
!
! finite-difference iteration: information in carl tape thesis2003, chap 5, (5.7)
      use propagationStartup
      implicit none
      ! local parameters
      logical:: looping

      ! loop for each delta location
      looping = .true.
      referenceRun = .true.
      do while( looping )
        ! prepare for simulation
        call prepareSimulation()

        ! do the time iteration
        call forwardIteration()

        ! save seismograms, determine phaseshifts and next delta location
        call processSolutions(looping)

      enddo !delta location

      end subroutine


!-----------------------------------------------------------------------
      subroutine prepareSimulation()
!-----------------------------------------------------------------------
! iterates through time and calculates the solutions at each time step
!
! finite-difference iteration: information in carl tape thesis2003, chap 5, (5.7)
      use propagationStartup; use parallel; use displacements; use verbosity
      implicit none
      ! local parameters
      integer:: ier

      ! make first a reference run without scatterer, then all others with scatterer
      ! assume that all processes have access to the same data directory
      if (PARALLELSEISMO) then
        if (referenceRun) then
          ! simulation without scatterer
          DELTA = .false.
          if (MAIN_PROCESS) then
            if (VERBOSE) print *,'running reference simulation without scatterer...'
          endif
        endif
      else
        if (referenceRun) then
          if (MAIN_PROCESS) then
            ! simulation without scatterer
            DELTA = .false.
            if (VERBOSE) print *,'running reference simulation without scatterer...'
          else
            ! only doing simulation with scatterer
            referenceRun = .false.
            DELTA = .true.

            ! wait for reference simulation
            call MPI_Barrier( MPI_COMM_WORLD, ier )
            if (ier /= 0) call stopProgram( "abort - program phaseshift MPI_Barrier failed    ")

            ! synchronize reference seismograms
            call syncReceiversRef()
          endif
        endif !referenceRun
      endif !parallelseismo

      ! benchmark
      if (MAIN_PROCESS) benchstart = MPI_WTIME()

      ! reset displacements
      displacement_old(:)=0.0_WP
      displacement(:)=0.0_WP
      newdisplacement(:)=0.0_WP

      ! precalculate the phase velocities for all grid points
      call constructPhaseVelocitySquare()

      end subroutine


!-----------------------------------------------------------------------
      subroutine processSolutions(looping)
!-----------------------------------------------------------------------
! determine phaseshifts and next location of scatterer
!
! finite-difference iteration: information in carl tape thesis2003, chap 5, (5.7)
      use propagationStartup; use cells; use phaseVelocityMap
      use parallel; use displacements
      use loop; use phaseBlockData; use griddomain; use verbosity
      implicit none
      logical,intent(inout):: looping
      ! local parameters
      logical:: available
      integer:: m,ier
      real(WP):: timerStart,timerEnd
      character(len=128):: datafile

      ! print seismogram to file
      if (.not. manyReceivers) call printSeismogram()

      ! benchmark output
      if (MAIN_PROCESS .and. VERBOSE) then
        benchend = MPI_WTIME()
        print *
        print *,'benchmark seconds:',benchend-benchstart
        print *
      endif


      if (.not. referenceRun) then
        ! see if the reference seismogram is available
        available = .false.
        if (.not. manyReceivers) then
          ! see if file is there
          datafile = trim(datadirectory)//'seismo.'//cphasetype(1:4)//'.withoutDelta.dat'
          inquire(file=trim(datafile),exist=available)
          if (.not. available) then
            ! wait and poll until it exists
            timerStart = MPI_WTIME()
            do while(.not. available )
              inquire(file=trim(datafile),exist=available)
              timerEnd = MPI_WTIME()
              if (abs(timerEnd - timerStart) > 120.0) then
                print *,'Error: did not find reference seismogram:',trim(datafile)
                print *,'       time out of process',rank
                print *,'shutting down...'
                call stopProgram( 'abort - processSolutions()   ')
              endif
            enddo
          endif
        endif
        !else
          !datafile = trim(datadirectory)//'seismo.'//cphasetype(1:4)//'.withoutDelta.360.dat'
        !endif

        ! both seismograms are available, get phaseshift
        if (manyKernels) then
          do m = 1,numofKernels
            currentKernel = m
            call calculatePhaseshifts()
          enddo
        else
          call calculatePhaseshifts()
        endif

        !move delta location
        if (DELTA) then
          if (PARALLELSEISMO) then
            call parallelFindnextLocation(looping)
          else
            call findnextLocation(looping)
          endif
        endif

        !check looping
        if (.not. MOVEDELTA  .or. .not. DELTA) looping = .false.

      else
        if (PARALLELSEISMO) then
          if (manyReceivers) then
            ! update reference seismograms of main process
            !print *,'collecting reference seismograms...',rank,referenceRun
            call collectReceiversSeismogram(referenceRun)
          endif

          !simulate reference run only once
          referenceRun = .false.
          DELTA = .true.
        else
          if (MAIN_PROCESS) then
            !simulate reference run only once
            referenceRun = .false.
            DELTA = .true.

            ! let all processes pass to simulate scattered seismograms
            call MPI_Barrier( MPI_COMM_WORLD, ier )
            if (ier /= 0) call stopProgram( "abort - processSolutions() MPI_Barrier failed    ")

            ! synchronize reference seismograms
            call syncReceiversRef()
          else
            print *,'Error: unknown condition, going to terminate this fuzzy program...'
            call stopProgram( 'abort - processSolutions()   ')
          endif
        endif !parallelseismo
      endif

      end subroutine


!-----------------------------------------------------------------------
      subroutine calculatePhaseshifts()
!-----------------------------------------------------------------------
! determines the phase shift between reference and scattered seismograms
!
! returns: new phaseshift array entry
      use cells; use loop; use phaseVelocityMap; use traveltime
      use filterType; use verbosity
      use parallel;use propagationStartup
      implicit none
      ! local parameters
      integer:: i,recvVertex,n,ier  !,filenameBaseLength,index
      character(len=128):: filename,fileReference
      real(WP):: vperturbation,phaseVelocityRef,realdeltaLat,realdeltaLon, &
                 correctedLat,correctedLon
      real(WP):: recvLat,recvLon,vsource(3),vreceiver(3),vdelta(3),mat(3,3),startTime
      real(WP),external:: getKernelValue

      ! in parallel simulation, only main process needs to calculate the phaseshift
      if (PARALLELSEISMO) then
        if (manyReceivers) then
          ! update the receivers seismograms of the main process (receivers
          ! may lie in the other domains as well)
          !print *,'collecting receivers seismograms...',rank,referenceRun
          call collectReceiversSeismogram(referenceRun)
        endif

        ! slaves can quit
        if (.not. MAIN_PROCESS) return
      endif

      ! get parameters of program initialization
      vperturbation = deltaPerturbation
      phaseVelocityRef = cphaseRef
      vsource(:) =  vertices(sourceVertex,:)
      if (manyKernels) then
        vreceiver(:) = vertices(kernelsReceivers(currentKernel,1),:)
      else
        vreceiver(:) = vertices(receiverVertex,:)
      endif

      ! file handling
      call prepareForKernel(filename,fileReference)

      ! get start time of reference seismogram
      call getStartingTime(fileReference,startTime)

      ! get rotation matrix to rotate into equatorial plane
      if (rotate_frame) call getInverseRotationMatrix(vsource,vreceiver,mat)

      ! console output only from main process
      if (.not. MAIN_PROCESS) beVerbose = .false.
      if (beVerbose) then
        print *
        print *,'phaseshift: '
        print *,'    velocity perturbation correction:',vperturbation
      endif

      ! calculate kernel-values
      if (.not. manyReceivers) then
        ! get time lag
        call getTimelag(t_lag,filename,fileReference,startTime)
        ! kernel value term
        kernel = getKernelValue(deltaVertex,heterogeneous,t_lag, &
                                phaseVelocityRef,arrivalTime,vperturbation)

        if (ANALYTICAL_CORRELATION) then
          ! bw_width must be set by earlier call to getTimelag()
          ! analytic time lag
          call getAnalyticalTimelag(t_lagAnalytic,filename,fileReference, &
                                    startTime,bw_width,BUTTERWORTH_POWER,bw_waveperiod,FILTERSEISMOGRAMS)
          ! kernel value term
          kernelAnalytic = getKernelValue(deltaVertex,heterogeneous,t_lagAnalytic, &
                                          phaseVelocityRef,arrivalTime,vperturbation)
        endif

        ! get actual vertex position
        call getSphericalCoord_Lat(deltaVertex,realdeltaLat,realdeltaLon)
        if (realdeltaLon > 180.0_WP) realdeltaLon= -(360.0_WP-realdeltaLon)

        call getSphericalCoord_Lat(receiverVertex,recvLat,recvLon)

        ! output
        write(80,'(2f8.2, 2e16.6, i10,5f12.6)') real(realdeltaLon), &
                  real(realdeltaLat),kernel,kernelAnalytic, &
                  receiverVertex, t_lag,t_lagAnalytic,vperturbation,recvLat,recvLon
        ! rotate into equatorial plane
        if (rotate_frame) then
          call getVector( realdeltaLat,realdeltaLon,vdelta(1),vdelta(2),vdelta(3) )
          call rotateVector(mat,vdelta,vdelta)
          call getSphericalCoordinates(vdelta,realdeltaLat,realdeltaLon)

          ! print to file
          !print '(2f8.2, f12.6, i10)',real(realdeltaLon),real(realdeltaLat),kernel, deltaVertex
          write(81,'(2f8.2, 2e16.6, i10,5f12.6)') real(realdeltaLon),real(realdeltaLat), &
                      kernel,kernelAnalytic, &
                      receiverVertex, t_lag,t_lagAnalytic,vperturbation,recvLat,recvLon
        endif
      else
        ! get actual perturbation vertex position
        call getSphericalCoord_Lat(deltaVertex,realdeltaLat,realdeltaLon)
        if (realdeltaLon > 180.0_WP) realdeltaLon= -(360.0_WP-realdeltaLon)
        !if (VERBOSE .and. MAIN_PROCESS) print *,'delta vertex:',deltaVertex,realdeltaLat,realdeltaLon

        ! get phase shifts for each receiver location
        ! we have receivers located all around the equator
        do i = 1,numofReceivers
          ! to have a order of increasing longitudes
          if (i == 1) then
            n = 1
          else
            n = numofReceivers-i+2
          endif

          ! determine filenames
          !write(recstr,'(i3.2)') int(n)
          !filename=filename(1:filenameBaseLength)//recstr//'.dat'
          !fileReference=fileReference(1:len_trim(fileReference)-7)//recstr//'.dat'

          ! get time lag
          call getTimelagRef(t_lag,n,startTime)

          ! check if NaN
          !if (isNaN(t_lag)) then
          !  print *,'Error: time lag has NaN value!',t_lag
          !  print *,'       receiver:',n
          !  print *,'       startTime:',startTime
          !  call stopProgram()
          !endif

          ! kernel value term
          kernel = getKernelValue(deltaVertex,heterogeneous,t_lag,phaseVelocityRef, &
                                arrivalTime,vperturbation)

          !if (ANALYTICAL_CORRELATION) then
          !  ! analytic time lag
          !  call getAnalyticalTimelag(t_lagAnalytic,filename,fileReference,startTime, &
          !                            bw_width,BUTTERWORTH_POWER,bw_waveperiod,FILTERSEISMOGRAMS)
          !  ! kernel value term
          !  kernelAnalytic = getKernelValue(deltaVertex,heterogeneous,t_lagAnalytic, &
          !                                  phaseVelocityRef,arrivalTime,vperturbation)
          !endif

          ! get receiver position
          if (manyKernels) then
            call getSphericalCoord_Lat(kernelsReceivers(currentKernel,n),recvLat,recvLon)
          else
            call getSphericalCoord_Lat(receivers(n),recvLat,recvLon)
          endif

          ! relativ position of delta towards receiver
          correctedLon = realdeltaLon - recvLon + 360.0_WP
          if (correctedLon >= 360.0_WP) correctedLon=correctedLon-360.0_WP
          correctedLat = realdeltaLat

          ! output
          if (manyKernels) then
            recvVertex=kernelsReceivers(currentKernel,n)
          else
            recvVertex=receivers(n)
          endif

          !print '(2f8.2,2f12.6,i10,5f12.6)',real(correctedLon),real(correctedLat), &
          !  kernel,kernelAnalytic,recvVertex,t_lag,t_lagAnalytic,vperturbation,recvLat,recvLon
          write(80,'(2f8.2, 2e16.6, i10,5f12.6)',iostat=ier) &
              real(correctedLon),real(correctedLat),kernel,kernelAnalytic, &
              recvVertex,t_lag,t_lagAnalytic,vperturbation,recvLat,recvLon
          if (ier /= 0) then
            print *,'Error: writing to file:',ier
            call stopProgram( 'abort - calculatePhaseshifts   ')
          endif

          ! rotate into equatorial plane
          if (rotate_frame) then
            call getVector( correctedLat,correctedLon,vdelta(1),vdelta(2),vdelta(3) )
            call rotateVector(mat,vdelta,vdelta)
            call getSphericalCoordinates(vdelta,correctedLat,correctedLon)
            ! print to file
            !print '(2f8.2, f12.6, i10)',real(realdeltaLon),real(realdeltaLat),kernel, deltaVertex
            write(81,'(2f8.2, 2e16.6, i10,5f12.6)') &
              real(correctedLon),real(correctedLat),kernel,kernelAnalytic, &
              recvVertex,t_lag,t_lagAnalytic,vperturbation,recvLat,recvLon
          endif
        enddo
      endif
      close(80)
      close(81)

      ! filter around 75 s and get kernels
      !call filterShort()

      end subroutine

!-----------------------------------------------------------------------
      subroutine prepareForKernel(filename,fileReference)
!-----------------------------------------------------------------------
! determines the filenames of reference and actual file and opens the ttkernel*.dat files
!
! returns: new phaseshift array entry
      use cells; use loop; use phaseVelocityMap; use traveltime
      use filterType; use verbosity
      use parallel;use propagationStartup
      implicit none
      character(len=128),intent(inout):: filename,fileReference
      ! local parameters
      character(len=3):: rankstr,kernelstr ! recstr
      character(len=6):: latstr,lonstr
      character(len=128):: ttkernelfile,ttkernelRotatedfile
      integer:: ier

      ! get filename for latitude/longitude
      if (.not. manyReceivers) then
        write(latstr,'(f6.1)') deltaLat
        write(lonstr,'(f6.1)') deltaLon
        latstr = trim(latstr)
        lonstr = trim(lonstr)
        cphasetype = trim(cphasetype)
        filename = trim(datadirectory)//'seismo.'//cphasetype(1:4)//'.'//latstr//'.'//lonstr//'.dat'
        filename = trim(filename)
        fileReference = trim(datadirectory)//'seismo.'//cphasetype(1:4)//'.withoutDelta.dat'
        fileReference = trim(fileReference)
      endif

      !if (manyReceivers) then
      !  fileReference=datadirectory(1:length)//'seismo.'//cphasetype(1:4)//'.withoutDelta. 01.dat'
      !  fileReference=trim(fileReference)
      !  filenameBaseLength=len_trim(filename)-3
      !endif

      ! open output - files for kernel values
      if (manyKernels) then
        write(rankstr,'(i3.3)') rank
        write(kernelstr,'(i3.3)') int(kernelStartDistance+currentKernel-1)

        ttkernelfile = trim(datadirectory)//'ttkernel'//kernelstr//'.rank'//rankstr//'.dat'
        open(80,file=trim(ttkernelfile),position='APPEND',iostat=ier)
        if (ier /= 0) then
          print *,'Error: could not open file ',trim(ttkernelfile)
          call stopProgram( 'abort - prepareForKernel     ')
        endif

        ttkernelRotatedfile = trim(datadirectory)//'ttkernel'//kernelstr//'.rot.rank'//rankstr//'.dat'
        open(81,file=trim(ttkernelRotatedfile),position='APPEND',iostat=ier)
        if (ier /= 0) then
          print *,'Error: could not open file ',trim(ttkernelRotatedfile)
          call stopProgram( 'abort - prepareForKernel     ')
        endif
      else
        write(rankstr,'(i3.3)') rank

        ttkernelfile = trim(datadirectory)//'ttkernel.rank'//rankstr//'.dat'
        open(80,file=trim(ttkernelfile),position='APPEND',iostat=ier)
        if (ier /= 0) then
          print *,'Error: could not open file',trim(ttkernelfile)
          call stopProgram( 'abort - prepareForKernel     ')
        endif

        ttkernelRotatedfile = trim(datadirectory)//'ttkernel.rot.rank'//rankstr//'.dat'
        open(81,file=trim(ttkernelRotatedfile),position='APPEND',iostat=ier)
        if (ier /= 0) then
          print *,'Error: could not open file ',trim(ttkernelRotatedfile)
          call stopProgram( 'abort - prepareForKernel     ')
        endif
      endif

      end subroutine

!-----------------------------------------------------------------------
      subroutine getStartingTime(fileReference,startTime)
!-----------------------------------------------------------------------
! determines the startTime and arrivalTime
!
! returns: startTime
      use cells; use loop; use phaseVelocityMap; use traveltime
      use filterType; use verbosity
      use parallel;use propagationStartup
      implicit none
      character(len=128),intent(in):: fileReference
      real(WP),intent(out):: startTime
      ! local parameters
      real(WP):: distance,seismo(2,lasttimestep-firsttimestep+1)

      ! initializes
      startTime = 0.0_WP

      ! start time
      if (manyReceivers) then
        if (manyKernels) then
          ! time
          seismo(1,:) = kernelsReceiversSeismogramRef(currentKernel,numofReceivers+1,:)
          ! displacement of first receiver
          seismo(2,:) = kernelsReceiversSeismogramRef(currentKernel,1,:)
          call getStartTimeSeismogram(seismo,lasttimestep-firsttimestep+1, &
                                      seismo(1,lasttimestep-firsttimestep+1),dt,startTime)
          !kernelsStartTime(currentKernel)=startTime
        else
          ! time
          seismo(1,:) = receiversSeismogramRef(numofReceivers+1,:)
          ! displacement of first receiver
          seismo(2,:) = receiversSeismogramRef(1,:)
          call getStartTimeSeismogram(seismo,lasttimestep-firsttimestep+1, &
                                      seismo(1,lasttimestep-firsttimestep+1),dt,startTime)
        endif
      else
        call getStartTime(fileReference,startTime)
      endif

      ! arrival time
      if (manyKernels) then
          ! calculate the analytical arrival time
          call greatCircleDistance(vertices(sourceVertex,:), &
                                   vertices(kernelsReceivers(currentKernel,1),:),distance)
          arrivalTime = distance*EARTHRADIUS/cphaseRef
          ! output
          if (MAIN_PROCESS .and. beVerbose) then
            call greatCircleDistance(vertices(sourceVertex,:), &
                                     vertices(kernelsReceivers(currentKernel,1),:),distance)
            print *,'kernel distance: ',distance*180.0/PI
            print *,'    arrival time:',arrivalTime
            print *,'    start time for fourier transformation:',startTime
            print *
          endif
      else
        ! calculate the analytical arrival time
        call greatCircleDistance(vertices(sourceVertex,:),vertices(receiverVertex,:), &
                                distance)
        arrivalTime = distance*EARTHRADIUS/cphaseRef
        ! output
        if (MAIN_PROCESS .and. beVerbose) then
          print *,'    arrival time:',arrivalTime
          print *,'    start time for fourier transformation:',startTime
          print *
        endif
      endif

      end subroutine


!-----------------------------------------------------------------------
      subroutine collectData()
!-----------------------------------------------------------------------
! collects kernel data from all process files
!
! returns: outputs to files 'ttkernel.dat' and 'ttkernel.rot.dat'
      use parallel; use propagationStartup;use traveltime;use cells;use verbosity
      implicit none
      ! local parameters
      integer:: index,n,i,ier,recvVertex
      character(len=128):: line,datafile,kernelfile  ! fileReference
      character(len=3):: rankstr,kernelstr  ! recstr
      real(WP):: lat,lon,t_lagfromFile,t_laganalyticfromFile,vperturbation,recvlon,recvlat
      real(WP):: kernelfromFile,kernelanalyticfromFile,distance

      ! only main process collects
      if (.not. MAIN_PROCESS) return

      ! collect data out of all rank-files
      if (VERBOSE) then
        print *
        print *,'collecting data...'
      endif

      ! phase shift data from all processes
      ! all data from each process should be available in same directory
      do index = 1,2
        ! open combined, total data file
        if (index == 1) then
          if (manyKernels) then
            write(kernelstr,'(i3.3)') int(kernelStartDistance+currentKernel-1)
            kernelfile = trim(datadirectory)//trim(cphasetype)//'_'//kernelstr//'.dat'
          else
            kernelfile = trim(datadirectory)//'ttkernel.dat'
          endif
          open(60,file=trim(kernelfile),iostat=ier)
          if (ier /= 0) call stopProgram('abort - collectData could not open ttkernel-file    ')
        else
          if (manyKernels) then
            write(kernelstr,'(i3.3)') int(kernelStartDistance+currentKernel-1)
            kernelfile = trim(datadirectory)//trim(cphasetype)//'_'//kernelstr//'.rot.dat'
          else
            kernelfile = trim(datadirectory)//'ttkernel.rot.dat'
          endif
          open(60,file=trim(kernelfile),iostat=ier)
          if (ier /= 0) call stopProgram('abort - collectData could not open ttkernel-file    ')
        endif

        if (manyKernels) then
          call greatCircleDistance(vertices(sourceVertex,:), &
                      vertices(kernelsReceivers(currentKernel,1),:),distance)
          distance = distance*180.0/PI
        else
          call greatCircleDistance(vertices(sourceVertex,:), &
                      vertices(receiverVertex,:),distance)
          distance = distance*180.0/PI
        endif

        if (index == 1) then
          write(60,*) '# total phase shift data ',distance
        else
          write(60,*) '# total phase shift data - rotated kernel',distance
        endif
        write(60,*) '# lon lat phaseshiftKernel relativePhaseshift absolutePhaseshift'

        n = 0
        do i = 0,nprocesses-1
          !open data file
          write(rankstr,'(i3.3)') i
          if (index == 1) then
            if (manyKernels) then
              datafile = trim(datadirectory)//'ttkernel'//kernelstr//'.rank'//rankstr//'.dat'
            else
              datafile = trim(datadirectory)//'ttkernel.rank'//rankstr//'.dat'
            endif
            open(70,file=trim(datafile),iostat=ier)
            if (ier /= 0) then
              print *,'  could not read data from ',trim(datafile)
              print *,'  continue with others...'
              continue
            endif
          else
            if (manyKernels) then
              datafile = trim(datadirectory)//'ttkernel'//kernelstr//'.rot.rank'//rankstr//'.dat'
            else
              datafile = trim(datadirectory)//'ttkernel.rot.rank'//rankstr//'.dat'
            endif
            open(70,file=trim(datafile),iostat=ier)
            if (ier /= 0) then
              print *,'  could not read data from ',trim(datafile)
              print *,'  continue with others...'
              continue
            endif
          endif

          ! skip first 2 lines
          read(70,*) line
          read(70,*) line

          ! collect data
          do while (ier == 0)
            read(70,'(2f8.2, 2e16.6, i10,5f12.6)',iostat=ier) &
              lon,lat,kernelfromFile,kernelAnalyticFromFile, &
              recvVertex, t_lagfromFile,t_lagAnalyticfromFile, &
              vperturbation,recvLat,recvLon
            if (ier == 0) then
              n = n+1
              !phaseshifts(n,1)=lon
              !phaseshifts(n,2)=lat
              !phaseshifts(n,3)=kernel
              write(80,'(2f8.2,e16.6,2e12.4)') lon,lat,kernelfromFile, &
                t_lagfromFile/arrivalTime,t_lagfromFile
            endif
          enddo
          close(70)
        enddo !nprocesses

        ! output the data
        !do i=1,n
        !  write(100,*) phaseshifts(i,1),phaseshifts(i,2),phaseshifts(i,3)
        !enddo
        close(60)
      enddo

      print *,'  wrote to ',kernelfile(1:(len_trim(kernelfile)-8))//'***.dat'
      print *,'  entries read: ',n
      if (n == 65160) print *,'  kernel full...'
      print *,'ended successfully!'

      ! clean up reference files
      ! for all processes, they should be available in same directory
      !if (manyReceivers) then
      !  fileReference=datadirectory(1:length)//'seismo.'//cphasetype(1:4)//'.withoutDelta. 01.dat'
      !  do i=1,size(receivers)
      !    ! determine filenames
      !    write(recstr,'(i3.2)') int(i)
      !    fileReference=fileReference(1:len_trim(fileReference)-7)//recstr//'.dat'
      !    open(300,file=fileReference,iostat=ier)
      !    if (ier == 0) close(300,status='DELETE')
      !  enddo
      !endif

      end subroutine
