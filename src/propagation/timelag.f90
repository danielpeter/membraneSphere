!-----------------------------------------------------------------------
      program timelag
!-----------------------------------------------------------------------
! calculates the time lag between two seismograms
!
! it reads in two (different) seismograms as input is given,
! then either filters or not the two and cross-correlates them.
! as result, the timelag (in sec) is printed to screen.
      use verbosity;
      use filterType, only: WindowSIZE
      use traveltime, only: t_lag
      use propagationStartup, only: lasttime,numofTimeSteps
      use parallel; use precisions
      implicit none
      ! local parameters
      real(WP):: startingTime,endingTime,amplification,analytict_lag
      character(len=128):: fileDelta,fileReference
      integer:: entries
      logical:: output_files = .true.  ! output debug files
      !-----------------------------------------------------------------------
      ! parameters
      ! most parameters concerning timelag calculation are set in file Timelag_Input
      ! (& default values in commonModules.f90)
      !-----------------------------------------------------------------------

      ! get input parameters
      ! console output
      print *,'Timelag'
      print *,'-----------------------------------------------------'
      call readInputParameters(fileDelta,fileReference,startingTime,endingTime)

      ! for debuging
      fileOutput = output_files
      beVerbose = output_files

      ! read in from startingTime
      if (VERBOSE) print *
      if (VERBOSE) print *,'determine file length...'
      call determineFileLengthTaped(fileReference,startingTime,endingTime,entries,lasttime)
      numofTimeSteps = entries

      ! fft window
      call determineFFTWindowsize(numofTimeSteps,WindowSIZE)

      ! determine filter bandwidth parameters (waveperiod and frequency)
      call determineFilterParameters()

      ! get time lag
      print *
      if (VERBOSE) print *,'calculating time lag...'
      call getTimelagTaped(t_lag,fileDelta,fileReference,startingTime,endingTime)


      ! time lag is the subsample position times the time step size
      ! if seismo lags seismoRef, i.e., is shifted to the right of it, then ans will
      ! show a peak at positive lags
      print *,'  numerical timelag          :',t_lag
      !print *,'  numerical phase anomaly (s):',t_lag/(2*PI/bw_waveperiod)
      ! compare with analytical formula
      if (ANALYTICAL_CORRELATION) then
        fileOutput = output_files
        print *
        if (VERBOSE) print *,'calculating analytical time lag...'
        call getAnalyticalTimelagTaped(analytict_lag,fileDelta,fileReference,startingTime,endingTime)
        print *,'  analytic timelag          :',analytict_lag
        !print *,'  analytic phase anomaly (s):',t_lag/(2*PI/bw_waveperiod)
      endif

      ! get time lag
      print *
      if (VERBOSE) print *,'downhill simplex time lag...'
      ! beVerbose = .true.
      fileOutput = output_files
      call getMinimized(t_lag,amplification,fileDelta,fileReference,startingTime,endingTime)

      print *,'  nonlinear timelag       :',t_lag
      print *,'  nonlinear amplification :',amplification
      print *

      end program


!-----------------------------------------------------------------------
      subroutine readInputParameters(fileDelta,fileReference,startingTime,endingTime)
!-----------------------------------------------------------------------
! read in parameters from file Timelag_Input
!
! input:
!     fileDelta,fileReference - file names of seismograms
!     startingTime                     - time to start fourier transformation window from
!
! returns: fileDelta,fileReference,startingTime
      use verbosity; use propagationStartup
      implicit none
      character(len=128),intent(out)::fileDelta,fileReference
      real(WP),intent(out):: startingTime
      real(WP),intent(out):: endingTime
      ! local parameters
      character(len=128):: inputName,tmp
      character(len=128):: line
      integer:: i,ier,length

      ! open input parameter file
      i = 0
      inputName = 'Timelag_Input'
      open(10,file=trim(inputName),status='old',iostat=ier)
      if (ier /= 0) then
        print *,'Error: opening file ',trim(inputName)
        stop 'abort - opening input'
      endif

      ! parse file for parameters
      do while( ier == 0)
        read(10,'(A128)',iostat=ier) line
        if (ier /= 0) exit

        length = len_trim(line)
        if (length == 0) then
          continue
        else
          if (line(1:1) == "%" .or. line(1:1) == " " .or. line(1:1) == "!") then
            continue
          else
            select case(line(1:5))

            ! start value to begin reading lines from
            case('START')
              read(line(35:len_trim(line)),*) startingTime
              if (Verbose)print *,'firsttime',startingTime
              FIRSTTIME = startingTime
            case('ENDIN')
              read(line(35:len_trim(line)),*) endingTime
              if (Verbose)print *,'endtime',endingTime
              LASTTIME = endingTime

            ! verbosity
            case('VERBO')
              read(line(35:len_trim(line)),*) beVerbose
              VERBOSE = beVerbose

            ! file names
            case('REFER')
              read(line(35:len_trim(line)),'(A128)') tmp
              fileReference = trim(tmp)
              i = i+1
              if (Verbose)print *,'reference file:',fileReference

            case('PERTU')
              read(line(35:len_trim(line)),'(A128)') tmp
              fileDelta = trim(tmp)
              i = i+1
              if (Verbose)print *,'perturbation file:',fileDelta

            ! file output directory
            case('DATAD')
              read(line(35:len_trim(line)),*) tmp
              datadirectory = trim(tmp)
              if (Verbose)print *,'data output directory : '//&
                              datadirectory(1:len_trim(datadirectory))

            ! wave type e.g. L150
            case('CPHAS')
              if (FILTERSEISMOGRAMS) then
                read(line(35:len_trim(line)),*) cphasetype
                cphasetype = trim(cphasetype)
                call determinePhaseRef(cphasetype,8,cphaseRef)
                if (Verbose)print *,'cphase    ', cphasetype, cphaseRef
              endif
            end select
          endif
        endif
      enddo
      close(10)

      ! check number of files
      if (i /= 2) then
        print *,'Error: reading input parameters'
        stop
      endif

      end subroutine


!-----------------------------------------------------------------------
      subroutine determineFileLengthTaped(fileName,start,ending,entries,lasttime)
!-----------------------------------------------------------------------
! reads in the seismogram
!
! input:
!     seismo    - data array
!     size        - data array size
!     fileName - data file name
!     start         - time position from which on data shall be read
!     timediff              - time step dt
!     entries       - (optional) number of entries
! returns: seismo array and dt
      use verbosity; use precisions
      implicit none
      character(len=128),intent(in):: fileName
      real(WP),intent(in):: start,ending
      integer,intent(out):: entries
      real(WP),intent(out):: lasttime
      ! local parameters
      integer:: ier,i
      character(len=128):: line
      real(WP):: time,displace,timediff,endtime
      integer:: index, offset,ZEROLINE,FILELINES !tests: ZEROLINE/0/ ! approximate time 0 s line

      !initialize
      ZEROLINE = 0
      FILELINES = 0

      ! open seismogram file
      if (beVerbose) then
        print *,'  file: ',trim(fileName)
      endif

      open(10,file=trim(fileName),status='old',iostat=ier)
      if (ier /= 0) then
        print *,'Error: opening file ',trim(fileName)
        stop 'abort - determineFileLength taped'
      endif

      ! parse file for line with zero time and number of lines
      index = 0
      offset = 0
      do while( ier == 0 )
        read(10,*,iostat=ier) line
        index = index+1
        line=trim(line)

        ! check for additional data in file
        if (line(1:3) == 'dis') offset=3
        if (line(1:3) == 'rec') offset=2

        ! check for line with time 0.0
        if (line(1:3) == '0.0') ZEROLINE= index

      enddo
      FILELINES = index-1
      if (start < 0.0) ZEROLINE= offset
      rewind(10)

      ! parse file for displacement values
      timediff = 0.0
      index = 0
      ier = 0
      do i = 1, FILELINES
        if (i >= ZEROLINE) then
          read(10, *, iostat=ier) time,displace !,sourceterm
          if (ier /= 0) then
            print *,'Error: reading input. last line ',i,ZEROLINE,index
            stop 'abort - determineFileLength taped'
          endif

          if (time-start >= 0.0 .and. ending-time >= 0.0) then ! starts with times from 'start time' on
            index = index+1
            endtime = time
          endif
        else
          !read until end is reached
          read(1, *, iostat=ier) time,displace !,sourceterm
          endtime = time
        endif
      enddo
      close(10)

      if (beVerbose) then
        print *,'    entry read:',index
        print *,'    last time readin:',endtime
      endif

      ! number of file entries read in
      entries = index

      ! end time in file
      lasttime = endtime

      end subroutine

!-----------------------------------------------------------------------
      subroutine stopProgram_fixedlength(textline)
!-----------------------------------------------------------------------
      use parallel
      implicit none
      character(len=128),intent(in):: textline
      ! local parameters
      integer:: endindex,ier
      logical:: flag

      ! console output
      endindex = index(textline,"  ")
      if (endindex < 1) endindex = 128
      print *,textline(1:endindex)

      ! on linux machines : i/o unit 6 is the stdout , on SGI 101
      flush(6)

      ! stop MPI
      call MPI_Initialized(flag,ier)
      if (flag .eqv. .true. .and. ier == 0) then
        ! note: MPI_ABORT does not return, it makes the program exit with an error code of 30
        call MPI_Abort(MPI_COMM_WORLD, 30, ier )
        call MPI_FINALIZE(ier)
      endif

      ! stop process
      stop 'abort - program '

      end subroutine

