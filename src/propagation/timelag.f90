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
      logical,parameter:: DEBUG_OUTPUT = .false.  ! output debug files
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

      ! timelag executable should be run as single process only
      MAIN_PROCESS = .true.

      ! for debuging
      fileOutput = DEBUG_OUTPUT
      beVerbose = DEBUG_OUTPUT

      ! read in from startingTime
      if (VERBOSE) then
        print *
        print *,'determine file length...'
      endif

      call determineFileLengthTaped(fileReference,startingTime,endingTime,entries,lasttime)
      numofTimeSteps = entries

      if (VERBOSE) then
        print *
        print *,'determine FFT parameters...'
      endif

      ! fft window size
      call determineFFTWindowsize(numofTimeSteps,WindowSIZE)

      ! determine filter bandwidth parameters (waveperiod and frequency)
      call determineFilterParameters()

      ! get time lag
      if (VERBOSE) then
        print *
        print *,'calculating time lag...'
      endif

      call getTimelagTaped(t_lag,fileDelta,fileReference,startingTime,endingTime)

      ! time lag is the subsample position times the time step size
      ! if seismo lags seismoRef, i.e., is shifted to the right of it, then ans will
      ! show a peak at positive lags
      if (VERBOSE) print *
      print *,'  numerical timelag          : ',t_lag
      !print *,'  numerical phase anomaly (s):',t_lag/(2*PI/bw_waveperiod)
      if (VERBOSE) print *

      ! compare with analytical formula
      if (ANALYTICAL_CORRELATION) then
        fileOutput = DEBUG_OUTPUT
        if (VERBOSE) then
          print *,'calculating analytical time lag...'
        endif

        call getAnalyticalTimelagTaped(analytict_lag,fileDelta,fileReference,startingTime,endingTime)

        if (VERBOSE) print *
        print *,'  analytic timelag           : ',analytict_lag
        if (VERBOSE) print *
        !print *,'  analytic phase anomaly (s):',t_lag/(2*PI/bw_waveperiod)
      endif

      ! for debugging
      fileOutput = DEBUG_OUTPUT
      beVerbose = DEBUG_OUTPUT

      ! get time lag
      if (VERBOSE) then
        print *,'downhill simplex time lag...'
      endif

      call getMinimized(t_lag,amplification,fileDelta,fileReference,startingTime,endingTime)

      if (VERBOSE) print *
      print *,'  nonlinear timelag          : ',t_lag
      print *,'  nonlinear amplification    : ',amplification
      if (VERBOSE) print *

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
      integer:: i,istart,iend,ier

      ! open input parameter file
      i = 0
      inputName = 'Timelag_Input'
      open(IIN,file=trim(inputName),status='old',iostat=ier)
      if (ier /= 0) then
        print *,'Error: opening file ',trim(inputName)
        stop 'Abort - opening input'
      endif

      ! parse file for parameters
      do while( ier == 0)
        read(IIN,'(A128)',iostat=ier) line
        if (ier /= 0) exit

        ! suppress leading white spaces, if any
        line = adjustl(line)

        ! suppress trailing carriage return (ASCII code 13) if any (e.g. if input text file coming from Windows/DOS)
        if (index(line,achar(13)) > 0) line = line(1:index(line,achar(13))-1)

        ! suppress trailing comments " .. ! my comment"
        if (index(line,'#') > 5) line = line(1:index(line,'#')-1)
        if (index(line,'!') > 5) line = line(1:index(line,'!')-1)

        line = trim(line)

        ! check line
        if (len_trim(line) == 0) cycle

        ! check if comment line
        if (line(1:1) == " " .or. line(1:1) == "!" .or. line(1:1) == "%") cycle

        ! get index of "=" sign
        istart = index(line,'=')
        if (istart == 0) cycle

        ! check if parameter name valid
        if (istart < 5) then
          print *,'Error: line with wrong format         : ***'//trim(line)//'***'
          print *,'       start index for parameter value: istart = ',istart
          call stopProgram('Abort - line with wrong format in input file    ')
        endif

        ! index range for parameter values
        istart = istart + 1
        iend = len_trim(line)

        !debug
        !print *,'debug: index range is/ie = ',is,'/',ie
        !print *,'debug: parameter string  = ***'//line(is:ie)//'***'

        select case(line(1:5))
        ! start value to begin reading lines from
        case('START')
          read(line(istart:iend),*) startingTime
          if (Verbose) print *,'firsttime        : ',startingTime
          FIRSTTIME = startingTime
        case('ENDIN')
          read(line(istart:iend),*) endingTime
          if (Verbose) print *,'endtime          : ',endingTime
          LASTTIME = endingTime

        ! verbosity
        case('VERBO')
          read(line(istart:iend),*) beVerbose
          VERBOSE = beVerbose

        ! file names
        case('REFER')
          read(line(istart:iend),'(A128)') tmp
          fileReference = trim(tmp)
          ! suppress leading white spaces, if any
          fileReference = adjustl(fileReference)
          i = i+1
          if (Verbose) print *,'reference file   : ',trim(fileReference)

        case('PERTU')
          read(line(istart:iend),'(A128)') tmp
          fileDelta = trim(tmp)
          ! suppress leading white spaces, if any
          fileDelta = adjustl(fileDelta)
          i = i+1
          if (Verbose) print *,'perturbation file: ',trim(fileDelta)

        ! file output directory
        case('DATAD')
          read(line(istart:iend),*) tmp
          datadirectory = trim(tmp)
          ! suppress leading white spaces, if any
          datadirectory = adjustl(datadirectory)
          if (Verbose) print *,'data output      : ',trim(datadirectory)

        ! wave type e.g. L150
        case('CPHAS')
          if (FILTER_SEISMOGRAMS) then
            read(line(istart:iend),*) cphasetype
            cphasetype = trim(cphasetype)
            ! suppress leading white spaces, if any
            cphasetype = adjustl(cphasetype)
            ! get corresponding phase velocity
            call determinePhaseRef(cphasetype,8,cphaseRef)
            if (Verbose) print *,'cphase         : ', cphasetype, cphaseRef
          endif
        end select
      enddo
      close(IIN)

      if (Verbose) print *,'-----------------------------------------------------'

      ! check number of files
      if (i /= 2) then
        print *,'Error: reading input parameters'
        stop 'Error reading input parameters'
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

      open(IIN,file=trim(fileName),status='old',iostat=ier)
      if (ier /= 0) then
        print *,'Error: opening file ',trim(fileName)
        stop 'Abort - determineFileLength taped'
      endif

      ! parse file for line with zero time and number of lines
      index = 0
      offset = 0
      do while( ier == 0 )
        read(IIN,*,iostat=ier) line
        index = index+1
        line = trim(line)

        ! check for additional data in file
        if (line(1:3) == 'dis') offset = 3
        if (line(1:3) == 'rec') offset = 2

        ! check for line with time 0.0
        if (line(1:3) == '0.0') ZEROLINE = index

      enddo
      FILELINES = index-1
      if (start < 0.0) ZEROLINE = offset

      rewind(IIN)

      ! parse file for displacement values
      timediff = 0.0
      index = 0
      ier = 0
      do i = 1, FILELINES
        if (i >= ZEROLINE) then
          read(IIN, *, iostat=ier) time,displace !,sourceterm
          if (ier /= 0) then
            print *,'Error: reading input. last line ',i,ZEROLINE,index
            stop 'Abort - determineFileLength taped'
          endif

          if (time-start >= 0.0 .and. ending-time >= 0.0) then ! starts with times from 'start time' on
            index = index+1
            endtime = time
          endif
        else
          !read until end is reached
          read(IIN, *, iostat=ier) time,displace !,sourceterm
          endtime = time
        endif
      enddo
      close(IIN)

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
      integer:: endindex

      ! console output
      endindex = index(textline,"  ")
      if (endindex < 1) endindex = 128
      print *,textline(1:endindex)

      ! on linux machines : i/o unit 6 is the stdout , on SGI 101
      flush(6)

      ! stop MPI
      call syncAbortMPI()

      ! stop process
      stop 'Abort - program '

      end subroutine

