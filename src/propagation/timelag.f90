!-----------------------------------------------------------------------
      program timelag
!-----------------------------------------------------------------------
! calculates the time lag between two seismograms
      use precision; use verbosity; use filterType;use traveltime; use propagationStartup
      use parallel
      implicit none
      double precision startingTime
      character*128 fileDelta,fileReference
      integer entries
            
      !-----------------------------------------------------------------------
      ! parameters      
      ! most parameters concerning timelag calculation are set in file Timelag_Input 
      ! (& default values in commonModules.f90)
      !-----------------------------------------------------------------------

      ! file output for debugging ('tmpTT*.dat')
      fileOutput       = .false.         

      ! this is a master process
      MASTER = .true.
      
      ! console output
      print*,'Timelag'
        
      ! get input parameters
      print*
      print*,'reading Input-Parameters...'
      call readInputParameters(fileDelta,fileReference,startingTime)
      
      ! read in from startingTime
      print*
      print*,'determine file length...'
      call determineFileLength(fileReference,startingTime,entries,lasttime)
      numofTimeSteps=entries

      ! fft window 
      call determineFFTWindowsize()

      ! determine filter bandwidth parameters (waveperiod and frequency)
      call determineFilterParameters()
      
      ! get time lag
      print*
      print*,'calculating time lag...'
      call getTimelag(t_lag,fileDelta,fileReference,startingTime)

      ! reset
      fileOutput=.true.
      beVerbose=.true.
        
      ! time lag is the subsample position times the time step size
      ! if seismo lags seismoRef, i.e., is shifted to the right of it, then ans will show a peak at positive lags
      if(BEVERBOSE) then
        print*,'timelag:',t_lag
        print*
      else
        print*,t_lag
      endif
      
      ! compare with analytical formula
      if( ANALYTICAL_CORRELATION ) then
        print*,'calculating analytical time lag...'
        call getAnalyticTimelag(t_lag,fileDelta,fileReference,startingTime,bw_width,BUTTERWORTH_POWER,bw_waveperiod,FILTERSEISMOGRAMS)
        print*,'analytic timelag:',t_lag
      endif
      
      end


!-----------------------------------------------------------------------
      subroutine readInputParameters(fileDelta,fileReference,startingTime)
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
      character*128:: line,fileDelta,fileReference
      character*128:: inputName,tmp
      integer:: i,ierror,length
      double precision:: startingTime
      
      ! open input parameter file
      i=0
      inputName = 'Timelag_Input' 
      inputName = trim(inputName)
      open(1,file=inputName,status='old',iostat=ierror)
      if( ierror .ne. 0) then
        print*,'error opening file Timelag_Input'
        stop 'abort - opening input'
      endif
      
      ! parse file for parameters
      do while( ierror .eq. 0)
        read(1,'(A128)',iostat=ierror) line
        if( ierror .ne. 0 ) exit
        
        length = len_trim(line)
        if( length .eq. 0 ) then
          continue
        else
          if( line(1:1) .eq. "%" .or. line(1:1) .eq. " " .or. line(1:1) .eq. "!" ) then 
            continue
          else
            select case( line(1:5) )

            ! start value to begin reading lines from
            case('START')
              read(line(35:len_trim(line)),*) startingTime
              if(beVerbose)print*,'firsttime',startingTime

            ! verbosity
            case('VERBO')
              read(line(35:len_trim(line)),*) beVerbose

            ! file names
            case('REFER')
              read(line(35:len_trim(line)),'(A128)') tmp            
              fileReference = trim(tmp)
              i=i+1
              if(beVerbose)print*,'reference file:',fileReference

            case('PERTU')
              read(line(35:len_trim(line)),'(A128)') tmp            
              fileDelta = trim(tmp)
              i=i+1
              if(beVerbose)print*,'perturbation file:',fileDelta
              
            ! file output directory  
            case('DATAD')
              read(line(35:len_trim(line)),*) tmp
              datadirectory = trim(tmp)
              if(beVerbose)print*,'data output directory : '//datadirectory(1:len_trim(datadirectory))

            ! wave type e.g. L150
            case('CPHAS')
              read(line(35:len_trim(line)),*) cphasetype
              cphasetype = trim(cphasetype)
              call determinePhaseRef()
              if(beVerbose)print*,'cphase    ', cphasetype, cphaseRef              
              
            end select
          endif
        endif
      enddo

      ! check number of files
      if( i .ne. 2) then
        print*,'error reading input parameters'
        stop
      endif

      end


!-----------------------------------------------------------------------
      subroutine determineFileLength(fileName,start,entries,lasttime)
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
      use verbosity
      implicit none
      integer:: size,ierror,i,entries
      character*128:: fileName,line
      double precision:: time,sourceterm,displace,timediff,lasttime,start,endtime
      integer:: index, offset,ZEROLINE,FILELINES !tests: ZEROLINE/0/ ! approximate time 0 s line

      !initialize
      ZEROLINE=0
      FILELINES=0

      ! open seismogram file
      fileName = trim(fileName)
      if(beVerbose) then 
        print*,'    ',fileName(1:len_trim(fileName))   
      endif
      open(1, file= fileName,status='old',iostat=ierror)
      if( ierror .ne. 0) then
        print*,'error opening file ',fileName
        stop 'abort - determineFileLength'
      endif
      
      ! parse file for line with zero time and number of lines
      ierror=0
      index=0
      offset=0
      do while( ierror.eq. 0 )
        read(1,*,iostat=ierror) line
        index=index+1
        line=trim(line)
        
        ! check for additional data in file
        if( line(1:3) .eq. 'dis') offset=3
        if( line(1:3) .eq. 'rec') offset=2
        
        ! check for line with time 0.0
        if( line(1:3) .eq. '0.0') ZEROLINE= index
        
      enddo
      FILELINES=index-1      
      if( start .lt. 0.0) ZEROLINE= offset
      rewind(1)
      
      ! parse file for displacement values      
      timediff=0.0
      index=0
      ierror=0
      do i=1, FILELINES
        if(i .ge. ZEROLINE ) then
          read(1, *, iostat=ierror) time,displace !,sourceterm
          if( ierror .ne. 0) then
            print*,'error reading input. last line ',i,ZEROLINE,index
            stop 'abort - determineFileLength'
          endif
          
          !debug
          if(DEBUG)print*,time,displace,sourceterm
          
          if( time - start .gt. 0.0) then ! starts with times from 'start time' on
            index=index+1
            endtime=time            
          endif
        else        
          !read until end is reached
          read(1, *, iostat=ierror) time,displace !,sourceterm
          endtime=time
        endif
      enddo
      close(1)
      
      if( beVerbose ) then
        print*,'    entry read:',index
        print*,'    last time readin:',endtime
      endif
      
      ! number of file entries read in
      entries=index
      
      ! end time in file
      lasttime=endtime
      
      end
