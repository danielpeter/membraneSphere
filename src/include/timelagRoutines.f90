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
      subroutine getTimelag(t_lag,fileDelta,fileReference,startingTime)
!-----------------------------------------------------------------------
! determines the time lag between two seismograms given 
! their filenames
! (both seismograms must have the same time step size)
!
! definition: fileDelta shifted to the right of fileReference, then timelag is positiv (i.e. waves arrive later in fileDelta)
!
! input:
!       t_lag                                 - time lag 
!       fileDelta,fileReference      - file names of seismograms
!       startingTime                     - window start time of fourier transformation for correlation
!
! returns: t_lag timelag in seconds
      use verbosity; use nrtype; use nrutil; use filterType
      use nr, ONLY : correl            
      implicit none
      character*128:: fileDelta,fileReference
      integer::i,entries,station,index,xcorrlength,zeropad,ierror
      real(WP):: seismo(2,WindowSIZE),seismoRef(2,WindowSIZE),t_lag,startingTime,lastseismotime,timedelta
            
      !add zero for non-periodic data, influences precision
      !zeropad=WindowSIZE 
      
      ! length of cross-correlation array      
      xcorrlength=WindowSIZE !or: xcorrlength=WindowSIZE+zeropad
      
      ! attention of seismo arrays: zero padding of length WindowSIZE at the end to get corresponding length of cross-correlation
      !allocate(seismo(2,xcorrlength),seismoRef(2,xcorrlength),stat=ierror)
      !if( ierror .ne. 0) call stopProgram( 'abort - getTimelag      ')
      
      !initialize
      seismo(:,:)=0.0_WP
      seismoRef(:,:)=0.0_WP
      
      ! read in from startingTime
      call readSeismo(seismoRef,WindowSIZE,fileReference,startingTime,timedelta,entries,lastseismotime)
      call readSeismo(seismo,WindowSIZE,fileDelta,startingTime,timedelta,entries,lastseismotime)
            
      ! get time lag
      call getTimelagSeismos(t_lag,seismo,seismoRef,entries,startingTime,timedelta)      

      end
      
!-----------------------------------------------------------------------
      subroutine getTimelagRef(t_lag,station,startingTime)
!-----------------------------------------------------------------------
! determines the time lag between two filtered seismograms given 
! their filenames
! (both seismograms must have the same time step size)
!
! definition: fileDelta shifted to the right of fileReference, then timelag is positiv (i.e. waves arrive later in fileDelta)
!
! input:
!       t_lag                                 - time lag 
!       station                              - index of reference station 
!       startingTime                     - window start time of fourier transformation for correlation
!       doFiltering                        - use filter before determining time lag
!
! returns: t_lag timelag in seconds
      use verbosity;use nrtype; use nrutil; use propagationStartup;use filterType
      use nr, ONLY : correl; use verbosity            
      implicit none
      integer::i,station,index,xcorrlength,zeropad,ierror
      real(WP):: seismo(2,WindowSIZE),seismoRef(2,WindowSIZE),t_lag,startingTime,endtime,time,timeRef,displace,displaceRef
            
      !add additional zero for non-periodic data, influences precision
      !zeropad=WindowSIZE 
      
      ! length of cross-correlation array      
      xcorrlength=WindowSIZE !or: xcorrlength=WindowSIZE+zeropad
      
      ! attention of seismo arrays: zero padding of length WindowSIZE at the end to get corresponding length of cross-correlation
      !allocate(seismo(2,xcorrlength),seismoRef(2,xcorrlength),stat=ierror)
      !if( ierror .ne. 0) call stopProgram('abort - getTimelagRef    ')
      
      !initialize
      seismo(:,:)=0.0_WP
      seismoRef(:,:)=0.0_WP
            
      ! read in from startingTime
      index=0
      do i=1,numofTimeSteps  
        if(index .lt. WindowSIZE) then
          if( manyKernels ) then
            time=kernelsReceiversSeismogram(currentKernel,numofReceivers+1,i)
            displace=kernelsReceiversSeismogram(currentKernel,station,i)
            timeRef=kernelsReceiversSeismogramRef(currentKernel,numofReceivers+1,i)
            displaceRef=kernelsReceiversSeismogramRef(currentKernel,station,i)
          else
            time=receiversSeismogram(size(receivers)+1,i)
            displace=receiversSeismogram(station,i)
            timeRef=receiversSeismogramRef(size(receivers)+1,i)
            displaceRef=receiversSeismogramRef(station,i)          
          endif
          !debug
          if(DEBUG) print*,time,displace
          
          if( time - startingTime .gt. 0.0_WP) then ! starts with times from 'start time' on
            index=index+1
            seismo(1,index) = time
            seismo(2,index) = displace
            seismoRef(1,index) = timeRef
            seismoRef(2,index) = displaceRef
            endtime=time            
          endif
        else        
          !read until end is reached
          time=receiversSeismogram(size(receivers)+1,i) !,sourceterm
          endtime=time
        endif
      enddo
      
      ! get time lag
      call getTimelagSeismos(t_lag,seismo,seismoRef,index,startingTime,dt)      
      end

!-----------------------------------------------------------------------
      subroutine getTimelagSeismos(t_lag,seismo1,seismo2,seismoLength,startingTime,timedelta)
!-----------------------------------------------------------------------
! determines the time lag between two seismograms
! (both seismograms must have the same time step size)
!
! definition: seismo1 shifted to the right of seismo2, then timelag is positiv (i.e. waves arrive later in seismo1)
!
! input:
!       t_lag                                 - time lag 
!       seismo1,seismo2             - seismograms (array 2 x seismoLength, first index is time, second displacement)
!       seismoLength                  - seismograms length (number of entries in seismogram)
!       startingTime                     - window start time of fourier transformation for correlation
!       doFiltering                        - use filter before determining time lag
!
! remember:
!       bw_width,BUTTERWORTH_POWER         - filter parameters
!       bw_waveperiod                - wave period in seconds of butterworth bandwidth frequency
! returns: t_lag timelag in seconds ( uses dt )
      use verbosity;use filterType; use propagationStartup
      use nrtype; use nrutil; use nr, ONLY : correl            
      implicit none
      logical:: doFiltering
      integer:: seismoLength
      integer:: indexmax,ierror,i,j,entries,station,index
      integer:: hannwindow,ileft,iright,xcorrlength,zeropad,centerfrequencyIndex
      real(WP):: t_lag, startingTime,timedelta
      real(WP):: seismo(2,WindowSIZE),seismoRef(2,WindowSIZE),seismo1(2,seismoLength),seismo2(2,seismoLength)
      real(WP):: crosscorrelation(WindowSIZE),seismoWindow(WindowSIZE),seismoRefWindow(WindowSIZE)
      real(WP):: sourceterm,correlation,max,sub, pointleft,pointright,samplingFreq,endtime
      real(WP):: getMaximum,hannfactor,hannA,time,timeRef,displace,displaceRef
      external:: getMaximum

      ! filtering necessary
      doFiltering=FILTERSEISMOGRAMS
            
      !add zero for non-periodic data, influences precision
      !zeropad=WindowSIZE 
      
      ! length of cross-correlation array
      xcorrlength=WindowSIZE !or:   xcorrlength=WindowSIZE+zeropad
      
      ! attention of seismo arrays: zero padding of length WindowSIZE at the end to get corresponding length of cross-correlation
      !allocate(seismo(2,xcorrlength),seismoRef(2,xcorrlength),crosscorrelation(xcorrlength),&
      !       seismoWindow(xcorrlength),seismoRefWindow(xcorrlength), stat=ierror)
      !if( ierror .ne. 0) call stopProgram( 'abort - getTimelagSeismos      ')
      
      !initialize
      seismo(:,:)=0.0_WP
      seismoRef(:,:)=0.0_WP
      seismoWindow(:)=0.0_WP
      seismoRefWindow(:)=0.0_WP
      
      ! read in seismograms from startingTime
      call captureSeismos(seismo1,seismo2,seismo,seismoRef,&
                  seismoLength,xcorrlength,startingTime,entries,endtime)
      
      if( beVerbose ) then
        print*      
        print*,'perturbed seismogram:'
        print*,'    dt:',timedelta
        print*,'    first entry:',seismo(1,1),seismo(2,1)
        print*,'    last entry:',seismo(1,entries),seismo(2,entries)
        print*,'    entry read:',entries
        print*,'    last time readin:',seismo(1,entries),'seismo:',endtime
        print*,'reference seismogram:'
        print*,'    first entry:',seismoRef(1,1),seismoRef(2,1)
        print*,'    last entry:',seismoRef(1,entries),seismoRef(2,entries)
        print*,'    entry read:',entries
        print*,'    last time readin:',seismoRef(1,entries)
        print*
      endif
                  
      ! debug output
      if( fileOutput) then
        open(10,file=datadirectory(1:len_trim(datadirectory))//'tmpgetTimelag_read.dat')
        do i=1,xcorrlength
          write(10,*) seismo(1,i),seismo(2,i)
        enddo
        close(10)
      endif
      
      if( doFiltering) then
        ! apply hanning window  to smooth seismograms ends
        call taperSeismos(seismo,seismoRef,xcorrlength,entries)
      
        ! smallest sampled frequency
        samplingFreq=1.0_WP/(xcorrlength*timedelta)
      
        ! position index of bandwidth frequency for given waveperiod
        centerfrequencyIndex=nint( 1.0_WP/((bw_waveperiod + MEMBRANECORRECTION)*samplingFreq) )           

        ! determine bandwidth integer depending on windowsize and requested bandwidth-frequency
        ! integer is a factor of 2 to have same number of frequencies filtered around the center frequency to the left and the right
        !bw_width= 2*abs(centerfrequencyIndex-nint(bw_frequency*dt*xcorrlength))
        bw_width= 2*nint(bw_frequency*timedelta*xcorrlength)
        
        ! console output
        if( beVerbose ) then
          print*
          print*,'filter parameter:'
          print*,'    sampled frequency:',samplingFreq
          print*,'    wave period:',bw_waveperiod+MEMBRANECORRECTION
          print*,'        initial period:',bw_waveperiod
          print*,'        membrane correction:',MEMBRANECORRECTION
          print*,'    center frequency:',1.0_WP/(bw_waveperiod+MEMBRANECORRECTION) 
          print*,'        index:',centerfrequencyIndex
          print*,'    bandwidth:', bw_width,' power:',BUTTERWORTH_POWER
          print*,'        as frequency:',bw_width*samplingFreq
          print*,'    fourier window start:',ARRIVAL_THRESHOLD
          print*,'    windowsize:',xcorrlength
          print*
        endif

        ! check if filter frequency is correct
        if( centerfrequencyIndex .lt. 1 .or. centerfrequencyIndex .gt. xcorrlength/2 ) then
          print*,'filter is incorrect!',centerfrequencyIndex,xcorrlength
          call stopProgram('abort filterSeismogram() - incorrect centerfrequencyIndex    ')
        endif
                
        ! filter seismograms in frequency domain ( e.g. with a butterworth filter around L150s frequency and bandwidth of 10*dfreq)
        call filterSeismos(seismo(2,:),seismoRef(2,:),xcorrlength,bw_width,BUTTERWORTH_POWER,centerfrequencyIndex)
        
        ! take window of seismogram data (no window function applied)
        ! avoid taking temporary-working data in second half of fourier-transformed seismogram
        !j=0
        !do i=1,seismoLength
        !  if( seismo(1,i) - startingTime .gt. 0.0) then      
        !    j=j+1
        !    seismoWindow(j)=seismo(2,i)
        !    seismoRefWindow(j)=seismoRef(2,i)
        !  endif
        !enddo
        seismoWindow(:)=seismo(2,:)
        seismoRefWindow(:)=seismoRef(2,:)        
      else
        ! just get initial seismograms
        seismoWindow(:)=seismo(2,:)
        seismoRefWindow(:)=seismoRef(2,:)
      endif

      ! normalize seismograms
      call normalizeArray(seismoWindow,xcorrlength)
      call normalizeArray(seismoRefWindow,xcorrlength)
      
      ! debug output
      if( fileOutput ) then
        open(10,file=datadirectory(1:len_trim(datadirectory))//'tmpgetTimelag_window.dat')
        j=0        
        do i=1,xcorrlength
          if( seismo(1,i) - startingTime .gt. 0.0) then
            if(j .eq. 0) then
              j=i
            endif
          endif
        enddo
        do i=1,xcorrlength
          write(10,*) seismo(1,1)+(j+i-2)*timedelta,seismoWindow(i),seismoRefWindow(i)
        enddo
        close(10)
      endif
      
      ! get cross-correlation (by numerical recipes) 
      crosscorrelation=correl(seismoWindow,seismoRefWindow)
      
      ! debug output
      if( fileOutput ) then
        open(1, file=datadirectory(1:len_trim(datadirectory))//'tmpgetTimelag_correlation.dat')
        !print*
        !print*,'autocorrelation'
        write(1,*) 'correlation'
        do i=1,xcorrlength
          !print*,i,ans(i)
          write(1,*) i,crosscorrelation(i)
        enddo
        close(1)
      endif
      
      ! calculate the subsample precision maximum position
      sub=getMaximum(crosscorrelation,size(crosscorrelation))
      
      !only be verbose once
      if( beVerbose ) then
        print*
        beVerbose=.false.
        fileOutput=.false.
      endif
      
      ! set timelag
      t_lag = sub*timedelta
      end


!-----------------------------------------------------------------------
      function getMaximum(crosscorrelation,length)
!-----------------------------------------------------------------------
! determines the maximum position with subsample precision
!
! input:
!       crosscorrelation                            - cross correlation values
!       length                                           - array size
!
! returns: subsample precision maximum
      use verbosity;use precision
      implicit none
      integer:: length,i,indexmax
      real(WP):: crosscorrelation(length),getMaximum,max,correlation
      real(WP):: pointleft,pointright,maxestimate
      
      !Search maximum correlation      
      max = abs(crosscorrelation(1))
      indexmax = 1
      do i=2,length
        correlation=abs(crosscorrelation(i))
        if( correlation.gt.max ) then
          max = correlation
          indexmax = i
        endif
      enddo
      max = crosscorrelation(indexmax)  ! could be negative maximum
      
      !if (indexmax.gt.WindowSIZE/2) then
      !  indexmax=indexmax-WindowSIZE
      !else
      !  indexmax=indexmax-1
      !endif
      if( beVerbose) then
        print*,'    maximum correlation:',max,indexmax
      endif
      
      !subsample precision
      ! first half (positive lag)
      if( indexmax .gt. 1 .and. indexmax .lt. length/2) then
        pointleft=crosscorrelation(indexmax-1)
        pointright=crosscorrelation(indexmax+1)
      endif
      if( indexmax .eq. 1 ) then
        pointleft=crosscorrelation(length)
        pointright=crosscorrelation(indexmax+1)
      endif
      if( indexmax .eq. length/2 ) then
        pointleft=crosscorrelation(indexmax-1)
        pointright=crosscorrelation(length/2+1)
      endif
      ! second half (negative lag)
      if( indexmax.gt.length/2+1.and.indexmax.lt.length)then
        pointleft=crosscorrelation(indexmax-1)
        pointright=crosscorrelation(indexmax+1)
      endif      
      if( indexmax .eq. length/2+1 ) then
        pointleft=crosscorrelation(length/2)
        pointright=crosscorrelation(indexmax+1)
      endif
      if( indexmax .eq. length ) then
        pointleft=crosscorrelation(indexmax-1)
        pointright=crosscorrelation(1)
      endif
            
      !debug
      if( beVerbose ) then
        print*,'        pointleft:',pointleft
        print*,'        pointright:',pointright
      endif

      ! respect wrap-around order
      if( indexmax.gt.length/2 ) then
        indexmax=indexmax-length-1
      else
        indexmax=indexmax-1
      endif      
      
      ! quadratic interpolation will give this maximum
      ! ( http://www-ccrma.stanford.edu/~jos/parshl/Peak_Detection_Steps_3.html )
      !getMaximum = indexmax - 0.5*(pointright-pointleft)/(pointleft - 2.0*max + pointright)
      getMaximum = indexmax + 0.5*(pointleft-pointright)/(pointleft - 2.0*max + pointright)
      if( beVerbose ) then
        maxestimate=max-0.25*(pointleft-pointright)*getMaximum
        print*,'    subsample position:',getMaximum
        print*,'        index:',indexmax
        print*,'        estimated maximum:',maxestimate        
      endif
      return
      end


!-----------------------------------------------------------------------
      subroutine readSeismo(seismo,size,fileName,start,timediff,entries,lastseismotime)
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
      use verbosity;use precision
      implicit none
      integer:: size,ierror,i,entries
      character*128:: fileName,line
      real(WP):: seismo(2,size),time,sourceterm,displace,timediff,lastseismotime,start,endtime
      integer:: index, offset,ZEROLINE,FILELINES !tests: ZEROLINE/0/ ! approximate time 0 s line

      !initialize
      ZEROLINE=0
      FILELINES=1100

      ! open seismogram file
      fileName = trim(fileName)
      if(beVerbose) then 
        print*
        print*,fileName   
        print*,'start:',start
      endif
      open(1, file= fileName,status='old',iostat=ierror)
      if( ierror .ne. 0) then
        print*,'error opening file ',fileName
        call stopProgram( 'abort - readSeismo   ')
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
      if( beVerbose ) then
        print*,'    file lines:',FILELINES
        print*,'    zero line:',ZEROLINE
        print*,'    offset:',offset
      endif
      rewind(1)
      
      ! parse file for displacement values      
      timediff=0.0
      index=0
      ierror=0
      do i=1, FILELINES
        if(i .ge. ZEROLINE .and. index .lt. size) then
          read(1, *, iostat=ierror) time,displace !,sourceterm
          if( ierror .ne. 0) then
            print*,'error reading input. last line ',i,ZEROLINE,index
            call stopProgram( 'abort - readSeismo   ')
          endif
          
          !debug
          if(DEBUG)print*,time,displace,sourceterm
          
          if( time - start .gt. 0.0) then ! starts with times from 'start time' on
            index=index+1
            seismo(1,index) = time
            seismo(2,index) = displace
            !if(beVerbose.and.(index.eq.1.or.index.eq.size))               &
            !&                print*,index,seismo(:,index)
            endtime=time
            
            !time step
            if( abs(timediff) .lt. 0.000001 .and. index .gt. 1) then
              timediff = seismo(1,index)-seismo(1,index-1)
            endif
          endif
        else        
          !read until end is reached
          if( index .ge. size) then
            !read values
            read(1, *, iostat=ierror) time,displace !,sourceterm
            endtime=time
          else
            !read a text line
            read(1, *, iostat=ierror) line
            if( ierror .ne. 0) then
              !print*,'error reading input. last line ',i
              exit
            endif                  
          endif
        endif
      enddo
      close(1)
      
      if( beVerbose ) then
        print*,'    timestep dt:',timediff
        print*,'    first entry:',seismo(1,1),seismo(2,1)
        print*,'    last entry:',seismo(1,index),seismo(2,index)
        print*,'    entry read:',index
        print*,'    last time readin:',seismo(1,index),'file:',endtime
      endif
      
      ! number of file entries read in
      entries=index
      
      ! end time in file
      lastseismotime=endtime
      
      end

!-----------------------------------------------------------------------
      subroutine getGridvalues()
!-----------------------------------------------------------------------      
! read in precalculated cell areas, edge lengths and center distances
      use propagationStartup;use cells
      implicit none
      integer::ierror
      
      ! read in vertices array
      call readData()      
      
      ! new arrays for precalculated cell attributes  
      allocate( cellAreas(numVertices),cellEdgesLength(numVertices,0:6), &
              cellCenterDistances(numVertices,0:6), stat=ierror )
      if( ierror .gt. 0 ) then
        print*,'error in allocating arrays for cell area,..'
        call stopProgram( 'abort - getGridvalues    ')
      endif      
                  
      ! read in cell areas
      call readPrecalculated()            
    
      end
 
      
!-----------------------------------------------------------------------
      subroutine getStartTime(filename,startTime)
!-----------------------------------------------------------------------   
! determines the phase arrival time based upon a threshold
! 
! input:
!       filename     - seismogram file
!       startTime   - start time of fourier transform window
!
! returns: arrivalTime (module traveltime) and startTime
      use traveltime; use verbosity; use propagationStartup
      implicit none
      character*128:: filename
      integer:: numEntries,entries1
      real(WP):: seismo(2,10000),startTime,endtime,timedelta,defaultStart
      
      !initialize
      numEntries=10000
      defaultStart=FIRSTTIME
      
      ! read in seismograms (hopefully complete)
      seismo(:,:)=0.0_WP
      call readSeismo(seismo,numEntries,filename,defaultStart,timedelta,entries1,endtime)
      seismo_timestep=timedelta
      ! determine startTime and arrivaltime
      call getStartTimeSeismogram(seismo,numEntries,endtime,timedelta,startTime)
      end


!-----------------------------------------------------------------------
      subroutine getStartTimeSeismogram(seismo,seismoLength,endtime,timedelta,startTime)
!-----------------------------------------------------------------------   
! determines the phase arrival time based upon a threshold
! 
! input:
!       seismo      - seismogram 
!       seismoLength - number of entries in seismogram
!       startTime   - start time of fourier transform window
!       endtime      - end time of seismogram
!       timedelta                - time step
!
! returns: arrivalTime (module traveltime) and startTime
      use traveltime;use verbosity;use propagationStartup; use filterType
      implicit none
      integer:: seismoLength,i
      real(WP),intent(in):: seismo(2,seismoLength),endtime,timedelta
      real(WP),intent(out)::startTime
      real(WP):: leftWindowTime,defaultStart
            
      !initialize
      defaultStart=FIRSTTIME
            
      ! get time when seismogram is above threshold
      arrivalTime=defaultStart
      do i=1,seismoLength
        if(seismo(2,i) .gt. ARRIVAL_THRESHOLD ) then
          arrivalTime = seismo(1,i)
          exit
        endif
      enddo
            
      ! adapt start time for fourier transform depending on arrival time and window size of fourier transformation
      ! to get at least all entries after the arrival time      
      if( seismoLength .gt. WindowSIZE) then
        ! count backward to have starting time of window
        leftWindowTime=endtime-timedelta*WindowSIZE
      
        if( leftWindowTime .lt. arrivalTime ) then
          startTime=leftWindowTime - timedelta/2.0_WP
        else
          if( beVerbose ) then
            print*,'window size not reaching end of seismogram...',arrivalTime,WindowSIZE
            print*
          endif
          ! start when wave arrives
          startTime=arrivalTime - timedelta/2.0_WP
        endif      
      else
        if( beVerbose ) then
          print*
          print*,'    defaultStart time:',defaultStart
        endif
        ! set to default start time
        startTime=defaultStart
      endif
      
      end



!-----------------------------------------------------------------------
      subroutine getAnalyticTimelag(t_lag,fileDelta,fileReference,startingTime,bandwidth,power,waveperiod,doFiltering)
!-----------------------------------------------------------------------
! determines the time lag between two filtered seismograms given 
! their filenames in an analytic approach
! (both seismograms must have the same time step size)
! input:
!       t_lag                             - time lag 
!       fileDelta,fileReference  - file names of seismograms
!       startingTime                      - window start time of fourier transformation for correlation
!       bandwidth,power - butterworth filter parameters
!       waveperiod                      - wave period in seconds of butterworth bandwidth frequency
! returns: t_lag timelag in seconds
      use verbosity; use filterType; use nrtype; use nrutil; use splineFunction
      use nr, ONLY : correl; use propagationStartup            
      implicit none
      character*128:: fileDelta,fileReference
      integer:: indexmax,ierror,i,j,entries
      real(WP):: t_lag,time,displace,sourceterm,correlation,max
      real(WP):: seismo(2,WindowSIZE),seismoRef(2,WindowSIZE)
      real(WP):: crosscorrelation(WindowSIZE),seismoWindow(WindowSIZE),seismoRefWindow(WindowSIZE)
      double complex:: ans(2*WindowSIZE)
      real(WP) startingTime,sub, pointleft,pointright,timedelta,samplingFreq,waveperiod
      integer xcorrlength,zeropad,bw_bfrequency,bandwidth,power,length
      real(DP):: seismo1(WindowSIZE),seismo2(WindowSIZE),seismoPerturbed(WindowSIZE)
      real(WP):: integral1,integral2,lastseismotime
      logical:: doFiltering
            
      !add zero for non-periodic data, influences precision
      !zeropad=WindowSIZE 
      
      ! length of cross-correlation array
      xcorrlength=WindowSIZE !or: xcorrlength=WindowSIZE+zeropad
      
      ! attention of seismo arrays: zero padding of length WindowSIZE at the end to get corresponding length of cross-correlation
      !allocate(seismo(2,xcorrlength),seismoRef(2,xcorrlength),crosscorrelation(xcorrlength),&
      !            seismoWindow(xcorrlength),seismoRefWindow(xcorrlength), &
      !            ans(2*(xcorrlength)),stat=ierror)
      !if( ierror .ne. 0) call stopProgram( 'abort - getAnalyticTimelag      ')
      
      !initialize
      seismo(:,:)=0.0_WP
      seismoRef(:,:)=0.0_WP
      
      if( beVerbose ) print*,'analytical: reading in ',fileOutput,doFiltering
      
      ! read in complete seismograms
      call readSeismo(seismoRef,WindowSIZE,fileReference,-60.0,timedelta,entries,lastseismotime)
      call readSeismo(seismo,WindowSIZE,fileDelta,-60.0,timedelta,entries,lastseismotime)
      
      
      
      if( fileOutput) then
        !open(10,file='tmpTTAnalyticread.'//fileDelta((len(trim(fileDelta))-12):(len(trim(fileDelta)))))
        open(10,file=datadirectory(1:len_trim(datadirectory))//'tmpgetAnalytic_read.dat')
        do i=1,xcorrlength
          write(10,*) seismo(1,i),seismo(2,i)
        enddo
        close(10)
      endif
      
      ! filter
      if( doFiltering) then
        samplingFreq=1.0_WP/(WindowSIZE*timedelta)
        
        ! position index of bandwidth frequency for given waveperiod
        bw_bfrequency=int( (1.0_WP/waveperiod)/samplingFreq )            

        if( beVerbose ) then
          print*
          print*,'butterworth parameter:'
          print*,'    sampled frequency:',samplingFreq
          print*,'    center frequency:',1.0_WP/waveperiod 
          print*,'        index:',bw_bfrequency
          print*,'    bandwidth:',bandwidth
          print*,'        as frequency:',bandwidth*samplingFreq
          print*,'    power:',power
          print*
        endif
        
        ! filter seismograms in frequency domain with a butterworth filter around L150s frequency and bandwidth of 10*dfreq
        call filterSeismos(seismo(2,1:WindowSIZE),seismoRef(2,1:WindowSIZE),WindowSIZE,bandwidth,power,bw_bfrequency)
      
      endif
      
      ! take window of seismogram data (no window weight function applied)
      seismoWindow(:)=0.0_WP
      seismoRefWindow(:)=0.0_WP
      j=1
      do i=1,WindowSIZE
        if( seismo(1,i) - startingTime .gt. 0.0) then      
          seismoWindow(j)=seismo(2,i)
          seismoRefWindow(j)=seismoRef(2,i)
          j=j+1
        endif
      enddo
      
      if( fileOutput ) then
        !open(10,file='tmpTTAnalyticwindow.'//fileDelta((len(trim(fileDelta))-15):(len(trim(fileDelta)))))
        open(10,file=datadirectory(1:len_trim(datadirectory))//'tmpgetAnalytic_window.dat')
        j=0        
        do i=1,WindowSIZE
          if( seismo(1,i) - startingTime .gt. 0.0) then
            if(j .eq. 0) then
              j=i
            endif
          endif
        enddo
        
        do i=1,xcorrlength
            write(10,*) seismo(1,1)+(j+i-2)*timedelta,seismoWindow(i),seismoRefWindow(i)
        enddo
        close(10)
      endif
      
      ! cross-correlation (by numerical recipes) 
      !call correl(seismoWindow,seismoRefWindow,WindowSIZE+zeropad,ans)
      length=size(seismo(2,:))
      !allocate(seismo1(length),seismo2(length),seismoPerturbed(length),stat=ierror)
      !if( ierror .ne. 0) call stopProgram( 'abort - getAnalyticTimelag at seismo1 allocation    ')
      seismo1(:)=0.0_WP
      seismo2(:)=0.0_WP
      seismoPerturbed(:)=0.0_WP
      
      ! number of file entries concerning our integral
      call readSeismo(seismo1,WindowSIZE,fileReference,startingTime,timedelta,entries,lastseismotime)
      if(beVerbose) print*,'entries to take:',entries
      seismo1(:)=0.0_WP 
      
      ! cut-off higher entries
      do i=1,length
        if( i .gt. entries) then
          seismoRef(2,i)=0.0_WP
          seismo(2,i)=0.0_WP
        endif
      enddo
      
      !--------------------------------------------------------------------------------------------------- 
      ! analytically derived expression for timelag
      ! get perturbed seismogram (only perturbations/differences between reference and delta seismogram)
      call getPerturbedSeismo(seismoRef(2,:),seismo(2,:),seismoPerturbed,length)
      
      ! for spline representation
      !allocate(X(length),Y(length),Q(3,length),F(3,length),stat=ierror)
      !if( ierror .ne. 0) then
      !  print*,'error allocating spline arrays'
      !  stop 'abort - getAnalyticTimelag'
      !endif

      ! get first time derivative of reference seismogram
      if( beVerbose ) print*,'    first deriv.. '      
      call getFirstDerivative(seismoRef(2,:),seismo1,timedelta,length)

      ! get second time derivative of reference seismogram
      if( beVerbose ) print*,'    second deriv.. '      
      call getSecondDerivative(seismoRef(2,:),seismo2,timedelta,length)

      ! calculate integrals
      if( beVerbose ) print*,'    first integral.. '            
      call getIntegral(seismo1,seismoPerturbed,integral1,timedelta,length)
      if( beVerbose ) print*,'    second integral.. '            
      call getIntegral(seismo2,seismoRef(2,:),integral2,timedelta,length)

      if(beVerbose) print*,'integral1:',integral1,'integral2:',integral2,'timestep:',timedelta
      t_lag=integral1/integral2
      
      end
      
      
!-----------------------------------------------------------------------
      subroutine getPerturbedSeismo(seismoRef,seismoDelta,seismoPerturbed,length)
!-----------------------------------------------------------------------
! calculates the difference between the reference and delta seismogram
! input:
!       seismoRef                             - reference seismogram 
!       seismoDelta                           - delta seismogram 
!       seismoPerturbed                   - perturbations
! returns: seismoPerturbed
      use precision         
      implicit none
      integer length,i
      real(WP):: seismoRef(length),seismoDelta(length),seismoPerturbed(length)
      
      length=size(seismoRef)      
      if( length .ne. size(seismoDelta)) then
        print*,' seismograms have different dimensions', length,size(seismoDelta)
        call stopProgram( 'abort - getPerturbedSeismo   ')
      endif
      
      do i=1,length
        seismoPerturbed(i)=seismoDelta(i)-seismoRef(i)
      enddo
      
      end
         
!-----------------------------------------------------------------------
      subroutine getFirstDerivative(seismo,seismoOut,timedelta,length)
!-----------------------------------------------------------------------
! calculates the first time derivative
! input:
!       seismo                             - reference seismogram 
!       seismo1                           - derivated seismogram 
!       timedelta                          - time step size
! returns: seismo1
      use precision;use splineFunction           
      implicit none
      integer:: length,i,ierror
      real(WP),intent(IN):: seismo(length)
      real(WP),intent(OUT):: seismoOut(length)
      real(WP):: timedelta
      double precision:: splineRepresentation,dfridr,err
      external:: dfridr,splineRepresentation
      
      ! check
      if( length .lt. 2 ) then
        print*,'derivative undefined!',length,timedelta
        call stopProgram('abort getFirstDerivative() - seismo too short!    ')      
      endif

      ! for spline representation
      allocate(X(length),Y(length),Q(3,length),F(3,length),stat=ierror)
      if( ierror .ne. 0) call stopProgram('getFirstDerivative() - error allocating spline arrays    ')
    
      !debug
      !if(DEBUG) print*,'first derivative:',length,size(seismo),size(seismoOut),size(X),size(Y)
        
      ! simplest finite difference
      !do i=1,length-1
      !    seismo1(i)=(seismo(i+1)-seismo(i))/timedelta
      !enddo
      !seismo1(length)=seismo1(length-1)
      
      ! get spline representation
      ! initialize
      i1=1
      i2=length
      do i=i1,i2    
        X(i)=dble((i-1)*timedelta)
        Y(i)=dble(seismo(i))
      enddo
      
      Q(:,:)=0.0d0
      F(:,:)=0.0d0
      call drspln(i1,i2,X,Y,Q,F)
      
      ! get time derivative
      do i=i1,i2
        if( WP .eq. 4 ) then        
          seismoOut(i)=real(dfridr(splineRepresentation,X(i),dble(2*timedelta),err))
        else
          seismoOut(i)=dfridr(splineRepresentation,X(i),dble(2*timedelta),err)
        endif
        if( err .gt. 0.1) then          
          print*,'derivative has big error:',i,seismo(i),seismoOut(i),err,i1,i2
        endif
      enddo      
      
      ! free memory
      deallocate(X,Y,Q,F)
      
      end
      

!-----------------------------------------------------------------------
      subroutine getSecondDerivative(seismo,seismoOut,timedelta,length)
!-----------------------------------------------------------------------
! calculates the first time derivative
! input:
!       seismo                             - reference seismogram 
!       seismo2                           - derivated seismogram 
!       timedelta                                      - time step size
! returns: seismo2
      use precision
      implicit none
      integer:: length
      real(WP),intent(IN):: seismo(length)
      real(WP),intent(OUT):: seismoOut(length)
      real(WP):: timedelta,seismo1(length)

      ! first time derivative
      call getFirstDerivative(seismo,seismo1,timedelta,length)
      ! second time
      call getFirstDerivative(seismo1,seismoOut,timedelta,length)      
      end

!-----------------------------------------------------------------------
      subroutine getIntegral(seismo1,seismo2,integral,timedelta,length)
!-----------------------------------------------------------------------
! calculates the integral of the product of two seismograms
! input:
!       seismo1,seismo2                     - reference seismograms
!       integral                                     - integral value
!       timedelta                                              - time step size
! returns: integral
      use precision;use nrtype; use nrutil; use splineFunction
      use nr            
      implicit none
      integer:: length,i,ierror
      real(WP),intent(IN):: seismo1(length),seismo2(length)
      real(WP):: timedelta,integral,val1,val2
      interface
        function spline_func(location)
        USE nrtype
        IMPLICIT NONE
        REAL(DP), DIMENSION(:),INTENT(IN) :: location
        REAL(DP), DIMENSION(size(location)):: spline_func
        END FUNCTION spline_func
      end interface
      
      ! initialize
      integral=0.0_WP
      
      !do i=1,length-1
      !  val1=seismo1(i)*seismo2(i)
      !  val2=seismo1(i+1)*seismo2(i+1)        
      !  integral=integral+val1*timedelta+(val2-val1)*timedelta*0.5_WP
      !enddo
      
      ! check with previous spline representation
      if( length .ne. size(Y) ) then
         print*,'different array lengths:',length,size(Y)
         call stopProgram( 'abort - getIntegral   ')
      endif

      ! for spline representation
      allocate(X(length),Y(length),Q(3,length),F(3,length),stat=ierror)
      if( ierror .ne. 0) call stopProgram('error allocating spline arrays    ')
            
      ! get spline representation
      ! initialize
      i1=1
      i2=length
      do i=i1,i2
        X(i)=dble((i-1)*timedelta)
        Y(i)=dble(seismo1(i)*seismo2(i))
      enddo
      Q(:,:)=0.0d0
      F(:,:)=0.0d0      
      call drspln(i1,i2,X,Y,Q,F)
      
      !get integral value
      if( WP .eq. 4 ) then
        integral=real(qsimp(spline_func, X(i1), X(i2)))
      else
        integral=qsimp(spline_func, X(i1), X(i2))
      endif
      !print*,'integral: ',i1,i2
      !print*,'    X:',X(i1),X(i2)
      !print*,'    Y:',Y(1),Y(2),Y(3)
      !print*,'    value=',integral
      !print*
      
      ! free memory
      deallocate(X,Y,Q,F)
      
      end

!-----------------------------------------------------------------------
	    double precision function splineRepresentation(location)
!-----------------------------------------------------------------------
      use nrtype; use splineFunction
      implicit none
      double precision, INTENT(IN) :: location
      double precision:: drsple
      external:: drsple
      
      ! calculate value for this representation
      splineRepresentation=drsple(i1,i2,X,Y,Q,location)      
      end function
            

!-----------------------------------------------------------------------   
      subroutine stopProgram(textline)
!-----------------------------------------------------------------------   
      use parallel
      implicit none
      character*128:: textline
      integer:: i,endindex,ierror
      logical:: flag
      
      ! console output
      endindex=index(textline,"  ")
      if( endindex .lt. 1 ) endindex = 128            
      print*,textline(1:endindex)
      
      ! stop MPI
      call MPI_Initialized(flag,ierror)
      if( flag == .true. .and. ierror == 0) then
        call MPI_Abort(MPI_COMM_WORLD, ierror )
        call MPI_FINALIZE(ierror)
      endif
      
      ! stop process
      stop 'abort - program '      
      end
      
