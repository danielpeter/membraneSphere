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
      subroutine determineFFTWindowsize(numofTimeSteps,WindowSIZE)
!-----------------------------------------------------------------------
! determines the length of the fft-seismograms (must be power of 2)
      !use propagationStartup; use verbosity; use filterType;
      use parallel, only: MASTER
      use verbosity, only: VERBOSE
      implicit none
      integer,intent(in):: numofTimeSteps
      integer,intent(out):: WindowSIZE
      real:: log2

      ! next higher number of power 2
      log2 = log(real(numofTimeSteps))/log(2.0)
      WindowSIZE = 2**(ceiling(log2))

      if (MASTER .and. VERBOSE) print *,'    fft window size =',WindowSIZE,'time steps=',numofTimeSteps
      end

!-----------------------------------------------------------------------
      subroutine determineFilterParameters()
!-----------------------------------------------------------------------
      use propagationStartup, only: cphasetype
      use filterType, only: bw_waveperiod,bw_frequency
      use precisions
      implicit none
      ! determine wave period
      if ( cphasetype == 'R35') bw_waveperiod=WAVEPERIOD_R35
      if ( cphasetype == 'R37') bw_waveperiod=WAVEPERIOD_R37
      if ( cphasetype == 'R40') bw_waveperiod=WAVEPERIOD_R40  ! wave period in seconds
      if ( cphasetype == 'R45') bw_waveperiod=WAVEPERIOD_R45
      if ( cphasetype == 'R50') bw_waveperiod=WAVEPERIOD_R50
      if ( cphasetype == 'R60') bw_waveperiod=WAVEPERIOD_R60
      if ( cphasetype == 'R75') bw_waveperiod=WAVEPERIOD_R75
      if ( cphasetype == 'R100') bw_waveperiod=WAVEPERIOD_R100
      if ( cphasetype == 'R150') bw_waveperiod=WAVEPERIOD_R150
      if ( cphasetype == 'R200') bw_waveperiod=WAVEPERIOD_R200
      if ( cphasetype == 'R250') bw_waveperiod=WAVEPERIOD_R250
      if ( cphasetype == 'R300') bw_waveperiod=WAVEPERIOD_R300

      if ( cphasetype == 'L35') bw_waveperiod=WAVEPERIOD_L35
      if ( cphasetype == 'L37') bw_waveperiod=WAVEPERIOD_L37
      if ( cphasetype == 'L40') bw_waveperiod=WAVEPERIOD_L40
      if ( cphasetype == 'L45') bw_waveperiod=WAVEPERIOD_L45
      if ( cphasetype == 'L50') bw_waveperiod=WAVEPERIOD_L50
      if ( cphasetype == 'L60') bw_waveperiod=WAVEPERIOD_L60
      if ( cphasetype == 'L75') bw_waveperiod=WAVEPERIOD_L75
      if ( cphasetype == 'L100') bw_waveperiod=WAVEPERIOD_L100
      if ( cphasetype == 'L150') bw_waveperiod=WAVEPERIOD_L150
      if ( cphasetype == 'L200') bw_waveperiod=WAVEPERIOD_L200
      if ( cphasetype == 'L250') bw_waveperiod=WAVEPERIOD_L250
      if ( cphasetype == 'L300') bw_waveperiod=WAVEPERIOD_L300

      ! get (upper) half-bandwidth frequency for filtering
      if ( BW_FIXFREQUENCY ) then
        bw_frequency = BW_HALFFREQUENCY
      else
        bw_frequency = BW_PERCENT/((BW_PERCENT+1.0)*bw_waveperiod)
      endif
      end


!-----------------------------------------------------------------------
      subroutine dofilterSeismogram(seismo,seismoLength)
!-----------------------------------------------------------------------
! filter given seismogram around the corner frequency
!
! returns: filtered seismo
      use verbosity; use filterType; use propagationStartup; use parallel
      implicit none
      integer:: seismoLength,entries,centerfrequencyindex,zeropad,i,xcorrlength
      real(WP),intent(INOUT):: seismo(2,seismoLength)
      real(WP):: startingTime,endtime,samplingFreq
      real(WP):: seismoTmp(2,WindowSIZE)

      !add zero for non-periodic data, influences precision
      !zeropad=WindowSIZE

      ! length of cross-correlation array
      xcorrlength = WindowSIZE !or: xcorrlength = WindowSIZE+zeropad

      ! attention of seismo arrays: zero padding of length WindowSIZE at the end to get corresponding length of cross-correlation
      !allocate(seismoTmp(2,xcorrlength),stat=ierror)
      !if ( ierror /= 0) call stopProgram( 'abort - dofilterSeismogram       ')

      ! determine startingTime
      endtime=seismo(1,seismoLength)
      call getStartTimeSeismogram(seismo,seismoLength,endtime,dt,startingTime)

      ! read in seismograms from startingTime
      seismoTmp(:,:)=0.0
      call captureSeismogram(seismo,seismoLength,seismoTmp,xcorrlength,startingTime,entries,endtime)

      ! console output
      if ( beVerbose ) then
        print *,'    captured seismogram:'
        print *,'      dt                  : ',dt
        print *,'      time first entry      : ',seismoTmp(1,1) !,seismoTmp(2,1)
        print *,'      time last entry      : ',seismoTmp(1,entries) !,seismoTmp(2,entries)
        print *,'      entries read        : ',entries
        print *,'      last time readin     : ',seismoTmp(1,entries) !,'seismo:',endtime
      endif

      ! apply hanning window  to smooth seismograms ends
      call taperSeismogram(seismoTmp,xcorrlength,entries,beVerbose)

      ! smallest sampled frequency
      samplingFreq=1.0_WP/(xcorrlength*dt)

      ! position index of bandwidth frequency for given waveperiod
      centerfrequencyIndex = nint( 1.0_WP/((bw_waveperiod + MEMBRANECORRECTION)*samplingFreq) )

      ! determine bandwidth integer depending on windowsize and requested bandwidth-frequency
      ! integer is a factor of 2 to have same number of frequencies filtered around the center frequency to the left and the right
      bw_width= 2*nint(bw_frequency*dt*xcorrlength)
      !bw_width= 2*abs(centerfrequencyIndex-nint(bw_frequency*dt*xcorrlength))

      ! console output
      if ( beVerbose ) then
        print *
        print *,'  filter parameter:'
        print *,'    sampling frequency    : ',samplingFreq
        print *,'    wave period           : ',bw_waveperiod+MEMBRANECORRECTION
        print *,'      initial period      : ',bw_waveperiod
        print *,'      membrane correction : ',MEMBRANECORRECTION
        print *,'    center frequency      : ',1.0_WP/(bw_waveperiod+MEMBRANECORRECTION)
        print *,'      index               : ',centerfrequencyIndex
        print *,'    bandwidth (full/half) : ',bw_width,bw_width/2
        print *,'        as frequency      : ',bw_width*samplingFreq,bw_width/2*samplingFreq
        if ( BUTTERWORTHFILTER) print *,'    power                 : ',BUTTERWORTH_POWER
        print *,'    windowsize            : ',xcorrlength
      endif

      ! check if filter frequency is correct
      if ( centerfrequencyIndex < 1 .or. centerfrequencyIndex > xcorrlength/2 ) then
        print *,'Error: filter is incorrect:',centerfrequencyIndex,xcorrlength
        call stopProgram('abort filterSeismogram() - incorrect centerfrequencyIndex    ')
      endif

      ! debug output
      if ( fileOutput .and. MASTER .and. beVerbose) then
        print *,'printing to file:',datadirectory(1:len_trim(datadirectory))//'Filter_input.dat'
        open(10,file=datadirectory(1:len_trim(datadirectory))//'Filter_input.dat')
        do i = 1,xcorrlength
          write(10,*) seismoTmp(1,i),seismoTmp(2,i)
        enddo
        close(10)
      endif

      ! filter seismograms in frequency domain
      call filterSeismogram(seismoTmp(2,:),xcorrlength,bw_width,BUTTERWORTH_POWER,centerfrequencyIndex)

      ! avoid taking temporary-working data in second half of fourier-transformed seismogram
      ! return as new filtered seismogram
      do i = 1,seismoLength
        seismo(2,i)=seismoTmp(2,i)
      enddo

      ! apply hanning window  to smooth seismograms ends again
      !call taperSeismogram(seismo,seismoLength,seismoLength,beVerbose)

      ! debug output
      if ( fileOutput .and. MASTER .and. beVerbose ) then
        print *,'printing to file:',datadirectory(1:len_trim(datadirectory))//'Filter_output.dat'
        open(10,file=datadirectory(1:len_trim(datadirectory))//'Filter_output.dat')
        do i = 1,seismoLength
          write(10,*) seismo(1,i),seismo(2,i)
        enddo
        close(10)
      endif

      ! be only once verbose
      beVerbose = .false.
      end


!-----------------------------------------------------------------------
      subroutine filterSeismos(data1,data2,n,bandwidth,power,centerIndex)
!-----------------------------------------------------------------------
! filters two input seismograms with a butterworth filter in the frequency domain
! around the L150s frequency and a bandwidth of 10*dfreq
!
! input:
!     data1,data2   - seismograms to filter
!     n                   - size of seismograms (must be power of 2)
!     bandwidth, power, centerIndex
!                           - parameters bandwidth, power and bandwidth frequency index for butterworth filtering
!
! returns: filtered data1, filtered data2
      use verbosity; use filterType
      USE nrtype; USE nrutil
      USE nr, only: four1
      implicit none
      integer,intent(in):: n,bandwidth,power,centerIndex
      real(WP),dimension(n),intent(inout):: data1,data2
      integer:: i
      real(WP),dimension(n):: butterworth
      double complex:: fftdata1(n),fftdata2(n),dn2,fftdataTest(n)

      ! call for each seismogram
      call filterSeismogram(data1,n,bandwidth,power,centerIndex)
      call filterSeismogram(data2,n,bandwidth,power,centerIndex)
      return
      end subroutine


!-----------------------------------------------------------------------
      subroutine filterSeismogram(seismodata,dataLength,bandwidth,power,centerIndex)
!-----------------------------------------------------------------------
! filters input seismogram (e.g. with a butterworth filter in the frequency domain
! around the L150s frequency and a bandwidth of 10*dfreq)
!
! input:
!     seismodata            - seismogram to filter
!     n                            - size of seismogram (must be power of 2)
!     bandwidth, power, centerIndex
!                           - parameters bandwidth, power and bandwidth frequency index for butterworth filtering
!
! returns: filtered seismodata
      use verbosity; use filterType; use nrtype; use nrutil;
      use nr, only: four1; use parallel; use propagationStartup
      implicit none
      integer,intent(in):: dataLength,bandwidth,power,centerIndex
      real(WP),dimension(dataLength),intent(inout):: seismodata
      integer:: i,firstEntry
      double complex::fftdata(dataLength),fftdataTest(dataLength),dn2,scale
      double precision::butterworth
      logical,parameter:: TEST = .false.

      !Normalization for inverse FFT
      dn2 = 1.0d0/dataLength

      ! construct complex array for fourier-transformation
      fftdata(:)=dcmplx(seismodata(:),0.0d0)

      ! test fourier transformation
      if ( TEST .and. fileOutput .and. MASTER ) then
        fftdataTest(:)=fftdata(:)

        ! debug output
        print *,'  printing to file:',datadirectory(1:len_trim(datadirectory))//'FilterTest_before.dat'
        open(10,file=datadirectory(1:len_trim(datadirectory))//'FilterTest_before.dat')
        do i = 1,dataLength
          write(10,*) i,real(fftdataTest(i)),aimag(fftdataTest(i))
        enddo
        close(10)
        ! fourier-transform
        call four1(fftdataTest,1)
        !debug output
        print *,'  printing to file:',datadirectory(1:len_trim(datadirectory))//'FilterTest_FFT.dat'
        open(10,file=datadirectory(1:len_trim(datadirectory))//'FilterTest_FFT.dat')
        do i = 1,dataLength
          write(10,*) i,real(fftdataTest(i)),aimag(fftdataTest(i))
        enddo
        close(10)
        ! normalize
        do i = 1,dataLength
          fftdataTest(i)=fftdataTest(i)*dn2  ! normalization
        enddo
        ! inverse fourier-transform
        call four1(fftdataTest,-1)
        ! debug output
        print *,'  printing to file:',datadirectory(1:len_trim(datadirectory))//'FilterTest_after.dat'
        open(10,file=datadirectory(1:len_trim(datadirectory))//'FilterTest_after.dat')
        do i = 1,dataLength
          write(10,*) i,real(fftdataTest(i)),aimag(fftdataTest(i))
        enddo
        close(10)
      endif

      ! fourier transform
      call four1(fftdata,1)

      ! debug output
      if ( TEST .and. fileOutput .and. MASTER ) then
        print *,'  printing to file:',datadirectory(1:len_trim(datadirectory))//'FilterSeismo_data.dat'
        open(10,file=datadirectory(1:len_trim(datadirectory))//'FilterSeismo_data.dat')
        do i = 1,dataLength
          write(10,*) i,real(fftdata(i)),aimag(fftdata(i))
        enddo
        close(10)
      endif

      ! apply filter
      if ( beVerbose) then
        if ( BUTTERWORTHFILTER) then
          print *,'    butterworth bandpass filter'
        else
          print *,'    bandpass filter'
        endif
      endif

      ! debug output
      if ( TEST .and. fileOutput .and. MASTER ) then
        print *,'  printing to file:',datadirectory(1:len_trim(datadirectory))//'FilterSeismo_datafiltered.dat'
        open(10,file=datadirectory(1:len_trim(datadirectory))//'FilterSeismo_datafiltered.dat')
        write(10,*) "# filterSeismogram"
        write(10,*) "# index realData imagData filter"
      endif

      do i = 1,dataLength
        ! butterworth bandpass filter
        if ( BUTTERWORTHFILTER) then
          butterworth = 1.0d0 - 1.0d0/(1.0d0+(dble(i*bandwidth)/((dble(i))**2-(dble(centerIndex))**2))**(2*power))
        else
          ! bandpass filter
          ! attention to storage of data in fftdata:
          ! fftdata(1) is only a real value representing a constant function ( 0 frequency), then fftdata(2) is first frequency part;
          ! ( see also http://www.parasitaere-kapazitaeten.net/Pd/fft_und_pd.htm )
          ! first "half" (2:n/2) are positive frequencies in increasing occurrence , (n:n/2) second half negative frequencies in decreasing manner
          ! filtering by bandpass has to filter both positive and negative frequencies
          ! ( see also http://www.library.cornell.edu/nr/bookcpdf/c12-2.pdf )
          if ( abs((centerIndex+1)-i) <= bandwidth/2 .or. abs(datalength-(centerIndex-1)-i) <= bandwidth/2) then
            butterworth = 1.0d0
          else
            butterworth = 0.0d0
          endif
        endif

        ! filter data
        fftdata(i)=butterworth*fftdata(i)

        ! use same amplitude for all frequencies which were considered
        !scale=fftdata(centerIndex)
        !fftdata(i)=dcmplx(butterworth*real(fftdata(i))/real(fftdata(i))*real(scale),butterworth*aimag(fftdata(i))/aimag(fftdata(i))*aimag(scale))

        ! debug output
        if ( TEST .and. fileOutput .and. MASTER ) then
          write(10,*) i,real(fftdata(i)),aimag(fftdata(i)),butterworth
        endif

        ! prepare for inverse FFT
        fftdata(i)=fftdata(i)*dn2  !with Normalization for inverse FFT
      enddo

      ! debug output
      if ( TEST .and. fileOutput .and. MASTER ) then
        close(10)
      endif

      ! applies inverse fourier transformation
      call four1(fftdata,-1)

      ! returns trace, as real part of the filtered data array
      if ( WP == 4 ) then
        seismodata(:)=real(fftdata(:))
      else
        seismodata(:)=dreal(fftdata(:))
      endif
      end subroutine


!-----------------------------------------------------------------------
      subroutine taperSeismos(seismo,seismoRef,seismoLengthRef,entries)
!-----------------------------------------------------------------------
! taper seismogram begining and end
!
! input:
!   seismo,seismoRef - seismograms to taper
!   seismoLengthRef   - length of seismograms
!  entries                     - non-zero entries of seismograms
!
! returns: tapered seismo & seismoRef
      use verbosity; use precisions
      implicit none
      integer::seismoLengthRef,entries,ileft,irigth,hannwindow
      real(WP)::seismo(2,seismoLengthRef),seismoRef(2,seismoLengthRef)
      real(WP)::hannA,hannfactor

      ! taper the two seismograms
      call taperSeismogram(seismo,seismoLengthRef,entries,.false.)
      call taperSeismogram(seismoRef,seismoLengthRef,entries,beVerbose)

      end

!-----------------------------------------------------------------------
      subroutine taperSeismogram(seismo,seismoLength,entries,verboseOutput)
!-----------------------------------------------------------------------
! taper seismogram begining and end
!
! input:
!   seismo                    - seismogram to taper
!   seismoLength        - length of seismogram
!  entries                     - non-zero entries of seismograms
!  verbose                    - be verbose
! returns: tapered seismo
      use filterType
      implicit none
      integer,intent(in):: seismoLength,entries
      real(WP),intent(inout):: seismo(2,seismoLength)
      logical,intent(in):: verboseOutput
      integer:: ileft,iright,hannwindow,i
      real(WP):: hannA,hannfactor

      ! apply hanning window  to smooth seismograms ends
      ! (info: http://www.teemach.com/FFTProp/FFTProperties/FFTProperties.htm )
      ileft = 0
      iright = entries+1
      hannwindow=int(entries*HANNWINDOW_PERCENT)  ! width is x% of actual entries of seismogram
      hannA = PI/hannwindow
      do i = 1,hannwindow
        hannfactor = 0.5*(1-cos(hannA*i))
        seismo(2,ileft+i)=seismo(2,ileft+i)*hannfactor
        seismo(2,iright-i)=seismo(2,iright-i)*hannfactor
      enddo
      if ( verboseOutput ) then
        print *,'    hanning window applied:',ileft+1,iright-1,'window width:',hannwindow
      endif
      end


!-----------------------------------------------------------------------
      subroutine captureSeismos(seismo1,seismo2,seismo,seismoRef, &
                  seismoLength,seismoLengthRef,startingTime,entries,endtime)
!-----------------------------------------------------------------------
! read in seismo and seismoRef with corresponding startingTime and window length
!
! input:
!  seismo1,seismo2  - given seismograms
!  seismoLength        - length of seismo1 and seismo2
!  startingTime          - read in from this time on
!
! returns: filled seismo & seismoRef, entries, endtime
      use precisions
      implicit none
      integer:: seismoLength,seismoLengthRef,entries
      real(WP):: startingTime,endtime
      real(WP),intent(in):: seismo1(2,seismoLength),seismo2(2,seismoLength)
      real(WP),intent(inout):: seismo(2,seismoLengthRef),seismoRef(2,seismoLengthRef)

      ! capture first seismogram
      call captureSeismogram(seismo1,seismoLength,seismo,seismoLengthRef,startingTime,entries,endtime)
      ! capture second seismogram
      call captureSeismogram(seismo2,seismoLength,seismoRef,seismoLengthRef,startingTime,entries,endtime)

      end

!-----------------------------------------------------------------------
      subroutine captureSeismogram(seismoData,seismoDataLength,seismoWindow, &
                            seismoWindowLength,startingTime,entries,endtime)
!-----------------------------------------------------------------------
! read in seismoData with corresponding startingTime and window length
!
! input:
!  seismoData                   - given seismogram
!  seismoDataLength        - length of seismo
!  startingTime          - read in from this time on
!
! returns: filled seismoWindow, entries, endtime
      use parallel; use filterType
      implicit none
      integer:: index,i,seismoDataLength,seismoWindowLength
      real(WP):: time,displace
      real(WP),intent(in):: seismoData(2,seismoDataLength),startingTime
      real(WP),intent(inout):: seismoWindow(2,seismoWindowLength)
      integer,intent(out):: entries
      real(WP),intent(out):: endtime

      index = 0
      do i = 1,seismoDataLength
        if (index < seismoWindowLength) then
          time=seismoData(1,i)
          displace=seismoData(2,i)

          if ( time - startingTime > 0.0) then ! starts with times from 'start time' on
            index = index+1
            seismoWindow(1,index) = time
            seismoWindow(2,index) = displace
            endtime = time
          endif
        else
          !read until end is reached
          time=seismoData(1,i) !,sourceterm
          endtime = time
        endif
      enddo
      entries = index

      end subroutine

