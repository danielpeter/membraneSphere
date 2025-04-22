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
  use parallel, only: MAIN_PROCESS
  use verbosity, only: VERBOSE
  implicit none
  integer,intent(in):: numofTimeSteps
  integer,intent(out):: WindowSIZE
  ! local parameters
  real:: log2

  ! next higher number of power 2
  log2 = log(real(numofTimeSteps))/log(2.0)
  WindowSIZE = 2**(ceiling(log2))

  ! debug - double window size to check influence on FFT
  !WindowSize = 2 * WindowSize

  if (MAIN_PROCESS .and. VERBOSE) then
    print *,'    FFT:'
    print *,'      window size         : ',WindowSIZE
    print *,'      number of time steps: ',numofTimeSteps
  endif

  end subroutine

!-----------------------------------------------------------------------
  subroutine determineFilterParameters()
!-----------------------------------------------------------------------
  use propagationStartup, only: cphasetype
  use filterType, only: bw_waveperiod,bw_frequency
  use precisions
  implicit none

  ! check if anything to do
  if (trim(cphasetype) == "") return

  ! determine wave period
  select case (trim(cphasetype))
  ! Rayleigh waves
  case ('R35')
    bw_waveperiod = WAVEPERIOD_R35
  case ('R37')
    bw_waveperiod = WAVEPERIOD_R37
  case ('R40')
    bw_waveperiod = WAVEPERIOD_R40  ! wave period in seconds
  case ('R45')
    bw_waveperiod = WAVEPERIOD_R45
  case ('R50')
    bw_waveperiod = WAVEPERIOD_R50
  case ('R60')
    bw_waveperiod = WAVEPERIOD_R60
  case ('R75')
    bw_waveperiod = WAVEPERIOD_R75
  case ('R100')
    bw_waveperiod = WAVEPERIOD_R100
  case ('R150')
    bw_waveperiod = WAVEPERIOD_R150
  case ('R200')
    bw_waveperiod = WAVEPERIOD_R200
  case ('R250')
    bw_waveperiod = WAVEPERIOD_R250
  case ('R300')
    bw_waveperiod = WAVEPERIOD_R300

  ! Love waves
  case ('L35')
    bw_waveperiod = WAVEPERIOD_L35
  case ('L37')
    bw_waveperiod = WAVEPERIOD_L37
  case ('L40')
    bw_waveperiod = WAVEPERIOD_L40
  case ('L45')
    bw_waveperiod = WAVEPERIOD_L45
  case ('L50')
    bw_waveperiod = WAVEPERIOD_L50
  case ('L60')
    bw_waveperiod = WAVEPERIOD_L60
  case ('L75')
    bw_waveperiod = WAVEPERIOD_L75
  case ('L100')
    bw_waveperiod = WAVEPERIOD_L100
  case ('L150')
    bw_waveperiod = WAVEPERIOD_L150
  case ('L200')
    bw_waveperiod = WAVEPERIOD_L200
  case ('L250')
    bw_waveperiod = WAVEPERIOD_L250
  case ('L300')
    bw_waveperiod = WAVEPERIOD_L300

  ! default
  case default
    print *,'Error: phase type not recognized: ',trim(cphasetype)
    call stopProgram('Invalid phase type in determineFilterParameters()    ')
  end select

  ! get (upper) half-bandwidth frequency for filtering
  if (BW_FIXFREQUENCY) then
    bw_frequency = BW_HALFFREQUENCY
  else
    bw_frequency = BW_PERCENT/((BW_PERCENT+1.0)*bw_waveperiod)
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine dofilterSeismogram(seismo,seismoLength)
!-----------------------------------------------------------------------
! filter given seismogram around the corner frequency
!
! returns: filtered seismo
  use verbosity; use filterType; use propagationStartup; use parallel
  use precisions
  implicit none
  integer,intent(in):: seismoLength
  real(WP),intent(inout):: seismo(2,seismoLength)
  ! local parameters
  integer:: entries,centerfrequencyindex,i,xcorrlength
  real(WP):: startingTime,endtime,samplingFreq
  real(WP):: seismoTmp(2,WindowSIZE)
  real(WP):: traceTmp(WindowSIZE)

  !add zero for non-periodic data, influences precision
  !zeropad=WindowSIZE

  ! length of cross-correlation array
  xcorrlength = WindowSIZE !or: xcorrlength = WindowSIZE+zeropad

  ! attention of seismo arrays: zero padding of length WindowSIZE at the end to get corresponding length of cross-correlation
  !allocate(seismoTmp(2,xcorrlength),stat=ier)
  !if (ier /= 0) call stopProgram( 'abort - dofilterSeismogram       ')

  ! determine startingTime
  endtime = seismo(1,seismoLength)
  call getStartTimeSeismogram(seismo,seismoLength,endtime,dt,startingTime)

  ! read in seismograms from startingTime
  seismoTmp(:,:) = 0.0
  call captureSeismogram(seismo,seismoLength,seismoTmp,xcorrlength,startingTime,entries,endtime)

  ! console output
  if (beVerbose) then
    print *,'    captured seismogram:'
    print *,'      dt                  : ',dt
    print *,'      time first entry    : ',seismoTmp(1,1) !,seismoTmp(2,1)
    print *,'      time last entry     : ',seismoTmp(1,entries) !,seismoTmp(2,entries)
    print *,'      entries read        : ',entries
    print *,'      last time readin    : ',seismoTmp(1,entries) !,'seismo:',endtime
  endif

  ! apply hanning window  to smooth seismograms ends
  call taperSeismogram(seismoTmp,xcorrlength,entries,beVerbose)

  ! smallest sampled frequency
  samplingFreq = 1.0_WP/(xcorrlength*dt)

  ! position index of bandwidth frequency for given waveperiod
  centerfrequencyIndex = nint( 1.0_WP/((bw_waveperiod + MEMBRANECORRECTION)*samplingFreq) )

  ! determine bandwidth integer depending on windowsize and requested bandwidth-frequency
  ! integer is a factor of 2 to have same number of frequencies filtered around the center frequency to the left and the right
  bw_width = 2*nint(bw_frequency*dt*xcorrlength)
  !bw_width = 2*abs(centerfrequencyIndex-nint(bw_frequency*dt*xcorrlength))

  ! console output
  if (beVerbose) then
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
    if (BUTTERWORTHFILTER) print *,'    power                 : ',BUTTERWORTH_POWER
    print *,'    windowsize            : ',xcorrlength
  endif

  ! check if filter frequency is correct
  if (centerfrequencyIndex < 1 .or. centerfrequencyIndex > xcorrlength/2) then
    print *,'Error: filter is incorrect:',centerfrequencyIndex,xcorrlength
    call stopProgram('abort filterSeismogram() - incorrect centerfrequencyIndex    ')
  endif

  ! debug output
  if (fileOutput .and. MAIN_PROCESS .and. beVerbose) then
    print *,'printing to file: ',trim(datadirectory)//'Filter_input.dat'
    open(IOUT,file=trim(datadirectory)//'Filter_input.dat')
    do i = 1,xcorrlength
      write(IOUT,*) seismoTmp(1,i),seismoTmp(2,i)
    enddo
    close(IOUT)
  endif

  ! seismo trace array
  traceTmp(:) = seismoTmp(2,:)

  ! filter seismograms in frequency domain
  call filterSeismogram(traceTmp,xcorrlength,bw_width,BUTTERWORTH_POWER,centerfrequencyIndex)

  ! avoid taking temporary-working data in second half of fourier-transformed seismogram
  ! return as new filtered seismogram
  do i = 1,seismoLength
    seismo(2,i) = traceTmp(i)
  enddo

  ! apply hanning window  to smooth seismograms ends again
  !call taperSeismogram(seismo,seismoLength,seismoLength,beVerbose)

  ! debug output
  if (fileOutput .and. MAIN_PROCESS .and. beVerbose) then
    print *,'printing to file: ',trim(datadirectory)//'Filter_output.dat'
    open(IOUT,file=trim(datadirectory)//'Filter_output.dat')
    do i = 1,seismoLength
      write(IOUT,*) seismo(1,i),seismo(2,i)
    enddo
    close(IOUT)
  endif

  ! be only once verbose
  beVerbose = .false.

  end subroutine


!-----------------------------------------------------------------------
  subroutine filterSeismos(inputdata1,inputdata2,n,bandwidth,power,centerIndex)
!-----------------------------------------------------------------------
! filters two input seismograms with a butterworth filter in the frequency domain
! around the L150s frequency and a bandwidth of 10*dfreq
!
! input:
!     inputdata1,inputdata2   - seismograms to filter
!     n                   - size of seismograms (must be power of 2)
!     bandwidth, power, centerIndex
!                           - parameters bandwidth, power and bandwidth frequency index for butterworth filtering
!
! returns: filtered data1, filtered data2
  use verbosity; use filterType
  implicit none
  integer,intent(in):: n,bandwidth,power,centerIndex
  real(WP),dimension(n),intent(inout):: inputdata1,inputdata2

  ! call for each seismogram
  call filterSeismogram(inputdata1,n,bandwidth,power,centerIndex)
  call filterSeismogram(inputdata2,n,bandwidth,power,centerIndex)

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
  use verbosity; use filterType
  use parallel; use propagationStartup; use precisions
  implicit none
  integer,intent(in):: dataLength,bandwidth,power,centerIndex
  real(WP),dimension(dataLength),intent(inout):: seismodata
  ! local parameters
  integer:: i,N
  complex(kind=8):: fftdata(dataLength),fftdataTest(dataLength) !,scale
  double precision:: dnorm
  double precision:: butterworth
  interface
    subroutine FFT_complex(data,N,isign)
      implicit none
      integer, intent(in) :: N,isign
      complex(kind=8), dimension(N), intent(inout) :: data
    end subroutine
  end interface
  logical,parameter:: DEBUG_TEST = .false.

  !Normalization for inverse FFT
  N = dataLength
  dnorm = 1.0d0 / dble(N)

  ! construct complex array for fourier-transformation
  fftdata(:) = cmplx(seismodata(:),0.0d0,kind=8)

  ! test fourier transformation
  if (DEBUG_TEST .and. fileOutput .and. MAIN_PROCESS) then
    print *,'  testing filter seismogram:'
    print *,'    size N = ',N

    fftdataTest(:) = fftdata(:)

    ! debug output
    print *,'    printing to file: ',trim(datadirectory)//'FilterTest_before.dat'
    open(IOUT,file=trim(datadirectory)//'FilterTest_before.dat')
    write(IOUT,*) "# format: #i #real-part #imaginary-part"
    do i = 1,dataLength
      write(IOUT,*) i,real(fftdataTest(i)),aimag(fftdataTest(i))
    enddo
    close(IOUT)

    ! fourier-transform
    call FFT_complex(fftdataTest,N,1)

    !debug output
    print *,'    printing to file: ',trim(datadirectory)//'FilterTest_FFT.dat'
    open(IOUT,file=trim(datadirectory)//'FilterTest_FFT.dat')
    write(IOUT,*) "# format: #i #real-part #imaginary-part"
    do i = 1,dataLength
      write(IOUT,*) i,real(fftdataTest(i)),aimag(fftdataTest(i))
    enddo
    close(IOUT)

    ! normalize
    do i = 1,dataLength
      fftdataTest(i) = fftdataTest(i) * dnorm  ! normalization
    enddo

    ! inverse fourier-transform
    call FFT_complex(fftdataTest,N,-1)

    ! debug output
    print *,'    printing to file: ',trim(datadirectory)//'FilterTest_after.dat'
    open(IOUT,file=trim(datadirectory)//'FilterTest_after.dat')
    write(IOUT,*) "# format: #i #real-part #imaginary-part"
    do i = 1,dataLength
      write(IOUT,*) i,real(fftdataTest(i)),aimag(fftdataTest(i))
    enddo
    close(IOUT)
  endif

  ! fourier transform
  call FFT_complex(fftdata,N,1)

  ! debug output
  if (DEBUG_TEST .and. fileOutput .and. MAIN_PROCESS) then
    print *,'    printing to file: ',trim(datadirectory)//'FilterSeismo_data.dat'
    open(IOUT,file=trim(datadirectory)//'FilterSeismo_data.dat')
    write(IOUT,*) "# format: #i #real-part #imaginary-part"
    do i = 1,dataLength
      write(IOUT,*) i,real(fftdata(i)),aimag(fftdata(i))
    enddo
    close(IOUT)
  endif

  ! apply filter
  if (beVerbose) then
    if (BUTTERWORTHFILTER) then
      print *,'    butterworth bandpass filter'
    else
      print *,'    bandpass filter'
    endif
  endif

  ! debug output
  if (DEBUG_TEST .and. fileOutput .and. MAIN_PROCESS) then
    print *,'    printing to file: ',trim(datadirectory)//'FilterSeismo_datafiltered.dat'
    open(IOUT,file=trim(datadirectory)//'FilterSeismo_datafiltered.dat')
    write(IOUT,*) "# filterSeismogram"
    write(IOUT,*) "# format: index realData imagData filter"
  endif

  do i = 1,dataLength
    ! butterworth bandpass filter
    if (BUTTERWORTHFILTER) then
      butterworth = 1.0d0 - 1.0d0/(1.0d0+(dble(i*bandwidth)/((dble(i))**2-(dble(centerIndex))**2))**(2*power))
    else
      ! bandpass filter
      ! attention to storage of data in fftdata:
      ! fftdata(1) is only a real value representing a constant function ( 0 frequency), then fftdata(2) is first frequency part;
      ! ( see also http://www.parasitaere-kapazitaeten.net/Pd/fft_und_pd.htm )
      ! first "half" (2:n/2) are positive frequencies in increasing occurrence , (n:n/2) second half negative frequencies in decreasing manner
      ! filtering by bandpass has to filter both positive and negative frequencies
      ! ( see also http://www.library.cornell.edu/nr/bookcpdf/c12-2.pdf )
      if (abs((centerIndex+1)-i) <= bandwidth/2 .or. abs(datalength-(centerIndex-1)-i) <= bandwidth/2) then
        butterworth = 1.0d0
      else
        butterworth = 0.0d0
      endif
    endif

    ! filter data
    fftdata(i) = butterworth*fftdata(i)

    ! use same amplitude for all frequencies which were considered
    !scale=fftdata(centerIndex)
    !fftdata(i)=cmplx(butterworth*real(fftdata(i))/real(fftdata(i))*real(scale), &
    !                 butterworth*aimag(fftdata(i))/aimag(fftdata(i))*aimag(scale),kind=8)

    ! debug output
    if (DEBUG_TEST .and. fileOutput .and. MAIN_PROCESS) then
      write(IOUT,*) i,real(fftdata(i)),aimag(fftdata(i)),butterworth
    endif

    ! prepare for inverse FFT
    fftdata(i) = fftdata(i) * dnorm  !with Normalization for inverse FFT
  enddo

  ! debug output
  if (DEBUG_TEST .and. fileOutput .and. MAIN_PROCESS) then
    close(IOUT)
    print *,'    test done'
  endif

  ! applies inverse fourier transformation
  call FFT_complex(fftdata,N,-1)

  ! returns trace, as real part of the filtered data array
  seismodata(:) = real(fftdata(:),kind=WP)

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
  integer,intent(in):: seismoLengthRef,entries
  real(WP),intent(inout):: seismo(2,seismoLengthRef),seismoRef(2,seismoLengthRef)

  ! taper the two seismograms
  call taperSeismogram(seismo,seismoLengthRef,entries,.false.)
  call taperSeismogram(seismoRef,seismoLengthRef,entries,beVerbose)

  end subroutine

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
  use precisions; use filterType
  implicit none
  integer,intent(in):: seismoLength,entries
  real(WP),intent(inout):: seismo(2,seismoLength)
  logical,intent(in):: verboseOutput
  ! local parameters
  integer:: ileft,iright,hannwindow,i
  real(WP):: hannA,hannfactor

  ! apply hanning window  to smooth seismograms ends
  ! (info: http://www.teemach.com/FFTProp/FFTProperties/FFTProperties.htm )
  ileft = 0
  iright = entries+1

  hannwindow = int(entries*HANNWINDOW_PERCENT)  ! width is x% of actual entries of seismogram
  hannA = PI/hannwindow

  do i = 1,hannwindow
    hannfactor = 0.5*(1-cos(hannA*i))
    seismo(2,ileft+i) = seismo(2,ileft+i)*hannfactor
    seismo(2,iright-i) = seismo(2,iright-i)*hannfactor
  enddo

  if (verboseOutput) then
    print *,'    hanning window applied:',ileft+1,iright-1,'window width:',hannwindow
  endif

  end subroutine


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
  integer,intent(in):: seismoLength,seismoLengthRef
  real(WP),intent(in):: seismo1(2,seismoLength),seismo2(2,seismoLength)
  real(WP),intent(inout):: seismo(2,seismoLengthRef),seismoRef(2,seismoLengthRef)
  real(WP),intent(in):: startingTime
  integer,intent(out):: entries
  real(WP),intent(out):: endtime

  ! capture first seismogram
  call captureSeismogram(seismo1,seismoLength,seismo,seismoLengthRef,startingTime,entries,endtime)
  ! capture second seismogram
  call captureSeismogram(seismo2,seismoLength,seismoRef,seismoLengthRef,startingTime,entries,endtime)

  end subroutine

!-----------------------------------------------------------------------
  subroutine captureSeismogram(seismoData,seismoDataLength,seismoWindow,seismoWindowLength, &
                               startingTime,entries,endtime)
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
  integer,intent(in):: seismoDataLength,seismoWindowLength
  real(WP),intent(in):: seismoData(2,seismoDataLength)
  real(WP),intent(inout):: seismoWindow(2,seismoWindowLength)
  real(WP),intent(in):: startingTime
  integer,intent(out):: entries
  real(WP),intent(out):: endtime
  ! local parameters
  integer:: index,i
  real(WP):: time,displace

  index = 0
  do i = 1,seismoDataLength
    if (index < seismoWindowLength) then
      time = seismoData(1,i)
      displace = seismoData(2,i)

      if (time - startingTime > 0.0) then ! starts with times from 'start time' on
        index = index+1
        seismoWindow(1,index) = time
        seismoWindow(2,index) = displace
        endtime = time
      endif
    else
      !read until end is reached
      time = seismoData(1,i) !,sourceterm
      endtime = time
    endif
  enddo
  entries = index

  end subroutine


!-----------------------------------------------------------------------
  subroutine FFT_complex(data,N,isign)
!-----------------------------------------------------------------------
! Performs a Fast Fourier Transform (FFT) on a complex data array.
! The input array size must be a power of 2.
!
! Arguments:
!   data      - Complex input/output array containing the signal.
!   isign     - Determines FFT direction: +1 for forward FFT, -1 for inverse FFT.
  implicit none
  integer, intent(in) :: N,isign
  complex(kind=8), dimension(N), intent(inout) :: data
  ! local parameters
  complex(kind=8), dimension(:,:), allocatable :: matrix_data,temp
  complex(kind=8), dimension(:), allocatable :: w_factors,w_prime
  double precision, dimension(:), allocatable :: theta
  integer :: nrow,ncol,j,ier
  double precision, parameter :: TWO_PI = 2.d0 * 3.1415926535897931d0

  ! Ensure data_size is a power of 2 (necessary for FFT)
  if (iand(N,N-1) /= 0) call stopProgram('FFT_complex: N must be a power of 2    ')

  ! Compute row and column sizes for 2D representation
  nrow = 2**ceiling(0.5d0 * log(real(N))/0.693147d0)
  ncol = N/nrow

  ! Allocate memory for matrices and auxiliary arrays
  allocate(matrix_data(nrow,ncol), &
           theta(nrow), &
           w_factors(nrow), &
           w_prime(nrow), &
           temp(ncol,nrow),stat=ier)
  if (ier /= 0) stop 'Error allocating matrix data for fft'

  ! Reshape 1D input array into a 2D matrix format
  matrix_data = reshape(data,shape(matrix_data))

  ! Perform FFT along matrix rows
  call myfourrow(matrix_data,isign)

  ! Compute phase angles for twiddle factors
  theta = arithmetic_sequence(0.d0,dble(isign),nrow) * TWO_PI / N

  ! Compute FFT twiddle factors
  w_prime = cmplx(-2.0d0*sin(0.5d0*theta)**2,sin(theta),kind=8)
  w_factors = cmplx(1.0d0,0.0d0)

  ! Apply twiddle factors to each column
  do j = 2,ncol
    w_factors = w_factors * w_prime + w_factors
    matrix_data(:,j) = matrix_data(:,j) * w_factors
  enddo

  ! Transpose the matrix for column-wise FFT processing
  temp = transpose(matrix_data)

  ! Perform FFT along matrix columns
  call myfourrow(temp,isign)

  ! Reshape back to 1D array format
  data = reshape(temp,shape(data))

  ! free arrays
  deallocate(matrix_data,w_factors,w_prime,theta,temp)

contains

  !-----------------------------------------------------------------------
  function arithmetic_sequence(first_term, step_size, num_terms)
  !-----------------------------------------------------------------------
  ! Generates an arithmetic sequence with a given first term, step size,
  ! and number of terms.
  !
  ! Arguments:
  !   first_term - The first term in the sequence.
  !   step_size  - The common difference between consecutive terms.
  !   num_terms  - The number of terms in the sequence.
  !
  ! returns:
  !   arithmetic_sequence - An array containing the arithmetic sequence.
    implicit none
    double precision, intent(in) :: first_term, step_size          ! Input: First term and step size
    integer, intent(in) :: num_terms                               ! Input: Number of terms in sequence
    double precision, dimension(num_terms) :: arithmetic_sequence  ! Output: Generated sequence
    ! local parameters
    integer :: index

    ! Initialize first term if sequence has elements
    if (num_terms > 0) arithmetic_sequence(1) = first_term

    ! use simple loop
    do index = 2, num_terms
      arithmetic_sequence(index) = arithmetic_sequence(index - 1) + step_size
    enddo

    return
  end function


  !-----------------------------------------------------------------------
  subroutine myfourrow(data, isign)
  !-----------------------------------------------------------------------
  ! Performs row-wise Fast Fourier Transform (FFT) on a 2D complex array.
  ! The input array size along the second dimension must be a power of 2.
  !
  ! Arguments:
  !   data      - Complex 2D array (modified in place with FFT results).
  !   isign     - Determines FFT direction: +1 for forward FFT, -1 for inverse FFT.
    implicit none
    ! Declare input/output variables
    complex(kind=8), dimension(:,:), intent(inout) :: data  ! 2D complex data array
    integer, intent(IN) :: isign                           ! FFT direction (+1 or -1)
    ! local parameters
    integer :: num_columns, index, step_size, swap_index
    integer :: sub_index, max_subsize, half_size
    ! FFT-related variables
    double precision :: rotation_angle
    complex(kind=8), dimension(size(data,1)) :: temp_values
    complex(kind=8):: twiddle_factor, twiddle_update
    complex(kind=8):: weight_factor
    double precision, PARAMETER :: PI = 3.1415926535897931d0

    ! Get the number of columns in the 2D array
    num_columns = size(data,2)

    ! Ensure number of columns is a power of 2 (required for FFT)
    if (iand(num_columns, num_columns - 1) /= 0) &
      call stopProgram('myfourrow num_columns must be a power of 2    ')

    ! Perform bit-reversal permutation
    half_size = num_columns / 2
    swap_index = half_size
    do index = 1, num_columns - 2
      if (swap_index > index) then
        call swap_complex_vectors(data(:, swap_index + 1), data(:, index + 1))
      endif

      sub_index = half_size
      do
        if (sub_index < 2 .or. swap_index < sub_index) exit
        swap_index = swap_index - sub_index
        sub_index = sub_index / 2
      enddo
      swap_index = swap_index + sub_index
    enddo

    ! Begin FFT computation using Danielson-Lanczos method
    max_subsize = 1
    do
      if (num_columns <= max_subsize) exit
      step_size = 2 * max_subsize

      ! Compute twiddle factors for FFT stage
      rotation_angle = PI / (isign * max_subsize)
      twiddle_update = cmplx(-2.0d0 * sin(0.5d0 * rotation_angle)**2, sin(rotation_angle), kind=8)
      twiddle_factor = cmplx(1.0d0, 0.0d0, kind=8)

      ! Apply FFT calculations to each subset of data
      do sub_index = 1, max_subsize
        weight_factor = twiddle_factor
        do index = sub_index, num_columns, step_size
          swap_index = index + max_subsize
          temp_values = weight_factor * data(:, swap_index)
          data(:, swap_index) = data(:, index) - temp_values
          data(:, index) = data(:, index) + temp_values
        enddo
        twiddle_factor = twiddle_factor * twiddle_update + twiddle_factor
      enddo

      max_subsize = step_size
    enddo

  end subroutine

  !-----------------------------------------------------------------------
  subroutine swap_complex_vectors(vector_a, vector_b)
  !-----------------------------------------------------------------------
  ! Swaps the contents of two complex single-precision vectors.
  !
  ! Arguments:
  !   vector_a - First complex vector (modified in place).
  !   vector_b - Second complex vector (modified in place).
  implicit none
  complex(kind=8), dimension(:), intent(inout) :: vector_a, vector_b  ! Input/output complex vectors
  ! local parameters
  complex(kind=8), dimension(size(vector_a)) :: temp_vector           ! Temporary storage for swapping

  ! Swap the contents of the two vectors
  temp_vector(:) = vector_a(:)
  vector_a(:) = vector_b(:)
  vector_b(:) = temp_vector(:)

  end subroutine

  end subroutine


!-----------------------------------------------------------------------
  subroutine FFT_real(data, N, isign, zdata)
!-----------------------------------------------------------------------
! Performs a real-to-complex or complex-to-real FFT using the
! Fast Fourier Transform (FFT) algorithm.
!
! Input:
!   data  - Real-valued input array of length N (must be a power of 2).
!   isign - +1 for forward transform, -1 for inverse transform.
!
! Output:
!   data  - Transformed data (in-place).
!   zdata - Optional complex array to store intermediate values.
!
! Dependencies:
!   - Uses `fft_complex` for FFT computation.
!   - Uses `nn_roots` for complex exponentials.
!
!-----------------------------------------------------------------------
  use precisions, only: WP
  implicit none
  ! Input/Output parameters
  integer, intent(in) :: N,isign
  real(WP), dimension(N), intent(inout) :: data
  complex(kind=8), dimension(N/2),intent(inout) :: zdata
  ! local parameters
  integer:: nhalf, nquarter
  complex(kind=8), dimension(N/4) :: w
  complex(kind=8), dimension(N/4-1) :: h1, h2
  complex(kind=8):: z
  double precision, parameter :: c1 = 0.5d0
  double precision :: c2
  interface
    subroutine FFT_complex(data,N,isign)
      implicit none
      integer, intent(in) :: N,isign
      complex(kind=8), dimension(N), intent(inout) :: data
    end subroutine
  end interface

  ! Ensure input size is a power of 2
  if (N == 0) call stopProgram('FFT_real: invalid N is zero    ')
  if (iand(N, N-1) /= 0) call stopProgram('FFT_real: N must be a power of 2    ')

  ! Define half and quarter lengths
  nhalf = N / 2
  nquarter = N / 4

  if (nhalf == 0) call stopProgram('FFT_real: invalid nhalf is zero    ')
  if (nquarter == 0) call stopProgram('FFT_real: invalid nquarter is zero    ')
  if (nhalf /= size(zdata)) call stopProgram('FFT_real: invalid zdata size (must be N/2)    ')

  ! initializes output arrays
  if (isign == 1) then
    ! forward FFT
    zdata(:) = cmplx(0.d0,0.d0,kind=8)
  else
    ! inverse FFT
    data(:) = 0.0_WP
  endif

  ! forward FFT
  if (isign == 1) zdata = cmplx(data(1:N-1:2), data(2:N:2), kind=8)

  ! Perform forward or inverse FFT
  if (isign == 1) then
    c2 = -0.5d0
    call FFT_complex(zdata,nhalf,+1)
  else
    c2 = 0.5d0
  endif

  ! Compute twiddle factors (roots of unity)
  w = nn_roots(sign(N, isign), nquarter)
  w = cmplx(-aimag(w), real(w), kind=8)

  ! Compute the Hermitian symmetry components
  h1 = c1 * (zdata(2:nquarter) + conjg(zdata(nhalf:nquarter+2:-1)))
  h2 = c2 * (zdata(2:nquarter) - conjg(zdata(nhalf:nquarter+2:-1)))

  ! Store transformed data
  zdata(2:nquarter) = h1 + w(2:nquarter) * h2
  zdata(nhalf:nquarter+2:-1) = conjg(h1 - w(2:nquarter) * h2)

  ! Special handling for zeroth frequency component
  z = zdata(1)
  if (isign == 1) then
    ! forward FFT
    zdata(1) = cmplx(real(z) + aimag(z), real(z) - aimag(z), kind=8)
  else
    ! inverse FFT
    zdata(1) = c1 * cmplx(real(z) + aimag(z), real(z) - aimag(z), kind=8)
    call FFT_complex(zdata,nhalf,-1)
  endif

  ! Copy data back to real array if necessary
  if (isign /= 1) then
    data(1:N-1:2) = real(zdata,kind=WP)
    data(2:N:2) = real(aimag(zdata),kind=WP)
  endif

contains

  !-----------------------------------------------------------------------
  function nn_roots(n_in, nn_in)
  !-----------------------------------------------------------------------
  ! Computes the nn-th roots of unity (complex exponential values)
  ! for a given integer n_in.
  !
  ! Input:
  !   n  - Base value defining the angular increment (e.g., TWOPI/n).
  !   nn - Number of roots to compute.
  !
  ! Output:
  !   nn_roots - Array of complex roots of unity.
  !
  ! The function constructs the roots iteratively, ensuring efficiency.
  !
  !-----------------------------------------------------------------------
    implicit none
    ! Input parameters
    integer, intent(in) :: n_in, nn_in
    ! Output array: Complex roots of unity
    complex(kind=8), dimension(nn_in) :: nn_roots
    ! local parameters
    integer :: k
    double precision :: theta
    double precision, parameter :: TWO_PI = 2.d0 * 3.1415926535897931d0

    ! initializes
    nn_roots(:) = cmplx(0.d0,0.d0)

    ! Initialize first root
    nn_roots(1) = 1.d0

    ! Compute angle increment
    theta = TWO_PI / n_in

    ! Iteratively compute the roots
    k = 1
    do
      if (k >= nn_in) exit

      ! Compute the k-th root
      nn_roots(k+1) = cmplx(cos(k * theta), sin(k * theta), kind=8)

      ! Fill in additional values using previous calculations
      nn_roots(k+2:min(2*k, nn_in)) = nn_roots(k+1) * nn_roots(2:min(k, nn_in-k))

      k = 2 * k  ! Double the index for next iteration
    enddo

    return
  end function

  end subroutine
