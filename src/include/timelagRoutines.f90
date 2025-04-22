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



!-----------------------------------------------------------------
module minimize_trace
!-----------------------------------------------------------------
  use precisions
  implicit none

  real(WP),allocatable,dimension(:,:):: seismoHom,seismoHet
  integer:: trace_length

contains

  !-----------------------------------------------------------------------
  real function costfunction(x)
  !-----------------------------------------------------------------------
  ! calculates the root-mean-square of the residual/difference between
  ! the shifted and reference trace
  !
  ! inputs:
  !     x     - array containing parameters t_lag and amplification to optimize
  !
  ! returns: rms of residuals of traces
    use precisions; use verbosity
    implicit none
    real,dimension(:),intent(in):: x
    ! local parameters
    real(WP):: seismoDiff(trace_length),seismoShift(trace_length)
    integer:: i,margin
    real(WP):: t_lag, amplification
    real(WP),external:: rootmeansquare
    ! rms inside the trace to avoid that the shifting produces artefacts
    logical,parameter:: REMOVE_MARGIN = .true.
    real,parameter:: MARGIN_PERCENT  = 0.1

    ! get parameters
    t_lag = x(1)
    amplification = x(2)
    if (amplification < 0.0) amplification = 0.0

    !debug
    !print *,'debug: costfunction: t_lag/amp = ',t_lag,amplification

    ! create shifted seismogram
    ! (move seismoHom, homogeneous trace, to seismoHet, the heterogeneous one)
    call sumSeismogram(t_lag,amplification,seismoHom,seismoHet,seismoShift,trace_length)

    ! get the rms of the difference between the two seismograms
    do i = 1,trace_length
      seismoDiff(i) = seismoHet(2,i) - seismoShift(i)
    enddo

    ! use a margin at beginning and end of seismogram
    if (REMOVE_MARGIN) then
      ! fixed margin
      margin = floor(trace_length * MARGIN_PERCENT)
      ! or depending on shift
      !margin = t_lag/abs(seismoHom(1,2)-seismoHom(1,1)) + 1
      if (margin > 1) then
        seismoDiff(1:margin) = 0.0
        seismoDiff(trace_length-margin:trace_length) = 0.0
      endif
    endif

    ! gets rms
    costfunction = rootmeansquare(seismoDiff,trace_length)
    return
  end function
end module


!-----------------------------------------------------------------
module module_spline
!du Numerical recipes
!-----------------------------------------------------------------
  implicit none
  public :: spline,splint
  private

contains

  !------------------------------------------------------------------
  subroutine spline(x,y,yp1,ypn,y2)
  !------------------------------------------------------------------
    implicit none
    double precision, dimension(:), intent(in) :: x,y
    double precision, intent(in) :: yp1,ypn
    double precision, dimension(:), intent(out) :: y2
    ! local parameters
    integer:: n
    double precision, dimension(size(x)) :: a,b,c,r

    n = size(x)
    if (n /= size(y) .or. n /= size(y2)) &
      call stopProgram('spline must have arrays with same size    ')

    c(1:n-1) = x(2:n)-x(1:n-1)
    r(1:n-1) = 6.0d0*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1) = r(2:n-1)-r(1:n-2)
    a(2:n-1) = c(1:n-2)
    b(2:n-1) = 2.0d0*(c(2:n-1)+a(2:n-1))
    b(1) = 1.0d0
    b(n) = 1.0d0
    if (yp1 > 0.99d30) then
       r(1) = 0.0d0
       c(1) = 0.0d0
    else
       r(1) = (3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       c(1) = 0.5d0
    endif
    if (ypn > 0.99d30) then
       r(n) = 0.0d0
       a(n) = 0.0d0
    else
       r(n) = (-3.0d0/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
       a(n) = 0.5d0
    endif
    call tridag_ser(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
  end subroutine

  !------------------------------------------------------------------
  double precision function splint(xa,ya,y2a,x)
  !------------------------------------------------------------------
    implicit none
    double precision, dimension(:), intent(in) :: xa,ya,y2a
    double precision, intent(in) :: x
    ! local parameters
    integer:: khi,klo,n,i
    double precision :: a,b,h

    n = size(xa)
    if (n /= size(ya) .or. n /= size(y2a)) &
      call stopProgram('splint: must have arrays with same size    ')

    klo = max(min(locate(xa,x),n-1),1)
    khi = klo+1
    h = xa(khi)-xa(klo)
    if (h == 0.0d0) then
          do i = 1,n
             print *,i,xa(i),khi,klo
          enddo
          call stopProgram('splint: bad xa input in splint    ')
    endif
    a = (xa(khi)-x)/h
    b = (x-xa(klo))/h
    splint = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi))*(h**2)/6.0d0
  end function

  !------------------------------------------------------------------
  subroutine tridag_ser(a,b,c,r,u)
  !------------------------------------------------------------------
    implicit none
    double precision, dimension(:), intent(in) :: a,b,c,r
    double precision, dimension(:), intent(out) :: u
    ! local parameters
    double precision, dimension(size(b)) :: gam
    integer:: n,j
    double precision :: bet

    n = size(a)+1
    if (n /= size(b) .or. n /= size(c)+1 .or. n /= size(r) .or. n /= size(u)) &
      call stopProgram('tridag_ser: must have same sized arrays    ')

    bet = b(1)
    if (bet == 0.0d0) &
      call stopProgram('tridag_ser: Error at code stage 1    ')

    u(1) = r(1)/bet
    do j = 2,n
       gam(j) = c(j-1)/bet
       bet = b(j)-a(j-1)*gam(j)
       if (bet == 0.0d0) &
            call stopProgram('tridag_ser: Error at code stage 2    ')
       u(j) = (r(j)-a(j-1)*u(j-1))/bet
    enddo
    do j = n-1,1,-1
       u(j) = u(j)-gam(j+1)*u(j+1)
    enddo
  end subroutine

  !-----------------------------------------------------------------
  integer function locate(xx,x)
  !-----------------------------------------------------------------
    implicit none
    double precision, dimension(:), intent(in) :: xx
    double precision, intent(in) :: x
    ! local parameters
    integer :: n,jl,jm,ju
    logical :: ascnd
    n = size(xx)
    ascnd = (xx(n) >= xx(1))
    jl = 0
    ju = n+1
    do
       if (ju-jl <= 1) exit
       jm = (ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    if (x == xx(1)) then
       locate = 1
    else if (x == xx(n)) then
       locate = n-1
    else
       locate = jl
    endif
  end function

end module


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
  use verbosity; use filterType
  implicit none
  real(WP),intent(out):: t_lag
  character(len=128),intent(in):: fileDelta,fileReference
  real(WP),intent(in):: startingTime
  ! local parameters
  integer:: entries,ier
  real(WP):: seismo(2,WindowSIZE),seismoRef(2,WindowSIZE),lastseismotime,timedelta
  real(WP),allocatable,dimension(:,:):: seismo1,seismo2

  !initialize
  seismo(:,:) = 0.0_WP
  seismoRef(:,:) = 0.0_WP

  ! read in from startingTime
  call readSeismo(seismoRef,WindowSIZE,fileReference,startingTime,timedelta,entries,lastseismotime)
  call readSeismo(seismo,WindowSIZE,fileDelta,startingTime,timedelta,entries,lastseismotime)

  allocate(seismo1(2,entries),seismo2(2,entries),stat=ier)
  if (ier /= 0) stop 'Error allocating seismo1/2'

  seismo1(:,:) = seismoRef(:,1:entries)
  seismo2(:,:) = seismo(:,1:entries)

  ! get time lag
  call getTimelagSeismos(t_lag,seismo2,seismo1,entries,timedelta)

  end subroutine


!-----------------------------------------------------------------------
  subroutine getTimelagTaped(t_lag,fileDelta,fileReference,startingTime,endingTime)
!-----------------------------------------------------------------------
! determines the time lag between two seismograms within a certain time window and given
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
  use verbosity; use filterType
  implicit none
  character(len=128),intent(in):: fileDelta,fileReference
  real(WP),intent(in):: startingTime,endingTime
  real(WP),intent(out):: t_lag
  ! local parameters
  integer:: entries,ier
  real(WP):: seismo(2,WindowSIZE),seismoRef(2,WindowSIZE),lastseismotime,timedelta
  real(WP),allocatable,dimension(:,:):: seismo1,seismo2

  !initialize
  seismo(:,:) = 0.0_WP
  seismoRef(:,:) = 0.0_WP

  ! read in from startingTime
  call readSeismoTaped(seismoRef,WindowSIZE,fileReference,startingTime,endingTime,timedelta,entries,lastseismotime)
  call readSeismoTaped(seismo,WindowSIZE,fileDelta,startingTime,endingTime,timedelta,entries,lastseismotime)

  allocate(seismo1(2,entries),seismo2(2,entries),stat=ier)
  if (ier /= 0) stop 'Error allocating seismo1/2'

  seismo1(:,:) = seismoRef(:,1:entries)
  seismo2(:,:) = seismo(:,1:entries)
  ! get time lag
  call getTimelagSeismos(t_lag,seismo2,seismo1,entries,timedelta)

  end subroutine


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
!
! returns: t_lag timelag in seconds
  use propagationStartup
  use verbosity; use filterType
  implicit none
  integer,intent(in):: station
  real(WP),intent(in):: startingTime
  real(WP),intent(out):: t_lag
  ! local parameters
  integer:: i,index,ier
  real(WP):: seismo(2,WindowSIZE),seismoRef(2,WindowSIZE),endtime,time,timeRef,displace,displaceRef
  real(WP),allocatable,dimension(:,:):: seismo1,seismo2

  !initialize
  seismo(:,:) = 0.0_WP
  seismoRef(:,:) = 0.0_WP

  ! read in from startingTime
  index = 0
  do i = 1,numofTimeSteps
    if (index < WindowSIZE) then
      if (manyKernels) then
        time = kernelsReceiversSeismogram(currentKernel,numofReceivers+1,i)
        displace = kernelsReceiversSeismogram(currentKernel,station,i)
        timeRef = kernelsReceiversSeismogramRef(currentKernel,numofReceivers+1,i)
        displaceRef = kernelsReceiversSeismogramRef(currentKernel,station,i)
      else
        time = receiversSeismogram(size(receivers)+1,i)
        displace = receiversSeismogram(station,i)
        timeRef = receiversSeismogramRef(size(receivers)+1,i)
        displaceRef = receiversSeismogramRef(station,i)
      endif

      if (time - startingTime > 0.0_WP) then ! starts with times from 'start time' on
        index = index+1
        seismo(1,index) = time
        seismo(2,index) = displace
        seismoRef(1,index) = timeRef
        seismoRef(2,index) = displaceRef
        endtime = time
      endif
    else
      !read until end is reached
      time = receiversSeismogram(size(receivers)+1,i) !,sourceterm
      endtime = time
    endif
  enddo

  allocate(seismo1(2,index),seismo2(2,index),stat=ier)
  if (ier /= 0) stop 'Error allocating seismo1/2'

  seismo1(:,:) = seismoRef(:,1:index)
  seismo2(:,:) = seismo(:,1:index)

  ! get time lag
  call getTimelagSeismos(t_lag,seismo2,seismo1,index,dt)

  end subroutine


!-----------------------------------------------------------------------
  subroutine getTimelagSeismos(t_lag,seismo1,seismo2,seismoLength,timedelta)
!-----------------------------------------------------------------------
! determines the time lag between two seismograms
! (both seismograms must have the same time step size)
!
! definition: seismo1 shifted to the right of seismo2, then timelag is positiv (i.e. waves arrive later in seismo1)
!
! input:
!       t_lag                   - time lag
!       seismo1,seismo2         - seismograms (array 2 x seismoLength, first index is time, second displacement)
!       seismoLength            - seismograms length (number of entries in seismogram)
!
! remember:
!       FILTER_SEISMOGRAMS      - use filter before determining time lag
! returns: t_lag timelag in seconds ( uses dt )
  use precisions
  use verbosity; use filterType
  use propagationStartup, only: datadirectory,dt
  implicit none
  real(WP),intent(out):: t_lag
  integer,intent(in):: seismoLength
  real(WP),intent(in):: seismo1(2,seismoLength),seismo2(2,seismoLength)
  real(WP),intent(in):: timedelta
  ! local parameters
  integer:: i,ier
  integer:: xcorrlength,zeropad
  real(WP):: seismo(2,WindowSIZE),seismoRef(2,WindowSIZE)
  real(WP):: sub
  real(WP):: timedeltaWindow,td
  real(WP),allocatable,dimension(:):: seismoWindow,seismoRefWindow,crosscorrelation

  real(WP),external:: getMaximum

  logical,parameter:: NORMALIZE         = .false.
  logical,parameter:: STRETCH           = .false.
  logical,parameter:: TAPER             = .false.
  logical,parameter:: TAPER_FILTERED    = .false.

  ! console
  if (beVerbose) then
    print *,'  getting timelag seismos by cross-correlation...'
    print *
    print *,'  perturbed seismogram:'
    print *,'    first entry: ',seismo1(1,1),seismo1(2,1)
    print *,'    last entry : ',seismo1(1,seismoLength),seismo1(2,seismoLength)
    print *,'    entry read : ',seismoLength
    print *,'  reference seismogram:'
    print *,'    first entry: ',seismo2(1,1),seismo2(2,1)
    print *,'    last entry : ',seismo2(1,seismoLength),seismo2(2,seismoLength)
    print *,'    entry read : ',seismoLength
    print *,'  time step: dt = ',timedelta
    print *
  endif

  ! set dt if not determined yet
  if (dt == 0.0) dt = timedelta

  ! debug output
  if (fileOutput) then
    print *,'  printing to file: ',trim(datadirectory)//'Timelag_read.dat'
    open(IOUT,file=trim(datadirectory)//'Timelag_read.dat')
    write(IOUT,*) "# format: #time #trace-perturbed #trace-reference"
    do i = 1,seismoLength
      write(IOUT,*) seismo1(1,i),seismo1(2,i),seismo2(2,i)
    enddo
    close(IOUT)
  endif

  ! copy arrays
  seismo(:,:) = 0.0_WP
  seismoRef(:,:) = 0.0_WP
  do i = 1,seismoLength
    seismo(:,i) = seismo1(:,i)
    seismoRef(:,i) = seismo2(:,i)
  enddo

  ! tapering ends
  if (TAPER) then
    if (beVerbose) print *,'    taper ends'
    call taperSeismogram(seismo(:,1:seismoLength),seismoLength,seismoLength,beVerbose)
    call taperSeismogram(seismoRef(:,1:seismoLength),seismoLength,seismoLength,beVerbose)
  endif

  ! filtering necessary
  if (FILTER_SEISMOGRAMS) then
    if (beVerbose) print *,'    filter seismograms'
    call dofilterSeismogram(seismo,WindowSIZE)
    call dofilterSeismogram(seismoRef,WindowSIZE)

    ! apply hanning window  to smooth seismograms ends
    if (TAPER_FILTERED) then
      if (beVerbose) print *,'    taper seismogram ends'
      call taperSeismogram(seismo,WindowSIZE,WindowSIZE,beVerbose)
      call taperSeismogram(seismoRef,WindowSIZE,WindowSIZE,beVerbose)
    endif
  endif

  ! length of cross-correlation array
  !initialize
  xcorrlength = WindowSIZE             ! WindowSize is power of 2 for FFT
  zeropad = WindowSIZE - seismoLength  ! zero padding

  ! increases array lengths to have a minimum of twice the actual seismo length
  if (zeropad < seismoLength) then
    xcorrlength = 2 * WindowSIZE
    zeropad = xcorrlength - seismoLength
  endif

  allocate(seismoWindow(xcorrlength), &
           seismoRefWindow(xcorrlength), &
           crosscorrelation(xcorrlength),stat=ier)
  if (ier /= 0) stop 'Error allocating correlation array'
  seismoWindow(:) = 0.0_WP
  seismoRefWindow(:) = 0.0_WP
  crosscorrelation(:) = 0.0_WP

  ! fills in the values into the window arrays
  if (STRETCH) then
    ! interpolated into window size
    td = (seismo1(1,seismoLength)-seismo1(1,1))/(seismoLength-1.0)
    timedeltaWindow = (seismo1(1,seismoLength)-seismo1(1,1))/(WindowSIZE-1.0)
    print *,'    stretched timedelta:',timedeltaWindow,WindowSIZE
    call resample(seismo1(2,1:seismoLength),seismoLength,td,seismoWindow,WindowSIZE, &
               timedeltaWindow )
    call resample(seismo2(2,1:seismoLength),seismoLength,td,seismoRefWindow,WindowSIZE, &
               timedeltaWindow)
  else
    ! window zeros are added
    if (WindowSIZE < seismoLength) stop 'Error windowsize/seismoLength'
    seismoWindow(1:WindowSIZE) = seismo(2,1:WindowSIZE)
    seismoRefWindow(1:WindowSIZE) = seismoRef(2,1:WindowSIZE)
  endif

  ! normalize seismograms
  if (NORMALIZE) then
    print *,'    normalize traces'
    call normalizeArray(seismoWindow,xcorrlength)
    call normalizeArray(seismoRefWindow,xcorrlength)
  endif

  ! debug output
  if (fileOutput) then
    print *,'  printing to file: ',trim(datadirectory)//'Timelag_window.dat'
    open(IOUT,file=trim(datadirectory)//'Timelag_window.dat')
    write(IOUT,*) "# format: #time #windowed-trace #windowed-reference-trace"
    do i = 1,xcorrlength
      if (STRETCH) then
        write(IOUT,*) (i-1)*timedeltaWindow+seismo2(1,1),seismoWindow(i),seismoRefWindow(i)
      else
        write(IOUT,*) (i-1)*timedelta+seismo2(1,1),seismoWindow(i),seismoRefWindow(i)
      endif
    enddo
    close(IOUT)
  endif

  ! gets cross-correlation
  call correlation_traces(seismoWindow,seismoRefWindow,crosscorrelation,xcorrlength)

  ! alternative
  !call time_correl(seismoWindow,seismoRefWindow,crosscorrelation,WindowSIZE)

  ! debug output
  if (fileOutput) then
    print *,'  printing to file: ',trim(datadirectory)//'Timelag_correlation.dat'
    open(IOUT, file=trim(datadirectory)//'Timelag_correlation.dat')
    write(IOUT,*) '# correlation'
    write(IOUT,*) "# format: #i #cross-correlation"
    do i = 1,xcorrlength
      write(IOUT,*) i,crosscorrelation(i)
    enddo
    close(IOUT)
  endif

  if (beVerbose) then
    print *,'  correlation array length = ',xcorrlength
  endif

  ! calculate the subsample precision maximum position
  sub = getMaximum(crosscorrelation,size(crosscorrelation))

  ! check with bounds: numerical recipes ( http://hebb.mit.edu/courses/9.29/2002/readings/c13-2.pdf )
  ! Just as in the case of convolution we have to consider end effects, since our
  ! data will not, in general, be periodic as intended by the correlation theorem. Here
  ! again, we can use zeropadding. If you are interested in the correlation for lags as
  ! large as
  ! +/- K, then you must append a buffer zone of K zeros at the end of both
  ! input datasets. If you want all possible lags from N datapoints (not a usual thing),
  ! then you will need to pad the data with an equal number of zeros; this is the extreme
  ! case.
  if (abs(sub) > zeropad) then
    print *,'Error: correlation returns lag:',sub
    print *,'       which is bigger than the buffer zone:',zeropad
    stop 'Error correlation'
  endif

  !only be verbose once
  if (beVerbose) then
    beVerbose = .false.
    fileOutput = .false.
  endif

  ! set timelag
  if (STRETCH) then
    ! windowed/streched traces
    t_lag = sub*timedeltaWindow
  else
    t_lag = sub*timedelta
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine time_correl(trace1,trace2,corr,nbsteps)
!-----------------------------------------------------------------------
  use precisions
  implicit none
  integer,intent(in):: nbsteps
  real(WP),intent(in):: trace1(nbsteps),trace2(nbsteps)
  real(WP),intent(out):: corr(nbsteps)
  ! local parameters
  integer:: m,n
  ! initiialize
  corr(:) = 0.0_WP
  ! cross-correlation in time domain where both traces are real
  do m = 1,nbsteps/2
    do n = 1,nbsteps-m+1
      corr(m) = corr(m) + trace1(n)*trace2(m-1+n)
    enddo
  enddo

  do m = nbsteps/2+1,nbsteps
    do n = 1,nbsteps-m+1
      corr(m) = corr(m) + trace1(n)*trace2(m-1+n)
    enddo
  enddo

  end subroutine


!-----------------------------------------------------------
  subroutine resample_alt(trace1,size1,dta,traceout,sizeout,dtb)
!-----------------------------------------------------------
  use precisions; use splinefunction
  implicit none
  integer,intent(in) :: size1,sizeout
  real(WP),dimension(size1), intent(in) :: trace1
  real(WP),intent(in) :: dta,dtb
  real(WP),dimension(sizeout), intent(out):: traceout
  ! local parameters
  integer:: i,length,ier
  double precision:: loc
  double precision,external:: cubicspline_eval

  ! spline arrays
  length = size1

  i1 = 1
  i2 = length

  if (allocated(X)) deallocate(X,Y,Q,F)
  allocate(X(i2),Y(i2),Q(3,i2),F(3,i2),stat=ier)
  if (ier /= 0) call stopProgram('resample error')

  ! set spline arrays
  do i = i1,i2
    X(i) = dble((i-1)*dta)
    Y(i) = dble(trace1(i))
  enddo
  call cubicspline_setup()

  ! interpolate values
  do i = 1,sizeout
     loc = (i-1)*dtb
     if (loc <= X(size1) .and. loc >= X(1)) then
        ! spline respresentation of trace at location loc
        traceout(i) = cubicspline_eval(loc)
     else
        traceout(i) = 0.0
     endif
  enddo
  ! free memory
  deallocate(X,Y,Q,F)

  end subroutine


!-----------------------------------------------------------
  subroutine resample(trace1,size1,dta,traceout,sizeout,dtb)
!-----------------------------------------------------------
  !use born_driver
  use precisions
  use module_spline, only: spline,splint
  implicit none
  integer,intent(in) :: size1,sizeout
  real(WP),dimension(size1), intent(in) :: trace1
  real(WP),intent(in) :: dta,dtb
  real(WP),dimension(sizeout), intent(out):: traceout
  ! local parameters
  double precision, dimension(size1):: ya
  real, dimension(sizeout):: yb
  double precision, dimension(size1) :: xa,y2
  integer :: i
  double precision :: yp1,ypn,x
  double precision,parameter:: EPS = 0.01

  ! fill arrays
  ya(:) = trace1(:)
  do i = 1,size1
     xa(i) = (i-1)*dta
  enddo

  yp1 = (ya(2)-ya(1)    )/dta
  ypn = (ya(size1)-ya(size1-1))/dta
  call spline(xa,ya,yp1,ypn,y2)

  ! get interpolated values
  do i = 1,sizeout
     x = (i-1)*dtb
     if (x <= xa(size1)+EPS .and. x >= xa(1)-EPS) then
        yb(i) = real(splint(xa,ya,y2,x))
     else
        yb(i) = 0.0
     endif
  enddo

  traceout(:) = yb(:)

  end subroutine


!-----------------------------------------------------------------------
  function getMaximum(crosscorrelation,length)
!-----------------------------------------------------------------------
! determines the maximum (positive) position with subsample precision
!
! input:
!       crosscorrelation                            - cross correlation values
!       length                                           - array size
!
! see remark from function correl():
! Computes the correlation of two real datasets data1 and data2 of length N (including any user-supplied zeropadding).
! N must be an integer power of 2. The answer is returned as the function correl, an array of length N.
! The answer is stored in wrap-around order, i.e., correlations at increasingly negative lags are in correl(N) on down to
! correl(N/2+1), while correlations at increasingly positive lags are in correl(1) (zero lag) on up to correl(N/2).
! Sign convention of this routine: if data1 lags data2, i.e., is shifted to the right of it, then correl will show
! a peak at positive lags.
!
! i found a problem when i take the absolute maximum of the cross-correlation. the
! traces could be shifted by PI, so that their "anti"-correlation was maximum.
! this routine now looks at the maximum positive correlation coefficient/value.
!
! returns: subsample precision maximum
  use verbosity;use precisions
  implicit none
  integer,intent(in):: length
  real(WP),intent(in):: crosscorrelation(length)
  real(WP):: getMaximum
  ! local parameters
  integer:: indexmax   !,i
  real(WP):: max  !,correlation
  real(WP):: pointleft,pointright,maxestimate,peaklocation
  !real(WP):: abscorr(length)
  integer,dimension(1):: indexmaxloc

  !Searches maximum correlation
  !max = abs(crosscorrelation(1))
  !indexmax = 1
  !do i=2,length
  !  correlation=abs(crosscorrelation(i))
  !  if (correlation>max) then
  !    max = correlation
  !    indexmax = i
  !  endif
  !enddo

  ! or uses intrinsic function to do so
  !!!! absolute maximum: see problem description from above
  !!!! indexmaxloc = maxloc( abs(crosscorrelation(:)) )
  ! positive maximum
  indexmaxloc = maxloc( crosscorrelation(:) )
  indexmax = indexmaxloc(1)

  max = crosscorrelation(indexmax)  ! could be negative maximum
  if (beVerbose) then
    print *,'    maximum correlation  :',max
    print *,'        index max        :',indexmax
  endif

  !subsample precision
  ! first half (positive lag)
  if (indexmax > 1 .and. indexmax <= length/2) then
    pointleft = crosscorrelation(indexmax-1)
    pointright = crosscorrelation(indexmax+1)
  endif
  if (indexmax == 1) then
    pointleft = crosscorrelation(length)
    pointright = crosscorrelation(indexmax+1)
  endif
  !if (indexmax == length/2) then
  !  pointleft=crosscorrelation(indexmax-1)
  !  pointright=crosscorrelation(length/2+1)
  !endif

  ! second half (negative lag)
  if (indexmax >= length/2+1 .and. indexmax < length) then
    pointleft = crosscorrelation(indexmax-1)
    pointright = crosscorrelation(indexmax+1)
  endif
  if (indexmax == length) then
    pointleft = crosscorrelation(indexmax-1)
    pointright = crosscorrelation(1)
  endif
  !if (indexmax == length/2+1) then
  !  pointleft=crosscorrelation(length/2)
  !  pointright=crosscorrelation(indexmax+1)
  !endif

  ! console output debug
  if (beVerbose) then
    print *,'        pointleft        : ',pointleft
    print *,'        pointright       : ',pointright
  endif

  ! respect wrap-around order of the returned crosscorrelation array by
  ! numerical recipes functions (see correl()  or
  ! see: http://hebb.mit.edu/courses/9.29/2002/readings/c13-2.pdf  )
  if (indexmax > length/2) then
    ! the correlation at lag -1 is in r_N-1, the last component; etc.
    indexmax = indexmax-length-1
  else
    ! index goes from 1 to N/2: correlation at zero lag is in r_0, the first component,
    ! the correlation at lag 1 is in r_1, the second component;
    indexmax = indexmax-1
  endif

  ! quadratic interpolation will give this maximum
  ! ( http://www-ccrma.stanford.edu/~jos/parshl/Peak_Detection_Steps_3.html )
  peaklocation = 0.5*(pointleft-pointright)/(pointleft - 2.0*max + pointright)
  getMaximum = indexmax + peaklocation

  ! console output debug
  if (beVerbose) then
    maxestimate = max-0.25*(pointleft-pointright)*getMaximum
    print *,'    subsample position   : ',getMaximum
    print *,'        index            : ',indexmax
    print *,'        estimated maximum: ',maxestimate
  endif

  return
  end function


!-----------------------------------------------------------------------
  subroutine readSeismo(seismo,size,fileName,start,timediff,entries,lastseismotime)
!-----------------------------------------------------------------------
! reads in the seismogram
!
! input:
!     seismo        - data array
!     size          - data array size
!     fileName      - data file name
!     start         - time position from which on data shall be read
!     timediff      - time step dt
!     entries       - (optional) number of entries
! returns: seismo array and dt
  use verbosity; use precisions
  implicit none
  integer,intent(in):: size
  character(len=128),intent(in):: fileName
  real(WP),intent(in):: start
  real(WP),intent(out):: seismo(2,size)
  real(WP),intent(out):: timediff,lastseismotime
  integer,intent(out):: entries
  ! local parameters
  integer:: ier,i
  character(len=128):: line
  real(WP):: time,displace,endtime
  integer:: index, offset,ZEROLINE,FILELINES !tests: ZEROLINE/0/ ! approximate time 0 s line

  !initialize
  ZEROLINE = 0
  FILELINES = 1100
  seismo(:,:) = 0.0_WP

  ! open seismogram file
  if (beVerbose) then
    print *,'  file : ',trim(fileName)
    print *,'  start: ',start
  endif

  open(IIN, file=trim(fileName),status='old',iostat=ier)
  if (ier /= 0) then
    print *,'Error: opening file ',trim(fileName)
    call stopProgram( 'abort - readSeismo   ')
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
  if (start < 0.0) ZEROLINE= offset
  if (beVerbose) then
    print *,'    file lines      : ',FILELINES
    print *,'    zero line       : ',ZEROLINE
    print *,'    offset          : ',offset
  endif
  rewind(IIN)

  ! parse file for displacement values
  timediff = 0.0
  index = 0
  ier = 0
  do i = 1, FILELINES
    if (i >= ZEROLINE .and. index < size) then
      read(IIN, *, iostat=ier) time,displace !,sourceterm
      if (ier /= 0) then
        print *,'Error: reading input. last line ',i,ZEROLINE,index
        call stopProgram( 'abort - readSeismo   ')
      endif

      if (time - start >= 0.0) then ! starts with times from 'start time' on
        index = index+1
        seismo(1,index) = time
        seismo(2,index) = displace
        !if (beVerbose.and.(index==1.or.index==size))               &
        !&                print *,index,seismo(:,index)
        endtime = time

        !time step
        if (abs(timediff) < 0.000001 .and. index > 1) then
          timediff = seismo(1,index)-seismo(1,index-1)
        endif
      endif
    else
      !read until end is reached
      if (index >= size) then
        !read values
        read(IIN, *, iostat=ier) time,displace !,sourceterm
        endtime = time
      else
        !read a text line
        read(IIN, *, iostat=ier) line
        if (ier /= 0) then
          !print *,'error reading input. last line ',i
          exit
        endif
      endif
    endif
  enddo
  close(IIN)

  if (beVerbose) then
    print *,'    timestep dt     : ',timediff
    print *,'    first entry     : ',seismo(1,1),seismo(2,1)
    print *,'    last entry      : ',seismo(1,index),seismo(2,index)
    print *,'    entry read      : ',index
    print *,'    last time       : read/file = ',seismo(1,index),'/',endtime
    print *
  endif

  ! number of file entries read in
  entries = index

  ! end time in file
  lastseismotime = endtime

  end subroutine


!-----------------------------------------------------------------------
  subroutine readSeismoTaped(seismo,size,fileName,start,ending,timediff,entries,lastseismotime)
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
  integer,intent(in):: size
  character(len=128),intent(in):: fileName
  real(WP),intent(in):: start,ending
  real(WP),intent(out):: seismo(2,size)
  real(WP),intent(out):: timediff,lastseismotime
  integer,intent(out):: entries
  ! local parameters
  integer:: ier,i
  character(len=128):: line
  real(WP):: time,displace,endtime
  integer:: index,offset,ZEROLINE,FILELINES !tests: ZEROLINE/0/ ! approximate time 0 s line

  !initialize
  ZEROLINE = 0
  FILELINES = 1100
  seismo(:,:) = 0.0_WP

  ! open seismogram file
  if (beVerbose) then
    print *,'  file : ',trim(fileName)
    print *,'  start: ',start
  endif

  open(IIN, file=trim(fileName),status='old',iostat=ier)
  if (ier /= 0) then
    print *,'Error: opening file ',trim(fileName)
    call stopProgram( 'abort - readSeismo   ')
  endif

  ! parse file for line with zero time and number of lines
  index = 0
  offset = 0
  do while( ier == 0 )
    read(IIN,*,iostat=ier) line
    index = index+1

    ! suppress leading white spaces, if any
    line = adjustl(line)
    line = trim(line)

    ! check for additional data in file
    if (line(1:3) == 'dis') offset = 3
    if (line(1:3) == 'rec') offset = 2

    ! check for line with time 0.0
    if (line(1:3) == '0.0') ZEROLINE = index

  enddo
  FILELINES = index-1
  if (start < 0.0) ZEROLINE= offset
  if (beVerbose) then
    print *,'    file lines      : ',FILELINES
    print *,'    zero line       : ',ZEROLINE
    print *,'    offset          : ',offset
  endif
  rewind(IIN)

  ! parse file for displacement values
  timediff = 0.0
  index = 0
  ier = 0
  do i = 1, FILELINES
    if (i >= ZEROLINE .and. index < size) then
      read(IIN, *, iostat=ier) time,displace !,sourceterm
      if (ier /= 0) then
        print *,'Error: reading input. last line ',i,ZEROLINE,index
        call stopProgram( 'abort - readSeismo   ')
      endif

      if (time - start >= 0.0 .and. ending - time >= 0.0) then ! starts with times from 'start time' on
        index = index+1
        seismo(1,index) = time
        seismo(2,index) = displace
        !if (beVerbose.and.(index==1.or.index==size))               &
        !&                print *,index,seismo(:,index)
        endtime = time

        !time step
        if (abs(timediff) < 0.000001 .and. index > 1) then
          timediff = seismo(1,index)-seismo(1,index-1)
        endif
      endif
    else
      !read until end is reached
      if (index >= size) then
        !read values
        read(IIN, *, iostat=ier) time,displace !,sourceterm
        endtime = time
      else
        !read a text line
        read(IIN, *, iostat=ier) line
        if (ier /= 0) then
          !print *,'error reading input. last line ',i
          exit
        endif
      endif
    endif
  enddo
  close(IIN)

  if (beVerbose) then
    print *,'    timestep dt     : ',timediff
    print *,'    first entry     : ',seismo(1,1),seismo(2,1)
    print *,'    last entry      : ',seismo(1,index),seismo(2,index)
    print *,'    entry read      : ',index
    print *,'    last time       : read/file = ',seismo(1,index),'/',endtime
    print *
  endif

  ! number of file entries read in
  entries = index

  ! end time in file
  lastseismotime = endtime

  end subroutine


!-----------------------------------------------------------------------
  subroutine getGridvalues()
!-----------------------------------------------------------------------
! read in precalculated cell areas, edge lengths and center distances
  use propagationStartup;use cells
  implicit none
  ! local parameters
  integer:: ier

  ! read in vertices array
  call readData(.true.)

  ! new arrays for precalculated cell attributes
  allocate(cellAreas(numVertices), &
           cellEdgesLength(numVertices,0:6), &
           cellCenterDistances(numVertices,0:6), stat=ier )
  if (ier /= 0) then
    print *,'Error: in allocating arrays for cell area,..'
    call stopProgram( 'abort - getGridvalues    ')
  endif

  if (CORRECT_RATIO) then
    allocate(cellFractions(numVertices,0:6),stat=ier)
    if (ier /= 0) call stopProgram('error allocating cellFractions')
  endif

  ! read in cell areas
  call readPrecalculated()

  end subroutine


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
  character(len=128),intent(in):: filename
  real(WP),intent(out):: startTime
  ! local parameters
  integer:: numEntries,entries1
  real(WP):: endtime,timedelta,defaultStart
  real(WP),dimension(:,:),allocatable:: seismo

  !initialize
  numEntries = 10000
  defaultStart = FIRSTTIME
  startTime = 0.0_WP

  ! read in seismograms (hopefully complete)
  allocate(seismo(2,numEntries))
  seismo(:,:) = 0.0_WP

  call readSeismo(seismo,numEntries,filename,defaultStart,timedelta,entries1,endtime)

  seismo_timestep = timedelta

  ! determine startTime and arrivaltime
  call getStartTimeSeismogram(seismo,numEntries,endtime,timedelta,startTime)

  ! free memory
  deallocate(seismo)

  end subroutine


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
  use precisions
  use traveltime; use verbosity; use propagationStartup; use filterType
  implicit none
  integer,intent(in):: seismoLength
  real(WP),intent(in):: seismo(2,seismoLength)
  real(WP),intent(in):: endtime,timedelta
  real(WP),intent(out):: startTime
  ! local parameters
  integer:: i
  real(WP):: leftWindowTime,defaultStart

  ! initialize
  defaultStart = FIRSTTIME
  startTime = 0.0_WP

  ! get time when seismogram is above threshold
  arrivalTime = defaultStart
  do i = 1,seismoLength
    if (seismo(2,i) > ARRIVAL_THRESHOLD) then
      arrivalTime = seismo(1,i)
      exit
    endif
  enddo

  ! adapt start time for fourier transform depending on arrival time and window size of fourier transformation
  ! to get at least all entries after the arrival time
  if (seismoLength > WindowSIZE) then
    ! count backward to have starting time of window
    leftWindowTime = endtime-timedelta*WindowSIZE

    if (leftWindowTime < arrivalTime) then
      startTime = leftWindowTime - timedelta/2.0_WP
    else
      if (beVerbose) then
        print *,'window size not reaching end of seismogram...',arrivalTime,WindowSIZE
        print *
      endif
      ! start when wave arrives
      startTime = arrivalTime - timedelta/2.0_WP
    endif
  else
    if (beVerbose) then
      print *,'  defaultStart time:',defaultStart
    endif
    ! set to default start time
    startTime = defaultStart
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine getPerturbedSeismo(seismoRef,seismoDelta,seismoPerturbed,length)
!-----------------------------------------------------------------------
! calculates the difference between the reference and delta seismogram
! input:
!       seismoRef                             - reference seismogram
!       seismoDelta                           - delta seismogram
!       seismoPerturbed                   - perturbations
! returns: seismoPerturbed
  use precisions
  implicit none
  integer,intent(in):: length
  real(WP),intent(in):: seismoRef(length),seismoDelta(length)
  real(WP),intent(out):: seismoPerturbed(length)
  ! local parameters
  integer:: i

  do i = 1,length
    seismoPerturbed(i) = seismoDelta(i) - seismoRef(i)
  enddo

  end subroutine


!-----------------------------------------------------------------------
  subroutine getFirstDerivative(seismo,seismoOut,timedelta,length)
!-----------------------------------------------------------------------
! calculates the first time derivative
! input:
!       seismo                             - reference seismogram
!       seismoOut                       - seismogram, derivated
!       timedelta                          - time step size
!       length                              - length of both seismograms seismo & seismoOut
! returns: first derivative in seismoOut
  use precisions; use splinefunction
  use parallel; use verbosity
  implicit none
  integer,intent(in):: length
  real(WP),intent(in):: seismo(length)
  real(WP),intent(out):: seismoOut(length)
  real(WP),intent(in):: timedelta
  ! local parameters
  integer:: i,ier
  double precision:: err
  !double precision,external:: splineRepresentation
  double precision,external:: cubicspline_derivative
  !interface
  !  double precision function splineRepresentation(location)
  !  use splinefunction, only: i1,i2,X,Y,Q
  !  implicit none
  !  double precision, INTENT(IN) :: location
  !  double precision,external:: cubicspline_eval
  !  end function splineRepresentation
  !end interface
  logical,parameter:: FINITE_DIFFERENCES = .true.

  ! check
  if (length < 2) then
    print *,'Error: derivative undefined!',length,timedelta
    call stopProgram('abort getFirstDerivative() - seismo too short!    ')
  endif

  ! simplest finite difference
  if (FINITE_DIFFERENCES) then
    do i = 2,length-1
      ! central differences scheme
      seismoOut(i) = ( seismo(i+1)-seismo(i-1) )/(2.0*timedelta)
    enddo
    seismoOut(1) = ( seismo(2)-seismo(1) )/timeDelta
    seismoOut(length) = ( seismo(length)-seismo(length-1) )/timeDelta
    return
  endif

  ! for spline representation
  i1 = 1
  i2 = length

  if (allocated(X)) deallocate(X,Y,Q,F)
  allocate(X(i2),Y(i2),Q(3,i2),F(3,i2),stat=ier)
  if (ier /= 0) call stopProgram('getFirstDerivative() - error allocating spline arrays    ')

  ! get spline representation
  do i = i1,i2
    X(i) = dble((i-1)*timedelta)
    Y(i) = dble(seismo(i))
  enddo
  call cubicspline_setup()

  ! get time derivative
  do i = i1,i2
    seismoOut(i) = real(cubicspline_derivative( X(i),dble(2*timedelta),err ),kind=WP)

    ! check
    if (MAIN_PROCESS .and. VERBOSE) then
      ! notifies in case of 10% error
      if (abs(seismoOut(i)) > 0.1) then
        if (abs(err/seismoOut(i)) > 0.1) then
          print *,'  derivative has big error:'
          print *,'    index/error               :  ',i,err
          print *,'    index range              : ',i1,i2
          print *,'    seismo/seismoDeriv :  ',seismo(i),seismoOut(i)
        endif
      endif
    endif
  enddo

  ! check
  if (MAIN_PROCESS .and. VERBOSE) then
    if (maxval(abs(seismoOut(:))) < 0.0001) then
      print *,'  check time series min/max:',minval(seismoOut(:)),maxval(seismoOut(:))
      print *,'    index range              : ',i1,i2
      print *,'    length:',length
    endif
  endif

  ! free memory
  deallocate(X,Y,Q,F)

  end subroutine


!-----------------------------------------------------------------------
  subroutine getSecondDerivative(seismo,seismoOut,timedelta,length)
!-----------------------------------------------------------------------
! calculates the first time derivative
! input:
!       seismo                             - reference seismogram
!       seismo2                           - derivated seismogram
!       timedelta                                      - time step size
! returns: seismo2
  use precisions
  implicit none
  integer,intent(in):: length
  real(WP),intent(in):: seismo(length)
  real(WP),intent(out):: seismoOut(length)
  real(WP),intent(in):: timedelta
  ! local parameters
  real(WP):: seismo1(length)

  ! first time derivative
  call getFirstDerivative(seismo,seismo1,timedelta,length)
  ! second time
  call getFirstDerivative(seismo1,seismoOut,timedelta,length)

  end subroutine


!-----------------------------------------------------------------------
  subroutine getIntegral(seismo1,seismo2,integral,timedelta,length)
!-----------------------------------------------------------------------
! calculates the integral of the product of two seismograms
! input:
!       seismo1,seismo2                     - reference seismograms
!       integral                                     - integral value
!       timedelta                                              - time step size
! returns: integral
  use precisions; use splinefunction
  use parallel
  implicit none
  integer,intent(in):: length
  real(WP),intent(in):: seismo1(length),seismo2(length)
  real(WP),intent(in):: timedelta
  real(WP),intent(out):: integral
  ! local parameters
  integer:: i,ier
  interface
    function myspline_func(location)
    implicit none
    double precision, dimension(:),intent(in) :: location
    double precision, dimension(size(location)):: myspline_func
    end function
  end interface
  double precision,external :: simpsons_integral_vecFunc
  logical,parameter:: FINITE_DIFFERENCES = .true.

  ! initialize
  integral = 0.0_WP

  ! simplest finite difference
  if (FINITE_DIFFERENCES) then
    do i = 1,length-1
      integral = integral + (seismo1(i)*seismo2(i)+seismo1(i+1)*seismo2(i+1))/2.0*timedelta
    enddo
    return
  endif

  ! for spline representation
  i1 = 1
  i2 = length

  if (.not. allocated(X)) then
    allocate(X(i2),Y(i2),Q(3,i2),F(3,i2),stat=ier)
    if (ier /= 0) call stopProgram('error allocating spline arrays    ')
  else
    ! check with previous spline representation
    if (i2 /= size(Y)) then
     print *,'Error: different array lengths:',length,size(Y)
     call stopProgram( 'abort - getIntegral   ')
    endif
  endif

  ! get spline representation
  do i = i1,i2
    X(i) = dble((i-1)*timedelta)
    Y(i) = dble(seismo1(i)*seismo2(i))
  enddo
  call cubicspline_setup()

  ! get integral value
  integral = real(simpsons_integral_vecFunc(myspline_func, X(i1), X(i2)),kind=WP)

  !debug
  !  val1 = 0.0
  !  do i=i1,i2
  !    val1 = val1 + Y(i)*timedelta
  !  enddo
  !  print *,'integral: ',i1,i2
  !  print *,'    X:',X(i1),X(i2)
  !  print *,'    Y:',Y(1),Y(2),Y(3)
  !  print *,'    value=',integral
  !  print *,'    approx value = ',val1
  !  print *

  ! free memory
  deallocate(X,Y,Q,F)

  end subroutine


!-----------------------------------------------------------------------
  subroutine getAnalyticalTimelag(t_lag,fileDelta,fileReference,startingTime)
!-----------------------------------------------------------------------
! determines the time lag between two filtered seismograms given
! their filenames in an analytic approach
! (both seismograms must have the same time step size)
! input:
!       t_lag                             - time lag
!       fileDelta,fileReference           - file names of seismograms
!       startingTime                      - window start time of fourier transformation for correlation
!       bandwidth,power                   - butterworth filter parameters
!       waveperiod                        - wave period in seconds of butterworth bandwidth frequency
! returns: t_lag timelag in seconds
  use verbosity, only: beVerbose
  use filterType, only: WindowSIZE
  use precisions
  implicit none
  real(WP),intent(out):: t_lag
  character(len=128),intent(in):: fileDelta,fileReference
  real(WP),intent(in):: startingTime
  ! local parameters
  integer:: entries,ier
  real(WP):: timedelta,lastseismotime
  real(WP):: seismo(2,WindowSIZE),seismoRef(2,WindowSIZE)
  real(WP),allocatable,dimension(:,:):: seismo1,seismo2

  ! read in complete seismograms
  if (beVerbose) print *,'analytical reading in... '
  call readSeismo(seismoRef,WindowSIZE,fileReference,startingTime,timedelta,entries,lastseismotime)
  call readSeismo(seismo,WindowSIZE,fileDelta,startingTime,timedelta,entries,lastseismotime)

  allocate(seismo1(2,entries),seismo2(2,entries),stat=ier)
  if (ier /= 0) stop 'Error allocating seismo1/2'

  seismo1(:,:) = seismoRef(:,1:entries)
  seismo2(:,:) = seismo(:,1:entries)

  ! return time lag calculated analytically
  call getAnalyticalTimelagSeismos(t_lag,seismo1,seismo2,entries,timedelta)

  end subroutine


!-----------------------------------------------------------------------
  subroutine getAnalyticalTimelagTaped(t_lag,fileDelta,fileReference,startingTime,endingTime)
!-----------------------------------------------------------------------
! determines the time lag between two filtered seismograms given
! their filenames in an analytic approach
! (both seismograms must have the same time step size)
! input:
!       t_lag                             - time lag
!       fileDelta,fileReference           - file names of seismograms
!       startingTime                      - window start time of fourier transformation for correlation
! returns: t_lag timelag in seconds
  use verbosity, only: beVerbose
  use filterType, only: WindowSIZE
  use precisions
  implicit none
  real(WP),intent(out):: t_lag
  character(len=128),intent(in):: fileDelta,fileReference
  real(WP),intent(in):: startingTime,endingTime
  ! local parameters
  integer:: entries,ier
  real(WP):: timedelta,lastseismotime
  real(WP):: seismo(2,WindowSIZE),seismoRef(2,WindowSIZE)
  real(WP),allocatable,dimension(:,:):: seismo1,seismo2

  ! read in complete seismograms
  if (beVerbose) print *,'analytical reading in... '
  call readSeismoTaped(seismoRef,WindowSIZE,fileReference,startingTime,endingTime,timedelta,entries,lastseismotime)
  call readSeismoTaped(seismo,WindowSIZE,fileDelta,startingTime,endingTime,timedelta,entries,lastseismotime)

  allocate(seismo1(2,entries),seismo2(2,entries),stat=ier)
  if (ier /= 0) stop 'Error allocating seismo1/2'

  seismo1(:,:) = seismoRef(:,1:entries)
  seismo2(:,:) = seismo(:,1:entries)

  ! return time lag calculated analytically
  call getAnalyticalTimelagSeismos(t_lag,seismo1,seismo2,entries,timedelta)

  end subroutine


!-----------------------------------------------------------------------
  subroutine getAnalyticalTimelagSeismos(t_lag,seismoRef,seismo,seismoLength,timedelta)
!-----------------------------------------------------------------------
! determines the time lag between two seismograms analytically
! (both seismograms must have the same time step size)
!
! definition: seismoRef shifted to the right of seismo, then timelag is positiv (i.e. waves arrive later in seismo1)
!
! input:
!       t_lag                        - time lag
!       seismoRef,seismo             - seismograms (array 2 x seismoLength, first index is time, second displacement)
!       seismoLength                 - seismograms length (number of entries in seismogram)
!
! remember:
!       FILTER_SEISMOGRAMS           - use filter before determining time lag
! returns: t_lag timelag in seconds ( uses dt )
  use precisions
  use verbosity, only: beVerbose,fileOutput
  use propagationStartup, only: datadirectory
  implicit none
  real(WP),intent(out):: t_lag
  integer,intent(in):: seismoLength
  real(WP),intent(inout):: seismoRef(2,seismoLength),seismo(2,seismoLength)
  real(WP),intent(in):: timedelta
  ! local parameters
  integer:: i
  real(WP):: seismo1(seismoLength),seismo2(seismoLength),seismoPerturbed(seismoLength)
  real(WP):: integral1,integral2

  logical,parameter:: TAPER             = .true.
  logical,parameter:: TAPER_FILTERED    = .false.

  ! initialize working traces
  seismo1(:) = 0.0_WP
  seismo2(:) = 0.0_WP
  seismoPerturbed(:) = 0.0_WP

  ! debug output
  if (fileOutput) then
    print *,'  printing to file: ',trim(datadirectory)//'TimelagAnalytic_read.dat'
    open(IOUT,file=trim(datadirectory)//'TimelagAnalytic_read.dat')
    do i = 1,seismoLength
      write(IOUT,*) seismo(1,i),seismo(2,i),seismoRef(2,i)
    enddo
    close(IOUT)
  endif

  ! tapering ends
  if (TAPER) then
    call taperSeismogram(seismo,seismoLength,seismoLength,beVerbose)
    call taperSeismogram(seismoRef,seismoLength,seismoLength,beVerbose)
  endif

  ! filter
  if (FILTER_SEISMOGRAMS) then
    call dofilterSeismogram(seismoRef,seismoLength)
    call dofilterSeismogram(seismo,seismoLength)

    ! apply hanning window  to smooth seismograms ends
    if (TAPER_FILTERED) then
      call taperSeismogram(seismo,seismoLength,seismoLength,beVerbose)
      call taperSeismogram(seismoRef,seismoLength,seismoLength,beVerbose)
    endif

    ! debug output
    if (fileOutput) then
      print *,'  printing to file: ',trim(datadirectory)//'TimelagAnalytic_filter.dat'
      open(IOUT,file=trim(datadirectory)//'TimelagAnalytic_filter.dat')
      do i = 1,seismoLength
        write(IOUT,*) seismo(1,i),seismo(2,i),seismoRef(2,i)
      enddo
      close(IOUT)
    endif
  endif

  ! analytically derived expression for timelag
  ! get perturbed seismogram (only perturbations/differences between reference and delta seismogram)
  call getPerturbedSeismo(seismoRef(2,:),seismo(2,:),seismoPerturbed,seismoLength)

  ! debug output
  if (fileOutput) then
    print *,'  printing to file: ',trim(datadirectory)//'TimelagAnalytic_perturbed.dat'
    open(IOUT,file=trim(datadirectory)//'TimelagAnalytic_perturbed.dat')
    do i = 1,seismoLength
      write(IOUT,*) seismo(1,i),seismoPerturbed(i)
    enddo
    close(IOUT)
  endif

  ! get first time derivative of reference seismogram
  if (beVerbose) print *,'    first deriv.. '
  call getFirstDerivative(seismoRef(2,:),seismo1,timedelta,seismoLength)

  ! debug output
  if (fileOutput) then
    print *,'  printing to file: ',trim(datadirectory)//'TimelagAnalytic_first.dat'
    open(IOUT,file=trim(datadirectory)//'TimelagAnalytic_first.dat')
    do i = 1,seismoLength
      write(IOUT,*) seismo(1,i),seismo1(i)
    enddo
    close(IOUT)
  endif

  ! get second time derivative of reference seismogram
  if (beVerbose) print *,'    second deriv.. '
  call getSecondDerivative(seismoRef(2,:),seismo2,timedelta,seismoLength)

  ! debug output
  if (fileOutput) then
    print *,'  printing to file: ',trim(datadirectory)//'TimelagAnalytic_second.dat'
    open(IOUT,file=trim(datadirectory)//'TimelagAnalytic_second.dat')
    do i = 1,seismoLength
      write(IOUT,*) seismo(1,i),seismo2(i)
    enddo
    close(IOUT)
  endif

  ! calculate integrals
  if (beVerbose) print *,'    first integral.. '
  call getIntegral(seismo1,seismoPerturbed,integral1,timedelta,seismoLength)

  if (beVerbose) print *,'    second integral.. '
  call getIntegral(seismo2,seismoRef(2,:),integral2,timedelta,seismoLength)

  ! time lag
  t_lag = integral1/integral2

  end subroutine


!-----------------------------------------------------------------------
  subroutine getMinimized(t_lag,amplification,fileDelta,fileReference,startingTime,endingTime)
!-----------------------------------------------------------------------
! determines the time lag and amplification by nonlinear minimization (downhill simplex)
! between two seismograms within a certain time window and given
! their filenames
! (both seismograms must have the same time step size)
!
! definition: fileDelta shifted to the right of fileReference,
!                  then timelag is positiv (i.e. waves arrive later in fileDelta)
!
! input:
!       t_lag                                 - time lag
!       fileDelta,fileReference      - file names of seismograms
!       startingTime                     - window start time of fourier transformation for correlation
!
! returns: t_lag timelag in seconds
  use verbosity; use filterType
  implicit none
  character(len=128),intent(in):: fileDelta,fileReference
  real(WP),intent(in):: startingTime,endingTime
  real(WP),intent(out):: t_lag,amplification
  ! local parameters
  integer:: entries,ier
  real(WP):: lastseismotime,timedelta
  real(WP):: seismo(2,WindowSIZE),seismoRef(2,WindowSIZE)
  real(WP),allocatable,dimension(:,:):: seismo1,seismo2
  integer:: seismoLength

  !initialize
  seismo(:,:) = 0.0_WP
  seismoRef(:,:) = 0.0_WP

  ! read in from startingTime
  call readSeismoTaped(seismoRef,WindowSIZE,fileReference,startingTime,endingTime,timedelta,entries,lastseismotime)
  call readSeismoTaped(seismo,WindowSIZE,fileDelta,startingTime,endingTime,timedelta,entries,lastseismotime)

  allocate(seismo1(2,entries),seismo2(2,entries),stat=ier)
  if (ier /= 0) stop 'Error allocating seismo1,seismo2'

  ! sets convention:
  ! seismo1 is the heterogeneous trace
  ! seismo2 is the homogeneous trace
  seismoLength = entries
  seismo1(:,:) = seismo(:,1:entries)
  seismo2(:,:) = seismoRef(:,1:entries)

  call getTimelagSimplex(t_lag,amplification,seismo1,seismo2,seismoLength)

  deallocate(seismo1,seismo2)

  end subroutine


!-----------------------------------------------------------------------
  subroutine getTimelagSimplex(t_lag,amplification,seismo1,seismo2,seismoLength)
!-----------------------------------------------------------------------
! determines the time lag and amplification by nonlinear minimization (downhill simplex)
! between two seismograms within a certain time window and given
! their filenames
! (both seismograms must have the same time step size)
!
! input:
!       t_lag                                 - time lag
!       seismo1,seismo2             - seismograms to compare
!
! sets convention:
! seismo1 is the heterogeneous trace
! seismo2 is the homogeneous trace
!
! returns: t_lag (timelag in seconds), amplification (factor to scale homogeneous seismogram with)
  use precisions
  use verbosity; use filterType
  use propagationStartup, only: datadirectory,dt
  implicit none
  real(WP),intent(out):: t_lag,amplification
  integer,intent(in):: seismoLength
  real(WP),intent(inout):: seismo1(2,seismoLength),seismo2(2,seismoLength)
  ! local parameters
  integer:: i
  real(WP):: seismoWindow(WindowSIZE),seismoRefWindow(WindowSIZE)
  real(WP),allocatable,dimension(:):: seismoShift

  logical,parameter:: NORMALIZE         = .false.
  logical,parameter:: STRETCH           = .false.
  logical,parameter:: TAPER             = .false.
  logical,parameter:: TAPER_FILTERED    = .false.

  ! console
  if (beVerbose) then
    print *,'  getting timelag seismos by downhill simplex minimization...'
    print *
    print *,'  perturbed seismogram:'
    print *,'    first entry: ',seismo1(1,1),seismo1(2,1)
    print *,'    last entry : ',seismo1(1,seismoLength),seismo1(2,seismoLength)
    print *,'    entry read : ',seismoLength
    print *,'  reference seismogram:'
    print *,'    first entry: ',seismo2(1,1),seismo2(2,1)
    print *,'    last entry : ',seismo2(1,seismoLength),seismo2(2,seismoLength)
    print *,'    entry read : ',seismoLength
    print *,'  time step: dt = ',dt
    print *
  endif

  ! checks if dt is set
  if (dt == 0.0) call stopProgram('Invalid dt in getTimelagSimplex()    ')

  ! debug output
  if (fileOutput) then
    print *,'  printing to file: ',trim(datadirectory)//'TimelagSimplex_read.dat'
    open(IOUT,file=trim(datadirectory)//'TimelagSimplex_read.dat')
    write(IOUT,*) "# format: #time #trace-perturbed #trace-reference"
    do i = 1,seismoLength
      write(IOUT,*) seismo1(1,i),seismo1(2,i),seismo2(2,i)
    enddo
    close(IOUT)
  endif

  ! tapering ends
  if (TAPER) then
    if (beVerbose) print *,'  taper ends'
    call taperSeismogram(seismo1,seismoLength,seismoLength,beVerbose)
    call taperSeismogram(seismo2,seismoLength,seismoLength,beVerbose)
  endif

  ! filtering necessary
  if (FILTER_SEISMOGRAMS) then
    if (beVerbose) print *,'  filter seismograms'
    call dofilterSeismogram(seismo1,seismoLength)
    call dofilterSeismogram(seismo2,seismoLength)

    ! apply hanning window to smooth seismograms ends
    if (TAPER_FILTERED) then
      if (beVerbose) print *,'  taper seismogram ends'
      call taperSeismogram(seismo1,seismoLength,seismoLength,beVerbose)
      call taperSeismogram(seismo2,seismoLength,seismoLength,beVerbose)
    endif
  endif

  !initialize
  seismoWindow(:) = 0.0_WP
  seismoRefWindow(:) = 0.0_WP
  ! window zeros are added
  if (WindowSIZE < seismoLength) stop 'Error windowsize/seismoLength'
  seismoWindow(1:seismoLength) = seismo1(2,1:seismoLength)
  seismoRefWindow(1:seismoLength) = seismo2(2,1:seismoLength)

  ! normalize seismograms
  if (NORMALIZE) then
    print *,'  normalize traces'
    call normalizeArray(seismo1(2,:),seismoLength)
    call normalizeArray(seismo2(2,:),seismoLength)
  endif

  ! debug output
  if (fileOutput) then
    print *,'  printing to file: ',trim(datadirectory)//'TimelagSimplex_window.dat'
    open(IOUT,file=trim(datadirectory)//'TimelagSimplex_window.dat')
    write(IOUT,*) "# format: #time #windowed-trace #windowed-reference-trace"
    do i = 1,seismoLength
      write(IOUT,*) seismo2(1,i),seismo1(2,i),seismo2(2,i)
    enddo
    close(IOUT)
  endif

  ! get minimized function (by numerical recipes)
  call minimizeTraces(t_lag,amplification,seismo1,seismo2,seismoLength)

  if (fileOutput) then
    print *,'  printing to file: ',trim(datadirectory)//'TimelagSimplex_shift.dat'
    allocate(seismoShift(seismoLength))
    call sumSeismogram(t_lag,amplification,seismo2,seismo1,seismoShift,seismoLength)
    open(IOUT,file=trim(datadirectory)//'TimelagSimplex_shift.dat')
    write(IOUT,*) '# simplex shift'
    write(IOUT,*) "# format: #time #seismo-shift"
    do i = 1,seismoLength
      write(IOUT,*) seismo1(1,i),seismoShift(i)
    enddo
    close(IOUT)
  endif

  !only be verbose once
  if (beVerbose) then
    beVerbose = .false.
    fileOutput = .false.
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine minimizeTraces(t_lag,amplification,seismo1,seismo2,seismoLength)
!-----------------------------------------------------------------------
! finds t_lag and amplification to match the seismo1 with seismo2
! moves the homogeneous trace seismo2 towards the heterogeneous trace seismo1
!
! inputs:
!     t_lag, amplification       - time lag and amplitude anomaly
!     seismo1,seismo2           - heterogeneous and reference trace
!     seismoLength                - size/length of traces
!
! returns: t_lag and amplification
  use precisions; use verbosity
  use minimize_trace, only: costfunction,seismoHom,seismoHet,trace_length
  implicit none
  integer,intent(in):: seismoLength
  real(WP),dimension(2,seismoLength):: seismo1,seismo2
  real(WP),intent(out):: t_lag,amplification
  ! local parameters
  integer:: i,ier
  integer:: iter
  real,dimension(3):: y
  real,dimension(3,2):: p
  real,dimension(2):: p_x
  real,parameter:: ftol = 1.e-3  ! tolerance
  interface
    subroutine minimize_downhill_simplex(p, y, ftol, func, iter)
      implicit none
      integer, intent(out) :: iter
      real, intent(in) :: ftol
      real, dimension(3), intent(inout) :: y
      real, dimension(3,2), intent(inout) :: p
      interface
        function func(x)
          implicit none
          real, dimension(:), intent(in) :: x
          real :: func
        end function func
      end interface
    end subroutine
  end interface
  logical,parameter:: DEBUG = .false.

  ! set for costfunction
  if (trace_length /= seismoLength) then
    if (allocated(seismoHom)) deallocate(seismoHom,seismoHet)
    allocate(seismoHom(2,seismoLength),seismoHet(2,seismoLength),stat=ier)
    if (ier /= 0) stop 'Error allocating minimizing traces'
  endif

  trace_length = seismoLength
  seismoHom(:,:) = seismo2(:,:)
  seismoHet(:,:) = seismo1(:,:)

  ! simplex startup values
  p(1,1) = 0.0  ! t_lag
  p(1,2) = 0.1  ! amp
  p(2,1) = -t_lag ! t_lag
  p(2,2) = 1.0  ! amp
  p(3,1) = t_lag ! t_lag
  p(3,2) = 1.0  ! amp
  do i = 1,3
    p_x(:) = p(i,:)
    y(i) = costfunction(p_x)
  enddo

  ! downhill simplex algorithm
  call minimize_downhill_simplex(p,y,ftol,costfunction,iter)

  if (DEBUG) then
    print *,'debug: minimize by downhill simplex: iterations = ',iter
    do i = 1,3
      print *,'    array ',i
      print *,'    rms                 :',y(i)
      print *,'    t_lag/amplification :',p(i,1),p(i,2)
    enddo
    print *
  endif

  ! set our two parameters:
  ! time shift ( positiv for seismo1 later arriving than seismo2)
  ! and amplification factor
  t_lag = p(1,1)
  amplification = p(1,2)

  end subroutine


!-----------------------------------------------------------------------
  subroutine sumSeismogram(arrival,amplification,seismo,seismoRef,seismoOut,length)
!-----------------------------------------------------------------------
! interpolates the values of the shifted (arrival) and amplified seismogram seismo
! at the times of seismoRef
!
! returns: seismoOut the interpolated values
  use precisions; use verbosity
  use module_spline, only: spline,splint
  use propagationStartup, only: dt
  implicit none
  real(WP),intent(in):: arrival,amplification
  integer,intent(in):: length
  real(WP),intent(in):: seismo(2,length),seismoRef(2,length)
  real(WP),intent(out):: seismoOut(length)
  ! local parameters
  double precision:: ya(length), xa(length),y2(length),yp1,ypn,loc,yb,dta
  double precision,parameter:: EPS = 0.01
  integer:: i

  ! initializes
  seismoOut(:) = 0.0_WP

  ! set spline arrays
  do i = 1,length
     xa(i) = seismo(1,i) + arrival
     ya(i) = seismo(2,i) * amplification
  enddo

  dta = dt

  ! check
  if (dta == 0.0) call stopProgram('Invalid dt in sumSeismogram()    ')

  yp1 = (ya(2) - ya(1)) / dta
  ypn = (ya(length) - ya(length-1)) / dta

  call spline(xa,ya,yp1,ypn,y2)

  ! interpolate values
  do i = 1,length
    loc = seismoRef(1,i)
    if (loc <= xa(length)+EPS .and. loc >= xa(1)-EPS) then
      yb = splint(xa,ya,y2,loc)
    else
      yb = 0.0
    endif

    seismoOut(i) = yb
  enddo

  end subroutine


!-----------------------------------------------------------------------
  subroutine minimize_downhill_simplex(p, y, ftol, func, iter)
!-----------------------------------------------------------------------
! Minimization of a function using the downhill simplex method of
! Nelder and Mead in N dimensions.
!
! Input:
!   p     - (N+1)xN matrix of simplex vertices.
!   y     - function values at simplex vertices.
!   ftol  - Convergence tolerance.
!   func  - function to minimize.
!
! Output:
!   p     - Updated simplex vertices.
!   y     - Updated function values at vertices.
!   iter  - Number of function evaluations taken.
  implicit none
  ! Input/Output parameters
  integer, intent(out) :: iter
  real, intent(in) :: ftol
  real, dimension(3), intent(inout) :: y
  real, dimension(3,2), intent(inout) :: p
  ! function interface
  interface
    function func(x)
      implicit none
      real, dimension(:), intent(in) :: x
      real :: func
    end function func
  end interface
  ! local parameters
  integer, parameter :: ITMAX = 500   ! Max iterations (reduced from 5000)
  real, parameter :: TINY = 1.0e-10   ! Small value to prevent division errors
  ! Internal variables
  integer:: ihi, ndim
  real, dimension(size(p, 2)) :: psum

  ! Initialize and run the private routine
  call minimize_func()

contains

  !-----------------------------------------------------------------------
  subroutine minimize_func()
  !-----------------------------------------------------------------------
  ! Internal routine that performs the actual minimization process.
  implicit none
  ! local parameters
  integer :: i, ilo, inhi, iloc(1)
  real :: rtol, ysave, ytry, ytmp

  ! Ensure proper array sizes
  ndim = size(p,2)
  if (ndim /= size(p,1)-1 .or. ndim /= size(y)-1) call stopProgram('minimize_func invalid ndim    ')

  iter = 0
  psum(:) = sum(p(:,:), dim=1)

  do  ! Main iteration loop
    ! Identify best (ilo), worst (ihi), and second worst (inhi) vertices
    iloc = minloc(y(:))
    ilo = iloc(1)

    iloc = maxloc(y(:))
    ihi = iloc(1)

    ytmp = y(ihi)
    y(ihi) = y(ilo)

    iloc = maxloc(y(:))
    inhi = iloc(1)

    y(ihi) = ytmp

    ! Convergence test: fractional range of function values
    rtol = 2.0 * abs(y(ihi) - y(ilo)) / (abs(y(ihi)) + abs(y(ilo)) + TINY)

    !debug
    !print *,'debug: minimize_func: iter = ',iter,' rtol/ftol = ',rtol,ftol,' yhi/ylo = ',y(ihi),y(ilo),ihi,ilo

    ! If convergence criterion met, exit
    if (rtol < ftol) then
      call swap_single(y(1), y(ilo))
      call swap_vectors(p(1,:), p(ilo,:))
      return
    endif

    ! Check for iteration limit
    if (iter >= ITMAX) call stopProgram('minimize_func: ITMAX exceeded    ')

    ! Reflect the simplex across the worst vertex
    ytry = extrapolate(-1.0)
    iter = iter + 1

    ! Check if new point is better than the best
    if (ytry <= y(ilo)) then
      ytry = extrapolate(2.0)  ! Expand simplex
      iter = iter + 1
    else if (ytry >= y(inhi)) then
      ! Contraction step if reflected point is still poor
      ysave = y(ihi)
      ytry = extrapolate(0.5)
      iter = iter + 1

      if (ytry >= ysave) then
        ! If contraction fails, shrink simplex around best vertex
        p(:,:) = 0.5 * (p(:,:) + spread(p(ilo,:), 1, size(p,1)))

        do i = 1, ndim+1
          if (i /= ilo) y(i) = func(p(i,:))
        enddo

        iter = iter + ndim
        psum(:) = sum(p(:,:), dim=1)
      endif
    endif
  enddo  ! End iteration loop

  end subroutine


  !-----------------------------------------------------------------------
  subroutine swap_single(a, b)
  !-----------------------------------------------------------------------
  ! Swaps two single entries
  !
  ! Arguments:
  !   a - First entry
  !   b - Second entry
  implicit none
  real, intent(inout) :: a, b  ! Input/output
  ! local parameters
  real:: temp                  ! Temporary storage for swapping

  ! Swap the contents of the two vectors
  temp = a
  a = b
  b = temp

  end subroutine


  !-----------------------------------------------------------------------
  subroutine swap_vectors(vector_a, vector_b)
  !-----------------------------------------------------------------------
  ! Swaps the contents of two vectors.
  !
  ! Arguments:
  !   vector_a - First complex vector (modified in place).
  !   vector_b - Second complex vector (modified in place).
  implicit none
  real, dimension(:), intent(inout) :: vector_a, vector_b  ! Input/output complex vectors
  ! local parameters
  real, dimension(size(vector_a)) :: temp_vector           ! Temporary storage for swapping

  ! Swap the contents of the two vectors
  temp_vector(:) = vector_a(:)
  vector_a(:) = vector_b(:)
  vector_b(:) = temp_vector(:)

  end subroutine


  !-----------------------------------------------------------------------
  real function extrapolate(fac)
  !-----------------------------------------------------------------------
  ! Extrapolates a new point through the face opposite the worst point
  ! and evaluates its function value.
  !
  ! Input:
  !   fac - Factor determining the nature of extrapolation:
  !         -1.0 for reflection,
  !          2.0 for expansion,
  !          0.5 for contraction.
  !
  ! Output:
  !   extrapolate - function value at new extrapolated point.
  implicit none
  real, intent(in) :: fac
  ! local parameters
  real :: fac1, fac2, ytry
  real, dimension(size(p,2)) :: ptry

  ! Compute extrapolated point
  fac1 = (1.0 - fac) / ndim
  fac2 = fac1 - fac
  ptry(:) = psum(:) * fac1 - p(ihi,:) * fac2

  ! Evaluate function at trial point
  ytry = func(ptry)

  ! If better than worst point, replace it
  if (ytry < y(ihi)) then
    y(ihi) = ytry
    psum(:) = psum(:) - p(ihi,:) + ptry(:)
    p(ihi,:) = ptry(:)
  endif

  extrapolate = ytry

  return
  end function

  end subroutine minimize_downhill_simplex

