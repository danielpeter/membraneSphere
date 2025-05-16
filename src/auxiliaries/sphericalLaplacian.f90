!=====================================================================
!
!       m e m b r a n e S p h e r e
!       --------------------------------------------------
!      (c) 2025
!
!=====================================================================
!
! Proper acknowledgement shall be made to the authors of the software in publications and
! presentations resulting from the use of this software:
!
! Peter, D., C. Tape, L. Boschi and J. H. Woodhouse, 2007. Surface wave tomography:
! global membrane waves and adjoint methods, Geophys. J. Int., , 171: p. 1098-1117.

!-----------------------------------------------------------------------
  program sphericalLaplacian
!-----------------------------------------------------------------------
!  sphericalLaplacian calculates the Laplacian on a sphere for a
!  hexagonal grid point given by the data in 'DvvertN.dat' and  'DvfaceN.dat'
  use cells
  use precisions, only: PRECALCULATED_CELLS
  use verbosity, only: VERBOSE
  implicit none

  ! console output
  print *,'-----------------------------------------------------------------'
  print *,"sphericalLaplacian"
  print *, 'Laplacian for hexagonal cells on a sphere'
  print *,'-----------------------------------------------------------------'

  ! safety check
  if (.not. STATION_CORRECTION) stop 'Please re-compile this executable with STATION_CORRECTION turned on in commonModules.f90'

  ! user input
  print *, 'Number of Subdivisions files to use: ( >= 0)'
  read(*,*) subdivisions

  print *
  print *,'subdivisions: ',subdivisions
  print *

  ! read vertices values from files
  print *,'reading in data files for',subdivisions,' subdivisions'
  ! be verbose
  VERBOSE = .true.

  ! read in triangular grid as well
  call allocateMesh()
  call readData(VERBOSE)
  call allocateData()
  if (PRECALCULATED_CELLS) call readPrecalculated()

  ! handles our model
  call manageComparison()

  end program

!-----------------------------------------------------------------------
  subroutine manageComparison()
!-----------------------------------------------------------------------
! manages the displacements at each grid point
  use cells
  implicit none
  integer:: ier
  real(WP),dimension(:),allocatable:: harmonicLaplacian,laplacian

  ! allocates arrays
  allocate(harmonicLaplacian(numVertices), &
           laplacian(numVertices),stat=ier)
  if (ier /= 0) stop 'Error allocating laplacian arrays'
  harmonicLaplacian(:) = 0.0_WP
  laplacian(:) = 0.0_WP

  ! initial displacements are real spherical harmonic function values
  ! (see: Tape (2003) )
  print *,'setup initial distribution...'
  call initialDistribution()

  ! numerical laplacian
  print *,'calculate numerical Laplacian...'
  call numericLaplacian(laplacian)

  ! analytical laplacian
  print *,'calculate analytical Laplacian...'
  call analyticLaplacian(harmonicLaplacian)

  ! calculate the error in the numerical approximation
  print *,'calculate errors...'
  call differenceError(harmonicLaplacian,laplacian)

  ! iteration time step
  !print *,'iterate numerically through time...'
  !call numericIterate()

  ! cleanup
  deallocate(harmonicLaplacian,laplacian)

  end subroutine


!-----------------------------------------------------------------------
  subroutine initialDistribution()
!-----------------------------------------------------------------------
! determines the starting distributions of displacements
! on all grid points
  use cells; use displacements
  use precisions, only: DEGREE_L,DEGREE_M
  use constants, only: IOUT
  implicit none
  integer::  n,i,k,ier
  real(WP):: vector(3), displace
  character(len=2):: div,degL, degM
  character(len=14):: append
  double precision,external:: sphericalHarmonics

  print *,'  using spherical harmonic:'
  print *,'    degree L = ',DEGREE_L
  print *,'    degree M = ',DEGREE_M

  ! start with a spherical harmonic value
  do n = 1,numVertices
    displacement(n) = sphericalHarmonics(n)
  enddo
  
  ! we use the real spherical harmonic value
  write(div,'(i0)') subdivisions
  write(degL,'(i0)') DEGREE_L
  write(degM,'(i0)') DEGREE_M

  append = trim(div)//'.L'//trim(degL)//'M'//trim(degM)//'.dat'

  ! file output
  print *,'  writing harmonics...'

  open(IOUT,file='OUTPUT/harmonics'//trim(append),iostat=ier)
  if (ier /= 0) stop 'Error opening harmonics file'
  write(IOUT,'(a,i2,a,i2)') '# Spherical harmonics - degree L = ',DEGREE_L,' / degree M = ',DEGREE_M
  write(IOUT,'(a)') '#'
  write(IOUT,'(a)') '#grid_vertex_value'
  do n = 1,numVertices
    write(IOUT,'(f12.4)') displacement(n)
  enddo
  close(IOUT)
  print *,'  written to: OUTPUT/sphericalharmonics'//trim(append)

  ! print out new file
  print *,'  writing face displacements...'

  open(IOUT,file='OUTPUT/harmonicsfaces'//trim(append),iostat=ier)
  if (ier /= 0) stop 'Error opening harmonicsfaces file'
  do n = 1, numTriangleFaces
    write(IOUT,*) '#',n
    do i = 1,3
      vector(:) = vertices(cellTriangleFace(n,i),:)
      displace = displacement(cellTriangleFace(n,i))
      displace = displace*displace
      call vectorscale(vector, displace, vector)
      write(IOUT,'(3f16.4)') (vector(k),k=1,3)
    enddo !i
    write(IOUT,*)
  enddo !n
  close(IOUT)
  print *,'  written to: OUTPUT/harmonicsfaces'//trim(append)
  print *

  end subroutine

!-----------------------------------------------------------------------
  subroutine numericLaplacian(laplacian)
!-----------------------------------------------------------------------
! calculates the laplacian of each grid point
  use cells
  use precisions, only: PRECALCULATED_CELLS
  implicit none
  real(WP),intent(out):: laplacian(numVertices)
  ! local parameters
  integer:: n
  real(WP):: D2
  real(WP),external:: precalc_discreteLaplacian,discreteLaplacian

  laplacian(:) = 0.0_WP

  ! propagate only corresponding vertices
  do n = 1, numVertices
    ! spherical Laplacian
    if (PRECALCULATED_CELLS) then
      D2 = precalc_discreteLaplacian(n) !uses displacement(..) field
    else
      D2 = discreteLaplacian(n)
    endif
    laplacian(n) = D2
  enddo

  end subroutine

!-----------------------------------------------------------------------
  subroutine analyticLaplacian(harmonicLaplacian)
!-----------------------------------------------------------------------
! calculates the laplacian of each grid point
  use cells; use displacements
  use precisions, only: DEGREE_L
  implicit none
  real(WP),intent(out):: harmonicLaplacian(numVertices)
  ! local parameters
  real(WP):: laplace
  integer:: n

  harmonicLaplacian(:) = 0.0_WP

  do n = 1, numVertices
    !calculate the analytical laplacian for the spherical harmonic
    laplace = - DEGREE_L * (DEGREE_L + 1) * displacement(n) / (EarthRadius**2)
    harmonicLaplacian(n) = laplace
  enddo

  end subroutine

!-----------------------------------------------------------------------
  subroutine differenceError(harmonicLaplacian,laplacian)
!-----------------------------------------------------------------------
! calculates the error between numerical and analytical
! laplacian
!
! uses: carl tape thesis2003, chap 4, 4.15
  use cells
  use precisions, only: DEGREE_L,DEGREE_M
  use constants, only: IOUT
  implicit none
  real(WP),intent(in):: harmonicLaplacian(numVertices),laplacian(numVertices)
  ! local parameters
  real(WP),dimension(:),allocatable:: numericalError
  integer:: n,k,ier
  real(WP):: maxAnalytical
  integer:: maxAnalyticalIndex
  real(WP):: maxNumerical
  integer:: maxNumericalIndex
  real(WP):: diff, diffMax
  integer:: diffMaxIndex
  double precision:: diffmean
  character(len=2):: div,degL, degM
  character(len=14):: append

  ! allocates array
  allocate(numericalError(numVertices),stat=ier)
  if (ier /= 0) stop 'Error allocating numericalError array'
  numericalError(:) = 0.0_WP

  ! find the maximum single value for the analytical laplacian on the sphere
  maxAnalytical = 0.0_WP
  maxNumerical = 0.0_WP
  do n = 1, numVertices
    !maximum analytical value
    if (abs(harmonicLaplacian(n)) >= abs(maxAnalytical)) then
      maxAnalytical = harmonicLaplacian(n)
      maxAnalyticalIndex = n
    endif
    !maximum numerical value
    if (abs(laplacian(n)) >= abs(maxNumerical)) then
      maxNumerical = laplacian(n)
      maxNumericalIndex = n
    endif
  enddo

  print *
  print *,'numerical approximation errors'
  print *,'  maximum analytical Laplacian    : ',maxAnalytical,'at vertex',maxAnalyticalIndex
  print *,'  maximum numerical  Laplacian    : ',maxNumerical,'at vertex',maxNumericalIndex

  ! check
  if (maxAnalytical == 0.0_WP) stop 'Error invalid maximum for analytical Laplacian, must be non-zero'

  ! calculate the error between numerical and analytical laplacian
  diffMax = 0.0_WP
  do n = 1, numVertices
    diff = (laplacian(n) - harmonicLaplacian(n)) / maxAnalytical
    numericalError(n) = diff

    !absolute
    !numericalError(n) = abs(diff)

    ! stats
    if (abs(diff) >= abs(diffMax)) then
      diffMax = diff
      diffMaxIndex = n
    endif
  enddo
  print *,'  maximum difference              : ',diffMax,'at vertex',diffMaxIndex

  ! calculate the mean error
  ! (compare with formula: Tape,2003,chap.4, 4.16)
  diffmean = 0.d0
  do n = 1, numVertices
    diffmean = diffmean + abs(numericalError(n))
  enddo
  diffmean = diffmean / dble(numVertices)

  print *
  print *,'  mean (absolute) difference error:',sngl(diffmean),'(D2error)'
  print *

  ! file output
  print *,'writing numerical approximation errors...'

  ! file extension use
  write(div,'(i0)') subdivisions
  write(degL,'(i0)') DEGREE_L
  write(degM,'(i0)') DEGREE_M

  append = trim(div)//'.L'//trim(degL)//'M'//trim(degM)//'.dat'

  open(IOUT,file='OUTPUT/numericalError'//trim(append),iostat=ier)
  if (ier /= 0) stop 'Error opening numericalError file'
  write(IOUT,'(a)') '# numerical approximation errors'
  write(IOUT,'(a)') '#'
  write(IOUT,'(a,i2)') '# grid level                        : ',subdivisions
  write(IOUT,'(a,i2,a,i2)') '# spherical harmonics (degree L & M): ',DEGREE_L,' / ',DEGREE_M
  write(IOUT,'(a)') '#'
  write(IOUT,'(a,es18.6)') '# mean (absolute) difference error  : ',sngl(diffmean)
  write(IOUT,'(a)') '#'
  write(IOUT,'(a)') '#vertex_x #vertex_y #vertex_z #numerical_error'
  do n = 1, numVertices
    !output
    write(IOUT,'(3f12.8,2x,e18.6)') ( vertices(n,k),k=1,3),numericalError(n)
  enddo
  close(IOUT)
  print *,'  written to: OUTPUT/numericalError'//trim(append)
  print *

  ! output main difference locations
  !do n = 1,numVertices
  ! if (abs((abs(numericalError(n))-abs(diffMax))) < 0.00000001) then
  !    print *,'  maximum error at point ',n, numericalError(n)
  !  endif
  !enddo
  !print *,'  error at(112): ', numericalError(112)

  end subroutine

!-----------------------------------------------------------------------
  subroutine numericIterate()
!-----------------------------------------------------------------------
! iterates through one time step with numerically
! calculation of the new displacement
!
! uses: carl tape thesis2003, chap 5, (5.7)
  use precisions; use cells
  use displacements
  use propagationStartup, only: cphaseRef,dt
  use parameter_inputs, only: WidthParameterMu
  implicit none
  real(WP):: displace_t, displace_tplus1, displace_tminus1
  real(WP):: forcing(numVertices)
  real(WP):: velocity
  real(WP),external:: discreteLaplacian
  real(WP):: colat, long, mu2
  integer::n,k,timestep
  character(len=4):: timestr
  real(WP):: vec(3),vecdisplace(3)
  real(WP):: time, exact, u_shape
  logical:: doExact
  real(WP),parameter:: sqrt2 = 1.414213562373095
  !drawing
  !real(WP):: tmpArray(MAXVERTICES*4)
  !external drawInit !$PRAGMA C( drawInit )
  !external drawArray !$PRAGMA C( drawArray )

  ! set initial conditions
  call readParameters()

  mu2 = WidthParameterMu**2
  velocity = cphaseRef ![km/s]

  ! dt
  call determineTimeStep()

  !debug
  print *,"time step dt =" ,dt
  timestep = 0
  time = timestep*dt

  write(timestr,'(i4.4)') timestep

  doExact = .true.

  if (doExact) then
    open(20,file='OUTPUT/simulation.'//timestr//'.error.dat')
    open(21,file='OUTPUT/simulation.'//timestr//'.exact.dat')
  endif
  do n = 1, numVertices
    !get spherical coordinates of vertex
    call getSphericalCoord(n, colat, long)

    !set initial conditions
    newdisplacement(n) = exp( - colat**2/(2.d0*mu2) ) / mu2
    displacement(n) = newdisplacement(n)
    displacement_old(n) = displacement(n)
    forcing(n) = 0.0

    !check if exact solution is right for time 0
    if (doExact) then
      exact = u_shape(colat,long, time)
      !print *,newdisplacement(n), exact
      if (abs(newdisplacement(n)-exact) > 0.001) then
        print *,"vertex ", n
        print *,"displacement ",newdisplacement(n)
        print *,"exact ", exact
      endif

      !write initial error output
      write(20,'(3f14.2,f14.10)') (vertices(n,k)*EarthRadius,k=1,3),newdisplacement(n)-exact
      write(21,'(3f14.2,f14.8)') (vertices(n,k)*EarthRadius,k=1,3),exact
    endif
  enddo
  if (doExact) then
    close(20)
    close(21)
  endif

  !write initial output
  open(101,file='OUTPUT/simulation.'//timestr//'.dat')
  open(102,file='OUTPUT/simulation.'//timestr//'.displace.dat')
  do n = 1, numVertices
    do k = 1,3
      vec(k) = vertices(n,k)
    enddo
    call vectorscale(vec,EarthRadius+newdisplacement(n),vecdisplace)
    call vectorscale(vec, EarthRadius, vec)
    write(101,'(3f14.2,f14.8)')(vec(k),k=1,3),newdisplacement(n)
    write(102,'(3f14.2)')(vecdisplace(k),k=1,3)

    !temporary array for drawing
    !tmpArray((n-1)*4+1) = vertices(n,1)
    !tmpArray((n-1)*4+2) = vertices(n,2)
    !tmpArray((n-1)*4+3) = vertices(n,3)
    !tmpArray((n-1)*4+4) = newdisplacement(n)
  enddo
  close(101)
  close(102)

  !call external c-function
  !do n=1, numVertices
  !  tmpArray(n)= vertices(n,1)
  !  tmpArray(n+1)= vertices(n,2)
  !  tmpArray(n+2)= vertices(n,3)
  !  tmpArray(n+3)= newdisplacement(n)
  !enddo
  !debug
  !print *,'initial temporary array for drawing ',(tmpArray(k),k=1,32)
  !call drawInit()
  !call drawArray(tmpArray, numVertices)

  !print *,'do time iteration  from 1 to 20'
  !do some steps forward in time
  do timestep = 1, 1000
    ! iterate one timestep
    !call numericLaplacian()
    do n = 1, numVertices
      !get old displacement
      displacement_old(n) = displacement(n)
      displace_tminus1 = displacement_old(n)

      !get actual displacement
      displacement(n) = newdisplacement(n)
      displace_t =displacement(n)

      ! calculate new displacement in time
      displace_tplus1 = 2.0*displace_t-displace_tminus1+(velocity**2)*(dt**2)*(discreteLaplacian(n)+forcing(n))
      newdisplacement(n) = displace_tplus1

      !debug
      !if (mod(timestep,50) == 0) then
      !  print *,'u_t-1 ', displace_tminus1
      !  print *,'u_t ', displace_t
      !  print *,'d2 ', discreteLaplacian(n)
      !  print *,'vertex ',n,': u_t+1 = ',displace_tplus1
      !endif
    enddo

    !write output
    if (mod(timestep,100) == 0) then
      write(timestr,'(i4.4)') timestep
      open(101,file='OUTPUT/simulation.'//timestr//'.dat')
      open(102,file='OUTPUT/simulation.'//timestr//'.displace.dat')
      !instant drawing uses file simulation.displace.dat
      !open(103, file='OUTPUT/simulation.displace.dat')
      do n = 1, numVertices
        !displacement_old(n) = displacement(n)
        !displacement(n) = newdisplacement(n)
        do k = 1,3
          vec(k) = vertices(n,k)
        enddo
        call vectorscale(vec, EarthRadius+newdisplacement(n),vecdisplace)
        call vectorscale(vec, EarthRadius, vec)
        write(101,'(3f14.2,f14.8)')(vec(k),k=1,3),newdisplacement(n)

        !temporary array for drawing
        !tmpArray((n-1)*4+1)= vertices(n,1)
        !tmpArray((n-1)*4+2)= vertices(n,2)
        !tmpArray((n-1)*4+3)= vertices(n,3)
        !tmpArray((n-1)*4+4)= newdisplacement(n)

        write(102,'(3f14.2)')(vecdisplace(k),k=1,3)
        !write(103,'(3f14.2)')(vecdisplace(k),k=1,3)
      enddo
      close(101)
      close(102)
      !close(103)
      ! let gnuplot draw actual displacements
      !call SYSTEM('./gnuplotdraw.sh')

      !call external c-function
      !do n=1, numVertices
      !  tmpArray(n)= vertices(n,1)
      !  tmpArray(n+1)= vertices(n,2)
      !  tmpArray(n+2)= vertices(n,3)
      !  tmpArray(n+3)= newdisplacement(n)
      !enddo
      !debug
      !print *,timestep,'temporary array for drawing '
      !print *,(tmpArray(k),k=1,4)
      !call drawArray(tmpArray, numVertices)

      !debug
      print *,"timestep ",timestep

      ! compare with analytic solution for a given initial displacement
      if (doExact) then
        print *,"comparision with analytic solution at time= ",timestep*dt
        open(21,file='OUTPUT/simulation.'//timestr//'.error.dat')
        open(22,file='OUTPUT/simulation.'//timestr//'.exact.dat')
        do n = 1, numVertices
         !get location and time of vertex
          call getSphericalCoord(n, colat, long)
          time = timestep*dt

          !calculate analytic solution
          exact = u_shape(colat, long, time)

          !output difference between numerical and exact solution for vertex
          write(21,'(3f14.2,f14.10)') (vertices(n,k)*EarthRadius, k=1,3),newdisplacement(n)-exact
          !print *,(vertices(n,k),k=1,3), (newdisplacement(n)-exact)

          write(22,'(3f14.2,f14.8)') (vertices(n,k)*EarthRadius, k=1,3),exact
        enddo
        close(21)
        close(22)
      endif !doExact
    endif

    !if (mod(timestep,100) == 0) print *,'timestep ',timestep
  enddo !timestep

  end subroutine

!-----------------------------------------------------------------------
  subroutine stopProgram(textline)
!-----------------------------------------------------------------------
! non-MPI version for this program
  use parallel
  implicit none
  character(len=*),intent(in):: textline
  ! local parameters
  integer:: endindex

  ! console output
  endindex = index(textline,"  ")
  if (endindex < 1) endindex = 128
  print *,textline(1:endindex)

  ! on linux machines : i/o unit 6 is the stdout , on SGI 101
  flush(6)

  ! stop process
  stop 'Abort - program'

  end subroutine
