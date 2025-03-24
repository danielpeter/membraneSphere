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

!-----------------------------------------------------------------------
  program propagation
!-----------------------------------------------------------------------
! iterates through one time step with numerically
! calculation of the new displacement
!
! uses: Carl Tape M.Sc. thesis, 2003, chap 5, (5.7)
  use propagationStartup
  use parallel
  use verbosity
  implicit none
  ! local parameters
  integer:: i,ier
  logical:: looping
  real(WP), external:: discreteLaplacian,precalc_discreteLaplacian
  real(WP), external:: syncWtime

  !-----------------------------------------------------------------------
  ! parameters
  ! most parameters concerning wave propagation set in file Parameter_Input
  ! (& default values in commonModules.f90)
  !-----------------------------------------------------------------------

  !print *,'Welcome to membrane waves on a spherical shell'
  !print *,'----------------------------------------------'
  !print *

  ! initialization of parameters and arrays
  call initialize()

  ! place all receiver stations
  if (manyReceivers) then
    ! useless for heterogeneous case
    if (HETEROGENEOUS) then
      print *,'multiple receiver stations and heterogeneous phase map not applicable!'
      stop 'Invalid multiple receivers for heterogeneous phase map'
    endif

    if (numofReceivers <= 0) then
      print *,'number of receiver stations invalid!',numofReceivers
      stop 'Invalid number of receiver stations'
    endif

    ! allocate seismograms for all receivers
    deallocate( seismogramReceiver )    ! only needed for single receiver setups
    allocate( receivers(numofReceivers), receiversSeismogram(numofReceivers+1,lasttimestep-firsttimestep+1),stat=ier)
    if (ier /= 0) then
      print *,'cannot allocate receivers array'
      stop 'Error allocating receivers array'
    endif

    ! assumes the source to be on north pole, the receiver station places along the same latitude
    do i = 1, numofReceivers
      ! vary longitude
      receiverLon = (i-1)*360.0/numofReceivers
      call findVertex(receiverLat,receiverLon,receivers(i))
    enddo
  endif

  if (MAIN_PROCESS .and. VERBOSE) print *,'starting time loop...'

  ! loop for each delta location
  looping = .true.
  do while( looping )
    ! benchmark
    if (MAIN_PROCESS) benchstart = syncWtime()

    ! forward simulation of membrane waves
    call forwardIteration()

    ! print seismogram to file
    call printSeismogram()

    ! move delta location
    if (DELTA) then
      if (PARALLELSEISMO) then
        call parallelFindnextLocation(looping)
      else
        call findnextLocation(looping)
      endif
    endif

    ! stop looping
    if (.not. MOVEDELTA .or. .not. DELTA) looping = .false.

    ! precalculate the phase velocities for all grid points (delta position may have changed)
    if (looping) call constructPhaseVelocitySquare()

    ! benchmark output
    if (MAIN_PROCESS .and. VERBOSE) then
      benchend = syncWtime()
      print *
      print *,'benchmark seconds ',benchend-benchstart
      print *
    endif
  enddo !delta location

  ! end parallelization
  call syncFinalizeMPI()

  end program
