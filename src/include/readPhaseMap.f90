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
  subroutine readPhaseVelocityMap(fileName)
!-----------------------------------------------------------------------
! reads in a phase velocity map file from directory /phasedata
! for corresponding grid level size e.g. 'phase.6.dat'
!
! input:
!     fileName -  phase velocity map file
!
! returns: phaseMap array filled up
  use phaseVelocityMap; use propagationStartup; use verbosity; use cells
  implicit none
  character(len=64),intent(in):: fileName
  ! local parameters
  integer:: i,ier,id
  real(WP):: phasepercent,lat,lon

  ! read phase map
  if (VERBOSE) print *,trim(fileName)

  !open file
  open(IIN, file=trim(fileName),status='old',iostat=ier)
  if (ier /= 0) then
    print *,'Error: opening file ',trim(fileName)
    call stopProgram( 'abort - readPhasemap   ')
  endif

  !fill in values
  do i = 1, numVertices
    read(IIN, *, iostat=ier) lat,lon,phaseMap(i)
    if (ier /= 0) then
      exit
    endif
    !print *,'vertex',lat,lon,phaseMap(i)
  enddo

  close(IIN)

  if (VERBOSE) print *,'number of phase velocities read in: ',i-1
  if (numVertices /= i-1) then
    print *,'Error: read entries do not match'
    call stopProgram( 'abort - readPhaseMap   ')
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine readPhaseMap()
!-----------------------------------------------------------------------
! reads in a phase velocity map file from directory /phasedata
! for corresponding grid level size e.g. 'phase.6.dat'
!
! returns: phaseMap array filled up
  use cells
  implicit none
  ! local parameters
  character(len=64):: fileName
  character(len=1):: divString

  ! parameters
  write(divString,'(i1)') subdivisions  ! grid level

  !filename of phase velocity map (values in absolute km/s)
  fileName = 'data/phasedata/phase.'//divString//'.L150.dat'

  ! read phase map
  call readPhaseVelocityMap(fileName)

  end subroutine
