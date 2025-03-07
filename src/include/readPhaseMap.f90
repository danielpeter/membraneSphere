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
      integer:: i,ioerror,id
      character*64:: fileName
      real(WP):: phasepercent,lat,lon

      ! read phase map
      fileName = trim(fileName)
      if (VERBOSE)print *,fileName

      !open file
      open(1, file= fileName,status='old',iostat=ioerror)
      if ( ioerror /= 0) then
        print *,'error opening file ',fileName
        call stopProgram( 'abort - readPhasemap   ')
      endif

      !fill in values
      do i = 1, numVertices
        read(1, *, iostat=ioerror) lat,lon,phaseMap(i)
        if ( ioerror /= 0) then
          exit
        endif
        !print *,'vertex',lat,lon,phaseMap(i)
      enddo

      close(1)
      if (VERBOSE)print *,'number of phase velocities read in: ',i-1
      if ( numVertices /= i-1 ) then
        print *,'read entries do not match'
        call stopProgram( 'abort - readPhaseMap   ')
      endif

      end

!-----------------------------------------------------------------------
      subroutine readPhaseMap()
!-----------------------------------------------------------------------
! reads in a phase velocity map file from directory /phasedata
! for corresponding grid level size e.g. 'phase.6.dat'
!
! returns: phaseMap array filled up
      use cells
      implicit none
      character*64:: fileName
      character*1:: divString

      ! parameters
      write(divString,'(i1)') subdivisions  ! grid level
      !filename of phase velocity map (values in absolute km/s)
      fileName = '../phasedata/phase.'//divString//'.L150.dat'

      ! read phase map
      call readPhaseVelocityMap(fileName)

      end
