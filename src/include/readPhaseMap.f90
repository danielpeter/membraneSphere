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
      subroutine readPhaseMap()
!-----------------------------------------------------------------------
! reads in a phase velocity map file from directory /phasedata
! for corresponding grid level size e.g. 'phase.6.dat'
!
! returns: phaseMap array filled up 
      use propagationStartup
      implicit none
      character*64:: fileName
      character*1:: divString
      
      ! parameters
      write(divString,'(i1)') subdivisions  ! grid level
      fileName = '../phasedata/phase.6.L150.dat'   !filename of phase velocity map (values in absolute km/s)
  
      ! read phase map
      call readPhaseVelocityMap(fileName)
      
      end

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
      use phaseVelocityMap; use propagationStartup; use verbosity      
      implicit none
      integer:: i,ioerror,id
      character*64:: fileName
      real(WP):: phasepercent,lat,lon
        
      ! read phase map
      fileName = trim(fileName)
      if(VERBOSE)print*,fileName   
      
      !open file      
      open(1, file= fileName,status='old',iostat=ioerror)
      if( ioerror .ne. 0) then
        print*,'error opening file ',fileName
        call stopProgram( 'abort - readPhasemap   ')
      endif
                  
      !fill in values
      do i=1, numVertices
        read(1, *, iostat=ioerror) lat,lon,phaseMap(i)
        if( ioerror .ne. 0) then
          exit
        endif
        !print*,'vertex',lat,lon,phaseMap(i)
      enddo    
      
      close(1)
      if(VERBOSE)print*,'number of phase velocities read in: ',i-1
      if( numVertices .ne. i-1 ) then
        print*,'read entries do not match'
        call stopProgram( 'abort - readPhaseMap   ')
      endif
      
      end
