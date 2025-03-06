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
      subroutine readPrecalculated()
!-----------------------------------------------------------------------
! reads in cellAreas, cellEdgesLength, cellCenterDistances
! (units are: cellArea [km2], length and distances [km] )
      use cells; use propagationStartup;use verbosity
      implicit none
      integer:: i,k,ioerror,count, iend
      character*1:: divString
      character*23:: fileName
           
      write(divString,'(i1)') subdivisions
! 100  format(f16.14,2x,f16.14,2x,f16.8)      
 
      ! read cellArea
      if(DEBUG)print*,'read in cellAreas ...'
      open(1, file=  '../griddata/cellAreas.'//divString//'.dat',status='old',iostat=ioerror)
      if( ioerror .ne. 0) call stopProgram( 'abort - readPrecalculated cellAreas    ')

      do i=1, numVertices
        read(1, *, iostat=ioerror) cellAreas(i)
        if( ioerror .ne. 0) then
          print*,'error at: ',i
          exit
        endif
      enddo     
      close(1)
      !debug
      if(DEBUG) print*,'cellArea',1,cellAreas(1)

      ! read cellEdgesLength
      if(DEBUG)print*,'read in cellEdgesLength ...'
      open(1, file=  '../griddata/cellEdgesLength.'//divString//'.dat',status='old',iostat=ioerror)
      if( ioerror .ne. 0) call stopProgram( 'abort - readPrecalculated cellEdgesLength   ')

      do i=1, numVertices
        read(1, *, iostat=ioerror) cellEdgesLength(i,0:6)
        if( ioerror .ne. 0) then
          print*,'error at: ',i
          exit
        endif
      enddo     
      close(1)
      !debug
      if(DEBUG)print*,'cellEdgesLength',1,(cellEdgesLength(1,k),k=0,6)

      ! read cellCenterDistances
      if(DEBUG)print*,'read in cellCenterDistances ...'
      open(1,file= '../griddata/cellCenterDistances.'//divString//'.dat',status='old',iostat=ioerror)
      if( ioerror .ne. 0) call stopProgram( 'abort - readPrecalculated cellEdgesLength   ')

      do i=1, numVertices
        read(1, *, iostat=ioerror) cellCenterDistances(i,0:6)
        if( ioerror .ne. 0) then
          print*,'error at: ',i
          exit
        endif
      enddo     
      close(1)
      !debug
      if(DEBUG)print*,'cellCenterDistances',1,(cellCenterDistances(1,k),k=0,6)
      
      if(DEBUG) print*
      
      end

