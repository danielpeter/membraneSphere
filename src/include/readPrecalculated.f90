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
      subroutine readPrecalculated()
!-----------------------------------------------------------------------
! reads in cellAreas, cellEdgesLength, cellCenterDistances
! (units are: cellArea [km2], length and distances [km] )
      use cells; use verbosity
      implicit none
      integer:: i,k,ioerror,count, iend
      character(len=1):: divString
      character(len=56):: fileName
      character(len=14):: ending

      !initialize
      if (RELAXED_GRID) then
        ending='.relaxed.dat'
      else
        ending='.dat'
      endif

      write(divString,'(i1)') subdivisions
! 100  format(f16.14,2x,f16.14,2x,f16.8)

      ! read cellArea
      open(1,file='../griddata/cellAreas.'//divString//trim(ending),status='old',iostat=ioerror)
      if ( ioerror /= 0) call stopProgram( 'abort - readPrecalculated cellAreas    ')

      do i = 1, numVertices
        read(1, *, iostat=ioerror) cellAreas(i)
        if ( ioerror /= 0) then
          print *,'error at: ',i
          exit
        endif
      enddo
      close(1)

      ! read cellEdgesLength
      open(1,file='../griddata/cellEdgesLength.'//divString//trim(ending),status='old',iostat=ioerror)
      if ( ioerror /= 0) call stopProgram( 'abort - readPrecalculated cellEdgesLength   ')

      do i = 1, numVertices
        read(1, *, iostat=ioerror) cellEdgesLength(i,0:6)
        if ( ioerror /= 0) then
          print *,'error at: ',i
          exit
        endif
      enddo
      close(1)

      ! read cellCenterDistances
      open(1,file='../griddata/cellCenterDistances.'//divString//trim(ending),status='old',iostat=ioerror)
      if ( ioerror /= 0) call stopProgram( 'abort - readPrecalculated cellEdgesLength   ')

      do i = 1, numVertices
        read(1, *, iostat=ioerror) cellCenterDistances(i,0:6)
        if ( ioerror /= 0) then
          print *,'error at: ',i
          exit
        endif
      enddo
      close(1)

      ! get heikes&randall ratios for each hexagonal cell edge
      if ( CORRECT_RATIO ) then
        open(1,file='../griddata/cellFractions.'//divString//trim(ending),status='old',iostat=ioerror)
        if ( ioerror /= 0) call stopProgram( 'abort - readPrecalculated cellFractions   ')

        do i = 1, numVertices
          read(1, *, iostat=ioerror) cellFractions(i,0:6)
          if ( ioerror /= 0) then
            print *,'error at: ',i
            exit
          endif
        enddo
        close(1)
        !debug
        !print *,'cellFractions',1,(cellFractions(1,k),k=0,6)
        !print *
      endif

      end

