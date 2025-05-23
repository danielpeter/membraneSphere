!=====================================================================
!
!       m e m b r a n e S p h e r e
!       --------------------------------------------------
!      (c) 2025
!
!=====================================================================

!-----------------------------------------------------------------------
  subroutine readPrecalculated()
!-----------------------------------------------------------------------
! reads in cellAreas, cellEdgesLength, cellCenterDistances
! (units are: cellArea [km2], length and distances [km] )
  use cells; use verbosity
  implicit none
  ! local parameters
  integer:: i,ier
  real(WP),dimension(0:6) :: tmp_vertex
  character(len=2):: strLev
  character(len=14):: ending
  character(len=128):: fileName

  if (VERBOSE) then
    print *,'  reading precalculated grid data:'
    print *,'    subdivisions            : ',subdivisions
  endif

  !initialize
  if (RELAXED_GRID) then
    ending = '.relaxed.dat'
  else
    ending = '.dat'
  endif

  write(strLev,'(i0)') subdivisions
  ! left adjust
  strLev = adjustl(strLev)

! 100  format(f16.14,2x,f16.14,2x,f16.8)

  ! read cellArea
  fileName = 'data/griddata/cellAreas.'//trim(strLev)//trim(ending)
  if (VERBOSE) print *,'    ',trim(fileName)
  open(IIN,file=trim(filename),status='old',iostat=ier)
  if (ier /= 0) call stopProgram( 'abort - readPrecalculated cellAreas    ')

  do i = 1, numVertices
    read(IIN, *, iostat=ier) cellAreas(i)
    if (ier /= 0) then
      print *,'error at: ',i
      exit
    endif
  enddo
  close(IIN)

  ! read cellEdgesLength
  fileName = 'data/griddata/cellEdgesLength.'//trim(strLev)//trim(ending)
  if (VERBOSE) print *,'    ',trim(fileName)
  open(IIN,file=trim(filename),status='old',iostat=ier)
  if (ier /= 0) call stopProgram( 'abort - readPrecalculated cellEdgesLength   ')

  do i = 1, numVertices
    read(IIN, *, iostat=ier) tmp_vertex
    if (ier /= 0) then
      print *,'error at: ',i
      exit
    endif
    cellEdgesLength(i,0:6) = tmp_vertex
  enddo
  close(IIN)

  ! read cellCenterDistances
  fileName = 'data/griddata/cellCenterDistances.'//trim(strLev)//trim(ending)
  if (VERBOSE) print *,'    ',trim(fileName)
  open(IIN,file=trim(filename),status='old',iostat=ier)
  if (ier /= 0) call stopProgram( 'abort - readPrecalculated cellEdgesLength   ')

  do i = 1, numVertices
    read(IIN, *, iostat=ier) tmp_vertex
    if (ier /= 0) then
      print *,'error at: ',i
      exit
    endif
    cellCenterDistances(i,0:6) = tmp_vertex
  enddo
  close(IIN)

  ! get heikes&randall ratios for each hexagonal cell edge
  if (CORRECT_RATIO) then
    fileName = 'data/griddata/cellFractions.'//trim(strLev)//trim(ending)
    if (VERBOSE) print *,'    ',trim(fileName)
    open(IIN,file=trim(filename),status='old',iostat=ier)
    if (ier /= 0) call stopProgram( 'abort - readPrecalculated cellFractions   ')

    do i = 1, numVertices
      read(IIN, *, iostat=ier) tmp_vertex
      if (ier /= 0) then
        print *,'error at: ',i
        exit
      endif
      cellFractions(i,0:6) = tmp_vertex
    enddo
    close(IIN)
    !debug
    !print *,'cellFractions',1,(cellFractions(1,k),k=0,6)
    !print *
  endif

  if (VERBOSE) print *

  end subroutine

