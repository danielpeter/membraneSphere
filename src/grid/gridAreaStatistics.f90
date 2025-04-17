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

! Description:
!
!  calculates the variations of the hexagonal face areas given by the hexagonal grid points.
!  data read in from files 'DvvertN.dat' and 'DvfaceN.dat'
!
!-----------------------------------------------------------------------
  program areas
!-----------------------------------------------------------------------
  use cells; use verbosity
  implicit none
  character(len=2) :: strLev
  logical, parameter :: STATFILE_OUTPUT = .true.
  logical, parameter :: PLOTFILE_OUTPUT = .false.

  print *,'hexagonal face areas variations'
  print *,'-------------------------------'

  ! console input
  print *, 'number of subdivisions to use: ( >= 0 )'
  read(*,*) subdivisions

  print *
  print *,'subdivisions: ',subdivisions
  print *

  !data input
  print *,'getting data ...'
  VERBOSE = .true.

  ! allocates arrays for mesh informations
  call allocateMesh()

  ! reads in sudivisions data files
  call readData(VERBOSE)

  ! file outputs
  write(strLev,'(i0)') subdivisions
  ! left adjust string number
  strLev = adjustl(strLev)

  if (STATFILE_OUTPUT) then
    !for print out of cells for gnuplot
    open(11,file='OUTPUT/cellAreas.'//trim(strLev)//'.dat')
    open(12,file='OUTPUT/cellFractionAverage.'//trim(strLev)//'.dat')

    open(13,file='OUTPUT/cellEdgesLength.'//trim(strLev)//'.dat')
    open(14,file='OUTPUT/cellCenterDistances.'//trim(strLev)//'.dat')
    open(15,file='OUTPUT/cellFractions.'//trim(strLev)//'.dat')
  endif

  if (PLOTFILE_OUTPUT) then
    open(21,file='OUTPUT/cellCenterLineMidpoints.'//trim(strLev)//'.dat')
    open(22,file='OUTPUT/cellEdgeMidpoints.'//trim(strLev)//'.dat')
    open(23,file='OUTPUT/cellCenterLines.'//trim(strLev)//'.dat')

    open(31,file='OUTPUT/cellCenterLineMidpoints.'//trim(strLev)//'.mat')
    open(32,file='OUTPUT/cellEdgeMidpoints.'//trim(strLev)//'.mat')
    open(33,file='OUTPUT/cellCenterLines.'//trim(strLev)//'.mat')
  endif

  ! area statistics calculations
  print *
  print *,'statistics'
  print *
  call dataVariation(STATFILE_OUTPUT,PLOTFILE_OUTPUT)

  if (STATFILE_OUTPUT) then
    !end print out
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    print *,'files written:'
    print *,'  OUTPUT/cellAreas.'//trim(strLev)//'.dat'
    print *,'  OUTPUT/cellFractionAverage.'//trim(strLev)//'.dat'
    print *,'  OUTPUT/cellEdgesLength.'//trim(strLev)//'.dat'
    print *,'  OUTPUT/cellCenterDistances.'//trim(strLev)//'.dat'
    print *,'  OUTPUT/cellFractions.'//trim(strLev)//'.dat'
    print *
  endif

  if (PLOTFILE_OUTPUT) then
    close(21)
    close(22)
    close(23)

    close(31)
    close(32)
    close(33)

    print *,'  OUTPUT/cellCenterLineMidpoints.'//trim(strLev)//'.dat'
    print *,'  OUTPUT/cellEdgeMidpoints.'//trim(strLev)//'.dat'
    print *,'  OUTPUT/cellCenterLines.'//trim(strLev)//'.dat'
    print *
  endif

  end program


!-----------------------------------------------------------------------
  subroutine stopProgram(textline)
!-----------------------------------------------------------------------
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

