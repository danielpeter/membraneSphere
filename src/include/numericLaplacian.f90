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
  function discreteLaplacian(vertex)
!-----------------------------------------------------------------------
! calculates the laplacian of a grid point given by index
!
! (formula used: Tape, 2003, Chap. 4, 4.13 )
  use displacements; use cells; use verbosity; use propagationStartup
  implicit none
  integer,intent(in):: vertex
  ! local parameters
  real(WP):: discreteLaplacian
  integer:: i,count
  real(WP):: centerDistances(0:6),edgesLength(0:6), &
             laplace,area,sumNeighbors, sumFactor

  !get spherical area of vertex's cell
  call getCellArea(vertex,area)
  !get arc distances from center to neighbors
  call getCellCenterDistances(vertex,centerDistances)
  !get cell edges length
  call getCellEdgesLength(vertex,edgesLength)

  !calculate laplacian approximation
  sumFactor = 0.0
  sumNeighbors = 0.0
  count = 0
  do i = 1, cellNeighbors(vertex,0)
    !factor for own displacement
    if (abs(centerDistances(i)*area) < 1.0e-9) then
      print *,'Error: invalid size in calculation for laplacian'
      print *,'       vertex:',vertex
      print *,'       neighbor:',i
      print *,'       centerDistance:',centerDistances(i)
      print *,'       area:',area
      call stopProgram( 'discreteLaplacian - centerDistance too small   ')
    endif

    sumFactor = sumFactor + edgesLength(i)/(centerDistances(i)*area)

    !influence by neighbors displacement
    sumNeighbors = sumNeighbors + (edgesLength(i)/(centerDistances(i)*area))* &
                                  displacement(cellNeighbors(vertex,i))


    !check
    count = count +1
    if (count > 6) call stopProgram( 'abort-discreteLaplacian neighbors   ')

    !debug
    !if (vertex == receiverVertex) then
    !  print *,'vertex ',vertex,' neighbor:',i, sumFactor, sumNeighbors
    !  print *,'edge ',edgesLength(i),' centers ',centerDistances(i)
    !  print *,'displacement ',cellNeighbors(vertex,i),displacement(cellNeighbors(vertex,i))
    !endif
  enddo

  laplace = sumNeighbors - sumFactor*displacement(vertex)
  discreteLaplacian = laplace

  return
  end function

!-----------------------------------------------------------------------
  function precalc_discreteLaplacian(vertex)
!-----------------------------------------------------------------------
! calculates the laplacian of a grid point given by index,
! uses precalculated values for cell areas, cell center distances and cell edge lengths
!
! (formula used: Tape, 2003, Chap. 4, 4.13 )
! needs arrays cellAreas(), cellEdgesLength(),cellCenterDistances(),
! cellNeighbors() and displacement()
  use cells; use displacements; use verbosity; use propagationStartup
  implicit none
  integer,intent(in):: vertex
  ! local parameters
  real(WP):: precalc_discreteLaplacian
  integer:: i,count !,index
  real(WP):: area,areareci,sumNeighbors,sumFactor
  real(WP):: ratioFactor,hr_ratio,alpha_i,alpha_zero !,alpha_minus1,alpha_plus1

  ! get spherical area of vertex's cell
  area = cellAreas(vertex)
  areareci = 1.0_WP/area

  !calculate laplacian approximation
  sumFactor = 0.0_WP
  sumNeighbors = 0.0_WP
  ratioFactor = 0.0_WP
  count = 0
  do i = 1, cellNeighbors(vertex,0)
    !factor for own displacement
    sumFactor = sumFactor + cellEdgesLength(vertex,i)/cellCenterDistances(vertex,i)

    !influence by neighbors displacement
    sumNeighbors = sumNeighbors + cellEdgesLength(vertex,i)/ &
          cellCenterDistances(vertex,i)*displacement(cellNeighbors(vertex,i))

    ! get term with heikes&randall ratio of distance between midpoints to length of edge
    if (CORRECT_RATIO) then

!          ! use the interpolated values at hexagonal grid corners  d/dn ( d/dx alpha)
!          index = i
!          if (index > 1 .and. index < cellNeighbors(vertex,0)) then
!            alpha_minus1 = displacement( cellNeighbors(vertex,index-1) )
!            alpha_plus1 = displacement( cellNeighbors( vertex,index+1) )
!          else if (index == 1) then
!            alpha_minus1 = displacement( cellNeighbors(vertex,cellNeighbors(vertex,0)) )
!            alpha_plus1 = displacement( cellNeighbors( vertex,index+1) )
!          else if (index == cellNeighbors(vertex,0)) then
!            alpha_minus1 = displacement( cellNeighbors(vertex,index-1) )
!            alpha_plus1 = displacement( cellNeighbors( vertex,1) )
!          endif
!
!          hr_ratio = cellFractions( vertex,i )
!          ratioFactor = ratioFactor + hr_ratio*cellEdgesLength(vertex,i)*cellEdgesLength(vertex,i)* &
!                (alpha_plus1 - alpha_minus1)/(cellCenterDistances(vertex,i)*3.0)

      ! use d/dx ( d/dn alpha) = 1/l_i * ( alpha_i - alpha_0 )/L_i
      alpha_i = displacement(cellNeighbors(vertex,i))
      alpha_zero = displacement(vertex)
      hr_ratio = cellFractions( vertex,i )
      ratioFactor = ratioFactor + hr_ratio*cellEdgesLength(vertex,i)* &
            (alpha_i - alpha_zero)/cellCenterDistances(vertex,i)
    endif
  enddo

  sumFactor    = sumFactor * displacement(vertex) * areareci
  sumNeighbors = sumNeighbors * areareci
  ratioFactor  = ratioFactor * areareci

  precalc_discreteLaplacian = sumNeighbors - sumFactor - ratioFactor

  !debug
  ! check if NaN
  !if (precalc_discreteLaplacian /= precalc_discreteLaplacian) then
  !  print *,'precalculated laplacian error:',precalc_discreteLaplacian,vertex, &
  !      sumFactor,sumNeighbors,area,areareci,count
  !  call stopProgram('precalc_discreteLaplacian error     ')
  !endif

  return
  end function


!-----------------------------------------------------------------------
  function precalc_backdiscreteLaplacian(vertex)
!-----------------------------------------------------------------------
! calculates the laplacian of a grid point given by index for backward calculation of previous forward displacement fields
!
! (formula used: Tape, 2003, Chap. 4, 4.13 )
! needs arrays cellAreas(), cellEdgesLength(),cellCenterDistances(),
! cellNeighbors() and backwarddisplacement()
  use cells; use propagationStartup; use adjointVariables; use verbosity
  implicit none
  integer,intent(in):: vertex
  ! local parameters
  integer:: i,count
  real(WP):: precalc_backdiscreteLaplacian,area,areareci,sumNeighbors,sumFactor

  !get spherical area of vertex's cell
  area = cellAreas(vertex)
  areareci = 1.0_WP/area

  !calculate laplacian approximation
  sumFactor = 0.0_WP
  sumNeighbors = 0.0_WP
  count = 0
  do i = 1, cellNeighbors(vertex,0)
    !factor for own displacement
    sumFactor = sumFactor+cellEdgesLength(vertex,i)/cellCenterDistances(vertex,i)

    !influence by neighbors displacement
    sumNeighbors = sumNeighbors+cellEdgesLength(vertex,i)/ &
          cellCenterDistances(vertex,i)*backwardDisplacement(cellNeighbors(vertex,i))

    !check
    count = count +1
  enddo

  sumFactor = sumFactor*backwardDisplacement(vertex)*areareci
  sumNeighbors = sumNeighbors*areareci
  precalc_backdiscreteLaplacian = sumNeighbors - sumFactor

  !debug
  ! check if Nan
  !if (precalc_backdiscreteLaplacian /= precalc_backdiscreteLaplacian) then
  !  print *,'precalculated backwardlaplacian error:',precalc_backdiscreteLaplacian, &
  !            vertex,sumFactor,sumNeighbors,area,areareci,count
  !  call stopProgram('precalc_backdiscreteLaplacian error     ')
  !endif

  return
  end function

