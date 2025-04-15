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
  subroutine forwardIteration()
!-----------------------------------------------------------------------
! iterates through time and calculates the solutions at each time step
!
! finite-difference iteration: information in carl tape thesis2003, chap 5, (5.7)
  use propagationStartup; use cells; use phaseVelocityMap; use parallel; use displacements
  use phaseBlockData; use griddomain; use adjointVariables; use verbosity
  implicit none
  ! local parameters
  real(WP):: u_t, u_tplus1, u_tminus1,forcing
  real(WP):: D2,time,newdisp
  integer:: n,timestep,vertex,index,jrec
  real(WP), external:: forceTerm2Source,forceTermExact
  real(WP), external:: discreteLaplacian,precalc_discreteLaplacian

  !initialize
  displacement(:)     = 0.0_WP
  displacement_old(:) = 0.0_WP
  newdisplacement(:)  = 0.0_WP

  ! time iteration of displacements
  index = 0
  do timestep = firsttimestep, lasttimestep
    ! model time
    time  = timestep*dt
    index = index+1

    ! swap displacement arrays
    displacement_old(:) = displacement(:)
    displacement(:)     = newdisplacement(:)

    ! propagate only corresponding vertices
    do n = 1, numDomainVertices
      ! choose vertex
      if (PARALLELSEISMO) then
        vertex = domainVertices(n)
      else
        vertex = n
      endif

      ! spherical Laplacian
      if (PRECALCULATED_CELLS) then
        D2 = precalc_discreteLaplacian(vertex) !uses displacement(..) field
      else
        D2 = discreteLaplacian(vertex)
      endif
      ! calculate new displacement in time
      u_t        = displacement(vertex)
      u_tminus1  = displacement_old(vertex)

      ! determines force
      if (PRESCRIBED_SOURCE) then
        ! get the value out of the prescribed force array
        ! prescribed array goes on till numDomainVertices, see indexing at initialization
        if (.not. sourceOnFile) then
          forcing = forceTermPrescribed(n,index)
        else
          jrec = (n-1)*numofTimeSteps + index
          read(sourceFileID,rec=jrec) forcing
        endif
      else
        if (STATION_CORRECTION) then
          forcing = forceTermExact(vertex,time,desiredSourceLat,desiredSourceLon)
        else
          forcing = forceTerm2Source(vertex,time,sourceVertex)
        endif
      endif

      ! gets phase velocity square
      cphase2 = phaseVelocitySquare(vertex)

      !propagation step
      u_tplus1 = u_t + u_t - u_tminus1 + cphase2*dt2*(D2+forcing)

      ! iterated displacements
      newdisplacement(vertex) = u_tplus1
    enddo

    ! synchronize new displacement arrays
    if (PARALLELSEISMO) then
      call syncNewdisplacement()
    endif

    ! fill displacement at receiver station into seismogram record of that receiver station
    ! note: newdisplacement is at time step t+dt
    if (manyReceivers) then
      call recordManySeismograms(index,time+dt,referenceRun)
    else
      ! set time in seismogram
      seismogramReceiver(1,index) = time+dt

      ! set displacement
      if (STATION_CORRECTION) then
        !interpolate new displacement for original receiver location
        call interpolateLinear(interpolation_distances,interpolation_triangleLengths, &
                              interpolation_corners,newdisp)
        seismogramReceiver(2,index) = newdisp
      else
        seismogramReceiver(2,index) = newdisplacement(receiverVertex)
      endif
    endif

    ! save displacements at each time step for future adjoint calculation
    if (Adjoint_Program) then
      call storeForwardDisplacements(timestep,index)
    endif

    ! file output for simulation snapshots
    if (SIMULATIONOUTPUT) then
      call printWavefield(timestep,time,index)
    endif
  enddo !timestep

  end subroutine


!-----------------------------------------------------------------------
  subroutine recordManySeismograms(index,time,refRun)
!-----------------------------------------------------------------------
! fill receivers seismograms with computed displacements
  use propagationStartup; use displacements
  implicit none
  integer,intent(in) :: index
  real(WP),intent(in):: time
  logical,intent(in):: refRun
  ! local parameters
  integer:: m,i

  ! record seismograms
  if (manyKernels) then
    if (refRun) then
      do m = 1,numofKernels
        ! fill in the time seismogram first
        kernelsReceiversSeismogramRef(m,numofReceivers+1,index) = time
        ! fill in the displacement seismograms for each receiver
        do i = 1,numofReceivers
          kernelsReceiversSeismogramRef(m,i,index) = newdisplacement(kernelsReceivers(m,i))
        enddo
      enddo
    else
      do m = 1,numofKernels
        ! fill in the time seismogram first
        kernelsReceiversSeismogram(m,numofReceivers+1,index) = time
        ! fill in the displacement seismograms for each receiver
        do i = 1,numofReceivers
          kernelsReceiversSeismogram(m,i,index) = newdisplacement(kernelsReceivers(m,i))
        enddo
      enddo
    endif
  else
    if (refRun) then
      ! fill in the time seismogram first
      receiversSeismogramRef(numofReceivers+1,index) = time
      ! fill in the displacement seismograms for each receiver
      do i = 1,numofReceivers
        receiversSeismogramRef(i,index) = newdisplacement(receivers(i))
      enddo
    else
      ! fill in the time seismogram first
      receiversSeismogram(numofReceivers+1,index) = time
      ! fill in the displacement seismograms for each receiver
      do i = 1,numofReceivers
        receiversSeismogram(i,index) = newdisplacement(receivers(i))
      enddo
    endif
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine parallelFindnextLocation(move)
!-----------------------------------------------------------------------
! increments delta location and finds new delta vertex
! (remember: different processes are used to compute the same simulation, not different simulations)
!
! input:
!     latDelta/lonDelta     - location of delta
!     vertexDelta             - nearest vertex
!     move                      - iteration flag
!
! return: latDelta,lonDelta,vertexDelta,move for new location
  use loop; use propagationStartup; use deltaSecondLocation
  implicit none
  logical,intent(inout):: move

  ! increment latitude
  if (abs(deltaLat-latitudeEnd) > deltaMoveIncrement/2) then
    deltaLat = deltaLat+deltaMoveIncrement
  else
    ! increment longitude
    ! if (int(latDelta) > 50) then
    deltaLat = latitudeStart
    deltaLon = deltaLon + deltaMoveIncrement
  endif

  ! get nearest vertex for new location
  call findVertex(deltaLat,deltaLon,deltaVertex)

  ! for 2 delta scatterers
  if (SECONDDELTA) call placeSecondDelta(deltaLat,deltaLon,deltaSecondLat, &
                                          deltaSecondLon,deltaSecondVertex,deltaSecondDistance)

  ! stop iteration for high longitude
  if (deltaLon > (longitudeEnd + deltaMoveIncrement/2)) move = .false.

  end subroutine


!-----------------------------------------------------------------------
  subroutine findnextLocation(move)
!-----------------------------------------------------------------------
! increments delta location for each process and finds new delta vertex
! (i.e. each process runs a different simulation)
! input:
!     latDelta/lonDelta     - location of delta
!     vertexDelta             - nearest vertex
!     move                      - iteration flag
!
! return: latDelta,lonDelta,vertexDelta,move for new location
  use loop; use parallel; use propagationStartup; use deltaSecondLocation
  implicit none
  logical,intent(inout):: move
  ! local parameters
  real(WP):: lat,lon,newlon
  logical:: bestfound

  ! increment latitude
  if ( (deltaLat - (latitudeEnd-nprocesses*deltaMoveIncrement)) < deltaMoveIncrement/2) then
    ! not yet reached the end, so increment
    deltaLat = deltaLat+nprocesses*deltaMoveIncrement
  else
    ! deltaLat is at the latitudeEnd or has already passed it
    if (manyReceivers) then
      ! moves only down the latitude at one chosen longitude, when end reached, it stops
      move = .false.
    else
      ! reset latitude to begining (take account of alternating which process will start from latitudeStart)
      deltaLat = latitudeStart + ( deltaLat-latitudeEnd+(nprocesses-1)*deltaMoveIncrement)
      ! increment longitude
      deltaLon = deltaLon+deltaMoveIncrement

      !stop iteration for high longitude
      if ( (deltaLon - longitudeEnd) > deltaMoveIncrement/2) move = .false.
    endif
  endif

  ! get nearest vertex for new location
  call findVertex(deltaLat,deltaLon,deltaVertex)

  ! avoid too big steps
  if (manyReceivers) then
    call getSphericalCoord_Lat(deltaVertex,lat,lon)
    newlon = deltaLon
    bestfound = .false.
    if (abs(lat-deltaLat) <= abs(deltaMoveIncrement/2)) bestfound = .true.

    do while( (.not. bestfound ) .and. (newlon < (deltaLon+360.0))  )
      ! search on the east for a better vertex
      newlon = newlon+deltaMoveIncrement*0.5
      call findVertex(deltaLat,newlon,deltaVertex)
      call getSphericalCoord_Lat(deltaVertex,lat,lon)
      if (abs(lat-deltaLat) <= abs(deltaMoveIncrement/2)) bestfound = .true.
    enddo
    ! ...we found the best vertex (best matching latitude)
    if (.not. bestfound) then
      ! take initial one
      call findVertex(deltaLat,deltaLon,deltaVertex)
      call getSphericalCoord_Lat(deltaVertex,lat,lon)
      ! output
      print *
      print *,'no a very good delta vertex found for position:',deltaLat,deltaLon
      print *,'   found only:',lat,lon
      print *
    endif

    ! for 2 delta scatterers
    if (SECONDDELTA) call placeSecondDelta(deltaLat,newlon,deltaSecondLat, &
                                            deltaSecondLon,deltaSecondVertex,deltaSecondDistance)
  else
    ! for 2 delta scatterers
    if (SECONDDELTA) call placeSecondDelta(deltaLat,deltaLon,deltaSecondLat, &
                                            deltaSecondLon,deltaSecondVertex,deltaSecondDistance)
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine constructPhaseVelocitySquare()
!-----------------------------------------------------------------------
! determines the phase velocities for all vertices and stores them in an array
! with corresponding indices
!
! returns: phaseVelocitySquare array
  use propagationStartup; use phaseVelocityMap; use parallel; use cells
  use verbosity
  implicit none
  ! local parameters
  integer:: i
  real(WP):: getPhaseSquare,lat,lon

  ! store to file
  if (MAIN_PROCESS) then
    if (VERBOSE) then
      print *,'  writing to file: ',trim(datadirectory)//'PhaseMap.dat'
    endif
    open(IOUT,file=trim(datadirectory)//'PhaseMap.dat')
    ! header info
    write(IOUT,*) '# Phase map'
    if (HETEROGENEOUS) then
      write(IOUT,*) '# heterogeneous phase velocity map'
    else
      write(IOUT,*) '# homogeneous phase velocity map'
    endif
    if (DELTA) then
      write(IOUT,*) '# added delta scatterer'
      if (SECONDDELTA) then
        write(IOUT,*) '# added second delta scatterer'
      endif
    endif
    write(IOUT,*) '# format:'
    ! format: #lon #lat #phase-velocity #vertexID
    write(IOUT,*) '#lon #lat #phase-velocity #vertexID'
  endif

  ! for all grid vertices
  phaseVelocitySquare(:) = 0.0_WP
  do i = 1,numVertices
    ! determine phase square
    phaseVelocitySquare(i) = getPhaseSquare(i)

    ! file output
    if (MAIN_PROCESS) then
      call getSphericalCoord_Lat(i,lat,lon)

      ! format: #lon #lat #phase-velocity #vertexID
      write(IOUT,*) lon,lat,sqrt(phaseVelocitySquare(i)),i
    endif
  enddo

  ! close files
  if (MAIN_PROCESS) then
    close(IOUT)
  endif

  end subroutine


!-----------------------------------------------------------------------
  function getPhaseSquare(n)
!-----------------------------------------------------------------------
! calculates the square of phase velocity depending on
! homogeneous, heterogeneous or delta phase maps
! for a certain vertex
!
! input:
!       n             - vertex index
!
! returns: phase velocity square
  use propagationStartup; use phaseVelocityMap; use parallel; use deltaSecondLocation
  use cells; use verbosity
  implicit none
  real(WP):: getPhaseSquare
  ! local parameters
  integer,intent(in):: n
  real(WP):: vectorV(3),vectorS(3),vectorSnd(3),referencePhaseVelocity
  real(WP):: distance,cphase,distanceSnd,lat,lon,csquare

  ! set reference velocity
  referencePhaseVelocity = cphaseRef

  ! homogeneous phase map (or delta phase map)
  csquare = cphaseRef*cphaseRef

  !heterogeneous phase velocities
  if (HETEROGENEOUS) then
    cphase = phaseMap(n)
    csquare = cphase*cphase

    !set as new reference for possible delta scatterers that will follow
    referencePhaseVelocity = cphase
  endif

  !delta phase map
  if (DELTA) then
    ! first and second scatterer vectors
    vectorV(:) = vertices(n,:)
    vectorS(:) = vertices(deltaVertex,:)
    call greatCircleDistance(vectorS,vectorV,distance)
    distance = EARTHRADIUS*distance
    if (SECONDDELTA) then
      vectorSnd(:) = vertices(deltaSecondVertex,:)
      call greatCircleDistance(vectorSnd,vectorV,distanceSnd)
      distanceSnd = EARTHRADIUS*distanceSnd
    endif

    ! determine which phase velocity to apply for first scatterer
    if (distance > DELTARADIUS) then
      cphase = referencePhaseVelocity
    else
      if (trim(DELTAfunction) == "plateau") then
        cphase = referencePhaseVelocity + deltaPerturbation
      else
        if (trim(DELTAfunction) == "Gaussian") then
          ! some experimental function which looks like a Gaussian for distances between 0 and 560 km
          cphase = (1.0-exp(-distance*distance/62500.0))*referencePhaseVelocity&
            + exp(-distance*distance/62500.0)*(referencePhaseVelocity + deltaPerturbation)
        endif
      endif
      ! console output
      if (MAIN_PROCESS .and. VERBOSE) then
            print *,'distance scatterer:',n,distance
            call getSphericalCoord_Lat(n,lat,lon)
            print *,'     delta location (lat/lon):',lat,lon
            print *,'     delta radius:',DELTARADIUS,'cphase:',cphase
      !    !endif
      endif
    endif !distance

    ! determine which phase velocity to apply for second scatterer
    if (SECONDDELTA) then
      if (distanceSnd <= DELTARADIUS) then
        if (trim(DELTAfunction) == "plateau") then
          cphase = referencePhaseVelocity + deltaPerturbation
        else
          if (trim(DELTAfunction) == "Gaussian") then
            ! some experimental function which looks like a Gaussian for distances between 0 and 560 km
            cphase = (1.0-exp(-distanceSnd*distanceSnd/62500.0))*referencePhaseVelocity &
              + exp(-distanceSnd*distanceSnd/62500.0)*(referencePhaseVelocity + deltaPerturbation)
          endif
        endif
        if (MAIN_PROCESS .and. VERBOSE) then
              print *,'distance second scatterer:',n,distanceSnd
              call getSphericalCoord_Lat(n,lat,lon)
              print *,'     delta location (lat/lon):',lat,lon
              print *,'     delta radius:',DELTARADIUS,'cphase:',cphase
              print *
        !    !endif
        endif
      endif
    endif !SECONDDELTA

    ! squared value for phase velocity
    csquare = cphase*cphase
  endif !DELTA

  ! return value
  getPhaseSquare = csquare

  return
  end function


!-----------------------------------------------------------------------
  subroutine printSeismogram()
!-----------------------------------------------------------------------
! plots all seismograms into corresponding files
! file named like:
!     'seismo.*L*.withoutDelta.*#*.dat'    when no delta scatterer was present, else
!     'seismo.*L*.*lat*.*lon*.dat'
! *L* is wave type, *#* is receiver id, *lat* is latitude, *lon* is longitude of delta scatterer
  use parallel;use griddomain;use propagationStartup; use verbosity
  implicit none
  ! local parameters
  integer:: i,n,domain,ier
  character(len=6):: latstr,lonstr
  character(len=3):: recstr
  character(len=128):: filename
  real(WP):: reclat,reclon
  integer, external:: getDomain

  write(latstr,'(f6.1)') deltaLat
  write(lonstr,'(f6.1)') deltaLon
  latstr = trim(latstr)
  lonstr = trim(lonstr)
  cphasetype = trim(cphasetype)

  if (DELTA) then
    filename = trim(datadirectory)//'seismo.'//trim(cphasetype)//'.'//latstr//'.'//lonstr//'.dat'
  else
    filename = trim(datadirectory)//'seismo.'//trim(cphasetype)//'.withoutDelta.dat'
  endif

  ! output
  if (MAIN_PROCESS .and. VERBOSE) print *,'  seismogram output:'

  ! just to remove unneccessary spaces
  filename = trim(filename)

  if (PARALLELSEISMO) then
    ! each process prints its own output to the current directory
    if (MAIN_PROCESS .and. VERBOSE) print *,'    parallel simulation seismos'

    ! displacement at receiver
    if (manyReceivers) then
      print *,'    printing to receiver files: '//filename(1:len_trim(filename)-4)

      do i = 1,size(receivers)
        ! get right seismogram
        domain = getDomain(receivers(i))
        if (myrank == domain) then
          write(recstr,'(i3.2)') int(i)
          ! get actual vertex position
          call getSphericalCoord_Lat(receivers(i),reclat,reclon)

          ! create filename
          if (DELTA) then
            filename = trim(datadirectory)//'seismo.'//trim(cphasetype)//'.'//latstr//'.'//lonstr//'.'//recstr//'.dat'
          else
            filename = trim(datadirectory)//'seismo.'//trim(cphasetype)//'.withoutDelta.'//recstr//'.dat'
          endif

          ! write to file
          open(200+i,file=trim(filename),iostat=ier)
          if (ier /= 0) call stopProgram('could not open '//trim(filename)//'   ')
          write(200+i,*) 'receiver:',i,real(reclat),real(reclon),receivers(i)
          do n = 1,size( receiversSeismogram(size(receivers)+1,:) )
            write(200+i,*) receiversSeismogram(size(receivers)+1,n),receiversSeismogram(i,n)
          enddo
          close(200+i)
        endif
      enddo
    else
      ! get right seismogram
      call syncReceiverSeismogram()

      ! output
      if (MAIN_PROCESS) then
        ! displacement at receiver
        print *,'    printing to file: ',trim(filename)
        open(200,file=trim(filename))
        do n = 1,size(seismogramReceiver(1,:))
          write(200,*) seismogramReceiver(1,n),seismogramReceiver(2,n)
        enddo
        close(200)
      endif
    endif
  else
    ! each process prints its own output to the current directory
    if (MAIN_PROCESS .and. VERBOSE) print *,'    distributed simulation seismos'

    ! displacement at receiver
    if (manyReceivers) then
      print *,'    printing to receiver files: ',filename(1:len_trim(filename)-4)

      do i = 1,size(receivers)
        write(recstr,'(i3.2)') int(i)
        ! get actual vertex position
        call getSphericalCoord_Lat(receivers(i),reclat,reclon)

        ! create filename
        if (DELTA) then
          filename = trim(datadirectory)//'seismo.'//trim(cphasetype)//'.'//latstr//'.'//lonstr//'.'//recstr//'.dat'
        else
          filename = trim(datadirectory)//'seismo.'//trim(cphasetype)//'.withoutDelta.'//recstr//'.dat'
        endif

        ! write to file
        open(200+i,file=trim(filename),iostat=ier)
        if (ier /= 0) call stopProgram('could not open '//trim(filename)//'   ')
        write(200+i,*) 'receiver:',i,real(reclat),real(reclon),receivers(i)
        do n = 1,size( receiversSeismogram(size(receivers)+1,:) )
          write(200+i,*) receiversSeismogram(size(receivers)+1,n),receiversSeismogram(i,n)
        enddo
        close(200+i)
      enddo
    else
      print *,'    printing to file: ',trim(filename)
      open(200,file=trim(filename))
      do n = 1,size(seismogramReceiver(1,:))
        write(200,*) seismogramReceiver(1,n),seismogramReceiver(2,n)
      enddo
      close(200)
    endif
  endif

  end subroutine


!-----------------------------------------------------------------------
  subroutine printWavefield(timestep,time,index)
!-----------------------------------------------------------------------
! prints complete newdisplacement wavefield to file
!
! returns: wavefiled in file #output-dir#/simulation.#time#.dat
  use precisions
  use cells, only: vertices,numVertices,numTriangleFaces,cellTriangleFace
  use parallel, only: MAIN_PROCESS
  use displacements, only: newdisplacement
  use propagationStartup, only: datadirectory
  implicit none
  integer,intent(in):: timestep,index
  real(WP),intent(in):: time
  ! local parameters
  integer:: n,k
  character(len=5):: timestr
  real(WP):: lat,lon
  logical,parameter:: FORMAT_VTK         = .true.  ! outputs as vtk file
  real,parameter:: SIMULATION_STARTTIME = -50.0

  ! only main process writes to files
  if (mod(timestep,SIMULATION_TIMESTEPPING) == 0 &
        .and. time >= SIMULATION_STARTTIME .and. MAIN_PROCESS) then
    write(timestr,'(i5.5)') index

    if (.not. FORMAT_VTK) then
      open(IOUT,file=trim(datadirectory)//'simulation.'//timestr//'.dat')
      do n = 1, numVertices
        ! #format: x, y, z, displacement  ( in Cartesian coordinates )
        !write(10,'(4f18.6)') (vertices(n,k),k=1,3),newdisplacement(n)

        ! #format: lon, lat, displacement  ( in spherical coordinates and degrees )
        call getSphericalCoordinates(vertices(n,:),lat,lon)
        write(IOUT,*) lon,lat,newdisplacement(n)

      enddo
      close(IOUT)
      print *,'    file written: ',trim(datadirectory)//'simulation.'//timestr//'.dat'
    else
      ! vtk file
      open(IOUT,file=trim(datadirectory)//'simulation.'//timestr//'.vtk')

      write(IOUT,'(a26)') "# vtk DataFile Version 3.1"
      write(IOUT,'(a14)') "membraneSphere"
      write(IOUT,'(a5)')  "ASCII"
      write(IOUT,'(a16)') "DATASET POLYDATA"
      !write(IOUT,'(a25)') "DATASET UNSTRUCTURED_GRID"
      write(IOUT,'(a6,i18,1x,a5)') "POINTS",numVertices,"float"
      do n = 1,numVertices
        write(IOUT,'(3(f16.8))') (vertices(n,k),k=1,3)
      enddo
      write(IOUT,*)

      write(IOUT,'(a8,i18,i18)') "POLYGONS",numTriangleFaces,numTriangleFaces*4
      !write(IOUT,'(a5,i,i)') "CELLS",numTriangleFaces,numTriangleFaces*4
      ! VTK starts indexing at 0
      ! (by that we have to shift the arrays by -1)
      do n = 1,numTriangleFaces
        write(IOUT,*) 3,cellTriangleFace(n,1)-1,cellTriangleFace(n,2)-1,cellTriangleFace(n,3)-1
      enddo
      write(IOUT,*)

      !write(IOUT,'(a10,1x,i)') "CELL_TYPES",numTriangleFaces
      !write(IOUT,*) (5,k=1,numTriangleFaces)
      !write(IOUT,*)

      write(IOUT,'(a10,i18)') "POINT_DATA",numVertices
      write(IOUT,'(a26)') "SCALARS displacement float"
      write(IOUT,'(a20)') "LOOKUP_TABLE default"
      do n = 1,numVertices
        write(IOUT,'(f16.8)') newdisplacement(n)
      enddo
      write(IOUT,*)
      close(IOUT)

      print *,'    file written: ',trim(datadirectory)//'simulation.'//timestr//'.vtk'
    endif
  endif

  end subroutine
