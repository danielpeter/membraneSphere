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
      subroutine readData(VERBOSE)
!-----------------------------------------------------------------------
! gets datas from corresponding data files in ../griddata directory
!
! needs to be set:
!   subdivisions - number of subfoldings
!
! fills up:
!   vertices, numVertices       - cell centers/triangle corners positions
!   cellCorners, numCorners - cell corners
!   cellNeighbors,numNeighbors - cell centers around referenced vertex
!   cellFace, numFaces        - cells
!   cellTriangleFace, numTriangleFaces - triangles
      !use propagationStartup;use cells;use verbosity;use adjointVariables
      use propagationStartup, only: SIMULATIONOUTPUT
      use cells
      implicit none
      logical,intent(in):: VERBOSE
      ! local parameters
      integer::i,k,ier,ilocal
      character(len=1):: divString
      character(len=56):: fileName
      character(len=14):: ending

      if (VERBOSE) then
        print *,'  reading grid data:'
        print *,'    subdivisions:',subdivisions
      endif

      !initialize
      if (RELAXED_GRID) then
        ending='.relaxed.dat'
      else
        ending='.dat'
      endif

      numVertices = 0
      write(divString,'(i1)') subdivisions

!100  format(f16.14,2x,f16.14,2x,f16.8)
!101  format(f16.10,2x,f16.10,2x,f16.10)

      ! read all voronoi cell centre vertices values from file
      !(voronoi cell centers are equal to triangle corners)
      fileName = 'data/griddata/Dtvert'//divString//ending
      if (VERBOSE) print *,'  ',trim(fileName)
      open(10,file=trim(fileName),status='old',iostat=ier)
      if (ier /= 0) then
        print *,'Error: could not open file ',trim(fileName)
        print *,'       Please check if the file exists...'
        call stopProgram( 'abort - readData File1   ')
      endif

      !print *,'getting voronoi cell centre vertices from',fileName
      numVertices = 0
      ilocal = 0
      do i = 1, MaxVertices
        read(10,*, iostat=ier) vertices(i,1), vertices(i,2),vertices(i,3)
        if (ier /= 0) exit
        numVertices = numVertices +1
        ! find locations
        !call getSphericalCoordinates(vertices(i,:),lat,lon)
        !if (abs(lon - 90) < 0.4 .or. abs(lon - 270) < 0.4) then
        !  ilocal = ilocal+1
        !  print *,'vertex ',ilocal,' :',lat,lon
        !endif
      enddo
      if (VERBOSE) print *,'   number of vertices: ',numVertices
      close(10)
      !debug
      !print *,'vertex',1,(vertices(1,k),k=1,3)

      !read in the corresponding cell corners
      fileName = 'data/griddata/Dvvert'//divString//ending
      if (VERBOSE) print *,'  ',trim(fileName)
      open(20,file=trim(fileName),status='old',iostat=ier)
      if (ier /= 0) call stopProgram( 'abort - readData File2   ')

      !print *,'getting voronoi cell corner vertices...'
      numCorners = 0
      do i = 1, MaxTriangles
        read(20,*,iostat=ier) (cellCorners(i,k),k=1,3)
        if (ier /= 0) exit
        numCorners = numCorners +1
      enddo
      if (VERBOSE) print *,'   number of corners: ',numCorners
      close(20)
      !debug
      !print *,'vertex',1,(cellCorners(1,k),k=1,3)

      !read in the corresponding cell neighbors indices
      fileName = 'data/griddata/Dnear'//divString//ending
      if (VERBOSE) print *,'  ',trim(fileName)
      open(30,file=trim(fileName),status='old',iostat=ier)
      if (ier /= 0) call stopProgram( 'abort - readData File3   ')

      !print *,'getting voronoi cell neighbors...'
      numNeighbors = 0
      do i = 1, MaxVertices
        read(30,*) (cellNeighbors(i,k),k=1,6)
        numNeighbors = numNeighbors +1
      enddo
      if (VERBOSE) print *,'   number of neighbors: ',numNeighbors
      close(30)

      !prepare cellNeighbors info
      do i = 1, numNeighbors
        if (cellNeighbors(i,1) == cellNeighbors(i,6)) then
          cellNeighbors(i,0) = 5
        else
          cellNeighbors(i,0) = 6
        endif
      enddo
      !debug
      !print *,'vertex',1,(cellNeighbors(1,k),k=0,6)

      !read in the corresponding cell face indices
      fileName = 'data/griddata/Dvface'//divString//ending
      if (VERBOSE) print *,'  ',trim(fileName)
      open(40,file=trim(fileName),status='old',iostat=ier)
      if (ier /= 0) call stopProgram( 'abort - readData File4   ')

      !print *,'getting voronoi cell face indices...'
      numFaces = 0
      do i = 1, MaxVertices
        read(40,*) (cellFace(i,k),k=1,6)
        numFaces = numFaces +1
      enddo
      if (VERBOSE) print *,'   number of faces: ',numFaces
      close(40)

      !prepare cellNeighbors info
      do i = 1, numFaces
        if (cellFace(i,1) == cellFace(i,6)) then
          cellFace(i,0) = 5
        else
          cellFace(i,0) = 6
        endif
      enddo
      !debug
      !print *,'cellFace',1,(cellFace(1,k),k=0,6)

      ! allocate triangle face array
      if (Station_Correction .or. SIMULATIONOUTPUT) then
        !read in the corresponding triangle face indices
        fileName = 'data/griddata/Dtface'//divString//ending
        if (VERBOSE) print *,'  ',trim(fileName)
        open(50,file=trim(fileName),status='old', iostat=ier)
        if (ier /= 0) call stopProgram( 'abort - readData Dtface file    ')

        numTriangleFaces = 0
        do i = 1, MaxTriangles
          if (numTriangleFaces == MaxTriangles) call stopProgram( 'abort-readData triangles    ')

          read(50,*,iostat=ier) (cellTriangleFace(i,k),k=1,3)
          if (ier /= 0) exit
          numTriangleFaces = numTriangleFaces +1
        enddo
        if (VERBOSE) print *,'   number of triangle faces: ',numTriangleFaces
        close(50)
      endif

      end subroutine



!-----------------------------------------------------------------------
      subroutine allocateMesh()
!-----------------------------------------------------------------------
! allocates the arrays which hold the mesh informations
      use cells; use parallel; use verbosity
      use propagationStartup, only: SIMULATIONOUTPUT
      implicit none
      ! local parameters
      integer:: ier

      ! sets maximum numbers of triangles and vertices
      ! (see Tape, 2003: table 4.3, p. 40)
      MaxVertices  = 30*(4**subdivisions) + 2
      MaxTriangles = 60*(4**subdivisions) ! ! from dedecv: 20*(4**(subdivisions+1) - 1)

      ! allocate arrays for mesh
      if (MAIN_PROCESS .and. VERBOSE) then
        print *,'  allocating vertices, cellNeighbors, cellFace and cellCorner arrays:'
        print *,'    max vertices  : ',MaxVertices
        print *,'    max triangles : ',MaxTriangles
        print *,'    size                : ',(17*MaxVertices*WP+MaxTriangles*WP)/1024./1024.,'Mb'
      endif
      allocate(vertices(MaxVertices,3),cellNeighbors(MaxVertices,0:6), &
            cellFace(MaxVertices,0:6),cellCorners(MaxTriangles,3),stat=ier )
      if (ier /= 0) call stopProgram('error in allocating arrays for cell face,..')

      if (Station_Correction .or. SIMULATIONOUTPUT) then
        if (MAIN_PROCESS .and. VERBOSE) then
          print *,'  allocating cellTriangleFace array:'
          print *,'    size                : ',3*MaxTriangles*WP/1024./1024.,'Mb'
        endif
        allocate( cellTriangleFace(MaxTriangles,3),stat=ier)
        if (ier /= 0) call stopProgram('error allocating cellTriangleFace')
      endif
      end subroutine


!-----------------------------------------------------------------------
      subroutine allocateData()
!-----------------------------------------------------------------------
! allocates the arrays which are needed for displacement informations
      use cells; use displacements; use phaseVelocityMap
      use parallel; use verbosity
      implicit none
      ! local parameters
      integer:: ier

      ! allocates displacement arrays
      if (MAIN_PROCESS .and. VERBOSE) then
        print *,'  allocating displacement, cell and phase velocity arrays:'
        print *,'  displacements array size  :  ',3*numVertices*WP/1024./1024.,'Mb'
        print *,'  cells array size          :  ',(2*numVertices*7*WP+numVertices*WP)/1024./1024.,'Mb'
        print *,'  phase velocity array size :  ',numVertices*WP/1024./1024.,'Mb'
        print *
      endif

      allocate(displacement(numVertices),displacement_old(numVertices), &
              newdisplacement(numVertices), stat=ier )
      if (ier /= 0) call stopProgram('error in allocating displacement arrays   ')
      displacement_old(:) = 0.0_WP
      displacement(:)     = 0.0_WP
      newdisplacement(:)  = 0.0_WP

      ! allocate new arrays for precalculated cell attributes
      allocate( cellAreas(numVertices),cellEdgesLength(numVertices,0:6),cellCenterDistances(numVertices,0:6), &
              stat=ier )
      if (ier /= 0) call stopProgram('error in allocating arrays for cell area,..   ')
      cellAreas(:)        = 0.0_WP
      cellEdgesLength(:,:)     = 0.0_WP
      cellCenterDistances(:,:) = 0.0_WP

      ! allocates phase velocity
      allocate( phaseVelocitySquare(numVertices),stat=ier )
      if (ier /= 0) call stopProgram('error allocating phaseVelocitySquare   ')
      phaseVelocitySquare(:)   = 0.0_WP

      ! allocates memory for array of heikes&randall ratios (calculated by ../geometry/areas)
      if (CORRECT_RATIO) then
        ! check
        if (.not. PRECALCULATED_CELLS) then
          print *,'Error: Laplacian calculation can not take account of midpoint ratios'
          call stopProgram('error correct_ratio   ')
        endif
        if (MAIN_PROCESS .and. VERBOSE) then
          print *,'  allocating cellFractions array:'
          print *,'  size: ',7*numVertices*WP/1024./1024.,'Mb'
          print *
        endif
        allocate( cellFractions(numVertices,0:6),stat=ier)
        if (ier /= 0) call stopProgram('error allocating cellFractions')
        ! initialize
        cellFractions(:,:) = 0.0_WP
      endif

      end subroutine
