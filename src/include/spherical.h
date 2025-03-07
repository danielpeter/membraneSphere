!
! spherical.h
!
! include file

  integer:: MAXLEVELS, MAXTRIANGLES, MAXVERTICES     ! maximum values for grid subfolding level,...
  logical:: VERBOSE
  parameter (MAXLEVELS        = 7)                          ! maximum grid level
  parameter (MAXTRIANGLES  = 20*(4**(MAXLEVELS+1) - 1)      ! maximum number of triangles of the spherical grid
  parameter (MAXVERTICES     = 30*(4**MAXLEVELS) + 2        ! maximum number of vertices of the spherical grid

    common/spherical/ &
    &                  vertices, &                         ! vertices coordinates of spherical grid (vertex_index,x,y,z), makes the triangular grid
    &                  cellCorners, &                      ! cell corner vertices coordinates (vertex_index,x,y,z), makes the hexagonal grid
    &                  cellFace,cellNeighbors,&            ! index arrays for a single hexagonal cell (references index of vertices() ); the neighbors of a single vertex;
    &                  cellTriangleFace,&                  ! the vertices of the triangle it belongs to
    &                  numTriangleFaces,           &       ! total number of triangles of the triangular spherical grid
    &                  numCorners,    &                    ! total number of corners (array cellCorners())
    &                  numFaces, numNeighbors,     &       ! total number of faces/cells (array cellFace()); total number of neighbors (array cellNeighbors() )
    &                  TimeParameterSigma,WidthParameterMu,& ! source parameters
    &                  numVertices, &                      ! total number of vertices (of array vertices());
    &                  subdivisions, &                     ! grid subfolding level
    &                  VERBOSE                             ! output verbose
      double precision vertices(MAXVERTICES,3)
      double precision cellCorners(MAXTRIANGLES,3)
      integer  numCorners
      integer cellNeighbors(MAXVERTICES,0:6), numNeighbors
      integer cellFace(MAXVERTICES,0:6), numFaces
      integer cellTriangleFace(MAXTRIANGLES,3), numTriangleFaces
      integer numVertices
      integer subdivisions
      double precision  TimeParameterSigma
      double precision  WidthParameterMu

      save /spherical/

