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
      subroutine readData()
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
      use propagationStartup;use cells;use verbosity;use adjointVariables
      implicit none
      integer::i,k,ioerror,count, iend
      character*1:: divString
      character*23:: fileName
           
      !initialize      
      numVertices = 0      
      write(divString,'(i1)') subdivisions
 100  format(f16.14,2x,f16.14,2x,f16.8)      
 
      ! read all voronoi cell centre vertices values from file 
      !(voronoi cell centers are equal to triangle corners)
      fileName = '../griddata/Dtvert'//divString//'.dat'      
      if( VERBOSE) print*,fileName
      open(1, file= fileName,status='old',iostat=ioerror)
      if( ioerror .ne. 0) call stopProgram( 'abort - readData File1   ')

      !print *,'getting voronoi cell centre vertices from',fileName
      numVertices=0
      do i=1, 30*(4**subdivisions) + 2
        read(1, 100, iostat=ioerror) vertices(i,1), vertices(i,2),vertices(i,3)
       if( ioerror .ne. 0) go to 11
        numVertices = numVertices +1
      enddo     
 11   continue
      if(VERBOSE)print*,'   number of vertices: ',numVertices
      close(1)
      !debug
      if(DEBUG) print*,'vertex',1,(vertices(1,k),k=1,3)
      
      !read in the corresponding cell corners
      fileName = '../griddata/Dvvert'//divString//'.dat'      
      if(VERBOSE)print*,fileName
      open(2,file=fileName,status='old',iostat=ioerror)
      if( ioerror .ne. 0) call stopProgram( 'abort - readData File2   ')

      !print *,'getting voronoi cell corner vertices...'
      numCorners = 0
      do i=1, 20*(4**(subdivisions+1) - 1)
        read(2, 100,iostat=ioerror) (cellCorners(i,k),k=1,3)
        if( ioerror .ne. 0) go to 12
        numCorners = numCorners +1
      enddo     
 12   continue      
      if(VERBOSE)print*,'   number of corners: ',numCorners
      close(2)      
      !debug
      if(DEBUG) print*,'vertex',1,(cellCorners(1,k),k=1,3)
      
      !read in the corresponding cell neighbors indices
      fileName = '../griddata/Dnear'//divString//'.dat'      
      if(VERBOSE)print*,fileName
      open(3,file=fileName,status='old',iostat=ioerror)
      if( ioerror .ne. 0) call stopProgram( 'abort - readData File3   ')

      !print *,'getting voronoi cell neighbors...'
      numNeighbors = 0
      do i=1, 30*(4**subdivisions) + 2
        read(3, '(6i8)') (cellNeighbors(i,k),k=1,6)
        numNeighbors = numNeighbors +1
      enddo     
      if(VERBOSE)print*,'   number of neighbors: ',numNeighbors
      close(3)            
 
      !prepare cellNeighbors info
      do i=1, numNeighbors
        if( cellNeighbors(i,1) .eq. cellNeighbors(i,6)) then
          cellNeighbors(i,0) = 5
        else
          cellNeighbors(i,0) = 6
        endif
      enddo
      !debug
      if(DEBUG) print*,'vertex',1,(cellNeighbors(1,k),k=0,6)
      
      !read in the corresponding cell face indices
      fileName = '../griddata/Dvface'//divString//'.dat'      
      if(VERBOSE)print*,fileName
      open(4,file=fileName,status='old',iostat=ioerror)
      if( ioerror .ne. 0) call stopProgram( 'abort - readData File4   ')

      !print *,'getting voronoi cell face indices...'
      numFaces = 0
      do i=1, 30*(4**subdivisions) + 2
        read(4, '(6i8)') (cellFace(i,k),k=1,6)
        numFaces = numFaces +1
      enddo     
      if(VERBOSE)print*,'   number of faces: ',numFaces
      close(4)                 

      !prepare cellNeighbors info
      do i=1, numFaces
        if( cellFace(i,1) .eq. cellFace(i,6)) then
          cellFace(i,0) = 5
        else
          cellFace(i,0) = 6
        endif
      enddo
      !debug
      if(DEBUG) print*,'cellFace',1,(cellFace(1,k),k=0,6)
                
      ! allocate triangle face array
      if( Station_Correction ) then                
        !read in the corresponding triangle face indices
        fileName = '../griddata/Dtface'//divString//'.dat'      
        if( VERBOSE ) print*,fileName
        open(5,file=fileName,status='old', iostat=ioerror)
        if( ioerror .ne. 0) call stopProgram( 'abort - readData Dtface file    ')

        numTriangleFaces = 0
        do i=1, 20*(4**(subdivisions+1) - 1) 
          if( numTriangleFaces .eq. MaxTriangles) call stopProgram( 'abort-readData triangles    ')
        
          read(5, '(6i8)',iostat=ioerror) (cellTriangleFace(i,k),k=1,3)
          if( ioerror .ne. 0) go to 13         
          numTriangleFaces = numTriangleFaces +1
        enddo     
 13     continue      
        if(VERBOSE) print*,'   number of triangle faces: ',numTriangleFaces
        close(5)                  
      endif
      
      end



