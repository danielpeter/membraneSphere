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
      subroutine syncParameters(rankid,nproc)
!-----------------------------------------------------------------------
! synchronizes parameters from 'Parameter_Input' 
      use propagationStartup;use parallel;use phaseBlockData
      use loop;use verbosity;use cells
      implicit none
      integer,intent(in):: rankid,nproc
      integer:: n,ierror
    
      ! check
      if( nproc < 2 ) return

      ! master process sends parameters to all other processes
      if(rankid .eq. 0) then
        ! send to slaves
        do n=1, nproc-1          
          call MPI_Send(VERBOSE,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ierror)
          call MPI_Send(subdivisions,1,MPI_INTEGER,n,n,MPI_COMM_WORLD,ierror)
          call MPI_Send(FIRSTTIME,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(LASTTIME,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(cphaseRef,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(HETEROGENEOUS,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(DELTA,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(DELTARADIUS,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(deltaLat,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(deltaLon,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(DELTAFUNCTION,8,MPI_CHARACTER,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(sourceLat,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(sourceLon,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(receiverLat,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(receiverLon,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(deltaMoveIncrement,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ierror)               
          call MPI_Send(manyReceivers,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(numofReceivers,1,MPI_INTEGER,n,n,MPI_COMM_WORLD,ierror)          
          call MPI_Send(manyKernels,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(kernelStartDistance,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(kernelEndDistance,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(importKernelsReceivers,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(SIMULATIONOUTPUT,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(cphasetype,8,MPI_CHARACTER,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(MOVEDELTA,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(SECONDDELTA,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(latitudeStart,1,MPI_INTEGER,n,n,MPI_COMM_WORLD,ierror)
          call MPI_Send(latitudeEnd,1,MPI_INTEGER,n,n,MPI_COMM_WORLD,ierror)
          call MPI_Send(longitudeEnd,1,MPI_INTEGER,n,n,MPI_COMM_WORLD,ierror)
          call MPI_Send(deltaPerturbation,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(datadirectory,len(datadirectory),MPI_CHARACTER,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(PARALLELSEISMO,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(phaseBlockFile,len(phaseBlockFile),MPI_CHARACTER,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(phaseBlockVelocityReference,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(heterogeneousDataFile,len(heterogeneousDataFile),MPI_CHARACTER,n,n,MPI_COMM_WORLD,ierror)     
          call MPI_Send(heterogeneousPixelsize,1,MPI_REAL,n,n,MPI_COMM_WORLD,ierror)
          call MPI_Send(gsh_maximum_expansion,1,MPI_INTEGER,n,n,MPI_COMM_WORLD,ierror)
        enddo
      else
          ! get parameters from master
          call MPI_Recv(VERBOSE,1,MPI_LOGICAL,0,rankid,MPI_COMM_WORLD,status,ierror)
          call MPI_Recv(subdivisions,1,MPI_INTEGER,0,rankid,MPI_COMM_WORLD,status,ierror)
          call MPI_Recv(FIRSTTIME,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(LASTTIME,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(cphaseRef,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(HETEROGENEOUS,1,MPI_LOGICAL,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(DELTA,1,MPI_LOGICAL,0,rank,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(DELTARADIUS,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(deltaLat,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(deltaLon,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(DELTAFUNCTION,8,MPI_CHARACTER,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(sourceLat,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(sourceLon,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(receiverLat,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(receiverLon,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(deltaMoveIncrement,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ierror)               
          call MPI_Recv(manyReceivers,1,MPI_LOGICAL,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(numofReceivers,1,MPI_INTEGER,0,rankid,MPI_COMM_WORLD,status,ierror)          
          call MPI_Recv(manyKernels,1,MPI_LOGICAL,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(kernelStartDistance,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(kernelEndDistance,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(importKernelsReceivers,1,MPI_LOGICAL,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(SIMULATIONOUTPUT,1,MPI_LOGICAL,0,rankid,MPI_COMM_WORLD,status,ierror)             
          call MPI_Recv(cphasetype,8,MPI_CHARACTER,0,rankid,MPI_COMM_WORLD,status,ierror)             
          call MPI_Recv(MOVEDELTA,1,MPI_LOGICAL,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(SECONDDELTA,1,MPI_LOGICAL,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(latitudeStart,1,MPI_INTEGER,0,rankid,MPI_COMM_WORLD,status,ierror)
          call MPI_Recv(latitudeEnd,1,MPI_INTEGER,0,rankid,MPI_COMM_WORLD,status,ierror)
          call MPI_Recv(longitudeEnd,1,MPI_INTEGER,0,rankid,MPI_COMM_WORLD,status,ierror)
          call MPI_Recv(deltaPerturbation,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(datadirectory,len(datadirectory),MPI_CHARACTER,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(PARALLELSEISMO,1,MPI_LOGICAL,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(phaseBlockFile,len(phaseBlockFile),MPI_CHARACTER,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(phaseBlockVelocityReference,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(heterogeneousDataFile,len(heterogeneousDataFile),MPI_CHARACTER,0,rankid,MPI_COMM_WORLD,status,ierror)     
          call MPI_Recv(heterogeneousPixelsize,1,MPI_REAL,0,rankid,MPI_COMM_WORLD,status,ierror)
          call MPI_Recv(gsh_maximum_expansion,1,MPI_INTEGER,0,rankid,MPI_COMM_WORLD,status,ierror)          
      endif      
      end
      
!-----------------------------------------------------------------------
      subroutine syncInitialData()
!-----------------------------------------------------------------------
! reads in grid point on master machine and synchronizes with
! its slave processes
      use cells; use parallel !; use propagationStartup; use adjointVariables; use verbosity
      use propagationStartup,only: SIMULATIONOUTPUT
      implicit none
      integer:: n,k,length,ierror

      ! read vertices values from files
      if( MASTER) then
        ! send to slaves
        do n=1, nprocesses-1          
          ! vertices
          length=MaxVertices*3
          tag = n
          call MPI_Send(vertices,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ierror)          
          if( ierror .ne. 0 ) then 
              print*,'syncInitial - vertices',n,' ierror',ierror
              call stopProgram( 'abort - syncInitialData    ')
          endif
          ! numVertices
          length=1
          call MPI_Send(numVertices,length,MPI_INTEGER,n,tag,MPI_COMM_WORLD,ierror)    
     
          ! cellNeighbors
          length=MaxVertices*7
          call MPI_Send(cellNeighbors,length,MPI_INTEGER,n,tag,MPI_COMM_WORLD,ierror)
          if( ierror .ne. 0 ) then 
              print*,'syncInitial - cellNeighbors',n,' ierror',ierror
              call stopProgram( 'abort - syncInitialData    ')
          endif
          ! numNeighbors
          length=1
          call MPI_Send(numNeighbors,length,MPI_INTEGER,n,tag,MPI_COMM_WORLD,ierror)
        
          ! cellCorners          
          length=MaxTriangles*3
          call MPI_Send(cellCorners,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ierror)
          if( ierror .ne. 0 ) then 
              print*,'syncInitial - cellCorners',n,' ierror',ierror
              call stopProgram( 'abort - syncInitialData    ')
          endif

          ! numCorners
          length=1
          call MPI_Send(numCorners,length,MPI_INTEGER,n,tag,MPI_COMM_WORLD,ierror)
          
          ! cellFace
          length=MaxVertices*7
          call MPI_Send(cellFace,length,MPI_INTEGER,n,tag,MPI_COMM_WORLD,ierror)
          if( ierror .ne. 0 ) then 
              print*,'syncInitial - cellFaces',n,' ierror',ierror
              call stopProgram( 'abort - syncInitialData    ')
          endif

          ! numFaces
          length=1
          call MPI_Send(numFaces,length,MPI_INTEGER,n,tag,MPI_COMM_WORLD,ierror)
          
          if( Station_Correction .or. SIMULATIONOUTPUT ) then
            ! cellTriangleFace
            length=MaxTriangles*3
            call MPI_Send(cellTriangleFace,length,MPI_INTEGER,n,tag,MPI_COMM_WORLD,ierror)
            if( ierror .ne. 0 ) then 
              print*,'syncInitial - cellTriangleFaces',n,' ierror',ierror
              call stopProgram( 'abort - syncInitialData    ')
            endif

            ! numFaces
            length=1
            call MPI_Send(numTriangleFaces,length,MPI_INTEGER,n,tag,MPI_COMM_WORLD,ierror)                    
          endif          
       enddo
      else
        ! receive from master       
        ! vertices
        !print*,'slave ',rank,' trying to receive vert. from master'
        length=MaxVertices*3
        tag = rank
        call MPI_RECV(vertices,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ierror)
        if( ierror .ne. 0 ) then 
              print*,'syncInitial - vertices',rank,' ierror',ierror
              call stopProgram( 'abort - syncInitialData    ')
        endif
        ! numVertices
        length=1
        call MPI_RECV(numVertices,length,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierror)     
        
        ! cellNeighbors
        length=MaxVertices*7
        call MPI_RECV(cellNeighbors,length,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierror)
        if( ierror .ne. 0 ) then 
              print*,'syncInitial - cellNeighbors',n,' ierror',ierror
              call stopProgram( 'abort - syncInitialData    ')
        endif
        ! numNeighbors
        length=1
        call MPI_RECV(numNeighbors,length,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierror)

        ! cellCorners
        length=MaxTriangles*3
        call MPI_RECV(cellCorners,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ierror)
        if( ierror .ne. 0 ) then 
              print*,'syncInitial - cellCorners',n,' ierror',ierror
              call stopProgram( 'abort - syncInitialData    ')
        endif
        ! numCorners
        length=1
        call MPI_RECV(numCorners,length,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierror)
        
        ! cellFace
        length=MaxVertices*7
        call MPI_RECV(cellFace,length,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierror)
        if( ierror .ne. 0 ) then 
              print*,'syncInitial - cellFace',n,' ierror',ierror
              call stopProgram( 'abort - syncInitialData    ')
        endif
        ! numFaces
        length=1
        call MPI_RECV(numFaces,length,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierror)     

        if( Station_Correction .or. SIMULATIONOUTPUT ) then
          ! cellTriangleFace
          length=MaxTriangles*3
          call MPI_RECV(cellTriangleFace,length,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierror)
          if( ierror .ne. 0 ) then 
            print*,'syncInitial - cellTriangleFace',n,' ierror',ierror
            call stopProgram( 'abort - syncInitialData    ')
          endif
          ! numFaces
          length=1
          call MPI_RECV(numTriangleFaces,length,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierror)             
        endif
       endif 
      end subroutine
      

!-----------------------------------------------------------------------
      subroutine syncDisplacement( myrank, numberofprocesses )
!-----------------------------------------------------------------------
! synchronizes displacement and displacement_old arrays with
! its slave processes
      use displacements;use propagationStartup
      use verbosity; use parallel; use cells
      implicit none
      integer:: myrank,length,n,numberofprocesses,ierror
      
      if( rank .eq. 0) then
        ! send to slaves
        do n=1, nprocesses-1
          ! displacement
          length=numVertices
          tag = n
          call MPI_Send(displacement,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ierror)
          ! displacement_old
          length=numVertices
          tag = n
          call MPI_Send(displacement_old,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ierror)        
        enddo
      else
        ! receive from master
        ! displacement
        length=numVertices
        tag = rank
        call MPI_RECV(displacement(1:numVertices),length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ierror)
        ! displacement_old
        length=numVertices
        tag = rank
        call MPI_RECV(displacement_old(1:numVertices),length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ierror)
      endif

      !debug
      !  if(rank .eq. 0) then
      !    print*
      !    print*,'master'
      !    do n=1340,1350
      !      print*,'  #',displacement(n)
      !    enddo
      !  else
      !    print*
      !    print*,'slave',rank
      !    do n=1340,1350
      !      print*,'  #',displacement(n)
      !    enddo      
      !  endif
      end
      
!-----------------------------------------------------------------------
      subroutine syncNewdisplacement()
!-----------------------------------------------------------------------
! synchronizes newdisplacement array 
! each process gets only the needed informations from the other ones
      use displacements;use griddomain;use parallel;use verbosity
      implicit none
      integer:: startLocationVertex, endLocationVertex,recStartVertex,recEndVertex
      integer:: dest,source,sendRange,recRange,n,k,neighbor,ierror
                            
      ! check if many processes
      if( nprocesses .eq. 1) return

      !vertex range
      sendRange = boundariesMaxRange
      recRange=sendRange
      if( sendRange .eq. 0 ) then
        ! no boundaries
        return
      endif
      
      !init
      !sendDisp(:)=0.0_WP
      !receiveDisp(:)=0.0_WP      
                        
      !get boundary information from all neighbors
      do n = 1, nprocesses - 1
        ! only talk to other processes/neighbors
        neighbor = domainNeighbors(rank,n)
        if( neighbor .ne. -1 .and. neighbor .ne. rank) then
          !fill send info (over a fixed range of values)
          do k=1,sendRange
            if( boundaries(rank,neighbor,k) .gt. 0) then
              sendDisp(k) = newdisplacement(boundaries(rank,neighbor,k))
            endif
          enddo

          ! talk to neighbor process  (destination to send and source to receive from are the same process)                  
          dest = neighbor
          source = neighbor
          call MPI_SendRecv(sendDisp,sendRange,MPI_CUSTOM,dest,rank,receiveDisp,recRange,MPI_CUSTOM,source, &
                           source,MPI_COMM_WORLD, status, ierror)
          if( ierror .ne. 0 ) then 
              print*,'syncNewdisplacement - p',rank,' ierror',ierror
          endif
          !these are the vertices we received
          do k=1,recRange
            ! check if entry is a valid vertex index in the neighbors boundary array
            if( boundaries(source,rank,k) .gt. 0) then
              ! add as new displacement value for this vertex
              newdisplacement(boundaries(source,rank,k))=receiveDisp(k)
            endif
          enddo
        endif
      enddo

      !debug
      !  print*,'      sync done, barrier to wait.',rank        
      !  ! wait until all processes reached this point
      !  call MPI_Barrier( MPI_COMM_WORLD, ierror )
      !  if( ierror .ne. 0) call stopProgram('abort - MPI_Barrier kernels failed    ')            
      end
      
!-----------------------------------------------------------------------
      subroutine syncNewdisplacement_alt( rank, nprocesses )
!-----------------------------------------------------------------------
! synchronizes newdisplacement array with
! its slave processes
      use displacements; use propagationStartup; use cells
      implicit none
      include 'mpif.h'
      integer:: rank,ierror,tag,n,nprocesses,status(MPI_STATUS_SIZE)
      integer:: startLocationVertex, endLocationVertex,recStartVertex,recEndVertex
      double precision, allocatable, dimension(:) :: tmpSend,tmpReceiv
      integer:: dest,source,range,vertexpart,sendRange,recRange,rangeMax
      
      !these are the vertices which are calculated by this process
      !startLocationVertex = rank*range+1
      !endLocationVertex = (rank+1)*range     
      ! determine the first and last vertice to send
      vertexpart=int(numVertices/nprocesses)
      startLocationVertex= rank*vertexpart + 1
      endLocationVertex = (rank+1)*vertexpart
      ! add remaining vertices to the last process
      if( rank+1 .eq. nprocesses) then
        endLocationVertex= endLocationVertex + mod(numVertices,nprocesses)
      endif

      !check vertex range
      sendRange = endLocationVertex - startLocationVertex + 1
      if( sendRange .lt. 1) call stopProgram( 'abort-syncNewdisplacement short range    ')  

      !allocate array for the data to send
      allocate( tmpSend(sendRange), stat=ierror)
      if( ierror .gt. 0 ) then
        print*,rank,' error in allocating send array'
        call stopProgram( 'abort - syncNewdisplacement_alt    ')
      endif
      tmpSend(1:sendRange) = newdisplacement(startLocationVertex:endLocationVertex)
      
      !allocate array for the data to receive
      rangeMax=0
      do n=0, nprocesses -1
        ! determine the first and last vertice to receive
        recStartVertex= n*int(numVertices/nprocesses) + 1
        recEndVertex = (n+1)*int(numVertices/nprocesses)
        if( n+1 .eq. nprocesses) then  ! add remaining vertices to the last process
          recEndVertex= recEndVertex + mod(numVertices,nprocesses)
        endif
        range=recEndVertex-recStartVertex+1
        if( range .gt. rangeMax) rangeMax=range
      enddo
      ! prepare array to receive data
      allocate( tmpReceiv(rangeMax), stat=ierror)
      if( ierror .gt. 0 ) then
        print*,rank,' error in allocating receiving array ',rangeMax
        call stopProgram( 'abort - syncNewdisplacement_alt    ')
      endif

      !print*,rank,'sendreceiving...'
      do n = 0, nprocesses - 1
        ! only talk to other processes
        if( n .ne. rank) then
          ! determine the first and last vertice to receive
          recStartVertex= n*int(numVertices/nprocesses) + 1
          recEndVertex = (n+1)*int(numVertices/nprocesses)
          if( n+1 .eq. nprocesses) then  ! add remaining vertices to the last process
            recEndVertex= recEndVertex + mod(numVertices,nprocesses)
          endif
          recRange=recEndVertex-recStartVertex+1

          ! talk to neighbor process                    
          dest = n
          source = n
          call MPI_SendRecv( tmpSend,sendRange, MPI_DOUBLE_PRECISION,   &
     &    dest,rank,tmpReceiv,recRange,MPI_DOUBLE_PRECISION,source,     &
     &    source,MPI_COMM_WORLD, status, ierror)
          if( ierror .ne. 0 ) then 
              print*,'syncNewdisplacement - p',rank,' ierror',ierror
          endif
          !these are the vertices we received
          !startLocationVertex = source*range+1
          !endLocationVertex = (source+1)*range
          newdisplacement(recStartVertex:recEndVertex)=tmpReceiv(1:recRange)
          
        endif
      enddo
      end
      
!-----------------------------------------------------------------------
      subroutine syncNewdisplacement_alt2(rank, nprocesses)
!-----------------------------------------------------------------------
! synchronizes newdisplacement array with
! its slave processes
      use displacements; use propagationStartup; use cells
      implicit none
      include 'mpif.h'
      integer:: rank, ierror, tag, length, n, nprocesses, i,status(MPI_STATUS_SIZE)
      integer:: startLocationVertex, endLocationVertex
      double precision:: tmp(numVertices/nprocesses)
      double precision:: tmpSend(numVertices/nprocesses)
      double precision:: tmpReceiv(numVertices/nprocesses)
      integer:: sendcount, sendtag, sendtype,dest, source,recvtag,comm,recvtype, recvcount,range
      
      range = numVertices/nprocesses
      if( range .lt. 1) call stopProgram( 'abort-syncNewdisplacement short range    ')                                                                                       
                                                                                                                                                                                                                                                
      if( rank .eq. 0) then
          ! receive from slaves
          !print*,'sync newdisplacement' !debug
          do n = 1, nprocesses-1
            !print*,'receiving tmpdisp. from slave',n !debug
            ! new displacements          

                        
!            length=numVertices/nprocesses
!            tag = n
!            call MPI_Recv(tmp,length,                       &
!     &                          MPI_DOUBLE_PRECISION,n,tag,             &
!     &                          MPI_COMM_WORLD,status,ierror)
            !print*,'master ', range
            
            tmpSend(1:range) = newdisplacement(1:range)
            tmpReceiv = 0.0
            
            sendcount = range
            sendtype = MPI_DOUBLE_PRECISION
            dest = n
            sendtag = 0
            recvcount = range
            recvtype = MPI_DOUBLE_PRECISION
            source = n
            recvtag = n
            comm = MPI_COMM_WORLD
            
            !print*,'master sendreceiving...'
            
            call MPI_SendRecv( tmpSend,sendcount, sendtype, dest,       &
     &        sendtag, tmpReceiv, recvcount, recvtype, source, recvtag,   &
     &        comm, status, ierror)
            if( ierror .ne. 0 ) then 
              print*,'syncNewdisplacement - master ierror',ierror
              call stopProgram( 'abort - syncNewdisplacement_alt2   ')              
            endif
       
            ! fill in corresponding values
            startLocationVertex = n*range + 1
            endLocationVertex = (n+1)*range
            !print*,'master received', startLocationVertex, endLocationVertex
            
            newdisplacement(startLocationVertex:endLocationVertex)=tmpReceiv(1:range)
            
            !do i=startLocationVertex, endLocationVertex
            !  newdisplacement(i) = tmpReceiv(i)
            !enddo
!            if( (numVertices/nprocesses) .ne.                           &
!     &          (endLocationVertex-startLocationVertex+1) ) then
!              stop 'syncNewdisplacement - wrong array size'
!            endif
              
!            newdisplacement(startLocationVertex:endLocationVertex) =                    &
!     &                   tmp(1:numVertices/nprocesses)
          enddo
     
          ! send updated whole array to slaves
          
!          do n=1, nprocesses-1
!            !print*,'sending newdisp. to slave',n !debug
!            length=numVertices
!            tag = n
!            call MPI_Send(newdisplacement,length,MPI_DOUBLE_PRECISION,  &
!     &                  n,tag,MPI_COMM_WORLD,ierror)        
!          enddo
      else
          ! send to master
          !print*,'slave',rank,'sending newdisp. to master' !debug
          ! new displacements
!          startLocationVertex = rank*numVertices/nprocesses + 1
!          endLocationVertex = (rank+1)*numVertices/nprocesses
!          if( (numVertices/nprocesses) .ne.                             &
!     &      (endLocationVertex-startLocationVertex+1) ) then
!              stop 'syncNewdisplacement - wrong array size'
!          endif
          
!          tmp(1:numVertices/nprocesses)=                                &
!     &              newdisplacement(startLocationVertex:endLocationVertex)
!          length=numVertices/nprocesses
!          tag = rank
!          call MPI_Send(tmp,length,MPI_DOUBLE_PRECISION,    &
!     &                0,tag,MPI_COMM_WORLD,ierror)
     
          ! receive from master new updated array
          !print*,'slave',rank,'receiving newdisp. from master'
!          length=numVertices
!          tag = rank
!          call MPI_Recv(newdisplacement,length,MPI_DOUBLE_PRECISION,    &
!     &                0,tag,MPI_COMM_WORLD,status,ierror)
     
            startLocationVertex = rank*range + 1
            endLocationVertex = (rank+1)*range     
            !print*,'slave ', startLocationVertex, endLocationVertex
            
            tmpSend(1:range) = newdisplacement(startLocationVertex:endLocationVertex)
            tmpReceiv = 0.0
            
            sendcount = range
            sendtype = MPI_DOUBLE_PRECISION
            dest = 0
            sendtag = rank
            recvcount = range
            recvtype = MPI_DOUBLE_PRECISION
            source = 0
            recvtag = 0
            comm = MPI_COMM_WORLD
            
            !print*,'slave', rank,' sendreceiving...'
            
            call MPI_SendRecv( tmpSend,sendcount, sendtype, dest,       &
     &        sendtag, tmpReceiv, recvcount, recvtype, source, recvtag,   &
     &        comm, status, ierror)
            if( ierror .ne. 0 ) then 
              print*,'syncNewdisplacement - slave', rank,'ierror',ierror
              call stopProgram( 'abort - syncNewdisplacement_alt2   ' )             
            endif
       
            ! fill in corresponding values
            !print*,'slave received', range
            newdisplacement(1:range)=tmpReceiv(1:range)
     
      endif
  
      !debug
      !if(rank .eq. 0) then
      !  print*
      !  print*,'new master'
      !  do l=1340,1350
      !    print*,'  m#',newdisplacement(l)
      !  enddo
      !open(11,file='syncNew0.dat')
      !  write(11,*) newdisplacement
      !  close(11)
      !else
      !  print*
      !  print*,'new slave',rank
      !  do l=1340,1350
      !    print*,'  #',newdisplacement(l)
      !  enddo        
      !  open(12,file='syncNew1.dat')
      !  write(12,*) newdisplacement
      ! close(12)
      !endif
        
      end subroutine

!-----------------------------------------------------------------------
      subroutine syncPrecalculated()
!-----------------------------------------------------------------------
! synchronize precalculated array cellAreas, cellEdgesLength, cellCenterDistances      
      use cells; use parallel; use propagationStartup; use verbosity
      implicit none
      integer:: length,n,first,lastsent,i,ierror
      real(WP), allocatable, dimension(:):: tmpExchange
      integer,parameter::SLICE=1000000 ! assuming MPI has a limited buffer size
        
      ! debug check arrays
      !  ! arrays
      !  if( .not. allocated(cellAreas) ) then
      !    print*,'cellAreas not allocated',rank
      !    call stopProgram("array not allocated    ")
      !  endif
      !  if( .not. allocated(cellEdgesLength) ) then
      !    print*,'cellEdgesLength not allocated',rank
      !    call stopProgram("array not allocated    ")
      !  endif
      !  if( .not. allocated(cellCenterDistances) ) then
      !    print*,'cellCenterDistances not allocated',rank
      !    call stopProgram("array not allocated    ")
      !  endif  
      !  ! barrier
      !  call MPI_Barrier(MPI_COMM_WORLD,ierror)
      !  if( ierror .ne. 0) call stopProgram('syncRoutines - MPI_Barrier kernels failed    ')              
            
      if( MASTER ) then
        !read in values
        call readPrecalculated()
        
        ! send to slaves
        do n=1, nprocesses-1                    
          tag = n        
          
          ! number of vertices
          length=numVertices
          
          ! split into chunks to avoid shared memory allocation problems of MPI
          if( length .gt. SLICE ) then
            lastsent=0
            do while(lastsent .lt. numVertices )
              first=lastsent+1
              lastsent=lastsent+SLICE
              if( lastsent .ge. numVertices ) lastsent=numVertices
              length=lastsent-first+1                            
              call MPI_Send(cellAreas(first:lastsent),length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ierror)              
            enddo
          else
            call MPI_Send(cellAreas,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ierror)          
          endif  

          length=numVertices
          ! split into chunks to avoid shared memory allocation problems of MPI
          if( length .gt. SLICE ) then
            lastsent=0
            do while(lastsent .lt. numVertices )
              first=lastsent+1
              lastsent=lastsent+SLICE
              if( lastsent .ge. numVertices ) lastsent=numVertices
              length=(lastsent-first+1)
              allocate(tmpExchange(length),stat=ierror)
              if( ierror .ne. 0 ) call stopProgram('could not allocate temporary exchange array!      ')
              
              do i=0,6
                tmpExchange(:)=cellEdgesLength(first:lastsent,i)              
                call MPI_Send(tmpExchange,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ierror)

                tmpExchange(:)=cellCenterDistances(first:lastsent,i)                                
                call MPI_Send(tmpExchange,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ierror)  
                
                if( CORRECT_RATIO ) then
                  tmpExchange(:)=cellFractions(first:lastsent,i)                
                  call MPI_Send(tmpExchange,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ierror)                    
                endif                
              enddo
              deallocate(tmpExchange)
            enddo
          else
            length=numVertices*7
            call MPI_Send(cellEdgesLength,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ierror)            
            call MPI_Send(cellCenterDistances,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ierror)       
            if( CORRECT_RATIO ) then
              call MPI_Send(cellFractions,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ierror)                   
            endif            
          endif  
       enddo       
      else
        ! receive from master
        tag = rank
        length=numVertices
        ! split into chunks to avoid shared memory allocation problems of MPI
        if( length .gt. SLICE ) then
          lastsent=0
          do while(lastsent .lt. numVertices )
            first=lastsent+1
            lastsent=lastsent+SLICE
            if( lastsent .ge. numVertices ) lastsent=numVertices
            length=lastsent-first+1
            call MPI_RECV(cellAreas(first:lastsent),length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ierror)
          enddo
        else
          call MPI_RECV(cellAreas,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ierror)
        endif          

        length=numVertices
        ! split into chunks to avoid shared memory allocation problems of MPI
        if( length .gt. SLICE ) then
          lastsent=0
          do while(lastsent .lt. numVertices )
            first=lastsent+1
            lastsent=lastsent+SLICE
            if( lastsent .ge. numVertices ) lastsent=numVertices
            length=(lastsent-first+1)
            allocate(tmpExchange(length),stat=ierror)
            if( ierror .ne. 0 ) call stopProgram('could not allocate temporary exchange array in slave!      ')
            
            do i=0,6            
              call MPI_RECV(tmpExchange,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ierror)
              cellEdgesLength(first:lastsent,i)=tmpExchange(:)
              
              call MPI_RECV(tmpExchange,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ierror)     
              cellCenterDistances(first:lastsent,i)=tmpExchange(:)              
              
              if( CORRECT_RATIO ) then
                call MPI_RECV(tmpExchange,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ierror)     
                cellFractions(first:lastsent,i)=tmpExchange(:)                            
              endif              
            enddo
            deallocate(tmpExchange)
          enddo
        else
          length=numVertices*7          
          call MPI_RECV(cellEdgesLength,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ierror)          
          call MPI_RECV(cellCenterDistances,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ierror)               
          if( CORRECT_RATIO ) then
            call MPI_RECV(cellFractions,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ierror)               
          endif
        endif          
      endif
      end subroutine

!-----------------------------------------------------------------------
      subroutine syncPhaseMap(myrank,numberofprocesses)
!-----------------------------------------------------------------------
! synchronize phase map
      use cells;use phaseVelocityMap;use propagationStartup;use verbosity;use parallel
      implicit none
      integer:: myrank,length,n,numberofprocesses,ierror
      
      if(rank .eq. 0) then        
        ! send to slaves
        do n=1, nprocesses-1          
          ! vertices
          length=numVertices
          tag = n
          
          call MPI_Send(phaseMap,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ierror)
       enddo
      else
        ! receive from master
        length=numVertices
        tag = rank        
        call MPI_RECV(phaseMap,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ierror)
      endif
      end


!-----------------------------------------------------------------------
      subroutine syncReceiversRef()
!-----------------------------------------------------------------------
! synchronize array receiversSeismogramRef     
      use cells; use parallel; use propagationStartup; use verbosity
      implicit none
      integer:: length,n,m,ierror
      
      if( MASTER ) then        
        ! send to slaves
        do n=1, nprocesses-1          
          ! vertices
          if( manyKernels ) then
            do m=1,numofKernels
              length=size(kernelsReceiversSeismogramRef(m,:,:))
              tag = n              
              call MPI_Send(kernelsReceiversSeismogramRef(m,:,:),length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ierror)
            enddo
          else
            length=size(receiversSeismogramRef)
            tag = n
            call MPI_Send(receiversSeismogramRef,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ierror)          
          endif
       enddo
      else
        if( manyKernels ) then
          do m=1,numofKernels
            ! receive from master
            length=size(kernelsReceiversSeismogramRef(m,:,:))
            tag = rank            
            call MPI_Recv(kernelsReceiversSeismogramRef(m,:,:),length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ierror)
          enddo
        else
          ! receive from master
          length=size(receiversSeismogramRef)
          tag = rank
          call MPI_Recv(receiversSeismogramRef,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ierror)
        endif
      endif
      end

!-----------------------------------------------------------------------
      subroutine collectReceiversSeismogram(doReference)
!-----------------------------------------------------------------------
! collects array receiversSeismogram from receivers not in the master domain
      use cells; use parallel; use propagationStartup;use griddomain; use verbosity
      implicit none
      integer:: i,countVertices,recRange,sendRange,n,m,ierror
      real(WP),allocatable,dimension(:,:)::tmpExchange
!      real(WP),allocatable,dimension(:,:,:)::tmpManyExchange
      logical:: doReference
      
      ! check if many processes
      if( nprocesses .eq. 1 ) return
      
      if( MASTER ) then        
        ! debug
        !print*,'collectReceiversSeismogram() -master',nprocesses,rank,MASTER,numofReceivers,doReference
      
        ! get from the other processes
        do n=1, nprocesses-1          
          !vertex range
          countVertices=0
          if( manyKernels ) then
          !  countVertices = numofReceivers
            continue
          else
            do i=1,numofReceivers
              if( vertexDomain(receivers(i)) .eq. n) then              
                countVertices=countVertices+1
              endif
            enddo          
          endif
          recRange = countVertices
          
          ! allocate array to receive
          if( manyKernels ) then
            !allocate(tmpManyExchange(numofKernels,recRange,numofTimeSteps),stat=ierror)                      
            !allocate(tmpExchange(recRange,numofTimeSteps),stat=ierror)            
            continue
          else
            allocate(tmpExchange(recRange,numofTimeSteps),stat=ierror)
            if( ierror .ne. 0 ) then
              print*,'master collects receivers seismogram:',numofKernels,recRange,numofTimeSteps
              call stopProgram('collectReceiversSeismogram() - error allocating receive array   ')
            endif            
          endif

          ! get info from slave
          tag = n
          if( manyKernels ) then
            !call MPI_Recv(tmpManyExchange,size(tmpManyExchange),MPI_CUSTOM,n,tag,MPI_COMM_WORLD,status,ierror)          
            do m=1,numofKernels
              countVertices = 0
              do i=1,numofReceivers
                if( vertexDomain(KernelsReceivers(m,i)) .eq. n) then              
                  countVertices=countVertices+1
                endif
              enddo                    
              recRange = countVertices
              
              allocate(tmpExchange(recRange,numofTimeSteps),stat=ierror)                          
              if( ierror .ne. 0 ) then
                print*,'master collects receivers seismogram:',m,numofKernels,recRange,numofTimeSteps
                call stopProgram('collectReceiversSeismogram() - error allocating receive array   ')
              endif
            
              call MPI_Recv(tmpExchange,size(tmpExchange),MPI_CUSTOM,n,tag,MPI_COMM_WORLD,status,ierror)      
              if( ierror .ne. 0 ) then 
                print*,'collectReceiversSeismogram - ',rank,n,ierror
                call stopProgram('collectReceiversSeismogram() - error mpi-recv   ')
              endif
              
              ! these are the vertices we received              
              countVertices = 0
              do i=1,numofReceivers
                ! extract only seismograms from receiver stations in that particular process domain                
                if( vertexDomain(kernelsReceivers(m,i)) .eq. n) then
                  countVertices = countVertices + 1
                  if( doReference) then
                    kernelsReceiversSeismogramRef(m,i,:)=tmpExchange(countVertices,:)
                  else
                    kernelsReceiversSeismogram(m,i,:)=tmpExchange(countVertices,:)                    
                  endif                
                endif
              enddo              
              
              deallocate( tmpExchange )
            enddo
          else
            call MPI_Recv(tmpExchange,size(tmpExchange),MPI_CUSTOM,n,tag,MPI_COMM_WORLD,status,ierror)
            if( ierror .ne. 0 ) then 
              print*,'collectReceiversSeismogram - ',rank,n,ierror
              call stopProgram('collectReceiversSeismogram() - error mpi-recv   ')
            endif            
          endif

          ! these are the vertices we received
          countVertices=0
          do i=1,numofReceivers
            ! check if entry is a valid vertex index in the neighbors boundary array
            if( manyKernels ) then
!              do m=1,numofKernels
!                ! extract only seismograms from receiver stations in that particular process domain
!                if( vertexDomain(kernelsReceivers(m,i)) .eq. n) then
!                  if( doReference) then
!                    kernelsReceiversSeismogramRef(m,i,:)=tmpManyExchange(m,i,:)
!                  else
!                    kernelsReceiversSeismogram(m,i,:)=tmpManyExchange(m,i,:)                    
!                  endif                
!                endif
!              enddo
              continue
            else
              if( vertexDomain(receivers(i)) .eq. n ) then
                countVertices=countVertices+1
                if( doReference ) then
                  receiversSeismogramRef(i,:)=tmpExchange(countVertices,:)              
                else
                  receiversSeismogram(i,:)=tmpExchange(countVertices,:)                              
                endif
              endif
            endif
          enddo
        enddo
      else
        ! debug
        !print*,'collectReceiversSeismogram() - slave',nprocesses,rank,MASTER,numofReceivers,doReference
        
        !vertex range
        countVertices=0
        if( manyKernels ) then
        !  countVertices=numofReceivers
          continue
        else
          do i=1,numofReceivers
            if( vertexDomain(receivers(i)) .eq. rank) then
              countVertices=countVertices+1
            endif
          enddo          
        endif
        sendRange = countVertices
        
        ! allocate array
        if( manyKernels ) then
          !allocate(tmpManyExchange(numofKernels,sendRange,numofTimeSteps),stat=ierror)
          !allocate(tmpExchange(sendRange,numofTimeSteps),stat=ierror)                  
          continue
        else
          allocate(tmpExchange(sendRange,numofTimeSteps),stat=ierror)        
          if( ierror .ne. 0 ) then
            print*,'slave sends receivers seismograms:',rank,numofKernels,sendRange,numofTimeSteps
            call stopProgram('collectReceiversSeismogram() - error allocating send array    ')
          endif          
        endif
          
        ! fill in corresponding values
        countVertices=0
        do i=1,numofReceivers
          if( manyKernels) then
!            ! fill in all seismograms
!            do m=1,numofKernels
!              if( doReference ) then
!                tmpManyExchange(m,i,:)=kernelsReceiversSeismogramRef(m,i,:)
!              else
!                tmpManyExchange(m,i,:)=kernelsReceiversSeismogram(m,i,:)              
!              endif            
!            enddo
            continue
          else
            ! fill in only seismogram from receivers in this process domain
            if( vertexDomain(receivers(i)) .eq. rank) then          
              countVertices=countVertices+1
              if( doReference ) then
                tmpExchange(countVertices,:)=receiversSeismogramRef(i,:)
              else
                tmpExchange(countVertices,:)=receiversSeismogram(i,:)            
              endif
            endif
          endif
        enddo
        
        ! send to master
        tag = rank
        if( manyKernels ) then
          !call MPI_Send(tmpManyExchange,size(tmpManyExchange),MPI_CUSTOM,0,tag,MPI_COMM_WORLD,ierror)
          do m=1,numofKernels
            countVertices = 0
            do i=1,numofReceivers
              if( vertexDomain(kernelsReceivers(m,i)) .eq. rank) then
                countVertices=countVertices+1
              endif
            enddo          
            sendRange = countVertices
            
            allocate(tmpExchange(sendRange,numofTimeSteps),stat=ierror)                  
            if( ierror .ne. 0 ) then
              print*,'slave sends receivers seismograms:',rank,m,numofKernels,sendRange,numofTimeSteps
              call stopProgram('collectReceiversSeismogram() - error allocating send array    ')
            endif          
          
            ! fill in only seismogram from receivers in this process domain
            countVertices = 0
            do i=1,numofReceivers            
              if( vertexDomain(kernelsReceivers(m,i)) .eq. rank) then          
                countVertices = countVertices + 1               
                if( doReference ) then
                  tmpExchange(countVertices,:)=kernelsReceiversSeismogramRef(m,i,:)
                else
                  tmpExchange(countVertices,:)=kernelsReceiversSeismogram(m,i,:)              
                endif            
              endif
            enddo
            ! send 
            call MPI_Send(tmpExchange,size(tmpExchange),MPI_CUSTOM,0,tag,MPI_COMM_WORLD,ierror)        
            if( ierror .ne. 0 ) then 
              print*,'collectReceiversSeismogram - ',rank,n,ierror
            endif
            
            deallocate( tmpExchange )
          enddo
        else
          call MPI_Send(tmpExchange,size(tmpExchange),MPI_CUSTOM,0,tag,MPI_COMM_WORLD,ierror)        
          if( ierror .ne. 0 ) then 
            print*,'collectReceiversSeismogram - ',rank,n,ierror
          endif          
        endif ! manykernels
      endif ! MASTER
      end

!-----------------------------------------------------------------------
      subroutine syncBackwardNewdisplacement()
!-----------------------------------------------------------------------
! synchronizes backwardnewdisplacement array 
! each process gets the needed (boundary) informations from the other ones
      use displacements;use griddomain;use parallel;use propagationStartup; use adjointVariables
      implicit none
      integer:: dest,source,sendRange,recRange,n,k,neighbor,ierror
      
      ! check if many processes
      if( nprocesses .eq. 1) return

      !vertex range
      sendRange = boundariesMaxRange
      recRange=sendRange
      if( sendRange .eq. 0 ) then
        ! no boundaries
        return
      endif
      
      !init
      sendDisp(:)=0.0_WP
      receiveDisp(:)=0.0_WP      
                        
      !get boundary information from all neighbors
      do n = 1, nprocesses - 1
        ! only talk to other processes/neighbors
        neighbor = domainNeighbors(rank,n)
        if( neighbor .ne. -1 .and. neighbor .ne. rank) then
          !fill send info (over a fixed range of values)
          do k=1,sendRange
            if( boundaries(rank,neighbor,k) .gt. 0) then
              sendDisp(k) = backwardNewdisplacement(boundaries(rank,neighbor,k))
            endif
          enddo

          ! talk to neighbor process  (destination to send and source to receive from are the same process)                  
          dest = neighbor
          source = neighbor
          call MPI_SendRecv( sendDisp,sendRange, MPI_CUSTOM,   &
     &                  dest,rank,receiveDisp,recRange,MPI_CUSTOM,source,     &
     &                  source,MPI_COMM_WORLD, status, ierror)
          if( ierror .ne. 0 ) then 
              print*,'syncbackwardNewdisplacement - p',rank,' ierror',ierror
          endif
          !these are the vertices we received
          do k=1,recRange
            ! check if entry is a valid vertex index in the neighbors boundary array
            if( boundaries(source,rank,k) .gt. 0) then
              ! add as new displacement value for this vertex
              backwardNewdisplacement(boundaries(source,rank,k))=receiveDisp(k)
            endif
          enddo
        endif
      enddo
      end


!-----------------------------------------------------------------------
      subroutine syncReceiverSeismogram()
!-----------------------------------------------------------------------
! let seismogram at receiver be available for all processes
      use parallel;use propagationStartup; use verbosity
      implicit none
      integer:: domain,i,length,ierror
      integer,external::getDomain
      
      ! check if different grid domains exist
      if( nprocesses < 2 ) return      
      
      ! get right seismogram at receiver station
      domain = getDomain(receiverVertex)
      if( domain .eq. rank ) then
        do i=0,nprocesses-1
          if( i .ne. domain ) then
            ! send seismogram 
            length = size(seismogramReceiver(:,:))
            tag = rank          
            call MPI_Send(seismogramReceiver,length,MPI_CUSTOM,i,tag,MPI_COMM_WORLD,ierror)
          endif
        enddo      
      else 
        ! receive from domain process
        length = size(seismogramReceiver(:,:))
        tag = domain          
        call MPI_Recv(seismogramReceiver,length,MPI_CUSTOM,domain,tag,MPI_COMM_WORLD,status,ierror)
      endif            

      end subroutine

      
!-----------------------------------------------------------------------
      subroutine collectAdjointKernel()
!-----------------------------------------------------------------------
! completes adjointKernel array for master
! master collects from all slave processes the needed informations
      use griddomain;use parallel;use propagationStartup 
      use adjointVariables; use cells
      implicit none
      integer:: sendRange,recRange,n,k,i,countVertices,ierror
      real(WP),allocatable,dimension(:)::tmpExchange
      
      ! check if many processes
      if( nprocesses < 2) return

      ! only receive from slaves
      if( MASTER) then
        do n = 1, nprocesses-1
          !vertex range
          countVertices=0
          do i=1,numVertices
            if( vertexDomain(i) .eq. n) then
              countVertices=countVertices+1
            endif
          enddo          
          recRange = countVertices
          
          ! allocate array to receive
          allocate(tmpExchange(recRange),stat=ierror)
          if( ierror .ne. 0 ) then
            print*,'master allocate error:',ierror,recRange,n
            call stopProgram('collectAdjointKernel() - error allocating receive array   ')
          endif
          
          ! get info from slave
          tag = n
          call MPI_Recv(tmpExchange,recRange,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,status,ierror)
          if( ierror .ne. 0 ) then 
              print*,'collectAdjointKernel - ',rank,n,ierror
              call stopProgram('collectAdjointKernel() - error mpi-recv   ')
          endif

          !these are the vertices we received
          countVertices=0
          do i=1,numVertices
            ! check if entry is a valid vertex index in the neighbors boundary array
            if( vertexDomain(i) .eq. n ) then
              countVertices=countVertices+1
              adjointKernel(i)=tmpExchange(countVertices)
            endif
          enddo
          
          ! free again
          deallocate(tmpExchange)

        enddo
      else  
        !vertex range
        sendRange = size(domainVertices(:))
        
        ! allocate array
        allocate(tmpExchange(sendRange),stat=ierror)
        if( ierror .ne. 0 ) call stopProgram('collectAdjointKernel() - error allocating send array    ')
          
        ! fill in corresponding values
        do i=1,sendRange
          tmpExchange(i) = adjointKernel(domainVertices(i))
        enddo
        
        ! send to master
        tag = rank
        call MPI_Send(tmpExchange,sendRange,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,ierror)
        if( ierror .ne. 0 ) then 
            print*,'error collectAdjointKernel - ',rank,n,ierror
        endif
      endif
      end subroutine
      
!-----------------------------------------------------------------------
      subroutine collectFullNewdisplacement()
!-----------------------------------------------------------------------
! completes newdisplacement array for master
! master collects from all slave processes the needed informations
      use displacements; use griddomain; use parallel
      use propagationStartup; use cells
      implicit none
      integer:: sendRange,recRange,n,i,countVertices,ierror
      
      ! check if many processes
      if( nprocesses .eq. 1) return
      
      !init
      sendDisp(:)=0.0_WP
                        
      ! only receive from slaves
      if( MASTER) then
        do n = 1, nprocesses-1
          !vertex range
          countVertices=0
          do i=1,numVertices
            if( vertexDomain(i) .eq. n) then
              countVertices=countVertices+1
            endif
          enddo          
          recRange = countVertices

          ! get info from slave
          receiveDisp(:)=0.0_WP                
          tag = n
          call MPI_Recv(receiveDisp,recRange,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,status,ierror)
          if( ierror .ne. 0 ) then 
              print*,'collectFullNewdisplacement - ',rank,n,ierror
          endif
          !these are the vertices we received
          countVertices=0
          do i=1,numVertices
            ! check if entry is a valid vertex index in the neighbors boundary array
            if( vertexDomain(i) .eq. n ) then
              countVertices=countVertices+1
              newdisplacement(i)=receiveDisp(countVertices)
            endif
          enddo
        enddo
      else  
        !vertex range
        sendRange = size(domainVertices(:))
        do i=1,sendRange
          sendDisp(i)=newdisplacement(domainVertices(i))
        enddo
        
        ! send to master
        tag = rank
        call MPI_Send(sendDisp,sendRange,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,ierror)
        if( ierror .ne. 0 ) then 
            print*,'collectFullNewdisplacement - ',rank,n,ierror
        endif
      endif
      end
      
!-----------------------------------------------------------------------
      subroutine syncInterpolationData()
!-----------------------------------------------------------------------
! synchronizes data related to interpolation with its slave processes
      use cells; use parallel; use propagationStartup; use verbosity
      implicit none
      integer:: n,k,length,ierror

      ! read vertices values from files
      if( MASTER) then
        ! send to slaves
        do n=1, nprocesses-1          
          ! distances
          length=3
          tag = n
          call MPI_Send(interpolation_distances,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ierror)          
          if( ierror .ne. 0 ) then 
              print*,'syncInterpolation - distances',n,' ierror',ierror
              call stopProgram( 'abort - syncInterpolation    ')
          endif
          ! side lengths
          length=3
          tag = n
          call MPI_Send(interpolation_triangleLengths,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ierror)          
          if( ierror .ne. 0 ) then 
              print*,'syncInterpolation - lengths',n,' ierror',ierror
              call stopProgram( 'abort - syncInterpolation    ')
          endif
          ! corner indices
          length=3
          tag = n
          call MPI_Send(interpolation_corners,length,MPI_INTEGER,n,tag,MPI_COMM_WORLD,ierror)          
          if( ierror .ne. 0 ) then 
              print*,'syncInterpolation - corners',n,' ierror',ierror
              call stopProgram( 'abort - syncInterpolation    ')
          endif
       enddo
      else
        ! receive from master       
        ! distances
        length=3
        tag = rank
        call MPI_RECV(interpolation_distances,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ierror)
        if( ierror .ne. 0 ) then 
              print*,'syncInterpolation - distances',rank,' ierror',ierror
              call stopProgram( 'abort - syncInterpolation    ')
        endif
        ! side lengths
        length=3
        tag = rank
        call MPI_RECV(interpolation_triangleLengths,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ierror)
        if( ierror .ne. 0 ) then 
              print*,'syncInterpolation - lengths',rank,' ierror',ierror
              call stopProgram( 'abort - syncInterpolation    ')
        endif
        ! corner indices
        length=3
        tag = rank
        call MPI_RECV(interpolation_corners,length,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierror)
        if( ierror .ne. 0 ) then 
              print*,'syncInterpolation - corners',rank,' ierror',ierror
              call stopProgram( 'abort - syncInterpolation    ')
        endif
      endif 
      end
      
!-----------------------------------------------------------------------
      subroutine syncPixelgridInput(rank,nprocesses,txtfile,filen,eq_incr,ifa)
!-----------------------------------------------------------------------     
! sync with other process
      implicit none
      include 'mpif.h'
      integer:: ierror,nprocesses,rank,i,length,status(MPI_STATUS_SIZE) 
      character:: filen*128,txtfile*128
      real:: eq_incr
      integer:: ifa

      ! string length
      length=128
      
      ! nothing to do for single processor job
      if( nprocesses .lt. 2 ) return
      
      ! synchronize input data
      if( rank .eq. 0 ) then
        ! MASTER sends
        do i=1, nprocesses-1
          call MPI_Send(txtfile,length,MPI_CHARACTER,i,i,MPI_COMM_WORLD,ierror)
          call MPI_Send(filen,length,MPI_CHARACTER,i,i,MPI_COMM_WORLD,ierror)
          call MPI_Send(eq_incr,1,MPI_REAL,i,i,MPI_COMM_WORLD,ierror)
          call MPI_Send(ifa,1,MPI_INTEGER,i,i,MPI_COMM_WORLD,ierror)     
        enddo
      else
        ! slave process receives information
        call MPI_Recv(txtfile,length,MPI_CHARACTER,0,rank,MPI_COMM_WORLD,status,ierror)
        call MPI_Recv(filen,length,MPI_CHARACTER,0,rank,MPI_COMM_WORLD,status,ierror)
        call MPI_Recv(eq_incr,1,MPI_REAL,0,rank,MPI_COMM_WORLD,status,ierror)
        call MPI_Recv(ifa,1,MPI_INTEGER,0,rank,MPI_COMM_WORLD,status,ierror)
      endif

      ! check
      if( ierror .ne. 0 ) then
        print*,'error in syncInput()',rank,ierror
        call MPI_Abort(MPI_COMM_WORLD, ierror )
        call MPI_FINALIZE(ierror)                      
        stop 'abort - syncInput'
      endif
      end


!-----------------------------------------------------------------------
      subroutine syncHits(rank,nprocesses,hitcount,nfree)
!-----------------------------------------------------------------------     
! sync with master process and print to file
      implicit none
      include 'mpif.h'
      integer:: ierror,nprocesses,rank,i,k,length,status(MPI_STATUS_SIZE) 
      integer:: nfree,datacount,jrecord,restprocesses    
      integer:: hitcount(nfree),hitcount_proc(nfree)  

      ! nothing to do for single processor job
      if( nprocesses .lt. 2 ) return

      !print*,'rank',rank,nfree,nz

      ! get row data and print to file
      if( rank .eq. 0 ) then        
        ! MASTER receives        
        do i=1,nprocesses-1
          ! get data
          !print*,'getting from',i
          call MPI_Recv(hitcount_proc,nfree,MPI_INTEGER,i,i,MPI_COMM_WORLD,status,ierror)          
          if( ierror .ne. 0 ) then 
            print*,'error hitcount_proc',rank,ierror
            call stopProgram("syncHits    ")
          endif          
               
          ! add to master's hitcounts
          do k=1,nfree
            hitcount(k)=hitcount(k)+hitcount_proc(k)          
          enddo
        enddo
      else
        ! slave process sends information
          call MPI_Send(hitcount,nfree,MPI_INTEGER,0,rank,MPI_COMM_WORLD,ierror)        
          if( ierror .ne. 0 ) then 
            print*,'error hitcount process',rank,ierror
            call stopProgram("syncHits    ")
          endif                    
      endif
      end subroutine
      
!-----------------------------------------------------------------------
      subroutine syncFlag(rank,nprocesses,flag)
!-----------------------------------------------------------------------     
! syncs a (logical) flag between all processes, 
! assumes that the flag by default is .false.
      implicit none
      include 'mpif.h'
      integer,intent(in):: rank,nprocesses      
      logical,intent(inout):: flag
      integer:: ierror,status(MPI_STATUS_SIZE) 
      logical:: sum_flag
      
      ! nothing to do for single processor job
      if( nprocesses < 2 ) return

      ! sees if a process has the flag set to true ( logical or )
      sum_flag = .false.
      call MPI_AllReduce(flag,sum_flag,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierror)      
      if( ierror /= 0 ) call stopProgram("syncFlag error  " )
      
      ! returns flag which should be the same for all processes now
      flag = sum_flag            

      end subroutine
      

!-----------------------------------------------------------------------
      subroutine syncProcesses()
!-----------------------------------------------------------------------
! syncs a (logical) flag between all processes,
! assumes that the flag by default is .false.
      implicit none
      include 'mpif.h'
      integer:: ierror

      ! wait until all processes reached this point
      call MPI_Barrier( MPI_COMM_WORLD, ierror )

      end subroutine

