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
      ! local parameters
      integer:: n,ier

      ! check
      if (nproc < 2) return

      ! main process sends parameters to all other processes
      if (rankid == 0) then
        ! send to secondary processes
        do n = 1, nproc-1
          call MPI_Send(VERBOSE,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(subdivisions,1,MPI_INTEGER,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(FIRSTTIME,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(LASTTIME,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(cphaseRef,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(HETEROGENEOUS,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(DELTA,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(DELTARADIUS,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(deltaLat,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(deltaLon,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(DELTAfunction,8,MPI_CHARACTER,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(sourceLat,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(sourceLon,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(receiverLat,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(receiverLon,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(deltaMoveIncrement,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(manyReceivers,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(numofReceivers,1,MPI_INTEGER,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(manyKernels,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(kernelStartDistance,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(kernelEndDistance,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(importKernelsReceivers,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(SIMULATIONOUTPUT,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(cphasetype,8,MPI_CHARACTER,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(MOVEDELTA,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(SECONDDELTA,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(latitudeStart,1,MPI_INTEGER,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(latitudeEnd,1,MPI_INTEGER,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(longitudeEnd,1,MPI_INTEGER,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(deltaPerturbation,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(datadirectory,len(datadirectory),MPI_CHARACTER,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(PARALLELSEISMO,1,MPI_LOGICAL,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(phaseBlockFile,len(phaseBlockFile),MPI_CHARACTER,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(phaseBlockVelocityReference,1,MPI_CUSTOM,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(heterogeneousDataFile,len(heterogeneousDataFile),MPI_CHARACTER,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(heterogeneousPixelsize,1,MPI_REAL,n,n,MPI_COMM_WORLD,ier)
          call MPI_Send(gsh_maximum_expansion,1,MPI_INTEGER,n,n,MPI_COMM_WORLD,ier)
        enddo
      else
          ! get parameters from main process
          call MPI_Recv(VERBOSE,1,MPI_LOGICAL,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(subdivisions,1,MPI_INTEGER,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(FIRSTTIME,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(LASTTIME,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(cphaseRef,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(HETEROGENEOUS,1,MPI_LOGICAL,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(DELTA,1,MPI_LOGICAL,0,rank,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(DELTARADIUS,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(deltaLat,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(deltaLon,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(DELTAfunction,8,MPI_CHARACTER,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(sourceLat,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(sourceLon,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(receiverLat,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(receiverLon,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(deltaMoveIncrement,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(manyReceivers,1,MPI_LOGICAL,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(numofReceivers,1,MPI_INTEGER,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(manyKernels,1,MPI_LOGICAL,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(kernelStartDistance,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(kernelEndDistance,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(importKernelsReceivers,1,MPI_LOGICAL,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(SIMULATIONOUTPUT,1,MPI_LOGICAL,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(cphasetype,8,MPI_CHARACTER,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(MOVEDELTA,1,MPI_LOGICAL,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(SECONDDELTA,1,MPI_LOGICAL,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(latitudeStart,1,MPI_INTEGER,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(latitudeEnd,1,MPI_INTEGER,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(longitudeEnd,1,MPI_INTEGER,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(deltaPerturbation,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(datadirectory,len(datadirectory),MPI_CHARACTER,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(PARALLELSEISMO,1,MPI_LOGICAL,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(phaseBlockFile,len(phaseBlockFile),MPI_CHARACTER,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(phaseBlockVelocityReference,1,MPI_CUSTOM,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(heterogeneousDataFile,len(heterogeneousDataFile),MPI_CHARACTER,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(heterogeneousPixelsize,1,MPI_REAL,0,rankid,MPI_COMM_WORLD,status,ier)
          call MPI_Recv(gsh_maximum_expansion,1,MPI_INTEGER,0,rankid,MPI_COMM_WORLD,status,ier)
      endif
      end subroutine

!-----------------------------------------------------------------------
      subroutine syncInitialData()
!-----------------------------------------------------------------------
! reads in grid point on main machine and synchronizes with
! its secondary processes
      use cells; use parallel !; use propagationStartup; use adjointVariables; use verbosity
      use propagationStartup, only: SIMULATIONOUTPUT
      implicit none
      ! local parameters
      integer:: n,length,ier

      ! read vertices values from files
      if (MAIN_PROCESS) then
        ! send to secondary processes
        do n = 1, nprocesses-1
          ! vertices
          length = MaxVertices*3
          tag = n
          call MPI_Send(vertices,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ier)
          if (ier /= 0) then
              print *,'Error: syncInitial - vertices',n,' error',ier
              call stopProgram( 'abort - syncInitialData    ')
          endif
          ! numVertices
          length = 1
          call MPI_Send(numVertices,length,MPI_INTEGER,n,tag,MPI_COMM_WORLD,ier)

          ! cellNeighbors
          length = MaxVertices*7
          call MPI_Send(cellNeighbors,length,MPI_INTEGER,n,tag,MPI_COMM_WORLD,ier)
          if (ier /= 0) then
              print *,'Error: syncInitial - cellNeighbors',n,' error',ier
              call stopProgram( 'abort - syncInitialData    ')
          endif
          ! numNeighbors
          length = 1
          call MPI_Send(numNeighbors,length,MPI_INTEGER,n,tag,MPI_COMM_WORLD,ier)

          ! cellCorners
          length = MaxTriangles*3
          call MPI_Send(cellCorners,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ier)
          if (ier /= 0) then
              print *,'Error: syncInitial - cellCorners',n,' error',ier
              call stopProgram( 'abort - syncInitialData    ')
          endif

          ! numCorners
          length = 1
          call MPI_Send(numCorners,length,MPI_INTEGER,n,tag,MPI_COMM_WORLD,ier)

          ! cellFace
          length = MaxVertices*7
          call MPI_Send(cellFace,length,MPI_INTEGER,n,tag,MPI_COMM_WORLD,ier)
          if (ier /= 0) then
              print *,'Error: syncInitial - cellFaces',n,' error',ier
              call stopProgram( 'abort - syncInitialData    ')
          endif

          ! numFaces
          length = 1
          call MPI_Send(numFaces,length,MPI_INTEGER,n,tag,MPI_COMM_WORLD,ier)

          if (Station_Correction .or. SIMULATIONOUTPUT) then
            ! cellTriangleFace
            length = MaxTriangles*3
            call MPI_Send(cellTriangleFace,length,MPI_INTEGER,n,tag,MPI_COMM_WORLD,ier)
            if (ier /= 0) then
              print *,'Error: syncInitial - cellTriangleFaces',n,' error',ier
              call stopProgram( 'abort - syncInitialData    ')
            endif

            ! numFaces
            length = 1
            call MPI_Send(numTriangleFaces,length,MPI_INTEGER,n,tag,MPI_COMM_WORLD,ier)
          endif
       enddo
      else
        ! receive from main process
        ! vertices
        !print *,'secondary ',rank,' trying to receive vert. from main process'
        length = MaxVertices*3
        tag = rank
        call MPI_RECV(vertices,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ier)
        if (ier /= 0) then
              print *,'Error: syncInitial - vertices',rank,' error',ier
              call stopProgram( 'abort - syncInitialData    ')
        endif
        ! numVertices
        length = 1
        call MPI_RECV(numVertices,length,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ier)

        ! cellNeighbors
        length = MaxVertices*7
        call MPI_RECV(cellNeighbors,length,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ier)
        if (ier /= 0) then
              print *,'Error: syncInitial - cellNeighbors',n,' error',ier
              call stopProgram( 'abort - syncInitialData    ')
        endif
        ! numNeighbors
        length = 1
        call MPI_RECV(numNeighbors,length,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ier)

        ! cellCorners
        length = MaxTriangles*3
        call MPI_RECV(cellCorners,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ier)
        if (ier /= 0) then
              print *,'Error: syncInitial - cellCorners',n,' error',ier
              call stopProgram( 'abort - syncInitialData    ')
        endif
        ! numCorners
        length = 1
        call MPI_RECV(numCorners,length,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ier)

        ! cellFace
        length = MaxVertices*7
        call MPI_RECV(cellFace,length,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ier)
        if (ier /= 0) then
              print *,'Error: syncInitial - cellFace',n,' error',ier
              call stopProgram( 'abort - syncInitialData    ')
        endif
        ! numFaces
        length = 1
        call MPI_RECV(numFaces,length,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ier)

        if (Station_Correction .or. SIMULATIONOUTPUT) then
          ! cellTriangleFace
          length = MaxTriangles*3
          call MPI_RECV(cellTriangleFace,length,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ier)
          if (ier /= 0) then
            print *,'Error: syncInitial - cellTriangleFace',n,' error',ier
            call stopProgram( 'abort - syncInitialData    ')
          endif
          ! numFaces
          length = 1
          call MPI_RECV(numTriangleFaces,length,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ier)
        endif
       endif
      end subroutine


!!-----------------------------------------------------------------------
!      subroutine syncDisplacement( myrank, numberofprocesses )
!!-----------------------------------------------------------------------
!! synchronizes displacement and displacement_old arrays with
!! its secondary processes
!      use displacements;use propagationStartup
!      use verbosity; use parallel; use cells
!      implicit none
!      integer,intent(in):: myrank,numberofprocesses
!      ! local parameters
!      integer:: length,n,ier
!
!      if (rank == 0) then
!        ! send to secondary processes
!        do n = 1, nprocesses-1
!          ! displacement
!          length = numVertices
!          tag = n
!          call MPI_Send(displacement,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ier)
!          ! displacement_old
!          length = numVertices
!          tag = n
!          call MPI_Send(displacement_old,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ier)
!        enddo
!      else
!        ! receive from main process
!        ! displacement
!        length = numVertices
!        tag = rank
!        call MPI_RECV(displacement(1:numVertices),length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ier)
!        ! displacement_old
!        length = numVertices
!        tag = rank
!        call MPI_RECV(displacement_old(1:numVertices),length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ier)
!      endif
!
!      !debug
!      !  if (rank == 0) then
!      !    print *
!      !    print *,'main process'
!      !    do n=1340,1350
!      !      print *,'  #',displacement(n)
!      !    enddo
!      !  else
!      !    print *
!      !    print *,'secondary',rank
!      !    do n=1340,1350
!      !      print *,'  #',displacement(n)
!      !    enddo
!      !  endif
!      end subroutine

!-----------------------------------------------------------------------
      subroutine syncNewdisplacement()
!-----------------------------------------------------------------------
! synchronizes newdisplacement array
! each process gets only the needed informations from the other ones
      use displacements;use griddomain;use parallel;use verbosity
      implicit none
      ! local parameters
      integer:: dest,source,sendRange,recRange,n,k,neighbor,ier

      ! check if many processes
      if (nprocesses == 1) return

      !vertex range
      sendRange = boundariesMaxRange
      recRange = sendRange
      if (sendRange == 0) then
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
        if (neighbor /= -1 .and. neighbor /= rank) then
          !fill send info (over a fixed range of values)
          do k = 1,sendRange
            if (boundaries(rank,neighbor,k) > 0) then
              sendDisp(k) = newdisplacement(boundaries(rank,neighbor,k))
            endif
          enddo

          ! talk to neighbor process  (destination to send and source to receive from are the same process)
          dest = neighbor
          source = neighbor
          call MPI_SendRecv(sendDisp,sendRange,MPI_CUSTOM,dest,rank,receiveDisp,recRange,MPI_CUSTOM,source, &
                           source,MPI_COMM_WORLD, status, ier)
          if (ier /= 0) then
              print *,'syncNewdisplacement - p',rank,' error',ier
          endif
          !these are the vertices we received
          do k = 1,recRange
            ! check if entry is a valid vertex index in the neighbors boundary array
            if (boundaries(source,rank,k) > 0) then
              ! add as new displacement value for this vertex
              newdisplacement(boundaries(source,rank,k))=receiveDisp(k)
            endif
          enddo
        endif
      enddo

      !debug
      !  print *,'      sync done, barrier to wait.',rank
      !  ! wait until all processes reached this point
      !  call MPI_Barrier( MPI_COMM_WORLD, ier )
      !  if (ier /= 0) call stopProgram('abort - MPI_Barrier kernels failed    ')

      end subroutine

!-----------------------------------------------------------------------
      subroutine syncNewdisplacement_alt( rank, nprocesses )
!-----------------------------------------------------------------------
! synchronizes newdisplacement array with
! its secondary processes
      use displacements; use propagationStartup; use cells
      use mpi
      implicit none
      integer,intent(in):: rank,nprocesses
      ! local parameters
      integer:: ier,n
      integer:: status(MPI_STATUS_SIZE)
      integer:: startLocationVertex, endLocationVertex,recStartVertex,recEndVertex
      double precision, allocatable, dimension(:) :: tmpSend,tmpReceiv
      integer:: dest,source,range,vertexpart,sendRange,recRange,rangeMax

      !these are the vertices which are calculated by this process
      !startLocationVertex = rank*range+1
      !endLocationVertex = (rank+1)*range
      ! determine the first and last vertice to send
      vertexpart = int(numVertices/nprocesses)
      startLocationVertex = rank*vertexpart + 1
      endLocationVertex = (rank+1)*vertexpart
      ! add remaining vertices to the last process
      if (rank+1 == nprocesses) then
        endLocationVertex = endLocationVertex + mod(numVertices,nprocesses)
      endif

      !check vertex range
      sendRange = endLocationVertex - startLocationVertex + 1
      if (sendRange < 1) call stopProgram( 'abort-syncNewdisplacement short range    ')

      !allocate array for the data to send
      allocate( tmpSend(sendRange), stat=ier)
      if (ier /= 0) then
        print *,'Error: rank ',rank,' error in allocating send array'
        call stopProgram( 'abort - syncNewdisplacement_alt    ')
      endif
      tmpSend(1:sendRange) = newdisplacement(startLocationVertex:endLocationVertex)

      !allocate array for the data to receive
      rangeMax = 0
      do n = 0, nprocesses -1
        ! determine the first and last vertice to receive
        recStartVertex = n*int(numVertices/nprocesses) + 1
        recEndVertex = (n+1)*int(numVertices/nprocesses)
        if (n+1 == nprocesses) then  ! add remaining vertices to the last process
          recEndVertex = recEndVertex + mod(numVertices,nprocesses)
        endif
        range = recEndVertex-recStartVertex+1
        if (range > rangeMax) rangeMax = range
      enddo
      ! prepare array to receive data
      allocate( tmpReceiv(rangeMax), stat=ier)
      if (ier /= 0) then
        print *,'Error: rank ',rank,' error in allocating receiving array ',rangeMax
        call stopProgram( 'abort - syncNewdisplacement_alt    ')
      endif

      !print *,rank,'sendreceiving...'
      do n = 0, nprocesses - 1
        ! only talk to other processes
        if (n /= rank) then
          ! determine the first and last vertice to receive
          recStartVertex = n*int(numVertices/nprocesses) + 1
          recEndVertex = (n+1)*int(numVertices/nprocesses)
          if (n+1 == nprocesses) then  ! add remaining vertices to the last process
            recEndVertex = recEndVertex + mod(numVertices,nprocesses)
          endif
          recRange = recEndVertex-recStartVertex+1

          ! talk to neighbor process
          dest = n
          source = n
          call MPI_SendRecv( tmpSend,sendRange, MPI_DOUBLE_PRECISION, &
                             dest,rank,tmpReceiv,recRange,MPI_DOUBLE_PRECISION,source, &
                             source,MPI_COMM_WORLD, status, ier)
          if (ier /= 0) then
              print *,'syncNewdisplacement - p',rank,' error',ier
          endif
          !these are the vertices we received
          !startLocationVertex = source*range+1
          !endLocationVertex = (source+1)*range
          newdisplacement(recStartVertex:recEndVertex) = tmpReceiv(1:recRange)

        endif
      enddo
      end subroutine

!!-----------------------------------------------------------------------
!      subroutine syncNewdisplacement_alt2(rank, nprocesses)
!!-----------------------------------------------------------------------
!! synchronizes newdisplacement array with
!! its secondary processes
!      use displacements; use propagationStartup; use cells
!      use mpi
!      implicit none
!      integer,intent(in):: rank,nprocesses
!      ! local parameters
!      integer:: ier,n
!      integer:: status(MPI_STATUS_SIZE)
!      integer:: startLocationVertex, endLocationVertex
!      !double precision:: tmp(numVertices/nprocesses)
!      double precision:: tmpSend(numVertices/nprocesses)
!      double precision:: tmpReceiv(numVertices/nprocesses)
!      integer:: sendcount, sendtag, sendtype,dest, source,recvtag,comm,recvtype, recvcount,range
!
!      range = numVertices/nprocesses
!      if (range < 1) call stopProgram( 'abort-syncNewdisplacement short range    ')
!
!      if (rank == 0) then
!          ! receive from secondary processes
!          !print *,'sync newdisplacement' !debug
!          do n = 1, nprocesses-1
!            !print *,'receiving tmpdisp. from secondary',n !debug
!            ! new displacements
!
!
!!            length=numVertices/nprocesses
!!            tag = n
!!            call MPI_Recv(tmp,length, &
!!     &                          MPI_DOUBLE_PRECISION,n,tag, &
!!     &                          MPI_COMM_WORLD,status,ier)
!            !print *,'main process ', range
!
!            tmpSend(1:range) = newdisplacement(1:range)
!            tmpReceiv = 0.0
!
!            sendcount = range
!            sendtype = MPI_DOUBLE_PRECISION
!            dest = n
!            sendtag = 0
!            recvcount = range
!            recvtype = MPI_DOUBLE_PRECISION
!            source = n
!            recvtag = n
!            comm = MPI_COMM_WORLD
!
!            !print *,'main process sendreceiving...'
!
!            call MPI_SendRecv( tmpSend,sendcount, sendtype, dest, &
!                               sendtag, tmpReceiv, recvcount, recvtype, source, recvtag, &
!                               comm, status, ier)
!            if (ier /= 0) then
!              print *,'Error: syncNewdisplacement - main process error',ier
!              call stopProgram( 'abort - syncNewdisplacement_alt2   ')
!            endif
!
!            ! fill in corresponding values
!            startLocationVertex = n*range + 1
!            endLocationVertex = (n+1)*range
!            !print *,'main process received', startLocationVertex, endLocationVertex
!
!            newdisplacement(startLocationVertex:endLocationVertex)=tmpReceiv(1:range)
!
!            !do i=startLocationVertex, endLocationVertex
!            !  newdisplacement(i) = tmpReceiv(i)
!            !enddo
!!            if ( (numVertices/nprocesses) /= &
!!     &          (endLocationVertex-startLocationVertex+1)) then
!!              stop 'syncNewdisplacement - wrong array size'
!!            endif
!
!!            newdisplacement(startLocationVertex:endLocationVertex) =                    &
!!     &                   tmp(1:numVertices/nprocesses)
!          enddo
!
!          ! send updated whole array to secondary processes
!
!!          do n=1, nprocesses-1
!!            !print *,'sending newdisp. to secondary',n !debug
!!            length=numVertices
!!            tag = n
!!            call MPI_Send(newdisplacement,length,MPI_DOUBLE_PRECISION, &
!!     &                  n,tag,MPI_COMM_WORLD,ier)
!!          enddo
!      else
!          ! send to main process
!          !print *,'secondary',rank,'sending newdisp. to main process' !debug
!          ! new displacements
!!          startLocationVertex = rank*numVertices/nprocesses + 1
!!          endLocationVertex = (rank+1)*numVertices/nprocesses
!!          if ( (numVertices/nprocesses) /= &
!!     &      (endLocationVertex-startLocationVertex+1)) then
!!              stop 'syncNewdisplacement - wrong array size'
!!          endif
!
!!          tmp(1:numVertices/nprocesses)=                                &
!!     &              newdisplacement(startLocationVertex:endLocationVertex)
!!          length=numVertices/nprocesses
!!          tag = rank
!!          call MPI_Send(tmp,length,MPI_DOUBLE_PRECISION, &
!!     &                0,tag,MPI_COMM_WORLD,ier)
!
!          ! receive from main process new updated array
!          !print *,'secondary',rank,'receiving newdisp. from main process'
!!          length=numVertices
!!          tag = rank
!!          call MPI_Recv(newdisplacement,length,MPI_DOUBLE_PRECISION, &
!!     &                0,tag,MPI_COMM_WORLD,status,ier)
!
!            startLocationVertex = rank*range + 1
!            endLocationVertex = (rank+1)*range
!            !print *,'secondary ', startLocationVertex, endLocationVertex
!
!            tmpSend(1:range) = newdisplacement(startLocationVertex:endLocationVertex)
!            tmpReceiv = 0.0
!
!            sendcount = range
!            sendtype = MPI_DOUBLE_PRECISION
!            dest = 0
!            sendtag = rank
!            recvcount = range
!            recvtype = MPI_DOUBLE_PRECISION
!            source = 0
!            recvtag = 0
!            comm = MPI_COMM_WORLD
!
!            !print *,'secondary', rank,' sendreceiving...'
!
!            call MPI_SendRecv( tmpSend,sendcount, sendtype, dest, &
!                               sendtag, tmpReceiv, recvcount, recvtype, source, recvtag, &
!                               comm, status, ier)
!            if (ier /= 0) then
!              print *,'Error: syncNewdisplacement - secondary', rank,'error',ier
!              call stopProgram( 'abort - syncNewdisplacement_alt2   ' )
!            endif
!
!            ! fill in corresponding values
!            !print *,'secondary received', range
!            newdisplacement(1:range)=tmpReceiv(1:range)
!
!      endif
!
!      !debug
!      !if (rank == 0) then
!      !  print *
!      !  print *,'new main process'
!      !  do l=1340,1350
!      !    print *,'  m#',newdisplacement(l)
!      !  enddo
!      !open(11,file='syncNew0.dat')
!      !  write(11,*) newdisplacement
!      !  close(11)
!      !else
!      !  print *
!      !  print *,'new secondary',rank
!      !  do l=1340,1350
!      !    print *,'  #',newdisplacement(l)
!      !  enddo
!      !  open(12,file='syncNew1.dat')
!      !  write(12,*) newdisplacement
!      ! close(12)
!      !endif
!
!      end subroutine

!-----------------------------------------------------------------------
      subroutine syncPrecalculated()
!-----------------------------------------------------------------------
! synchronize precalculated array cellAreas, cellEdgesLength, cellCenterDistances
      use cells; use parallel; use propagationStartup; use verbosity
      implicit none
      ! local parameters
      integer:: length,n,first,lastsent,i,ier
      real(WP), allocatable, dimension(:):: tmpExchange
      integer,parameter::SLICE = 1000000 ! assuming MPI has a limited buffer size

      ! debug check arrays
      !  ! arrays
      !  if (.not. allocated(cellAreas)) then
      !    print *,'Error: cellAreas not allocated',rank
      !    call stopProgram("array not allocated    ")
      !  endif
      !  if (.not. allocated(cellEdgesLength)) then
      !    print *,'Error: cellEdgesLength not allocated',rank
      !    call stopProgram("array not allocated    ")
      !  endif
      !  if (.not. allocated(cellCenterDistances)) then
      !    print *,'Error: cellCenterDistances not allocated',rank
      !    call stopProgram("array not allocated    ")
      !  endif
      !  ! barrier
      !  call MPI_Barrier(MPI_COMM_WORLD,ier)
      !  if (ier /= 0) call stopProgram('syncRoutines - MPI_Barrier kernels failed    ')

      if (MAIN_PROCESS) then
        !read in values
        call readPrecalculated()

        ! send to secondary processes
        do n = 1, nprocesses-1
          tag = n

          ! number of vertices
          length = numVertices

          ! split into chunks to avoid shared memory allocation problems of MPI
          if (length > SLICE) then
            lastsent = 0
            do while(lastsent < numVertices )
              first = lastsent+1
              lastsent = lastsent+SLICE
              if (lastsent >= numVertices) lastsent=numVertices
              length = lastsent-first+1
              call MPI_Send(cellAreas(first:lastsent),length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ier)
            enddo
          else
            call MPI_Send(cellAreas,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ier)
          endif

          length = numVertices
          ! split into chunks to avoid shared memory allocation problems of MPI
          if (length > SLICE) then
            lastsent = 0
            do while(lastsent < numVertices )
              first = lastsent+1
              lastsent = lastsent+SLICE
              if (lastsent >= numVertices) lastsent=numVertices
              length=(lastsent-first+1)
              allocate(tmpExchange(length),stat=ier)
              if (ier /= 0) call stopProgram('could not allocate temporary exchange array!      ')

              do i = 0,6
                tmpExchange(:) = cellEdgesLength(first:lastsent,i)
                call MPI_Send(tmpExchange,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ier)

                tmpExchange(:) = cellCenterDistances(first:lastsent,i)
                call MPI_Send(tmpExchange,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ier)

                if (CORRECT_RATIO) then
                  tmpExchange(:) = cellFractions(first:lastsent,i)
                  call MPI_Send(tmpExchange,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ier)
                endif
              enddo
              deallocate(tmpExchange)
            enddo
          else
            length = numVertices*7
            call MPI_Send(cellEdgesLength,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ier)
            call MPI_Send(cellCenterDistances,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ier)
            if (CORRECT_RATIO) then
              call MPI_Send(cellFractions,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ier)
            endif
          endif
       enddo
      else
        ! receive from main process
        tag = rank
        length = numVertices
        ! split into chunks to avoid shared memory allocation problems of MPI
        if (length > SLICE) then
          lastsent = 0
          do while(lastsent < numVertices )
            first = lastsent+1
            lastsent = lastsent+SLICE
            if (lastsent >= numVertices) lastsent=numVertices
            length = lastsent-first+1
            call MPI_RECV(cellAreas(first:lastsent),length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ier)
          enddo
        else
          call MPI_RECV(cellAreas,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ier)
        endif

        length = numVertices
        ! split into chunks to avoid shared memory allocation problems of MPI
        if (length > SLICE) then
          lastsent = 0
          do while(lastsent < numVertices )
            first = lastsent+1
            lastsent = lastsent+SLICE
            if (lastsent >= numVertices) lastsent=numVertices
            length=(lastsent-first+1)
            allocate(tmpExchange(length),stat=ier)
            if (ier /= 0) call stopProgram('could not allocate temporary exchange array in secondary!      ')

            do i = 0,6
              call MPI_RECV(tmpExchange,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ier)
              cellEdgesLength(first:lastsent,i) = tmpExchange(:)

              call MPI_RECV(tmpExchange,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ier)
              cellCenterDistances(first:lastsent,i) = tmpExchange(:)

              if (CORRECT_RATIO) then
                call MPI_RECV(tmpExchange,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ier)
                cellFractions(first:lastsent,i) = tmpExchange(:)
              endif
            enddo
            deallocate(tmpExchange)
          enddo
        else
          length = numVertices*7
          call MPI_RECV(cellEdgesLength,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ier)
          call MPI_RECV(cellCenterDistances,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ier)
          if (CORRECT_RATIO) then
            call MPI_RECV(cellFractions,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ier)
          endif
        endif
      endif
      end subroutine

!-----------------------------------------------------------------------
      subroutine syncPhaseMap()
!-----------------------------------------------------------------------
! synchronize phase map
      use cells;use phaseVelocityMap;use propagationStartup;use verbosity;use parallel
      implicit none
      ! local parameters
      integer:: length,n,ier

      if (rank == 0) then
        ! send to secondary processes
        do n = 1, nprocesses-1
          ! vertices
          length = numVertices
          tag = n

          call MPI_Send(phaseMap,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ier)
       enddo
      else
        ! receive from main process
        length = numVertices
        tag = rank
        call MPI_RECV(phaseMap,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ier)
      endif
      end subroutine


!-----------------------------------------------------------------------
      subroutine syncReceiversRef()
!-----------------------------------------------------------------------
! synchronize array receiversSeismogramRef
      use cells; use parallel; use propagationStartup; use verbosity
      implicit none
      ! local parameters
      integer:: length,n,m,ier

      if (MAIN_PROCESS) then
        ! send to secondary processes
        do n = 1, nprocesses-1
          ! vertices
          if (manyKernels) then
            do m = 1,numofKernels
              length = size(kernelsReceiversSeismogramRef(m,:,:))
              tag = n
              call MPI_Send(kernelsReceiversSeismogramRef(m,:,:),length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ier)
            enddo
          else
            length=size(receiversSeismogramRef)
            tag = n
            call MPI_Send(receiversSeismogramRef,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ier)
          endif
       enddo
      else
        if (manyKernels) then
          do m = 1,numofKernels
            ! receive from main process
            length = size(kernelsReceiversSeismogramRef(m,:,:))
            tag = rank
            call MPI_Recv(kernelsReceiversSeismogramRef(m,:,:),length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ier)
          enddo
        else
          ! receive from main process
          length=size(receiversSeismogramRef)
          tag = rank
          call MPI_Recv(receiversSeismogramRef,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ier)
        endif
      endif
      end subroutine

!-----------------------------------------------------------------------
      subroutine collectReceiversSeismogram(doReference)
!-----------------------------------------------------------------------
! collects array receiversSeismogram from receivers not in the main process domain
      use cells; use parallel; use propagationStartup;use griddomain; use verbosity
      implicit none
      logical,intent(in):: doReference
      ! local parameters
      integer:: i,countVertices,recRange,sendRange,n,m,ier
      real(WP),allocatable,dimension(:,:)::tmpExchange
!      real(WP),allocatable,dimension(:,:,:)::tmpManyExchange


      ! check if many processes
      if (nprocesses == 1) return

      if (MAIN_PROCESS) then
        ! debug
        !print *,'collectReceiversSeismogram() -main process',nprocesses,rank,MAIN_PROCESS,numofReceivers,doReference

        ! get from the other processes
        do n = 1, nprocesses-1
          !vertex range
          countVertices = 0
          if (manyKernels) then
          !  countVertices = numofReceivers
            continue
          else
            do i = 1,numofReceivers
              if (vertexDomain(receivers(i)) == n) then
                countVertices = countVertices+1
              endif
            enddo
          endif
          recRange = countVertices

          ! allocate array to receive
          if (manyKernels) then
            !allocate(tmpManyExchange(numofKernels,recRange,numofTimeSteps),stat=ier)
            !allocate(tmpExchange(recRange,numofTimeSteps),stat=ier)
            continue
          else
            allocate(tmpExchange(recRange,numofTimeSteps),stat=ier)
            if (ier /= 0) then
              print *,'Error: main process collects receivers seismogram:',numofKernels,recRange,numofTimeSteps
              call stopProgram('collectReceiversSeismogram() - error allocating receive array   ')
            endif
          endif

          ! get info from secondary
          tag = n
          if (manyKernels) then
            !call MPI_Recv(tmpManyExchange,size(tmpManyExchange),MPI_CUSTOM,n,tag,MPI_COMM_WORLD,status,ier)
            do m = 1,numofKernels
              countVertices = 0
              do i = 1,numofReceivers
                if (vertexDomain(KernelsReceivers(m,i)) == n) then
                  countVertices = countVertices+1
                endif
              enddo
              recRange = countVertices

              allocate(tmpExchange(recRange,numofTimeSteps),stat=ier)
              if (ier /= 0) then
                print *,'Error: main process collects receivers seismogram:',m,numofKernels,recRange,numofTimeSteps
                call stopProgram('collectReceiversSeismogram() - error allocating receive array   ')
              endif

              call MPI_Recv(tmpExchange,size(tmpExchange),MPI_CUSTOM,n,tag,MPI_COMM_WORLD,status,ier)
              if (ier /= 0) then
                print *,'Error: collectReceiversSeismogram - ',rank,n,ier
                call stopProgram('collectReceiversSeismogram() - error mpi-recv   ')
              endif

              ! these are the vertices we received
              countVertices = 0
              do i = 1,numofReceivers
                ! extract only seismograms from receiver stations in that particular process domain
                if (vertexDomain(kernelsReceivers(m,i)) == n) then
                  countVertices = countVertices + 1
                  if (doReference) then
                    kernelsReceiversSeismogramRef(m,i,:)=tmpExchange(countVertices,:)
                  else
                    kernelsReceiversSeismogram(m,i,:)=tmpExchange(countVertices,:)
                  endif
                endif
              enddo

              deallocate( tmpExchange )
            enddo
          else
            call MPI_Recv(tmpExchange,size(tmpExchange),MPI_CUSTOM,n,tag,MPI_COMM_WORLD,status,ier)
            if (ier /= 0) then
              print *,'Error: collectReceiversSeismogram - ',rank,n,ier
              call stopProgram('collectReceiversSeismogram() - error mpi-recv   ')
            endif
          endif

          ! these are the vertices we received
          countVertices = 0
          do i = 1,numofReceivers
            ! check if entry is a valid vertex index in the neighbors boundary array
            if (manyKernels) then
!              do m=1,numofKernels
!                ! extract only seismograms from receiver stations in that particular process domain
!                if (vertexDomain(kernelsReceivers(m,i)) == n) then
!                  if (doReference) then
!                    kernelsReceiversSeismogramRef(m,i,:)=tmpManyExchange(m,i,:)
!                  else
!                    kernelsReceiversSeismogram(m,i,:)=tmpManyExchange(m,i,:)
!                  endif
!                endif
!              enddo
              continue
            else
              if (vertexDomain(receivers(i)) == n) then
                countVertices = countVertices+1
                if (doReference) then
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
        !print *,'collectReceiversSeismogram() - secondary',nprocesses,rank,MAIN_PROCESS,numofReceivers,doReference

        !vertex range
        countVertices = 0
        if (manyKernels) then
        !  countVertices=numofReceivers
          continue
        else
          do i = 1,numofReceivers
            if (vertexDomain(receivers(i)) == rank) then
              countVertices = countVertices+1
            endif
          enddo
        endif
        sendRange = countVertices

        ! allocate array
        if (manyKernels) then
          !allocate(tmpManyExchange(numofKernels,sendRange,numofTimeSteps),stat=ier)
          !allocate(tmpExchange(sendRange,numofTimeSteps),stat=ier)
          continue
        else
          allocate(tmpExchange(sendRange,numofTimeSteps),stat=ier)
          if (ier /= 0) then
            print *,'Error: secondary sends receivers seismograms:',rank,numofKernels,sendRange,numofTimeSteps
            call stopProgram('collectReceiversSeismogram() - error allocating send array    ')
          endif
        endif

        ! fill in corresponding values
        countVertices = 0
        do i = 1,numofReceivers
          if (manyKernels) then
!            ! fill in all seismograms
!            do m=1,numofKernels
!              if (doReference) then
!                tmpManyExchange(m,i,:)=kernelsReceiversSeismogramRef(m,i,:)
!              else
!                tmpManyExchange(m,i,:)=kernelsReceiversSeismogram(m,i,:)
!              endif
!            enddo
            continue
          else
            ! fill in only seismogram from receivers in this process domain
            if (vertexDomain(receivers(i)) == rank) then
              countVertices = countVertices+1
              if (doReference) then
                tmpExchange(countVertices,:)=receiversSeismogramRef(i,:)
              else
                tmpExchange(countVertices,:)=receiversSeismogram(i,:)
              endif
            endif
          endif
        enddo

        ! send to main process
        tag = rank
        if (manyKernels) then
          !call MPI_Send(tmpManyExchange,size(tmpManyExchange),MPI_CUSTOM,0,tag,MPI_COMM_WORLD,ier)
          do m = 1,numofKernels
            countVertices = 0
            do i = 1,numofReceivers
              if (vertexDomain(kernelsReceivers(m,i)) == rank) then
                countVertices = countVertices+1
              endif
            enddo
            sendRange = countVertices

            allocate(tmpExchange(sendRange,numofTimeSteps),stat=ier)
            if (ier /= 0) then
              print *,'Error: secondary sends receivers seismograms:',rank,m,numofKernels,sendRange,numofTimeSteps
              call stopProgram('collectReceiversSeismogram() - error allocating send array    ')
            endif

            ! fill in only seismogram from receivers in this process domain
            countVertices = 0
            do i = 1,numofReceivers
              if (vertexDomain(kernelsReceivers(m,i)) == rank) then
                countVertices = countVertices + 1
                if (doReference) then
                  tmpExchange(countVertices,:)=kernelsReceiversSeismogramRef(m,i,:)
                else
                  tmpExchange(countVertices,:)=kernelsReceiversSeismogram(m,i,:)
                endif
              endif
            enddo
            ! send
            call MPI_Send(tmpExchange,size(tmpExchange),MPI_CUSTOM,0,tag,MPI_COMM_WORLD,ier)
            if (ier /= 0) then
              print *,'collectReceiversSeismogram - ',rank,n,ier
            endif

            deallocate( tmpExchange )
          enddo
        else
          call MPI_Send(tmpExchange,size(tmpExchange),MPI_CUSTOM,0,tag,MPI_COMM_WORLD,ier)
          if (ier /= 0) then
            print *,'collectReceiversSeismogram - ',rank,n,ier
          endif
        endif ! manykernels
      endif ! MAIN_PROCESS
      end subroutine

!-----------------------------------------------------------------------
      subroutine syncBackwardNewdisplacement()
!-----------------------------------------------------------------------
! synchronizes backwardnewdisplacement array
! each process gets the needed (boundary) informations from the other ones
      use displacements;use griddomain;use parallel;use propagationStartup; use adjointVariables
      implicit none
      ! local parameters
      integer:: dest,source,sendRange,recRange,n,k,neighbor,ier

      ! check if many processes
      if (nprocesses == 1) return

      !vertex range
      sendRange = boundariesMaxRange
      recRange = sendRange
      if (sendRange == 0) then
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
        if (neighbor /= -1 .and. neighbor /= rank) then
          !fill send info (over a fixed range of values)
          do k = 1,sendRange
            if (boundaries(rank,neighbor,k) > 0) then
              sendDisp(k) = backwardNewdisplacement(boundaries(rank,neighbor,k))
            endif
          enddo

          ! talk to neighbor process  (destination to send and source to receive from are the same process)
          dest = neighbor
          source = neighbor
          call MPI_SendRecv( sendDisp,sendRange, MPI_CUSTOM, &
                             dest,rank,receiveDisp,recRange,MPI_CUSTOM,source, &
                             source,MPI_COMM_WORLD, status, ier)
          if (ier /= 0) then
              print *,'syncbackwardNewdisplacement - p',rank,' error',ier
          endif
          !these are the vertices we received
          do k = 1,recRange
            ! check if entry is a valid vertex index in the neighbors boundary array
            if (boundaries(source,rank,k) > 0) then
              ! add as new displacement value for this vertex
              backwardNewdisplacement(boundaries(source,rank,k))=receiveDisp(k)
            endif
          enddo
        endif
      enddo
      end subroutine


!-----------------------------------------------------------------------
      subroutine syncReceiverSeismogram()
!-----------------------------------------------------------------------
! let seismogram at receiver be available for all processes
      use parallel;use propagationStartup; use verbosity
      implicit none
      ! local parameters
      integer:: domain,i,length,ier
      integer,external::getDomain

      ! check if different grid domains exist
      if (nprocesses < 2) return

      ! get right seismogram at receiver station
      domain = getDomain(receiverVertex)
      if (domain == rank) then
        do i = 0,nprocesses-1
          if (i /= domain) then
            ! send seismogram
            length = size(seismogramReceiver(:,:))
            tag = rank
            call MPI_Send(seismogramReceiver,length,MPI_CUSTOM,i,tag,MPI_COMM_WORLD,ier)
          endif
        enddo
      else
        ! receive from domain process
        length = size(seismogramReceiver(:,:))
        tag = domain
        call MPI_Recv(seismogramReceiver,length,MPI_CUSTOM,domain,tag,MPI_COMM_WORLD,status,ier)
      endif

      end subroutine


!-----------------------------------------------------------------------
      subroutine collectAdjointKernel()
!-----------------------------------------------------------------------
! completes adjointKernel array for main process
! main process collects from all secondary processes the needed informations
      use griddomain;use parallel;use propagationStartup
      use adjointVariables; use cells
      implicit none
      ! local parameters
      integer:: sendRange,recRange,n,i,countVertices,ier
      real(WP),allocatable,dimension(:)::tmpExchange

      ! check if many processes
      if (nprocesses < 2) return

      ! only receive from secondary processes
      if (MAIN_PROCESS) then
        do n = 1, nprocesses-1
          !vertex range
          countVertices = 0
          do i = 1,numVertices
            if (vertexDomain(i) == n) then
              countVertices = countVertices+1
            endif
          enddo
          recRange = countVertices

          ! allocate array to receive
          allocate(tmpExchange(recRange),stat=ier)
          if (ier /= 0) then
            print *,'Error: main process allocate error:',ier,recRange,n
            call stopProgram('collectAdjointKernel() - error allocating receive array   ')
          endif

          ! get info from secondary
          tag = n
          call MPI_Recv(tmpExchange,recRange,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,status,ier)
          if (ier /= 0) then
              print *,'Error: collectAdjointKernel - ',rank,n,ier
              call stopProgram('collectAdjointKernel() - error mpi-recv   ')
          endif

          !these are the vertices we received
          countVertices = 0
          do i = 1,numVertices
            ! check if entry is a valid vertex index in the neighbors boundary array
            if (vertexDomain(i) == n) then
              countVertices = countVertices+1
              adjointKernel(i) = tmpExchange(countVertices)
            endif
          enddo

          ! free again
          deallocate(tmpExchange)

        enddo
      else
        !vertex range
        sendRange = size(domainVertices(:))

        ! allocate array
        allocate(tmpExchange(sendRange),stat=ier)
        if (ier /= 0) call stopProgram('collectAdjointKernel() - error allocating send array    ')

        ! fill in corresponding values
        do i = 1,sendRange
          tmpExchange(i) = adjointKernel(domainVertices(i))
        enddo

        ! send to main process
        tag = rank
        call MPI_Send(tmpExchange,sendRange,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,ier)
        if (ier /= 0) then
            print *,'error collectAdjointKernel - ',rank,n,ier
        endif
      endif
      end subroutine

!-----------------------------------------------------------------------
      subroutine collectFullNewdisplacement()
!-----------------------------------------------------------------------
! completes newdisplacement array for main process
! main process collects from all secondary processes the needed informations
      use displacements; use griddomain; use parallel
      use propagationStartup; use cells
      implicit none
      ! local parameters
      integer:: sendRange,recRange,n,i,countVertices,ier

      ! check if many processes
      if (nprocesses == 1) return

      !init
      sendDisp(:)=0.0_WP

      ! only receive from secondary processes
      if (MAIN_PROCESS) then
        do n = 1, nprocesses-1
          !vertex range
          countVertices = 0
          do i = 1,numVertices
            if (vertexDomain(i) == n) then
              countVertices = countVertices+1
            endif
          enddo
          recRange = countVertices

          ! get info from secondary
          receiveDisp(:)=0.0_WP
          tag = n
          call MPI_Recv(receiveDisp,recRange,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,status,ier)
          if (ier /= 0) then
              print *,'collectFullNewdisplacement - ',rank,n,ier
          endif
          !these are the vertices we received
          countVertices = 0
          do i = 1,numVertices
            ! check if entry is a valid vertex index in the neighbors boundary array
            if (vertexDomain(i) == n) then
              countVertices = countVertices+1
              newdisplacement(i) = receiveDisp(countVertices)
            endif
          enddo
        enddo
      else
        !vertex range
        sendRange = size(domainVertices(:))
        do i = 1,sendRange
          sendDisp(i) = newdisplacement(domainVertices(i))
        enddo

        ! send to main process
        tag = rank
        call MPI_Send(sendDisp,sendRange,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,ier)
        if (ier /= 0) then
            print *,'collectFullNewdisplacement - ',rank,n,ier
        endif
      endif
      end subroutine

!-----------------------------------------------------------------------
      subroutine syncInterpolationData()
!-----------------------------------------------------------------------
! synchronizes data related to interpolation with its secondary processes
      use cells; use parallel; use propagationStartup; use verbosity
      implicit none
      ! local parameters
      integer:: n,length,ier

      ! read vertices values from files
      if (MAIN_PROCESS) then
        ! send to secondary processes
        do n = 1, nprocesses-1
          ! distances
          length = 3
          tag = n
          call MPI_Send(interpolation_distances,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ier)
          if (ier /= 0) then
              print *,'Error: syncInterpolation - distances',n,' error',ier
              call stopProgram( 'abort - syncInterpolation    ')
          endif
          ! side lengths
          length = 3
          tag = n
          call MPI_Send(interpolation_triangleLengths,length,MPI_CUSTOM,n,tag,MPI_COMM_WORLD,ier)
          if (ier /= 0) then
              print *,'Error: syncInterpolation - lengths',n,' error',ier
              call stopProgram( 'abort - syncInterpolation    ')
          endif
          ! corner indices
          length = 3
          tag = n
          call MPI_Send(interpolation_corners,length,MPI_INTEGER,n,tag,MPI_COMM_WORLD,ier)
          if (ier /= 0) then
              print *,'Error: syncInterpolation - corners',n,' error',ier
              call stopProgram( 'abort - syncInterpolation    ')
          endif
       enddo
      else
        ! receive from main process
        ! distances
        length = 3
        tag = rank
        call MPI_RECV(interpolation_distances,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ier)
        if (ier /= 0) then
              print *,'Error: syncInterpolation - distances',rank,' error',ier
              call stopProgram( 'abort - syncInterpolation    ')
        endif
        ! side lengths
        length = 3
        tag = rank
        call MPI_RECV(interpolation_triangleLengths,length,MPI_CUSTOM,0,tag,MPI_COMM_WORLD,status,ier)
        if (ier /= 0) then
              print *,'Error: syncInterpolation - lengths',rank,' error',ier
              call stopProgram( 'abort - syncInterpolation    ')
        endif
        ! corner indices
        length = 3
        tag = rank
        call MPI_RECV(interpolation_corners,length,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ier)
        if (ier /= 0) then
              print *,'Error: syncInterpolation - corners',rank,' error',ier
              call stopProgram( 'abort - syncInterpolation    ')
        endif
      endif
      end subroutine

!-----------------------------------------------------------------------
      subroutine syncPixelgridInput(rank,nprocesses,txtfile,filen,eq_incr,ifa)
!-----------------------------------------------------------------------
! sync with other process
      use mpi
      implicit none
      integer,intent(in):: rank,nprocesses
      character(len=128),intent(inout):: txtfile,filen
      real,intent(inout):: eq_incr
      integer,intent(inout):: ifa
      ! local parameters
      integer:: ier,i,length
      integer:: status(MPI_STATUS_SIZE)

      ! string length
      length = 128

      ! nothing to do for single processor job
      if (nprocesses < 2) return

      ! synchronize input data
      if (rank == 0) then
        ! main process sends
        do i = 1, nprocesses-1
          call MPI_Send(txtfile,length,MPI_CHARACTER,i,i,MPI_COMM_WORLD,ier)
          call MPI_Send(filen,length,MPI_CHARACTER,i,i,MPI_COMM_WORLD,ier)
          call MPI_Send(eq_incr,1,MPI_REAL,i,i,MPI_COMM_WORLD,ier)
          call MPI_Send(ifa,1,MPI_INTEGER,i,i,MPI_COMM_WORLD,ier)
        enddo
      else
        ! secondary process receives information
        call MPI_Recv(txtfile,length,MPI_CHARACTER,0,rank,MPI_COMM_WORLD,status,ier)
        call MPI_Recv(filen,length,MPI_CHARACTER,0,rank,MPI_COMM_WORLD,status,ier)
        call MPI_Recv(eq_incr,1,MPI_REAL,0,rank,MPI_COMM_WORLD,status,ier)
        call MPI_Recv(ifa,1,MPI_INTEGER,0,rank,MPI_COMM_WORLD,status,ier)
      endif

      ! check
      if (ier /= 0) then
        print *,'error in syncInput()',rank,ier
        ! note: MPI_ABORT does not return, it makes the program exit with an error code of 30
        call MPI_Abort(MPI_COMM_WORLD, 30, ier )
        call MPI_FINALIZE(ier)
        stop 'Abort - syncInput'
      endif
      end subroutine


!-----------------------------------------------------------------------
      subroutine syncHits(rank,nprocesses,hitcount,nfree)
!-----------------------------------------------------------------------
! sync with main process and print to file
      use mpi
      implicit none
      integer,intent(in):: rank,nprocesses
      integer,intent(inout):: nfree
      integer,intent(inout):: hitcount(nfree)
      ! local parameters
      integer:: ier,i,k,status(MPI_STATUS_SIZE)
      integer:: hitcount_proc(nfree)

      ! nothing to do for single processor job
      if (nprocesses < 2) return

      !print *,'rank',rank,nfree,nz

      ! get row data and print to file
      if (rank == 0) then
        ! main process receives
        do i = 1,nprocesses-1
          ! get data
          !print *,'getting from',i
          call MPI_Recv(hitcount_proc,nfree,MPI_INTEGER,i,i,MPI_COMM_WORLD,status,ier)
          if (ier /= 0) then
            print *,'Error: hitcount_proc',rank,ier
            call stopProgram("syncHits    ")
          endif

          ! add to main process hitcounts
          do k = 1,nfree
            hitcount(k) = hitcount(k)+hitcount_proc(k)
          enddo
        enddo
      else
        ! secondary process sends information
          call MPI_Send(hitcount,nfree,MPI_INTEGER,0,rank,MPI_COMM_WORLD,ier)
          if (ier /= 0) then
            print *,'Error: hitcount process',rank,ier
            call stopProgram("syncHits    ")
          endif
      endif
      end subroutine

!-----------------------------------------------------------------------
      subroutine syncFlag(rank,nprocesses,flag)
!-----------------------------------------------------------------------
! syncs a (logical) flag between all processes,
! assumes that the flag by default is .false.
      use mpi
      implicit none
      integer,intent(in):: rank,nprocesses
      logical,intent(inout):: flag
      ! local parameters
      integer:: ier
      logical:: sum_flag

      ! nothing to do for single processor job
      if (nprocesses < 2) return

      ! sees if a process has the flag set to true ( logical or )
      sum_flag = .false.
      call MPI_AllReduce(flag,sum_flag,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ier)
      if (ier /= 0) then
        print *,'Error: rank ',rank,' failed in allReduce'
        call stopProgram("syncFlag error  " )
      endif

      ! returns flag which should be the same for all processes now
      flag = sum_flag

      end subroutine


!-----------------------------------------------------------------------
      subroutine syncProcesses()
!-----------------------------------------------------------------------
! syncs a (logical) flag between all processes,
! assumes that the flag by default is .false.
      use mpi
      implicit none
      ! local parameters
      integer:: ier

      ! wait until all processes reached this point
      call MPI_Barrier( MPI_COMM_WORLD, ier )

      end subroutine

