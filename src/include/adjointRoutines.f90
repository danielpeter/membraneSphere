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
      subroutine allocateAdjointArrays()
!-----------------------------------------------------------------------
! allocates dynamically the arrays:   adjointKernel, adjointSource, wavefieldForward
      use propagationStartup;use parallel;use adjointVariables;use verbosity
      implicit none
      real(WP)::memory
      integer:: ierror
      
      ! allocate arrays
      if( allocated(adjointKernel) ) deallocate(adjointKernel)
      allocate(adjointKernel(numVertices), stat=ierror)
      if( ierror .ne. 0 ) call stopProgram('error allocating adjoint kernel arrays    ')
      ! initialize array
      adjointKernel(:) = 0.0_WP

      ! velocity seismogram
      if( allocated(adjointSource) ) deallocate(adjointSource)
      allocate(adjointSource(2,numofTimeSteps),stat=ierror)
      if( ierror .ne. 0 ) call stopProgram('error allocating adjoint source arrays   ')
      adjointSource(:,:)= 0.0_WP
      
      ! machine memory holds for 2 GB RAM: 
      !   level 6: numVertices=122'882, numofTimeSteps~500, double precision 8 byte -> needs ~ 500 MB per wavefield, still o.k.
      !   level 7: numVertices=491'522, numofTimeSteps~100, dp 8 byte -> needs ~ 3.8 GB ! per wavefield, too big
      !if( subdivisions .gt. 6 ) then
      !  storeAsFile=.true.
      !endif

      ! for faster computation try to store values in arrays than files
      storeAsFile = .false.   
      ! try to allocate wavefield arrays
      if( allocated(wavefieldForward) ) deallocate(wavefieldForward)
      allocate(wavefieldForward(numDomainVertices,numofTimeSteps),stat=ierror)
      if( ierror .ne. 0 ) then
        ! required harddisk memory shouldn't be bigger than 10 GB
        memory=numVertices*numofTimeSteps*8
        if( DEBUG ) print*,'required memory:',memory
        if( memory .gt. 10.0e9 ) call stopProgram('too much memory use to store data!   ')
      
        storeAsFile = .true.
        ! debug
        if( DEBUG) print*,'   allocating forward wavefield failed...',rank          
      else
        !initialize
        wavefieldForward(:,:)=0.0_WP
        !wavefieldAdjoint(:,:)=0.0_WP
      endif
      
      ! storage of adjoint wavefields
      if( storeAsFile ) then
        ! free up space if they were allocated successfully
        if(allocated(wavefieldForward)) deallocate(wavefieldForward)
        
        ! console output
        if( MASTER .and. VERBOSE) then
          print*
          print*,'adjoint wavefields will be stored as data-files, slower performance expected.'
          print*
        endif
      endif
      end

!-----------------------------------------------------------------------
      subroutine storeForwardDisplacements(timestep,index)
!-----------------------------------------------------------------------
      use adjointVariables;use displacements;use parallel;use propagationStartup
      use verbosity;use griddomain
      implicit none
      integer:: timestep,index,vertex,n,ierror
      character*8:: timestepstr
      character*3:: rankstr
      
      !debug
      if(DEBUG) print*,'    storing forward field...',rank,numDomainVertices,index
    
      if( storeAsFile ) then
        ! collect displacements
        !if( PARALLELSEISMO) then
        !  call collectFullNewdisplacement()
        !endif               
        !if( MASTER ) then     
        !endif
        
        ! output to file for each process
        write(timestepstr,'(i8.7)')timestep
        write(rankstr,'(i3.3)') rank            
        open(10,file=datadirectory(1:len_trim(datadirectory))//'wavefield_'//timestepstr//'.rank'//rankstr//'.dat',&
              access='direct',recl=8,iostat=ierror)
        !open(10,file=datadirectory(1:len_trim(datadirectory))//'wavefield_'//timestepstr//'.rank'//rankstr//'.dat',iostat=ierror)
        if( ierror .ne. 0 ) call stopProgram('could not open file: wavefield_'//timestepstr//'.rank'//rankstr//'.dat ....     ')
        ! fill in displacements
        !do i=1,numVertices
        !  call getSphericalCoord_Lat(i,lat,lon)
        !  write(10,*) lon,lat,newdisplacement(i)
        !enddo
        do n=1, numDomainVertices
          ! choose vertex
          if( PARALLELSEISMO ) then
            vertex=domainVertices(n)
          else
            vertex=n
          endif
          !call getSphericalCoord_Lat(vertex,lat,lon)
          !write(10,*) lon,lat,newdisplacement(vertex)
          write(10,rec=n) newdisplacement(vertex)
        enddo
        close(10)
      else
        ! store in wavefield array
        do n=1, numDomainVertices
          ! choose vertex
          if( PARALLELSEISMO ) then
            vertex=domainVertices(n)
          else
            vertex=n
          endif
          wavefieldForward(n,index) = newdisplacement(vertex)
        enddo            
      endif
      end


!-----------------------------------------------------------------------
      subroutine getAdjointSource()
!-----------------------------------------------------------------------
! determine adjoint source
! (Tromp et al., 2005, sec 4.1 eq (45) )
      use propagationStartup; use parallel; use splineFunction
      use verbosity; use filterType; use adjointVariables
      implicit none
      real(WP)::seismo(numofTimeSteps),seismoDisplacement(numofTimeSteps)
      integer::n,i,ierror
      
      ! console output
      if( MASTER .and. VERBOSE) print*,'getting adjoint source...'            

      ! synchronize seismogram at receiver location
      if( PARALLELSEISMO ) call syncReceiverSeismogram()
                  
      ! for spline representation
      !allocate(X(numofTimeSteps),Y(numofTimeSteps),Q(3,numofTimeSteps),F(3,numofTimeSteps),stat=ierror)
      !if( ierror .ne. 0) call stopProgram('getAdjointSource() - error allocating spline arrays    ')
      
      ! filter around corner frequency
      !if( FILTERSEISMOGRAMS ) then
      !  if( MASTER .and. VERBOSE) print*,'    using filtered receiver seismogram...' 
      !  if( .not. MASTER ) beVerbose=.false.     
      !  call dofilterSeismogram(seismogramReceiver,numofTimeSteps)
      !endif
            
      ! get receiver's velocity seismogram (first time derivative of recorded seismogram at receiver)
      seismo(:)=0.0_WP
      seismoDisplacement(:)=seismogramReceiver(2,:)
      call getFirstDerivative(seismoDisplacement,seismo,dt,numofTimeSteps)

      ! file output
      if( MASTER .and. VERBOSE) then
        print*,'    printing to file: '//datadirectory(1:len_trim(datadirectory))//'seismo.velocity.dat'
        open(10,file=datadirectory(1:len_trim(datadirectory))//'seismo.velocity.dat')
        do n=1,numofTimeSteps
          write(10,*) seismogramReceiver(1,n),seismo(n)
        enddo
        close(10)
      endif
      
      ! source is the velocity seismogram reversed in time
      do i=1,numofTimeSteps
        adjointSource(1,i)=seismogramReceiver(1,numofTimeSteps-i+1)
        adjointSource(2,i)=seismo(numofTimeSteps-i+1)
      enddo
      
      ! apply hanning window  to smooth adjoint source ends
      call taperSeismogram(adjointSource,numofTimeSteps,numofTimeSteps,.false.)
      
      ! normalization factor (Tromp et al., 2005, sec. 4.1 eq. (42) )          
      ! get second time derivative
      seismo(:)=0.0
      call getSecondDerivative(seismoDisplacement,seismo,dt,numofTimeSteps)
      
      ! file output
      if( MASTER .and. VERBOSE ) then
        print*,'    printing to file: '//datadirectory(1:len_trim(datadirectory))//'seismo.acceleration.dat'
        open(10,file=datadirectory(1:len_trim(datadirectory))//'seismo.acceleration.dat')
        do n=1,numofTimeSteps
          write(10,*) seismogramReceiver(1,n),seismo(n)
        enddo
        close(10)
      endif
      
      ! normalization factor is the time integral
      normFactor=0.0
      call getIntegral(seismo,seismoDisplacement,normFactor,dt,numofTimeSteps)

      ! check for division by zero
      if( abs(normFactor) .lt. 1e-6 ) call stopProgram('normalization factor zero!   ')

      ! normalize source
      adjointSource(2,:)=adjointSource(2,:)/normFactor

      ! adjoint source localization is at the receiver
      adjointSourceVertex=receiverVertex      
      
      ! console & file output
      if( MASTER .and. VERBOSE ) then
        ! store adjoint source seismogram
        print*,'    printing to file: '//datadirectory(1:len_trim(datadirectory))//'seismo.adjointSource.dat'
        open(10,file=datadirectory(1:len_trim(datadirectory))//'seismo.adjointSource.dat')
        do n=1,numofTimeSteps
          write(10,*) adjointSource(1,n),adjointSource(2,n)
        enddo
        close(10)
        ! normalization output
        print*,'    normalization factor: ',normFactor
      endif
                        
      ! precalculate the second time derivatives of all vertices belonging to this process
      if( ADJOINT_ONTHEFLY .and. PRECALCULATE_DERIVATIVES ) call precalculateSecondDerivatives()      
      end
      
!-----------------------------------------------------------------------
      subroutine precalculateSecondDerivatives()
!-----------------------------------------------------------------------
! calculates the second derivatives of the forward wavefield seismogram at each vertex
!
! returns: stores the second derivatives in the wavefieldForward array
      use propagationStartup; use parallel; use griddomain; use filterType; use adjointVariables;use verbosity
      implicit none
      integer::n,vertex,index
      real(WP)::seismo(2,numofTimeSteps),seismoTmp(numofTimeSteps)
      
      ! console output
      if( MASTER .and. VERBOSE) print*,'precalculate time derivatives...'
            
      ! every vertex has its own seismogram to derive
      seismo(1,:)=seismogramReceiver(1,:) ! time is the same as from the receiver station seismogram
      if( MASTER .and. VERBOSE .and. FILTERSEISMOGRAMS ) print*,'    using filtered forward seismograms...'
      
      do n=1,numDomainVertices
        ! get cell vertex
        if( PARALLELSEISMO ) then
          vertex=domainVertices(n)
        else
          vertex=n
        endif
                
        ! get corresponding seismogram
        if( storeAsFile ) then
          call getForwardWave(n,seismoTmp)
        else
          seismoTmp(:)=wavefieldForward(n,:)
        endif
        seismo(2,:)=seismoTmp(:)            
                        
        ! filter seismogram        
        if( FILTERSEISMOGRAMS ) then
          call dofilterSeismogram(seismo,numofTimeSteps)        
          seismoTmp(:)=seismo(2,:)
        endif

        ! compute second time derivative (by spline)                        
        call getSecondDerivative(seismo(2,:),seismoTmp,dt,numofTimeSteps)
        
        ! compute second time derivative (by central-difference time scheme)
        !do index=2,numofTimeSteps-1
        !  seismoSecondDerivative(vertex,index)=(seismoTmp(2,index+1) &
        !                - 2.0_WP*seismoTmp(2,index) &
        !                + seismoTmp(2,index-1))/dt2
        !enddo        
        
        ! store as new seismogram
        if( storeAsFile ) then
          call setForwardWave(n,seismoTmp)
        else
          wavefieldForward(n,:)=seismoTmp(:)
        endif
        
      enddo
      end
            
                        
!-----------------------------------------------------------------------
      subroutine backwardIteration()
!-----------------------------------------------------------------------
! (back-propagation) simlation with the adjoint source
      use propagationStartup; use cells; use phaseVelocityMap; use parallel; use displacements;use verbosity
      use loop; use phaseBlockData; use griddomain; use splineFunction; use adjointVariables; use precision
      implicit none
      real(WP):: u_t, u_tplus1, u_tminus1,forcing,forcingRef,forceAdjointSource,arrivalTime,vertexCellArea
      real(WP):: D2,time,discreteLaplacian, precalc_discreteLaplacian,precalc_backdiscreteLaplacian,distance
      integer:: n,i,k,timestep,vertex,index,getDomain,vertexIndex,ierror
      external:: forceAdjointSource,discreteLaplacian,precalc_discreteLaplacian,precalc_backdiscreteLaplacian,getDomain
      character*8:: timestepstr
      character*3:: rankstr
      character*4:: timestr
      
      ! benchmark
      if( MASTER .and. VERBOSE ) benchstart = MPI_WTIME()      
      
      if( ADJOINT_ONTHEFLY ) then
        ! simultaneous backward displacement (same calculation as forward one but reversed) initial start
        !backwardDisplacement_old(:) = newdisplacement(:)
        !backwardDisplacement(:) = newdisplacement(:)
        !backwardNewDisplacement(:) = displacement(:)      

        ! look for vertex in middle of source/receiver
        call findVertex(sourceLat+nint((receiverLat-sourceLat)/2.0),sourceLon+nint((receiverLon-sourceLon)/2.0),midpointVertex)

        if( MASTER .and. VERBOSE ) print*,'   storing file:',datadirectory(1:len_trim(datadirectory))//'seismo.integral_source.dat'
        open(adjSourceFileID,file=datadirectory(1:len_trim(datadirectory))//'seismo.integral_source.dat')
        if( MASTER .and. VERBOSE ) print*,'   storing file:',datadirectory(1:len_trim(datadirectory))//'seismo.integral_rec.dat'
        open(adjRecFileID,file=datadirectory(1:len_trim(datadirectory))//'seismo.integral_rec.dat')
        if( MASTER .and. VERBOSE ) print*,'   storing file:',datadirectory(1:len_trim(datadirectory))//'seismo.integral_midpoint.dat'
        open(adjMidpointFileID,file=datadirectory(1:len_trim(datadirectory))//'seismo.integral_midpoint.dat')
      else
        if( .not. storeAsFile ) then      
          ! allocate adjoint wavefield 
          if( .not. allocated(wavefieldAdjoint) ) then
            allocate(wavefieldAdjoint(numDomainVertices,numofTimeSteps),stat=ierror)
            if( ierror .ne. 0 ) then
              call stopProgram('error allocating adjoint wavefield    ')
            endif
            !initialize
            wavefieldAdjoint(:,:)=0.0_WP     
          endif
        endif   
      endif
            
      ! reset displacements for adjoint propagation
      displacement_old(:)=0.0_WP
      displacement(:)=0.0_WP
      newdisplacement(:)=0.0_WP      
              
      ! time iteration of displacements
      index=0
      do timestep=lasttimestep,firsttimestep,-1
        ! model time
        time = timestep*dt      
        !debug 
        if(DEBUG) then
          print*,'    time:',rank,time
        endif             
        
        ! time index
        index=index+1
        
        ! swap displacement arrays
        displacement_old(:) = displacement(:)
        displacement(:) = newdisplacement(:)
        
        !if( ADJOINT_ONTHEFLY ) then
        !  backwardDisplacement_old(:) = backwardDisplacement(:)
        !  backwardDisplacement(:) = backwardNewdisplacement(:)
        !endif
                                                            
        ! propagate only corresponding vertices
        do n=1, numDomainVertices
          ! choose vertex
          if( PARALLELSEISMO ) then
            vertex=domainVertices(n)
          else
            vertex=n
          endif
          
          ! spherical Laplacian    
          if( PRECALCULATED_CELLS ) then      
            D2 = precalc_discreteLaplacian(vertex) !uses displacement(..) field
          else
            D2 = discreteLaplacian(vertex)
          endif
          
          ! calculate new displacement in time
          u_t = displacement(vertex)
          u_tminus1 = displacement_old(vertex)
          
          forcing = forceAdjointSource(vertex,index)

          ! calculate phase velocity square 
          cphase2 = phaseVelocitySquare(vertex)
                           
          !propagation step                                                     
          u_tplus1=u_t+u_t-u_tminus1+cphase2*dt2*(D2+forcing)
                        
          ! iterated displacements
          newdisplacement(vertex) = u_tplus1
          
          if( ADJOINT_ONTHEFLY ) then
          !  ! simulate backward iteration in parallel
          !  ! spherical Laplacian          
          !  D2 = precalc_backdiscreteLaplacian(vertex) !uses displacement(..) field
          !  
          !  ! calculate new displacement in time
          !  u_t = backwardDisplacement(vertex)
          !  u_tminus1 = backwardDisplacement_old(vertex)            
          !  
          !  ! calculate phase velocity square 
          !  cphase2 = phaseVelocitySquare(vertex)
          !                   
          !  !propagation step ( no force applies for backward calculation of previous forward field)
          !  u_tplus1=u_t+u_t-u_tminus1+cphase2*dt2*(D2)
          !                
          !  ! iterated displacements
          !  backwardNewdisplacement(vertex) = u_tplus1
          !  
            ! compute adjoint kernel value step by step 
            call getAdjointKernel(vertex,timestep)
          endif
        enddo

        !print*,'synchronizing...',rank,nprocesses,index,time
        if( PARALLELSEISMO ) then
          ! synchronize new displacement arrays
          call syncNewdisplacement()
          ! synchronize new displacement arrays
          !if( ADJOINT_ONTHEFLY ) then
          !  call syncBackwardNewdisplacement(rank,nprocesses)
          !endif
        endif

        !file output for simulation snapshots
        if(SIMULATIONOUTPUT) then
          if( mod(timestep,10) .eq. 0 .and. time .ge. 0.0 .and. MASTER ) then
            write(timestr,'(i4.4)') timestep
            open(10,file=datadirectory(1:len_trim(datadirectory))//'simulationAdjoint.'//timestr//'.dat')        
            do n=1, numVertices
              write(10,'(4f18.6)')(vertices(n,k),k=1,3),newdisplacement(n)
            enddo
            close(10) 
          endif          
        endif

        if( .not. ADJOINT_ONTHEFLY ) then
          !debug
          if(DEBUG) print*,'    storing adjoint field...',rank,numDomainVertices,index

          ! save displacements at each time step for future adjoint calculation        
          if( storeAsFile ) then        
            ! output to file for each process
            write(timestepstr,'(i8.7)')timestep
            write(rankstr,'(i3.3)') rank            
            open(10,file=datadirectory(1:len_trim(datadirectory))//'wavefieldAdj_'//timestepstr//'.rank'//rankstr//'.dat',&
                  access='direct',recl=8,iostat=ierror)
            !open(10,file=datadirectory(1:len_trim(datadirectory))//'wavefield_'//timestepstr//'.rank'//rankstr//'.dat',iostat=ierror)
            if( ierror .ne. 0 ) call stopProgram('could not open file wavefieldAdj_'//timestepstr//'.rank'//rankstr//'.dat ....     ')
            ! fill in displacements
            do n=1, numDomainVertices
              ! choose vertex
              if( PARALLELSEISMO ) then
                vertex=domainVertices(n)
              else
                vertex=n
              endif
              write(10,rec=n) newdisplacement(vertex)
            enddo
            close(10)
          else
            ! store in adjoint wavefield array
            do n=1, numDomainVertices
              ! choose vertex
              if( PARALLELSEISMO ) then
                vertex=domainVertices(n)
              else
                vertex=n
              endif
              wavefieldAdjoint(n,index) = newdisplacement(vertex)
            enddo            
          endif           
        endif ! .not. ADJOINT_ONTHEFLY
      enddo !timestep

      ! on the fly calculation needs to be scaled finally
      if( ADJOINT_ONTHEFLY ) then
        ! close output files
        close(adjSourceFileID)
        close(adjRecFileID)
        close(adjMidpointFileID)
      
        ! cell area in [rad^2]
        vertexCellArea= cellAreas(receiverVertex)/(EARTHRADIUS*EARTHRADIUS)
      
        ! calculate the reference travel time
        call greatCircleDistance(vertices(sourceVertex,:),vertices(receiverVertex,:),distance)
        arrivalTime=distance*EARTHRADIUS/cphaseRef      

        ! scale each entry of the adjoint kernel        
        adjointKernel(:)= adjointKernel(:)/(arrivalTime*vertexCellArea)      
      endif

      !debug file output for seismogram at source again
      if(DEBUG) then
        if( rank .eq. getDomain(sourceVertex) .and. .not. storeAsFile .and. .not. ADJOINT_ONTHEFLY ) then
          ! index for source vertex location in the wavefield files
          if( PARALLELSEISMO ) then
            do n=1,numDomainVertices
              if( sourceVertex .eq. domainVertices(n) ) then
                vertexIndex=n
                exit
              endif
            enddo
          else
            vertexIndex=vertex
          endif
          ! write to file  
          print*,'    printing to file: '//datadirectory(1:len_trim(datadirectory))//'seismo.atSource.dat'          
          open(10,file=datadirectory(1:len_trim(datadirectory))//'seismo.atSource.dat')
          do n=1,numofTimeSteps
            write(10,*) (lasttimestep-n+1)*dt+dt,wavefieldAdjoint(vertexIndex,n)
          enddo
          close(10)
        endif
      endif
      
      ! free memory to be able to have some to collect data
      !deallocate(cellEdgesLength,cellCenterDistances,cellNeighbors)
      !deallocate(displacement,displacement_old,newdisplacement)
      !if( PARALLELSEISMO ) deallocate(boundaries,sendDisp,receiveDisp,domainNeighbors)    
      !if( ADJOINT_ONTHEFLY ) deallocate(backwardDisplacement,backwardDisplacement_old,backwardNewdisplacement)
      
      ! benchmark output
      if( MASTER .and. VERBOSE ) then
        benchend = MPI_WTIME()
        print*,'    benchmark seconds:',benchend-benchstart
        print*
      endif            
      end
      
            
!-----------------------------------------------------------------------
      subroutine getAdjointKernel(vertex,timestep)
!-----------------------------------------------------------------------
! determine sensitivity kernel value when calculating "on the fly"
!
! return: stores values in the adjointKernel array
      use propagationStartup; use displacements; use phaseVelocityMap; use traveltime 
      use adjointVariables; use parallel; use griddomain; use cells; use verbosity
      implicit none
      integer:: vertex,timestep,n,index,vertexIndex
      real(WP):: seismo(3),seismoAdjoint(3),seismoTmp(numofTimeSteps),timewindow
      real(WP):: kernelVal,val1,val2,for1,for2,for3,for4,derivativeActual,derivativeNext,kernelfactor
      logical:: bySpline=.false.
      
      ! current time step index
      index= timestep - firsttimestep  + 1

      ! time window for integration
      if( WINDOWEDINTEGRATION ) then
        if( timestep*dt .lt. 0.0 .and. timestep*dt .gt. arrivalTime ) then
          return
        endif
      endif
      
      ! index for vertex location in the wavefield files
      if( PARALLELSEISMO ) then
        do n=1,numDomainVertices
          if( vertex .eq. domainVertices(n) ) then
            vertexIndex=n
            exit
          endif
        enddo
      else
        vertexIndex=vertex
      endif
      
      ! calculate kernel value with actual available displacement values
      if( bySpline ) then
        ! time integral calculated as a sum of discrete rectangles with size dt
        !if( (index .gt. numofTimeSteps-2) .or. (index .lt. 1) ) then 
        !  return        
        !endif
        
        ! get seismogram
        !if( storeAsFile ) then
        !  call getForwardWave(vertexIndex,seismoTmp)          
        !else
        !  seismoTmp(:)=wavefieldForward(vertexIndex,:)
        !endif        
        !seismo(1)=seismoTmp(index)
        !seismo(2)=seismoTmp(index+1)
        !seismo(3)=seismoTmp(index+2)      
        
        seismo(1)=backwardNewdisplacement(vertex)
        seismo(2)=backwardDisplacement(vertex)
        seismo(3)=backwardDisplacement_old(vertex)
        
        seismoAdjoint(1)=newdisplacement(vertex)
        seismoAdjoint(2)=displacement(vertex)
        seismoAdjoint(3)=displacement_old(vertex)
        
        ! time derivative of seismo      
        if( .not. PRECALCULATE_DERIVATIVES ) then     
          if( seismo(1) .eq. 0.0 .and. seismo(2) .eq. 0.0 .and. seismo(3) .eq. 0.0) then
            continue
          else
            call getSecondDerivative(seismo,seismo,dt,3)
          endif
        endif
        
        ! kernel value is the time integral
        kernelVal=0.0
        val1=seismo(1)*seismoAdjoint(1)
        val2=seismo(2)*seismoAdjoint(2)        
        kernelVal=0.5_WP*(val1+val2)*dt         
      else
        ! time integral calculated as a sum of discrete rectangles with size dt
        !if( (index .gt. numofTimeSteps-1).or.(index .lt. 2 ) ) then
        !  return
        !endif        
        
        ! get seismogram
        !seismo(1)=backwardNewdisplacement(vertex)
        !seismo(2)=backwardDisplacement(vertex)
        !seismo(3)=backwardDisplacement_old(vertex)
        ! get adjoint seismogram
        !seismoAdjoint(1)=newdisplacement(vertex)
        seismoAdjoint(2)=displacement(vertex)
        !seismoAdjoint(3)=displacement_old(vertex)
        
        ! time derivative of seismo (by central-differences)
        !val1=(seismo(1)+seismo(3)-2.0_WP*seismo(2))/dt2
        !kernelval=val1*seismoAdjoint(2)*dt
        
        ! get seismogram
        if( storeAsFile ) then
          call getForwardWave(vertexIndex,seismoTmp)          
        else
          seismoTmp(:)=wavefieldForward(vertexIndex,:)
        endif
                
        ! compute second time derivative (by central-difference time scheme)
        if( PRECALCULATE_DERIVATIVES ) then
          val1=seismoTmp(index)*newdisplacement(vertex)
          val2=seismoTmp(index+1)*displacement(vertex)  
        else
          for1=seismoTmp(index-1)
          for2=seismoTmp(index)
          for3=seismoTmp(index+1)
        !  for4=seismoTmp(index+2)
        
        !  !seismoSecondDerivative(vertex,index)=(for1 + for3 - 2.0_WP*for2)/dt2
        !  !seismoSecondDerivative(vertex,index+1)=(for2 + for4 - 2.0_WP*for3)/dt2
        !  derivativeActual=(for1 + for3 - 2.0_WP*for2)/dt2
        !  derivativeNext=(for2 + for4 - 2.0_WP*for3)/dt2
        
          ! kernel value increment
        !  !val1= seismoSecondDerivative(vertex,index)*newdisplacement(vertex)
        !  !val2= seismoSecondDerivative(vertex,index+1)*displacement(vertex)
        !  val1= derivativeActual*newdisplacement(vertex)
        !  val2= derivativeNext*displacement(vertex)          
          val1=(for1 + for3 - 2.0_WP*for2)/dt2
        endif        
        !kernelVal=0.5_WP*(val1+val2)*dt
        kernelVal=val1*seismoAdjoint(2)*dt
        
      endif

      ! calculate kernel factor
      kernelfactor=2.0_WP/(phaseVelocitySquare(vertex))      
      kernelVal=kernelfactor*kernelVal

      ! sum up to build adjointKernel value
      adjointKernel(vertex)=adjointKernel(vertex)+kernelVal
      
      !debug
      !if(DEBUG) print*,'adjoint:',adjointKernel(vertex),kernelVal,vertex

      ! output to files (timestep: correct +1 by back-propagation and to account for values stored in displacement() and not newdisplacement() )
      if( vertex .eq. sourceVertex ) then
          write(adjSourceFileID,'(6e16.4e3)') (timestep+1)*dt,seismo(2),(timestep+1)*dt,seismoAdjoint(2),adjointKernel(vertex),kernelVal/kernelFactor
      endif          
      
      if( vertex .eq. receiverVertex ) then
          write(adjRecFileID,'(6e16.4e3)') (timestep+1)*dt,seismo(2),(timestep+1)*dt,seismoAdjoint(2),adjointKernel(vertex),kernelVal/kernelFactor
      endif          

      if( vertex .eq. midpointVertex ) then
          write(adjMidpointFileID,'(6e16.4e3)') (timestep+1)*dt,seismo(2),(timestep+1)*dt,seismoAdjoint(2),adjointKernel(vertex),kernelVal/kernelFactor
      endif                    
      
      end   
      
!-----------------------------------------------------------------------
      subroutine frechetKernel()
!-----------------------------------------------------------------------
! determine adjoint kernel value when calculating NOT "on the fly"
!
! return: stores values in the adjointKernel() array
      use propagationStartup; use parallel; use griddomain; use splineFunction 
      use phaseVelocityMap; use filterType; use verbosity; use cells; use adjointVariables
      implicit none
      logical:: bySpline=.false.
      integer:: timestep,index,indexadjoint,vertex,n,i
      real(WP):: kernelVal,val1,val2,derivativeActual,derivativeNext,kernelfactor
      real(WP):: arrivalTime,distance,vertexCellArea,vertexCellAreaRad,timewindow,for1,for2,for3
      real(WP)::seismo(2,numofTimeSteps),seismoAdjoint(2,numofTimeSteps),seismoTmp(numofTimeSteps)      
      
      ! console output
      if( MASTER .and. VERBOSE) print*,'calculating kernel values...'      
            
      ! initialize with time (dealing with newdisplacements means at time steps t+dt)
      do index=1,numofTimeSteps
        seismo(1,index)=(firsttimestep+index-1)*dt+dt
        seismoAdjoint(1,index)=(lasttimestep-index+1)*dt+dt
      enddo

      ! calculate the reference travel time 
      ! attention: for a heterogeneous background it takes here the PREM value as well.
      !                this should be considered when inverting and using these kernels.
      call greatCircleDistance(vertices(sourceVertex,:),vertices(receiverVertex,:),distance)
      arrivalTime=distance*EARTHRADIUS/cphaseRef      
      
      ! cell area 
      vertexCellArea=cellAreas(receiverVertex)
      ! convert cell area into [rad^2]
      vertexCellAreaRad=vertexCellArea/(EARTHRADIUS*EARTHRADIUS)
            
      ! console output
      if( MASTER .and. VERBOSE) then
        print*,'    reference travel time [s]:',arrivalTime      
        print*,'    receiver cell area [km2]:',vertexCellArea !for [rad2]: vertexCellArea/(EARTHRADIUS*EARTHRADIUS)      
        print*,'time integration:'
        print*,'    starting seismogram at:',seismo(1,1)
        print*,'    ending seismogram at:',seismo(1,numofTimeSteps)      
        if( WINDOWEDINTEGRATION ) print*,'    windowing between:',0.0,arrivalTime
        print*

      endif  
            
      ! look for vertex in middle of source/receiver
      call findVertex(sourceLat+nint((receiverLat-sourceLat)/2.0),sourceLon+nint((receiverLon-sourceLon)/2.0),midpointVertex)
      
      ! determine kernel value at each grid point for the corresponding processor domain
      do n=1, numDomainVertices        
        ! choose vertex
        if( PARALLELSEISMO ) then
          vertex=domainVertices(n)
        else
          vertex=n
        endif

        ! get forward seismogram
        if( storeAsFile ) then
          call getForwardWave(n,seismoTmp)          
          seismo(2,:)=seismoTmp(:)
        else
          seismo(2,:)=wavefieldForward(n,:)
        endif

        ! get adjoint seismogram
        if( storeAsFile ) then
          call getAdjointWave(n,seismoTmp)  
          seismoAdjoint(2,:)=seismoTmp(:)        
        else
          seismoAdjoint(2,:)=wavefieldAdjoint(n,:)
        endif
                        
        ! filter (and taper) the seismograms 
        if( FILTERSEISMOGRAMS ) then
          ! reset the verbosity
          if( MASTER .and. VERBOSE) then
            beverbose = .true.
          else
            beverbose = .false.
          endif        
          
          ! filter
          if( beverbose ) print*,'   filtering seismograms...'
          call dofilterSeismogram(seismo,numofTimeSteps)
          call dofilterSeismogram(seismoAdjoint,numofTimeSteps)
          beverbose=.false.
        endif

        ! kernel value
        kernelVal=0.0_WP   
        seismoTmp(:)=0.0_WP     
        if( bySpline ) then        
          ! compute second derivative of forward seismogram       
          call getSecondDerivative(seismo(2,:),seismo(2,:),dt,numofTimeSteps)
                    
          ! kernel value is the time integral  
          do index=1,numofTimeSteps-1
            ! time window for integration
            timewindow=1.0_WP            
            if( WINDOWEDINTEGRATION ) then
              if( seismo(1,index) .lt. 0.0 .and. seismo(1,index) .gt. arrivalTime ) then
                timewindow=0.0_WP
              endif
            endif
          
            derivativeActual=seismo(2,index)
            derivativeNext=seismo(2,index+1)
  
            ! kernel value is the time integral (by a sum of discrete rectangles with size dt)                
            val1=derivativeActual*seismoAdjoint(2,numofTimeSteps-index+1)
            val2=derivativeNext*seismoAdjoint(2,numofTimeSteps-index+1-1)                    
            
            kernelVal=kernelVal+0.5_WP*(val1+val2)*dt*timewindow
            
            ! check & store temporary
            if( index .gt. 1 ) then
              if( val1 .ne. seismoTmp(index) ) print*,'kernel values:',val1,seismoTmp(index)
            endif
            seismoTmp(index)=val1
            seismoTmp(index+1)=val2
          enddo                            
          
          ! integral by spline representation: little bit better accuracy (influence on ~ 4. digit; neglectable), little bit slower
          !reverse adjoint
          do index=1,numofTimeSteps/2
            val1=seismoAdjoint(2,index)
            val2=seismoAdjoint(2,numofTimeSteps-index+1)
            seismoAdjoint(2,index)=val2
            seismoAdjoint(2,numofTimeSteps-index+1)=val1
          enddo          
          call getIntegral(seismo(2,:),seismoAdjoint(2,:),kernelVal,dt,numofTimeSteps)          
        else          
          ! central differences - scheme          
          do index=1,numofTimeSteps-1          
            ! time window for integration
            timewindow=1.0_WP            
            if( WINDOWEDINTEGRATION ) then
              if( seismo(1,index) .lt. 0.0 .and. seismo(1,index) .gt. arrivalTime ) then
                timewindow=0.0_WP
              endif
            endif
            
            ! compute second derivative of forward seismogram                
            if( index .eq. 1 ) continue
            
            for1=seismo(2,index-1)
            for2=seismo(2,index)
            for3=seismo(2,index+1)

            val1=(for1 + for3 - 2.0_WP*for2)/dt2
            
            ! adjoint wavefield index
            ! careful: the index of the adjoint wavefield corresponds to T-t, therefore (numofTimeSteps - index) for this second derivative
            indexadjoint=numofTimeSteps-index            
            
            kernelVal=kernelVal+val1*seismoAdjoint(2,indexadjoint)*dt
            
            !seismoSecondDerivative(vertex,index)=(seismo(2,index+1)-2.0_WP*seismo(2,index)+seismo(2,index-1))/dt2
            !seismoSecondDerivative(vertex,index+1)=(seismo(2,index+2)-2.0_WP*seismo(2,index+1)+seismo(2,index))/dt2
            !if( index .eq. 1 ) then
            !  derivativeActual=(seismo(2,index+1)-2.0_WP*seismo(2,index)+seismo(2,index))/dt2
            !  !derivativeNext=(seismo(2,index+2)-2.0_WP*seismo(2,index+1)+seismo(2,index))/dt2
            !else if( index .eq. numofTimeSteps-1 ) then
            !  derivativeActual=(seismo(2,index+1)-2.0_WP*seismo(2,index)+seismo(2,index-1))/dt2
            !  !derivativeNext=(seismo(2,index+1)-2.0_WP*seismo(2,index+1)+seismo(2,index))/dt2                        
            !else            
            !  derivativeActual=(seismo(2,index+1)-2.0_WP*seismo(2,index)+seismo(2,index-1))/dt2
            !  !derivativeNext=(seismo(2,index+2)-2.0_WP*seismo(2,index+1)+seismo(2,index))/dt2
            !endif
  
            ! kernel value is the time integral (by a sum of discrete trapezoids with size dt)                
            !val1=seismoAdjoint(2,index)*seismoSecondDerivative(vertex,index)
            !val2=seismoAdjoint(2,index+1)*seismoSecondDerivative(vertex,index+1)                    
            !val1=derivativeActual*seismoAdjoint(2,numofTimeSteps-index+1)
            !val2=derivativeNext*seismoAdjoint(2,numofTimeSteps-index+1-1)                                
            !kernelVal=kernelVal+0.5_WP*(val1+val2)*dt*timewindow
            
            ! check & store temporary
            !if( index .gt. 1 ) then
            !  if( val1 .ne. seismoTmp(index) ) print*,'kernel values:',val1,seismoTmp(index)
            !endif
            !seismoTmp(index)=val1
            !seismoTmp(index+1)=val2            
          enddo
        endif
                                                                        
        ! calculate kernel factor for relative phase shift kernel (units in radians; sign convection depending on time lag definition)
        kernelfactor=2.0_WP/(phaseVelocitySquare(vertex)*arrivalTime*vertexCellAreaRad)
        kernelVal=kernelfactor*kernelVal
        
        ! store in array
        adjointKernel(vertex)=kernelVal
        
        ! console output &  output to files
        if( vertex .eq. sourceVertex .and. VERBOSE ) print*,'    kernel values source:',vertex,kernelVal,kernelfactor,kernelVal/kernelfactor               
        if( vertex .eq. sourceVertex .and. VERBOSE ) then
          if( MASTER .and. VERBOSE ) print*,'   storing file:',datadirectory(1:len_trim(datadirectory))//'seismo.integral_source.dat'
          open(10,file=datadirectory(1:len_trim(datadirectory))//'seismo.integral_source.dat')
          do i=1,numofTimeSteps
            write(10,'(6e16.4e3)') seismo(1,i),seismo(2,i),seismoAdjoint(1,i),seismoAdjoint(2,i),seismo(2,i)*seismoAdjoint(2,i),seismoTmp(i)
          enddo
          close(10)
        endif          
        if( vertex .eq. receiverVertex .and. VERBOSE ) print*,'    kernel values receiver:',vertex,kernelVal,kernelfactor,kernelVal/kernelfactor                
        if( vertex .eq. receiverVertex .and. VERBOSE ) then
          if( MASTER .and. VERBOSE ) print*,'   storing file:',datadirectory(1:len_trim(datadirectory))//'seismo.integral_rec.dat'
          open(10,file=datadirectory(1:len_trim(datadirectory))//'seismo.integral_rec.dat')
          do i=1,numofTimeSteps
            write(10,'(6e16.4e3)') seismo(1,i),seismo(2,i),seismoAdjoint(1,i),seismoAdjoint(2,i),seismo(2,i)*seismoAdjoint(2,i),seismoTmp(i)
          enddo
          close(10)
        endif          
        if( vertex .eq. midpointVertex .and. VERBOSE ) print*,'    kernel values midpoint (',sourceLat+nint((receiverLat-sourceLat)/2.0),'/',sourceLon+nint((receiverLon-sourceLon)/2.0),'):',vertex,kernelVal,kernelfactor,kernelVal/kernelfactor                
        if( vertex .eq. midpointVertex .and. VERBOSE ) then
          if( MASTER .and. VERBOSE ) print*,'   storing file:',datadirectory(1:len_trim(datadirectory))//'seismo.integral_midpoint.dat'
          open(10,file=datadirectory(1:len_trim(datadirectory))//'seismo.integral_midpoint.dat')
          do i=1,numofTimeSteps
            write(10,'(6e16.4e3)') seismo(1,i),seismo(2,i),seismoAdjoint(1,i),seismoAdjoint(2,i),seismo(2,i)*seismoAdjoint(2,i),seismoTmp(i)
          enddo
          close(10)
        endif                    
      enddo            
      end
      
!-----------------------------------------------------------------------
      subroutine storeAdjointKernel()
!-----------------------------------------------------------------------
! print the adjointKernel()-array values to a file 'adjointkernel.dat'
      use propagationStartup; use parallel; use adjointVariables;use verbosity; use cells
      implicit none
      real(WP):: lat,lon,sum_kern
      integer:: i,ierror
      character*128:: kernelfile

      ! console output
      if( MASTER .and. VERBOSE) print*,'writing values to kernel file...'
            
      ! get complete adjoint kernel array in case we run a single simulation on parallel processors      
      if( PARALLELSEISMO ) call collectAdjointKernel()

      ! remove obsolete files
      if( storeAsFile ) call cleanupWaveFiles()   
      
      ! slaves are done
      if( .not. MASTER) return
      
      ! open kernel file
      kernelfile=datadirectory(1:len_trim(datadirectory))//adjointKernelName(1:len_trim(adjointKernelName))
      open(10,file=kernelfile,iostat=ierror)
      if( ierror .ne. 0) then
        print*,'could not open '//kernelfile
        call stopProgram( 'abort - storeAdjointKernel()    ')
      endif
      
      ! file header
      write(10,*) '# adjoint method - sensitivity kernel'
      write(10,*) '# lon lat kernel vertexID'
      
      !debug
      if(DEBUG) open(30,file=datadirectory(1:len_trim(datadirectory))//'adjointKernel.lon45.dat')
      
      ! write values to file
      sum_kern=0.0_WP
      do i=1,numVertices
        call getSphericalCoord_Lat(i,lat,lon)
        write(10,'(2f8.2,e18.6e3,i12)') lon,lat,adjointKernel(i),i
        
        ! debug a trace
        if(DEBUG) then
          if( abs(lon-45._WP) .lt. 0.5 ) then
            write(30,'(2f8.2,e18.6e3,i12)') lon,lat,adjointKernel(i),i          
          endif        
        endif
        
        ! summate values
        sum_kern=sum_kern+adjointKernel(i)*cellAreas(i)/(EARTHRADIUS*EARTHRADIUS)
      enddo      
      close(10)
      !debug
      if(DEBUG) close(30)
      if(DEBUG) print*,'    trace at longitude 45 written to: ',datadirectory(1:len_trim(datadirectory))//'adjointKernel.lon45.dat'
      
      if( VERBOSE) then
        print*,'    stored kernel values in file: ',kernelfile(1:len_trim(kernelfile))
        print*,'    integrated over sphere: ',sum_kern
        print*
      endif
                     
      end
      
!-----------------------------------------------------------------------
      subroutine getForwardWave(index,seismo)
!-----------------------------------------------------------------------
! get forward seismogram for at given vertex location from the wavefield storage-files
      use propagationStartup; use parallel
      real(WP),intent(OUT)::seismo(numofTimeSteps)
      character*8:: timestepstr
      character*3:: rankstr
      integer:: i,timestep,index,ierror
      
      ! input
      write(rankstr,'(i3.3)') rank 
      i=0                 
      do timestep= firsttimestep, lasttimestep
        i=i+1
        write(timestepstr,'(i8.7)')timestep      
        !if(rank .eq. 1) print*,'    open timestep:',i,timestep
        open(10,file=datadirectory(1:len_trim(datadirectory))//'wavefield_'//timestepstr//'.rank'//rankstr//'.dat',&
            access='direct',recl=8,iostat=ierror)
        !open(10,file=datadirectory(1:len_trim(datadirectory))//'wavefield_'//timestepstr//'.rank'//rankstr//'.dat',iostat=ierror)
        if( ierror .ne. 0 ) call stopProgram('could not open for input, file wavefield_'//timestepstr//'.rank'//rankstr//'.dat ...    ')
        ! fill in displacements
        !if(rank .eq. 1) print*,'    read timestep:',i,timestep        
        read(10,rec=index,iostat=ierror) seismo(i)
        if( ierror .ne. 0 ) then
          print*,'read error:',timestep,index,ierror
          call stopProgram('could not read input from file wavefield_'//timestepstr//'.rank'//rankstr//'.dat ...    ')
        endif
        !if(rank .eq. 1) print*,'    close timestep:',i,timestep        
        close(10)      
      enddo
      
      end
      
      
!-----------------------------------------------------------------------
      subroutine setForwardWave(index,seismo)
!-----------------------------------------------------------------------
! set forward seismogram for at given vertex location to the wavefield storage-files
      use propagationStartup; use parallel
      real(WP),intent(IN)::seismo(numofTimeSteps)
      character*8:: timestepstr
      character*3:: rankstr
      integer:: i,timestep,index,ierror
      
      ! output
      write(rankstr,'(i3.3)') rank 
      i=0                 
      do timestep= firsttimestep, lasttimestep
        i=i+1
        write(timestepstr,'(i8.7)')timestep      
        open(10,file=datadirectory(1:len_trim(datadirectory))//'wavefield_'//timestepstr//'.rank'//rankstr//'.dat',&
            access='direct',recl=8,iostat=ierror)
        !open(10,file=datadirectory(1:len_trim(datadirectory))//'wavefield_'//timestepstr//'.rank'//rankstr//'.dat',iostat=ierror)
        if( ierror .ne. 0 ) call stopProgram('could not open for output, file wavefield_'//timestepstr//'.rank'//rankstr//'.dat ...    ')
        ! fill in displacements
        write(10,rec=index,iostat=ierror) seismo(i)
        if( ierror .ne. 0 ) then
          print*,'write error:',timestep,index,ierror
          call stopProgram('could not write input from file wavefield_'//timestepstr//'.rank'//rankstr//'.dat ...    ')
        endif
        
        close(10)      
      enddo
      
      end
      
!-----------------------------------------------------------------------
      subroutine getAdjointWave(index,seismo)
!-----------------------------------------------------------------------
! get adjoint seismogram for at given vertex location from the wavefield storage-files
      use propagationStartup; use parallel
      real(WP),intent(OUT)::seismo(numofTimeSteps)
      character*8:: timestepstr
      character*3:: rankstr
      integer:: i,timestep,index,ierror
      
      ! input
      write(rankstr,'(i3.3)') rank 
      i=numofTimeSteps                 
      do timestep= firsttimestep, lasttimestep
        write(timestepstr,'(i8.7)')timestep      
        open(10,file=datadirectory(1:len_trim(datadirectory))//'wavefieldAdj_'//timestepstr//'.rank'//rankstr//'.dat',&
            access='direct',recl=8,iostat=ierror)
        !open(10,file=datadirectory(1:len_trim(datadirectory))//'wavefield_'//timestepstr//'.rank'//rankstr//'.dat',iostat=ierror)
        if( ierror .ne. 0 ) call stopProgram('could not open for input, file wavefieldAdj_'//timestepstr//'.rank'//rankstr//'.dat ...    ')
        ! fill in displacements
        read(10,rec=index,iostat=ierror) seismo(i)
        if( ierror .ne. 0 ) then
          print*,'read error:',timestep,index,ierror
          call stopProgram('could not read input from file wavefieldAdj_'//timestepstr//'.rank'//rankstr//'.dat ...    ')
        endif        
        close(10)      
        i=i-1  ! adjoint wavefield: stores T-t meaning that firsttimestep value (t=0) is at seismo(T=numofTimeSteps) and lasttimestep at seismo(1)
      enddo
      
      end
      
!-----------------------------------------------------------------------
      subroutine cleanupWaveFiles()
!-----------------------------------------------------------------------
! delete wavefield files
      use propagationStartup; use adjointVariables; use parallel
      character*8:: timestepstr
      character*3:: rankstr
      integer:: timestep,ierror

      ! get rid of files
      write(rankstr,'(i3.3)') rank 
      do timestep= firsttimestep, lasttimestep
        write(timestepstr,'(i8.7)')timestep      
        ! forward wavefield files
        open(10,file=datadirectory(1:len_trim(datadirectory))//'wavefield_'//timestepstr//'.rank'//rankstr//'.dat',&
            access='direct',recl=8,iostat=ierror)
        !open(10,file=datadirectory(1:len_trim(datadirectory))//'wavefield_'//timestepstr//'.rank'//rankstr//'.dat',iostat=ierror)
        if( ierror .eq. 0 ) close(10,status='DELETE')      
        ! adjoint wavefield files
        open(20,file=datadirectory(1:len_trim(datadirectory))//'wavefieldAdj_'//timestepstr//'.rank'//rankstr//'.dat',&
            access='direct',recl=8,iostat=ierror)
        !open(10,file=datadirectory(1:len_trim(datadirectory))//'wavefield_'//timestepstr//'.rank'//rankstr//'.dat',iostat=ierror)
        if( ierror .eq. 0 ) close(20,status='DELETE')      
      enddo
      end
      