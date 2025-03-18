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
!
! Proper acknowledgement shall be made to the authors of the software in publications and
! presentations resulting from the use of this software:
!
! Peter, D., C. Tape, L. Boschi and J. H. Woodhouse, 2007. Surface wave tomography:
! global membrane waves and adjoint methods, Geophys. J. Int., , 171: p. 1098-1117.

!-----------------------------------------------------------------------
      program adjointMethod
!-----------------------------------------------------------------------
! performs a simulation, where no scatterer is present, and calculates the adjoint source.
! then make a time-reversed simulation and computes the kernel values.
!
! adjoint method: information in tromp et al. (2005)
! finite-difference iteration: information in carl tape thesis (2003), chap 5, (5.7)
      use propagationStartup; use cells; use phaseVelocityMap
      use parallel; use displacements
      use loop; use phaseBlockData; use griddomain; use adjointVariables;use verbosity
      implicit none
      ! local parameters
      integer:: kernel,ierror
      real(WP):: window_start_org,window_end_org

      !-----------------------------------------------------------------------
      ! parameters
      ! most parameters concerning wave propagation set in file Parameter_Input
      ! (& commonModules.f90)
      ! initialize parameters
      Adjoint_Program                 = .true.   ! we calculate kernels by adjoint method
      Set_Antipode_Time               = .true.   ! simulation time ends at reference antipode time (overrides LASTTIME )
      !-----------------------------------------------------------------------
      ! machine memory holds for 2 GB RAM:
      !   level 6: numVertices=122'882, numofTimeSteps~500,
      !                 double precision 8 byte -> needs ~ 500 MB per wavefield, still o.k.
      !   level 7: numVertices=491'522, numofTimeSteps~100,
      !                 dp 8 byte -> needs ~ 3.8 GB ! per wavefield, too big

      ! initialization of parameters and arrays
      call initialize()

      ! wait until all processes reached this point
      call MPI_Barrier( MPI_COMM_WORLD, ierror )
      if (ierror /= 0) call stopProgram('abort - MPI_Barrier kernels failed    ')

      ! prepare for simulation
      if (MAIN_PROCESS .and. VERBOSE) then
        print *
        print *,'running reference simulation...'
        print *
      endif

      ! benchmark
      if (MAIN_PROCESS) benchstart = MPI_WTIME()

      ! determine how many kernels
      if (kernelIteration) then
        ! only for homogeneous background earth
        if (HETEROGENEOUS) then
          print *,'Error: many kernels can only be calculated for a homogeneous background earth.'
          print *,'       the receivers are fixed on the equator.'
          call stopProgram( 'abort - adjointMethod')
        endif

        ! determine number of kernels
        numofKernels=int(kernelEndDistance-kernelStartDistance+1)
        if (numofKernels <= 0) then
          print *,'Error: kernels cannot be found correctly:',kernelStartDistance,kernelEndDistance
          call stopProgram( 'abort - adjointMethod')
        endif
        if (MAIN_PROCESS .and. VERBOSE) then
          print *
          print *,'iteration for number of kernels:',numofKernels
          print *,'    start/end distance :',kernelStartDistance,kernelEndDistance
          print *
        endif
        window_start_org = WINDOW_START
        window_end_org   = WINDOW_END
      else
        numofKernels = 1
      endif

      ! propagates membrane waves
      do kernel = 1,numofKernels
        if (kernelIteration) then
          ! resets time integration window
          WINDOW_START = window_start_org
          WINDOW_END = window_end_org

          ! sets new receiver position
          call prepareKernel(kernel)
        endif

        ! do the time iteration
        if (MAIN_PROCESS .and. VERBOSE) print *,'    forward simulation...'

        call forwardIteration()

        ! save seismogram at receiver
        if (.not. kernelIteration) call printSeismogram()

        ! benchmark output
        if (MAIN_PROCESS .and. VERBOSE) then
          benchend = MPI_WTIME()
          print *,'    benchmark seconds:',benchend-benchstart
          print *
        endif

        ! determine adjoint source
        call getAdjointSource()

        ! prepare for simulation
        if (MAIN_PROCESS .and. VERBOSE) then
          print *
          print *,'running adjoint simulation...'
          print *
        endif

        ! do back-in-time iteration
        call backwardIteration()

        ! wait until all processes reached this point
        !call MPI_Barrier( MPI_COMM_WORLD, ierror )
        !if (ierror /= 0) call stopProgram('abort - MPI_Barrier iterations failed    ')

        ! compute kernel
        if (.not. ADJOINT_ONTHEFLY) call frechetKernel()

        ! wait until all processes reached this point
        !call MPI_Barrier( MPI_COMM_WORLD, ierror )
        !if (ierror /= 0) call stopProgram('abort - MPI_Barrier kernels failed    ')

        ! output to kernel file
        call storeAdjointKernel()

        ! benchmark
        if (MAIN_PROCESS .and. VERBOSE) then
          benchAllEnd = MPI_WTIME()
          print *
          print *,'running time: ',int((benchAllEnd-benchAllStart)/60.0),'min ', &
                  mod((benchAllEnd-benchAllStart),60.0),'sec'
          print *
        endif

      enddo !kernel

      if (MAIN_PROCESS) print *,'done.'

      ! wait until all processes reached this point
      call MPI_Barrier( MPI_COMM_WORLD, ierror )
      if (ierror /= 0) call stopProgram('abort - final MPI_Barrier failed    ')

      ! end parallelization
      call MPI_FINALIZE(ierror)
      if (ierror /= 0) call stopProgram('abort - finalize failed    ')

      end program

