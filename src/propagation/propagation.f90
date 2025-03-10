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
      program propagation
!-----------------------------------------------------------------------
! iterates through one time step with numerically
! calculation of the new displacement
!
! uses: Carl Tape M.Sc. thesis, 2003, chap 5, (5.7)
      use propagationStartup; use cells; use phaseVelocityMap; use parallel
      use displacements; use griddomain; use loop; use phaseBlockData; use verbosity
      implicit none
      double precision:: u_t, u_tplus1, u_tminus1,forcing,initialShape,forcingRef
      double precision:: D2,time,discreteLaplacian, precalc_discreteLaplacian
      double precision:: forceTerm2Source, initialShapeTerm
      double precision:: area,cphaseDefault,vectorV(3),lat,lon
      integer::         n,k,i,timestep,vertex,index,ierror
      character(len=1)::     rankstr, dstr
      character(len=4)::     timestr
      logical:: looping
      external:: initialShapeTerm,discreteLaplacian,precalc_discreteLaplacian

      !-----------------------------------------------------------------------
      ! parameters
      ! most parameters concerning wave propagation set in file Parameter_Input
      ! (& default values in commonModules.f90)
      !-----------------------------------------------------------------------

      !print *,'Welcome to membrane waves on a spherical shell'
      !print *,'----------------------------------------------'
      !print *

      ! initialization of parameters and arrays
      call initialize()

      ! place all receiver stations
      if ( manyReceivers ) then
        ! useless for heterogeneous case
        if ( HETEROGENEOUS ) then
          print *,'multiple receiver stations and heterogeneous phase map not applicable!'
          stop
        endif

        if ( numofReceivers <= 0 ) then
          print *,'number of receiver stations invalid!',numofReceivers
          stop
        endif

        ! allocate seismograms for all receivers
        deallocate( seismogramReceiver )    ! only needed for single receiver setups
        allocate( receivers(numofReceivers), receiversSeismogram(numofReceivers+1,lasttimestep-firsttimestep+1),stat=ierror)
        if ( ierror > 0 ) then
          print *,'cannot allocate receivers array'
          stop
        endif

        ! assumes the source to be on north pole, the receiver station places along the same latitude
        do i = 1, numofReceivers
          ! vary longitude
          receiverLon=(i-1)*360.0/numofReceivers
          call findVertex(receiverLat,receiverLon,receivers(i))
        enddo
      endif

      if ( MASTER .and. VERBOSE ) print *,'starting time loop...'

      ! loop for each delta location
      looping = .true.
      do while( looping )
        ! benchmark
        if ( MASTER ) benchstart = MPI_WTIME()

        ! forward simulation of membrane waves
        call forwardIteration()

        ! print seismogram to file
        call printSeismogram()

        ! move delta location
        if ( DELTA ) then
          if ( PARALLELSEISMO ) then
            call parallelFindnextLocation(looping)
          else
            call findnextLocation(looping)
          endif
        endif

        ! stop looping
        if (.not. MOVEDELTA .or. .not. DELTA ) looping = .false.

        ! precalculate the phase velocities for all grid points (delta position may have changed)
        if (looping) call constructPhaseVelocitySquare()

        ! benchmark output
        if ( MASTER .and. VERBOSE) then
          benchend = MPI_WTIME()
          print *
          print *,'benchmark seconds:',benchend-benchstart
          print *
        endif
      enddo !delta location

      ! wait until all processes reached this point
      call MPI_Barrier( MPI_COMM_WORLD, ierror )
      if ( ierror /= 0) call stopProgram('abort - final MPI_Barrier failed    ')

      ! end parallelization
      call MPI_FINALIZE(ierror)
      if ( ierror /= 0) call stopProgram('abort - finalize failed    ')

      end
