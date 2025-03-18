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
! modules needed by e.g. program phaseshift

      ! working precision & parameters
      ! these values are not set in 'Parameter_Input'
      !------------------------------------------------------------------------------------------
      module precisions
      !------------------------------------------------------------------------------------------
        ! user can modify these settings
        implicit none
        ! working precision in bytes: real = 4,double precision = 8
        integer, parameter :: WP                   = 4

        ! source parameters
        ! fixed_sourceparameter: use fixed source parameters, important when a higher grid level
        !                       is choosen but the same source should be simulated
        ! adapt_source: used to change the time parameter sigma in a empirical way to have a maximum
        !                         energy band at the considered period
        ! TimeParameterSigma: level 4: 164.6128 ; L150: 41.1795264700414 ; L75: 42.9967494824469
        !                                      TimeParameterSigma = 42.0216_WP
        !                                      good other values: q4: 204.5, q6: 51.125, 42.0216(L150)
        ! WidthParameterMu: level 4: 7.1343616E-02 ; L150: 1.783786819208144E-002 ;
        !                                  L75: 1.783786819208144E-002
        !                                  WidthParameterMu = 0.017968_WP
        !                                  good other values q4:0.0713 (20 degree),
        !                                                   q6: 0.017826 (5 degree),0.017968 (factor 8x)
        ! theta_wid: factor is taken to calculate source parameters mu
        !                     when no FIXED_SOURCEPARAMETER
        !                     (see Tape, thesis 2003, eq. 3.21)
        !                     for Love L150 - subdivision 8: fW8.0 have main period at 139.6 s ,
        !                                              subdivision 7: fW10.0 main period at 206.9 s
        !                     (check with waveparameters executable)
        ! default:
        !     sigma = 60.0
        !     mu     = 4.e-2
        ! other comparison: 180. 15.e-2 (bit better...), Carl 204.5/theta_wid=5
        real(WP) :: TimeParameterSigma              = 60.0 ! with adapt_source: L150 30.0_WP  ! L100: 30.0  , L50: 10.0 , L150: 60.0
        real(WP) :: WidthParameterMu                = 0.04 ! with adapt source: L150 0.015_WP ! L100: 2.e-2 , L50: 1.e-2, L150: 0.04
        logical,parameter :: FIXED_SOURCEPARAMETER  = .true.
        logical,parameter:: ADAPT_SOURCE            = .false.
        real(WP),parameter :: THETA_WID             = 5.0_WP ! (in degrees)

        ! filter parameters
        ! FILTERINITIALSOURCE: filter initial prescribed source too
        ! FILTERSEISMOGRAMS: use filtering when calculating timelag
        ! bw_fixfrequency: use the fixed halfband-frequency or a percentage of the wave period
        !
        ! half-bandwidth frequency: use same bandwidth-half as for analytic comparision
        !                                           (compare: Kernels/BORN/find_kernels.f, parameter deltau)
        !                                           paper suggests:
        !                                           Spetzler et al. 2002: delta_f_half = 2.5 mHz
        !                                           Ekstrom et al. 1997: delta_f_half ~ 2.27 mHz (derived from eq. 13 and fig. 2)
        !
        ! bandwidth percent: 10% of corner period used as filter width on both sides
        !                    ( e.g. Love 150s gets +- 15 s filter ), 0.6059 for half-bandwidth of 2.5 mHz
        ! butterworthfilter: just use bandpass filter or butterworth bandpass filter
        ! butterworth_power: butterworth power exponent
        ! arrival_threshold: threshold value for arrival time picking also determines startTime
        !                             for fourier transformation
        ! hannwindow_percent: hanning window width to apply at the ends of the seismogram
        !                                     in percent of the seismogram length
        ! membranecorrection: [s] to correct wave period in filtering routines;
        !                                   -10.0 for grid level 6 (empirical value; compare with single
        !                                   frequency kernel cross-sections)
        ! ANALYTICAL_CORRELATION: calculates also the analytically derived timelag kernel value
        logical,parameter :: FILTERINITIALSOURCE    = .true. !.true. for inversions/kernels
        logical,parameter :: FILTERSEISMOGRAMS      = .false.
        logical,parameter :: BW_FIXFREQUENCY        = .true.
        real(WP),parameter :: BW_HALFFREQUENCY     = 2.5e-3  ! default: 2.5e-3; Ekstrom: 2.27e-3
        real(WP),parameter :: BW_PERCENT           = 0.0_WP
        logical,parameter :: BUTTERWORTHFILTER      = .false.
        integer,parameter :: BUTTERWORTH_POWER      = 1
        real(WP),parameter :: ARRIVAL_THRESHOLD    = 1.0_WP
        real(WP),parameter :: HANNWINDOW_PERCENT   = 0.1
        real(WP),parameter :: MEMBRANECORRECTION   = 0.0
        logical,parameter :: ANALYTICAL_CORRELATION = .false.
        ! filter the source around a new period
        logical:: new_period                        = .false.
        real(WP):: newperiod                       = 150.919  ! (s)


        ! adjoint method
        ! time window of signal
        ! WINDOW_START: for major-arcs R150: 4700.;  R200: 4400.; R250: 4100.; R300: 3800.
        ! 1. orbit:         0  to   4000   for L150
        ! 2. orbit:   4000  to   8000
        ! 3. orbit:   8000  to  13500
        ! 4. orbit:  13500 to  17000
        logical,parameter:: WINDOWEDINTEGRATION      = .false.
        real(WP):: WINDOW_START                      = 0.0 ! in seconds
        real(WP):: WINDOW_END                        = 4300.0 ! in seconds


        ! routine optimization
        ! precalculated_cells: speeding up computation by using precalculated cell areas
        !                                 and distances
        ! prescribedsource: precalculates source values and stores in an array for faster
        !                              computation of finite-difference iteration
        ! relaxed_grid: takes cell infos from grid files ***.relaxed.dat
        logical,parameter::PRECALCULATED_CELLS       = .true.
        logical:: PRESCRIBEDSOURCE                   = .true.
        logical,parameter:: RELAXED_GRID             = .false.
        logical,parameter:: CORRECT_RATIO            = .false.

        ! analytical spherical harmonics degrees
        ! degree's l & m : spherical harmonic degrees for comparision with laplacian accuracy
        ! degreeAccuracy : legendre accuracy for analytic solution of seismogram with source 2
        !                              see plot fig. (3.3)B, Tape, p. 32
        !                             ( 38 for 180/0.15; 80 for 60/0.04; >120 for 10/0.02 )
        integer, parameter :: DEGREE_L              = 6
        integer, parameter :: DEGREE_M              = 1
        integer, parameter :: DEGREEACCURACY        = 120

        ! inversion
        ! compatible_refgridXrefgrid: set equal to 1 to make the grid
        !                           compatible with a refgridXrefgrid grid
        !   0 for L0035.jwkb.lsqr.mod.pcn-.2-.2 or L150.crust.2degree.gsh.modelVector.*.pcn.blk
        !   1 for L150.crust
        ! grid pixelsize: inversion grid consists of approximately equal sized
        ! pixels with this size (in  degrees)
        integer,parameter:: compatible_refgridXrefgrid  = 0

        ! source/station location
        ! station_correction: tries to use exact locations of source and receiver.
        !                               interpolates linearly for the receiver displacement.
        logical:: Station_Correction                    = .true.

        ! use overtime : for the adaption of the simulation time for different source/receiver setups
        !                         otherwise the antipode time is taken when calculating a new simulation time
        !                         (e.g. for many-kernels calculation)
        ! simulation overtime %: how much percent more time shall be given for the overtime
        logical,parameter:: USEOVERTIME                 = .false.
        real(WP),parameter:: SIMULATIONOVERTIMEPERCENT = 0.8_WP

        ! checkerboard map parameters
        ! instead of reading in a heterogeneous phase velocity map, it creates
        ! a checkerboard map with corresponding parameters
        logical,parameter :: DO_CHECKERBOARD          = .false.
        real(WP),parameter :: MAP_PERCENT_AMPLITUDE  = 2.0    !given in percent, i.e. 5% = 5.0
        integer,parameter :: MAP_DEGREE_L            = 13 ! 9 ! 20
        integer,parameter :: MAP_DEGREE_M            = 7 ! 5 ! 10

        ! simulation output
        integer,parameter :: SIMULATION_TIMESTEPPING = 4  ! takes every fifth timestep for output

        !----------------------------------------------------------------------------------------
        ! PREM values
        ! radius measured in km; phase velocities measured in km/s
        ! high periods 450/600: determined with spline representation for wave period
        ! (see: Kernels/BORN/find_kernels.f)
        real(WP), parameter :: EARTHRADIUS         = 6371.0_WP
        real(WP), parameter :: EARTHRADIUS_SQUARED = 40589641.0_WP

        ! rayleigh waves
        real(WP), parameter :: PHASEVELOCITY_R35 = 3.91253_WP
        real(WP), parameter :: PHASEVELOCITY_R37 = 3.91926_WP
        !real(WP), parameter :: PHASEVELOCITY_R40 = 3.9284_WP  ! TW 3.9284
        real(WP), parameter :: PHASEVELOCITY_R40 = 3.92799_WP  ! ETL 3.92799
        real(WP), parameter :: PHASEVELOCITY_R45 = 3.94044_WP
        real(WP), parameter :: PHASEVELOCITY_R50 = 3.95146_WP
        real(WP), parameter :: PHASEVELOCITY_R60 = 3.97339_WP
        real(WP), parameter :: PHASEVELOCITY_R75 = 4.01077_WP
        real(WP), parameter :: PHASEVELOCITY_R100 = 4.08957_WP
        real(WP), parameter :: PHASEVELOCITY_R150 = 4.30508_WP
        real(WP), parameter :: PHASEVELOCITY_R200 = 4.57676_WP
        real(WP), parameter :: PHASEVELOCITY_R250 = 4.91844_WP
        real(WP), parameter :: PHASEVELOCITY_R300 = 5.29076_WP
        real(WP), parameter :: WAVEPERIOD_R35    = 35.0987_WP
        real(WP), parameter :: WAVEPERIOD_R37    = 37.0734_WP
        !real(WP), parameter :: WAVEPERIOD_R40    = 40.1970_WP  ! TW: 40.1970
        real(WP), parameter :: WAVEPERIOD_R40    = 40.0432_WP  ! ETL: 40.0432
        real(WP), parameter :: WAVEPERIOD_R45    = 45.0501_WP
        real(WP), parameter :: WAVEPERIOD_R50    = 50.0271_WP
        real(WP), parameter :: WAVEPERIOD_R60    = 60.1467_WP
        real(WP), parameter :: WAVEPERIOD_R75    = 75.3258_WP
        real(WP), parameter :: WAVEPERIOD_R100    = 100.3930_WP
        real(WP), parameter :: WAVEPERIOD_R150    = 151.1930_WP ! minor arc?
        !real(WP), parameter :: WAVEPERIOD_R150    = 150.0_WP ! major arc?
        real(WP), parameter :: WAVEPERIOD_R200    = 200.0_WP
        real(WP), parameter :: WAVEPERIOD_R250    = 250.0_WP
        real(WP), parameter :: WAVEPERIOD_R300    = 300.0_WP

        ! love waves
        real(WP), parameter :: PHASEVELOCITY_L35 = 4.39987_WP
        real(WP), parameter :: PHASEVELOCITY_L37 = 4.41795_WP
        real(WP), parameter :: PHASEVELOCITY_L40 = 4.44109_WP
        real(WP), parameter :: PHASEVELOCITY_L45 = 4.47005_WP
        real(WP), parameter :: PHASEVELOCITY_L50 = 4.49336_WP
        real(WP), parameter :: PHASEVELOCITY_L60 = 4.53031_WP
        real(WP), parameter :: PHASEVELOCITY_L75 = 4.57474_WP
        real(WP), parameter :: PHASEVELOCITY_L100= 4.64431_WP
        ! Trampert & Woodhouse 1995/1996
        real(WP), parameter :: PHASEVELOCITY_L150= 4.78619_WP ! TW: 4.78619
        ! Ekstrom et al. 1997
        !real(WP), parameter :: PHASEVELOCITY_L150= 4.77915_WP ! ETL: 4.77915 first arrival
        real(WP), parameter :: PHASEVELOCITY_L200= 4.91928_WP
        real(WP), parameter :: PHASEVELOCITY_L250= 5.07097_WP
        real(WP), parameter :: PHASEVELOCITY_L300= 5.22906_WP
        real(WP), parameter :: PHASEVELOCITY_L450= 5.74979_WP
        real(WP), parameter :: PHASEVELOCITY_L600= 6.36841_WP
        real(WP), parameter :: WAVEPERIOD_L35    = 35.0599_WP
        real(WP), parameter :: WAVEPERIOD_L37    = 37.0585_WP
        real(WP), parameter :: WAVEPERIOD_L40    = 40.1497_WP
        real(WP), parameter :: WAVEPERIOD_L45    = 45.1143_WP
        real(WP), parameter :: WAVEPERIOD_L50    = 50.1901_WP
        real(WP), parameter :: WAVEPERIOD_L60    = 60.3145_WP
        real(WP), parameter :: WAVEPERIOD_L75    = 75.1095_WP
        real(WP), parameter :: WAVEPERIOD_L100   = 100.8090_WP
        ! Trampert & Woodhouse 1995/1996
        real(WP), parameter :: WAVEPERIOD_L150   = 153.462_WP ! TW: 153.462
        ! Ekstrom et al. 1997
        !real(WP), parameter :: WAVEPERIOD_L150   = 150.9190_WP ! ETL: 150.9190
        real(WP), parameter :: WAVEPERIOD_L200    = 200.0_WP
        real(WP), parameter :: WAVEPERIOD_L250    = 250.0_WP
        real(WP), parameter :: WAVEPERIOD_L300    = 300.0_WP

        ! useful parameters
        real(WP), parameter :: PI                 = 3.1415926535897931_WP
        real(WP), parameter :: RAD2DEGREE         = 180.0_WP/PI
        real(WP), parameter :: DEGREE2RAD         = PI/180.0_WP

        ! file I/O
        integer, parameter :: IIN = 40
        integer, parameter :: IOUT = 41
      end module precisions

!-----------------------------------------------------------------------
      module phaseVelocityMap
!-----------------------------------------------------------------------
      ! phaseVelocityMap module
      ! used for phase velocity maps
        use precisions, only: WP,IIN,IOUT, &
          PRECALCULATED_CELLS, &
          WINDOW_START,WINDOW_END, &
          MAP_DEGREE_L,MAP_DEGREE_M,MAP_PERCENT_AMPLITUDE
        implicit none
        integer:: numPhaseEntries
        real(WP), allocatable, dimension(:):: phaseMap
        real(WP), allocatable, dimension(:):: phaseVelocitySquare
      end module

!-----------------------------------------------------------------------
      module verbosity
!-----------------------------------------------------------------------
        implicit none
        logical:: VERBOSE     = .true.         !for console verbosity
        logical:: beVerbose   = .false.         !for getTimelag verbosity
        logical:: fileOutput  = .false.        !for debuging: outputs to file 'tmp*.dat'
      end module

!-----------------------------------------------------------------------
      module adjointVariables
!-----------------------------------------------------------------------
      ! adjoint module
      ! used for adjoint routines
        use precisions, only: WP,PRECALCULATED_CELLS,DO_CHECKERBOARD,MAP_DEGREE_L,MAP_DEGREE_M
        implicit none
        ! adjoint parameters
        ! onthefly: instead of integration at the end, after each time step
        ! startatzero: start integration from zero, not for whole seismogram
        ! precalculate derivatives: calculate derivatives after forward iteration has finished,
        !                                         instead of just before integration
        ! windowed integration: do integration between zero and arrivaltime for seismogram
        logical,parameter:: ADJOINT_ONTHEFLY            = .false.
        logical,parameter:: ADJOINT_STARTATZERO         = .false.
        logical,parameter:: PRECALCULATE_DERIVATIVES    = .false.
        logical,parameter:: BYSPLINE                    = .false.

        ! adjoint variables
        logical:: Adjoint_Program                        = .false.
        logical:: Adjoint_InversionProgram               = .false.
        logical:: kernelIteration,storeAsFile,Set_Antipode_Time
        integer:: adjointSourceVertex
        character(len=64):: adjointKernelName
        real(WP):: normFactor
        real(WP),allocatable,dimension(:,:):: adjointSource
        real(WP),allocatable,dimension(:):: adjointKernel
        real(WP),allocatable,dimension(:):: backwardDisplacement, &
                                           backwardDisplacement_old, &
                                           backwardNewdisplacement
        real(WP),allocatable,dimension(:,:):: wavefieldForward, wavefieldAdjoint, &
                                            seismoSecondDerivative
        integer,parameter:: adjSourceFileID = 101,adjRecFileID = 102,adjMidpointFileID = 103
      end module

!-----------------------------------------------------------------------
      module cells
!-----------------------------------------------------------------------
      ! cells module
      ! used for precalculated values of cell areas etc.
        use precisions, only: WP,EARTHRADIUS,EARTHRADIUS_SQUARED,PI,DEGREE2RAD,IIN, &
          CORRECT_RATIO,RELAXED_GRID,Station_Correction
        implicit none
        integer:: subdivisions,MaxTriangles,MaxVertices
        integer,allocatable, dimension(:,:):: cellFace,cellNeighbors,cellTriangleFace
        integer:: numFaces,numNeighbors,numTriangleFaces,numCorners,numVertices, &
                 numDomainVertices
        real(WP), allocatable, dimension(:):: cellAreas
        real(WP), allocatable, dimension(:,:):: vertices,cellEdgesLength, &
                                                cellCenterDistances,cellCorners,cellFractions
        real(WP):: interpolation_distances(3),interpolation_triangleLengths(3)
        integer:: interpolation_corners(3),interpolation_triangleIndex
      end module cells

!-----------------------------------------------------------------------
      module propagationStartup
!-----------------------------------------------------------------------
      ! propagationStartup module
      ! used for initialization startup process
        use precisions, only: WP,IOUT,PI,EARTHRADIUS, &
          SIMULATIONOVERTIMEPERCENT,USEOVERTIME,WINDOWEDINTEGRATION, &
          FILTERINITIALSOURCE,PRESCRIBEDSOURCE, &
          FILTERSEISMOGRAMS, &
          Station_Correction,timeParameterSigma
        implicit none
        logical:: Phaseshift_Program = .false.
        logical:: HetPhaseshift_Program = .false.
        logical:: HETEROGENEOUS,DELTA,SECONDDELTA,MOVEDELTA
        logical:: SIMULATIONOUTPUT
        logical:: manyReceivers,manyKernels
        logical:: importKernelsReceivers,referenceRun
        integer:: startVertex,endVertex,firsttimestep,lasttimestep, &
                 numofTimeSteps,midpointVertex
        integer:: sourceVertex,receiverVertex,deltaVertex,originSourceVertex, &
                  originReceiverVertex
        integer:: numofReceivers,numofKernels,currentKernel
        integer,allocatable,dimension(:):: receivers
        integer,allocatable,dimension(:,:):: kernelsReceivers
        character(len=8):: DELTAfunction,cphasetype
        character(len=128):: datadirectory
        real(WP):: FIRSTTIME,LASTTIME,DELTARADIUS
        real(WP):: sourceLat,sourceLon,receiverLat,receiverLon,cphaseRef, &
                   desiredSourceLat,desiredSourceLon,muSquare,muTwo
        real(WP):: deltaLat,deltaLon,deltaPerturbation,deltaMoveIncrement, &
                   desiredReceiverLat,desiredReceiverLon
        real(WP):: cphase2,dt,dt2,averageCellDistance
        real(WP):: kernelStartDistance,kernelEndDistance,distanceQuit
        real(WP):: benchAllStart,benchAllEnd,benchstart,benchend
        real(WP),allocatable,dimension(:,:):: seismogramReceiver, &
                                             receiversSeismogram, &
                                             receiversSeismogramRef
        real(WP),allocatable,dimension(:,:,:)::kernelsReceiversSeismogram, &
                                              kernelsReceiversSeismogramRef
        real(WP),allocatable,dimension(:,:)::forceTermPrescribed
        logical:: sourceOnFile = .false.
        integer:: sourceFileID = 201
        integer:: SH_lmx,SH_ncoef
        real,allocatable,dimension(:):: SH_coef
        ! rotate source&receiver locations /
        ! heterogeneous phase map to equatorial plane (by default)
        logical:: rotate_frame = .false.
      end module

!-----------------------------------------------------------------------
      module deltaSecondLocation
!-----------------------------------------------------------------------
      ! additional delta location
        use precisions, only: WP
        implicit none
        integer:: deltaSecondVertex
        real(WP):: deltaSecondLat,deltaSecondLon
        real(WP),parameter:: deltaSecondDistance = 3.0_WP
      end module

!-----------------------------------------------------------------------
      module parallel
!-----------------------------------------------------------------------
      ! parallel module
      ! used for MPI parallelization
        use mpi
        implicit none
        logical:: PARALLELSEISMO
        logical:: MAIN_PROCESS = .false.
        integer:: nprocesses,rank,tag,MPI_CUSTOM
        integer:: status(MPI_STATUS_SIZE)
      end module

!-----------------------------------------------------------------------
      module displacements
!-----------------------------------------------------------------------
      ! displacement arrays
        use precisions, only: WP
        implicit none
        real(WP),allocatable,dimension(:):: displacement
        real(WP),allocatable,dimension(:):: displacement_old
        real(WP),allocatable,dimension(:):: newdisplacement
      end module

!-----------------------------------------------------------------------
      module griddomain
!-----------------------------------------------------------------------
      ! grid domain array
        use precisions, only: WP,PI
        implicit none
        integer:: boundariesMaxRange
        integer,allocatable,dimension(:)::domainVertices
        integer,allocatable,dimension(:)::vertexDomain
        integer,allocatable,dimension(:,:,:)::boundaries
        integer,allocatable,dimension(:,:)::domainNeighbors
        real(WP),allocatable,dimension(:)::sendDisp,receiveDisp
      end

!-----------------------------------------------------------------------
      module loop
!-----------------------------------------------------------------------
      ! looping delta location
        implicit none
        integer:: latitudeStart, latitudeEnd
        integer:: longitudeEnd
      end

!-----------------------------------------------------------------------
      module phaseBlockData
!-----------------------------------------------------------------------
      ! phase module
        use precisions, only: WP,IIN,compatible_refgridXrefgrid
        implicit none
        integer:: numBlocks
        real(WP), allocatable, dimension(:):: phaseBlock
        real(WP):: phaseBlockVelocityReference
        character(len=128):: phaseBlockFile,heterogeneousDataFile,heterogeneousOutput
        real:: heterogeneousPixelsize
        integer:: gsh_maximum_expansion
      end module

!-----------------------------------------------------------------------
      module traveltime
!-----------------------------------------------------------------------
      ! traveltime common parameters
        use precisions, only: WP
        implicit none
        real(WP):: t_lag,arrivalTime,vertexCellArea,t_lagAnalytic,seismo_timestep
        real(WP):: phasevelocity,kernel,kernelAnalytic
      end module

!-----------------------------------------------------------------------
      module filterType
!-----------------------------------------------------------------------
        use precisions
        implicit none
        ! window size for fft
        integer:: WindowSIZE
        ! full bandwidth (in frequency domain)
        integer:: bw_width
        ! half-bandwidth frequency, wave period used for bandwidth filtering
        real(WP):: bw_frequency,bw_waveperiod
      end module

!-----------------------------------------------------------------------
      module splinefunction
!-----------------------------------------------------------------------
        implicit none
        integer:: i1,i2
        double precision, allocatable,dimension(:):: X,Y    ! dimension (i2)
        double precision, allocatable,dimension(:,:):: Q,F  ! dimension (3,i2)
      end module

!-----------------------------------------------------------------------
      module heterogeneousMatrix
!-----------------------------------------------------------------------
        implicit none
        real:: Amin,Amax
      end

