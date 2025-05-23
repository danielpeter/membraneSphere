!=====================================================================
!
!       m e m b r a n e S p h e r e
!       --------------------------------------------------
!      (c) 2025
!
!=====================================================================
!
! modules needed by e.g. program phaseshift

! working precision & parameters

!------------------------------------------------------------------------------------------
module constants
!------------------------------------------------------------------------------------------
  implicit none
  ! defines working precision (WP) in bytes: real = 4,double precision = 8
  include "constants.h"
end module constants

! these values are not set in 'Parameter_Input'
!------------------------------------------------------------------------------------------
module precisions
!------------------------------------------------------------------------------------------
  ! user can modify these settings
  use constants
  implicit none

  ! source parameters
  ! FIXED_SOURCEPARAMETER  - use fixed source parameters, important when a higher grid level
  !                          is choosen but the same source should be simulated
  ! ADAPT_SOURCE           - used to change the time parameter sigma in a empirical way to have a maximum
  !                          energy band at the considered period
  ! THETA_WIDTH            - factor is taken to calculate source parameters mu when no FIXED_SOURCEPARAMETER
  !                          (see Tape, thesis 2003, eq. 3.21)
  !                          for Love L150 - subdivision 8: fW8.0 have main period at 139.6 s ,
  !                                          subdivision 7: fW10.0 main period at 206.9 s
  !                          (check with waveparameters executable)
  logical,parameter :: FIXED_SOURCEPARAMETER  = .true.
  logical,parameter :: ADAPT_SOURCE           = .false.
  real(WP),parameter :: THETA_WIDTH           = 5.0_WP ! (in degrees)

  ! filter parameters
  ! FILTER_SEISMOGRAMS     - use filtering when calculating timelag
  !
  ! BW_FIXFREQUENCY        - bandwidth filter use the fixed halfband-frequency or a percentage of the wave period
  ! BW_HALFFREQUENCY       - use same bandwidth-half as for analytic comparision
  !                          (compare: Kernels/BORN/find_kernels.f, parameter deltau)
  !                          paper suggests:
  !                            Spetzler et al. 2002: delta_f_half = 2.5 mHz
  !                            Ekstrom et al. 1997 : delta_f_half ~ 2.27 mHz (derived from eq. 13 and fig. 2)
  ! BW_PERCENT             - 10% of corner period used as filter width on both sides
  !                          (e.g., Love 150s gets +- 15 s filter ), 0.6059 for half-bandwidth of 2.5 mHz
  !
  ! BUTTERWORTHFILTER      - use butterworth bandpass filter or just bandpass filter
  ! BUTTERWORTH_POWER      - butterworth power exponent
  !
  ! ARRIVAL_THRESHOLD      - threshold value for arrival time picking also determines startTime
  !                          for fourier transformation
  ! HANNWINDOW_PERCENT     - Hanning window width to apply at the ends of the seismogram
  !                          in percent of the seismogram length
  ! MEMBRANECORRECTION     - in [s] to correct wave period in filtering routines
  !                          for grid level 6 == -10.0 (empirical value; compare with single frequency kernel cross-sections)
  logical,parameter :: FILTER_SEISMOGRAMS     = .false.

  logical,parameter :: BW_FIXFREQUENCY        = .true.
  real(WP),parameter :: BW_HALFFREQUENCY      = 2.5e-3  ! default: 2.5e-3; Ekstrom: 2.27e-3
  real(WP),parameter :: BW_PERCENT            = 0.0_WP

  logical,parameter :: BUTTERWORTHFILTER      = .false.
  integer,parameter :: BUTTERWORTH_POWER      = 1

  real(WP),parameter :: ARRIVAL_THRESHOLD     = 1.0_WP
  real(WP),parameter :: HANNWINDOW_PERCENT    = 0.1_WP
  real(WP),parameter :: MEMBRANECORRECTION    = 0.0_WP

  ! calculates also the analytically derived timelag kernel value
  logical,parameter :: ANALYTICAL_CORRELATION = .false.

  ! filter the source around a new period (needs FILTER_INITIALSOURCE turned on)
  logical,parameter :: FILTER_AT_NEWPERIOD    = .false.
  real(WP),parameter :: BW_NEWPERIOD          = 150.919  ! (s)

  ! routine optimization
  ! PRECALCULATED_CELLS    - speeding up computation by using precalculated cell areas and distances
  ! PRESCRIBED_SOURCE      - precalculates source values and stores in an array for faster
  !                          computation of finite-difference iteration
  ! RELAXED_GRID           - takes cell infos from grid files ***.relaxed.dat
  logical,parameter :: PRECALCULATED_CELLS     = .true.
  logical,parameter :: PRESCRIBED_SOURCE       = .true.
  logical,parameter :: RELAXED_GRID            = .false.
  logical,parameter :: CORRECT_RATIO           = .false.

  ! analytical spherical harmonics degrees
  ! DEGREE_L & DEGREE_M    - spherical harmonic degrees for comparision with laplacian accuracy
  ! DEGREEACCURACY         - legendre accuracy for analytic solution of seismogram with source 2
  !                           see plot fig. (3.3)B, Tape, p. 32
  !                           (38 for 180/0.15; 80 for 60/0.04; >120 for 10/0.02)
  integer, parameter :: DEGREE_L              = 6
  integer, parameter :: DEGREE_M              = 1
  integer, parameter :: DEGREEACCURACY        = 120

  ! inversion
  ! compatible_refgridXrefgrid - set equal to 1 to make the grid compatible with a refgridXrefgrid grid
  !                              0 for L0035.jwkb.lsqr.mod.pcn-.2-.2 or L150.crust.2degree.gsh.modelVector.*.pcn.blk
  !                              1 for L150.crust
  ! grid pixelsize: inversion grid consists of approximately equal sized
  ! pixels with this size (in  degrees)
  integer,parameter :: compatible_refgridXrefgrid = 0

  ! source/station location
  ! STATION_CORRECTION     - tries to use exact locations of source and receiver.
  !                          interpolates linearly for the receiver displacement.
  logical,parameter :: STATION_CORRECTION         = .true.

  ! USE_OVERTIME           - for the adaption of the simulation time for different source/receiver setups
  !                          otherwise the antipode time is taken when calculating a new simulation time
  !                          (e.g. for many-kernels calculation)
  ! OVERTIME_PERCENT       - how much percent more time shall be given for the overtime
  logical,parameter :: USE_OVERTIME               = .false.
  real(WP),parameter :: OVERTIME_PERCENT          = 0.8_WP

  ! checkerboard map parameters
  ! instead of reading in a heterogeneous phase velocity map, it creates
  ! a checkerboard map with corresponding parameters
  logical,parameter :: DO_CHECKERBOARD           = .false.
  integer,parameter :: MAP_DEGREE_L              = 13       ! or: 9, 13, 20
  integer,parameter :: MAP_DEGREE_M              = 7        ! or: 5, 7, 10
  real(WP),parameter :: MAP_PERCENT_AMPLITUDE    = 2.0      ! given in percent, i.e. 5% = 5.0

  ! simulation output
  integer,parameter :: SIMULATION_TIMESTEPPING   = 1        ! saves every N timestep a file output (4 == every fourth)
  real(WP),parameter :: SIMULATION_STARTTIME     = -1000.0  ! starts outputting files after this time
  logical,parameter :: SIMULATION_OUTPUT_VTK     = .true.   ! outputs as vtk file, otherwise as .dat file

  ! rotate source&receiver locations / heterogeneous phase map to equatorial plane (by default)
  logical,parameter :: ROTATE_FRAME              = .false.

  !----------------------------------------------------------------------------------------
  ! Surface waves
  ! phase velocities measured in km/s
  ! high periods 450/600: determined with spline representation for wave period
  ! (see: Kernels/BORN/find_kernels.f)

  ! Rayleigh waves
  ! reference phase velocities
  real(WP), parameter :: PHASEVELOCITY_R35  = 3.91253_WP
  real(WP), parameter :: PHASEVELOCITY_R37  = 3.91926_WP
  !real(WP), parameter :: PHASEVELOCITY_R40  = 3.9284_WP  ! TW 3.9284
  real(WP), parameter :: PHASEVELOCITY_R40  = 3.92799_WP  ! ETL 3.92799
  real(WP), parameter :: PHASEVELOCITY_R45  = 3.94044_WP
  real(WP), parameter :: PHASEVELOCITY_R50  = 3.95146_WP
  real(WP), parameter :: PHASEVELOCITY_R60  = 3.97339_WP
  real(WP), parameter :: PHASEVELOCITY_R75  = 4.01077_WP
  real(WP), parameter :: PHASEVELOCITY_R80  = 4.01860_WP  ! TW
  real(WP), parameter :: PHASEVELOCITY_R100 = 4.08957_WP
  real(WP), parameter :: PHASEVELOCITY_R150 = 4.30508_WP
  real(WP), parameter :: PHASEVELOCITY_R200 = 4.57676_WP
  real(WP), parameter :: PHASEVELOCITY_R250 = 4.91844_WP
  real(WP), parameter :: PHASEVELOCITY_R300 = 5.29076_WP
  ! reference periods
  real(WP), parameter :: WAVEPERIOD_R35     = 35.0987_WP
  real(WP), parameter :: WAVEPERIOD_R37     = 37.0734_WP
  !real(WP), parameter :: WAVEPERIOD_R40     = 40.1970_WP  ! TW: 40.1970
  real(WP), parameter :: WAVEPERIOD_R40     = 40.0432_WP  ! ETL: 40.0432
  real(WP), parameter :: WAVEPERIOD_R45     = 45.0501_WP
  real(WP), parameter :: WAVEPERIOD_R50     = 50.0271_WP
  real(WP), parameter :: WAVEPERIOD_R60     = 60.1467_WP
  real(WP), parameter :: WAVEPERIOD_R75     = 75.3258_WP
  real(WP), parameter :: WAVEPERIOD_R80     = 78.1290_WP  ! TW
  real(WP), parameter :: WAVEPERIOD_R100    = 100.3930_WP
  real(WP), parameter :: WAVEPERIOD_R150    = 151.1930_WP ! minor arc?
  !real(WP), parameter :: WAVEPERIOD_R150    = 150.0_WP ! major arc?
  real(WP), parameter :: WAVEPERIOD_R200    = 200.0_WP
  real(WP), parameter :: WAVEPERIOD_R250    = 250.0_WP
  real(WP), parameter :: WAVEPERIOD_R300    = 300.0_WP

  ! Love waves
  ! reference phase velocities
  real(WP), parameter :: PHASEVELOCITY_L35  = 4.39987_WP
  real(WP), parameter :: PHASEVELOCITY_L37  = 4.41795_WP
  real(WP), parameter :: PHASEVELOCITY_L40  = 4.44109_WP
  real(WP), parameter :: PHASEVELOCITY_L45  = 4.47005_WP
  real(WP), parameter :: PHASEVELOCITY_L50  = 4.49336_WP
  real(WP), parameter :: PHASEVELOCITY_L60  = 4.53031_WP
  real(WP), parameter :: PHASEVELOCITY_L75  = 4.57474_WP
  real(WP), parameter :: PHASEVELOCITY_L80  = 4.58561_WP
  real(WP), parameter :: PHASEVELOCITY_L100 = 4.64431_WP
  ! Trampert & Woodhouse 1995/1996
  real(WP), parameter :: PHASEVELOCITY_L150 = 4.78619_WP ! TW: 4.78619
  ! Ekstrom et al. 1997
  !real(WP), parameter :: PHASEVELOCITY_L150 = 4.77915_WP ! ETL: 4.77915 first arrival
  real(WP), parameter :: PHASEVELOCITY_L200 = 4.91928_WP
  real(WP), parameter :: PHASEVELOCITY_L250 = 5.07097_WP
  real(WP), parameter :: PHASEVELOCITY_L300 = 5.22906_WP
  real(WP), parameter :: PHASEVELOCITY_L450 = 5.74979_WP
  real(WP), parameter :: PHASEVELOCITY_L600 = 6.36841_WP
  ! reference periods
  real(WP), parameter :: WAVEPERIOD_L35     = 35.0599_WP
  real(WP), parameter :: WAVEPERIOD_L37     = 37.0585_WP
  real(WP), parameter :: WAVEPERIOD_L40     = 40.1497_WP
  real(WP), parameter :: WAVEPERIOD_L45     = 45.1143_WP
  real(WP), parameter :: WAVEPERIOD_L50     = 50.1901_WP
  real(WP), parameter :: WAVEPERIOD_L60     = 60.3145_WP
  real(WP), parameter :: WAVEPERIOD_L75     = 75.1095_WP
  real(WP), parameter :: WAVEPERIOD_L80     = 79.0020_WP ! TW
  real(WP), parameter :: WAVEPERIOD_L100    = 100.8090_WP
  ! Trampert & Woodhouse 1995/1996
  real(WP), parameter :: WAVEPERIOD_L150    = 153.462_WP ! TW: 153.462
  ! Ekstrom et al. 1997
  !real(WP), parameter :: WAVEPERIOD_L150    = 150.9190_WP ! ETL: 150.9190
  real(WP), parameter :: WAVEPERIOD_L200    = 200.0_WP
  real(WP), parameter :: WAVEPERIOD_L250    = 250.0_WP
  real(WP), parameter :: WAVEPERIOD_L300    = 300.0_WP
end module precisions

!-----------------------------------------------------------------------
module parameter_inputs
!-----------------------------------------------------------------------
  use constants, only: WP
  implicit none
  ! Parameter_file inputs
  ! Physical model
  integer:: subdivisions           = 0        ! LEVEL
  real(WP):: FIRSTTIME             = 0.0
  real(WP):: LASTTIME              = 0.0
  character(len=8):: cphasetype    = ""
  real(WP):: cphaseRef

  ! Geometry
  real(WP):: sourceLat,sourceLon
  real(WP):: receiverLat,receiverLon
  logical:: manyReceivers          = .false.
  integer:: numofReceivers         = 0
  logical:: manyKernels            = .false.
  real(WP):: kernelStartDistance,kernelEndDistance
  logical:: importKernelsReceivers = .false.

  ! Source
  ! source default:
  !     sigma = 60.0
  !     mu    = 4.e-2
  ! other comparison: 180. 15.e-2 (bit better...), Carl 204.5/theta_width=5
  real(WP):: TimeParameterSigma    = 60.0 ! with adapt_source: L150 30.0_WP; L100: 30.0; L50: 10.0
  real(WP):: WidthParameterMu      = 0.04 ! with adapt source: L150 0.015_WP; L100: 0.02; L50: 0.01
  ! FILTER_INITIALSOURCE   - filter initial prescribed source
  !                          (applied only when PRESCRIBED_SOURCE is turned on; useful when calculating kernels
  !                           to filter out spurious frequency artefacts)
  logical:: FILTER_INITIALSOURCE   = .false.  !.true. for inversions/kernels
  ! adjoint source
  ! ADJOINT_TAPERSIGNAL       - apply tapering to adjoint source signal
  logical:: ADJOINT_TAPERSIGNAL    = .true.
  ! WINDOWED_INTEGRATION       - do integration between zero and arrivaltime for seismogram
  !
  ! time window of signal
  ! WINDOW_START: for major-arcs R150: 4700.;  R200: 4400.; R250: 4100.; R300: 3800.
  ! 1. orbit:      0  to   4000   for L150
  ! 2. orbit:   4000  to   8000
  ! 3. orbit:   8000  to  13500
  ! 4. orbit:  13500  to  17000
  logical:: WINDOWED_INTEGRATION   = .false.
  real(WP):: WINDOW_START          = 0.0_WP    ! in seconds
  real(WP):: WINDOW_END            = 4300.0_WP ! in seconds
  ! ADJOINT_ANTIPODE_TIME     - simulation time ends at reference antipode time (overrides LASTTIME )
  logical:: ADJOINT_ANTIPODE_TIME  = .false.

  ! Scatterer
  logical:: DELTA                  = .false.
  real(WP):: DELTARADIUS           = 0.0
  character(len=8):: DELTAfunction = ""
  real(WP):: deltaPerturbation
  real(WP):: deltaLat,deltaLon
  logical:: MOVEDELTA              = .false.
  integer:: latitudeStart, latitudeEnd
  integer:: longitudeEnd
  real(WP):: deltaMoveIncrement
  logical:: SECONDDELTA            = .false.

  ! heterogeneous model
  logical:: HETEROGENEOUS          = .false.
  character(len=128):: phaseBlockFile
  real:: heterogeneousPixelsize
  integer:: gsh_maximum_expansion
  real(WP):: phaseBlockVelocityReference
  character(len=128):: heterogeneousDataFile,heterogeneousOutput

  ! file output & misc
  character(len=128):: datadirectory    = ""
  character(len=64):: adjointKernelName = ""
  logical:: VERBOSE                = .true.          !for console verbosity
  logical:: SIMULATIONOUTPUT       = .false.
  logical:: PARALLELSEISMO         = .false.
end module

!-----------------------------------------------------------------------
module phaseVelocityMap
!-----------------------------------------------------------------------
! phaseVelocityMap module
! used for phase velocity maps
  use constants, only: WP,IIN,IOUT
  use precisions, only: &
    PRECALCULATED_CELLS, &
    MAP_DEGREE_L,MAP_DEGREE_M,MAP_PERCENT_AMPLITUDE
  implicit none
  integer:: numPhaseEntries
  real(WP), allocatable, dimension(:):: phaseMap
  real(WP), allocatable, dimension(:):: phaseVelocitySquare
end module

!-----------------------------------------------------------------------
module verbosity
!-----------------------------------------------------------------------
  use parameter_inputs, only: VERBOSE
  implicit none
  logical:: beVerbose   = .false.         !for getTimelag verbosity
  logical:: fileOutput  = .false.         !for debuging: outputs to files, e.g. 'seismo.displacement.dat',..
end module

!-----------------------------------------------------------------------
module adjointVariables
!-----------------------------------------------------------------------
! adjoint module
! used for adjoint routines
  use constants, only: WP
  use precisions, only: PRECALCULATED_CELLS,DO_CHECKERBOARD,MAP_DEGREE_L,MAP_DEGREE_M
  use parameter_inputs, only: adjointKernelName
  implicit none
  ! adjoint parameters
  ! ADJOINT_ONTHEFLY          - instead of integration at the end, after each time step
  ! ADJOINT_STARTATZERO       - start integration from zero, not for whole seismogram
  ! PRECALCULATE_DERIVATIVES  - calculate derivatives after forward iteration has finished,
  !                             instead of just before integration
  ! ADJOINT_INTEGRALBYSPLINE  - calculates kernel integral by spline representation
  logical,parameter:: ADJOINT_ONTHEFLY            = .false.
  logical,parameter:: ADJOINT_STARTATZERO         = .false.
  logical,parameter:: PRECALCULATE_DERIVATIVES    = .false.
  logical,parameter:: ADJOINT_INTEGRALBYSPLINE    = .false.

  ! adjoint variables
  logical:: Adjoint_Program                       = .false.
  logical:: Adjoint_InversionProgram              = .false.
  logical:: kernelIteration,storeAsFile
  integer:: adjointSourceVertex
  real(WP),allocatable,dimension(:,:):: adjointSource
  real(WP),allocatable,dimension(:):: adjointKernel
  real(WP),allocatable,dimension(:):: backwardDisplacement, backwardDisplacement_old, backwardNewdisplacement
  real(WP),allocatable,dimension(:,:):: wavefieldForward, wavefieldAdjoint, seismoSecondDerivative
  integer,parameter:: adjSourceFileID = 301,adjRecFileID = 302,adjMidpointFileID = 303
end module

!-----------------------------------------------------------------------
module cells
!-----------------------------------------------------------------------
! cells module
! used for precalculated values of cell areas etc.
  use constants, only: WP,EARTHRADIUS,EARTHRADIUS_SQUARED,PI,DEGREE2RAD,IIN
  use precisions, only: CORRECT_RATIO,RELAXED_GRID,STATION_CORRECTION
  use parameter_inputs, only: subdivisions
  implicit none
  integer:: MaxTriangles,MaxVertices
  integer,allocatable, dimension(:,:):: cellFace,cellNeighbors,cellTriangleFace
  integer:: numFaces,numNeighbors,numTriangleFaces,numCorners,numVertices,numDomainVertices

  real(WP), allocatable, dimension(:):: cellAreas
  real(WP), allocatable, dimension(:,:):: vertices,cellEdgesLength, &
                                          cellCenterDistances,cellCorners,cellFractions
  real(WP):: interpolation_distances(3),interpolation_triangleLengths(3)
  integer:: interpolation_corners(3),interpolation_triangleIndex
  real(WP):: averageCellDistance
end module cells

!-----------------------------------------------------------------------
module propagationStartup
!-----------------------------------------------------------------------
! propagationStartup module
! used for initialization startup process
  use constants, only: WP,IIN,IOUT,PI,EARTHRADIUS
  use precisions, only: USE_OVERTIME,OVERTIME_PERCENT, &
    PRESCRIBED_SOURCE,FILTER_SEISMOGRAMS,STATION_CORRECTION,ROTATE_FRAME
  use parameter_inputs
  implicit none
  logical:: Phaseshift_Program     = .false.
  logical:: HetPhaseshift_Program  = .false.
  logical:: referenceRun           = .false.

  ! simulation parameters
  integer:: startVertex,endVertex,firsttimestep,lasttimestep, &
            numofTimeSteps,midpointVertex
  integer:: sourceVertex,receiverVertex,deltaVertex,originSourceVertex, &
            originReceiverVertex
  integer:: numofKernels,currentKernel
  integer,allocatable,dimension(:):: receivers
  integer,allocatable,dimension(:,:):: kernelsReceivers
  real(WP):: desiredSourceLat,desiredSourceLon,muSquare,muTwo
  real(WP):: desiredReceiverLat,desiredReceiverLon
  real(WP):: dt = 0.0_WP
  real(WP):: dt2,cphase2
  real(WP):: distanceQuit
  real(WP):: benchAllStart,benchAllEnd,benchstart,benchend
  real(WP),allocatable,dimension(:,:):: seismogramReceiver, &
                                        receiversSeismogram, &
                                        receiversSeismogramRef
  real(WP),allocatable,dimension(:,:,:):: kernelsReceiversSeismogram, &
                                          kernelsReceiversSeismogramRef
  real(WP),allocatable,dimension(:,:):: forceTermPrescribed
  logical:: sourceOnFile = .false.
  integer:: sourceFileID = 201
  integer:: SH_lmx,SH_ncoef
  real,allocatable,dimension(:):: SH_coef
end module

!-----------------------------------------------------------------------
module deltaSecondLocation
!-----------------------------------------------------------------------
! additional delta location
  use constants, only: WP
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
  use parameter_inputs, only: PARALLELSEISMO
  implicit none
  logical:: MAIN_PROCESS   = .false.
  integer:: nprocesses
  integer:: myrank
  integer:: MPI_CUSTOM
end module

!-----------------------------------------------------------------------
module displacements
!-----------------------------------------------------------------------
! displacement arrays
  use constants, only: WP
  implicit none
  real(WP),allocatable,dimension(:):: displacement
  real(WP),allocatable,dimension(:):: displacement_old
  real(WP),allocatable,dimension(:):: newdisplacement
end module

!-----------------------------------------------------------------------
module griddomain
!-----------------------------------------------------------------------
! grid domain array
  use constants, only: WP,PI
  implicit none
  integer:: boundariesMaxRange
  integer,allocatable,dimension(:):: domainVertices
  integer,allocatable,dimension(:):: vertexDomain
  integer,allocatable,dimension(:,:,:):: boundaries
  integer,allocatable,dimension(:,:):: domainNeighbors
  real(WP),allocatable,dimension(:):: sendDisp,receiveDisp
end

!-----------------------------------------------------------------------
module loop
!-----------------------------------------------------------------------
! looping delta location
  use parameter_inputs, only: latitudeStart,latitudeEnd,longitudeEnd
  implicit none
end

!-----------------------------------------------------------------------
module phaseBlockData
!-----------------------------------------------------------------------
! phase module
  use constants, only: WP,IIN
  use precisions, only: compatible_refgridXrefgrid
  use parameter_inputs, only: phaseBlockFile,heterogeneousPixelsize,phaseBlockVelocityReference, &
    heterogeneousDataFile,heterogeneousOutput,gsh_maximum_expansion
  implicit none
  integer:: numBlocks
  real(WP), allocatable, dimension(:):: phaseBlock
end module

!-----------------------------------------------------------------------
module traveltime
!-----------------------------------------------------------------------
! traveltime common parameters
  use constants, only: WP
  implicit none
  real(WP):: t_lag,arrivalTime,vertexCellArea,t_lagAnalytic,seismo_timestep
  real(WP):: phasevelocity,kernel,kernelAnalytic
end module

!-----------------------------------------------------------------------
module filterType
!-----------------------------------------------------------------------
  use constants, only: WP
  implicit none
  ! window size for fft
  integer:: WindowSIZE = 0
  ! full bandwidth (in frequency domain)
  integer:: bw_width = 0
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

