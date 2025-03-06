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
!
! modules needed by e.g. program propagation

      ! working precision & parameters
      ! these values are not set in 'Parameter_Input'      
      module precision
        !-------------------------------------------------------------------------------------------------------------
        ! user can modify these settings
        ! working precision in bytes: real = 4,double precision = 8 
        integer, parameter :: WP                   = 4            
  
        ! source parameters
        ! fixed_sourceparameter: use fixed source parameters, important when a higher grid level 
        !                       is choosen but the same source should be simulated
        ! TimeParameterSigma: level 4: 164.6128 ; L150: 41.1795264700414 ; L75: 42.9967494824469
        !                                      TimeParameterSigma = 42.0216_WP       
        !                                      good other values: q4: 204.5, q6: 51.125, 42.0216(L150)
        ! WidthParameterMu: level 4: 7.1343616E-02 ; L150: 1.783786819208144E-002 ; 
        !                                  L75: 1.783786819208144E-002
        !                                  WidthParameterMu = 0.017968_WP         
        !                                  good other values q4:0.0713 (20 degree), 
        !                                                       q6: 0.017826 (5 degree),0.017968 (factor 8x)
        ! factor_time: these factors are taken to calculate source parameters sigma/mu 
        !                     when no FIXED_SOURCEPARAMETER
        !                     for Love L150 - subdivision 8: fT6.0,fW8.0 have main period at 139.6 s ,
        !                                              subdivision 7: fT4.0,fW10.0 main period at 206.9 s
        !                     (check with waveparameters executable)                                
        logical,parameter :: FIXED_SOURCEPARAMETER  = .true. 
        real(WP) :: TimeParameterSigma            = 60.0_WP  
        real(WP) :: WidthParameterMu              = 4.0E-02  
        real(WP),parameter :: FACTOR_TIME          = 4.0_WP     
        real(WP),parameter :: FACTOR_WIDTH         = 8.0_WP          

        ! filter parameters
        ! FILTERINITIALSOURCE: filter initial prescribed source too
        ! FILTERSEISMOGRAMS: use filtering when calculating timelag
        ! bw_fixfrequency: use the fixed halfband-frequency or a percentage of the wave period
        ! half-bandwidth frequency: use same bandwidth-half as for analytic comparision 
        !                                           (compare: Kernels/BORN/find_kernels.f, parameter deltau)
        ! bandwidth percent: 10% of corner period used as filter width on both sides 
        !                                ( e.g. Love 150s gets +- 15 s filter ), 0.6059 for half-bandwidth of 2.5 mHz        
        ! butterworthfilter: just use bandpass filter or butterworth bandpass filter
        ! butterworth_power: butterworth power exponent
        ! arrival_threshold: threshold value for arrival time picking also determines startTime 
        !                             for fourier transformation        
        ! hannwindow_percent: hanning window width to apply at the ends of the seismogram 
        !                                     in percent of the seismogram length            
        ! membranecorrection: [s] to correct wave period in filtering routines; 
        !                                   -10.0 for grid level 6 (empirical value; compare with single 
        !                                   frequency kernel cross-sections)        
        ! ANALYTICAL_CORRELATION: calculates also the analytically derived timelag & 
        !                                                   kernel value
        logical,parameter :: FILTERINITIALSOURCE    = .true.    
        logical,parameter :: FILTERSEISMOGRAMS      = .false.         
        logical,parameter :: BW_FIXFREQUENCY        = .true.     
        real(WP),parameter :: BW_HALFFREQUENCY     = 2.5e-3 
        real(WP),parameter :: BW_PERCENT           = 0.0_WP     
        logical,parameter :: BUTTERWORTHFILTER      = .false.   
        integer,parameter :: BUTTERWORTH_POWER      = 1         
        real(WP),parameter :: ARRIVAL_THRESHOLD    = 1.0_WP    
        real(WP),parameter :: HANNWINDOW_PERCENT   = 0.1        
        real(WP),parameter :: MEMBRANECORRECTION   = -0.0      
        logical,parameter :: ANALYTICAL_CORRELATION = .false.   
        
        ! routine optimization
        ! precalculated_cells: speeding up computation by using precalculated cell areas 
        !                                 and distances
        ! prescribedsource: precalculates source values and stores in an array for faster 
        !                              computation of finite-difference iteration 
        logical,parameter::PRECALCULATED_CELLS      = .true.         
        logical:: PRESCRIBEDSOURCE                 = .true.              
        
        ! analytical spherical harmonics degrees
        ! degree's l & m : spherical harmonic degrees for comparision with laplacian accuracy
        ! degreeAccuracy : legendre accuracy for analytic solution of seismogram with source 2
        !                              default: 80; see plot fig. (3.3)B, Tape, p. 32
        integer, parameter :: DEGREE_L              = 6
        integer, parameter :: DEGREE_M              = 1        
        integer, parameter :: DEGREEACCURACY        = 100

        ! inversion
        ! compatible_refgridXrefgrid: set equal to 1 to make the grid compatible with a refgridXrefgrid grid
        !   0 for L0035.jwkb.lsqr.mod.pcn-.2-.2 or L150.crust.2degree.gsh.modelVector.*.pcn.blk
        !   1 for L150.crust
        ! grid pixelsize: inversion grid consists of approximately equal sized pixels with this size (in  degrees)
        integer,parameter:: compatible_refgridXrefgrid  = 0

        ! source/station location
        ! station_correction: tries to use exact locations of source and receiver.
        !                               interpolates linearly for the receiver displacement.
        logical,parameter:: Station_Correction                = .false.

        ! use overtime : for the adaption of the simulation time for different source/receiver setups
        !                         (e.g. when more kernels either uses an overtime or antipode time)
        ! simulation overtime %: how much percent more time shall be given for the overtime
        logical,parameter:: USEOVERTIME                 = .false.    
        real(WP),parameter:: SIMULATIONOVERTIMEPERCENT = 0.8_WP

        
        !-------------------------------------------------------------------------------------------------------------
        ! user shouldn't need to modify the following ones
        ! PREM values
        ! radius measured in km; phase velocities measured in km/s
        ! high periods 450/600: determined with spline representation for wave period 
        ! (see: Kernels/BORN/find_kernels.f)   
        real(WP), parameter :: EARTHRADIUS = 6371.0_WP        
        real(WP), parameter :: PHASEVELOCITY_R40 = 3.92799_WP !rayleigh wave 40s
        real(WP), parameter :: PHASEVELOCITY_R50 = 3.95146_WP    
        real(WP), parameter :: PHASEVELOCITY_R60 = 3.97339_WP                
        real(WP), parameter :: PHASEVELOCITY_R75 = 4.01077_WP                      
        real(WP), parameter :: PHASEVELOCITY_L35 = 4.39987_WP !love wave 35s 
        real(WP), parameter :: PHASEVELOCITY_L50 = 4.49336_WP    
        real(WP), parameter :: PHASEVELOCITY_L75 = 4.57474_WP    
        real(WP), parameter :: PHASEVELOCITY_L100= 4.64431_WP    
        real(WP), parameter :: PHASEVELOCITY_L150= 4.77915_WP ! first arrival phase velocity             
        real(WP), parameter :: PHASEVELOCITY_L200= 4.91928_WP          
        real(WP), parameter :: PHASEVELOCITY_L250= 5.07097_WP          
        real(WP), parameter :: PHASEVELOCITY_L300= 5.22906_WP    
        real(WP), parameter :: PHASEVELOCITY_L450= 5.74979_WP 
        real(WP), parameter :: PHASEVELOCITY_L600= 6.36841_WP 
        real(WP), parameter :: PI = 3.1415926535897931_WP       
        !-------------------------------------------------------------------------------------------------------------        
      end module precision      
!-----------------------------------------------------------------------
      module verbosity
        logical:: VERBOSE     = .true.         !for console verbosity
        logical:: beVerbose   = .true.         !for getTimelag verbosity
        logical:: fileOutput  = .false.        !for debuging: outputs to file 'tmp*.dat'
        logical:: DEBUG       = .false.        !for debugging: code very verbose!
      end module
!-----------------------------------------------------------------------
      ! adjoint module
      ! used for adjoint routines
      module adjointVariables
        use precision
        ! adjoint parameters
        ! onthefly: instead of integration at the end, after each time step 
        ! startatzero: start integration from zero, not for whole seismogram
        ! precalculate derivatives: calculate derivatives after forward iteration has finished, 
        !                                         instead of just before integration       
        ! windowed integration: do integration between zero and arrivaltime for seismogram
        logical,parameter:: ADJOINT_ONTHEFLY            = .false. 
        logical,parameter:: ADJOINT_STARTATZERO         = .false. 
        logical,parameter:: PRECALCULATE_DERIVATIVES    = .false. 
        logical,parameter:: WINDOWEDINTEGRATION         = .false. 
        
        ! adjoint variables
        logical:: Adjoint_Program                   = .false.
        logical:: Adjoint_InversionProgram          = .false.
        logical:: kernelIteration,storeAsFile
        integer:: adjointSourceVertex     
        character*64:: adjointKernelName
        real(WP):: normFactor   
        real(WP),allocatable,dimension(:,:):: adjointSource
        real(WP),allocatable,dimension(:):: adjointKernel
        real(WP),allocatable,dimension(:):: backwardDisplacement,backwardDisplacement_old,backwardNewdisplacement
        real(WP),allocatable,dimension(:,:):: wavefieldForward, wavefieldAdjoint,seismoSecondDerivative      
        integer,parameter:: adjSourceFileID=101,adjRecFileID=102,adjMidpointFileID=103
      end module
!-----------------------------------------------------------------------
      ! cells module
      ! used for precalculated values of cell areas etc.
      module cells  
        use precision
        integer:: MaxTriangles,MaxVertices
        integer,allocatable, dimension(:,:):: cellFace,cellNeighbors,cellTriangleFace
        integer:: numFaces,numNeighbors,numTriangleFaces,numCorners
        real(WP), allocatable, dimension(:)::cellAreas
        real(WP), allocatable, dimension(:,:)::vertices,cellEdgesLength,cellCenterDistances,cellCorners
        real(WP):: interpolation_distances(3),interpolation_triangleLengths(3)
        integer:: interpolation_corners(3),interpolation_triangleIndex
      end module cells    
!-----------------------------------------------------------------------
      ! phaseVelocityMap module
      ! used for phase velocity maps
      module phaseVelocityMap
        use precision
        integer:: numPhaseEntries
        real(WP), allocatable, dimension(:):: phaseMap
        real(WP), allocatable, dimension(:):: phaseVelocitySquare
      end module
!-----------------------------------------------------------------------
      ! propagationStartup module
      ! used for initialization startup process
      module propagationStartup
        use precision
        logical:: Phaseshift_Program = .false.
        logical:: HETEROGENEOUS,DELTA,SIMULATIONOUTPUT,SECONDDELTA,MOVEDELTA    
        logical:: manyReceivers,manyKernels
        logical:: importKernelsReceivers,referenceRun
        integer:: startVertex,endVertex,firsttimestep,lasttimestep,numofTimeSteps,midpointVertex
        integer:: sourceVertex,receiverVertex,deltaVertex,originSourceVertex,originReceiverVertex
        integer:: numVertices,numDomainVertices,subdivisions,numofReceivers,numofKernels,currentKernel   
        integer,allocatable,dimension(:):: receivers
        integer,allocatable,dimension(:,:)::kernelsReceivers
        character*8:: DELTAFUNCTION,cphasetype
        character*128:: datadirectory
        real(WP):: FIRSTTIME,LASTTIME,DELTARADIUS
        real(WP):: sourceLat,sourceLon,receiverLat,receiverLon,cphaseRef,desiredSourceLat,desiredSourceLon
        real(WP):: deltaLat,deltaLon,deltaPerturbation,deltaMoveIncrement,desiredReceiverLat,desiredReceiverLon 
        real(WP):: cphase2,dt,dt2,averageCellDistance
        real(WP):: kernelStartDistance,kernelEndDistance,distanceQuit
        real(WP):: benchAllStart,benchAllEnd,benchstart,benchend
        real(WP),allocatable,dimension(:,:):: seismogramReceiver,receiversSeismogram,receiversSeismogramRef                
        real(WP),allocatable,dimension(:,:,:)::kernelsReceiversSeismogram,kernelsReceiversSeismogramRef
        real(WP),allocatable,dimension(:,:)::forceTermPrescribed
        logical:: sourceOnFile = .false.
        integer:: sourceFileID = 201 
        integer::SH_lmx,SH_ncoef
        real,allocatable,dimension(:)::SH_coef
        logical:: rotate_frame = .true.  ! rotate source&receiver locations/heterogeneous phase map to equatorial plane (by default)
      end module
!-----------------------------------------------------------------------
      ! additional delta location
      module deltaSecondLocation
        use precision
        integer:: deltaSecondVertex
        real(WP):: deltaSecondLat,deltaSecondLon
        real(WP),parameter:: deltaSecondDistance = 3.0
      end module
!-----------------------------------------------------------------------
      ! parallel module
      ! used for MPI parallelization
      module parallel
        include 'mpif.h'
        logical:: PARALLELSEISMO,MASTER
        integer:: nprocesses,rank,tag,MPI_CUSTOM
        integer:: status(MPI_STATUS_SIZE)
      end module
!-----------------------------------------------------------------------
      ! displacement arrays
      module displacements
        use precision
        real(WP),allocatable,dimension(:):: displacement
        real(WP),allocatable,dimension(:):: displacement_old
        real(WP),allocatable,dimension(:):: newdisplacement
      end module
!-----------------------------------------------------------------------
      ! grid domain array
      module griddomain
        use precision
        integer boundariesMaxRange
        integer,allocatable,dimension(:)::domainVertices
        integer,allocatable,dimension(:)::vertexDomain
        integer,allocatable,dimension(:,:,:)::boundaries
        integer,allocatable,dimension(:,:)::domainNeighbors
        real(WP),allocatable,dimension(:)::sendDisp,receiveDisp
      end
!-----------------------------------------------------------------------
      ! looping delta location
      module loop
        integer:: latitudeStart, latitudeEnd
        integer:: longitudeEnd
      end
!-----------------------------------------------------------------------      
      ! phase module
      module phaseBlockData
        use precision
        integer:: numBlocks
        real(WP), allocatable, dimension(:):: phaseBlock
        real(WP):: phaseBlockVelocityReference
        character*128:: phaseBlockFile,heterogeneousDataFile,heterogeneousOutput
        real:: heterogeneousPixelsize        
      end module                  
!-----------------------------------------------------------------------
      ! traveltime common parameters
      module traveltime
        use precision
        real(WP):: t_lag,arrivalTime,vertexCellArea,t_lagAnalytic,seismo_timestep
        real(WP):: phasevelocity,kernel,kernelAnalytic
      end module
!-----------------------------------------------------------------------      
      module filterType
        use precision
        integer:: WindowSIZE                    ! window size for fft
        integer:: bw_width                      ! full bandwidth (in frequency domain)
        real(WP):: bw_frequency,bw_waveperiod  ! half-bandwidth frequency, wave period used for bandwidth filtering
      end module      
!-----------------------------------------------------------------------      
      module splineFunction
        integer:: i1,i2
        double precision, allocatable,dimension(:):: X,Y
        double precision, allocatable,dimension(:,:):: Q,F
      end module          
!-----------------------------------------------------------------------      
      module heterogeneousMatrix
        real:: Amin,Amax
      end       
