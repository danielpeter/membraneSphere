!=====================================================================
!
!       m e m b r a n e S p h e r e
!       --------------------------------------------------
!      (c) 2025
!
!=====================================================================

!-----------------------------------------------------------------------
  subroutine readParameters()
!-----------------------------------------------------------------------
! read in parameters from file 'Parameter_Input'
  use constants, only: IIN
  use parameter_inputs
  use adjointVariables, only: Adjoint_Program,Adjoint_InversionProgram
  use propagationStartup, only: HetPhaseshift_Program
  use precisions, only: DO_CHECKERBOARD,MAP_DEGREE_L,MAP_DEGREE_M
  implicit none
  ! local parameters
  integer:: ier,tmpInteger,istart,iend
  character(len=35):: inputName,parName
  character(len=150):: line
  character(len=128):: tmp

  ! input parameter
  inputName = 'Parameter_Input'
  inputName = trim(inputName)
  !if (VERBOSE) print *,inputName
  open(IIN, file=trim(inputName),status='old',iostat=ier)
  if (ier /= 0) then
    print *,'Error: opening file ',trim(inputName)
    call stopProgram('Abort - readParameters() opening input    ')
  endif

  do while( ier == 0)
    read(IIN,'(a150)', iostat=ier) line
    if (ier /= 0) exit

    ! suppress leading white spaces, if any
    line = adjustl(line)

    ! suppress trailing carriage return (ASCII code 13) if any (e.g. if input text file coming from Windows/DOS)
    if (index(line,achar(13)) > 0) line = line(1:index(line,achar(13))-1)

    ! suppress trailing comments " .. ! my comment"
    if (index(line,'#') > 5) line = line(1:index(line,'#')-1)
    if (index(line,'!') > 5) line = line(1:index(line,'!')-1)

    line = trim(line)

    ! check line
    if (len_trim(line) == 0) cycle

    ! check if comment line
    if (line(1:1) == " " .or. line(1:1) == "!" .or. line(1:1) == "%") cycle

    ! get index of "=" sign
    istart = index(line,'=')
    if (istart == 0) cycle

    ! check if parameter name valid
    if (istart < 5) then
      print *,'Error: line with wrong format         : ***'//trim(line)//'***'
      print *,'       start index for parameter value: istart = ',istart
      call stopProgram('Abort - line with wrong format in input file    ')
    endif

    ! parameter name
    parName = line(1:istart-1)
    parName = trim(adjustl(parName))

    ! index range for parameter values
    istart = istart + 1
    iend = len_trim(line)

    !debug
    !print *,'debug: index range istart/iend = ',istart,'/',iend
    !print *,'debug: parameter string  = ***'//line(istart:iend)//'***'

    ! read in parameter
    select case(trim(parName))
    ! Physical model
    case('LEVEL')
      read(line(istart:iend),*) subdivisions
    case('FIRSTTIME')
      read(line(istart:iend),*) FIRSTTIME
    case('LASTTIME')
      read(line(istart:iend),*) LASTTIME
    case('CPHASE')
      read(line(istart:iend),*) cphasetype
      ! trim & suppress leading white spaces, if any
      cphasetype = trim(adjustl(cphasetype))
      call determinePhaseRef(cphasetype,8,cphaseRef)

    ! Geometry
    case('SOURCE')
      read(line(istart:iend),*) sourceLat,sourceLon
    case('RECEIVER')
      read(line(istart:iend),*) receiverLat,receiverLon
    case('MANYRECEIVERS')
      read(line(istart:iend),*) manyReceivers
    case('MANYNUMOFRECEIVERS')
      read(line(istart:iend),*) numofReceivers
    case('MANYKERNELS')
      read(line(istart:iend),*) manyKernels
    case('KRNEPI')
      read(line(istart:iend),*) kernelStartDistance,kernelEndDistance
    case('IMPORTKERNELSRECEIVERS')
      read(line(istart:iend),*) importKernelsReceivers

    ! Source
    case('SOURCE_TIME_SIGMA')
      read(line(istart:iend),*) TimeParameterSigma
    case('SOURCE_WIDTH_MU')
      read(line(istart:iend),*) WidthParameterMu
    case('FILTER_INITIALSOURCE')
      read(line(istart:iend),*) FILTER_INITIALSOURCE
    ! adjoint source
    case('ADJOINT_TAPERSIGNAL')
      read(line(istart:iend),*) ADJOINT_TAPERSIGNAL
    case('WINDOWED_INTEGRATION')
      read(line(istart:iend),*) WINDOWED_INTEGRATION
    case('WINDOW_START')
      read(line(istart:iend),*) WINDOW_START
    case('WINDOW_END')
      read(line(istart:iend),*) WINDOW_END
    case('ADJOINT_ANTIPODE_TIME')
      read(line(istart:iend),*) ADJOINT_ANTIPODE_TIME

    ! Scatterer
    case('DELTA')
      read(line(istart:iend),*) DELTA
    case('DRADIUS')
      read(line(istart:iend),*) DELTARADIUS
    case('DTYPE')
      read(line(istart:iend),*) DELTAfunction
      DELTAfunction = trim(DELTAfunction)
    case('DPERTURBATION')
      read(line(istart:iend),*) deltaPerturbation
    case('DLOCATION')
      read(line(istart:iend),*) deltaLat,deltaLon
    case('MOVEDELTA')
      read(line(istart:iend),*) MOVEDELTA
    case('DLATSTART')
      read(line(istart:iend),*) latitudeStart
    case('DLATEND')
      read(line(istart:iend),*) latitudeEnd
    case('DLONEND')
      read(line(istart:iend),*) longitudeEnd
    case('DINCREMENT')
      read(line(istart:iend),*) deltaMoveIncrement
    case('SECONDDELTA')
      read(line(istart:iend),*) SECONDDELTA

    ! heterogeneous model
    case('HETEROGENEOUS')
      read(line(istart:iend),*) HETEROGENEOUS
    case('BLKFILE')
      read(line(istart:iend),*) phaseBlockFile
      ! trim & suppress leading white spaces, if any
      phaseBlockFile = trim(adjustl(phaseBlockFile))
    case('BLK/INV_PIXELSIZE')
      read(line(istart:iend),*) heterogeneousPixelsize
      if (.not. DO_CHECKERBOARD) then
        if (phaseBlockFile(len_trim(phaseBlockFile)-2:len_trim(phaseBlockFile)) == 'gsh') then
          gsh_maximum_expansion = heterogeneousPixelsize
          heterogeneousPixelsize = 0
        endif
      endif
    case('BLKVELOCITYREFERENCE')
      read(line(istart:iend),*) tmp
      ! trim & suppress leading white spaces, if any
      tmp = trim(adjustl(tmp))
      call determinePhaseRef(tmp,128,phaseBlockVelocityReference)
    case('INV_DATA')
      read(line(istart:iend),*) heterogeneousDataFile
      ! trim & suppress leading white spaces, if any
      heterogeneousDataFile = trim(adjustl(heterogeneousDataFile))
    case('INV_OUTPUT')
      read(line(istart:iend),*) heterogeneousOutput
      ! trim & suppress leading white spaces, if any
      heterogeneousOutput = trim(adjustl(heterogeneousOutput))
    case('INV_VOXELSIZE')
      read(line(istart:iend),*) tmpInteger
      if (phaseBlockFile(len_trim(phaseBlockFile)-2:len_trim(phaseBlockFile)) /= 'gsh') then
        if (tmpInteger /= floor(heterogeneousPixelsize+0.1)) then
          print *,'Error: inversion grid size and block phase velocity map are different', &
                  tmpInteger,floor(heterogeneousPixelsize+0.1)
          call stopProgram('abort - readParameters ')
        endif
      endif
      heterogeneousPixelsize = tmpInteger

    ! file output & misc
    case('DATADIRECTORY')
      read(line(istart:iend),*) datadirectory
      ! trim & suppress leading white spaces, if any
      datadirectory = trim(adjustl(datadirectory))
    case('ADJOINTKERNEL')
      read(line(istart:iend),*) adjointKernelName
      ! trim & suppress leading white spaces, if any
      adjointKernelName = trim(adjustl(adjointKernelName))
    case('VERBOSE')
      read(line(istart:iend),*) VERBOSE
    case('SIMULATIONOUTPUT')
      read(line(istart:iend),*) SIMULATIONOUTPUT
    case('PARALLELSEISMO')
      read(line(istart:iend),*) PARALLELSEISMO
    end select
  enddo
  close(IIN)

  ! info output
  if (VERBOSE) then
    print *,'level         ',subdivisions
    print *,'firsttimestep ',FIRSTTIME
    print *,'lasttimestep  ',LASTTIME
    print *,'cphase        ', cphasetype, cphaseRef
    print *,'source        ',sourceLat,sourceLon
    print *,'receiver      ',receiverLat,receiverLon
    print *,'many receiver stations', manyReceivers
    if (manyReceivers) print *,'  receiver stations ',numofReceivers
    print *,'many kernels          ', manyKernels
    if (manyKernels) then
      print *,'  epicentral distances ',kernelStartDistance,kernelEndDistance
      print *,'  importing receivers  ', importKernelsReceivers
    endif
    print *,'delta ', DELTA
    if (DELTA) then
      print *,'  delta radius       ',DELTARADIUS
      print *,'  delta type         ',DELTAfunction
      print *,'  delta perturbation ',deltaPerturbation
      print *,'  delta location     ',deltaLat,deltaLon
      print *,'  move delta ', MOVEDELTA
      if (MOVEDELTA) then
        print *,'    delta start latitude ',latitudeStart
        print *,'    delta end latitude   ',latitudeEnd
        print *,'    delta end longitude  ',longitudeEnd
        print *,'    delta move increment ',deltaMoveIncrement
      endif
      print *,'  second delta ',SECONDDELTA
    endif
    print *,'heterogeneous ', HETEROGENEOUS
    if (HETEROGENEOUS) then
      if (DO_CHECKERBOARD) then
        print *,'  checkerboard used  : L=',MAP_DEGREE_L,'M=',MAP_DEGREE_M
      else
        print *,'  phaseBlockFile used: ',trim(phaseBlockFile)
        if (phaseBlockFile(len_trim(phaseBlockFile)-2:len_trim(phaseBlockFile)) == 'gsh') then
          print *,'  maximum degree expansion ',gsh_maximum_expansion
        else
          print *,'  grid pixel size ',heterogeneousPixelsize
        endif
      endif
      print *,'  phaseBlock VelocityReference',phaseBlockVelocityReference
    endif
    if (HetPhaseshift_Program .or. Adjoint_InversionProgram) then
      print *,'  data input                : ',trim(heterogeneousDataFile)
      print *,'  output files              : ',trim(heterogeneousOutput)
      print *,'  inversion Grid pixel size : ',heterogeneousPixelsize
    endif
    print *,'data output directory : ',trim(datadirectory)
    if (Adjoint_Program) print *,'adjoint kernel name : ',trim(adjointKernelName)
    print *,'verbose ',VERBOSE
    print *,'simulationoutput ',SIMULATIONOUTPUT
    print *,'parallelize single simulation ',PARALLELSEISMO
  endif

  ! check parameters
  call checkParameters()

  end subroutine


!-----------------------------------------------------------------------
  subroutine determinePhaseRef(ctype,length,cphase)
!-----------------------------------------------------------------------
! from phase type (e.g. L150) determine the reference phase speed
  use precisions
  implicit none
  integer,intent(in):: length
  character(len=length),intent(in):: ctype
  real(WP),intent(out):: cphase

  ! initialize
  cphase = -1.0

  ! gets phase speed
  select case(trim(ctype))

  ! rayleigh waves
  case('R35')
    cphase = PHASEVELOCITY_R35
  case('R37')
    cphase = PHASEVELOCITY_R37
  case('R40')
    cphase = PHASEVELOCITY_R40
  case('R45')
    cphase = PHASEVELOCITY_R45
  case('R50')
    cphase = PHASEVELOCITY_R50
  case('R60')
    cphase = PHASEVELOCITY_R60
  case('R75')
    cphase = PHASEVELOCITY_R75
  case('R80')
    cphase = PHASEVELOCITY_R80
  case('R100')
    cphase = PHASEVELOCITY_R100
  case('R150')
    cphase = PHASEVELOCITY_R150
  case('R200')
    cphase = PHASEVELOCITY_R200
  case('R250')
    cphase = PHASEVELOCITY_R250
  case('R300')
    cphase = PHASEVELOCITY_R300

  ! love waves
  case('L35')
    cphase = PHASEVELOCITY_L35
  case('L37')
    cphase = PHASEVELOCITY_L37
  case('L40')
    cphase = PHASEVELOCITY_L40
  case('L45')
    cphase = PHASEVELOCITY_L45
  case('L50')
    cphase = PHASEVELOCITY_L50
  case('L60')
    cphase = PHASEVELOCITY_L60
  case('L75')
    cphase = PHASEVELOCITY_L75
  case('L80')
    cphase = PHASEVELOCITY_L80
  case('L100')
    cphase = PHASEVELOCITY_L100
  case('L150')
    cphase = PHASEVELOCITY_L150
  case('L200')
    cphase = PHASEVELOCITY_L200
  case('L250')
    cphase = PHASEVELOCITY_L250
  case('L300')
    cphase = PHASEVELOCITY_L300
  case('L450')
    cphase = PHASEVELOCITY_L450
  case('L600')
    cphase = PHASEVELOCITY_L600

  case default
    print *,'Error: phase speed not recognized:',trim(ctype)
    call stopProgram('phase speed not recognized    ')
  end select

  if (cphase < 0.0) then
    print *,'Error: phase speed not recognized:',trim(ctype)
    call stopProgram('abort - readParameters')
  endif

  end subroutine

!-----------------------------------------------------------------------
  subroutine checkParameters()
!-----------------------------------------------------------------------
! checks parameters read in from file 'Parameter_Input'
  use parameter_inputs
  implicit none

  ! for now, only checks subdivision level

  ! checks integer overflow
  !   - 32-bit integer < 2,147,483,647
  !   - 64-bit integer < 9,223,372,036,854,775,807
  !
  ! for 32-bit integers: MAXTRIANGLES limit is reached at MAXLEVELS == 12
  !                      that is, for levels == 13 MAXTRIANGLES will overflow
  if (subdivisions > 12) then
    print *,'Error: the maximum of supported subdivision levels is equal to 12.'
    print *,'       This is due to integer overflow (of the number of triangles) beyond these levels.'
    print *
    print *,'Please choose a subdivision level <= 12 in your Parameter_Input file, exiting...'
    print *
    call stopProgram('Invalid subdivisions requested, maximum of subdivisions == 12')
  endif

  end subroutine
