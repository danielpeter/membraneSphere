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
  subroutine readParameters()
!-----------------------------------------------------------------------
! read in parameters from file 'Parameter_Input'
  use propagationStartup;use loop;use phaseBlockData;use parallel
  use verbosity; use adjointVariables; use cells
  implicit none
  ! local parameters
  integer:: ier,tmpInteger,istart,iend
  character(len=32):: inputName
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

    ! index range for parameter values
    istart = istart + 1
    iend = len_trim(line)

    !debug
    !print *,'debug: index range istart/iend = ',istart,'/',iend
    !print *,'debug: parameter string  = ***'//line(istart:iend)//'***'

    ! read in parameter
    select case(line(1:5))
    case('VERBO')
      read(line(istart:iend),*) VERBOSE
      if (VERBOSE) print *,'verbose output'
    case('LEVEL')
      read(line(istart:iend),*) subdivisions
      if (VERBOSE) print *,'level',subdivisions
    case('FIRST')
      read(line(istart:iend),*) FIRSTTIME
      if (VERBOSE) print *,'firsttimestep',FIRSTTIME
    case('LASTT')
      read(line(istart:iend),*) LASTTIME
      if (VERBOSE) print *,'lasttimestep',LASTTIME
    case('CPHAS')
      read(line(istart:iend),*) cphasetype
      ! suppress leading white spaces, if any
      cphasetype = adjustl(cphasetype)
      ! trim
      cphasetype = trim(cphasetype)
      call determinePhaseRef(cphasetype,8,cphaseRef)
      if (VERBOSE) print *,'cphase    ', cphasetype, cphaseRef
    case('SOURC')
      read(line(istart:iend),*) sourceLat,sourceLon
      if (VERBOSE) print *,'source',sourceLat,sourceLon
    case('RECEI')
      read(line(istart:iend),*) receiverLat,receiverLon
      if (VERBOSE) print *,'receiver',receiverLat,receiverLon
    case('MANYR')
      read(line(istart:iend),*) manyReceivers
      if (VERBOSE) print *,'many receiver stations', manyReceivers
    case('MANYN')
      read(line(istart:iend),*) numofReceivers
      if (VERBOSE .and. manyReceivers) print *,'    receiver stations',numofReceivers
    case('MANYK')
      read(line(istart:iend),*) manyKernels
      if (VERBOSE) print *,'many kernels', manyKernels
    case('KRNEP')
      read(line(istart:iend),*) kernelStartDistance,kernelEndDistance
      if (VERBOSE .and. manyKernels) print *,'    epicentral distances',kernelStartDistance,kernelEndDistance
    case('IMPOR')
      read(line(istart:iend),*) importKernelsReceivers
      if (VERBOSE .and. manyKernels) print *,'    importing receivers', importKernelsReceivers
    case('DELTA')
      read(line(istart:iend),*) DELTA
      if (VERBOSE) print *,'delta', DELTA
    case('MOVED')
      read(line(istart:iend),*) MOVEDELTA
      if (VERBOSE .and. DELTA) print *,'  move delta', MOVEDELTA
    case('DRADI')
      read(line(istart:iend),*) DELTARADIUS
      if (VERBOSE .and. DELTA) print *,'  delta radius',DELTARADIUS
    case('DTYPE')
      read(line(istart:iend),*) DELTAfunction
      DELTAfunction=trim(DELTAfunction)
      if (VERBOSE .and. DELTA) print *,'  delta type ',DELTAfunction
    case('DLOCA')
      read(line(istart:iend),*) deltaLat,deltaLon
      if (VERBOSE .and. DELTA) print *,'  delta location',deltaLat,deltaLon
    case('DPERT')
      read(line(istart:iend),*) deltaPerturbation
      if (VERBOSE .and. DELTA) print *,'  delta perturbation',deltaPerturbation
    case('DLATS')
      read(line(istart:iend),*) latitudeStart
      if (VERBOSE .and. DELTA) print *,'  delta start latitude',latitudeStart
    case('DLATE')
      read(line(istart:iend),*) latitudeEnd
      if (VERBOSE .and. DELTA) print *,'  delta end latitude',latitudeEnd
    case('DLONE')
      read(line(istart:iend),*) longitudeEnd
      if (VERBOSE .and. DELTA) print *,'  delta end longitude',longitudeEnd
    case('DINCR')
      read(line(istart:iend),*) deltaMoveIncrement
      if (VERBOSE .and. DELTA) print *,'  delta move increment',deltaMoveIncrement
    case('SECON')
      read(line(istart:iend),*) SECONDDELTA
      if (VERBOSE .and. DELTA) print *,'  second delta:',SECONDDELTA
    case('HETER')
      read(line(istart:iend),*) HETEROGENEOUS
      if (VERBOSE) print *,'heterogeneous', HETEROGENEOUS
    case('BLKFI')
      read(line(istart:iend),*) line
      phaseBlockFile = trim(line)
      ! suppress leading white spaces, if any
      phaseBlockFile = adjustl(phaseBlockFile)
      if (DO_CHECKERBOARD) then
        if (VERBOSE .and. HETEROGENEOUS) &
          print *,'  checkerboard used: L=',MAP_DEGREE_L,'M=',MAP_DEGREE_M
      else
        if (VERBOSE .and. HETEROGENEOUS) &
          print *,'  phaseBlockFile used: ',trim(phaseBlockFile)
      endif
    case('BLK/I')
      read(line(istart:iend),*) heterogeneousPixelsize
      if (.not. DO_CHECKERBOARD) then
        if (phaseBlockFile(len_trim(phaseBlockFile)-2:len_trim(phaseBlockFile)) == 'gsh') then
          gsh_maximum_expansion = heterogeneousPixelsize
          heterogeneousPixelsize = 0
          if (VERBOSE .and. HETEROGENEOUS) print *,'    maximum degree expansion ',gsh_maximum_expansion
        else
          if (VERBOSE .and. HETEROGENEOUS) print *,'  Grid pixel size ',heterogeneousPixelsize
        endif
      endif
    case('BLKVE')
      read(line(istart:iend),*) tmp
      tmp = trim(tmp)
      ! suppress leading white spaces, if any
      tmp = adjustl(tmp)
      call determinePhaseRef(tmp,128,phaseBlockVelocityReference)
      if (VERBOSE .and. HETEROGENEOUS) print *,'  phaseBlock VelocityReference',phaseBlockVelocityReference
    case('INV_D')
      read(line(istart:iend),*) line
      heterogeneousDataFile = trim(line)
      ! suppress leading white spaces, if any
      heterogeneousDataFile = adjustl(heterogeneousDataFile)
      if (VERBOSE .and. (HetPhaseshift_Program .or. Adjoint_InversionProgram)) &
        print *,'  Data input: ',trim(heterogeneousDataFile)
    case('INV_O')
      read(line(istart:iend),*) line
      heterogeneousOutput = trim(line)
      ! suppress leading white spaces, if any
      heterogeneousOutput = adjustl(heterogeneousOutput)
      if (VERBOSE .and. (HetPhaseshift_Program .or. Adjoint_InversionProgram)) &
        print *,'  Output files: ',trim(heterogeneousOutput)
    case('INV_V')
      read(line(istart:iend),*) tmpInteger
      if (phaseBlockFile(len_trim(phaseBlockFile)-2:len_trim(phaseBlockFile)) /= 'gsh') then
        if (tmpInteger /= floor(heterogeneousPixelsize+0.1)) then
          print *,'Error: inversion grid size and block phase velocity map are different', &
                  tmpInteger,floor(heterogeneousPixelsize+0.1)
          call stopProgram('abort - readParameters ')
        endif
      endif
      heterogeneousPixelsize = tmpInteger
      if (VERBOSE .and. (HetPhaseshift_Program .or. Adjoint_InversionProgram)) &
        print *,'  Inversion Grid pixel size ',heterogeneousPixelsize
    case('SIMUL')
      read(line(istart:iend),*) SIMULATIONOUTPUT
      if (VERBOSE) print *,'simulationoutput',SIMULATIONOUTPUT
    case('DATAD')
      read(line(istart:iend),*) tmp
      datadirectory = trim(tmp)
      ! suppress leading white spaces, if any
      datadirectory = adjustl(datadirectory)
      if (VERBOSE) print *,'data output directory : ',trim(datadirectory)
    case('ADJOI')
      read(line(istart:iend),*) tmp
      adjointKernelName = trim(tmp)
      ! suppress leading white spaces, if any
      adjointKernelName = adjustl(adjointKernelName)
      if (VERBOSE .and. Adjoint_Program) &
        print *,'adjoint kernel name : ',trim(adjointKernelName)
    case('PARAL')
      read(line(istart:iend),*) PARALLELSEISMO
      if (VERBOSE) print *,'parallelize single simulation',PARALLELSEISMO
    end select
  enddo

  !if (VERBOSE) print *,'delta location(lat/lon):',deltaLat,deltaLon
  close(IIN)

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
