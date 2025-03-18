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
      integer:: length,ierror,tmpInteger
      character(len=32):: inputName
      character(len=150):: line
      character(len=128):: tmp

      ! input parameter
      inputName = 'Parameter_Input'
      inputName = trim(inputName)
      !if (VERBOSE) print *,inputName
      open(IIN, file=trim(inputName),status='old',iostat=ierror)
      if (ierror /= 0) then
        print *,'Error: opening file ',trim(inputName)
        call stopProgram( 'abort - readParameters() opening input    ')
      endif

      do while( ierror == 0)
        read(IIN,'(a150)', iostat=ierror) line
        if (ierror /= 0) exit

        length = len_trim(line)
        if (length == 0) then
          continue
        else
          if (line(1:1) == " " .or. line(1:1) == "!") then
            continue
          else
            select case(line(1:5))
            case('VERBO')
              read(line(35:len_trim(line)),*) VERBOSE
              if (VERBOSE) print *,'verbose output'
            case('LEVEL')
              read(line(35:len_trim(line)),*) subdivisions
              if (VERBOSE) print *,'level',subdivisions
            case('FIRST')
              read(line(35:len_trim(line)),*) FIRSTTIME
              if (VERBOSE) print *,'firsttimestep',FIRSTTIME
            case('LASTT')
              read(line(35:len_trim(line)),*) LASTTIME
              if (VERBOSE) print *,'lasttimestep',LASTTIME
            case('CPHAS')
              read(line(35:len_trim(line)),*) cphasetype
              cphasetype=trim(cphasetype)
              call determinePhaseRef(cphasetype,8,cphaseRef)
              if (VERBOSE) print *,'cphase    ', cphasetype, cphaseRef
            case('SOURC')
              read(line(35:len_trim(line)),*) sourceLat,sourceLon
              if (VERBOSE) print *,'source',sourceLat,sourceLon
            case('RECEI')
              read(line(35:len_trim(line)),*) receiverLat,receiverLon
              if (VERBOSE) print *,'receiver',receiverLat,receiverLon
            case('MANYR')
              read(line(35:len_trim(line)),*) manyReceivers
              if (VERBOSE) print *,'many receiver stations', manyReceivers
            case('MANYN')
              read(line(35:len_trim(line)),*) numofReceivers
              if (VERBOSE .and. manyReceivers) print *,'    receiver stations',numofReceivers
            case('MANYK')
              read(line(35:len_trim(line)),*) manyKernels
              if (VERBOSE) print *,'many kernels', manyKernels
            case('KRNEP')
              read(line(35:len_trim(line)),*) kernelStartDistance,kernelEndDistance
              if (VERBOSE .and. manyKernels) print *,'    epicentral distances',kernelStartDistance,kernelEndDistance
            case('IMPOR')
              read(line(35:len_trim(line)),*) importKernelsReceivers
              if (VERBOSE .and. manyKernels) print *,'    importing receivers', importKernelsReceivers
            case('DELTA')
              read(line(35:len_trim(line)),*)DELTA
              if (VERBOSE) print *,'delta', DELTA
            case('MOVED')
              read(line(35:len_trim(line)),*) MOVEDELTA
              if (VERBOSE .and. DELTA) print *,'  move delta', MOVEDELTA
            case('DRADI')
              read(line(35:len_trim(line)),*) DELTARADIUS
              if (VERBOSE .and. DELTA) print *,'  delta radius',DELTARADIUS
            case('DTYPE')
              read(line(35:len_trim(line)),*) DELTAfunction
              DELTAfunction=trim(DELTAfunction)
              if (VERBOSE .and. DELTA) print *,'  delta type ',DELTAfunction
            case('DLOCA')
              read(line(35:len_trim(line)),*) deltaLat,deltaLon
              if (VERBOSE .and. DELTA) print *,'  delta location',deltaLat,deltaLon
            case('DPERT')
              read(line(35:len_trim(line)),*) deltaPerturbation
              if (VERBOSE .and. DELTA) print *,'  delta perturbation',deltaPerturbation
            case('DLATS')
              read(line(35:len_trim(line)),*) latitudeStart
              if (VERBOSE .and. DELTA) print *,'  delta start latitude',latitudeStart
            case('DLATE')
              read(line(35:len_trim(line)),*) latitudeEnd
              if (VERBOSE .and. DELTA) print *,'  delta end latitude',latitudeEnd
            case('DLONE')
              read(line(35:len_trim(line)),*) longitudeEnd
              if (VERBOSE .and. DELTA) print *,'  delta end longitude',longitudeEnd
            case('DINCR')
              read(line(35:len_trim(line)),*) deltaMoveIncrement
              if (VERBOSE .and. DELTA) print *,'  delta move increment',deltaMoveIncrement
            case('SECON')
              read(line(35:len_trim(line)),*)SECONDDELTA
              if (VERBOSE .and. DELTA) print *,'  second delta:',SECONDDELTA
            case('HETER')
              read(line(35:len_trim(line)),*) HETEROGENEOUS
              if (VERBOSE) print *,'heterogeneous', HETEROGENEOUS
            case('BLKFI')
              read(line(35:len_trim(line)),*) line
              phaseBlockFile = trim(line)
              if (DO_CHECKERBOARD) then
                if (VERBOSE .and. HETEROGENEOUS) &
                  print *,'   checkerboard - L',MAP_DEGREE_L,'M',MAP_DEGREE_M
              else
                if (VERBOSE .and. HETEROGENEOUS) &
                  print *,'  phaseBlockFile used: ',trim(phaseBlockFile)
              endif
            case('BLK/I')
              read(line(35:len_trim(line)),*) heterogeneousPixelsize
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
              read(line(35:len_trim(line)),*) tmp
              tmp = trim(tmp)
              call determinePhaseRef(tmp,128,phaseBlockVelocityReference)
              if (VERBOSE .and. HETEROGENEOUS) print *,'  phaseBlock VelocityReference',phaseBlockVelocityReference
            case('INV_D')
              read(line(35:len_trim(line)),*) line
              heterogeneousDataFile = trim(line)
              if (VERBOSE .and. (HetPhaseshift_Program .or. Adjoint_InversionProgram)) &
                print *,'  Data input: ',trim(heterogeneousDataFile)
            case('INV_O')
              read(line(35:len_trim(line)),*) line
              heterogeneousOutput = trim(line)
              if (VERBOSE .and. (HetPhaseshift_Program .or. Adjoint_InversionProgram)) &
                print *,'  Output files: ',trim(heterogeneousOutput)
            case('INV_V')
              read(line(35:len_trim(line)),*) tmpInteger
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
              read(line(35:len_trim(line)),*) SIMULATIONOUTPUT
              if (VERBOSE) print *,'simulationoutput',SIMULATIONOUTPUT
            case('DATAD')
              read(line(35:len_trim(line)),*) tmp
              datadirectory = trim(tmp)
              if (VERBOSE) print *,'data output directory : ',trim(datadirectory)
            case('ADJOI')
              read(line(35:len_trim(line)),*) tmp
              adjointKernelName = trim(tmp)
              if (VERBOSE .and. Adjoint_Program) &
                print *,'adjoint kernel name : ',trim(adjointKernelName)
            case('PARAL')
              read(line(35:len_trim(line)),*)PARALLELSEISMO
              if (VERBOSE) print *,'parallelize single simulation',PARALLELSEISMO
            end select
          endif
        endif
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
