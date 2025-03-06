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
      subroutine readParameters()
!-----------------------------------------------------------------------
! read in parameters from file 'Parameter_Input'
      use propagationStartup;use loop;use phaseBlockData;use parallel; use verbosity; use adjointVariables
      implicit none
      integer:: length,ierror
      character*32:: inputName 
      character*150:: line
      character*128:: tmp
            
      ! input parameter
      inputName = 'Parameter_Input' 
      inputName = trim(inputName)
      !if(VERBOSE) print*,inputName   
      open(1, file= inputName,status='old',iostat=ierror)
      if( ierror .ne. 0) then
        print*,'error opening file ',inputName
        call stopProgram( 'abort - readParameters() opening input    ')
      endif
      
      do while( ierror .eq. 0)
        read(1,'(a150)', iostat=ierror) line
        if( ierror .ne. 0 ) exit
        
        length = len_trim(line)
        if( length .eq. 0 ) then
          continue
        else
          if(line(1:1) .eq. " " .or. line(1:1) .eq. "!" ) then 
            continue
          else
            select case( line(1:5))
            case('VERBO')
              read(line(35:len_trim(line)),*) VERBOSE
              if(VERBOSE)print*,'verbose output'         
            case('LEVEL')
              read(line(35:len_trim(line)),*) subdivisions
              if(VERBOSE)print*,'level',subdivisions
            case('FIRST')
              read(line(35:len_trim(line)),*) FIRSTTIME
              if(VERBOSE)print*,'firsttimestep',FIRSTTIME
            case('LASTT')
              read(line(35:len_trim(line)),*) LASTTIME
              if(VERBOSE)print*,'lasttimestep',LASTTIME
            case('CPHAS')
              read(line(35:len_trim(line)),*) cphasetype
              cphasetype=trim(cphasetype)
              call determinePhaseRef()
              if(VERBOSE)print*,'cphase    ', cphasetype, cphaseRef
            case('SOURC')
              read(line(35:len_trim(line)),*) sourceLat,sourceLon
              if(VERBOSE)print*,'source',sourceLat,sourceLon
            case('RECEI')
              read(line(35:len_trim(line)),*) receiverLat,receiverLon
              if(VERBOSE)print*,'receiver',receiverLat,receiverLon              
            case('MANYR')
              read(line(35:len_trim(line)),*) manyReceivers
              if(VERBOSE)print*,'many receiver stations', manyReceivers
            case('MANYN')
              read(line(35:len_trim(line)),*) numofReceivers
              if(VERBOSE)print*,'    receiver stations',numofReceivers
            case('MANYK')
              read(line(35:len_trim(line)),*) manyKernels
              if(VERBOSE)print*,'many kernels', manyKernels
            case('KRNEP')
              read(line(35:len_trim(line)),*) kernelStartDistance,kernelEndDistance
              if(VERBOSE)print*,'    epicentral distances',kernelStartDistance,kernelEndDistance
            case('IMPOR')
              read(line(35:len_trim(line)),*) importKernelsReceivers
              if(VERBOSE)print*,'    importing receivers', importKernelsReceivers
            case('DELTA')
              read(line(35:len_trim(line)),*)DELTA
              if(VERBOSE)print*,'delta', DELTA
            case('MOVED')
              read(line(35:len_trim(line)),*) MOVEDELTA
              if(VERBOSE)print*,'  move delta', MOVEDELTA
            case('DRADI')
              read(line(35:len_trim(line)),*) DELTARADIUS
              if(VERBOSE)print*,'  delta radius',DELTARADIUS
            case('DTYPE')
              read(line(35:len_trim(line)),*) DELTAFUNCTION
              DELTAFUNCTION=trim(DELTAFUNCTION)
              if(VERBOSE)print*,'  delta type ',DELTAFUNCTION
            case('DLOCA')
              read(line(35:len_trim(line)),*) deltaLat,deltaLon
              if(VERBOSE)print*,'  delta location',deltaLat,deltaLon
            case('DPERT')
              read(line(35:len_trim(line)),*) deltaPerturbation
              if(VERBOSE)print*,'  delta perturbation',deltaPerturbation
            case('SECON')
              read(line(35:len_trim(line)),*)SECONDDELTA
              if(VERBOSE)print*,'  second delta:',SECONDDELTA
            case('DLATS')
              read(line(35:len_trim(line)),*) latitudeStart
              if(VERBOSE)print*,'  delta start latitude',latitudeStart
            case('DLATE')
              read(line(35:len_trim(line)),*) latitudeEnd
              if(VERBOSE)print*,'  delta end latitude',latitudeEnd
            case('DLONE')
              read(line(35:len_trim(line)),*) longitudeEnd
              if(VERBOSE)print*,'  delta end longitude',longitudeEnd
            case('DINCR')
              read(line(35:len_trim(line)),*) deltaMoveIncrement
              if(VERBOSE)print*,'  delta move increment',deltaMoveIncrement              
            case('HETER')
              read(line(35:len_trim(line)),*) HETEROGENEOUS
              if(VERBOSE)print*,'heterogeneous', HETEROGENEOUS
            case('BLKFI')
              read(line(35:len_trim(line)),*) line
              phaseBlockFile = trim(line)              
              if(VERBOSE)print*,'  phaseBlockFile used: ', phaseBlockFile(1:len_trim(phaseBlockFile))
            case('BLKVE')
              read(line(35:len_trim(line)),*) tmp
              tmp = trim(tmp)
              if( tmp .eq. 'R40')phaseBlockVelocityReference=PHASEVELOCITY_R40
              if( tmp .eq. 'R50')phaseBlockVelocityReference=PHASEVELOCITY_R50
              if( tmp .eq. 'R60')phaseBlockVelocityReference=PHASEVELOCITY_R60
              if( tmp .eq. 'R75')phaseBlockVelocityReference=PHASEVELOCITY_R75
              if( tmp .eq. 'L35')phaseBlockVelocityReference=PHASEVELOCITY_L35
              if( tmp .eq. 'L50')phaseBlockVelocityReference=PHASEVELOCITY_L50              
              if( tmp .eq. 'L75')phaseBlockVelocityReference=PHASEVELOCITY_L75
              if( tmp .eq. 'L100')phaseBlockVelocityReference=PHASEVELOCITY_L100
              if( tmp .eq. 'L150')phaseBlockVelocityReference=PHASEVELOCITY_L150
              if( tmp .eq. 'L200')phaseBlockVelocityReference=PHASEVELOCITY_L200
              if( tmp .eq. 'L250')phaseBlockVelocityReference=PHASEVELOCITY_L250
              if( tmp .eq. 'L300')phaseBlockVelocityReference=PHASEVELOCITY_L300
              if( tmp .eq. 'L450')phaseBlockVelocityReference=PHASEVELOCITY_L450
              if( tmp .eq. 'L600')phaseBlockVelocityReference=PHASEVELOCITY_L600
              if(VERBOSE)print*,'  phaseBlock VelocityReference',phaseBlockVelocityReference
            case('INV_D')
              read(line(35:len_trim(line)),*) line
              heterogeneousDataFile = trim(line)              
              if(VERBOSE)print*,'  Data input: ',heterogeneousDataFile(1:len_trim(heterogeneousDataFile))
            case('INV_O')
              read(line(35:len_trim(line)),*) line
              heterogeneousOutput = trim(line)              
              if(VERBOSE)print*,'  Output files: ',heterogeneousOutput(1:len_trim(heterogeneousOutput))
            case('BLK/I')
              read(line(35:len_trim(line)),*) heterogeneousPixelsize
              if(VERBOSE)print*,'  Grid pixel size ',heterogeneousPixelsize
            case('SIMUL')
              read(line(35:len_trim(line)),*) SIMULATIONOUTPUT
              if(VERBOSE)print*,'simulationoutput',SIMULATIONOUTPUT
            case('DATAD')
              read(line(35:len_trim(line)),*) tmp
              datadirectory = trim(tmp)
              if(VERBOSE)print*,'data output directory : '//datadirectory(1:len_trim(datadirectory))
            case('ADJOI')
              read(line(35:len_trim(line)),*) tmp
              adjointKernelName = trim(tmp)
              if(VERBOSE)print*,'adjoint kernel name : '//adjointKernelName(1:len_trim(adjointKernelName))
            case('PARAL')
              read(line(35:len_trim(line)),*)PARALLELSEISMO
              if(VERBOSE)print*,'parallelize single simulation',PARALLELSEISMO              
            end select
          endif
        endif
      end do
    
      !if(VERBOSE) print*,'delta location(lat/lon):',deltaLat,deltaLon
      close(1)
      
      end

