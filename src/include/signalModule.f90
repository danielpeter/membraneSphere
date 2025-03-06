!=====================================================================
!
!       m e m b r a n e S p h e r e  1 . 1
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
      module  exception
      ! exception handling, using signal caught in c-routines
        integer:: dataLoopIndex     
        public:: finterrupt2, finterrupt15,setupSignalHandling         
      contains
        subroutine setupSignalHandling()  
          implicit none
          external def_interrupt2, def_interrupt15 !$PRAGMA C(def_interrupt2, def_interrupt15)              
          call def_interrupt2()
          ! call def_interrupt11() 
          call def_interrupt15() 
        end subroutine
      
        subroutine finterrupt2() 
          call terminate(2)
        end subroutine

        subroutine finterrupt15() 
          call terminate(15)
        end subroutine
      
        subroutine terminate(signalnumber)
          use parallel
          implicit none
          integer:: i,signalnumber
          ! console output
          print*,'stopping process:',rank,' by signal:',signalnumber
          print*,'    last data inverted:',dataLoopIndex
          print*
          ! wait until all processes reached this point
          !call MPI_Barrier( MPI_COMM_WORLD, ierror )
          !if( ierror .ne. 0) call stopProgram('terminate MPI_Barrier failed    ')      
          
          ! stop program
          !call stopProgram("terminate     ")
          stop
        end subroutine      
      end module      
