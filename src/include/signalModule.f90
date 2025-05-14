!=====================================================================
!
!       m e m b r a n e S p h e r e
!       --------------------------------------------------
!      (c) 2025
!
!=====================================================================

!-----------------------------------------------------------------------
module  exception
!-----------------------------------------------------------------------
! exception handling, using signal caught in c-routines
  implicit none

  integer:: dataLoopIndex
  public:: finterrupt2, finterrupt15, setupSignalHandling

contains

  subroutine setupSignalHandling()
    implicit none
    !external def_interrupt2, def_interrupt15 !$PRAGMA C(def_interrupt2, def_interrupt15)
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
    print *,'stopping process: rank ',myrank,' by signal:',signalnumber
    print *,'    last data inverted:',dataLoopIndex
    print *
    ! wait until all processes reached this point
    !call syncProcesses()

    ! stop program
    !call stopProgram("terminate     ")
    stop
  end subroutine
end module
