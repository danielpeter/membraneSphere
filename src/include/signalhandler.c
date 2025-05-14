/*
!=====================================================================
!
!       m e m b r a n e S p h e r e
!       --------------------------------------------------
!      (c) 2025
!
!=====================================================================
*/

/* signal handling (easiest within c-codes)

when jobs are in a queue and the walltime limit is reached, they receive a warning signal
to exit gracefully.
here, we catch the SIGINT interrupt signal and call a fortran routine to stop the program

*/
#include <signal.h>
#include <stdio.h>


#pragma linkage(finterrupt2,FORTRAN)
#pragma linkage(finterrupt15,FORTRAN)

#define finterrupt2 F77_NAME(finterrupt2,finterrupt2)
#define finterrupt15 F77_NAME(finterrupt15,finterrupt15)

#define def_interrupt2 C_NAME(def_interrupt2,def_interrupt2)
#define def_interrupt15 C_NAME(def_interrupt15,def_interrupt15)


/* Portland */
//#define F77_NAME(name,NAME) __exception_MOD_##name
//#define C_NAME(name,NAME) name

/* intel */
#define F77_NAME(name,NAME) exception_mp_##name##_
#define C_NAME(name,NAME) name##_

/* others */
//#define F77_NAME(name,NAME) name##_

// fortran routines which will be called when signal is caught
//extern void finterrupt11_();
//extern void exception_mp_finterrupt2_();
//extern void exception_mp_finterrupt15_();
extern void finterrupt2();
extern void finterrupt15();

// set new SIGINT interrupt routine
void def_interrupt2()
{
  //printf("setting new interrupt handler \n");
  sigset(2,finterrupt2);
}
// SIGTERM termination
void def_interrupt15()
{
  sigset(15,finterrupt15);
}
/*
// SIGSEGV segmentation fault
void def_interrupt11_()
{
  sigset(11,finterrupt11_);
}
*/



