dnl                                               -*- Autoconf -*-
dnl
dnl Process this file with autoconf to produce the 'configure' script:
dnl
dnl autoreconf -i
dnl
dnl You must have recent versions of Autoconf and Automake installed.

############################################################

AC_PREREQ([2.72])
AC_INIT([membraneSphere],[m4_normalize(m4_include(VERSION))],[see the doc],[membraneSphere])
AC_CONFIG_SRCDIR([src/include/commonModules.f90])

############################################################

AS_BOX([setting up compilation flags])
AC_CANONICAL_HOST

# 'configure' options

###
### precision
###

AC_ARG_ENABLE([double-precision],
    [AS_HELP_STRING([--enable-double-precision],
        [solver in double precision @<:@default=no@:>@])],
    [want_double_precision="$enableval"],
    [want_double_precision=no])

AS_IF([test x"$want_double_precision" = xno],[
    WP=4
],[
    WP=8
])
AC_SUBST([WP])


###
### debugging
###

AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
        [build with debugging options enabled @<:@default=no@:>@])],
    [want_debug="$enableval"],
    [want_debug=no])


############################################################

# Checks for programs.

# a courtesy to the installed base of users
if test x"$FC" = x && test x"$F90" != x; then
    FC="$F90"
fi
if test x"$MPIFC" = x && test x"$MPIF90" != x; then
    MPIFC="$MPIF90"
fi
# note: if not set, default FCFLAGS (and CFLAGS) will be set by autoconf/configure script (mostly -g -O2) by AC_PROG_FC
if test x"$FCFLAGS" = x ; then
    FCFLAGS=""
fi
FLAGS_CHECK="$FCFLAGS"

AC_PROG_FC
export FC
export MPIFC
FFLAGS="$FCFLAGS"

###
### compiler flags / debugging
###

# First find the "real" compiler
# e.g., Cray provides a Programming Environment that wraps around compilers,
# but the command-line options are still compiler-specific.
case $FC in
    ftn|*/ftn)
        case $PE_ENV in
            CRAY)
                my_FC=crayftn
                ;;
            GNU)
                my_FC=gfortran
                ;;
            INTEL)
                my_FC=ifort
                ;;
            *)
                # Unrecognized (this will not set any default flags)
                my_FC="$FC"
                ;;
        esac
        ;;
    *)
        my_FC="$FC"
        ;;
esac
case $my_FC in
    ftn|*/ftn|crayftn|*/crayftn)
        # Cray Fortran
        OPT_FFLAGS="-O3 -Onoaggress -Oipa0 -hfp2 -Ovector3 -Oscalar3 -Ocache2 -Ounroll2 -Ofusion2"
        DEBUG_FFLAGS="-eC -eD -ec -en -eI -ea -g -G0"
        ;;
    ifort|*/ifort)
        # Intel ifort Fortran90 for Linux
        OPT_FFLAGS="-O3 -check nobounds -assume byterecl"
        DEBUG_FFLAGS="-check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv -gen-interfaces -warn all,noexternal"
        ;;
    gfortran|*/gfortran|f95|*/f95)
        # GNU gfortran
        OPT_FFLAGS="-O3"
        DEBUG_FFLAGS="-O0 -g -ggdb -fbacktrace -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow -std=f2008 -Wall -Wno-conversion"
        ;;
esac

# note: if not set, default FCFLAGS (and CFLAGS) would be set by autoconf/configure script (mostly -g -O2).
#       we try to guess the flags based on the above compiler settings.
AS_IF([test x"$FLAGS_CHECK" = x && test x"${OPT_FFLAGS}" != x],[
    AS_IF([test x"$want_debug" = xno],[
        # optimization flags
        AC_MSG_NOTICE([guessed optimization flags ${OPT_FFLAGS}])
        FCFLAGS="${OPT_FFLAGS}"
    ],[
        # debug flags
        AC_MSG_NOTICE([guessed debugging flags ${DEBUG_FFLAGS}])
        FCFLAGS="${DEBUG_FFLAGS}"
    ])
],[
    AC_MSG_NOTICE([with flags: $FLAGS_CHECK])
])

AC_SUBST(srcdir)

# Fortran, use .f90 as default extension
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)
AC_FC_FREEFORM()

# check for C compiler
AC_PROG_CC

############################################################
# influential environment variables

AC_ARG_VAR(MPIFC, [MPI Fortran compiler command])
AC_ARG_VAR(MPILIBS, [extra libraries for linking MPI programs])

#AC_ARG_VAR(FCFLAGS, [Fortran compiler flags])

############################################################

# Checks for typedefs, structures, and compiler characteristics.
# check fortran modules flag
AS_BOX([module extensions])

AC_FC_MODULE_EXTENSION
AC_SUBST([FC_MODEXT])
AC_FC_MODULE_FLAG([], [])
AC_FC_MODULE_OUTPUT_FLAG([
  FCFLAGS_f90="$FC_MODOUT./obj $FC_MODINC./obj $FC_MODINC. $FCFLAGS_f90"
  FC_MODDIR=./obj
], [
  FC_MODDIR=.
])
AC_SUBST([FC_MODDIR])


###
### MPI
###

AS_BOX([MPI])
AC_MSG_NOTICE([MPI is needed])

# tests
AS_IF([test x"$MPIFC" = x],[
    MPIFC=mpif90
])

# checking for Fortran
AC_LANG_PUSH(Fortran)
AC_MSG_CHECKING([for MPI modules])
FCFLAGS="$FCFLAGS"
# adds include directory
if test "x$MPI_INC" != "x"; then
  MPI_INCLUDES="-I$MPI_INC"
  FCFLAGS+=" $MPI_INCLUDES"
fi
# compilation test
ac_compile='$MPIFC -c $FCFLAGS conftest.$ac_ext >&5'
AC_COMPILE_IFELSE([
  AC_LANG_PROGRAM([], [[
  use mpi
  ]])
], [
  AC_MSG_RESULT(yes)
], [
  AC_MSG_RESULT(no)
  AC_MSG_ERROR([MPI module not found; is MPI built with Fortran support for this compiler?])
])
AC_LANG_POP([Fortran])
AC_SUBST([MPI_INCLUDES])

############################################################

AS_BOX([setting up default simulation setup])

# Output results.
AC_CONFIG_FILES([
    Makefile
    src/include/constants.h
])

AC_CONFIG_COMMANDS([bin], [AS_MKDIR_P(bin)])
AC_CONFIG_COMMANDS([obj], [AS_MKDIR_P(obj)])

AC_OUTPUT


dnl FYI, this is not a "result", but AS_ECHO is not available in 2.61, and it
dnl basically works like a portable echo.

AC_MSG_RESULT([])
AS_BOX([$PACKAGE_NAME $PACKAGE_VERSION])
AC_MSG_RESULT([])
AC_MSG_RESULT([$0 has completed and set up a default configuration to build.])
AC_MSG_RESULT([])
AC_MSG_RESULT([You may wish to modify the following files before running a simulation:])
AC_MSG_RESULT([AS_HELP_STRING([Parameter_Input],
                              [Set parameters affecting the simulation.])])
AC_MSG_RESULT([AS_HELP_STRING([Timelag_Input],
                              [Set parameters for timelag executable.])])
AC_MSG_RESULT([])

dnl end of configure.ac
