##
## compilation directory
##
S := ${S_TOP}/src/propagation

################################################################
##
## compilation targets
##
################################################################

propagation_TARGETS = \
	$E/adjointMethod \
	$E/heterogeneousPhaseshift \
	$E/propagation \
	$E/timelag \
	$(EMPTY_MACRO)

propagation_OBJECTS = \
	$(OBJS_ADJOINT) \
	$(OBJS_HETSHIFT) \
	$(OBJS_PROPAGATION) \
	$(OBJS_TIMELAG) \
	$(EMPTY_MACRO)


################################################################
##
## object files
##
################################################################

OBJS_COMMON = \
	$O/cellfunctions.o \
	$O/commonModules.module.o \
	$O/drsple.o \
	$O/drspln.o \
	$O/filterRoutines.o \
	$O/forcingTerms.o \
	$O/gridfunctions.o \
	$O/initialize.o \
	$O/kernelValue.o \
	$O/numericalRecipes_correl.o \
	$O/numericalRecipes_nr.o \
	$O/numericLaplacian.o \
	$O/phasemapRoutines.o  \
	$O/phaseshiftRoutines.o \
	$O/propagationRoutines.o \
	$O/readData.o \
	$O/readParameters.o \
	$O/readPrecalculated.o \
	$O/sphericalHarmonics.o \
	$O/syncRoutines.o \
	$O/timelagRoutines.o \
	$O/vectorfunctions.o \
	$(EMPTY_MACRO)

OBJS_PROPAGATION = $(OBJS_COMMON) \
	$O/adjointRoutines.o \
	$O/propagation.o \
	$(EMPTY_MACRO)

OBJS_TIMELAG = $(OBJS_COMMON) \
	$O/timelag.o \
	$(EMPTY_MACRO)

OBJS_ADJOINT = $(OBJS_COMMON) \
	$O/adjointRoutines.o \
	$O/adjointMethod.o \
	$(EMPTY_MACRO)

OBJS_HETSHIFT = $(OBJS_COMMON) \
	$O/adjointRoutines.o \
	$O/heterogeneousPhaseshift.o \
	$(EMPTY_MACRO)

	
################################################################
##
## targets
##
################################################################

adjointMethod: $E/adjointMethod
heterogeneousPhaseshift: $E/heterogeneousPhaseshift
propagation: $E/propagation
timelag: $E/timelag


$E/propagation:$(OBJS_PROPAGATION)
	$(F90) $(LDFLAGS) -o $@ $+

$E/timelag: $(OBJS_TIMELAG)
	$(F90) $(LDFLAGS) -o $@ $+

$E/adjointMethod: $(OBJS_ADJOINT)
	$(F90) $(LDFLAGS) -o $@ $+

$E/heterogeneousPhaseshift: $(OBJS_HETSHIFT)
	$(F90) $(LDFLAGS) -o $@ $+


################################################################
##
## file compilation
##
################################################################

#.SUFFIXES: $(SUFFIXES) .f90

##
## rule for each .o file below
##

$O/%.o: src/propagation/%.f90 $O/commonModules.module.o
	$(F90) $(FFLAGS) -c -o $@ $<

