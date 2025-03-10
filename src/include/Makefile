##
## compilation directory
##
S := ${S_TOP}/src/include

################################################################

##
## object files
##

include_OBJECTS = \
	$O/adjointRoutines.o \
	$O/analyticSolutions.o \
	$O/cellfunctions.o \
	$O/commonModules.module.o \
	$O/correlationRoutines.module.o \
	$O/filterRoutines.o \
	$O/forcingTerms.o \
	$O/gridfunctions.o \
	$O/initialize.o \
	$O/kernelValue.o \
	$O/numericalRecipes_correl.o \
	$O/numericalRecipes_nr.o \
	$O/numericalRecipes_splines.o \
	$O/numericLaplacian.o \
	$O/phasemapRoutines.o  \
	$O/phaseshiftRoutines.o \
	$O/propagationRoutines.o \
	$O/readData.o \
	$O/readParameters.o \
	$O/readPrecalculated.o \
	$O/signalModule.o \
	$O/sphericalHarmonics.o \
	$O/splineRoutines.o \
	$O/syncRoutines.o \
	$O/timelagRoutines.o \
	$O/vectorfunctions.o \
	$(EMPTY_MACRO)



################################################################
##
## file compilation
##
################################################################

## dependencies
$O/numericalRecipes_correl.o: $O/numericalRecipes_nr.o
$O/numericalRecipes_splines.o: $O/numericalRecipes_nr.o
$O/filterRoutines.o: $O/numericalRecipes_nr.o
$O/timelagRoutines.o: $O/numericalRecipes_nr.o

#.SUFFIXES: $(SUFFIXES) .f90

##
## rule for each .o file below
##

$O/%.module.o: src/include/%.f90
	$(F90) $(FFLAGS) -c -o $@ $<

$O/%.o: src/include/%.f
	$(F77) $(FFLAGS) -c -o $@ $<

$O/%.o: src/include/%.f90 $O/commonModules.module.o
	$(F90) $(FFLAGS) -c -o $@ $<

$O/%.o: src/include/%.c
	$(CC) $(CFLAGS) -c -o $@ $<


