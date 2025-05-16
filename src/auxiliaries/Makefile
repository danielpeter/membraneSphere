################################################################
# Makefile
DIR = auxiliaries

################################################################
##
## compilation targets
##
################################################################

auxiliaries_TARGETS = \
	$E/sphericalLaplacian \
	$(EMPTY_MACRO)

auxiliaries_OBJECTS = \
	$(OBJS_SPHERICALLAPLACIAN) \
	$(EMPTY_MACRO)


################################################################
##
## object files
##
################################################################

OBJS_AUXILIARIESCOMMON = \
	$O/analyticSolutions.o \
	$O/cellfunctions.o \
	$O/commonModules.module.o \
	$O/correlationRoutines.o \
	$O/filterRoutines.o \
	$O/gridfunctions.o \
	$O/numericLaplacian.o \
	$O/readData.o \
	$O/readParameters.o \
	$O/readPrecalculated.o \
	$O/sphericalHarmonics.o \
	$O/splineRoutines.o \
	$O/timelagRoutines.o \
	$O/vectorfunctions.o \
	$(EMPTY_MACRO)

OBJS_SPHERICALLAPLACIAN = $(OBJS_AUXILIARIESCOMMON) \
	$O/sphericalLaplacian.o \
	$(EMPTY_MACRO)


################################################################
##
## targets
##
################################################################

sphericalLaplacian: $E/sphericalLaplacian

auxiliaries: $(auxiliaries_TARGETS)

$E/sphericalLaplacian: $(OBJS_SPHERICALLAPLACIAN)
	@echo ""
	@echo "building sphericalLaplacian"
	@echo ""
	$(F90) $(LDFLAGS) -o $@ $+


.PHONY: auxiliaries

################################################################
##
## file compilation
##
################################################################

##
## compilation directory
##
S := ${S_TOP}/src/auxiliaries
$(auxiliaries_OBJECTS): S := ${S_TOP}/src/auxiliaries

#.SUFFIXES: $(SUFFIXES) .f90

##
## rule for each .o file below
##

$O/%.o: $S/%.f90 $O/commonModules.module.o
	$(F90) $(FFLAGS) -c -o $@ $<

