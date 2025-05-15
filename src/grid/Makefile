################################################################
# Makefile
DIR = grid

################################################################
##
## compilation targets
##
################################################################

grid_TARGETS = \
	$E/gridConstruction \
	$E/gridAreaStatistics \
	$(EMPTY_MACRO)

grid_OBJECTS = \
	$(OBJS_GRIDCONSTRUCTION) \
	$(OBJS_GRIDAREASTATISTICS) \
	$(EMPTY_MACRO)


################################################################
##
## object files
##
################################################################

OBJS_GRIDCOMMON = \
	$O/cellfunctions.o \
	$O/commonModules.module.o \
	$O/gridfunctions.o \
	$O/readData.o \
	$O/vectorfunctions.o \
	$(EMPTY_MACRO)

OBJS_GRIDCONSTRUCTION = \
	$O/gridConstruction.o \
	$(EMPTY_MACRO)

OBJS_GRIDAREASTATISTICS = $(OBJS_GRIDCOMMON) \
	$O/gridAreaStatistics.o \
	$(EMPTY_MACRO)
	
################################################################
##
## targets
##
################################################################

gridConstruction: $E/gridConstruction
gridAreaStatistics: $E/gridAreaStatistics

grid: $(grid_TARGETS)

$E/gridConstruction: $(OBJS_GRIDCONSTRUCTION)
	@echo ""
	@echo "building gridConstruction"
	@echo ""
	$(F90) $(LDFLAGS) -o $@ $+

$E/gridAreaStatistics: $(OBJS_GRIDAREASTATISTICS)
	@echo ""
	@echo "building gridAreaStatistics"
	@echo ""
	$(F90) $(LDFLAGS) -o $@ $+

.PHONY: grid

################################################################
##
## file compilation
##
################################################################

##
## compilation directory
##
S := ${S_TOP}/src/grid
$(grid_OBJECTS): S := ${S_TOP}/src/grid

#.SUFFIXES: $(SUFFIXES) .f90

##
## rule for each .o file below
##

$O/%.o: $S/%.f90 $O/commonModules.module.o
	$(F90) $(FFLAGS) -c -o $@ $<

