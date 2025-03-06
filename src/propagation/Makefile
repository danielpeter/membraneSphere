#############################################################################################
# compilation flags
#############################################################################################

# compilation directories
MOD = mod/
LOBJ = obj/


# intel fortran optimization 
#(optional: -assume byterecl: open direct-access file assumes by default 4 byte unit, 
#                 so double precision 8 byte has recl=2; the flag resets this behaviour)
FFLAGS = -O3 -module $(MOD) -I/opt/mpich/include/
LDFLAGS =
F90 = mpif90

# portland fortran optimization (fastsse same as axW on intel) for gonzales cluster 
#FFLAGS = -Mnobounds -Mneginfo -Knoieee -fast -fastsse -I../include -I/usr/lib/mpi/include  -qmoddir=mod -Imod 
#LDFLAGS = -L/usr/lib/mpi/lib -lmpi -lmpifarg -lg2c 
#F90 = /opt/pgi/linux86-64/5.2/bin/pgf90

# optimization for mac
#FFLAGS =  -I../include -I/opt/mpich/src/fortran/include -qsuppress=cmpmsg -qmoddir=mod -Imod 
#LDFLAGS = 
#F90 = mpif90

#############################################################################################
# compilation targets
#############################################################################################

all: propagation timelag adjointMethod  

#############################################################################################
# source files
#############################################################################################

SRC = ../include/commonModules.f90 ../include/initialize.f90 ../include/phaseshiftRoutines.f90 \
	../include/propagationRoutines.f90 ../include/gridfunctions.f90 readParameters.f90 \
	../include/kernelValue.f90 ../include/numericalRecipes_nr.f90 ../include/numericalRecipes_correl.f90 \
	../include/filterRoutines.f90 ../include/timelagRoutines.f90 	../include/readData.f90 \
	../include/vectorfunctions.f90 ../include/readPhaseMap.f90 ../include/syncRoutines.f90 \
	../include/forcingTerms.f90 ../include/readPrecalculated.f90 ../include/numericLaplacian.f90 \
	../include/cellfunctions.f90 ../include/phasemapRoutines.f90 phasedata.f90 heterogeneousInversion.f90  

#############################################################################################
# object files
#############################################################################################

OBJS_PROPAGATION = $(LOBJ)commonModules.o $(LOBJ)numericalRecipes_nr.o $(LOBJ)numericalRecipes_correl.o \
	$(LOBJ)initialize.o $(LOBJ)propagationRoutines.o $(LOBJ)timelagRoutines.o \
	$(LOBJ)adjointRoutines.o $(LOBJ)gridfunctions.o $(LOBJ)phasemapRoutines.o \
	$(LOBJ)filterRoutines.o $(LOBJ)vectorfunctions.o $(LOBJ)forcingTerms.o \
	$(LOBJ)readParameters.o $(LOBJ)syncRoutines.o $(LOBJ)readData.o \
	$(LOBJ)readPrecalculated.o $(LOBJ)numericLaplacian.o $(LOBJ)cellfunctions.o \
	$(LOBJ)propagation.o 

OBJS_TIMELAG = $(LOBJ)commonModules.o $(LOBJ)numericalRecipes_nr.o $(LOBJ)numericalRecipes_correl.o \
	$(LOBJ)filterRoutines.o $(LOBJ)timelagRoutines.o $(LOBJ)vectorfunctions.o \
	$(LOBJ)readData.o $(LOBJ)readPrecalculated.o $(LOBJ)timelag.o 

OBJS_COMMON = $(LOBJ)commonModules.o $(LOBJ)initialize.o $(LOBJ)phaseshiftRoutines.o $(LOBJ)propagationRoutines.o \
	$(LOBJ)gridfunctions.o $(LOBJ)readParameters.o $(LOBJ)kernelValue.o $(LOBJ)numericalRecipes_nr.o \
	$(LOBJ)timelagRoutines.o $(LOBJ)filterRoutines.o $(LOBJ)numericalRecipes_correl.o $(LOBJ)readData.o \
	$(LOBJ)vectorfunctions.o $(LOBJ)readPhaseMap.o $(LOBJ)syncRoutines.o $(LOBJ)forcingTerms.o \
	$(LOBJ)readPrecalculated.o $(LOBJ)numericLaplacian.o $(LOBJ)cellfunctions.o $(LOBJ)phasemapRoutines.o 

OBJS_ADJOINT = $(OBJS_COMMON) $(LOBJ)adjointRoutines.o $(LOBJ)adjointMethod.o 
	
#############################################################################################
# targets
#############################################################################################

propagation:$(OBJS_PROPAGATION)
	$(F90) $(LDFLAGS) $(OBJS_PROPAGATION) -o propagation

timelag: $(OBJS_TIMELAG)
	$(F90) $(LDFLAGS) $(OBJS_TIMELAG) -o timelag

adjointMethod: $(OBJS_ADJOINT)
	$(F90) $(LDFLAGS) $(OBJS_ADJOINT) -o adjointMethod
	
clean:
	rm -f propagation timelag adjointMethod $(LOBJ)/*.o *.o *.mod $(MOD)/*.mod


#############################################################################################
# file compilation
#############################################################################################

.SUFFIXES: $(SUFFIXES) .f90

$(LOBJ)%.o: %.f90 	
	$(F90) $(FFLAGS) -c $*.f90
	mv -f $*.o $(LOBJ)

$(LOBJ)%.o: ../include/%.f90 	
	$(F90) $(FFLAGS) -c ../include/$*.f90
	mv -f $*.o $(LOBJ)

$(LOBJ)%.o: ../include/%.c 	
	cc -c ../include/$*.c
	mv -f $*.o $(LOBJ)
	