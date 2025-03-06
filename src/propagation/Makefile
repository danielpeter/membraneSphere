################################################################
# compilation flags
################################################################
## compilation directories
MOD = mod/
LOBJ = obj/

## intel fortran optimization 
##(optional: -assume byterecl: open direct-access file assumes by default 4 byte unit, 
##                 so double precision 8 byte has recl=2; the flag resets this behaviour)
FFLAGS = -O3 -module $(MOD) -I/opt/mpich/include/ -assume byterecl
LDFLAGS =
F90 = mpif90
F77 = mpif77

## gfortran optimization 
#FFLAGS = -O3 -J$(MOD) -I/opt/mpich/include/ -ffree-line-length-none
#LDFLAGS =
#F90 = mpif90

## brutus optimization by using openmpi (be careful that you should first load
## the correct module, otherwise it will not work, see infos about how modules
## work on brutus at http://www.brutus.ethz.ch/wiki/index.php/Working_with_modules
#FFLAGS = -O3 -module $(MOD) -I/cluster/apps/openmpi/1.2.6/x86_64/intel/include -assume byterecl 
#LDFLAGS= -lm
#F90 = mpif90
#F77 = mpif77

## sun compilers cluster-one.earth.ox.ac.uk
#FFLAGS = -O3 -M$(MOD) -e
#LDFLAGS = -L/opt/SUNWhpc/HPC6.0/lib -lmpi -lfmpi -lsocket -lnsl -laio
#F90 = mpf90
#F77 = mpf77

## optimization for mac
#FFLAGS =  -I../include -I/opt/mpich/src/fortran/include -qsuppress=cmpmsg -qmoddir=mod -Imod 
#LDFLAGS = 
#F90 = mpif90

#ifeq ($(findstring Linux,$(shell uname -a)),Linux)
#FFLAGS = -O3 -axW -Imod -qmoddir=mod   
#endif

## portland fortran optimization (fastsse same as axW on intel) for gonzales cluster 
#FFLAGS = -Mnobounds -Mneginfo -Knoieee -fast -fastsse \
#			-I../include -I/usr/lib/mpi/include  -module mod -Imod #-O3 -fastsse -Mextend
#LDFLAGS = -L/usr/lib/mpi/lib -lmpi -lmpifarg -lg2c  #-static-libcxa
#F90 = /opt/pgi/linux86-64/5.2/bin/pgf90

################################################################
# compilation targets
################################################################

all: propagation timelag adjointMethod heterogeneousPhaseshift 

################################################################
# source files
################################################################

SRC = ../include/commonModules.f90 ../include/initialize.f90 ../include/phaseshiftRoutines.f90 \
	../include/propagationRoutines.f90 ../include/gridfunctions.f90 readParameters.f90 \
	../include/kernelValue.f90 ../include/numericalRecipes_nr.f90 ../include/numericalRecipes_correl.f90 \
	../include/filterRoutines.f90 ../include/timelagRoutines.f90 	../include/readData.f90 \
	../include/vectorfunctions.f90 ../include/readPhaseMap.f90 ../include/syncRoutines.f90 \
	../include/forcingTerms.f90 ../include/readPrecalculated.f90 ../include/numericLaplacian.f90 \
	../include/cellfunctions.f90 ../include/phasemapRoutines.f90 phasedata.f90 heterogeneousInversion.f90  \
	../include/signalhandler.c ../include/signalModule.f90 ../include/sphericalHarmonics.f90 \
	heterogeneousPhaseshift.f90 phaseshift.f90 ../include/drsple.f ../include/drspln.f

################################################################
# object files
################################################################

OBJS_PROPAGATION = $(LOBJ)commonModules.o \
	$(LOBJ)numericalRecipes_nr.o $(LOBJ)numericalRecipes_correl.o \
	$(LOBJ)initialize.o $(LOBJ)propagationRoutines.o $(LOBJ)timelagRoutines.o \
	$(LOBJ)adjointRoutines.o $(LOBJ)gridfunctions.o $(LOBJ)phasemapRoutines.o \
	$(LOBJ)filterRoutines.o $(LOBJ)vectorfunctions.o $(LOBJ)forcingTerms.o \
	$(LOBJ)readParameters.o $(LOBJ)syncRoutines.o $(LOBJ)readData.o \
	$(LOBJ)readPrecalculated.o $(LOBJ)numericLaplacian.o $(LOBJ)cellfunctions.o \
	$(LOBJ)kernelValue.o $(LOBJ)sphericalHarmonics.o $(LOBJ)drsple.o $(LOBJ)drspln.o $(LOBJ)propagation.o 

OBJS_TIMELAG = $(LOBJ)commonModules.o \
	$(LOBJ)numericalRecipes_nr.o $(LOBJ)numericalRecipes_correl.o \
	$(LOBJ)filterRoutines.o $(LOBJ)timelagRoutines.o $(LOBJ)vectorfunctions.o \
	$(LOBJ)readData.o $(LOBJ)readPrecalculated.o $(LOBJ)readParameters.o $(LOBJ)drsple.o $(LOBJ)drspln.o $(LOBJ)timelag.o 

OBJS_COMMON = $(LOBJ)commonModules.o $(LOBJ)initialize.o \
	$(LOBJ)phaseshiftRoutines.o $(LOBJ)propagationRoutines.o \
	$(LOBJ)gridfunctions.o $(LOBJ)readParameters.o \
	$(LOBJ)kernelValue.o $(LOBJ)numericalRecipes_nr.o \
	$(LOBJ)timelagRoutines.o $(LOBJ)filterRoutines.o \
	$(LOBJ)numericalRecipes_correl.o $(LOBJ)readData.o \
	$(LOBJ)vectorfunctions.o \
	$(LOBJ)syncRoutines.o $(LOBJ)forcingTerms.o \
	$(LOBJ)readPrecalculated.o $(LOBJ)numericLaplacian.o \
	$(LOBJ)cellfunctions.o $(LOBJ)phasemapRoutines.o  \
	$(LOBJ)sphericalHarmonics.o $(LOBJ)drsple.o $(LOBJ)drspln.o  

OBJS_ADJOINT = $(OBJS_COMMON) $(LOBJ)adjointRoutines.o $(LOBJ)adjointMethod.o 

OBJS_HETSHIFT = $(OBJS_COMMON)  \
	$(LOBJ)adjointRoutines.o $(LOBJ)heterogeneousPhaseshift.o   

	
################################################################
# targets
################################################################

propagation:$(OBJS_PROPAGATION)
	$(F90) $(LDFLAGS) $(OBJS_PROPAGATION) -o propagation

timelag: $(OBJS_TIMELAG)
	$(F90) $(LDFLAGS) $(OBJS_TIMELAG) -o timelag

adjointMethod: $(OBJS_ADJOINT)
	$(F90) $(LDFLAGS) $(OBJS_ADJOINT) -o adjointMethod

heterogeneousPhaseshift: $(OBJS_HETSHIFT)
	$(F90) $(LDFLAGS) $(OBJS_HETSHIFT) -o heterogeneousPhaseshift	

	
clean:
	rm -f propagation timelag adjointMethod heterogeneousPhaseshift \
	$(LOBJ)/*.o *.o *.mod $(MOD)/*.mod 


################################################################
# file compilation
################################################################

.SUFFIXES: $(SUFFIXES) .f90

$(LOBJ)%.o: %.f
	$(F77) $(FFLAGS) -c $*.f
	mv -f $*.o $(LOBJ)

$(LOBJ)%.o: ../include/%.f
	$(F77) $(FFLAGS) -c ../include/$*.f
	mv -f $*.o $(LOBJ)

$(LOBJ)%.o: %.f90 	
	$(F90) $(FFLAGS) -c $*.f90
	mv -f $*.o $(LOBJ)

$(LOBJ)%.o: ../include/%.f90 	
	$(F90) $(FFLAGS) -c ../include/$*.f90
	mv -f $*.o $(LOBJ)

$(LOBJ)%.o: ../include/%.c 	
	cc -c ../include/$*.c
	mv -f $*.o $(LOBJ)

