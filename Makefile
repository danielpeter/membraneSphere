#######################################
#
# membraneSphere
#
#######################################

## gfortran optimization 
#FFLAGS = -O3    # -Jobj/ -I/opt/mpisrc/include/ -ffree-line-length-none
FFLAGS = -O0 -g -ggdb -fbacktrace -fbounds-check -frange-check -Werror -Jobj/
LDFLAGS =
F90 = mpif90
F77 = mpif77 #gfortran
CFLAGS =
CC = gcc

## intel fortran optimization 
##(optional: -assume byterecl: open direct-access file assumes by default 4 byte unit, 
##                 so double precision 8 byte has recl=2; the flag resets this behaviour)
#FFLAGS = -O3 -module obj/ -I/opt/mpisrc/include/ -assume byterecl
#LDFLAGS =
#F90 = mpif90
#F77 = mpif77
#CFLAGS =
#CC = cc

## brutus optimization by using openmpi (be careful that you should first load
## the correct module, otherwise it will not work, see infos about how modules
## work on brutus at http://www.brutus.ethz.ch/wiki/index.php/Working_with_modules
#FFLAGS = -O3 -module obj/ -I/cluster/apps/openmpi/1.2.6/x86_64/intsrc/include -assume byterecl
#LDFLAGS= -lm
#F90 = mpif90
#F77 = mpif77
#CC = cc

## sun compilers cluster-one.earth.ox.ac.uk
#FFLAGS = -O3 -Mobj/ -e
#LDFLAGS = -L/opt/SUNWhpc/HPC6.0/lib -lmpi -lfmpi -lsocket -lnsl -laio
#F90 = mpf90
#F77 = mpf77
#CC = cc

## optimization for mac
#FFLAGS =  -Isrc/include -I/opt/mpich/src/fortrsrc/include -qsuppress=cmpmsg -qmoddir=obj/ -Iobj/
#LDFLAGS =
#F90 = mpif90
#CC = cc

#ifeq ($(findstring Linux,$(shell uname -a)),Linux)
#  FFLAGS = -O3 -axW -Iobj/ -qmoddir=obj/
#endif

## portland fortran optimization (fastsse same as axW on intel) for gonzales cluster 
#FFLAGS = -Mnobounds -Mneginfo -Knoieee -fast -fastsse \
#         -Isrc/include -I/usr/lib/msrc/include  -module obj/ -Iobj/ #-O3 -fastsse -Mextend
#LDFLAGS = -L/usr/lib/mpi/lib -lmpi -lmpifarg -lg2c  #-static-libcxa
#F90 = /opt/pgi/linux86-64/5.2/bin/pgf90
#CC = cc


#######################################
##
## directories
##
#######################################

## compilation directories
O = ./obj
E = ./bin
# S_TOP : source file root directory
S_TOP = .

#######################################
##
## targets
##
#######################################

# code subdirectories
SUBDIRS = \
	include \
	propagation \
	$(EMPTY_MACRO)

# default targets for the pure Fortran version
DEFAULT = \
	adjointMethod \
	propagation \
	timelag \
	$(EMPTY_MACRO)

default: $(DEFAULT)

all: default

clean:
	@echo "cleaning"
	-rm -f $(foreach dir, $(SUBDIRS), $($(dir)_OBJECTS) $($(dir)_MODULES) $($(dir)_TARGETS))
	-rm -rf $E/* $O/*

help:
	@echo "usage: make [executable]"
	@echo ""
	@echo "supported main executables:"
	@echo "    adjointMethod"
	@echo "    propagation"
	@echo "    timelag"
	@echo ""

.PHONY: all default clean help

#######################################

# Get dependencies and rules for building stuff
include $(patsubst %, ${S_TOP}/src/%/Makefile, $(SUBDIRS))

#######################################

##
## Shortcuts
##

# Shortcut for: <prog>/<prog> -> bin/<prog>
define target_shortcut
$(patsubst $E/%, %, $(1)): $(1)
.PHONY: $(patsubst $E/%, %, $(1))
endef

# Shortcut for: dir -> src/dir/<targets in here>
define shortcut
$(1): $($(1)_TARGETS)
.PHONY: $(1)
$$(foreach target, $$(filter $E/%,$$($(1)_TARGETS)), $$(eval $$(call target_shortcut,$$(target))))
endef

$(foreach dir, $(SUBDIRS), $(eval $(call shortcut,$(dir))))
