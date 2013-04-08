program = diffusion

templates = $(wildcard templates/*.o)
xml_fort = xml-fortran/xmlparse.o \
           xml-fortran/read_xml_primitives.o \
           xml-fortran/write_xml_primitives.o
slatec = $(wildcard slatec/*.o)

#===============================================================================
# Object Files
#===============================================================================

include OBJECTS

#===============================================================================
# User Options
#===============================================================================

COMPILER = gnu
DEBUG    = no
PROFILE  = no
OPTIMIZE = no 

#===============================================================================
# Add git SHA-1 hash
#===============================================================================

GIT_SHA1 = $(shell git log -1 | head -n 1 | awk '{print $$2}')

#===============================================================================
# GNU Fortran compiler options
#===============================================================================

ifeq ($(COMPILER),gnu)
  F90 = gfortran
  F90FLAGS := -cpp -fbacktrace -DNO_F2008
  LDFLAGS =

  # Debugging
  ifeq ($(DEBUG),yes)
    F90FLAGS += -g -Wall -pedantic -std=f2008 -fbounds-check \
                -ffpe-trap=invalid,zero,overflow,underflow
    LDFLAGS  += -g
  endif

  # Profiling
  ifeq ($(PROFILE),yes)
    F90FLAGS += -pg
    LDFLAGS  += -pg
  endif

  # Optimization
  ifeq ($(OPTIMIZE),yes)
    F90FLAGS += -O3
  endif
endif

#===============================================================================
# Slatec
#===============================================================================

LDFLAGS += -L slatec -lslatec

#===============================================================================
# Targets
#===============================================================================

all: xml-fortran slatec $(program)
xml-fortran:
	cd xml-fortran; make F90=$(F90) F90FLAGS="$(F90FLAGS)"
	cd templates; make F90=$(F90) F90FLAGS="$(F90FLAGS)"
slatec:
	cd slatec; make MACHINE=$(MACHINE) F90=$(F90) F90FLAGS="$(F90FLAGS)"
$(program): $(objects)
	$(F90) $(objects) $(templates) $(xml_fort) -o $@ $(LDFLAGS)
distclean: clean
	cd xml-fortran; make clean
	cd templates; make clean
	cd slatec; make clean
clean:
	@rm -f *.o *.mod $(program)
neat:
	@rm -f *.o *.mod

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .F90 .o
.PHONY: all xml-fortran clean neat distclean slatec 

%.o: %.F90
	$(F90) $(F90FLAGS) -DGIT_SHA1="\"$(GIT_SHA1)\"" -Ixml-fortran -Itemplates -Islatec -c $<

#===============================================================================
# Dependencies
#===============================================================================

include DEPENDENCIES
