program = libslatec.a
source  = $(wildcard *.f)
objects = $(source:.f=.o)

#===============================================================================
# Targets
#===============================================================================

RANLIB = ranlib

all: $(program)

$(program): $(objects)
	ar ruv $(program) $(objects)
	$(RANLIB) $(program)

clean:
	@rm -f *.o $(program)

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .f .o
.PHONY: all clean

%.o: %.f
	$(F90) $(F90FLAGS) -c $<
