# Makefile to compile examples
# ----------------------------
# J.-M. Brankart, November 2018
# 
include ../build/Makefile.macro
#
# Add EnsDAM library to compile examples:
LIB_ENSDAM = ../lib/libensdam.a
LIB = -L../lib -lensdam $(LDFLAGS)

# Define source and executables
src1 = $(wildcard *.f90)
exe1 = $(src1:.f90=.x)

src2 = $(wildcard *.F90)
exe2 = $(src2:.F90=.x)


# Compilation rules
%.x : %.f90 $(LIB_ENSDAM)
	$(FC) $(FPPFLAGS) $(FFLAGS) -o $@ $< $(LIB)

%.x : %.F90 $(LIB_ENSDAM)
	$(FC) $(FPPFLAGS) $(FFLAGS) -o $@ $< $(LIB)

examples: $(exe1) $(exe2)

clean:
	rm -f *.x
