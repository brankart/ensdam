#!/bin/bash

# location of ensdam fortran library
libf="../build/src/libensdam.a"

# location of ensdam library with C-callable wrappers
libc="../build/python/src/libcensdam.a"

# Get all objects from these libraries
ar x $libc
ar x $libf

# prepare a shared library that is usable by julia
gfortran -shared -o ensdam.so *.o -lsvml -lifcore

# Clean obkects files
rm -f *.F90.o *.f90.o

