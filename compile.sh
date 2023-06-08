#!/bin/sh

# Make sure to have access to cmake, fortran compiler, python
# and then build the python module as follows.

# To know the default compilation settings: cmake -LAH
# Edit the flags to mmodify the default

flags = ""

# Uncomment to explicitly specify the python3 exectuable
# flags="$flags -DPython3_EXECUTABLE=/gpfslocalsup/pub/anaconda-py3/2021.05/envs/python-3.10.4/bin/python3"

# Uncomment to explicitly specify the Fortran compiler
# flags="$flags -DCMAKE_Fortran_COMPILER=mpiifort"

# Uncomment to use MPI
# flags="$flags -DCMAKE_Fortran_FLAGS=-DMPI"

cmake $flags -B build -S .

cmake --build build --target wheel && pip install --user build/python/dist/pyensdam-0.1.1-py3-none-any.whl --force-reinstall
