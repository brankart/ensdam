#!/etc/bash

# Make sure to have access to cmake, fortran compiler, python
# and then build the python module with:

cmake -B build -S .

cmake --build build --target wheel && pip install --user build/python/dist/pyensdam-0.1.1-py3-none-any.whl --force-reinstall
