import pyensdam as edam
import numpy as np

# Example illustrating the module: interpolation
# ==============================================

# The module contains the functions:
# - locate1D : locate positions in 1D grid and compute interpolation weights
# - interp1D : apply interpolation on 1D input field
# - define2D : define 2D grid
# - locate2D : locate positions in 2D grid and compute interpolation weights
# - interp2D : apply interpolation on 2D input field

# Print documentation
print(edam.interpolation.__doc__)
print(edam.interpolation.locate1D.__doc__)
print(edam.interpolation.interp1D.__doc__)
print(edam.interpolation.define2D.__doc__)
print(edam.interpolation.locate2D.__doc__)
print(edam.interpolation.interp2D.__doc__)

