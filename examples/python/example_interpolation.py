import pyensdam as edam

# Example illustrating the module: interpolation
# ==============================================

# The module contains the functions:
# - locate1D : locate positions in 1D grid and compute interpolation weights
# - interp1D : apply interpolation on 1D input field
# - define2D : define 2D grid
# - locate2D : locate positions in 2D grid and compute interpolation weights
# - interp2D : apply interpolation on 2D input field

print('----------------------')
print('1. Interpolation in 1D')
print('----------------------')

# Generate a 1D random ensemble using the ensdam.random functions
# See example_random.py for more details

# Definition of the spectrum (discretization of a continuous spectrum)
f0 = 0.1 ; df = 0.1 ; nf = 10
spct_freq = np.arange(f0, f0 + nf * df, df, dtype=np.float64)
spct_power = np.ones(10)

# Initialization of the sampler
edam.random.field1d_init(spct_freq,spct_power)

# Definition of the oroginal grid (here global on the sphere)
x0 = 0. ; dx = 0.5 ; nx = 101
x = np.arange(x0, x0 + nx * dx, dx, dtype=np.float64)

