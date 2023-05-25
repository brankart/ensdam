import pyensdam as edam
import numpy as np
from netCDF4 import Dataset

# Example illustrating the module: random
# =======================================

# The module contains the functions:
# - random.seed: Seed random number generator
# - random.seed_save: Save seed in restart file
# - random.seed_load: Load seed from restart file
# - random.check: Check random number generator
# - random.swap: Random array swapping
# - random.uniform: Draw random numbers with uniform distribution
# - random.normal: Draw random numbers with normal distribution
# - random.exp: Draw random numbers with exponential distribution
# - random.gamma: Draw random numbers with gamma distribution
# - random.beta: Draw random numbers with beta distribution
# - random.truncated_exp: Draw random numbers with truncated exponential distribution
# - random.truncated_normal: Draw random numbers with truncated normal distribution
# - random.truncated_normal_vec: Draw random vectors with truncated normal distribution
# - random.field1d_init: Initialization for the sampling of 1D random fields
# - random.field2d_init: Initialization for the sampling of 2D random fields
# - random.field1d_sample: Sample 1D random fields with given spectrum
# - random.field2d_sample: Sample 2D random fields with given spectrum

# Seed the random number generator
# --------------------------------

print('--------------------------------')
print('Seed the random number generator')
print('--------------------------------')

# Perform short check of the random number generator (silent, but stop if failed)
edam.random.check('short')
# Perform long check of the random number generator (with message, and stop if failed)
edam.random.check('long')  # comment out to save time

# Use default seed
edam.random.seed(0)

# Save current state in restart file
edam.random.seed_save()

# Use seed index 8
edam.random.seed(8)

# Load seed from restart file (i.e. go back to default seed)
edam.random.seed_load()

# Generate random numbers with various distributions
# --------------------------------------------------

print('--------------------------------------------------')
print('Generate random numbers with various distributions')
print('--------------------------------------------------')

zran=edam.random.uniform([5])
print ('  Uniform random number:',zran)

zran=edam.random.normal([5])
print ('  Normal random number:',zran)

zran=edam.random.exp([5])
print ('  Exponential random number:',zran)

zran=edam.random.gamma(1.,[5])
print ('  Gamma random number:',zran)

zran=edam.random.beta(1.,1.,[5])
print ('  Beta random number:',zran)

zran=edam.random.truncated_exp(1.,[5])
print ('  Truncated exponential random number:',zran)

zran=edam.random.truncated_normal(1.,2.,[5])
print ('  Truncated normal random number:',zran)

# Randomly swap elements of input array
# -------------------------------------

print('-------------------------------------')
print('Randomly swap elements of input array')
print('-------------------------------------')

a = np.arange(10,dtype=np.intc)
print ('  Ordered array:',a)
edam.random.swap(a)
print ('  Randomly swapped array:',a)

# Sample random vector with multivariate truncated Gaussian distribution
# ----------------------------------------------------------------------

print('----------------------------------------------------------------------')
print('Sample random vector with multivariate truncated Gaussian distribution')
print('----------------------------------------------------------------------')

# Define constraint (Ax <= b, here: -x-y <= 0 i.e. x+y >= 0)
# Just add rows to A and B to impose more inequality constraints
A = np.ones((1,2)) - 2  # A = [ -1, -1 ]
b = np.zeros(1)         # b = 0
# Sample random vectors with truncated Gaussian distribution
# All iterates of the Gibbs sampler are output
# -> drop first draws and subsample the sequence of draws as appropriate
sample = edam.random.truncated_normal_vec(100,A,b);
for i in range(9,99,10):
  print ('  Truncated normal random vector (x+y>0):',sample[:,i])

# Sample 1D random field with given continuous spectrum
# -----------------------------------------------------

print('-----------------------------------------------------')
print('Sample 1D random field with given continuous spectrum')
print('-----------------------------------------------------')

# Definition of the spectrum (discretization of a continuous spectrum)
f0 = 0.1 ; df = 0.1 ; nf = 10
spct_freq = np.arange(f0, f0 + nf * df, df, dtype=np.float64)
spct_power = np.ones(10)

# Initialization of the sampler
edam.random.field1d_init(spct_freq,spct_power)

# Definition of the output grid (here global on the sphere)
x0 = 0. ; dx = 0.5 ; nx = 100
x = np.arange(x0, x0 + nx * dx, dx, dtype=np.float64)

# Generate random field
nharm = 100 # number of harmonics to superpose (sampled from the continuous spectrum)
field1d = edam.random.field1d_sample(x,nharm)

print('  Random field 1d stored in file random_field_1d.nc')
nc_field1d = Dataset('random_field_1d.nc', 'w', format='NETCDF4')
nc_field1d.createDimension('x', nx)
nc_field1d_x = nc_field1d.createVariable('x', 'f4', ('x'))
nc_field1d_field = nc_field1d.createVariable('random_field_1d', 'f4', ('x'))
nc_field1d_x[:] = x
nc_field1d_field[:] = field1d
nc_field1d.close()

# Sample 2D random field with given continuous (isotropic) spectrum
# -----------------------------------------------------------------

print('-----------------------------------------------------------------')
print('Sample 2D random field with given continuous (isotropic) spectrum')
print('-----------------------------------------------------------------')

# Initialization of the sampler (with same power spectrum as before)
edam.random.field2d_init(spct_freq,spct_power)

# Definition of the output grid (here global on the sphere)
xx, yy = np.meshgrid(x, x)

# Generate random field
nharm = 1000 # number of harmonics to superpose (sampled from the continuous spectrum)
field2d = edam.random.field2d_sample(xx,yy,nharm)

print('  Random field 2d stored in file random_field_2d.nc')
nc_field2d = Dataset('random_field_2d.nc', 'w', format='NETCDF4')
nc_field2d.createDimension('x', nx)
nc_field2d.createDimension('y', nx)
random_field_2d = nc_field2d.createVariable('random_field_2d', 'f4', ('x','y'))
random_field_2d[:,:] = field2d
nc_field2d.close()

# Sample 2D random field on the sphere with given spectrum
# in the basis of the spherical harmonics
# --------------------------------------------------------

print('--------------------------------------------------------')
print('Sample 2D random field on the sphere with given spectrum')
print('--------------------------------------------------------')

# Callback function defining the requested power spectrum
# as a function of degree l and order m of the spherical harmonics
lmin=0 ; lmax=30
def pow_spectrum(l,m):
  global lmax

  # Power spectrum
  lc = 6 # characteristic degree of the spherical harmonics
  power = 1. / ( 1. + l*l/(lc*lc) )

  # Normalize the spectrum
  norm = 0.
  for ll in range(0,lmax+1):
    norm = norm + 1. / ( 1. + ll*ll/(lc*lc) )
  power = power / norm

  # Scale to account for the multiplicity of each degree
  power = power / ( 1. + 2. * l )

  return power

# Definition of the output grid (here global on the sphere)
lon0 = 0. ; dlon = 1 ; nlon=360
lat0 = -90. ; dlat = 1 ; nlat=181
lon = np.arange(lon0, lon0 + nlon * dlon, dlon, dtype=np.double)
lat = np.arange(lat0, lat0 + nlat * dlat, dlat, dtype=np.double)
lon2d, lat2d = np.meshgrid(lon, lat)
lon1d = np.ravel(lon2d)
lat1d = np.ravel(lat2d)

# Generate random field
field2s_1d = edam.random.field2s_sample(lon1d,lat1d,pow_spectrum,lmin,lmax)
field2s = np.reshape(field2s_1d, lon2d.shape)

# Output in NetCDF file
print('  Random field on the sphere stored in file random_field_2s.nc')
nc_field2s = Dataset('random_field_2s.nc', 'w', format='NETCDF4')
nc_field2s.createDimension('lon', nlon)
nc_field2s.createDimension('lat', nlat)
nc_field2s_lon = nc_field2s.createVariable('lon', lon2d.dtype, ('lon'))
nc_field2s_lat = nc_field2s.createVariable('lat', lat2d.dtype, ('lat'))
nc_field2s_field = nc_field2s.createVariable('random_field_2s', field2s.dtype, ('lat','lon'))
nc_field2s_lon[:] = lon
nc_field2s_lat[:] = lat
nc_field2s_field[:,:] = field2s
nc_field2s.close()

