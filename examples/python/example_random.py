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

# The module contains the parameters:
# - 

# Notes:
# -

#edam.random.seed_save()

#edam.random.seed(0)

#edam.random.seed(1)

#edam.random.seed_save()

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

a = np.arange(10,dtype=np.intc)
print ('  Ordered array:',a)
edam.random.swap(a)
print ('  Randomly swapped array:',a)

# Sample random vector with multivariate truncated Gaussian distribution
# ----------------------------------------------------------------------

# Define constraint (-x-y <= 0)
A = np.ones((1,2)) - 2  # A = [ -1, -1 ]
b = np.zeros(1)         # b = 0
# Sample random vectors with truncated Gaussian distribution
sample = edam.random.truncated_normal_vec(100,A,b);
for i in range(9,99,10):
  print ('  Truncated normal random vector (x+y>0):',sample[:,i])

# Sample 1D random field with given continuous spectrum
# -----------------------------------------------------

f0 = 0.1 ; df = 0.1 ; nf = 10
spct_freq = np.arange(f0, f0 + nf * df, df, dtype=np.float64)
spct_power = np.ones(10)
edam.random.field1d_init(spct_freq,spct_power)

x0 = 0. ; dx = 0.1 ; nx = 100
x = np.arange(x0, x0 + nx * dx, dx, dtype=np.float64)
nharm = 100 # number of harmonics to superpose

field1d = edam.random.field1d_sample(x,nharm)
print('  Random field 1d',field1d)

nc_field1d = Dataset('random_field_1d.nc', 'w', format='NETCDF4')
nc_field1d.createDimension('x', nx)
nc_field1d_x = nc_field1d.createVariable('x', 'f4', ('x'))
nc_field1d_field = nc_field1d.createVariable('random_field_1d', 'f4', ('x'))
nc_field1d_x[:] = x
nc_field1d_field[:] = field1d
nc_field1d.close()

# Sample 2D random field with given continuous (isotropic) spectrum
# -----------------------------------------------------------------

edam.random.field2d_init(spct_freq,spct_power)

xx, yy = np.meshgrid(x, x)
field2d = edam.random.field2d_sample(xx,yy,nharm)

nc_field2d = Dataset('random_field_2d.nc', 'w', format='NETCDF4')
nc_field2d.createDimension('x', nx)
nc_field2d.createDimension('y', nx)
random_field_2d = nc_field2d.createVariable('random_field_2d', 'f4', ('x','y'))
random_field_2d[:,:] = field2d
nc_field2d.close()


# Sample 2D random field on the sphere with given spectrum
# in the basis of the spherical harmonics
# --------------------------------------------------------

lon0 = 0. ; dlon = 1 ; nlon=360
lat0 = -90. ; dlat = 1 ; nlat=181
lon = np.arange(lon0, lon0 + nlon * dlon, dlon, dtype=np.double)
lat = np.arange(lat0, lat0 + nlat * dlat, dlat, dtype=np.double)
lon2d, lat2d = np.meshgrid(lon, lat)
lon1d = np.ravel(lon2d)
lat1d = np.ravel(lat2d)

lmin=0 ; lmax=20
field2s_1d = lon1d
#field2s_1d = edam.random.field2s_sample(lon1d,lat1d,pow_spectrum,lmin,lmax

field2s = np.reshape(field2s_1d, lon2d.shape)

nc_field2s = Dataset('random_field_2s.nc', 'w', format='NETCDF4')
nc_field2s.createDimension('lon', nlon)
nc_field2s.createDimension('lat', nlat)
nc_field2s_lon = nc_field2s.createVariable('lon', lon2d.dtype, ('lat','lon'))
nc_field2s_lat = nc_field2s.createVariable('lat', lat2d.dtype, ('lat','lon'))
nc_field2s_field = nc_field2s.createVariable('random_field_2s', field2s.dtype, ('lat','lon'))
#nc_field2s_lon[:,:] = lon2d
#nc_field2s_lat[:,:] = lat2d
nc_field2s_field[:,:] = field2s
nc_field2s.close()

