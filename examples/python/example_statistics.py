import pyensdam as edam
import numpy as np
from netCDF4 import Dataset

# Example illustrating the module: statistics
# ===========================================

# The module contains 4 functions:
# - statistics.meanstd : compute ensemble mean and standard deviation
# - statistics.correlation : compute ensemble correlation
# - statistics.covariance : compute ensemble covariance
# - statistics.representer : compute ensemble representer

print('-------------------------------')
print('Statistics of a scalar ensemble')
print('-------------------------------')

# Ensemble size
m=5

# Generate a 5-member scalar ensemble of normal random numbers
ens = edam.random.normal([m])
print ('  Input scalar ensemble:',ens)

# Compute ensemble mean only
mean = edam.statistics.meanstd(ens,std=False)
print ('  Ensemble mean:',mean)

# Compute ensemble mean and standard deviation
mean,std = edam.statistics.meanstd(ens)
print ('  Ensemble mean, std:',mean,std)

# Generate weight of ensemble members
weight = edam.random.gamma(1.,[m])
print ('  Input ensemble weights:',weight)

# Compute ensemble mean only
mean = edam.statistics.meanstd(ens,weight=weight,std=False)
print ('  Weighted ensemble mean:',mean)

# Compute ensemble mean and standard deviation
mean,std = edam.statistics.meanstd(ens,weight=weight)
print ('  Weighted ensemble mean, std:',mean,std)

print('-------------------------------------')
print('Statistics of a multivariate ensemble')
print('-------------------------------------')

# Ensemble size
m=5

# Generate random ensemble using the ensdam.random functions
# See example_random.py for more details

# Callback function defining the power spectrum
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

# Generate random ensemble as an example
for i in range(m):
  print('  Generating ensemble member:',i)
  # Generate new random member
  member_1d = edam.random.field2s_sample(lon1d,lat1d,pow_spectrum,lmin,lmax)
  if i==0:
    # initialize ensemble with first member
    ens = np.array([member_1d])
  else:
    # include new member to ensemble
    ens.resize((ens.shape[0] + 1, ens.shape[1]))
    ens[-1] = member_1d

# Refrence ensemble with which computing ensemble correlation
ensref = np.ascontiguousarray(ens[:,nlon*nlat//2])
#print('contiguous:',ensref.flags['C_CONTIGUOUS'])  # False

# Compute ensemble mean and std
print('  Computing ensemble mean and std')
mean_1d, std_1d = edam.statistics.meanstd(ens)
mean = np.reshape(mean_1d, lon2d.shape)
std  = np.reshape(std_1d,  lon2d.shape)

# Compute ensemble correlation
print('  Computing ensemble correlation')
correl_1d = edam.statistics.correlation(ens,ensref)
correl = np.reshape(correl_1d,  lon2d.shape)

# Output in NetCDF file
print('  Statistics stored in file statistics_2s.nc')
nc_field2s = Dataset('statistics_2s.nc', 'w', format='NETCDF4')
nc_field2s.createDimension('lon', nlon)
nc_field2s.createDimension('lat', nlat)
nc_field2s_lon = nc_field2s.createVariable('lon', lon2d.dtype, ('lon'))
nc_field2s_lat = nc_field2s.createVariable('lat', lat2d.dtype, ('lat'))
nc_field2s_mean = nc_field2s.createVariable('mean', mean.dtype, ('lat','lon'))
nc_field2s_std = nc_field2s.createVariable('std', mean.dtype, ('lat','lon'))
nc_field2s_correl = nc_field2s.createVariable('correl', mean.dtype, ('lat','lon'))
nc_field2s_lon[:] = lon
nc_field2s_lat[:] = lat
nc_field2s_mean[:,:] = mean
nc_field2s_mean[:,:] = std
nc_field2s_correl[:,:] = correl
nc_field2s.close()

