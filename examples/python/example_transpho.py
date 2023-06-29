import pyensdam as edam
import numpy as np
from netCDF4 import Dataset

# Example illustrating the module: transpho
# =========================================

# The module contains the functions:
# - transpho.forward : forward transformation (compute spectrum)
# - transpho.backward : backward transformation
# - transpho.mesh_area : compute area of grid cells on the sphere

# Callback function defining the power spectrum of the example field
# as a function of degree l and order m of the spherical harmonics
def pow_spectrum(l,m):
  global lmax, lc

  # Power spectrum
  power = 1. / ( 1. + l*l/(lc*lc) )

  # Normalize the spectrum
  norm = 0.
  for ll in range(0,lmax+1):
    norm = norm + 1. / ( 1. + ll*ll/(lc*lc) )
  power = power / norm

  # Scale to account for the multiplicity of each degree
  power = power / ( 1. + 2. * l )

  return power

# Definition of the output grid on the sphere
lon0 = 0.5   ; dlon = 1 ; nlon=360
lat0 = -89.5 ; dlat = 1 ; nlat=180
lon = np.arange(lon0, lon0 + nlon * dlon, dlon, dtype=np.double)
lat = np.arange(lat0, lat0 + nlat * dlat, dlat, dtype=np.double)
lon2d, lat2d = np.meshgrid(lon, lat)
lon1d = np.ravel(lon2d)
lat1d = np.ravel(lat2d)

print('---------------------------------------------------------')
print('1. Preliminary step: generate example field on the sphere')
print('---------------------------------------------------------')
# Generate random field using the ensdam.random functions
# See example_random.py for more details

lmin=0 ; lmax=30 ; lc=6
field = edam.random.field2s_sample(lon1d,lat1d,pow_spectrum,lmin,lmax)

print('-------------------------------------')
print('2. Compute the spectrum of this field')
print('-------------------------------------')

# Compute area associated to each grid point
lon0 = 0.   ; dlon = 1 ; nlone=361
lat0 = -90. ; dlat = 1 ; nlate=181
lone = np.arange(lon0, lon0 + nlone * dlon, dlon, dtype=np.double)
late = np.arange(lat0, lat0 + nlate * dlat, dlat, dtype=np.double)
lon2de, lat2de = np.meshgrid(lone, late)
area = edam.transpho.mesh_area(lon2de,lat2de)
area1d = np.ravel(area)

# Perform transformation
edam.transpho.lmax=30
edam.transpho.latres=0.5
spectrum = edam.transpho.forward(field*area1d,lon1d,lat1d)

print('----------------------------------')
print('3. Perform backward transformation')
print('----------------------------------')

# Transform back with the full spectrum
field_back = edam.transpho.backward(spectrum,lon1d,lat1d)
# Transform back keeping only spherical harmonics up to degree 10
#field_back = edam.transpho.backward(spectrum,lon1d,lat1d,l0=0,l1=10)

# Output in NetCDF file
print('  Example of the transpho module stored in: transpho_field.nc')
nc_field2s = Dataset('transpho_field.nc', 'w', format='NETCDF4')
nc_field2s.createDimension('lon', nlon)
nc_field2s.createDimension('lat', nlat)
nc_field2s_lon = nc_field2s.createVariable('lon', lon2d.dtype, ('lon'))
nc_field2s_lat = nc_field2s.createVariable('lat', lat2d.dtype, ('lat'))
nc_field2s_field = nc_field2s.createVariable('field', field.dtype, ('lat','lon'))
nc_field2s_field_back = nc_field2s.createVariable('field_back', field.dtype, ('lat','lon'))
nc_field2s_area = nc_field2s.createVariable('area', field.dtype, ('lat','lon'))
nc_field2s_lon[:] = lon
nc_field2s_lat[:] = lat
field2s = np.reshape(field, lon2d.shape)
nc_field2s_field[:,:] = field2s
field2s = np.reshape(field_back, lon2d.shape)
nc_field2s_field_back[:,:] = field2s
nc_field2s_area[:,:] = area
nc_field2s.close()

