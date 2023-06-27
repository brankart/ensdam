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
lon0 = 0. ; dlon = 1 ; nlon=360
lat0 = -90. ; dlat = 1 ; nlat=181
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

edam.transpho.lmax=30
edam.transpho.latres=0.5
spectrum = edam.transpho.forward(field,lon1d,lat1d)

print('shape',spectrum.shape)

print('----------------------------------')
print('3. Perform backward transformation')
print('----------------------------------')

field_back = edam.transpho.backward(spectrum,lon1d,lat1d)


print('max difference',np.max(np.abs(field-field_back)))
