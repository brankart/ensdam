import pyensdam_mpi as edam
import numpy as np
from netCDF4 import Dataset
from mpi4py import MPI

# Example illustrating the module: augment
# ========================================

# The module contains the functions:
# - sample_mcmc : resample input ensemble with MCMC sampler,
#                 using covariance localization

# Get information from MPI
mpicomm = MPI.COMM_WORLD
mpirank = mpicomm.Get_rank()
mpisize = mpicomm.Get_size()

# Callback function defining the power spectrum of ensemble members
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

# Definition of the output grid (here global on the sphere)
lon0 = 0. ; dlon = 1 ; nlon=360
lat0 = -90. ; dlat = 1 ; nlat=181
lon = np.arange(lon0, lon0 + nlon * dlon, dlon, dtype=np.double)
lat = np.arange(lat0, lat0 + nlat * dlat, dlat, dtype=np.double)
lon2d, lat2d = np.meshgrid(lon, lat)
lon1dglo = np.ravel(lon2d)
lat1dglo = np.ravel(lat2d)

# Subsample grid points for current processor
lon1d = np.ascontiguousarray(lon1dglo[mpirank::mpisize])
lat1d = np.ascontiguousarray(lat1dglo[mpirank::mpisize])

if mpirank == 0:
  print('-----------------------------------------------')
  print('1. Preliminary step: generate original ensemble')
  print('-----------------------------------------------')

# Generate random ensemble using the ensdam.random functions
# See example_random.py for more details

# Ensemble size
m=50

# Generate random ensemble as an example
lmin=0 ; lmax=30 ; lc=6
for i in range(m):
  if mpirank == 0:
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

# Generate the multiscale ensemble as an array of ensembles
ensmulti = np.array([ens])
del ens

if mpirank == 0:
  print('---------------------------------------------------------')
  print('2. Preliminary step: generate large-scale random patterns')
  print('---------------------------------------------------------')

# Generate large-scale random pattern for covariance localization
lmin=0 ; lmax=6 ; lc=1.2
for i in range(m):
  if mpirank == 0:
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

# Renormalize the large-scale patterns to unit std (to describe a localizing correlation matrix)
mean, std = edam.statistics.meanstd(ens)
for i in range(m):
  ens[i,:] = ( ens[i,:] - mean[:] ) / std[:]

# Include large-scale patterns in multiscale ensemble
ensmulti.resize((ensmulti.shape[0] + 1, ensmulti.shape[1], ensmulti.shape[2]))
ensmulti[-1] = ens
del ens

# Output multiscale ensemble in NetCDF file
if mpirank == 0:
  print('  First member of multiscale ensemble stored in file: augment_example.nc')
  nc_ens = Dataset('augment_example.nc', 'w', format='NETCDF4')
  nc_ens.createDimension('lon', nlon)
  nc_ens.createDimension('lat', nlat)
  nc_ens_lon = nc_ens.createVariable('lon', lon.dtype, ('lon'))
  nc_ens_lat = nc_ens.createVariable('lat', lat.dtype, ('lat'))
  nc_ens_ens0 = nc_ens.createVariable('ens0', ensmulti.dtype, ('lat','lon'))
  nc_ens_ens1 = nc_ens.createVariable('ens1', ensmulti.dtype, ('lat','lon'))
  nc_ens_lon[:] = lon
  nc_ens_lat[:] = lat

# Reconstruct global member from the contribution of every processor
mem1dglo = np.zeros_like(lon1dglo)
mem1dglo[mpirank::mpisize] = ensmulti[0,0,:]
mpicomm.Allreduce(MPI.IN_PLACE, mem1dglo, op=MPI.SUM)

# Save scale 0 of member 0
if mpirank == 0:
  mem2d = np.reshape(mem1dglo, lon2d.shape)
  nc_ens_ens0[:,:] = mem2d

# Reconstruct global member from the contribution of every processor
mem1dglo = np.zeros_like(lon1dglo)
mem1dglo[mpirank::mpisize] = ensmulti[1,0,:]
mpicomm.Allreduce(MPI.IN_PLACE, mem1dglo, op=MPI.SUM)

# Save scale 1 of member 0
if mpirank == 0:
  mem2d = np.reshape(mem1dglo, lon2d.shape)
  nc_ens_ens1[:,:] = mem2d

if mpirank == 0:
  print('---------------------------------------------')
  print('3. Generate new members with the MCMC sampler')
  print('---------------------------------------------')

naug = 1 # number of augmented members requested
maxchain = 1000 # number of iteration in the MCMC sampler
multiplicity = np.array([1, 4],np.intc) # mutliplicity of each scale in the Schur products
augens = edam.augment.sample_mcmc(ensmulti,multiplicity,naug,maxchain)

# Reconstruct global member from the contribution of every processor
mem1dglo = np.zeros_like(lon1dglo)
mem1dglo[mpirank::mpisize] = augens[0,:]
mpicomm.Allreduce(MPI.IN_PLACE, mem1dglo, op=MPI.SUM)

# Output augmented member in NetCDF file
# Looking at the new member, we can see that the local correlation structure
# is preserved, as in the original ensemble, bu the global correlation structure is lost
if mpirank == 0:
  print('  First member of augmented ensemble stored in file: augment_example.nc')
  nc_ens_aug = nc_ens.createVariable('aug', augens.dtype, ('lat','lon'))
  mem2d = np.reshape(mem1dglo, lon2d.shape)
  nc_ens_aug[:,:] = mem2d

# Close NetCDF file
if mpirank == 0:
  nc_ens.close()

# Finalize MPI
MPI.Finalize()
