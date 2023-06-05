import pyensdam as edam
import numpy as np
from netCDF4 import Dataset

# Example illustrating the module: anamorphosis
# =============================================

# The module contains the functions:
# - anamorphosis.quantiles : compute ensemble quantiles
# - anamorphosis.forward : forward anamorphosis transformation
# - anamorphosis.backward : backward anamorphosis transformation
# - anamorphosis.forward_obs : forward anamorphosis transformation of observations
# - anamorphosis.forward_obs_sym : forward anamorphosis transformation of observations (symmetric)

# The module contains the parameters:
# - anamorphosis.target: target probability distribution to use (default=normal, uniform, gamma, beta, custom)
# - anamorphosis.quaref : quantiles of the target probability distribution to use (if target == 'custom'))
# - anamorphosis.obstype: probability distribution of observations (default=normal, gamma, beta)

# Notes:
# -

print('-----------------------------------------------')
print('1. Preliminary step: generate original ensemble')
print('-----------------------------------------------')

# Generate a 1D random ensemble using the ensdam.random functions
# See example_random.py for more details

# Definition of the spectrum (discretization of a continuous spectrum)
f0 = 0.1 ; df = 0.1 ; nf = 10
spct_freq = np.arange(f0, f0 + nf * df, df, dtype=np.float64)
spct_power = np.ones(10)

# Initialization of the sampler
edam.random.field1d_init(spct_freq,spct_power)

# Definition of the output grid (here global on the sphere)
x0 = 0. ; dx = 0.5 ; nx = 101
x = np.arange(x0, x0 + nx * dx, dx, dtype=np.float64)

# Generate random ensemble as an example
nens = 100      # ensemble size
nharm = 200 # number of harmonics to superpose in each member (sampled from the continuous spectrum)
for i in range(nens):
  print('  Generating ensemble member:',i)
  # Generate new random member
  member_1d = edam.random.field1d_sample(x,nharm)
  # Transform the Gaussian distribution into a beta distribution
  edam.probability.type = 'normal'
  cdf = edam.probability.cdf(member_1d)
  edam.probability.type = 'beta'
  edam.probability.a = 2.
  edam.probability.b = 2.
  member_1d = edam.probability.invcdf(cdf)
  if i==0:
    # initialize ensemble with first member
    ens = np.array([member_1d])
  else:
    # include new member to ensemble
    ens.resize((ens.shape[0] + 1, ens.shape[1]))
    ens[-1] = member_1d

print('  Original ensemble stored in file anamorphosis_example.nc')
nc_field1d = Dataset('anamorphosis_example.nc', 'w', format='NETCDF4')
nc_field1d.createDimension('x', nx)
nc_field1d.createDimension('nens', nens)
nc_field1d_x = nc_field1d.createVariable('x', 'f4', ('x'))
nc_field1d_field = nc_field1d.createVariable('original_ensemble', 'f4', ('nens','x'))
nc_field1d_x[:] = x
nc_field1d_field[:,:] = ens

print('-------------------------------------------------')
print('2. Compute the quantiles of the original ensemble')
print('-------------------------------------------------')

# Definition of the quantiles
qua0 = 0. ; dqua = 0.05 ; nqua = 21
quadef = np.arange(qua0, qua0 + nqua * dqua, dqua, dtype=np.float64)

# Compute the quantiles of the original ensemble
qua = edam.anamorphosis.quantiles(ens,quadef)
# By default, this also prepares the quantiles of the target distribution (anamorphosis.quaref).
# If customized by the user (anamorphosis.target='custom'), quaref must have the same dimension as quadef.

print('  Ensemble quantiles stored in file anamorphosis_example.nc')
nc_field1d.createDimension('nqua', nqua)
nc_field1d_quantiles = nc_field1d.createVariable('ensemble_quantiles', 'f4', ('nqua','x'))
nc_field1d_quantiles[:,:] = qua

print('--------------------------------------------------------')
print('3. Perform forward anamorphosis of the original ensemble')
print('--------------------------------------------------------')

edam.anamorphosis.forward(ens,qua)

print('  Transformed ensemble stored in file anamorphosis_example.nc')
nc_field1d_transformed = nc_field1d.createVariable('transformed_ensemble', 'f4', ('nens','x'))
nc_field1d_transformed[:,:] = ens

print('------------------------------------------------------------')
print('4. Perform backward anamorphosis of the transformed ensemble')
print('------------------------------------------------------------')

edam.anamorphosis.backward(ens,qua)

print('  Restored ensemble stored in file anamorphosis_example.nc')
nc_field1d_restored = nc_field1d.createVariable('restored_ensemble', 'f4', ('nens','x'))
nc_field1d_restored[:,:] = ens

nc_field1d.close()
