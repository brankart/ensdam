
import ensdam
import numpy as np

# EnsDAM functions aliases
gran = ensdam.stochtools.kiss_gaussian
ens_meanstd = ensdam.ensstat.ensemble_meanstd
crps_score = ensdam.ensscores.crps_score
rcrv_score = ensdam.ensscores.rcrv_score

# Parameters of the example
m = 100      # Size of the ensemble
n = 1000     # Size of the state vector
sigma = 0.3  # Observation error standard deviation

# Sample prior ensemble from N(0,I) distribution
prior_ensemble = np.empty([n,m])
for i, j in np.ndindex(prior_ensemble.shape):
  prior_ensemble[i,j] = gran()

# Sample reference truth from the same distribution
reference_truth = np.empty(n)
for i in np.ndindex(reference_truth.shape):
  reference_truth[i] = gran()

# Generate observations by adding N(0,sigma) perturbations to the reference truth
observations = np.empty(n)
for i in np.ndindex(observations.shape):
  observations[i] = gran()

observations = reference_truth + sigma * observations

# Compute the posterior ensemble by conditioning the prior ensemble on observations
# i) Generate perttubations to observations for each ensemble member
posterior_ensemble = np.empty([n,m])
for i, j in np.ndindex(prior_ensemble.shape):
  posterior_ensemble[i,j] = observations[i] + sigma * gran()
# ii) Compute innovation for each ensemble member
posterior_ensemble = posterior_ensemble - prior_ensemble
# iii) Compute mean and variance of prior ensemble
ens_mean = np.empty(n)
ens_var = np.empty(n)
ens_mean,ens_var = ens_meanstd(prior_ensemble)
ens_var = ens_var * ens_var
# iii) Multiply by Kalman gain
for i, j in np.ndindex(prior_ensemble.shape):
  posterior_ensemble[i,j] = posterior_ensemble[i,j] * ens_var[i] / ( ens_var[i] + sigma*sigma )
# iv) Add prior ensemble
posterior_ensemble = posterior_ensemble + prior_ensemble

# Compute CRPS score, using reference truth as verification data
crps,crps_reliability,crps_resolution = crps_score(prior_ensemble,reference_truth)
print ('Prior CRPS reliability and resolution:    ',crps_reliability,crps_resolution)

crps,crps_reliability,crps_resolution = crps_score(posterior_ensemble,reference_truth)
print ('Posterior CRPS reliability and resolution:',crps_reliability,crps_resolution)

# Compute RCRV score, using reference truth as verification data
rcrv_bias,rcrv_spread = rcrv_score(prior_ensemble,reference_truth)
print ('Prior RCRV bias and spread:    ',rcrv_bias,rcrv_spread)

rcrv_bias,rcrv_spread = rcrv_score(posterior_ensemble,reference_truth)
print ('Posterior RCRV bias and spread:',rcrv_bias,rcrv_spread)

