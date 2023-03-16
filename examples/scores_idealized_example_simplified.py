
import ensdam
import numpy as np

# EnsDAM functions aliases
crps_score = ensdam.ensscores.crps_score
rcrv_score = ensdam.ensscores.rcrv_score
compute_ranks = ensdam.ensscores.compute_ranks

# Parameters of the example
m = 20       # Size of the ensemble
n = 1000     # Size of the state vector

mu = 0.0     # Ensemble mean
sigma = 1.0  # Error standard deviation

# Sample prior ensemble from N(0,I) distribution
prior_ensemble = np.random.normal(mu, sigma,(n,m))

# Sample reference truth from the same distribution
reference_truth = np.random.normal(mu, sigma,n)

# Compute CRPS score, using reference truth as verification data
crps,crps_reliability,crps_resolution = crps_score(prior_ensemble,reference_truth)
print ('Prior CRPS reliability and resolution:    ',crps_reliability,crps_resolution)

# Compute RCRV score, using reference truth as verification data
rcrv_bias,rcrv_spread = rcrv_score(prior_ensemble,reference_truth)
print ('Prior RCRV bias and spread:    ',rcrv_bias,rcrv_spread)

# Compute rank histogram
ranks,rank_histogram = compute_ranks(prior_ensemble,reference_truth)
print ('Prior rank histogram:',rank_histogram)

