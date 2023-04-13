import pyensdam as edam
import numpy as np

# Example illustrating the module: scores
# =======================================

# The module contains 5 functions:
# - rank_histogram : compute rank histograms (with option to output ranks)
# - crps : compute CRPS score (total, reliability, resolution)
# - rcrv : compute RCRV score (bias, spread)
# - optimality : compute OPTIMALITY score
# - entropy : compute ENTROPY score (with option to compute entropy components)

# The module contains 6 parameters:
# - crps_missing_value : missing value for CRPS score
# - rcrv_missing_value : missing value for RCRV score
# - rcrv_with_anamorphosis : apply anamorphosis rather than center-reduction in RCRV score
# - rcrv_number_of_quantiles : number of quantiles used in the anamorphosis transformation
# - optimality_missing_value : missing value for OPTIMALITY score
# - entropy_base : basis for the logarithm in entropy computations

# Notes:
# - CRPS, RCRV and OPTIMALITY scores have the option to partition the input data
#   and compute the score separately for each element of the partition.

# Print documentation
print(edam.scores.__doc__)
print(edam.scores.rank_histogram.__doc__)
print(edam.scores.crps.__doc__)
print(edam.scores.rcrv.__doc__)
print(edam.scores.optimality.__doc__)
print(edam.scores.entropy.__doc__)

# Parameters of the example
n = 10     # number of variables
m = 20     # size of the ensemble

mu = 0.0     # Ensemble mean
sigma = 1.0  # Ensemble standard deviation

# Example ensemble sampled from N(0,I) distribution
prior_ensemble = np.random.normal(mu, sigma,(n,m))

# Example reference truth sampled from the same distribution
reference_truth = np.random.normal(mu, sigma,n)

# Example partition of the data
partition = np.zeros(n,np.intc)
half=np.intc(n/2)
partition[0:half]=1  # 1st half belongs to component 1 of the partition
partition[half:n]=2  # 2nd half belongs to component 2 of the partition

# Compute CRPS score, using reference truth as verification data
edam.scores.crps_missing_value=-1.
crps,crps_reliability,crps_resolution = edam.scores.crps(prior_ensemble,reference_truth)
print ('CRPS:                                     ',crps)
print ('Prior CRPS reliability and resolution:    ',crps_reliability,crps_resolution)

crps,crps_reliability,crps_resolution = edam.scores.crps(prior_ensemble,reference_truth,partition=partition)
print ('CRPS:                                     ',crps)
print ('Prior CRPS reliability and resolution:    ',crps_reliability,crps_resolution)

partition[1:half]=2
crps,crps_reliability,crps_resolution = edam.scores.crps(prior_ensemble,reference_truth,partition=partition)
print ('CRPS:                                     ',crps)
print ('Prior CRPS reliability and resolution:    ',crps_reliability,crps_resolution)

# Compute RCRV score, using reference truth as verification data
rcrv_bias,rcrv_spread = edam.scores.rcrv(prior_ensemble,reference_truth)
print ('Prior RCRV bias and spread:    ',rcrv_bias,rcrv_spread)

# Compute rank histogram
rank_histogram = edam.scores.rank_histogram(prior_ensemble,reference_truth)
print ('Prior rank histogram:',rank_histogram)

rank_histogram, ranks = edam.scores.rank_histogram(prior_ensemble,reference_truth,histogram_only=False)
print ('Prior ranks :',ranks[0:15],'.....')

# Compute optimality score
optimality_score = edam.scores.optimality(prior_ensemble,reference_truth,sigma)
print ('Optimality score:',optimality_score)
print ('obstype',edam.scores.obstype)

edam.scores.obstype='normal'
#edam.scores.obssigma = np.zeros(n)+1.
std_obs = np.zeros(n)+1.
optimality_score = edam.scores.optimality(prior_ensemble,reference_truth,std_obs)
print ('Optimality score:',optimality_score)
print ('obstype',edam.scores.obstype)

sigma=2.
optimality_score = edam.scores.optimality(prior_ensemble,reference_truth,sigma)
print ('Optimality score:',optimality_score)
print ('obstype',edam.scores.obstype)

std_obs = np.zeros(n)+1.
optimality_score = edam.scores.optimality(prior_ensemble,reference_truth,std_obs)
print ('Optimality score:',optimality_score)
print ('obstype',edam.scores.obstype)

# Callback function for entropy score
def binary_event_outcomes(member):

  outcome = np.zeros(2,dtype=np.intc)

  if np.sum(np.square(member))/(member.size-1) < 1. :
    outcome[0] = 1
  else:
    outcome[0] = 2

  if np.amax(np.abs(member)) < 3.4 :
    outcome[1] = 1
  else:
    outcome[1] = 2

  return outcome

pref = np.zeros((2,2))+0.5
score_only=False
entropy_score =  edam.scores.entropy(prior_ensemble,pref,binary_event_outcomes,score_only=score_only)
print ('Entropy score:',entropy_score)
#if not score_only :
  #print ('Relative entropy:',entropy_score.relative_entropy)
  #print ('Cross entropy:',entropy_score.cross_entropy)
  #print ('Entropy:',entropy_score.entropy)


