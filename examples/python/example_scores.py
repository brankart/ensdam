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

# Parameters of the example
n = 1000     # number of variables
m = 20       # size of the ensemble

mu = 0.0     # Ensemble mean
sigma = 1.0  # Ensemble standard deviation

# Sample reference truth from N(0,I) distribution
reference_truth = np.random.normal(mu,sigma,n)

# Sample ensemble from N(0,I) distribution
ensemble = np.random.normal(mu,sigma,(m,n))

# Example partition of the data
partition = np.zeros(n,np.intc)
half=np.intc(n/2)
partition[0:half]=0  # 1st half belongs to component 0 of the partition
partition[half:n]=1  # 2nd half belongs to component 1 of the partition

# 1. Computation of rank histograms
print('---------------------------------')
print('1. Computation of rank histograms')
print('---------------------------------')

# Compute rank histogram
rank_histogram = edam.scores.rank_histogram(ensemble,reference_truth)
print ('  Rank histogram:',rank_histogram)

# Compute rank histogram, and output ranks
rank_histogram, ranks = edam.scores.rank_histogram(ensemble,reference_truth,histogram_only=False)
print ('  Ranks :',ranks[0:15],'...')

# 2. Computation of CRPS score
print('----------------------------')
print('2. Computation of CRPS score')
print('----------------------------')

# Example of change of the missing value (default=-9999.)
edam.scores.crps_missing_value=-1.

# Compute CRPS score, using reference truth as verification data
print('2.a Without partition of the data')
crps,crps_reliability,crps_resolution = edam.scores.crps(ensemble,reference_truth)
print ('  CRPS total:       ',crps)
print ('  CRPS reliability: ',crps_reliability)
print ('  CRPS resolution:  ',crps_resolution)

# Compute CRPS score, using reference truth as verification data, with partition of the data
print('2.b With partition of the data')
crps,crps_reliability,crps_resolution = edam.scores.crps(ensemble,reference_truth,partition=partition)
print ('  CRPS total (component 0):       ',crps[0])
print ('  CRPS reliability (component 0): ',crps_reliability[0])
print ('  CRPS resolution (component 0):  ',crps_resolution[0])
print ('  CRPS total (component 1):       ',crps[1])
print ('  CRPS reliability (component 1): ',crps_reliability[1])
print ('  CRPS resolution (component 1):  ',crps_resolution[1])

# 3. Computation of RCRV score
print('----------------------------')
print('3. Computation of RCRV score')
print('----------------------------')

# Example of change of the missing value (default=-9999.)
edam.scores.rcrv_missing_value=-1.

# Compute RCRV score, using reference truth as verification data
print('3.a Without partition of the data')
rcrv_bias,rcrv_spread = edam.scores.rcrv(ensemble,reference_truth)
print ('  RCRV bias:   ',rcrv_bias)
print ('  RCRV spread: ',rcrv_spread)

# Compute RCRV score, using reference truth as verification data, with partition of the data
print('3.b With partition of the data')
rcrv_bias,rcrv_spread = edam.scores.rcrv(ensemble,reference_truth,partition=partition)
print ('  RCRV bias (component 0):       ',rcrv_bias[0])
print ('  RCRV spread (component 0): ',rcrv_spread[0])
print ('  RCRV bias (component 1):       ',rcrv_bias[1])
print ('  RCRV spread (component 1): ',rcrv_spread[1])

print('3.c With anamorphosis')
edam.scores.rcrv_with_anamorphosis=True
edam.scores.rcrv_number_of_quantiles=21
rcrv_bias,rcrv_spread = edam.scores.rcrv(ensemble,reference_truth,partition=partition)
print ('  RCRV bias (component 0):       ',rcrv_bias[0])
print ('  RCRV spread (component 0): ',rcrv_spread[0])
print ('  RCRV bias (component 1):       ',rcrv_bias[1])
print ('  RCRV spread (component 1): ',rcrv_spread[1])


# 4. Computation of OPTIMALITY score
print('----------------------------------')
print('4. Computation of OPTIMALITY score')
print('----------------------------------')

# Redefine ensemble with reference truth as expected value
# and reduced standard deviation sigma
sigma=0.3
improved_ensemble = np.zeros((m,n))
for i in range(m):
   improved_ensemble[i,:] = np.random.normal(reference_truth,sigma,(n))

# Example of change of the missing value (default=-9999.)
edam.scores.optimality_missing_value=-1.

# Define type of observation error (normal, lognormal, gamma or beta)
edam.scores.obstype='normal'
print('  Type of observation error:',edam.scores.obstype)

# Compute optimality score
print('4.a Without partition of the data')
optimality_score = edam.scores.optimality(improved_ensemble,reference_truth,sigma)
print ('  Optimality score:',optimality_score)

std_obs = np.zeros(n)+sigma
optimality_score = edam.scores.optimality(improved_ensemble,reference_truth,std_obs)
print ('  Optimality score (vector sigma):',optimality_score)

print('4.b With partition of the data')
optimality_score = edam.scores.optimality(improved_ensemble,reference_truth,sigma,partition=partition)
print ('  Optimality score (component 0):',optimality_score[0])
print ('  Optimality score (component 1):',optimality_score[1])

# 5. Computation of ENTROPY score
print('-------------------------------')
print('5. Computation of ENTROPY score')
print('-------------------------------')

# Callback function for entropy score
def binary_event_outcomes(member):

  outcome = np.zeros(2,dtype=np.intc)

  # Event 1: is the mean square of the ensemble member larger than 1 ?
  if np.sum(np.square(member))/member.size < 1. :
    outcome[0] = 1
  else:
    outcome[0] = 2

  # Event 2: is the maximum of the ensemble member larger than 3.4 ?
  if np.amax(np.abs(member)) < 3.4 :
    outcome[1] = 1
  else:
    outcome[1] = 2

  return outcome

# Set reference probability distribution to equal probability for yes and no
pref = np.zeros((2,2))+0.5

# Compute entropy score
print('5.a Compute entropy score')
score_only=True
entropy_score =  edam.scores.entropy(ensemble,pref,binary_event_outcomes,score_only=score_only)
print ('  Entropy score for non-informative ensemble (event 1):',entropy_score[0])
print ('  Entropy score for non-informative ensemble (event 2):',entropy_score[1])

entropy_score =  edam.scores.entropy(improved_ensemble,pref,binary_event_outcomes,score_only=score_only)
print ('  Entropy score for improved ensemble (event 1):',entropy_score[0])
print ('  Entropy score for improved ensemble (event 2):',entropy_score[1])

print('5.b Compute other entropy components')
score_only=False
entropy_score,relative_entropy,cross_entropy,entropy =  edam.scores.entropy(ensemble,pref,binary_event_outcomes,score_only=score_only)
print ('  Entropy score for non-informative ensemble (event 1):',entropy_score[0])
print ('  Entropy score for non-informative ensemble (event 2):',entropy_score[1])
print ('  Relative entropy for non-informative ensemble (event 1):',relative_entropy[0])
print ('  Relative entropy for non-informative ensemble (event 2):',relative_entropy[1])
print ('  Cross entropy for non-informative ensemble (event 1):',cross_entropy[0])
print ('  Cross entropy for non-informative ensemble (event 2):',cross_entropy[1])
print ('  Entropy for non-informative ensemble (event 1):',entropy[0])
print ('  Entropy for non-informative ensemble (event 2):',entropy[1])
print ()

entropy_score,relative_entropy,cross_entropy,entropy =  edam.scores.entropy(improved_ensemble,pref,binary_event_outcomes,score_only=score_only)
print ('  Entropy score for improved ensemble (event 1):',entropy_score[0])
print ('  Entropy score for improved ensemble (event 2):',entropy_score[1])
print ('  Relative entropy for improved ensemble (event 1):',relative_entropy[0])
print ('  Relative entropy for improved ensemble (event 2):',relative_entropy[1])
print ('  Cross entropy for improved ensemble (event 1):',cross_entropy[0])
print ('  Cross entropy for improved ensemble (event 2):',cross_entropy[1])
print ('  Entropy for improved ensemble (event 1):',entropy[0])
print ('  Entropy for improved ensemble (event 2):',entropy[1])

