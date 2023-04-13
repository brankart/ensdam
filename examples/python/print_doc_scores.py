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

