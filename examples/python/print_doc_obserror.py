import pyensdam as edam
import numpy as np

# Example illustrating the module: obserror
# =========================================

# The module contains 3 functions:
# - logpdf : compute the logarithm of the observation error pdf,
#            it is used to compute the observation cost function
# - cdf    : compute the observation error cdf,
#            it is needed in the computation of the optimality score
# - sample : sample the observation error probability distribution,
#            it is used to perturb a reference state with observation errors

# The module contains 1 parameter:
# - obstype : type of observation error probability distribution to use (default=normal), 
#             supported distributions: normal, lognormal, gamma and beta

# Notes:
# - When applied to observation vectors, observation errors are assumed independent,
#   only marginal distributions for each component are used.
#   In logpdf, the contributions of the vector components are summed.
# - The random number generator used by this module is in the module: stochtools,
#   it can be seeded using functions provided there, this example uses the default seed.

# Print documentation
print(edam.obserror.__doc__)
print(edam.obserror.logpdf.__doc__)
print(edam.obserror.cdf.__doc__)
print(edam.obserror.sample.__doc__)

