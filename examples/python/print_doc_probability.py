import pyensdam as edam
import numpy as np

# Example illustrating the module: probability
# ============================================

# The module contains 4 functions:
# - pdf    : compute the probability density function
# - logpdf : compute the logartihm of the probability density function
# - cdf    : compute the cumulative distribution function
# - invcdf : compute the iinverse cumulative distribution function

# The module contains 1 parameter:
# - type : type of probability distribution (normal, gamma, beta)
# - k: shape parameter of the gamma distribution
# - a: parameter alpha of the beta distribution
# - b: parameter beta of the beta distribution

# Notes:
# - The logpdf of the gamma distribution is computed
#   by dropping a additive function of the k parameter to save cost

# Print documentation
print(edam.probability.__doc__)
print(edam.probability.pdf.__doc__)
print(edam.probability.logpdf.__doc__)
print(edam.probability.cdf.__doc__)
print(edam.probability.invcdf.__doc__)

