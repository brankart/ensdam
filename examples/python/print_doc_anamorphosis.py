import pyensdam as edam
import numpy as np

# Example illustrating the module: anamorphosis
# =============================================

# The module contains the functions:
# - anamorphosis.quantiles : compute ensemble quantiles
# - anamorphosis.forward : forward anamorphosis transformation
# - anamorphosis.backward : backward anamorphosis transformation
# - anamorphosis.forward_obs : forward anamorphosis transformation of observations
# - anamorphosis.forward_obs_sym : forward anamorphosis transformation of observations (symmetric)

# The module contains the parameters:
# - anamorphosis.target: target probability distribution to use (default=normal, uniform, gamma, beta)
# - anamorphosis.obstype: probability distribution of observations (default=normal, gamma, beta)

# Notes:
# -

# Print documentation
print(edam.anamorphosis.__doc__)
print(edam.anamorphosis.forward.__doc__)
print(edam.anamorphosis.backward.__doc__)
print(edam.anamorphosis.forward_obs.__doc__)
print(edam.anamorphosis.forward_obs_sym.__doc__)

