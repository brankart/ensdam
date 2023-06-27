import pyensdam as edam

# Example illustrating the module: update
# =======================================

# The module contains the functions:
# - update.sample_mcmc : apply new condition (e.g. observations) on input ensemble,
#                        using MCMC sampler with covariance localization
#
# Module parameters:
# - update.chain_index : current chain index
# - update.zero_start : start from zero (T) or from restart ensemble (F)
# - update.control_print : number of iteration between control prints
# - update.convergence_check : number of iteration between convergence checks
# - update.convergence_stop : stop at convergence (T) or perform full requested iterations (F)



# Print documentation
print(edam.update.__doc__)
print(edam.update.sample_mcmc.__doc__)
