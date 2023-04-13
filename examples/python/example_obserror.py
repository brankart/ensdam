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

# 1. Example with a normal distribution
print('-------------------------------------')
print('1. Example with a normal distribution')
print('-------------------------------------')
edam.obserror.obstype = 'normal'
print('Type of observation error:',edam.obserror.obstype)


# 1.a Scalar observation
print('1.a Scalar observation')
y = 1.       # Value of the observation yo
x = 0.       # Reference value in the conditional distribution p(yo|x)
             # This can be the true value or a current estimate
sigma = 1.   # Observation error standard deviation
print('  y=',y,'x=',x,'sigma=',sigma)

logpdf = edam.obserror.logpdf(y, x, sigma)
print ('  logpdf :',logpdf)

cdf = edam.obserror.cdf(y, x, sigma)
print ('  cdf :',cdf)

sample = edam.obserror.sample(x,sigma)
print ('  sample :',sample)

sample = edam.obserror.sample(x,sigma)
print ('  new sample, with same x=0 :',sample)

xmod=1
sample = edam.obserror.sample(xmod,sigma,reuse_last_rank=True)
print ('  new sample, with modified x=1, but same rank :',sample)

# 1.b Vector observation, with homogeneous observation error statistics
print('1.b Vector observation, with homogeneous observation error statistics')
n = 5    # Size of the observation vector
xvector = np.zeros(n)   # Reference value
yvector = np.random.normal(x, sigma, n) # Observations
print('  yvector=',yvector)

vlogpdf = edam.obserror.logpdf(yvector, xvector, sigma)
print ('  logpdf :',vlogpdf)

vcdf = edam.obserror.cdf(yvector, xvector, sigma)
print ('  cdf :',vcdf)

vsample = edam.obserror.sample(xvector,sigma)
print ('  sample :',vsample)

# 1.c Vector observation, with nonhomogeneous observation error statistics
print('1.c Vector observation, with homogeneous observation error statistics')
svector = xvector + 1
print ('  sigma vector :',svector)

vlogpdf = edam.obserror.logpdf(yvector, xvector, svector)
print ('  logpdf :',vlogpdf)

vcdf = edam.obserror.cdf(yvector, xvector, svector)
print ('  cdf :',vcdf)

vsample = edam.obserror.sample(xvector,svector)
print ('  sample :',vsample)

# 2. Example with a lognormal distribution
print('----------------------------------------')
print('2. Example with a lognormal distribution')
print('----------------------------------------')
edam.obserror.obstype = 'lognormal'
print('Type of observation error:',edam.obserror.obstype)

# 2.a Scalar observation
print('2.a Scalar observation')
y = 2.       # Value of the observation yo
x = 1.       # Reference value in the conditional distribution p(yo|x)
             # This can be the true value or a current estimate
sigma = 1.   # Observation error standard deviation
print('  y=',y,'x=',x,'sigma=',sigma)

logpdf = edam.obserror.logpdf(y, x, sigma)
print ('  logpdf :',logpdf)

cdf = edam.obserror.cdf(y, x, sigma)
print ('  cdf :',cdf)

sample = edam.obserror.sample(x,sigma)
print ('  sample :',sample)

# 2.b Vector observation, with homogeneous observation error statistics
print('2.b Vector observation, with homogeneous observation error statistics')
n = 5    # Size of the observation vector
xvector = np.ones(n)   # Reference value
yvector = yvector * yvector
print('  yvector=',yvector)

vlogpdf = edam.obserror.logpdf(yvector, xvector, sigma)
print ('  logpdf :',vlogpdf)

vcdf = edam.obserror.cdf(yvector, xvector, sigma)
print ('  cdf :',vcdf)

vsample = edam.obserror.sample(xvector,sigma)
print ('  sample :',vsample)

# 2.c Vector observation, with nonhomogeneous observation error statistics
print('2.c Vector observation, with homogeneous observation error statistics')
svector = xvector
print ('  sigma vector :',svector)

vlogpdf = edam.obserror.logpdf(yvector, xvector, svector)
print ('  logpdf :',vlogpdf)

vcdf = edam.obserror.cdf(yvector, xvector, svector)
print ('  cdf :',vcdf)

vsample = edam.obserror.sample(xvector,svector)
print ('  sample :',vsample)

# 3. Example with a gamma distribution
print('------------------------------------')
print('3. Example with a gamma distribution')
print('------------------------------------')
edam.obserror.obstype = 'gamma'
print('Type of observation error:',edam.obserror.obstype)

# 3.a Scalar observation
print('3.a Scalar observation')
y = 2.       # Value of the observation yo
x = 1.       # Reference value in the conditional distribution p(yo|x)
             # This can be the true value or a current estimate
sigma = 1.   # Observation error standard deviation
print('  y=',y,'x=',x,'sigma=',sigma)

logpdf = edam.obserror.logpdf(y, x, sigma)
print ('  logpdf :',logpdf)

cdf = edam.obserror.cdf(y, x, sigma)
print ('  cdf :',cdf)

sample = edam.obserror.sample(x,sigma)
print ('  sample :',sample)

# 3.b Vector observation, with homogeneous observation error statistics
print('3.b Vector observation, with homogeneous observation error statistics')
n = 5    # Size of the observation vector
xvector = np.ones(n)   # Reference value
print('  yvector=',yvector)

vlogpdf = edam.obserror.logpdf(yvector, xvector, sigma)
print ('  logpdf :',vlogpdf)

vcdf = edam.obserror.cdf(yvector, xvector, sigma)
print ('  cdf :',vcdf)

vsample = edam.obserror.sample(xvector,sigma)
print ('  sample :',vsample)

# 3.c Vector observation, with nonhomogeneous observation error statistics
print('3.c Vector observation, with homogeneous observation error statistics')
svector = xvector
print ('  sigma vector :',svector)

vlogpdf = edam.obserror.logpdf(yvector, xvector, svector)
print ('  logpdf :',vlogpdf)

vcdf = edam.obserror.cdf(yvector, xvector, svector)
print ('  cdf :',vcdf)

vsample = edam.obserror.sample(xvector,svector)
print ('  sample :',vsample)

# 4. Example with a beta distribution
print('-----------------------------------')
print('4. Example with a beta distribution')
print('-----------------------------------')
edam.obserror.obstype = 'beta'
print('Type of observation error:',edam.obserror.obstype)

# 4.a Scalar observation
print('4.a Scalar observation')
y = 0.7      # Value of the observation yo
x = 0.5      # Reference value in the conditional distribution p(yo|x)
             # This can be the true value or a current estimate
sigma = 0.2  # Observation error standard deviation
print('  y=',y,'x=',x,'sigma=',sigma)

logpdf = edam.obserror.logpdf(y, x, sigma)
print ('  logpdf :',logpdf)

cdf = edam.obserror.cdf(y, x, sigma)
print ('  cdf :',cdf)

sample = edam.obserror.sample(x,sigma)
print ('  sample :',sample)

# 4.b Vector observation, with homogeneous observation error statistics
print('4.b Vector observation, with homogeneous observation error statistics')
n = 5    # Size of the observation vector
xvector = 0.5 * np.ones(n)   # Reference value
yvector = 0.7 * np.ones(n)
print('  yvector=',yvector)

vlogpdf = edam.obserror.logpdf(yvector, xvector, sigma)
print ('  logpdf :',vlogpdf)

vcdf = edam.obserror.cdf(yvector, xvector, sigma)
print ('  cdf :',vcdf)

vsample = edam.obserror.sample(xvector,sigma)
print ('  sample :',vsample)

# 4.c Vector observation, with nonhomogeneous observation error statistics
print('4.c Vector observation, with homogeneous observation error statistics')
svector = 0.2 * np.ones(n)
print ('  sigma vector :',svector)

vlogpdf = edam.obserror.logpdf(yvector, xvector, svector)
print ('  logpdf :',vlogpdf)

vcdf = edam.obserror.cdf(yvector, xvector, svector)
print ('  cdf :',vcdf)

vsample = edam.obserror.sample(xvector,svector)
print ('  sample :',vsample)

