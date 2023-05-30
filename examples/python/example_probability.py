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

# 1. Example with a normal distribution
print('-------------------------------------')
print('1. Example with a normal distribution')
print('-------------------------------------')
edam.probability.type = 'normal'
print('Type of observation error:',edam.probability.type)

# 1.a Scalar input
print('1.a Scalar input')
x = 1.       # Value of the random variable
print('  x=',x)

pdf = edam.probability.pdf(x)
print ('  pdf :',pdf)
logpdf = edam.probability.logpdf(x)
print ('  logpdf :',logpdf)
cdf = edam.probability.cdf(x)
print ('  cdf :',cdf)
invcdf = edam.probability.invcdf(cdf)
print ('  invcdf :',invcdf)

# 1.b Vector input
print('1.b Vector input')
x0 = 0. ; dx = 0.5 ; nx = 5
x = np.arange(x0, x0 + nx * dx, dx, dtype=np.float64)
print('  x=',x)

pdf = edam.probability.pdf(x)
print ('  pdf :',pdf)
logpdf = edam.probability.logpdf(x)
print ('  logpdf :',logpdf)
cdf = edam.probability.cdf(x)
print ('  cdf :',cdf)
invcdf = edam.probability.invcdf(cdf)
print ('  invcdf :',invcdf)


# 2. Example with a normal distribution
print('------------------------------------')
print('2. Example with a gamma distribution')
print('------------------------------------')
edam.probability.type = 'gamma'
edam.probability.k = 1.
edam.probability.theta = 1.
print('Type of observation error:',edam.probability.type)

# 2.a Scalar input
print('2.a Scalar input')
x = 1.       # Value of the random variable
print('  x=',x)

pdf = edam.probability.pdf(x)
print ('  pdf :',pdf)
logpdf = edam.probability.logpdf(x)
print ('  logpdf :',logpdf)
cdf = edam.probability.cdf(x)
print ('  cdf :',cdf)
invcdf = edam.probability.invcdf(cdf)
print ('  invcdf :',invcdf)

# 2.b Vector input
print('2.b Vector input')
x0 = 0.1 ; dx = 0.5 ; nx = 5
x = np.arange(x0, x0 + nx * dx, dx, dtype=np.float64)
print('  x=',x)

pdf = edam.probability.pdf(x)
print ('  pdf :',pdf)
logpdf = edam.probability.logpdf(x)
print ('  logpdf :',logpdf)
cdf = edam.probability.cdf(x)
print ('  cdf :',cdf)
invcdf = edam.probability.invcdf(cdf)
print ('  invcdf :',invcdf)

# 3. Example with a beta distribution
print('-----------------------------------')
print('3. Example with a beta distribution')
print('-----------------------------------')
edam.probability.type = 'beta'
edam.probability.k = 2.
edam.probability.theta = 2.
print('Type of observation error:',edam.probability.type)

# 3.a Scalar input
print('3.a Scalar input')
x = 0.5       # Value of the random variable
print('  x=',x)

pdf = edam.probability.pdf(x)
print ('  pdf :',pdf)
logpdf = edam.probability.logpdf(x)
print ('  logpdf :',logpdf)
cdf = edam.probability.cdf(x)
print ('  cdf :',cdf)
invcdf = edam.probability.invcdf(cdf)
print ('  invcdf :',invcdf)

# 3.b Vector input
print('3.b Vector input')
x0 = 0.1 ; dx = 0.2 ; nx = 5
x = np.arange(x0, x0 + nx * dx, dx, dtype=np.float64)
print('  x=',x)

pdf = edam.probability.pdf(x)
print ('  pdf :',pdf)
logpdf = edam.probability.logpdf(x)
print ('  logpdf :',logpdf)
cdf = edam.probability.cdf(x)
print ('  cdf :',cdf)
invcdf = edam.probability.invcdf(cdf)
print ('  invcdf :',invcdf)

