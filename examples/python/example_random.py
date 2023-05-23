import pyensdam as edam
import numpy as np

# Example illustrating the module: random
# =======================================

# The module contains the functions:
# - random.seed: Seed random number generator
# - random.seed_save: Save seed in restart file
# - random.seed_load: Load seed from restart file
# - random.check: Check random number generator
# - random.swap: Random array swapping
# - random.uniform: Draw random numbers with uniform distribution
# - random.normal: Draw random numbers with normal distribution
# - random.exp: Draw random numbers with exponential distribution
# - random.gamma: Draw random numbers with gamma distribution
# - random.beta: Draw random numbers with beta distribution
# - random.truncated_exp: Draw random numbers with truncated exponential distribution
# - random.truncated_normal: Draw random numbers with truncated normal distribution
# - random.truncated_normal_vec: Draw random vectors with truncated normal distribution
# - random.field1d_init: Initialization for the sampling of 1D random fields
# - random.field2d_init: Initialization for the sampling of 2D random fields
# - random.field1d_sample: Sample 1D random fields with given spectrum
# - random.field2d_sample: Sample 2D random fields with given spectrum

# The module contains the parameters:
# - 

# Notes:
# -

#edam.random.seed_save()

#edam.random.seed(0)

#edam.random.seed(1)

#edam.random.seed_save()

zran=edam.random.uniform([5])
print ('  Uniform random number:',zran)

zran=edam.random.normal([5])
print ('  Normal random number:',zran)

zran=edam.random.exp([5])
print ('  Exponential random number:',zran)

zran=edam.random.gamma(1.,[5])
print ('  Gamma random number:',zran)

zran=edam.random.beta(1.,1.,[5])
print ('  Beta random number:',zran)

zran=edam.random.truncated_exp(1.,[5])
print ('  Truncated exponential random number:',zran)

zran=edam.random.truncated_normal(1.,2.,[5])
print ('  Truncated normal random number:',zran)

a = np.arange(10,dtype=np.intc)
print ('  Ordered array:',a)
edam.random.swap(a)
print ('  Randomly swapped array:',a)

# Define constraint (-x-y <= 0)
A = np.ones((1,2)) - 2  # A = [ -1, -1 ]
b = np.zeros(1)         # b = 0
# Sample random vectors with truncated Gaussian distribution
sample = edam.random.truncated_normal_vec(100,A,b);
for i in range(9,99,10):
  print ('  Truncated normal random vector (x+y>0):',sample[:,i])

