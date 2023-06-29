using PyCall
@pyimport pyensdam as edam
@pyimport numpy as np

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


# Seed the random number generator
# --------------------------------

println("--------------------------------")
println("Seed the random number generator")
println("--------------------------------")

# Perform short check of the random number generator (silent, but stop if failed)
edam.random.check("short")
# Perform long check of the random number generator (with message, and stop if failed)
edam.random.check("long")  # comment out to save time

# Use default seed
edam.random.seed(0)

# Save current state in restart file
edam.random.seed_save()

# Use seed index 8
edam.random.seed(8)

# Load seed from restart file (i.e. go back to default seed)
edam.random.seed_load()

# Generate random numbers with various distributions
# --------------------------------------------------

println("--------------------------------------------------")
println("Generate random numbers with various distributions")
println("--------------------------------------------------")

zran=edam.random.uniform([5])

println("  Uniform random number: $zran")

zran=edam.random.normal([5])
println("  Normal random number: $zran")

zran=edam.random.exp([5])
println("  Exponential random number: $zran")

zran=edam.random.gamma(1.,[5])
println("  Gamma random number: $zran")

zran=edam.random.beta(1.,1.,[5])
println("  Beta random number: $zran")

zran=edam.random.truncated_exp(1.,[5])
println("  Truncated exponential random number: $zran")

zran=edam.random.truncated_normal(1.,2.,[5])
println("  Truncated normal random number: $zran")

# Randomly swap elements of input array
# -------------------------------------

println("-------------------------------------")
println("Randomly swap elements of input array")
println("-------------------------------------")

a = collect(Int32, 0:9)
println("  Ordered array: $a")
edam.random.swap(a)
println("  Randomly swapped array: $a")

# Sample random vector with multivariate truncated Gaussian distribution
# ----------------------------------------------------------------------

println("----------------------------------------------------------------------")
println("Sample random vector with multivariate truncated Gaussian distribution")
println("----------------------------------------------------------------------")

# Define constraint (Ax <= b, here: -x-y <= 0 i.e. x+y >= 0)
# Just add rows to A and B to impose more inequality constraints
A = - np.ones((1,2))    # A = [ -1, -1 ]
b = np.zeros(1)         # b = 0
# Sample random vectors with truncated Gaussian distribution
# All iterates of the Gibbs sampler are output
# -> drop first draws and subsample the sequence of draws as appropriate
sample = edam.random.truncated_normal_vec(100,A,b);
for i in range(9, step=10, stop=99)
  println("  Truncated normal random vector (x+y>0):",sample[:,i])
end

# Sample 1D random field with given continuous spectrum
# -----------------------------------------------------

println("-----------------------------------------------------")
println("Sample 1D random field with given continuous spectrum")
println("-----------------------------------------------------")

# Definition of the spectrum (discretization of a continuous spectrum)
f0 = 0.1 ; df = 0.1 ; nf = 10
spct_freq = np.arange(f0, f0 + nf * df, df, dtype=np.double)
spct_power = np.ones(10)

# Initialization of the sampler
edam.random.field1d_init(spct_freq,spct_power)

# Definition of the output grid
x0 = 0. ; dx = 0.5 ; nx = 100
x = np.arange(x0, x0 + nx * dx, dx, dtype=np.double)

# Generate random field
nharm = 100 # number of harmonics to superpose (sampled from the continuous spectrum)
field1d = edam.random.field1d_sample(x,nharm)

println("  1D random field:",field1d[:])
