push!(LOAD_PATH, "../../julia/")
import scores

using Random

# Example illustrating the module: scores
# =======================================

# Parameters of the example
n = 1000    # number of variables
m = 50      # size of the ensemble

mu = 0.0     # Ensemble mean
sigma = 1.0  # Ensemble standard deviation

# Sample reference truth from N(0,I) distribution
reference_truth = randn(n) .* sigma .+ mu

# Sample ensemble from N(0,I) distribution
ensemble = randn(n, m) .* sigma .+ mu

# 1. Computation of rank histograms
println("---------------------------------")
println("1. Computation of rank histograms")
println("---------------------------------")

# Compute rank histogram
ranks, rank_histogram = scores.rank_histogram(ensemble,reference_truth)
println("Rank histogram:",rank_histogram)

# 2. Computation of CRPS score
println("----------------------------")
println("2. Computation of CRPS score")
println("----------------------------")

crps,crps_reliability,crps_resolution = scores.crps_score_global(ensemble,reference_truth)
println("  CRPS total:       ",crps)
println("  CRPS reliability: ",crps_reliability)
println("  CRPS resolution:  ",crps_resolution)

# 3. Computation of RCRV score
println("----------------------------")
println("3. Computation of RCRV score")
println("----------------------------")

rcrv_bias,rcrv_spread = scores.rcrv_score_global(ensemble,reference_truth)
println("  RCRV bias:   ",rcrv_bias)
println("  RCRV spread: ",rcrv_spread)

