"""
module scores: Ensemble probablistic scores
===========================================

Available functions:
 -  scores.rank_histogram: Compute rank histogram
 -  scores.crps : Compute CRPS score (total, reliability, resolution)
 -  scores.rcrv : Compute RCRV score (bias, spread)

"""
module scores

const libdir = dirname(pathof(@__MODULE__))
const lib = joinpath(libdir,"ensdam")

# Computation of rank histograms
# ------------------------------

# Function to compute ranks and rank histograms
"""
    ranks_histogram(ens::Array{Cdouble, 2}, verif::Array{Cdouble, 1})

Compute ranks and rank histogram

Arguments:
- `ens::Array{Cdouble, 2}`  : ensemble simulation (nvar,nens)
- `verif::Array{Cdouble, 1}`: verification data (nvar)

Returns:
- `ranks`: Array of ranks.
- `rank_histogram`: Histogram of ranks.
"""
function rank_histogram(ens::Array{Cdouble,2}, verif::Array{Cdouble,1})
    n=size(ens,1)
    m=size(ens,2)
    if !(n==length(verif)) println("rank_histogram: incoherent sizes, not calling library"); end

    ranks = zeros(Cint, n)            # Pre-allocate array for output
    rank_histogram = zeros(Cint, m+1) # Pre-allocate array for output

    j_compute_ranks_histogram(Cint(n), Cint(m), ens, verif, ranks, rank_histogram)

    (ranks, rank_histogram)
end

# Interface to corresponding C-callable fortran function
function j_compute_ranks_histogram(nstate::Cint, nens::Cint, ens::Array{Cdouble, 2}, verif::Array{Cdouble, 1},
                                   ranks::Array{Cint, 1}, rank_histogram::Array{Cint, 1})
    ccall((:c_compute_ranks_histogram,lib), Cvoid,
          (Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}),
          nstate, nens, ens, verif, ranks, rank_histogram)
end

# Computation of CRPS score
# -------------------------

"""
    crps_score_global(ens::Array{Cdouble, 2}, verif::Array{Cdouble, 1})

Calculate CRPS (Continuous Ranked Probability Score) score, reliability, and resolution.

Arguments:
- `ens::Array{Cdouble, 2}`  : ensemble simulation (nvar,nens)
- `verif::Array{Cdouble, 1}`: verification data (nvar)

Returns a tuple:
- `crps::Float64`: CRPS score.
- `reliability::Float64`: Reliability.
- `resolution::Float64`: Resolution.
"""
function crps_score_global(ens::Array{Cdouble,2}, verif::Array{Cdouble,1})
    n=size(ens,1)
    m=size(ens,2)
    if !(n==length(verif)) println("crps_score_global: incoherent sizes, not calling library"); end

    crps = Ref{Cdouble}()
    reliability = Ref{Cdouble}()
    resolution = Ref{Cdouble}()

    ccall((:c_crps_score_global,lib), Cvoid,
          (Cint, Cint, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          Cint(n), Cint(m), crps, reliability, resolution, ens, verif)

    return crps[], reliability[], resolution[]
end

# Computation of RCRV score
# -------------------------
"""
    rcrv_score_global(ens::Array{Cdouble, 2}, verif::Array{Cdouble, 1})

Compute RCRV score  (bias and spread components)

This function computes the bias and spread of ensemble using the ensemble data (`ens`) and verification data (`verif`).

Arguments:
- `ens::Array{Cdouble, 2}`  : ensemble simulation (nvar,nens)
- `verif::Array{Cdouble, 1}`: verification data (nvar)

Returns a tuple:
- `ens_bias::Float64`: bias component of RCRV score
- `ens_spread::Float64`: spread component of RCRV score
"""
function rcrv_score_global(ens::Array{Cdouble, 2}, verif::Array{Cdouble, 1})
    n = size(ens, 1)
    m = size(ens, 2)

    ens_bias = Ref{Cdouble}()
    ens_spread = Ref{Cdouble}()

    ccall((:c_rcrv_score_global,lib), Cvoid,
          (Cint, Cint, Ref{Cdouble}, Ref{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          Cint(n), Cint(m), ens_bias, ens_spread, ens, verif)

    return ens_bias[], ens_spread[]
end

# end of module
end
