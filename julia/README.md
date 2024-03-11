
## Julia interface to EnsDAM

### Warning: This interface is still experimental.

First prepare the shared library with C-callable fortran routines
to be used by the Julia interface.

Only a few simple functions have been implemented julia (see the list below).
To import then in julia:

```julia
push!(LOAD_PATH, "./")
import ensdam
```

### List of available functions

In the probabilty module (scores.jl):

```
rank_histogram : compute rank histogram
crps_score_global : compute CRPS score
rcrv_score_global : compute RCRV score
```

In the probabilty module (probability.jl):

```
kiss_uniform : sample from uniform distribution
kiss_gaussian : sample from Gaussian distribution
kiss_gamma : sample from gamma distribution
kiss_beta : sample from beta distribution
cdf_gaussian : compute cdf of the Gaussian distribution N(0,1)
cdf_gamma : compute cdf of the Gamma distribution Gamma(k,1)
cdf_beta : compute cdf of the Beta distribution Beta(a,b)
pdf_gaussian : compute pdf of the Gaussian distribution N(0,1)
pdf_gamma : compute pdf of the Gamma distribution Gamma(k,1)
pdf_beta : compute pdf of the Beta distribution Beta(a,b)
invcdf_gaussian : compute inverse cdf of the Gaussian distribution N(0,1)
invcdf_gamma : compute inverse cdf of the Gamma distribution Gamma(k,1)
invcdf_beta : compute inverse cdf of the Beta distribution Beta(a,b)
logpdf_gaussian : compute the log of the pdf of the Gaussian distribution N(0,1)
logpdf_gamma : compute the log of the pdf of the Gamma distribution Gamma(k,1)
logpdf_beta : compute the log of the pdf of the Beta distribution Beta(a,b)
```
