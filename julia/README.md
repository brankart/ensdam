
## Julia interface to EnsDAM

### Warning: This interface is still experimental.

Only a few simple functions have been tested from julia (see the list below).

To access to more functions, another possibility
is to import the python modules in julia:

import Pkg; Pkg.add("PyCall")

using PyCall

@pyimport ensdam

### List of available functions in ensdam.jl

```
ensdam
o       ensdam.cdf_gaussian : compute cdf of the Gaussian distribution N(0,1)
o       ensdam.cdf_gamma : compute cdf of the Gamma distribution Gamma(k,1)
o       ensdam.cdf_beta : compute cdf of the Beta distribution Beta(a,b)
o       ensdam.pdf_gaussian : compute pdf of the Gaussian distribution N(0,1)
o       ensdam.pdf_gamma : compute pdf of the Gamma distribution Gamma(k,1)
o       ensdam.pdf_beta : compute pdf of the Beta distribution Beta(a,b)
o       ensdam.invcdf_gaussian : compute inverse cdf of the Gaussian distribution N(0,1)
o       ensdam.invcdf_gamma : compute inverse cdf of the Gamma distribution Gamma(k,1)
o       ensdam.invcdf_beta : compute inverse cdf of the Beta distribution Beta(a,b)
o       ensdam.logpdf_gaussian : compute the log of the pdf of the Gaussian distribution N(0,1)
o       ensdam.logpdf_gamma : compute the log of the pdf of the Gamma distribution Gamma(k,1)
o       ensdam.logpdf_beta : compute the log of the pdf of the Beta distribution Beta(a,b)
```
