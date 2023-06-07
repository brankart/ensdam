## EnsDAM : Ensemble Data Assimilation Modules

EnsDAM is a collection of FORTRAN modules
that can be useful to ensemble data assimilation systems.

EnsDAM is distributed under the terms of the CeCILL free software license agreement.
See LICENSE.txt for more information. 

EnsDAM code is written in Fortran.
Interface to python is in development.

### Installation of the EnsDAM python module

Requirements: cmake, Fortran compiler, python

Edit and run the script: compile.bash

Run the examples in the 'example/python' directory.

The pyensdam module currently includes:
pyensdam.scores and pyensdam.obserror.

### Installation of the EnsDAM Fortran library only

To compile the Fortran library only and the Fortran examples,
go to the 'tools' directory and follow the compilation guidelines.

If everything goes well, the EnsDAM library 'libensdam.a'
has been created in the 'lib' directory, the module files (.mod)
have been created in the 'include' directory' and the example
exexutables (.x) have been created in the 'examples/fortran' directory.

Run the examples in the 'example/fortran' directory.

### List of available EnsDAM modules

See documentation (in the doc directory) for more details.

```
EnsAnam : ensemble anamorphosis transformation
  module ensdam_anaqua : computation of ensemble quantiles
  module ensdam_anatra : application of the anamorphosis transformation
  module ensdam_anaobs : transformation of observations
  module ensdam_anautil : utilities for anamorphosis transformation
EnsAugm : ensemble augmentation
  module ensdam_ensaugm : sample augmented ensemble
  module ensdam_schurprod : compute Schur products
EnsScores : ensemble scores
  module ensdam_score_crps : compute CRPS score
  module ensdam_score_rcrv : compute RCRV score
  module ensdam_score_optimality : compute optimality score
  module ensdam_score_entropy : compute entropy score
EnsStat : ensemble statistics
  module ensdam_meanstd : compute/update ensemble mean and standard deviation
  module ensdam_covariance : compute/update ensemble covariance/correlation
EnsUpdate : ensemble observational update
  module ensdam_mcmc_update : MCMC ensemble observational update
InterpTools : interpolation tools
  module ensdam_interp : locate/interpolate in 1D/2D grids (cartesian or spherical coordinates)
ObsError : observation error
  module ensdam_obserror : compute/sample observation error probability distribution
StochTools : stochastic tools
  module ensdam_storng : sample from uniform/normal/gamma/beta distributions
  module ensdam_stotge : sample from truncated exponential/normal distributions
  module ensdam_stoutil : compute pdf/cdf/logpdf/invcdf of normal/gamma/beta distributions
  module ensdam_stoanam : transform normal number into gamma/beta/normal_product number
  module ensdam_storfg : generate 2D random field with specified power spectrum (on the sphere)
TranSpHO : scale separation (by projection on the spherical harmonics)
  module ensdam_spharea : compute spherical areas
  module ensdam_sphylm : project on spherical harmonics / transform back on the sphere
```

