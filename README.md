## EnsDAM : Ensemble Data Assimilation Modules

EnsDAM is a collection of FORTRAN modules
that can be useful to ensemble data assimilation systems.

EnsDAM is distributed under the terms of the CeCILL free software license agreement.
See LICENSE.txt for more information. 

### Installation of EnsDAM

You need a FORTRAN-90 compiler and the NetCDF library (with f90 support) installed.

To compile the library and the examples :

- create a 'make.macro' file corresponding to your compiler in the 'macro' directory.
  This is the Makefile configurable part, which specifies
  your compiler options and where to find the NetCDF library.

- edit the file 'build/Makefile.template' to include this 'make.macro' file (first line below the title)

- compile (library and examples) with:

```bash
cd build
make
make examples
```

- if everything goes well, the EnsDAM library 'libensdam.a'
  has been created in the 'lib' directory, the module files (*.mod)
  have been created in the 'include' directory' and the example
  exexutables (*.x) have been created in the 'examples' directory.

 To check the installation :

 - run the examples in the 'example' directory:

```bash
cd examples
./random_field_on_the_sphere.x
./low_pass_filter_on_the_sphere.x
./observation_regression_on_the_sphere.x
```

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
EnsScores :
  module ensdam_score_crps : compute CRPS score
  module ensdam_score_rcrv : compute RCRV score
  module ensdam_score_optimality : compute optimality score
  module ensdam_score_entropy : compute entropy score
EnsStat : ensemble statistics
  module ensdam_meanstd : compute/update ensemble mean and standard deviation
  module ensdam_covariance : compute/update ensemble covariance/correlation
EnsUpdate : ensemble observational update
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

