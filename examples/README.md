## EnsDAM : examples

The main example is `mcmc_ensemble_update.F90`,
which corresponds to the example application described in the paper:

*Brankart J.-M., 2019: 
Implicitly localized ensemble observational update
to cope with nonlocal/nonlinear data constraints
in large-size inverse problems,
submitted to Frontiers in Applied Mathematics and Statistics.*

This main example illustrates many of the EnsDAM modules.
The other examples are mainly dedicated
to illustrating the TranSpHO modules.

### Example: mcmc_ensemble_update

The parameter file is `mcmc_ensemble_update_parameters.h90`,
(pointing to `mcmc_ensemble_update_obs100err20loc06.h90`).
This corresponds to the parameters used in the reference example
described in the paper.

Two other parameter files are provides:
`mcmc_ensemble_update_obs100err20loc06x2.h90` and
`mcmc_ensemble_update_obs100err20loc06x4.h90`,
which correspond to multiplying the grid resolution
by 2 and 4 respectively.

All experiments described in the paper
can be reproduced by changing only the parameter file.

To run the experiments in parallel (with MPI), the number of processors must be tuned
to have at least one observation in the domain of each processor.
With this constraint, the reference experiment can be run
with 64 processors, and the higher resolution grids,
with 256 and 1024 processors respectively.
More processors can be used if the number of observations
is increased in the parameter file.

### Example: random_field_on_the_sphere

Generation of a random field with specified spectrum on the sphere.

This illustrates the routine *gen_field_2s* from the StochTools modules.

### Example: low_pass_filter_on_the_sphere

Apply low pass filter to a field defined on the sphere.

This illustrates the routines *proj_ylm* and *back_ylm* from the TranSpHO modules.

### Example: observation_regression_on_the_sphere

Perform regression of observations in the basis of the spherical harmonics.

This illustrates the routine *regr_ylm* from the TranSpHO modules.
