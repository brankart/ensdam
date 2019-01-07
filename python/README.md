
Python interface to EnsDAM

Warning: This interface is still experimental.
Many modules have not been tested from python (see the list below).

To generate the Python interface (with f2py) :

 - edit the 'make.macro' file corresponding to your compiler in the 'macro' directory.
   This is the Makefile configurable part, which specifies options to pass to f2py.

 - edit the Makefile to include this 'make.macro' file (first line below the title)

 - compile with "make" (gmake)

 - if everything goes well, the EnsDAM shared librares should
   have been created in the 'lib' directory,
   ready to be imported in the python interface module (ensdam.py).

To import EnsDAM in python:

 - import ensdam

List of available EnsDAM modules
  x=interface not yet available in ensdam.py, try importing directly the .so file generated bu f2py
  u=interface is still untested
  p=interface is only partially tested
  b=interface has bugs
  o=interface has been tested, to use with care

ensdam
    ensdam.ensanam : ensemble anamorphosis transformation
x       ensdam.ensanam.ens_quantiles : computation of ensemble quantiles
x       ensdam.ensanam.ana_forward : forward transformation
x       ensdam.ensanam.ana_backward : backward transformation
x       ensdam.ensanam.ana_obs : transformation of observations
x       ensdam.ensanam.ana_obs_sym : transformation of observations (symmetric case)
x       ensdam.ensanam.ana_util_quaref : compute the quantiles of the target distribution
    ensdam.ensaugm : ensemble augmentation
x       ensdam.ensaugm.sample_augmented_ensemble : sample augmented ensemble
x       ensdam.ensaugm.newproduct : sample new Schur product
x       ensdam.ensaugm.getproduct : get specified Schur product
    ensdam.ensscores : ensemble scores
x       ensdam.ensscores.crps_score : compute CRPS score (with option to partition the data)
x       ensdam.ensscores.crps_cumul : accumulate data to prepare the final computation of the CRPS score
x       ensdam.ensscores.crps_final : compute final score from accumulated data
x       ensdam.ensscores.rcrv_score : compute RCRV score (with option to partition the data)
x       ensdam.ensscores.rcrv_cumul : accumulate data to prepare the final computation of the RCRV score
x       ensdam.ensscores.optimality_score : compute optimality score (with option to partition the data)
x       ensdam.ensscores.optimality_cumul : accumulate data to prepare the final computation of the RCRV score
x       ensdam.ensscores.events_score : compute ensemble entropy score for the required events
x       ensdam.ensscores.events_relative_entropy : compute relative entropy between ensemble distribution and reference distribution
x       ensdam.ensscores.events_cross_entropy : compute cross entropy between ensemble distribution and reference distribution
x       ensdam.ensscores.events_entropy : compute entropy of ensemble distribution for the required events
x       ensdam.ensscores.events_probability : compute events marginal probability distributions from the ensemble
    ensdam.ensstat : ensemble statistics
p       ensdam.ensstat.ensemble_meanstd : compute mean and standard deviation from input ensemble
p       ensdam.ensstat.update_meanstd : update mean and standard deviation with one additional input member
u       ensdam.ensstat.ensemble_correlation : compute correlation from input ensemble
u       ensdam.ensstat.ensemble_representer : compute representer from input ensemble
u       ensdam.ensstat.ensemble_covariance : compute covariance from input ensemble
u       ensdam.ensstat.update_meancov : update mean and covariance with one additional input member
    ensdam.ensupdate ! ensemble observational update
    ensdam.interptools ! interpolation tools
x       ensdam.interptools.grid1D_locate : locate data point in 1D grid
x       ensdam.interptools.grid1D_interp : compute interpolation weights in 1D grid
x       ensdam.interptools.grid2D_init : initialize 2D grid
x       ensdam.interptools.grid2D_locate : locate data point in 2D grid
x       ensdam.interptools.grid2D_interp : compute interpolation weights in 2D grid
    ensdam.obserror ! observation error
x       ensdam.obserror.obserror_logpdf : compute logarithm of observation error probability density function
x       ensdam.obserror.obserror_cdf : compute cumulate distribution function for observation errors
x       ensdam.obserror.obserror_sample : sample probability distribution of observation errors
    ensdam.stochtools ! stochastic tools
o       ensdam.stochtools.kiss : 64-bit KISS random number generator (period ~ 2^250)
o       ensdam.stochtools.kiss_seed : define seeds for KISS random number generator
o       ensdam.stochtools.kiss_save : save current state of KISS (for future restart)
o       ensdam.stochtools.kiss_load : load the saved state of KISS
o       ensdam.stochtools.kiss_reset : reset the default seeds
o       ensdam.stochtools.kiss_check : check the KISS pseudo-random sequence
o       ensdam.stochtools.kiss_uniform : real random numbers with uniform distribution in [0,1]
o       ensdam.stochtools.kiss_gaussian : real random numbers with Gaussian distribution N(0,1)
o       ensdam.stochtools.kiss_gamma : real random numbers with Gamma distribution Gamma(k,1)
o       ensdam.stochtools.kiss_beta : real random numbers with Beta distribution Beta(a,b)
o       ensdam.stochtools.kiss_sample : select a random sample from a set of integers
o       ensdam.stochtools.ran_te : sample random number with truncated exponential distribution
o       ensdam.stochtools.ran_tg : sample random number with truncated Gaussian distribution
o       ensdam.stochtools.ranv_tg : sample random vector with truncated Gaussian distribution
o       ensdam.stochtools.cdf_gaussian : compute cdf of the Gaussian distribution N(0,1)
o       ensdam.stochtools.cdf_gamma : compute cdf of the Gamma distribution Gamma(k,1)
o       ensdam.stochtools.cdf_beta : compute cdf of the Beta distribution Beta(a,b)
o       ensdam.stochtools.pdf_gaussian : compute pdf of the Gaussian distribution N(0,1)
o       ensdam.stochtools.pdf_gamma : compute pdf of the Gamma distribution Gamma(k,1)
o       ensdam.stochtools.pdf_beta : compute pdf of the Beta distribution Beta(a,b)
o       ensdam.stochtools.invcdf_gaussian : compute inverse cdf of the Gaussian distribution N(0,1)
o       ensdam.stochtools.invcdf_gamma : compute inverse cdf of the Gamma distribution Gamma(k,1)
o       ensdam.stochtools.invcdf_beta : compute inverse cdf of the Beta distribution Beta(a,b)
o       ensdam.stochtools.logpdf_gaussian : compute the log of the pdf of the Gaussian distribution N(0,1)
o       ensdam.stochtools.logpdf_gamma : compute the log of the pdf of the Gamma distribution Gamma(k,1)
o       ensdam.stochtools.logpdf_beta : compute the log of the pdf of the Beta distribution Beta(a,b)
o       ensdam.stochtools.gau_to_gam : transform Gaussian number into gamma number
o       ensdam.stochtools.gau_to_beta : transform Gaussian number into beta number
o       ensdam.stochtools.gprod_to_gau : transform Gaussian product into Gaussian number
u       ensdam.stochtools.gen_field_2s : generate 2D random field with specified power spectrum (on the sphere)
    ensdam.transpho ! scale separation (by projection on the spherical harmonics)
x       ensdam.transpho.init_ylm : initialize computation of spherical harmonics
x       ensdam.transpho.init_regr_ylm : initialize regression of observations
x       ensdam.transpho.proj_ylm : project on spherical harmonics
x       ensdam.transpho.back_ylm : transform back on the sphere
x       ensdam.transpho.regr_ylm : regression of observations on spherical harmonics
x       ensdam.transpho.disp_ylm : output one single spherical harmonics
x       ensdam.transpho.ylm : evaluate spherical harmonics
