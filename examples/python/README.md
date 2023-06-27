
pyensdam: Ensemble data assimilation modules
============================================

Avalaible modules:
 - pyensdam.anamorphosis : ensemble anamorphosis transformation
 - pyensdam.augment : ensemble augmentation with MCMC sampler
 - pyensdam.interpolation : interpolation in 1D and 2D grids
 - pyensdam.obserror : observation error (normal, lognormal, gamma, beta)
 - pyensdam.probability : probability distribution (normal, lognormal, gamma, beta)
 - pyensdam.random : random field generator
 - pyensdam.scores : ensemble probabilistic scores
 - pyensdam.statistics : ensemble statistics (mean, std, correlation, representer)
 - pyensdam.transpho : transformation in the basis of the spherical harmonics
 - pyensdam.update : ensemble observational update, with an MCMC sampler



pyensdam.anamorphosis: ensemble anamorphosis
============================================

Available functions:
 -  anamorphosis.quantiles : compute ensemble quantiles
 -  anamorphosis.forward : forward anamorphosis transformation
 -  anamorphosis.backward : backward anamorphosis transformation
 -  anamorphosis.forward_obs : forward anamorphosis transformation of observations
 -  anamorphosis.forward_obs_sym : forward anamorphosis transformation of observations (symmetric)

Module parameters:
 -  anamorphosis.target : target probability distribution to use (default=normal, uniform, gamma, beta, custom)
 -  anamorphosis.quaref : quantiles of the target probability distribution to use (if target == 'custom'))
 -  anamorphosis.obstype : probability distribution of observations (default=normal, gamma, beta)


forward(var,qua,[rank])

       Forward anamorphosis transformation

       Inputs
       ------
       var [scalar, rank-1 or rank-2 double array] : variable(s) to transform (nvar) or (nens,nvar)
       qua [rank-1 or rank-2 double array] : ensemble quantiles (nqua) or (nqua,nvar)
       rank [scalar or rank-1 double array] : rank to use to resolve probabiity concentration

       Returns
       -------
       var [scalar, rank-1 or rank-2 double array] : transformed input data (in-place)

    
backward(var,qua)

       Backward anamorphosis transformation

       Inputs
       ------
       var [scalar, rank-1 or rank-2 double array] : variable(s) to transform (nvar) or (nens,nvar)
       qua [rank-1 or rank-2 double array] : ensemble quantiles (nqua) or (nqua,nvar)

       Returns
       -------
       var [scalar, rank-1 or rank-2 double array] : transformed input data (in-place)

    
anaobs = forward_obs(nsmp,obs,obs_std,obsens,quadef)

       Forward anamorphosis transformation of observations

       Inputs
       ------
       nsmp [integer] : size of transformed observation sample
       obs [rank-1 double array] : observations to transform (nobs)
       obs_std [rank-1 double array] : observation error standard deviation
       obsens [rank-2 double array] : ensemble equivalent to observations (nens,nobs)
       quadef [rank-1 double array] : definition of the quantiles used in anamorphosis (nqua)

       Returns
       -------
       anaobs [rank-2 double array] : sample of transformed observations (nsmp,nobs)

    
anaobs = forward_obs_sym(nsmp,obs,obs_std,obsqua)

       Forward anamorphosis transformation of observations (symmetric)

       Inputs
       ------
       nsmp [integer] : size of transformed observation sample
       obs [rank-1 double array] : observations to transform (nobs)
       obs_std [rank-1 double array] : observation error standard deviation
       obsqua [rank-1 double array] : quantiles of the ensemble equivalent to observations (nqua,nobs)

       Returns
       -------
       anaobs [rank-2 double array] : sample of transformed observations (nsmp,nobs)

    

pyensdam.augment: ensemble augmentation
=======================================

Available functions:
 -  augment.sample_mcmc : resample input ensemble with MCMC sampler,
                          using covariance localization


augens = sample_mcmc(ens,multiplicity,naug,maxchain)

       Resample input ensemble with MCMC sampler, using covariance localization

       Inputs
       ------
       ens [rank-3 double array] : multiscale input ensemble (nscl,nens,nvar)
       multiplicity [rank-1 integer array] : multiplicity of each scale in the Schur products
       naug [integer] : number of members requested for the augmented ensemble
       maxchain [integer] : number of iteration of the MCMC sampler

       Returns
       -------
       augens [rank-2 double array] : augmented ensemble (naug,nvar)
    

pyensdam.interpolation: interpolation tools
===========================================

Available functions:
 -  interpolation.locate1D : locate positions in 1D grid and compute interpolation weights
 -  interpolation.interp1D : apply interpolation on 1D input field
 -  interpolation.define2D : define 2D grid
 -  interpolation.locate2D : locate positions in 2D grid and compute interpolation weights
 -  interpolation.interp2D : apply interpolation on 2D input field
 -  interpolation.unmask2D : unmask input 2D array

Available parameters:
 -  unmask_spval   : special value to unmask
 -  unmask_max     : maximum number of iteration for unmasking
 -  unmask_window  : maximum averaging window in unmasking
 -  unmask_k_ew    : grid periodicity in unmasking
 -  unmask_damping : damp extrapolation as a function of distance to coast (0 = no damping)


location, weight = locate1D(grid,x)

       Locate positions in 1D grid and compute interpolation weights

       Inputs
       ------
       grid  [rank-1 double array] : definition of the grid locations (in ascending order (ngrid)
       x [double array] : array of positions to locate in the grid

       Returns
       -------
       location [integer array] : index of the grid cell where postions are located (same shape as x)
       weight [double array] : interpolation weight to use (same shape as x)

    
field_interpolated = interp1D(field,location,weight)

       Apply interpolation on 1D input field

       Inputs
       ------
       field [rank-1 double array] : field from which interpolating (same dimensions are as grid)
       location [integer array] : index of the grid cell where interpolated postions are located
       weight [double array] : interpolation weight to use (same shape as location)

       Returns
       -------
       field_interpolated [double array] : interpolation weight to use (same shape as location)

    
define2D(xgrid,ygrid,[grid_type])

       Define 2D grid

       Inputs
       ------
       xgrid  [rank-2 double array] : definition of the x-coordinates
       ygrid  [rank-2 double array] : definition of the x-coordinates
       grid_type: typr of coordinates (default='cartesian', spherical)

    
location, weight = locate2D(x,y)

       Locate positions in 2D grid and compute interpolation weights

       Inputs
       ------
       x [double array] : array of positions to locate in the grid (x-coordinate)
       y [double array] : array of positions to locate in the grid (y-coordinate)

       Returns
       -------
       location [integer array] : index of the grid cell where postions are located
                                  (same shape as x and y, for a 2-dimension vector)
       weight [double array] : interpolation weight to use
                               (same shape as x and y, for a '2 by 2' matrix)

    
field_interpolated = interp2D(field,location,weight)

       Apply interpolation on 2D input field

       Inputs
       ------
       field [rank-2 double array] : field from which interpolating (same dimensions are as grid)
       location [integer array] : index of the grid cell where interpolated postions are located
       weight [double array] : interpolation weight to use (same shape as location)

       Returns
       -------
       field_interpolated [double array] : interpolation weight to use (same shape as location)

    
unmask2D(field)

       Unmask 2D input field

       Inputs
       ------
       field [rank-2 double array] : field to unmask (nj,ni)

       Returns
       -------
       field [rank-2 double array] : unmask field (in place)

    

pyensdam.obserror: Operations related to observation error
==========================================================

Available functions:
 -  obserror.logpdf : Compute the logarithm of the observation error pdf
 -  obserror.cdf : Compute the observation error cdf
 -  obserror.sample : Sample the observation error probability distribution

Module parameters:
 -  obserror.obstype : Type of observation error
                       (normal, lognormal, gamma, beta)

Notes:
 - When applied to observation vectors, observation errors are assumed independent,
   only marginal distributions for each component are used.
   In logpdf, the contributions of the vector components are summed.
 - The random number generator used by this module is in the module: stochtools,
   it can be seeded using functions provided there.


logpdf = logpdf(y,x,sigma)

       Compute the logarithm of the observation error pdf
       Supported distributions (obstype module parameter): normal, lognormal, gamma, beta

       Inputs
       ------
       y [double scalar or vector] : value of the observation
       x [double scalar or vector] : expected value of the distribution
       sigma [double scalar or vector] : spread of the distribution
         * the meaning of sigma is not the same for all distributions:
            -normal: sigma = standard deviation
            -lognormal: sigma = standard deviation / expected value
            -gamma: sigma = standard deviation / expected value
            -beta: sigma = max standard deviation (for x=0.5)
         * sigma can be scalar when x and y are vectors
       Returns
       -------
       logpdf : [double scalar] : logarithm of the observation error pdf
    
cdf = cdf(y,x,sigma)

       Compute the observation error cdf
       Supported distributions (obstype module parameter): normal, lognormal, gamma, beta

       Inputs
       ------
       y [double scalar or vector] : value of the observation
       x [double scalar or vector] : expected value of the distribution
       sigma [double scalar or vector] : spread of the distribution
         * the meaning of sigma is not the same for all distributions:
            -normal: sigma = standard deviation
            -lognormal: sigma = standard deviation / expected value
            -gamma: sigma = standard deviation / expected value
            -beta: sigma = max standard deviation (for x=0.5)
         * sigma can be scalar when x and y are vectors
       Returns
       -------
       cdf : [double scalar or vector] : value of the cdf for each observation
    
sample = sample(x,sigma)

       Sample the observation error probability distribution
       Supported distributions (obstype module parameter): normal, lognormal, gamma, beta

       Inputs
       ----------
       x [double scalar or vector] : expected value of the distribution
       sigma [double scalar or vector] : spread of the distribution
         * the meaning of sigma is not the same for all distributions:
            -normal: sigma = standard deviation
            -lognormal: sigma = standard deviation / expected value
            -gamma: sigma = standard deviation / expected value
            -beta: sigma = max standard deviation (for x=0.5)
         * sigma can be scalar when x and y are vectors
       Returns
       -------
       sample : [double scalar or vector] : draw from the observation error pdf
       

pyensdam.probability: probability distribution
==============================================

Available functions:
 -  probability.pdf: compute the probability density function
 -  probability.logpdf: compute the logartihm of the probability density function
 -  probability.cdf: compute the cumulative distribution function
 -  probability.invcdf: compute the iinverse cumulative distribution function

Module parameters:
 -  probability.type: type of probability distribution (normal, gamma, beta)
 -  probability.k: shape parameter of the gamma distribution
 -  probability.theta: scale parameter of the gamma distribution
 -  probability.a: parameter alpha of the beta distribution
 -  probability.b: parameter beta of the beta distribution


pdf=pdf(x)

       Compute the probability density function

       Input
       -----
       x [double]: value of the random variable

       Output
       ------
       pdf [double]: value of the pdf function

    
logpdf=logpdf(x)

       Compute the log of probability density function

       Input
       -----
       x [double]: value of the random variable

       Output
       ------
       logpdf [double]: value of the log of the pdf function

    
cdf=cdf(x)

       Compute the cumulative distribution function

       Input
       -----
       x [double]: value of the random variable

       Output
       ------
       cdf [double]: value of the cumulative distribution function

    
invcdf=invcdf(rank)

       Compute the inverse cumulative distribution function

       Input
       -----
       rank [double]: value of the random variable

       Output
       ------
       invcdf [double]: value of the inverse cumulative distribution function

    

pyensdam.random: generate random numbers and fields
===================================================

Available functions:
 -  random.seed: Seed random number generator
 -  random.seed_save: Save seed in restart file
 -  random.seed_load: Load seed from restart file
 -  random.check: Check random number generator
 -  random.swap: Random array swapping
 -  random.uniform: Draw random numbers with uniform distribution
 -  random.normal: Draw random numbers with normal distribution
 -  random.exp: Draw random numbers with exponential distribution
 -  random.gamma: Draw random numbers with gamma distribution
 -  random.beta: Draw random numbers with beta distribution
 -  random.truncated_exp: Draw random numbers with truncated exponential distribution
 -  random.truncated_normal: Draw random numbers with truncated normal distribution
 -  random.truncated_normal_vec: Draw random vectors with truncated normal distribution
 -  random.field1d_init: Initialization for the sampling of 1D random fields
 -  random.field2d_init: Initialization for the sampling of 2D random fields
 -  random.field1d_sample: Sample 1D random fields with given spectrum
 -  random.field2d_sample: Sample 2D random fields with given spectrum
 -  random.field2s_sample: Sample 2D random fields with given spectrum on the sphere

Module parameters:


seed(seed_idx)

       Seed random number generator

       Inputs
       ------
       seed_idx[integer] : index of seed to use
    
seed_save()

       Save current state of random number generator in restart file (.kiss_restart)
    
seed_load()

       Load seed from restart file (.kiss_restart)
    
check(check_type)

       Check random number generator

       Inputs
       ------
       check_type[string] : type of check ('short' or 'long')

    
swap(a,[k])

       Random array swapping

       Inputs
       ------
       a [integer vector] : input array to swap
       k [integer]: number of elements to sample from array (default=size of a)

    
zran=uniform([shape])

       Draw random numbers with uniform distribution

       Input
       -----
       shape: shape of numpy array to produce (None-> double scalar)

       Output
       ------
       zran [double]: scalar or array with independent random numbers

    
zran=normal([shape])

       Draw random numbers with normal distribution

       Input
       -----
       shape: shape of numpy array to produce (None-> double scalar)

       Output
       ------
       zran [double]: scalar or array with independent random numbers

    
zran=exp([shape])

       Draw random numbers with exponential distribution

       Input
       -----
       shape: shape of numpy array to produce (None-> double scalar)

       Output
       ------
       zran [double]: scalar or array with independent random numbers

    
zran=gamma(k,[shape])

       Draw random numbers with gamma distribution

       Input
       -----
       k [double]: parameter of gamma distribution
       shape: shape of numpy array to produce (None-> double scalar)

       Output
       ------
       zran [double]: scalar or array with independent random numbers

    
zran=beta(a,b,[shape])

       Draw random numbers with beta distribution

       Input
       -----
       a [double]: parameter of beta distribution
       b [double]: parameter of beta distribution
       shape: shape of numpy array to produce (None-> double scalar)

       Output
       ------
       zran [double]: scalar or array with independent random numbers

    
zran=truncated_exp(a,[shape])

       Draw random numbers with truncated exponential distribution

       Input
       -----
       a [double]: minimum of truncated exponential distribution
       shape: shape of numpy array to produce (None-> double scalar)

       Output
       ------
       zran [double]: scalar or array with independent random numbers

    
zran=truncated_normal([shape])

       Draw random numbers with truncated normal distribution

       Input
       -----
       a,b [double]: bounds of truncated normal distribution
       shape: shape of numpy array to produce (None-> double scalar)

       Output
       ------
       zran [double]: scalar or array with independent random numbers

    
sample = truncated_normal_vec(nsmp,A,b)

       Draw random vectors x with truncated normal distribution
       with N(0,I) refrerence Gaussian and contsraint Ax <= b

       Input
       -----
       nsmp [integer]: size of the sample to produce
       A [double]: matrix defining the constraint (ncstr,nvar)
       b [double]: vector defining the constraint (ncstr)

       Output
       ------
       sample [double]: sample of random vectors (nvar,nsmp)

    
field1d_init(spct_freq,spct_power)

       Initialize the sampling of 1D random fields

       Input
       -----
       spct_freq [double array]: list of frequencies
       spct_power [double array]: spectrum power at each frewuency

    
field2d_init(spct_freq,spct_power)

       Initialize the sampling of 2D random fields

       Input
       -----
       spct_freq [double array]: list of frequencies
       spct_power [double array]: spectrum power at each frewuency

    
field = field1d_sample(x,nharm)

       Sample 1D random fields with given spectrum

       Input
       -----
       x [double array]: grid of the output random field (1D)
       nharm [integer]: number of harmonics to sample from the spectrum and superpose

       Output
       ------
       field [double array]: random field with required spectrum (1D)

    
field = field2d_sample(x,y,nharm)

       Sample 2D random fields with given spectrum

       Input
       -----
       x [double array]: x coordinate of the output random field (2D)
       y [double array]: y coordinate of the output random field (2D)
       nharm [integer]: number of harmonics to sample from the spectrum and superpose

       Output
       ------
       field [double array]: random field with required spectrum (2D)

    
field = field2s_sample(lon,lat,pow_spectrum,lmin,lmax)

       Sample 2D random fields with given spectrum on the sphere

       Input
       -----
       lon [double array]: longitudes of the output random field (1D)
       lat [double array]: latitudes of the output random field (1D)
       pow_spectrum [callback routine]: power spectrum in the basis of the spherical harmonics
       lmin [integer] : minimum degree of the spherical harmonics
       lmax [integer] : maximum degree of the spherical harmonics

       Output
       ------
       field [double array]: random field with required spectrum (1D)

    

pyensdam.scores: Ensemble probablistic scores
=============================================

Available functions:
 -  scores.rank_histogram: Compute rank histogram
 -  scores.crps : Compute CRPS score (total, reliability, resolution)
 -  scores.rcrv : Compute RCRV score (bias, spread)
 -  scores.optimality : Compute OPTIMALITY score
 -  scores.entropy : Compute ENTROPY score (with option to compute entropy components)

Module parameters:
 -  scores.crps_missing_value: Missing value for CRPS score
 -  scores.rcrv_missing_value : missing value for RCRV score
 -  scores.rcrv_with_anamorphosis : apply anamorphosis rather than center-reduction in RCRV score
 -  scores.rcrv_number_of_quantiles : number of quantiles used in the anamorphosis transformation
 -  scores.optimality_missing_value : missing value for OPTIMALITY score
 -  scores.entropy_base : basis for the logarithm in entropy computations

Notes:
 - CRPS, RCRV and OPTIMALITY scores have the option to partition the input data
   and compute the score separately for each element of the partition.


rank_histogram,[ranks] = rank_histogram(ens,verif,[histogram_only=True])

       Compute ranks and rank histogram

       Inputs
       ------
       ens [rank-2 double array] : ensemble simulation (nens,nvar)
       verif [rank-1 double array] : verification data (nvar)
       histogram_only : compute only histogram (true or false)

       Returns
       -------
       rank_histogram : histogram of ranks
       ...and optionnally:
       ranks[rank-1 intc array] : ranks of observations in the input ensemble
    
crps, reliability, resolution = crps(ens,verif,[partition])

       Compute CRPS score (with option to partition the data)

       Inputs
       ------
       ens [rank-2 double array] : ensemble simulation (nens,nvar)
       verif [rank-1 double array] : verification data (nvar)
       partition [rank-1 intc array] : partition (nvar)

       Returns
       -------
       [global, or with as many component as elements in the partition]
       crps : total CRPS score
       reliability : reliability component of CRPS
       resolution : resolution component of CRPS
    
bias,spread = crps(ens,verif,[partition])

       Compute RCRV score (with option to partition the data)

       Inputs
       ------
       ens [rank-2 double array] : ensemble simulation (nens,nvar)
       verif [rank-1 double array] : verification data (nvar)
       partition [rank-1 intc array] : partition (nvar)

       Returns
       -------
       [global, or with as many component as elements in the partition]
       bias : bias component of RCRV score
       spread : spread component of RCRV score
    
optimality = optimality(ens,obs,[partition])

       Compute OPTIMALITY score (with option to partition the data)

       Observation error type ans spread are specified using attributes:
        obstype (normal, lognormal, gamma or beta) and
        obssigma (scalar or vector), see obserror module for definition

       Inputs
       ------
       ens [rank-2 double array] : ensemble simulation (nens,nvar)
       obs [rank-1 double array] : observations (nvar)
       obs_std [scalar or rank-1 double array] : observation error std (nvar)
       partition [rank-1 intc array] : partition (nvar)

       Returns
       -------
       [global, or with as many component as elements in the partition]
       optimality : optimality score
    
score = entropy(ens,pref,events_outcome,[score_only=True]))

       Compute ENTROPY score (with option to compute entropy components)

       Inputs
       ------
       ens [rank-2 double array] : ensemble simulation (nens,nvar)
       pref [rank-2 double array] : reference probability distribution (noutcomes,nevents)
       events_outcome [rank-1 intc array] : callback function

       Returns
       -------
       score [rank-1 double array] : entropy score (nevents)
       ...and optionnally:
       relative entropy [rank-1 double array (nevents)]
       cross entropy [rank-1 double array (nevents)]
       entropy [rank-1 double array (nevents)]

       Call-back function:
       -------------------
          def events_outcome(member): return outcome
          Required arguments:
            member [rank-1 double array] : ensemble member (nvar)
          Return:
            outcome [rank-1 intc array] : outcome for each event (nevents)
    

pyensdam.statistics: ensemble statistics
========================================

Available functions:
 -  statistics.meanstd : compute ensemble mean and standard deviation
 -  statistics.correlation : compute ensemble correlation
 -  statistics.covariance : compute ensemble covariance
 -  statistics.representer : compute ensemble representer


 mean,[std] = meanstd(ens,[weight],std=True)

        Compute ensemble mean and standard deviation

        Inputs
        ------
        ens [rank-2 double array] : ensemble simulation (nens,nvar)
        weight [rank-1 double array] : weight of ensemble members (nens)
        std : return also standard deviation (default=True)

        Outputs
        -------
        mean [rank-1 double array] : ensemble mean (nvar)
        std [rank-1 double array] : ensemble standard deviation (nvar)

    
 correlation = correlation(ens,ensref,[weight])

        Compute ensemble correlation
        with reference scalar ensemble

        Inputs
        ------
        ens [rank-2 double array] : ensemble simulation (nens,nvar)
        ensref [rank-1 double array] : reference scalar ensemble (nens)
        weight [rank-1 double array] : weight of ensemble members (nens)

        Outputs
        -------
        correlation [rank-1 double array] : correlation (nvar)

    
 covariance = covariance(ens,ensref,[weight])

        Compute ensemble covariance
        with reference scalar ensemble

        Inputs
        ------
        ens [rank-2 double array] : ensemble simulation (nens,nvar)
        ensref [rank-1 double array] : reference scalar ensemble (nens)
        weight [rank-1 double array] : weight of ensemble members (nens)

        Outputs
        -------
        covariance [rank-1 double array] : covariance (nvar)

    
 representer = representer(ens,ensref,[weight])

        Compute ensemble representer
        with reference scalar ensemble

        Inputs
        ------
        ens [rank-2 double array] : ensemble simulation (nens,nvar)
        ensref [rank-1 double array] : reference scalar ensemble (nens)
        weight [rank-1 double array] : weight of ensemble members (nens)

        Outputs
        -------
        representer [rank-1 double array] : representer (nvar)

    

pyensdam.transpho: transformation in the basis of the spherical harmonics
=========================================================================

Available functions:
 - transpho.forward : forward transformation (compute spectrum)
 - transpho.backward : backward transformation
 - transpho.mesh_area : compute area of grid cells on the sphere

Module parameters:


spectrum = forward(field,lon,lat)

       Forward transformation (compute spectrum)

       Inputs
       ------
       field [rank-1 double array] : field to transform (nvar)
       lon [rank-1 double array] : longitude of grid points (nvar)
       lat [rank-1 double array] : latitude of grid points (nvar)

       Returns
       -------
       spectrum [rank-2 double array] : spectrum in the basis of the spherical harmonics

    
field = forward(spectrum,lon,lat,[l0,l1])

       Backward transformation (from spectrum)

       Inputs
       ------
       spectrum [rank-2 double array] : spectrum in the basis of the spherical harmonics
       lon [rank-1 double array] : longitude of grid points (nvar)
       lat [rank-1 double array] : latitude of grid points (nvar)
       l0 [integer] : minimum degree to use in the backward transformation
       l1 [integer] : maximum degree to use in the backward transformation

       Returns
       -------
       field [rank-1 double array] : reconstructed field (nvar)

    
area = mesh_area(lon,lat)

       Compute the area of grid cells on the sphere

       Inputs
       ------
       lon [rank-2 double array] : longitude of grid points (nx,ny)
       lat [rank-2 double array] : latitude of grid points (nx,ny)

       Returns
       -------
       area [rank-2 double array] : area of grid cells (nx,ny)

    

pyensdam.update: ensemble update with addtional condtitions
===========================================================

Available functions:
 -  update.sample_mcmc : apply new condition (e.g. observations) on input ensemble,
                         using MCMC sampler with covariance localization

Module parameters:

 -  update.chain_index : current chain index
 -  update.zero_start : start from zero (T) or from restart ensemble (F)
 -  update.control_print : number of iteration between control prints
 -  update.convergence_check : number of iteration between convergence checks
 -  update.convergence_stop : stop at convergence (T) or perform full requested iterations (F)


upens, [upxens] = sample_mcmc(ens,multiplicity,nup,maxchain,my_jo,[my_test,xens])

 -     Apply new condition (e.g. observations) on input ensemble,
       using MCMC sampler with covariance localization

       Option xens is to include extra variables, not involved in the evaluation of the condition

       Inputs
       ------
       ens [rank-3 double array] : multiscale input ensemble (nscl,nens,nvar)
       multiplicity [rank-1 integer array] : multiplicity of each scale in the Schur products
       nup [integer] : number of members requested for the updated ensemble
       maxchain [integer] : number of iteration of the MCMC sampler
       my_jo [double function] : callback function with condition (cost function, taking member state as input)
       my_test [integer function] : callback function to test convergence and diagnose current properties of the updated ensemble
       xens [rank-3 double array] : extra variables in multiscale input ensemble (nscl,nens,nextra)

       Returns
       -------
       upens [rank-2 double array] : updated ensemble (nup,nvar)
       upxens [rank-2 double array] : updated ensemble for extra variables (nup,nextra), if xens is provided
    
