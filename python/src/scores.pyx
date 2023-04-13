# cython: language_level=3
# cython: profile=True
"""
ensdam.scores: Ensemble probablistic scores
===========================================

Available functions:
    scores.rank_histogram: Compute rank histogram
    scores.crps : Compute CRPS score (total, reliability, resolution)
    scores.rcrv : Compute RCRV score (bias, spread)
    scores.optimality : Compute OPTIMALITY score
    scores.entropy : Compute ENTROPY score (with option to compute entropy components)

Module parameters:
    scores.crps_missing_value: Missing value for CRPS score
    scores.rcrv_missing_value : missing value for RCRV score
    scores.rcrv_with_anamorphosis : apply anamorphosis rather than center-reduction in RCRV score
    scores.rcrv_number_of_quantiles : number of quantiles used in the anamorphosis transformation
    scores.optimality_missing_value : missing value for OPTIMALITY score
    scores.entropy_base : basis for the logarithm in entropy computations

Notes:
 - CRPS, RCRV and OPTIMALITY scores have the option to partition the input data
   and compute the score separately for each element of the partition.

"""

cimport cython
import tempfile
import atexit
import os

cimport numpy
import numpy
import pyensdam.obserror as obserror

# Definition of external C-callable routines
cdef extern void c_compute_ranks(int nstate, int nens, double* ens, double* verif, int* ranks) nogil
cdef extern void c_compute_ranks_histogram(int nstate, int nens, double* ens, double* verif, int* ranks, int* rank_histogram) nogil

cdef extern void c_crps_score_global(int nstate, int nens, double* crps, double* reliability, double* resolution, double* ens, double* verif) nogil
cdef extern void c_crps_score_partition(int nstate, int nens, int npartition, double* crps, double* reliability, double* resolution, double* ens, double* verif, int* partition) nogil

cdef extern void c_rcrv_score_global(int nstate, int nens, double* ens_bias, double* ens_spread, double* ens, double* verif) nogil
cdef extern void c_rcrv_score_partition(int nstate, int nens, int npartition, double* ens_bias, double* ens_spread, double* ens, double* verif, int* partition) nogil

cdef extern void c_optimality_score_global(int nstate, int nens, double* ens_optimality, double* ens, double* obs) nogil
cdef extern void c_optimality_score_global_gaussian(int nstate, int nens, double* ens_optimality, double* ens, double* obs, double* std_obs) nogil
cdef extern void c_optimality_score_partition(int nstate, int nens, int npartition, double* ens_optimality, double* ens, double* obs, int* partition) nogil
cdef extern void c_optimality_score_partition_gaussian(int nstate, int nens, int npartition, double* ens_optimality, double* ens, double* obs, double* std_obs, int* partition) nogil
cdef extern void c_associate_cdf_callback(void* cdf_callback_) nogil

cdef extern void c_events_score(int nstate, int nens, int nevents, int noutcomes, double* score, double* ens, double* pref)
cdef extern void c_events_relative_entropy(int nstate, int nens, int nevents, int noutcomes, double* relative_entropy, double* ens, double* pref)
cdef extern void c_events_cross_entropy(int nstate, int nens, int nevents, int noutcomes, double* cross_entropy, double* entropy, double* ens, double* pref)
cdef extern void c_events_entropy(int nstate, int nens, int nevents, int noutcomes, double* entropy, double* ens)
cdef extern void c_events_probability(int nstate, int nens, int nevents, int noutcomes, double* pens, double* ens)
cdef extern void c_associate_events_callback(void* events_callback_) nogil

# Routines to exchange the module global variables
cdef extern void c_get_crps_missing_value(double* var) nogil
cdef extern void c_set_crps_missing_value(double* var) nogil
cdef extern void c_get_rcrv_missing_value(double* var) nogil
cdef extern void c_set_rcrv_missing_value(double* var) nogil
cdef extern void c_get_rcrv_with_anamorphosis(int* var) nogil
cdef extern void c_set_rcrv_with_anamorphosis(int* var) nogil
cdef extern void c_get_rcrv_number_of_quantiles(int* var) nogil
cdef extern void c_set_rcrv_number_of_quantiles(int* var) nogil
cdef extern void c_get_optimality_missing_value(double* var) nogil
cdef extern void c_set_optimality_missing_value(double* var) nogil
cdef extern void c_get_score_entropy_base(double* var) nogil
cdef extern void c_set_score_entropy_base(double* var) nogil

# Define callback routines needed in C-callable wrapper
cdef double cdf_callback(double* o, double* y, int* obs_idx) with gil:
    cdef double o_ = numpy.asarray(<double[:1:1]> o)
    cdef double y_ = numpy.asarray(<double[:1:1]> y)
    cdef int i_ = numpy.asarray(<int[:1:1]> obs_idx)
    cdef double output
    obserror.obstype=obstype
    if numpy.isscalar(obssigma):
        output=obserror.cdf(o_,y_,obssigma)
    else:
        output=obserror.cdf(o_,y_,obssigma[i_-1])
    return output

cdef void events_callback(int nstate, int nevents, double* member, int* outcome) with gil:
    member_ = numpy.asarray(<double[:nstate:1]> member)
    result = glob_events_outcome(member_)
    for i in range(nevents):
       outcome[i] = result[i]
    return

# Interface global variables of Fortran module into attributes of this module
cdef class __module_variable:

  # Missing value for CRPS score
  property crps_missing_value:
    def __get__(self):
      cdef double var
      c_get_crps_missing_value(&var)
      return var
    def __set__(self, double var):
      c_set_crps_missing_value(&var)

  # Missing value for RCRV score
  property rcrv_missing_value:
    def __get__(self):
      cdef double var
      c_get_rcrv_missing_value(&var)
      return var
    def __set__(self, double var):
      c_set_rcrv_missing_value(&var)

  # Use anamorphosis to compute RCRV score
  property rcrv_with_anamorphosis:
    def __get__(self):
      cdef int var
      c_get_rcrv_with_anamorphosis(&var)
      return var
    def __set__(self, int var):
      c_set_rcrv_with_anamorphosis(&var)

  # Number of quantiles to perform anamorphosis in the computation of the RCRV score
  property rcrv_number_of_quantiles:
    def __get__(self):
      cdef int var
      c_get_rcrv_number_of_quantiles(&var)
      return var
    def __set__(self, int var):
      c_set_rcrv_number_of_quantiles(&var)

  # Missing value for OPTIMALITY score
  property optimality_missing_value:
    def __get__(self):
      cdef double var
      c_get_optimality_missing_value(&var)
      return var
    def __set__(self, double var):
      c_set_optimality_missing_value(&var)

  # Basis of logarithm in entropy computations
  property entropy_base:
    def __get__(self):
      cdef double var
      c_get_score_entropy_base(&var)
      return var
    def __set__(self, double var):
      c_set_score_entropy_base(&var)

attr = __module_variable()

# Get default values of module attributes from Fortran module
crps_missing_value = attr.crps_missing_value
rcrv_missing_value = attr.rcrv_missing_value
rcrv_with_anamorphosis = attr.rcrv_with_anamorphosis
rcrv_number_of_quantiles = attr.rcrv_number_of_quantiles
optimality_missing_value = attr.optimality_missing_value
entropy_base = attr.entropy_base
# Get default values of module attributes from other pyensdam modules
obstype=obserror.obstype
# Set default values of additional attributes
obssigma=1.
glob_events_outcome=None

# Public function to compute ranks and RANK HISTOGRAMS
def rank_histogram(ens,verif,histogram_only=True):
    """rank_histogram,[ranks] = rank_histogram(ens,verif,[histogram_only=True])

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
    """
    # Apply appropriate function
    ranks, rank_histogram = compute_ranks(ens,verif,histogram=True)

    # Return required output
    if histogram_only:
        return rank_histogram
    else:
        return rank_histogram,ranks

# Interfaces to corresponding FORTRAN functions
def compute_ranks(double[:,::1] ens not None, double[::1] verif not None, histogram=False):
    ranks = numpy.zeros((<int>ens.shape[1]), dtype=numpy.intc)
    cdef int[::1] ranks_ = ranks
    rank_histogram = numpy.zeros((<int>ens.shape[0]+1), dtype=numpy.intc)
    cdef int[::1] rank_histogram_ = rank_histogram

    if histogram:
        c_compute_ranks_histogram(<int>ens.shape[1], <int>ens.shape[0], &ens[0,0], &verif[0], &ranks_[0], &rank_histogram_[0])
        return ranks, rank_histogram
    else:
        c_compute_ranks(<int>ens.shape[1], <int>ens.shape[0], &ens[0,0], &verif[0], &ranks_[0])
        return ranks

# Public function to compute the CRPS score
def crps(ens,verif,partition=False):
    """crps, reliability, resolution = crps(ens,verif,[partition])

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
    """
    # Update Fortran module public variables
    attr.crps_missing_value = crps_missing_value

    # Apply appropriate function
    if numpy.isscalar(partition):
        crps,reliability,resolution = crps_score_global(ens, verif)
    else:
        partition = partition + 1
        crps,reliability,resolution = crps_score_partition(ens, verif, partition)

    return crps, reliability, resolution

# Interfaces to corresponding FORTRAN functions
def crps_score_global(double[:,::1] ens not None, double[::1] verif not None):
    cdef double crps
    cdef double reliability
    cdef double resolution

    c_crps_score_global(<int>ens.shape[1], <int>ens.shape[0], &crps, &reliability, &resolution, &ens[0,0], &verif[0])
    return crps, reliability, resolution

def crps_score_partition(double[:,::1] ens not None, double[::1] verif not None, int[::1] partition not None):
    cdef int npartition=numpy.amax(partition)

    if numpy.amin(partition) != 1:
        raise ValueError("Invalid partition in scores")

    crps = numpy.zeros((npartition), dtype=numpy.double)
    reliability = numpy.zeros((npartition), dtype=numpy.double)
    resolution = numpy.zeros((npartition), dtype=numpy.double)
    cdef double[::1] crps_ = crps
    cdef double[::1] reliability_ = reliability
    cdef double[::1] resolution_ = resolution

    c_crps_score_partition(<int>ens.shape[1], <int>ens.shape[0], <int>npartition, &crps_[0], &reliability_[0], &resolution_[0], &ens[0,0], &verif[0], &partition[0])
    return crps, reliability, resolution

# Public function to compute the RCRV score
def rcrv(ens,verif,partition=False):
    """bias,spread = crps(ens,verif,[partition])

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
    """
    # Update Fortran module public variables
    attr.rcrv_missing_value = rcrv_missing_value
    attr.rcrv_with_anamorphosis = rcrv_with_anamorphosis
    attr.rcrv_number_of_quantiles = rcrv_number_of_quantiles

    # Apply appropriate function
    if numpy.isscalar(partition):
        bias,spread = rcrv_score_global(ens, verif)
    else:
        partition = partition + 1
        bias,spread = rcrv_score_partition(ens, verif, partition)

    return bias, spread

# Interfaces to corresponding FORTRAN functions
def rcrv_score_global(double[:,::1] ens not None, double[::1] verif not None):
    cdef double ens_bias
    cdef double ens_spread

    c_rcrv_score_global(<int>ens.shape[1], <int>ens.shape[0], &ens_bias, &ens_spread, &ens[0,0], &verif[0])
    return ens_bias, ens_spread

def rcrv_score_partition(double[:,::1] ens not None, double[::1] verif not None, int[::1] partition not None):
    cdef int npartition=numpy.amax(partition)

    if numpy.amin(partition) != 1:
        raise ValueError("Invalid partition in scores")

    ens_bias = numpy.zeros((npartition), dtype=numpy.double)
    ens_spread = numpy.zeros((npartition), dtype=numpy.double)
    cdef double[::1] ens_bias_ = ens_bias
    cdef double[::1] ens_spread_ = ens_spread

    c_rcrv_score_partition(<int>ens.shape[1], <int>ens.shape[0], <int>npartition, &ens_bias_[0], &ens_spread_[0], &ens[0,0], &verif[0], &partition[0])
    return ens_bias, ens_spread

# Public function to compute the OPTIMALITY score
def optimality(ens,obs,obs_std,partition=False):
    """optimality = optimality(ens,obs,[partition])

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
    """
    # Update Fortran module public variables
    attr.optimality_missing_value = optimality_missing_value
    # Update module global variables
    obssigma = obs_std
    # Associate callback routines in C-callable wrapper
    c_associate_cdf_callback(&cdf_callback)

    # Apply appropriate function
    if obstype == 'normal' or obstype == 'gaussian' :
       if numpy.isscalar(obssigma):
          obssigma_ = numpy.zeros((ens.shape[1]), dtype=numpy.double) + obssigma
       else:
          obssigma_ = obssigma

       if numpy.isscalar(partition):
          optimality = optimality_score_global_gaussian(ens, obs, obssigma_)
       else:
          partition = partition + 1
          optimality = optimality_score_partition_gaussian(ens, obs, obssigma_, partition)
    else:
       if numpy.isscalar(partition):
          optimality = optimality_score_global(ens, obs)
       else:
          partition = partition + 1
          optimality = optimality_score_partition(ens, obs, partition)

    return optimality

# Interfaces to corresponding FORTRAN functions
def optimality_score_global(double[:,::1] ens not None, double[::1] obs not None):
    cdef double ens_optimality

    c_optimality_score_global(<int>ens.shape[1], <int>ens.shape[0], &ens_optimality, &ens[0,0], &obs[0])
    return ens_optimality

def optimality_score_global_gaussian(double[:,::1] ens not None, double[::1] obs not None, double[::1] std_obs not None):
    cdef double ens_optimality

    c_optimality_score_global_gaussian(<int>ens.shape[1], <int>ens.shape[0], &ens_optimality, &ens[0,0], &obs[0], &std_obs[0])
    return ens_optimality

def optimality_score_partition(double[:,::1] ens not None, double[::1] obs not None, int[::1] partition not None):
    cdef int npartition=numpy.amax(partition)

    if numpy.amin(partition) != 1:
        raise ValueError("Invalid partition in scores")

    ens_optimality = numpy.zeros((npartition), dtype=numpy.double)
    cdef double[::1] ens_optimality_ = ens_optimality

    c_optimality_score_partition(<int>ens.shape[1], <int>ens.shape[0], <int>npartition, &ens_optimality_[0], &ens[0,0], &obs[0], &partition[0])
    return ens_optimality

def optimality_score_partition_gaussian(double[:,::1] ens not None, double[::1] obs not None, double[::1] std_obs not None, int[::1] partition not None):
    cdef int npartition=numpy.amax(partition)

    if numpy.amin(partition) != 1:
        raise ValueError("Invalid partition in scores")

    ens_optimality = numpy.zeros((npartition), dtype=numpy.double)
    cdef double[::1] ens_optimality_ = ens_optimality

    c_optimality_score_partition_gaussian(<int>ens.shape[1], <int>ens.shape[0], <int>npartition, &ens_optimality_[0], &ens[0,0], &obs[0], &std_obs[0], &partition[0])
    return ens_optimality

# Public function to compute the ENTROPY score
def entropy(ens,pref,events_outcome,score_only=True):
    """score = entropy(ens,pref,events_outcome,[score_only=True]))

       Compute ENTROPY score (with option to compute entropy components)

       Inputs
       ------
       ens [rank-2 double array] : ensemble simulation (nens,nvar)
       pref [rank-2 double array] : reference probability distribution (nevents,noutcomes)
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
    """
    # Update Fortran module public variables
    attr.entropy_base = entropy_base
    # Associate callback routines in C-callable wrapper
    c_associate_events_callback(&events_callback)
    # Associate callback function to global name
    global glob_events_outcome
    glob_events_outcome = events_outcome

    # Apply appropriate function
    score = events_score(ens,pref)

    if not score_only :
        relative_entropy = events_relative_entropy(ens,pref)
        cross_entropy, entropy = events_cross_entropy(ens,pref)
        return score, relative_entropy, cross_entropy, entropy
    else: 
        return score

# Interfaces to corresponding FORTRAN functions
def events_score(double[:,::1] ens not None, double[:,::1] pref not None):

    score = numpy.zeros((<int>pref.shape[0]), dtype=numpy.double)
    cdef double[::1] score_ = score

    c_events_score(<int>ens.shape[1], <int>ens.shape[0], <int>pref.shape[0], <int>pref.shape[1], &score_[0], &ens[0,0], &pref[0,0])
    return score

def events_relative_entropy(double[:,::1] ens not None, double[:,::1] pref not None):

    relative_entropy = numpy.zeros((<int>pref.shape[0]), dtype=numpy.double)
    cdef double[::1] relative_entropy_ = relative_entropy

    c_events_relative_entropy(<int>ens.shape[1], <int>ens.shape[0], <int>pref.shape[0], <int>pref.shape[1], &relative_entropy_[0], &ens[0,0], &pref[0,0])
    return relative_entropy

def events_cross_entropy(double[:,::1] ens not None, double[:,::1] pref not None):

    cross_entropy = numpy.zeros((<int>pref.shape[0]), dtype=numpy.double)
    entropy = numpy.zeros((<int>pref.shape[0]), dtype=numpy.double)
    cdef double[::1] cross_entropy_ = cross_entropy
    cdef double[::1] entropy_ = entropy

    c_events_cross_entropy(<int>ens.shape[1], <int>ens.shape[0], <int>pref.shape[0], <int>pref.shape[1], &cross_entropy_[0], &entropy_[0], &ens[0,0], &pref[0,0])
    return cross_entropy, entropy

def events_entropy(double[:,::1] ens not None, nevents not None, noutcomes not None):

    entropy = numpy.zeros((nevents), dtype=numpy.double)
    cdef double[::1] entropy_ = entropy

    c_events_entropy(<int>ens.shape[1], <int>ens.shape[0], <int>nevents, <int>noutcomes, &entropy_[0], &ens[0,0])
    return entropy

def events_probability(double[:,::1] ens not None, nevents not None, noutcomes not None):

    pens = numpy.zeros((nevents,noutcomes), dtype=numpy.double)
    cdef double[:,::1] pens_ = pens

    c_events_probability(<int>ens.shape[1], <int>ens.shape[0], <int>nevents, <int>noutcomes, &pens_[0,0], &ens[0,0])
    return pens
