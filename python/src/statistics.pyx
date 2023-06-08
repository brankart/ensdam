# cython: language_level=3
# cython: profile=True
"""
pyensdam.statistics: ensemble statistics
========================================

Available functions:
 -  statistics.meanstd : compute ensemble mean and standard deviation
 -  statistics.correlation : compute ensemble correlation
 -  statistics.covariance : compute ensemble covariance
 -  statistics.representer : compute ensemble representer

"""

cimport cython
import tempfile
import atexit
import os

cimport numpy
import numpy

# Definition of external C-callable routines
cdef extern void c_ensemble_meanstd_vector(int nstate, int nens, double* ens, double* mean, double* std, double* weight, int* argcase)
cdef extern void c_ensemble_meanstd_variable(int nens, double* ens, double* mean, double* std, double* weight, int* argcase)
cdef extern void c_ensemble_correlation(int nstate,int nens, double* ens, double* ensref, double* correl, double* weight, int* argcase)
cdef extern void c_ensemble_representer(int nstate,int nens, double* ens, double* ensref, double* representer, double* weight, int* argcase)
cdef extern void c_ensemble_covariance(int nstate,int nens, double* ens, double* ensref, double* covariance, double* weight, int* argcase)

# Public function to compute ensemble mean and standard deviation
def meanstd(ens,weight=False,std=True):
    """ mean,[std] = meanstd(ens,[weight],std=True)

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

    """
    if std == True:
      if ens.ndim == 1:
        if not numpy.isscalar(weight):
          mean,std = meanstd_variable_weight(ens,weight)
        else:
          mean,std = meanstd_variable(ens)
      else:
        if not numpy.isscalar(weight):
          mean,std = meanstd_vector_weight(ens,weight)
        else:
          mean,std = meanstd_vector(ens)

      return mean,std
    else:
      if ens.ndim == 1:
        if not numpy.isscalar(weight):
          mean = mean_variable_weight(ens,weight)
        else:
          mean = mean_variable(ens)
      else:
        if not numpy.isscalar(weight):
          mean = mean_vector_weight(ens,weight)
        else:
          mean = mean_vector(ens)

      return mean

# Interfaces to corresponding FORTRAN functions
def meanstd_vector(double[:,::1] ens not None):
    cdef int argcase = 1
    mean = numpy.zeros((ens.shape[1]), dtype=numpy.double)
    std  = numpy.zeros((ens.shape[1]), dtype=numpy.double)
    cdef double[::1] mean_ = mean
    cdef double[::1] std_  = std
    c_ensemble_meanstd_vector(<int>ens.shape[1],<int>ens.shape[0],&ens[0,0],&mean_[0],&std_[0],&ens[0,0],&argcase)
    return mean, std

def meanstd_variable(double[::1] ens not None):
    cdef int argcase = 1
    cdef double mean
    cdef double std
    c_ensemble_meanstd_variable(<int>ens.shape[0],&ens[0],&mean,&std,&ens[0],&argcase)
    return mean, std

def mean_vector(double[:,::1] ens not None):
    cdef int argcase = 0
    mean = numpy.zeros((ens.shape[1]), dtype=numpy.double)
    cdef double[::1] mean_ = mean
    c_ensemble_meanstd_vector(<int>ens.shape[1],<int>ens.shape[0],&ens[0,0],&mean_[0],&mean_[0],&ens[0,0],&argcase)
    return mean

def mean_variable(double[::1] ens not None):
    cdef int argcase = 0
    cdef double mean
    c_ensemble_meanstd_variable(<int>ens.shape[0],&ens[0],&mean,&mean,&ens[0],&argcase)
    return mean

def meanstd_vector_weight(double[:,::1] ens not None, double[:,::1] weight not None):
    cdef int argcase = 3
    mean = numpy.zeros((ens.shape[1]), dtype=numpy.double)
    std  = numpy.zeros((ens.shape[1]), dtype=numpy.double)
    cdef double[::1] mean_ = mean
    cdef double[::1] std_  = std
    c_ensemble_meanstd_vector(<int>ens.shape[1],<int>ens.shape[0],&ens[0,0],&mean_[0],&std_[0],&weight[0,0],&argcase)
    return mean, std

def meanstd_variable_weight(double[::1] ens not None, double[::1] weight not None):
    cdef int argcase = 3
    cdef double mean
    cdef double std
    c_ensemble_meanstd_variable(<int>ens.shape[0],&ens[0],&mean,&std,&weight[0],&argcase)
    return mean, std

def mean_vector_weight(double[:,::1] ens not None, double[:,::1] weight not None):
    cdef int argcase = 2
    mean = numpy.zeros((ens.shape[1]), dtype=numpy.double)
    cdef double[::1] mean_ = mean
    c_ensemble_meanstd_vector(<int>ens.shape[1],<int>ens.shape[0],&ens[0,0],&mean_[0],&mean_[0],&weight[0,0],&argcase)
    return mean

def mean_variable_weight(double[::1] ens not None, double[::1] weight not None):
    cdef int argcase = 2
    cdef double mean
    c_ensemble_meanstd_variable(<int>ens.shape[0],&ens[0],&mean,&mean,&weight[0],&argcase)
    return mean

# Public function to compute ensemble correlation
def correlation(ens,ensref,weight=None):
    """ correlation = correlation(ens,ensref,[weight])

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

    """
    if weight == None:
      correlation = correlation_noweight(ens,ensref)
    else:
      correlation = correlation_weight(ens,ensref,weight)

    return correlation

# Interfaces to corresponding FORTRAN functions
def correlation_noweight(double[:,::1] ens not None, double[::1] ensref not None):
    cdef int argcase = 0
    correlation = numpy.zeros((ens.shape[1]), dtype=numpy.double)
    cdef double[::1] correlation_ = correlation
    c_ensemble_correlation(<int>ens.shape[1],<int>ens.shape[0],&ens[0,0],&ensref[0],&correlation_[0],&ens[0,0],&argcase)
    return correlation

def correlation_weight(double[:,::1] ens not None, double[::1] ensref not None, double[:,::1] weight not None):
    cdef int argcase = 1
    correlation = numpy.zeros((ens.shape[1]), dtype=numpy.double)
    cdef double[::1] correlation_ = correlation
    c_ensemble_correlation(<int>ens.shape[1],<int>ens.shape[0],&ens[0,0],&ensref[0],&correlation_[0],&weight[0,0],&argcase)
    return correlation

# Public function to compute ensemble covariance
def covariance(ens,ensref,weight=None):
    """ covariance = covariance(ens,ensref,[weight])

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

    """
    if weight == None:
      covariance = covariance_noweight(ens,ensref)
    else:
      covariance = covariance_weight(ens,ensref,weight)

    return covariance

# Interfaces to corresponding FORTRAN functions
def covariance_noweight(double[:,::1] ens not None, double[::1] ensref not None):
    cdef int argcase = 0
    covariance = numpy.zeros((ens.shape[1]), dtype=numpy.double)
    cdef double[::1] covariance_ = covariance
    c_ensemble_covariance(<int>ens.shape[1],<int>ens.shape[0],&ens[0,0],&ensref[0],&covariance_[0],&ens[0,0],&argcase)
    return covariance

def covariance_weight(double[:,::1] ens not None, double[::1] ensref not None, double[:,::1] weight not None):
    cdef int argcase = 1
    covariance = numpy.zeros((ens.shape[1]), dtype=numpy.double)
    cdef double[::1] covariance_ = covariance
    c_ensemble_covariance(<int>ens.shape[1],<int>ens.shape[0],&ens[0,0],&ensref[0],&covariance_[0],&weight[0,0],&argcase)
    return covariance

# Public function to compute ensemble representer
def representer(ens,ensref,weight=None):
    """ representer = representer(ens,ensref,[weight])

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

    """
    if weight == None:
      representerrepresenter = representer_noweight(ens,ensref)
    else:
      representer = representer_weight(ens,ensref,weight)

    return representer

# Interfaces to corresponding FORTRAN functions
def representer_noweight(double[:,::1] ens not None, double[::1] ensref not None):
    cdef int argcase = 0
    representer = numpy.zeros((ens.shape[1]), dtype=numpy.double)
    cdef double[::1] representer_ = representer
    c_ensemble_representer(<int>ens.shape[1],<int>ens.shape[0],&ens[0,0],&ensref[0],&representer_[0],&ens[0,0],&argcase)
    return representer

def representer_weight(double[:,::1] ens not None, double[::1] ensref not None, double[:,::1] weight not None):
    cdef int argcase = 1
    representer = numpy.zeros((ens.shape[1]), dtype=numpy.double)
    cdef double[::1] representer_ = representer
    c_ensemble_representer(<int>ens.shape[1],<int>ens.shape[0],&ens[0,0],&ensref[0],&representer_[0],&weight[0,0],&argcase)
    return representer
