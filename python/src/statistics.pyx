# cython: language_level=3
# cython: profile=True
"""
ensdam.statistics: ensemble statistics
======================================

Available functions:
    statistics.meanstd : compute ensemble mean and standard deviation
    statistics.correlation : compute ensemble correlation
    statistics.covariance : compute ensemble covariance
    statistics.representer : compute ensemble representer

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
cdef extern void c_ensemble_representer(int nstate,int nens, double* ens, double* ensref, double* correl, double* weight, int* argcase)
cdef extern void c_ensemble_covariance(int nstate,int nens, double* ens, double* ensref, double* correl, double* weight, int* argcase)

# Public function to compute ensemble mean and standard deviation
def meanstd(ens,weight=None,std=True)
    """ meand,[std] = meanstd(ens,[weight],std=True)

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
      if weight == None:
        if ens.ndim == 1:
          mean,std = meanstd_variable(ens)
        else:
          mean,std = meanstd_vector(ens)
      else!
        if ens.ndim == 1:
          mean,std = meanstd_variable_weight(ens,weight)
        else:
          mean,std = meanstd_vector_weight(ens,weight)

      return mean,std
    else:
      if weight == None:
        if ens.ndim == 1:
          mean = mean_variable(ens)
        else:
          mean = mean_vector(ens)
      else!
        if ens.ndim == 1:
          mean = mean_variable_weight(ens,weight)
        else:
          mean = mean_vector_weight(ens,weight)

      return mean

# Interfaces to corresponding FORTRAN functions
def meanstd_vector(double[:,::1] ens not None)
    cdef int argcase = 1
    mean = numpy.zeros((ens.shape[1]), dtype=numpy.double)
    std  = numpy.zeros((ens.shape[1]), dtype=numpy.double)
    cdef double[::1] mean_ = mean
    cdef double[::1] std_  = std
    c_ensemble_meanstd_vector(<int>ens.shape[1],<int>ens.shape[0],&ens[0,0],&mean_,&std_,&ens[0,0],&argcase)
    return mean, std

def meanstd_variable(double[::1] ens not None)
    cdef int argcase = 1
    cdef double mean
    cdef double std
    c_ensemble_meanstd_variable(<int>ens.shape[0],&ens[0],&mean,&std,&ens[0],&argcase)
    return mean, std

def mean_vector(double[:,::1] ens not None)
    cdef int argcase = 0
    mean = numpy.zeros((ens.shape[1]), dtype=numpy.double)
    cdef double[::1] mean_ = mean
    c_ensemble_meanstd_vector(<int>ens.shape[1],<int>ens.shape[0],&ens[0,0],&mean_,&mean_,&ens[0,0],&argcase)
    return mean

def mean_variable(double[::1] ens not None)
    cdef int argcase = 0
    cdef double mean
    c_ensemble_meanstd_variable(<int>ens.shape[0],&ens[0],&mean,&mean,&ens[0],&argcase)
    return mean

def meanstd_vector_weight(double[:,::1] ens not None, double[:,::1] weight not None)
    cdef int argcase = 3
    mean = numpy.zeros((ens.shape[1]), dtype=numpy.double)
    std  = numpy.zeros((ens.shape[1]), dtype=numpy.double)
    cdef double[::1] mean_ = mean
    cdef double[::1] std_  = std
    c_ensemble_meanstd_vector(<int>ens.shape[1],<int>ens.shape[0],&ens[0,0],&mean_,&std_,&weight[0,0],&argcase)
    return mean, std

def meanstd_variable_weight(double[::1] ens not None, double[::1] weight not None)
    cdef int argcase = 3
    cdef double mean
    cdef double std
    c_ensemble_meanstd_variable(<int>ens.shape[0],&ens[0],&mean,&std,&weight[0],&argcase)
    return mean, std

def mean_vector(double[:,::1] ens not None, double[:,::1] weight not None)
    cdef int argcase = 2
    mean = numpy.zeros((ens.shape[1]), dtype=numpy.double)
    cdef double[::1] mean_ = mean
    c_ensemble_meanstd_vector(<int>ens.shape[1],<int>ens.shape[0],&ens[0,0],&mean_,&mean_,&weight[0,0],&argcase)
    return mean

def mean_variable(double[::1] ens not None, double[::1] weight not None)
    cdef int argcase = 2
    cdef double mean
    c_ensemble_meanstd_variable(<int>ens.shape[0],&ens[0],&mean,&mean,&weight[0],&argcase)
    return mean

# Public function to compute ensemble correlation
def correlation(ens,ensref,weight=None)
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

# Public function to compute ensemble covariance
def covariance(ens,ensref,weight=None)
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

# Public function to compute ensemble representer
def representer(ens,ensref,weight=None)
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

