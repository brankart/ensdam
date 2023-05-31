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
def meanstd(ens,weight=None)
    """ meand,[std] = meanstd(ens,[weight])

        Compute ensemble mean and standard deviation

        Inputs
        ------
        ens [rank-2 double array] : ensemble simulation (nens,nvar)
        weight [rank-1 double array] : weight of ensemble members (nens)

        Outputs
        -------
        mean [rank-1 double array] : ensemble mean (nvar)
        std [rank-1 double array] : ensemble standard deviation (nvar)

    """

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

