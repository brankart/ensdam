# cython: language_level=3
# cython: profile=True
"""
ensdam.obserror: Operations related to observation error
========================================================

Module attributes:
    obserror.obstype : Type of observation error
                       (normal, lognormal, gamma, beta)

Available functions:
    obserror.logpdf : Compute the logarithm of the observation error pdf
    obserror.cdf : Compute the observation error cdf
    obserror.sample : Sample the observation error probability distribution

"""

cimport cython
import tempfile
import atexit
import os

cimport numpy
import numpy

from libc.stdint cimport int8_t, int16_t, int32_t, int64_t

# Definition of external C-callable routines
cdef extern void c_get_obstype(char *var) nogil
cdef extern void c_set_obstype(int32_t *len1, char *var) nogil

cdef extern double c_obserror_logpdf_vector(int nobs, double* y, double* x, double* sigma) nogil
cdef extern double c_obserror_logpdf_vector_homogeneous(int nobs, double* y, double* x, double* sigma) nogil
cdef extern double c_obserror_logpdf_variable(double* y, double* x, double* sigma) nogil

cdef extern void c_obserror_cdf_vector(int nobs, double* y, double* x, double* sigma, double* cdf) nogil
cdef extern void c_obserror_cdf_vector_homogeneous(int nobs, double* y, double* x, double* sigma, double* cdf) nogil
cdef extern void c_obserror_cdf_variable(double* y, double* x, double* sigma, double* cdf) nogil

cdef extern void c_obserror_sample_vector(int nobs, double* x, double* sigma, double* sample) nogil
cdef extern void c_obserror_sample_vector_homogeneous(int nobs, double* x, double* sigma, double* sample, double* rank, int* uniform_rank, int* argcase) nogil
cdef extern void c_obserror_sample_variable(double* x, double* sigma, double* sample, double* rank, int* reuse_last_rank, int* argcase) nogil


# Utility to convert python string into C string
cdef pystring2cstring(str pystring):
  # add a null c char, and convert to bytes
  cdef bytes cstring = (pystring+'\0').encode('utf-8')
  return cstring

# Interface global variables of Fortran module into attributes of this module  
cdef class __module_variable:

  # Observation type (normal, lognormal, gamma or beta)
  property obstype:
    def __get__(self):
      cdef char var[80+1]
      c_get_obstype(var)
      return var.decode('utf-8')
    def __set__(self, str var):
      cdef int32_t len1 = len(var)
      cdef bytes var_b = pystring2cstring(var)
      cdef char *var_c = var_b
      c_set_obstype(&len1, var_c)

attr = __module_variable()

# Get default values of module attributes from Fortran module
obstype = attr.obstype

# Public function to compute the logarithm of the observation error pdf
def logpdf(y, x, sigma):
    """logpdf = logpdf(y,x,sigma)

       Compute the logarithm of the observation error pdf
       Supported distributions (obstype module attribute): normal, lognormal, gamma, beta

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
    """  
    cdef double result

    # Update Fortran module public variables
    attr.obstype = obstype

    # Apply appropriate function
    if numpy.isscalar(y):
        result = logpdf_variable(y, x, sigma)
    elif y.ndim == 1:
        if numpy.isscalar(sigma):
            result = logpdf_vector_homogeneous(y, x, sigma)
        elif sigma.ndim == 1:
            result = logpdf_vector(y, x, sigma)
        else:
            raise ValueError("Invalid array dimensions in logpdf")

    else:
        raise ValueError("Invalid array dimensions in logpdf")
    return result

# Interfaces to corresponding FORTRAN functions
def logpdf_variable(y not None, x not None, sigma not None):
    cdef double x_ = x
    cdef double y_ = y
    cdef double sigma_ = sigma
    cdef double result
    result = c_obserror_logpdf_variable(&y_, &x_, &sigma_)
    return result

def logpdf_vector(double[::1] y not None, double[::1] x not None, double[::1] sigma not None):
    cdef double result
    result = c_obserror_logpdf_vector(<int>y.shape[0], &y[0], &x[0], &sigma[0])
    return result

def logpdf_vector_homogeneous(double[::1] y not None, double[::1] x not None, sigma not None):
    cdef double result
    cdef double sigma_ = sigma
    result = c_obserror_logpdf_vector_homogeneous(<int>y.shape[0], &y[0], &x[0], &sigma_)
    return result


# Public function to compute the observation error cdf
def cdf(y, x, sigma):
    """cdf = cdf(y,x,sigma)

       Compute the observation error cdf
       Supported distributions (obstype module attribute): normal, lognormal, gamma, beta

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
    """
    # Update Fortran module public variables
    attr.obstype = obstype

    # Apply appropriate function
    if numpy.isscalar(y):
        result = cdf_variable(y, x, sigma)
    elif y.ndim == 1:
        if numpy.isscalar(sigma):
            result = cdf_vector_homogeneous(y, x, sigma)
        elif sigma.ndim == 1:
            result = cdf_vector(y, x, sigma)
        else:
            raise ValueError("Invalid array dimensions in cdf")

    else:
        raise ValueError("Invalid array dimensions in cdf")
    return result

# Interfaces to corresponding FORTRAN functions
def cdf_variable(y not None, x not None, sigma not None):
    cdef double x_ = x
    cdef double y_ = y
    cdef double sigma_ = sigma
    cdef double result
    c_obserror_cdf_variable(&y_, &x_, &sigma_, &result)
    return result

def cdf_vector(double[::1] y not None, double[::1] x not None, double[::1] sigma not None):
    result = numpy.zeros((<int>y.shape[0]), dtype=numpy.double)
    cdef double[::1] result_ = result
    c_obserror_cdf_vector(<int>y.shape[0], &y[0], &x[0], &sigma[0], &result_[0])
    return result

def cdf_vector_homogeneous(double[::1] y not None, double[::1] x not None, sigma not None):
    cdef double sigma_ = sigma
    result = numpy.zeros((<int>y.shape[0]), dtype=numpy.double)
    cdef double[::1] result_ = result
    c_obserror_cdf_vector_homogeneous(<int>y.shape[0], &y[0], &x[0], &sigma_, &result_[0])
    return result


# Public function to sample the observation error probability distribution
def sample(x, sigma, reuse_last_rank=False):
    """sample = sample(x,sigma)

       Sample the observation error probability distribution
       Supported distributions (obstype module attribute): normal, lognormal, gamma, beta

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
       """
    cdef int reuse_last_rank_ = 0
  
    if reuse_last_rank:
        reuse_last_rank_ = 1
    
    # Update Fortran module public variables
    attr.obstype = obstype

    # Apply appropriate function
    if numpy.isscalar(x):
        result = sample_variable(x, sigma, reuse_last_rank_)
    elif x.ndim == 1:
        if numpy.isscalar(sigma):
            result = sample_vector_homogeneous(x, sigma)
        elif sigma.ndim == 1:
            result = sample_vector(x, sigma)
        else:
            raise ValueError("Invalid array dimensions in sample")

    else:
        raise ValueError("Invalid array dimensions in sample")
    return result

# Interfaces to corresponding FORTRAN functions
def sample_variable(x not None, sigma not None, reuse_last_rank):
    cdef double x_ = x
    cdef double sigma_ = sigma
    cdef double result
    cdef double rank
    cdef int reuse_last_rank_ = reuse_last_rank
    cdef int argcase = 0
    if reuse_last_rank==1:
        argcase = 2

    c_obserror_sample_variable(&x_, &sigma_, &result, &rank, &reuse_last_rank_, &argcase)

    return result

def sample_vector(double[::1] x not None, double[::1] sigma not None):
    result = numpy.zeros((<int>x.shape[0]), dtype=numpy.double)
    cdef double[::1] result_ = result
    c_obserror_sample_vector(<int>x.shape[0], &x[0], &sigma[0], &result_[0])
    return result

def sample_vector_homogeneous(double[::1] x not None, sigma not None):
    cdef double sigma_ = sigma
    cdef double rank
    cdef int uniform_rank
    cdef int argcase = 0
    result = numpy.zeros((<int>x.shape[0]), dtype=numpy.double)
    cdef double[::1] result_ = result
    c_obserror_sample_vector_homogeneous(<int>x.shape[0], &x[0], &sigma_, &result_[0], &rank, &uniform_rank, &argcase)
    return result

