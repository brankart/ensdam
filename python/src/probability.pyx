# cython: language_level=3
# cython: profile=True
"""
ensdam.probability: probability distribution
============================================

Available functions:
    probability.pdf: compute the probability density function
    probability.logpdf: compute the logartihm of the probability density function
    probability.cdf: compute the cumulative distribution function
    probability.invcdf: compute the iinverse cumulative distribution function

Module parameters:
    probability.type: type of probability distribution (normal, gamma, beta)
    probability.k: shape parameter of the gamma distribution
    probability.theta: scale parameter of the gamma distribution
    probability.a: parameter alpha of the beta distribution
    probability.b: parameter beta of the beta distribution

"""

cimport cython
import tempfile
import atexit
import os

cimport numpy
import numpy

# Definition of external C-callable routines
cdef extern void c_cdf_gaussian(double* x,double* fun) nogil
cdef extern void c_pdf_gaussian(double* x,double* fun) nogil
cdef extern void c_logpdf_gaussian(double* x,double* fun) nogil
cdef extern void c_invcdf_gaussian(double* rank,double* fun) nogil
cdef extern void c_cdf_gamma(double* a,double* x,double* fun) nogil
cdef extern void c_pdf_gamma(double* a,double* x,double* fun) nogil
cdef extern void c_logpdf_gamma(double* k,double* theta,double* x,double* fun) nogil
cdef extern void c_invcdf_gamma(double* a,double* rank,double* fun) nogil
cdef extern void c_cdf_beta(double* a,double* b,double* x,double* fun) nogil
cdef extern void c_pdf_beta(double* a,double* b,double* x,double* fun) nogil
cdef extern void c_logpdf_beta(double* a,double* b,double* x,double* fun) nogil
cdef extern void c_invcdf_beta(double* a,double* b,double* rank,double* fun) nogil

# Set default values of additional attributes
type='normal'
k=1.
theta=1.
a=2.
b=2.

# Public function to compute the probability density function
def pdf(x):
    """pdf=pdf(x)

       Compute the probability density function

       Input
       -----
       x [double]: value of the random variable

       Output
       ------
       pdf [double]: value of the pdf function

    """
    if numpy.isscalar(x):
      fun = pdf_any(x)
    else:
      fun = numpy.zeros_like(x)
      for i in range(numpy.prod(x.shape)):
        indices = numpy.unravel_index(i,x.shape)
        fun[indices] = pdf_any(x[indices])

    return fun

def pdf_any(x):
    cdef double fun_, k_, theta_, a_, b_
    cdef double x_ = x

    if type == 'normal' :
      c_pdf_gaussian(&x_,&fun_)
    elif type == 'gamma' :
      k_ = k ; theta_ = theta
      x_ = x_ / theta_
      c_pdf_gamma(&k_,&x_,&fun_)
      fun_ = fun_ / theta_
    elif type == 'beta' :
      a_ = a ; b_ = b
      c_pdf_beta(&a_,&b_,&x_,&fun_)
    else:
      raise ValueError("Invalid distribution type in probability")

    return fun_

# Public function to compute the log of probability density function
def logpdf(x):
    """logpdf=logpdf(x)

       Compute the log of probability density function

       Input
       -----
       x [double]: value of the random variable

       Output
       ------
       logpdf [double]: value of the log of the pdf function

    """
    if numpy.isscalar(x):
      fun = logpdf_any(x)
    else:
      fun = numpy.zeros_like(x)
      for i in range(numpy.prod(x.shape)):
        indices = numpy.unravel_index(i,x.shape)
        fun[indices] = logpdf_any(x[indices])

    return fun

def logpdf_any(x):
    cdef double fun_, k_, theta_, a_, b_
    cdef double x_ = x

    if type == 'normal' :
      c_logpdf_gaussian(&x_,&fun_)
    elif type == 'gamma' :
      k_ = k ; theta_ = theta
      c_logpdf_gamma(&k_,&theta_,&x_,&fun_)
    elif type == 'beta' :
      a_ = a ; b_ = b
      c_logpdf_beta(&a_,&b_,&x_,&fun_)
    else:
      raise ValueError("Invalid distribution type in probability")

    return fun_

# Public function to compute the cumulative distribution function
def cdf(x):
    """cdf=cdf(x)

       Compute the cumulative distribution function

       Input
       -----
       x [double]: value of the random variable

       Output
       ------
       cdf [double]: value of the cumulative distribution function

    """
    if numpy.isscalar(x):
      fun = cdf_any(x)
    else:
      fun = numpy.zeros_like(x)
      for i in range(numpy.prod(x.shape)):
        indices = numpy.unravel_index(i,x.shape)
        fun[indices] = cdf_any(x[indices])

    return fun

def cdf_any(x):
    cdef double fun_, k_, theta_, a_, b_
    cdef double x_ = x

    if type == 'normal' :
      c_cdf_gaussian(&x_,&fun_)
    elif type == 'gamma' :
      k_ = k ; theta_ = theta
      x_ = x_ / theta_
      c_cdf_gamma(&k_,&x_,&fun_)
    elif type == 'beta' :
      a_ = a ; b_ = b
      c_cdf_beta(&a_,&b_,&x_,&fun_)
    else:
      raise ValueError("Invalid distribution type in probability")

    return fun_

# Public function to compute the inverse cumulative distribution function
def invcdf(rank):
    """invcdf=invcdf(rank)

       Compute the inverse cumulative distribution function

       Input
       -----
       rank [double]: value of the random variable

       Output
       ------
       invcdf [double]: value of the inverse cumulative distribution function

    """
    if numpy.isscalar(rank):
      fun = invcdf_any(rank)
    else:
      fun = numpy.zeros_like(rank)
      for i in range(numpy.prod(rank.shape)):
        indices = numpy.unravel_index(i,rank.shape)
        fun[indices] = invcdf_any(rank[indices])

    return fun

def invcdf_any(rank):
    cdef double fun_, k_, theta_, a_, b_
    cdef double rank_ = rank

    if type == 'normal' :
      c_invcdf_gaussian(&rank_,&fun_)
    elif type == 'gamma' :
      k_ = k ; theta_ = theta
      c_invcdf_gamma(&k_,&rank_,&fun_)
      fun_ = fun_ * theta_
    elif type == 'beta' :
      a_ = a ; b_ = b
      c_invcdf_beta(&a_,&b_,&rank_,&fun_)
    else:
      raise ValueError("Invalid distribution type in probability")

    return fun_

