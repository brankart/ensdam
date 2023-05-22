# cython: language_level=3
# cython: profile=True
"""
ensdam.random: generate random numbers and fields
=================================================

Available functions:
    random.seed: Seed random number generator
    random.seed_save: Save seed in restart file
    random.seed_load: Load seed from restart file
    random.check: Check random number generator
    random.swap: Random array swapping
    random.uniform: Draw random numbers with uniform distribution
    random.normal: Draw random numbers with normal distribution
    random.exp: Draw random numbers with exponential distribution
    random.gamma: Draw random numbers with gamma distribution
    random.beta: Draw random numbers with beta distribution
    random.truncated_exp: Draw random numbers with truncated exponential distribution
    random.truncated_normal: Draw random numbers with truncated normal distribution
    random.truncated_normal_vec: Draw random vectors with truncated normal distribution
    random.field1d_init: Initialization for the sampling of 1D random fields
    random.field2d_init: Initialization for the sampling of 2D random fields
    random.field1d_sample: Sample 1D random fields with given spectrum
    random.field2d_sample: Sample 2D random fields with given spectrum

Module parameters:

"""

cimport cython
import tempfile
import atexit
import os

cimport numpy
import numpy

from libc.stdint cimport int8_t, int16_t, int32_t, int64_t

# Definition of external C-callable routines
cdef extern void c_kiss(long iran) nogil
cdef extern void c_kiss_seed(long ix,long iy,long iz,long iw) nogil
cdef extern void c_kiss_save() nogil
cdef extern void c_kiss_load() nogil
cdef extern void c_kiss_reset() nogil
cdef extern void c_kiss_check(int32_t *len_check_type, char *check_type) nogil
cdef extern void c_kiss_uniform(double* uran) nogil
cdef extern void c_kiss_gaussian(double* gran) nogil
cdef extern void c_kiss_gamma(double* gamr,double* k) nogil
cdef extern void c_kiss_beta(double* betar,double* a,double* b) nogil
cdef extern void c_kiss_sample(int* a,int n,int* k) nogil
cdef extern void c_ran_te(double* teran,double* a) nogil
cdef extern void c_ran_tg(int nsmpl,double* tgsmpl,double* aa,double* bb) nogil
cdef extern void c_ranv_tg(int nvar,int ncstr,int nsmpl,double* tgvsmpl,double* matArm,double* vecbm) nogil
cdef extern void c_def_spect_init(int nfreq,int nspct1d,int nspct2d,int nspct2s) nogil
cdef extern void c_def_spect_power(int spct_type,int spct_idx,int nspct,double* spct_freq,double* spct_power) nogil
cdef extern void c_def_sample_size(int nsmp1d,int nsmp2d,int nsmp2s,int nseed) nogil
cdef extern void c_sample_freq_1d(int spct_idx) nogil
cdef extern void c_sample_freq_2d(int spct_idx) nogil
cdef extern void c_gen_field_1d(int spct_idx,int nx,ranfield,x) nogil
cdef extern void c_gen_field_2d(int spct_idx,int nx,int ny,ranfield,x,y) nogil
cdef extern void c_gen_field_2s(int ngrid,ranfield,lon,lat,int lmin,int lmax) nogil

# Utility to convert python string into C string
cdef pystring2cstring(str pystring):
    # add a null c char, and convert to bytes
    cdef bytes cstring = (pystring+'\0').encode('utf-8')
    return cstring

# Define callback routines needed in C-callable wrapper
cdef double pow_spectrum_callback(int* l,int* m):
    cdef int l_ = numpy.asarray(<int[:1:1]> l) 
    cdef int m_ = numpy.asarray(<int[:1:1]> m) 
    cdef double output
    output = glob_pow_spectrum(l_,m_)
    return output

# Set default values of additional attributes
glob_pow_spectrum=None

# Public function to seed random number generator
def seed(seed_idx not None):
    """seed(seed_idx)

       Seed random number generator

       Inputs
       ------
       seed_idx[integer] : index of seed to use
    """
    cdef long seed1=0, seed2=0, seed3=0, seed4=0

    c_kiss_reset()

    for i in range(1,seed_idx):
       c_kiss(seed1)
       c_kiss(seed2)
       c_kiss(seed3)
       c_kiss(seed4)

    c_kiss_seed(seed1,seed2,seed3,seed4)

# Public function to save seed in restart file
def seed_save():
    """seed_save()

       Save current state of random number generator in restart file
    """
    c_kiss_save()

# Public function to load seed from restart file
def seed_load():
    """seed_load()

       Load seed from restart file
    """
    c_kiss_load()

# Public function to check random number generator
def check(check_type not None):
    """check(check_type)

       Check random number generator

       Inputs
       ------
       check_type[string] : type of check ('short' or 'long')

    """
    cdef int32_t len_check_type = len(check_type)
    cdef bytes b_check_type = pystring2cstring(check_type)
    cdef char *c_check_type = b_check_type

    c_kiss_check(&len_check_type, c_check_type)

# Public function to perform random swap of input array
def swap(int[::1] a not None,k):
    """swap()

       Random array swapping

       Inputs
       ------
       a [integer vector] : input array to swap
       k [integer]: number of elements to sample from array (default=size of a)

    """
    cdef int k_

    if k == None:
      k_ = a.shape[0]
    else:
      k_ = k
    
    c_kiss_sample(&a[0],<int>a.shape[0],&k_)

# Public function to draw random numbers with uniform distribution
def uniform(zran):
    """uniform(zran)

       Draw random numbers with uniform distribution

       Input/Output
       ------------
       zran[double] : scalar/array to fill with uniform random numbers

    """
    cdef double zran_

    if numpy.isscalar(zran):
      zran_ = zran 
      c_kiss_uniform(&zran_)
    else:
      for zran1 in zran:
        zran_ = zran1
        c_kiss_uniform(&zran_)

# Public function to draw random numbers with normal distribution
def normal():
    """normal()

       Draw random numbers with normal distribution
    """

# Public function to draw random numbers with exponential distribution
def exp():
    """exp()

       Draw random numbers with exponential distribution
    """

# Public function to draw random numbers with gamma distribution
def gamma():
    """gamma()

       Draw random numbers with gamma distribution
    """

# Public function to draw random numbers with beta distribution
def beta():
    """beta()

       Draw random numbers with beta distribution
    """

# Public function to draw random numbers with truncated exponential distribution
def truncated_exp():
    """truncated_exp()

       Draw random numbers with truncated exponential distribution
    """

# Public function to draw random numbers with truncated normal distribution
def truncated_normal():
    """truncated_normal()

       Draw random numbers with truncated normal distribution
    """

# Public function to draw random vectors with truncated normal distribution
def truncated_normal_vec():
    """truncated_normal_vec()

       Draw random numbers with truncated normal distribution
    """

# Public function to initialize the sampling of 1D random fields
def field1d_init():
    """field1d_init()

       Initialize the sampling of 1D random fields
    """

# Public function to initialize the sampling of 2D random fields
def field2d_init():
    """field2d_init()

       Initialize the sampling of 2D random fields
    """

# Public function to sample 1D random fields with given spectrum
def field1d_sample():
    """field1d_sample()

       Sample 1D random fields with given spectrum
    """

# Public function to sample 2D random fields with given spectrum
def field2d_sample():
    """field2d_sample()

       Sample 2D random fields with given spectrum
    """

