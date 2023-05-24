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
    random.field2s_sample: Sample 2D random fields with given spectrum on the sphere

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
cdef extern void c_kiss(long* iran) nogil
cdef extern void c_kiss_seed(long ix,long iy,long iz,long iw) nogil
cdef extern void c_kiss_save() nogil
cdef extern void c_kiss_load() nogil
cdef extern void c_kiss_reset() nogil
cdef extern void c_kiss_check(int32_t *len_check_type, char *check_type) nogil
cdef extern void c_kiss_uniform(double* uran) nogil
cdef extern void c_kiss_gaussian(double* gran) nogil
cdef extern void c_kiss_gamma(double* gamr,double* k) nogil
cdef extern void c_kiss_beta(double* betar,double* a,double* b) nogil
cdef extern void c_kiss_sample(int* a,int n,int k) nogil
cdef extern void c_ran_te(double* teran,double* a) nogil
cdef extern void c_ran_tg(int nsmpl,double* tgsmpl,double* aa,double* bb) nogil
cdef extern void c_ranv_tg(int nvar,int ncstr,int nsmpl,double* tgvsmpl,double* matArm,double* vecbm) nogil
cdef extern void c_def_spect_init(int nfreq,int nspct1d,int nspct2d,int nspct2s) nogil
cdef extern void c_def_spect_power(int spct_type,int spct_idx,int nspct,double* spct_freq,double* spct_power) nogil
cdef extern void c_def_sample_size(int nsmp1d,int nsmp2d,int nsmp2s) nogil
cdef extern void c_sample_freq_1d(int spct_idx) nogil
cdef extern void c_sample_freq_2d(int spct_idx) nogil
cdef extern void c_gen_field_1d(int spct_idx,int nx,double* ranfield,double* x) nogil
cdef extern void c_gen_field_2d(int spct_idx,int nx,int ny,double* ranfield,double* x,double* y) nogil
cdef extern void c_gen_field_2s(int ngrid,double* ranfield,double* lon,double* lat,int lmin,int lmax) nogil
cdef extern void c_associate_pow_spectrum_callback(void* pow_spectrum_callback_) nogil

# Utility to convert python string into C string
cdef pystring2cstring(str pystring):
    # add a null c char, and convert to bytes
    cdef bytes cstring = (pystring+'\0').encode('utf-8')
    return cstring

# Define callback routines needed in C-callable wrapper
cdef double pow_spectrum_callback(int* l,int* m) with gil:
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

    for i in range(0,seed_idx):
       c_kiss(&seed1)
       c_kiss(&seed2)
       c_kiss(&seed3)
       c_kiss(&seed4)

    if seed_idx>0:
       c_kiss_seed(seed1,seed2,seed3,seed4)

# Public function to save seed in restart file
def seed_save():
    """seed_save()

       Save current state of random number generator in restart file (.kiss_restart)
    """
    c_kiss_save()

# Public function to load seed from restart file
def seed_load():
    """seed_load()

       Load seed from restart file (.kiss_restart)
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
def swap(int[::1] a,k=None):
    """swap(a,[k])

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
    
    c_kiss_sample(&a[0],<int>a.shape[0],k_)

# Public function to draw random numbers with uniform distribution
def uniform(shape=None):
    """zran=uniform([shape])

       Draw random numbers with uniform distribution

       Input
       -----
       shape: shape of numpy array to produce (None-> double scalar)

       Output
       ------
       zran [double]: scalar or array with independent random numbers

    """
    cdef double zran_

    if shape==None:
      c_kiss_uniform(&zran_)
      return zran_
    else:
      zran = numpy.zeros(shape,dtype=numpy.double)
      with numpy.nditer(zran,op_flags=['readwrite']) as it:
        for x in it:
          c_kiss_uniform(&zran_)
          x[...] = zran_
      return zran

# Public function to draw random numbers with normal distribution
def normal(shape=None):
    """zran=normal([shape])

       Draw random numbers with normal distribution

       Input
       -----
       shape: shape of numpy array to produce (None-> double scalar)

       Output
       ------
       zran [double]: scalar or array with independent random numbers

    """
    cdef double zran_

    if shape==None:
      c_kiss_gaussian(&zran_)
      return zran_
    else:
      zran = numpy.zeros(shape,dtype=numpy.double)
      with numpy.nditer(zran,op_flags=['readwrite']) as it:
        for x in it:
          c_kiss_gaussian(&zran_)
          x[...] = zran_
      return zran

# Public function to draw random numbers with exponential distribution
def exp(shape=None):
    """zran=exp([shape])

       Draw random numbers with exponential distribution

       Input
       -----
       shape: shape of numpy array to produce (None-> double scalar)

       Output
       ------
       zran [double]: scalar or array with independent random numbers

    """
    cdef double zran_
    cdef double a_ = 0.

    if shape==None:
      c_ran_te(&zran_,&a_)
      return zran_
    else:
      zran = numpy.zeros(shape,dtype=numpy.double)
      with numpy.nditer(zran,op_flags=['readwrite']) as it:
        for x in it:
          c_ran_te(&zran_,&a_)
          x[...] = zran_
      return zran

# Public function to draw random numbers with gamma distribution
def gamma(k,shape=None):
    """zran=gamma(k,[shape])

       Draw random numbers with gamma distribution

       Input
       -----
       k [double]: parameter of gamma distribution
       shape: shape of numpy array to produce (None-> double scalar)

       Output
       ------
       zran [double]: scalar or array with independent random numbers

    """
    cdef double zran_
    cdef double k_ = k

    if shape==None:
      c_kiss_gamma(&zran_,&k_)
      return zran_
    else:
      zran = numpy.zeros(shape,dtype=numpy.double)
      with numpy.nditer(zran,op_flags=['readwrite']) as it:
        for x in it:
          c_kiss_gamma(&zran_,&k_)
          x[...] = zran_
      return zran

# Public function to draw random numbers with beta distribution
def beta(a,b,shape=None):
    """zran=beta(a,b,[shape])

       Draw random numbers with beta distribution

       Input
       -----
       a [double]: parameter of beta distribution
       b [double]: parameter of beta distribution
       shape: shape of numpy array to produce (None-> double scalar)

       Output
       ------
       zran [double]: scalar or array with independent random numbers

    """
    cdef double zran_
    cdef double a_ = a
    cdef double b_ = b

    if shape==None:
      c_kiss_beta(&zran_,&a_,&b_)
      return zran_
    else:
      zran = numpy.zeros(shape,dtype=numpy.double)
      with numpy.nditer(zran,op_flags=['readwrite']) as it:
        for x in it:
          c_kiss_beta(&zran_,&a_,&b_)
          x[...] = zran_
      return zran

# Public function to draw random numbers with truncated exponential distribution
def truncated_exp(a,shape=None):
    """zran=truncated_exp(a,[shape])

       Draw random numbers with truncated exponential distribution

       Input
       -----
       a [double]: minimum of truncated exponential distribution
       shape: shape of numpy array to produce (None-> double scalar)

       Output
       ------
       zran [double]: scalar or array with independent random numbers

    """
    cdef double zran_
    cdef double a_ = a

    if shape==None:
      c_ran_te(&zran_,&a_)
      return zran_
    else:
      zran = numpy.zeros(shape,dtype=numpy.double)
      with numpy.nditer(zran,op_flags=['readwrite']) as it:
        for x in it:
          c_ran_te(&zran_,&a_)
          x[...] = zran_
      return zran

# Public function to draw random numbers with truncated normal distribution
def truncated_normal(a,b,shape=None):
    """zran=truncated_normal([shape])

       Draw random numbers with truncated normal distribution

       Input
       -----
       a,b [double]: bounds of truncated normal distribution
       shape: shape of numpy array to produce (None-> double scalar)

       Output
       ------
       zran [double]: scalar or array with independent random numbers

    """
    cdef double zran_
    cdef double a_ = a
    cdef double b_ = b
    cdef int n = 1

    if shape==None:
      c_ran_tg(n,&zran_,&a_,&b_)
      return zran_
    else:
      zran = numpy.zeros(shape,dtype=numpy.double)
      with numpy.nditer(zran,op_flags=['readwrite']) as it:
        for x in it:
          c_ran_tg(n,&zran_,&a_,&b_)
          x[...] = zran_
      return zran

# Public function to draw random vectors with truncated normal distribution
def truncated_normal_vec(nsmp,double[:,::1] A,double[::1] b):
    """sample = truncated_normal_vec(nsmp,A,b)

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

    """
    sample = numpy.zeros((A.shape[1],nsmp), dtype=numpy.double)
    cdef double[:,::1] sample_ = sample
    cdef int nsmp_ = nsmp

    c_ranv_tg(<int>A.shape[1], <int>A.shape[0], nsmp_, &sample_[0,0], &A[0,0], &b[0])

    return sample

# Public function to initialize the sampling of 1D random fields
def field1d_init(double[::1] spct_freq, double[::1] spct_power):
    """field1d_init(spct_freq,spct_power)

       Initialize the sampling of 1D random fields

       Input
       -----
       spct_freq [double array]: list of frequencies
       spct_power [double array]: spectrum power at each frewuency

    """
    cdef int nspct1d=1,nspct2d=0,nspct2s=0,spct_type=1,spct_idx=1

    c_def_spect_init(<int>spct_freq.shape[0],nspct1d,nspct2d,nspct2s)
    c_def_spect_power(spct_type,spct_idx,<int>spct_freq.shape[0],&spct_freq[0],&spct_power[0])

# Public function to initialize the sampling of 2D random fields
def field2d_init(double[::1] spct_freq, double[::1] spct_power):
    """field2d_init(spct_freq,spct_power)

       Initialize the sampling of 2D random fields

       Input
       -----
       spct_freq [double array]: list of frequencies
       spct_power [double array]: spectrum power at each frewuency

    """
    cdef int nspct1d=0,nspct2d=1,nspct2s=0,spct_type=2,spct_idx=1

    c_def_spect_init(<int>spct_freq.shape[0],nspct1d,nspct2d,nspct2s)
    c_def_spect_power(spct_type,spct_idx,<int>spct_freq.shape[0],&spct_freq[0],&spct_power[0])

# Public function to sample 1D random fields with given spectrum
def field1d_sample(double[::1] x,nharm):
    """field = field1d_sample(x,nharm)

       Sample 1D random fields with given spectrum

       Input
       -----
       x [double array]: grid of the output random field (1D)
       nharm [integer]: number of harmonics to sample from the spectrum and superpose

       Output
       ------
       field [double array]: random field with required spectrum (1D)

    """
    cdef int nsmp2d=0,nsmp2s=0,spct_idx=1
    field = numpy.zeros((<int>x.shape[0]), dtype=numpy.double)
    cdef double[::1] field_ = field

    c_def_sample_size(<int>nharm,nsmp2d,nsmp2s)
    c_sample_freq_1d(spct_idx)
    c_gen_field_1d(spct_idx,<int>x.shape[0],&field_[0],&x[0])

    return field

# Public function to sample 2D random fields with given spectrum
def field2d_sample(double[:,::1] x, double[:,::1] y,nharm):
    """field = field2d_sample(x,y,nharm)

       Sample 2D random fields with given spectrum

       Input
       -----
       x [double array]: x coordinate of the output random field (2D)
       y [double array]: y coordinate of the output random field (2D)
       nharm [integer]: number of harmonics to sample from the spectrum and superpose

       Output
       ------
       field [double array]: random field with required spectrum (2D)

    """
    cdef int nsmp1d=0,nsmp2s=0,spct_idx=1
    field = numpy.zeros((<int>x.shape[0],<int>x.shape[1]), dtype=numpy.double)
    cdef double[:,::1] field_ = field

    c_def_sample_size(nsmp1d,<int>nharm,nsmp2s)
    c_sample_freq_2d(spct_idx)
    c_gen_field_2d(spct_idx,<int>x.shape[0],<int>x.shape[1],&field_[0,0],&x[0,0],&y[0,0])

    return field

# Public function to sample 2D random fields with given spectrum on the sphere
def field2s_sample(double[::1] lon, double[::1] lat,pow_spectrum,lmin,lmax):
    """field = field2s_sample(lon,lat,pow_spectrum,lmin,lmax)

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

    """
    field = numpy.zeros((<int>lon.shape[0]), dtype=numpy.double)
    cdef double[::1] field_ = field

    # Associate callback routines in C-callable wrapper
    c_associate_pow_spectrum_callback(&pow_spectrum_callback)
    # Associate callback function to global name
    global glob_pow_spectrum
    glob_pow_spectrum = pow_spectrum

    c_gen_field_2s(<int>lon.shape[0],&field_[0],&lon[0],&lat[0],<int>lmin,<int>lmax)

    return field

