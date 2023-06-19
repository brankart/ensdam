# cython: language_level=3
# cython: profile=True
"""
pyensdam.augment: ensemble augmentation
=======================================

Available functions:
 -  augment.sample_mcmc : resample input ensemble with MCMC sampler,
                          using covariance localization

"""

cimport cython
import tempfile
import atexit
import os

cimport numpy
import numpy

# Definition of external C-callable routines
cdef extern void c_sample_augmented_ensemble(int nvar,int nens,int nres,int naug,int* maxchain,double* augens,double* ens,int* multiplicity)

# Routines to exchange the module global variables
cdef extern void c_get_augment_chain_index(int* var) nogil
cdef extern void c_set_augment_chain_index(int* var) nogil

# Interface global variables of Fortran module into attributes of this module
cdef class __module_variable:
  # Current chain index
  property chain_index:
    def __get__(self):
      cdef int var
      c_get_augment_chain_index(&var)
      return var
    def __set__(self, int var):
      c_set_augment_chain_index(&var)

attr = __module_variable()

# Get default values of module attributes from Fortran module
chain_index = attr.chain_index

# Public function to augment input ensemble with MCMC sampler
def sample_mcmc(double[:,:,::1] ens not None,int[::1] multiplicity not None,naug not None, maxchain not None):
    """augens = sample_mcmc(ens,multiplicity,naug,maxchain)

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
    """
    # Update Fortran module public variables
    attr.chain_index = chain_index

    cdef int maxchain_ = maxchain
    augens = numpy.zeros((naug,<int>ens.shape[2]), dtype=numpy.double)
    cdef double[:,::1] augens_ = augens

    c_sample_augmented_ensemble(<int>ens.shape[2],<int>ens.shape[1],<int>ens.shape[0],<int>augens.shape[0],\
                                &maxchain_,&augens_[0,0],&ens[0,0,0],&multiplicity[0])

    return augens
