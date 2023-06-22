# cython: language_level=3
# cython: profile=True
"""
pyensdam.update: ensemble update with addtional condtitions
===========================================================

Available functions:
 -  update.sample_mcmc : apply new condition (e.g. observations) on input ensemble,
                         using MCMC sampler with covariance localization

Module parameters:

 -  update.chain_index :
 -  update.zero_start : start from zero (T) or from restart ensemble (F)
 -  update.control_print : number of iteration between control prints
 -  update.convergence_check : number of iteration between convergence checks
 -  update.convergence_stop : stop at convergence (T) or perform full requested iterations (F)

"""

cimport cython
import tempfile
import atexit
import os

cimport numpy
import numpy

# Definition of external C-callable routines
cdef extern void c_mcmc_iteration(int nvar,int nens,int nres,int nup,int nextra,int* maxchain,double* upens,double* ens,int* multiplicity,double* upxens,double* xens,int* argcase)

cdef extern void c_associate_my_jo_callback(void* my_jo_in) 
cdef extern void c_associate_my_test_callback(void* my_test_in)

# Routines to exchange the module global variables
cdef extern void c_get_mcmc_index(int* var) nogil
cdef extern void c_set_mcmc_index(int* var) nogil
cdef extern void c_get_mcmc_zero_start(int* var) nogil
cdef extern void c_set_mcmc_zero_start(int* var) nogil
cdef extern void c_get_mcmc_control_print(int* var) nogil
cdef extern void c_set_mcmc_control_print(int* var) nogil
cdef extern void c_get_mcmc_convergence_check(int* var) nogil
cdef extern void c_set_mcmc_convergence_check(int* var) nogil
cdef extern void c_get_mcmc_convergence_stop(int* var) nogil
cdef extern void c_set_mcmc_convergence_stop(int* var) nogil
cdef extern void c_get_mcmc_proposal(int* var) nogil
cdef extern void c_set_mcmc_proposal(int* var) nogil
cdef extern void c_get_mcmc_proposal_std(double* var) nogil
cdef extern void c_set_mcmc_proposal_std(double* var) nogil
cdef extern void c_get_mcmc_schedule(double* var) nogil
cdef extern void c_set_mcmc_schedule(double* var) nogil

# Define callback routines needed in C-callable wrapper
cdef double my_jo_callback(int nvar, double* v) with gil:
  cdef double my_jo
  v_ = numpy.asarray(<double[:nvar:1]>v)
  my_jo = glob_my_jo(v_)
  return my_jo

cdef int my_test_callback(int nvar, int nens, int nextra, double* upens, double* upxens) with gil:
  cdef int my_test

  if nextra == 0:
    upens_ = numpy.asarray(<double[:nens:1,:nvar:1]>upens)
    my_test = glob_my_test(upens_)
  else:
    upens_ = numpy.asarray(<double[:nens:1,:nvar:1]>upens)
    upxens_ = numpy.asarray(<double[:nens:1,:nextra:1]>upxens)
    my_test = glob_my_test_extra(upens_,upxens_)

  return my_test

# Interface global variables of Fortran module into attributes of this module
cdef class __module_variable:
  # Current chain index
  property chain_index:
    def __get__(self):
      cdef int var
      c_get_mcmc_index(&var)
      return var
    def __set__(self, int var):
      c_set_mcmc_index(&var)
  # Start from zero or from a restart ensemble
  property zero_start:
    def __get__(self):
      cdef int var
      c_get_mcmc_zero_start(&var)
      return var
    def __set__(self, int var):
      c_set_mcmc_zero_start(&var)
  # Number of iteration between control prints
  property control_print:
    def __get__(self):
      cdef int var
      c_get_mcmc_control_print(&var)
      return var
    def __set__(self, int var):
      c_set_mcmc_control_print(&var)
  # Number of iteration between convergence checks
  property convergence_check:
    def __get__(self):
      cdef int var
      c_get_mcmc_convergence_check(&var)
      return var
    def __set__(self, int var):
      c_set_mcmc_convergence_check(&var)
  # Stop at convergence or not
  property convergence_stop:
    def __get__(self):
      cdef int var
      c_get_mcmc_convergence_stop(&var)
      return var
    def __set__(self, int var):
      c_set_mcmc_convergence_stop(&var)
  # Input ensemble is a proposal distribution
  property mcmc_proposal:
    def __get__(self):
      cdef int var
      c_get_mcmc_proposal(&var)
      return var
    def __set__(self, int var):
      c_set_mcmc_proposal(&var)
  # Proposal distribution std
  property mcmc_proposal_std:
    def __get__(self):
      cdef double var
      c_get_mcmc_proposal_std(&var)
      return var
    def __set__(self, double var):
      c_set_mcmc_proposal_std(&var)
  # MCMC schedule
  property mcmc_schedule:
    def __get__(self):
      cdef double var
      c_get_mcmc_schedule(&var)
      return var
    def __set__(self, double var):
      c_set_mcmc_schedule(&var)

attr = __module_variable()

# Get default values of module attributes from Fortran module
chain_index = attr.chain_index
zero_start = attr.zero_start
control_print = attr.control_print
convergence_check = attr.convergence_check
convergence_stop = attr.convergence_stop
mcmc_proposal = attr.mcmc_proposal
mcmc_proposal_std = attr.mcmc_proposal_std
mcmc_schedule = attr.mcmc_schedule

# Set default values of additional attributes
glob_my_jo=None
glob_my_test=None
glob_my_test_extra=None

# Public function to update input ensemble with MCMC sampler
def sample_mcmc(ens not None, multiplicity not None, nup not None, maxchain not None, my_jo not None, my_test=None, xens=False):
    """upens, [upxens] = sample_mcmc(ens,multiplicity,nup,maxchain,my_jo,[my_test,xens])

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
    """
    # Update Fortran module public variables
    attr.chain_index = chain_index
    attr.zero_start = zero_start
    attr.control_print = control_print
    attr.convergence_check = convergence_check
    attr.convergence_stop = convergence_stop
    attr.mcmc_proposal = mcmc_proposal
    attr.mcmc_proposal_std = mcmc_proposal_std
    attr.mcmc_schedule = mcmc_schedule

    # Associate callback routine in C-callable wrapper
    c_associate_my_jo_callback(&my_jo_callback)
    if my_test != None:
      c_associate_my_test_callback(&my_test_callback)

    # Associate callback functions to global name
    global glob_my_jo, glob_my_test
    glob_my_jo = my_jo
    glob_my_test = my_test

    # Compute updated ensemble
    if numpy.isscalar(xens):
      # No extra variables
      glob_my_test = my_test
      upens = sample_mcmc_noextra(ens,multiplicity,nup,maxchain)
      return upens
    else:
      # With extra variables
      glob_my_test_extra = my_test
      upens, upxens = sample_mcmc_extra(ens,multiplicity,nup,maxchain,xens)
      return upens, upxens

# Ensemble update without extra variables
def sample_mcmc_noextra(double[:,:,::1] ens not None,int[::1] multiplicity not None,nup not None, maxchain not None):
    cdef int argcase = 0
    cdef int nextra = 1
    cdef int maxchain_ = maxchain
    upens = numpy.zeros((nup,<int>ens.shape[2]), dtype=numpy.double)
    cdef double[:,::1] upens_ = upens

    c_mcmc_iteration(<int>ens.shape[2],<int>ens.shape[1],<int>ens.shape[0],<int>upens.shape[0],nextra,&maxchain_,\
                     &upens_[0,0],&ens[0,0,0],&multiplicity[0],&upens_[0,0],&ens[0,0,0],&argcase)

    return upens

# Ensemble update with extra variables
def sample_mcmc_extra(double[:,:,::1] ens not None,int[::1] multiplicity not None,nup not None,maxchain not None,double[:,:,::1] xens):
    cdef int argcase = 1
    cdef int nextra = xens.shape[2]
    cdef int maxchain_ = maxchain
    upens  = numpy.zeros((nup,<int>ens.shape[2]),  dtype=numpy.double)
    cdef double[:,::1] upens_ = upens
    upxens = numpy.zeros((nup,<int>xens.shape[2]), dtype=numpy.double)
    cdef double[:,::1] upxens_ = upxens

    c_mcmc_iteration(<int>ens.shape[2],<int>ens.shape[1],<int>ens.shape[0],<int>upens.shape[0],nextra,&maxchain_,\
                     &upens_[0,0],&ens[0,0,0],&multiplicity[0],&upxens_[0,0],&xens[0,0,0],&argcase)

    return upens, upxens

