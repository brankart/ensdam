# cython: language_level=3
# cython: profile=True
"""
ensdam.anamorphosis: ensemble anamorphosis
==========================================

Available functions:
    anamorphosis.quantiles : compute ensemble quantiles
    anamorphosis.forward : forward anamorphosis transformation
    anamorphosis.backward : backward anamorphosis transformation
    anamorphosis.forward_obs : forward anamorphosis transformation of observations
    anamorphosis.forward_obs_sym : forward anamorphosis transformation of observations (symmetric)

Module parameters:
    anamorphosis.target : target probability distribution to use (default=normal, uniform, gamma, beta)
    anamorphosis.obstype : probability distribution of observations (default=normal, gamma, beta)

"""

cimport cython
import tempfile
import atexit
import os

cimport numpy
import numpy

import pyensdam.obserror as obserror

# Definition of external C-callable routines
cdef extern void c_ens_quantiles_vector(int nvar,int nens,int nqua,double* qua,double* ens,double* quadef,double* enswei,double* ensweiloc,int* argcase)
cdef extern void c_ens_quantiles_variable(int nens,int nqua,double* qua,double* ens,double* quadef,double* enswei,int* argcase)
cdef extern void c_ana_forward_ensemble(int nvar,int nens,int nqua,double* ens,double* qua,double* quaref,double* rank,int* argcase)
cdef extern void c_ana_forward_vector(int nvar,int nqua,double* vct,double* qua,double* quaref,double* rank,int* argcase)
cdef extern void c_ana_forward_variable(int nqua,double* var,double* qua,double* quaref,double* rank,int* argcase)
cdef extern void c_ana_backward_ensemble(int nvar,int nens,int nqua,double* ens,double* qua,double* quaref)
cdef extern void c_ana_backward_vector(int nvar,int nqua,double* vct,double* qua,double* quaref)
cdef extern void c_ana_backward_variable(int nqua,double* var,double* qua,double* quaref)
cdef extern void c_ana_obs(int nobs,int nens,int nqua,int nsmp,double* anaobs,double* obsens,double* obs,double* obs_std,double* quadef,double* quaref)
cdef extern void c_ana_obs_sym(int nobs,int nqua,int nsmp,double* anaobs,double* obs,double* obs_std,double* obsqua,double* quaref)

# Get default values of module attributes from other pyensdam modules
obstype=obserror.obstype
# Set default values of additional attributes
target='normal'   # type of target distribution
quaref=None       # quantiles of the target distribution

# Public function to compute ensemble quantiles
def quantiles(ens,quadef,weight=False,local_weight=False):
    """qua = quantiles(ens,quadef,[weight,local_weight])

       Compute ensemble quantiles

       Inputs
       ------
       ens [rank-1 or rank-2 double array] : ensemble simulation (nens) or (nens,nvar)
       quadef [rank-1 double array] : definition of the quantiles (nqua)
       weight [rank-1 double array] : weight of ensemble members (nens)
       local_weight [rank-2 double array] : local weight of ensemble members (nens,nvar)

       Returns
       -------
       qua [rank-1 or rank-2 double array] : ensemble quantiles (nqua) or (nqua,nvar)

    """
    if ens.ndim == 1:
      if not numpy.isscalar(weight):
        qua = quantiles_variable_weight(ens,quadef,weight)
      else:
        qua = quantiles_variable(ens,quadef)
    elif ens.ndim == 2:
      if not numpy.isscalar(local_weight):
        qua = quantiles_vector_locweight(ens,quadef,local_weight)
      elif not numpy.isscalar(weight):
        qua = quantiles_vector_weight(ens,quadef,weight)
      else:
        qua = quantiles_vector(ens,quadef)

# Interfaces to corresponding FORTRAN functions
def quantiles_variable(double[::1] ens,double[::1] quadef):
    cdef int argcase = 0
    qua  = numpy.zeros((quadef.shape[0]), dtype=numpy.double)
    cdef double[::1] qua_ = qua
    c_ens_quantiles_variable(<int>ens.shape[0],<int>quadef.shape[0],&qua_[0],&ens[0],&quadef[0],&ens[0],&argcase)
    return qua

def quantiles_variable_weight(double[::1] ens,double[::1] quadef,double[::1] weight):
    cdef int argcase = 1
    qua  = numpy.zeros((quadef.shape[0]), dtype=numpy.double)
    cdef double[::1] qua_ = qua
    c_ens_quantiles_variable(<int>ens.shape[0],<int>quadef.shape[0],&qua_[0],&ens[0],&quadef[0],&weight[0],&argcase)
    return qua

def quantiles_vector(double[:,::1] ens,double[::1] quadef):
    cdef int argcase = 0
    qua  = numpy.zeros((ens.shape[1],quadef.shape[0]), dtype=numpy.double)
    cdef double[:,::1] qua_ = qua
    c_ens_quantiles_vector(<int>ens.shape[1],<int>ens.shape[0],<int>quadef.shape[0],&qua_[0,0],&ens[0,0],&quadef[0],&ens[0,0],&ens[0,0],&argcase)
    return qua

def quantiles_vector_weight(double[:,::1] ens,double[::1] quadef,double[::1] weight):
    cdef int argcase = 1
    qua  = numpy.zeros((ens.shape[1],quadef.shape[0]), dtype=numpy.double)
    cdef double[:,::1] qua_ = qua
    c_ens_quantiles_vector(<int>ens.shape[1],<int>ens.shape[0],<int>quadef.shape[0],&qua_[0,0],&ens[0,0],&quadef[0],&weight[0],&ens[0,0],&argcase)
    return qua

def quantiles_vector_locweight(double[:,::1] ens,double[::1] quadef,double[:,::1] locweight):
    cdef int argcase = 2
    qua  = numpy.zeros((ens.shape[1],quadef.shape[0]), dtype=numpy.double)
    cdef double[:,::1] qua_ = qua
    c_ens_quantiles_vector(<int>ens.shape[1],<int>ens.shape[0],<int>quadef.shape[0],&qua_[0,0],&ens[0,0],&quadef[0],&ens[0,0],&locweight[0,0],&argcase)
    return qua

# Public function to perform forward anamorphosis transformation
def forward(var,qua,rank=None):
    """forward(var,qua,[rank])

       Forward anamorphosis transformation

       Inputs
       ------
       var [scalar, rank-1 or rank-2 double array] : variable(s) to transform (nvar) or (nens,nvar)
       qua [rank-1 or rank-2 double array] : ensemble quantiles (nqua) or (nqua,nvar)
       rank [scalar or rank-1 double array] : rank to use to resolve probabiity concentration

       Returns
       -------
       var [scalar, rank-1 or rank-2 double array] : transformed input data (in-place)

    """
    global quaref
    if numpy.isscalar(quaref):
      raise ValueError("Quantiles of target distribution undefined in anamorphosis")

    if numpy.isscalar(var):
      if rank == None:
        forward_variable(var,qua)
      else:
        forward_variable_rank(var,qua,rank)
    elif var.ndim == 1:
      if not numpy.isscalar(rank):
        forward_vector(var,qua)
      else:
        forward_vector_rank(var,qua,rank)
    elif var.ndim == 2:
      if not numpy.isscalar(rank):
        forward_ensemble(var,qua)
      else:
        forward_ensemble_rank(var,qua,rank)

# Interfaces to corresponding FORTRAN functions
def forward_variable(double var,double[::1] qua):
    cdef int argcase = 0
    cdef double rank
    global quaref
    cdef double[::1] quaref_ = quaref
    c_ana_forward_variable(<int>qua.shape[0],&var,&qua[0],&quaref_[0],&rank,&argcase)

def forward_variable_rank(double var,double[::1] qua,double rank):
    cdef int argcase = 1
    global quaref
    cdef double[::1] quaref_ = quaref
    c_ana_forward_variable(<int>qua.shape[0],&var,&qua[0],&quaref_[0],&rank,&argcase)

def forward_vector(double[::1] var,double[:,::1] qua):
    cdef int argcase = 0
    cdef double rank
    global quaref
    cdef double[::1] quaref_ = quaref
    c_ana_forward_vector(<int>qua.shape[1],<int>qua.shape[0],&var[0],&qua[0,0],&quaref_[0],&rank,&argcase)

def forward_vector_rank(double[::1] var,double[:,::1] qua, double rank):
    cdef int argcase = 1
    global quaref
    cdef double[::1] quaref_ = quaref
    c_ana_forward_vector(<int>qua.shape[1],<int>qua.shape[0],&var[0],&qua[0,0],&quaref_[0],&rank,&argcase)

def forward_ensemble(double[:,::1] var,double[:,::1] qua):
    cdef int argcase = 0
    global quaref
    cdef double[::1] quaref_ = quaref
    c_ana_forward_ensemble(<int>var.shape[1],<int>var.shape[0],<int>qua.shape[0],&var[0,0],&qua[0,0],&quaref_[0],&var[0,0],&argcase)

def forward_ensemble_rank(double[:,::1] var,double[:,::1] qua,double[::1] rank):
    cdef int argcase = 1
    global quaref
    cdef double[::1] quaref_ = quaref
    c_ana_forward_ensemble(<int>var.shape[1],<int>var.shape[0],<int>qua.shape[0],&var[0,0],&qua[0,0],&quaref_[0],&rank[0],&argcase)

# Public function to perform backward anamorphosis transformation
def backward(var,qua):
    """backward(var,qua)

       Backward anamorphosis transformation

       Inputs
       ------
       var [scalar, rank-1 or rank-2 double array] : variable(s) to transform (nvar) or (nens,nvar)
       qua [rank-1 or rank-2 double array] : ensemble quantiles (nqua) or (nqua,nvar)

       Returns
       -------
       var [scalar, rank-1 or rank-2 double array] : transformed input data (in-place)

    """
    global quaref
    if numpy.isscalar(quaref):
      raise ValueError("Quantiles of target distribution undefined in anamorphosis")

    if numpy.isscalar(var):
      backward_variable(var,qua)
    elif var.ndim == 1:
      backward_vector(var,qua)
    elif var.ndim == 2:
      backward_ensemble(var,qua)

# Interfaces to corresponding FORTRAN functions
def backward_variable(double var,double[::1] qua):
    global quaref
    cdef double[::1] quaref_ = quaref
    c_ana_backward_variable(<int>qua.shape[0],&var,&qua[0],&quaref_[0])

def backward_vector(double[::1] var,double[:,::1] qua):
    global quaref
    cdef double[::1] quaref_ = quaref
    c_ana_backward_vector(<int>qua.shape[1],<int>qua.shape[0],&var[0],&qua[0,0],&quaref_[0])

def backward_ensemble(double[:,::1] var,double[:,::1] qua):
    global quaref
    cdef double[::1] quaref_ = quaref
    c_ana_backward_ensemble(<int>var.shape[1],<int>var.shape[0],<int>qua.shape[0],&var[0,0],&qua[0,0],&quaref_[0])

# Public function to perform forward anamorphosis transformation of observations
def forward_obs(nsmp,double[::1] obs,double[::1] obs_std,double[:,::1] obsens,double[::1] quadef):
    """anaobs = forward_obs(nsmp,obs,obs_std,obsens,quadef)

       Forward anamorphosis transformation of observations

       Inputs
       ------
       nsmp [integer] : size of transformed observation sample
       obs [rank-1 double array] : observations to transform (nobs)
       obs_std [rank-1 double array] : observation error standard deviation
       obsens [rank-2 double array] : ensemble equivalent to observations (nens,nobs)
       quadef [rank-1 double array] : definition of the quantiles used in anamorphosis (nqua)

       Returns
       -------
       anaobs [rank-2 double array] : sample of transformed observations (nsmp,nobs)

    """
    cdef int nsmp_ = nsmp
    anaobs = numpy.zeros((nsmp,obs.shape[0]), dtype=numpy.double)
    cdef double[:,::1] anaobs_ = anaobs
    global quaref
    cdef double[::1] quaref_ = quaref
    c_ana_obs(<int>obs.shape[0],<int>obsens.shape[0],<int>quadef.shape[0],nsmp_,&anaobs_[0,0],&obsens[0,0],&obs[0],&obs_std[0],&quadef[0],&quaref_[0])
    return anaobs

# Public function to perform forward anamorphosis transformation of observations (symmetric)
def forward_obs_sym(nsmp,double[::1] obs,double[::1] obs_std,double[:,::1] obsqua):
    """anaobs = forward_obs_sym(nsmp,obs,obs_std,obsqua)

       Forward anamorphosis transformation of observations (symmetric)

       Inputs
       ------
       nsmp [integer] : size of transformed observation sample
       obs [rank-1 double array] : observations to transform (nobs)
       obs_std [rank-1 double array] : observation error standard deviation
       obsqua [rank-1 double array] : quantiles of the ensemble equivalent to observations (nqua,nobs)

       Returns
       -------
       anaobs [rank-2 double array] : sample of transformed observations (nsmp,nobs)

    """
    cdef int nsmp_ = nsmp
    anaobs = numpy.zeros((nsmp,obs.shape[0]), dtype=numpy.double)
    cdef double[:,::1] anaobs_ = anaobs
    global quaref
    cdef double[::1] quaref_ = quaref
    c_ana_obs_sym(<int>obsqua.shape[1],<int>obsqua.shape[0],nsmp_,&anaobs_[0,0],&obs[0],&obs_std[0],&obsqua[0,0],&quaref_[0])
    return anaobs
