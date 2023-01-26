#---------------------------------------------------------------------
# Copyright: CNRS - Universite Grenoble Alpes
#
# Contributors : Jean-Michel Brankart
#
# Jean-Michel.Brankart@univ-grenoble-alpes.fr
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
#---------------------------------------------------------------------
"""
EnsDAM: Ensamble Data Assimilation modules

Available modules:
    ensdam.ensanam : ensemble anamorphosis transformation
    ensdam.ensaugm : ensemble augmentation
    ensdam.ensscores : ensemble scores
    ensdam.ensstat : ensemble statistics
    ensdam.ensupdate : ensemble observational update
    ensdam.interptools : interpolation tools
    ensdam.obserror : observation error
    ensdam.stochtools : stochastic tools
    ensdam.transpho : scale separation (by projection on the spherical harmonics)
"""
import sys
sys.path.append("./lib")

# Import EnsDAM modules (prepared by f2py)
import anaobs
import anaqua
import anatra
import anautil
import covariance
import ensaugm
import interp
import obserror
import meanstd
import schurprod
import score_crps
import score_entropy
import score_optimality
import score_rcrv
import spharea
import sphylm
import stoanam
import stogprod
import storfg
import storng
import stotge
import stoutil

# Import ctypes and numpy
from ctypes import *
import numpy

# Get MPI public variables from the FORTRAN modules
mpi=False
if mpi:
    mpi_comm_ensaugm = ensaugm.ensdam_ensaugm.mpi_comm_ensaugm
    mpi_comm_meanstd = meanstd.ensdam_meanstd.mpi_comm_meanstd
    mpi_comm_score_crps = score_crps.ensdam_score_crps.mpi_comm_score_crps
    mpi_comm_score_entropy = score_entropy.ensdam_score_entropy.mpi_comm_score_entropy
    mpi_comm_score_rcrv = score_rcrv.ensdam_score_rcrv.mpi_comm_score_rcrv
    mpi_comm_sphylm = sphylm.ensdam_sphylm.mpi_comm_sphylm
    mpi_comm_storfg = storfg.ensdam_storfg.mpi_comm_storfg

# Get default public variables from the FORTRAN modules
# From anautil
anautil_reference_cdf = anautil.ensdam_anautil.anautil_reference_cdf
anautil_a = anautil.ensdam_anautil.anautil_a
anautil_b = anautil.ensdam_anautil.anautil_b
# From ensaugm
ensaugm_chain_index  = ensaugm.ensdam_ensaugm.ensaugm_chain_index
ensaugm_with_renormalization = ensaugm.ensdam_ensaugm.ensaugm_with_renormalization
# From schurprod
schurprod_gpmin = schurprod.ensdam_schurprod.schurprod_gpmin
schurprod_gpmax = schurprod.ensdam_schurprod.schurprod_gpmax
schurprod_precompute = schurprod.ensdam_schurprod.schurprod_precompute
schurprod_tablesize = schurprod.ensdam_schurprod.schurprod_tablesize
# From score_crps
crps_missing_value = score_crps.ensdam_score_crps.crps_missing_value
# From score_entropy
score_entropy_base = score_entropy.ensdam_score_entropy.score_entropy_base
# From score_rcrv
rcrv_number_of_quantiles = score_rcrv.ensdam_score_rcrv.rcrv_number_of_quantiles
rcrv_with_anamorphosis = score_rcrv.ensdam_score_rcrv.rcrv_with_anamorphosis
rcrv_missing_value = score_rcrv.ensdam_score_rcrv.rcrv_missing_value
# From sphylm
regr_epsilon = sphylm.ensdam_sphylm.regr_epsilon
regr_overlap = sphylm.ensdam_sphylm.regr_overlap
regr_maxiter = sphylm.ensdam_sphylm.regr_maxiter
regr_type = sphylm.ensdam_sphylm.regr_type
regr_rho = sphylm.ensdam_sphylm.regr_rho
regr_maxbloc = sphylm.ensdam_sphylm.regr_maxbloc
external_vector_decomposition = sphylm.ensdam_sphylm.external_vector_decomposition
# From storfg:
storfg_ylm_resolution = storfg.ensdam_storfg.storfg_ylm_resolution
# From stoutil:
nominal_accuracy = stoutil.ensdam_stoutil.nominal_accuracy
accuracy = stoutil.ensdam_stoutil.accuracy
maxiter = stoutil.ensdam_stoutil.maxiter

# API general functions
def _iarrayint(o):
     #array = numpy.ascontiguousarray(o, numpy.int32)
     array = numpy.asfortranarray(o, numpy.int32)
     return  array

def _iarraydouble(o):
     #array = numpy.ascontiguousarray(o, numpy.float64)
     array = numpy.asfortranarray(o, numpy.float64)
     return  array

# API for EnsDAM FORTRAN modules
class ensanam:
    """
    The purpose of EnsAnam is to provide tools
    to apply anamorphosis transformation to ensemble simulations.
    The objective is to transform the marginal distribution of all variables
    to the same distribution defined by the user
    (usually a normalized Gaussian distribution).

    Available variables:
    anautil_reference_cdf,anautil_a,anautil_b

    Available functions:
    ens_quantiles(),ana_forward(),ana_backward(),
    ana_obs(),ana_obs_sym(),ana_util_quaref()
    """

    @staticmethod
    def ens_quantiles(ens,quadef,enswei=(),ensweiloc=()):
        """
        qua = ens_quantiles(ens,quadef,[enswei,ensweiloc])

        Compute quantiles from input ensemble

        Parameters
        ----------
        ens : input rank-1 or rank-2 array('d')
        quadef : input rank-1 array('d')

        Other Parameters
        ----------------
        enswei : input rank-1 array('d')
        ensweiloc : input rank-1 or rank-2 array('d')

        Returns
        -------
        qua : rank-1 or rank-2 array('d')
        """
        if ens.ndim == 1:
            m = ens.shape[0] # ensemble size
            if enswei == ():
                enswei=numpy.ones([m])
            api_ens = _iarraydouble(ens)
            api_quadef = _iarraydouble(quadef)
            api_enswei = _iarraydouble(enswei)
            qua = anaqua.ensdam_anaqua.ens_quantiles_variable(api_ens,api_quadef,api_enswei)
        elif ens.ndim == 2:
            n = ens.shape[0] # number of variables
            m = ens.shape[1] # ensemble size
            if enswei == ():
                enswei=numpy.ones([m])
            if ensweiloc == ():
                enswei=numpy.ones([n,m])
            api_ens = _iarraydouble(ens)
            api_quadef = _iarraydouble(quadef)
            api_enswei = _iarraydouble(enswei)
            api_ensweiloc = _iarraydouble(ensweiloc)
            qua = anaqua.ensdam_anaqua.ens_quantiles_vector(api_ens,api_quadef,api_enswei,api_ensweiloc)
        else:
            raise ValueError("Invalid array dimension: ",ens.ndim)
        return qua

    @staticmethod
    def ana_forward(ens,qua,quaref,rank):
        """
        ana_forward_ensemble(ens,qua,quaref,[rank])

        Perform forward anamorphosis transformation

        Parameters
        ----------
        ens : in/output rank-2 array('d') with bounds (f2py_ens_d0,f2py_ens_d1)
        qua : input rank-2 array('d') with bounds (size(ens,1),f2py_qua_d1)
        quaref : input rank-1 array('d') with bounds (size(qua,2))
        rank : input rank-1 array('d') with bounds (f2py_rank_bn_d0)
        """
        api_qua = _iarraydouble(qua)
        api_quaref = _iarraydouble(quaref)
        if numpy.isscalar(ens):
            anatra.ensdam_anatra.ana_forward_variable(ens,api_qua,api_quaref,rank)
        elif ens.ndim == 1:
            m = ens.shape[0] # ensemble size
            api_ens = _iarraydouble(ens)
            anatra.ensdam_anatra.ana_forward_ensemble(api_ens,api_qua,api_quaref,rank)
            ens = api_ens
        elif ens.ndim == 2:
            n = ens.shape[0] # number of variables
            m = ens.shape[1] # ensemble size
            api_ens = _iarraydouble(ens)
            api_rank = _iarraydouble(rank)
            anatra.ensdam_anatra.ana_forward_ensemble(api_ens,api_qua,api_quaref,rank)
            ens = api_ens
        else:
            raise ValueError("Invalid array dimension: ",ens.ndim)

    @staticmethod
    def ana_backward(ens,qua,quaref):
        """
        ana_backward_ensemble(ens,qua,quaref)

        Perform backward anamorphosis transformation

        Parameters
        ----------
        ens : in/output rank-2 array('d') with bounds (f2py_ens_d0,f2py_ens_d1)
        qua : input rank-2 array('d') with bounds (size(ens,1),f2py_qua_d1)
        quaref : input rank-1 array('d') with bounds (size(qua,2))
        """
        api_qua = _iarraydouble(qua)
        api_quaref = _iarraydouble(quaref)
        if numpy.isscalar(ens):
            anatra.ensdam_anatra.ana_backward_variable(ens,api_qua,api_quaref)
        elif ens.ndim == 1:
            m = ens.shape[0] # ensemble size
            api_ens = _iarraydouble(ens)
            anatra.ensdam_anatra.ana_backward_ensemble(api_ens,api_qua,api_quaref)
            ens = api_ens
        elif ens.ndim == 2:
            n = ens.shape[0] # number of variables
            m = ens.shape[1] # ensemble size
            api_ens = _iarraydouble(ens)
            anatra.ensdam_anatra.ana_backward_ensemble(api_ens,api_qua,api_quaref)
            ens = api_ens
        else:
            raise ValueError("Invalid array dimension: ",ens.ndim)

    @staticmethod
    def ana_obs(obsens,obs,obserror,quadef,quaref):
        """
        anaobs = ana_obs(obsens,obs,obserror,quadef,quaref)

        Transform observation probability distribution

        Parameters
        ----------
        obsens : input rank-2 array('d') with bounds (size(anaobs,1),f2py_obsens_d1)
        obs : input rank-1 array('d') with bounds (size(anaobs,1))
        obserror : input rank-1 array('d') with bounds (size(anaobs,1))
        quadef : input rank-1 array('d') with bounds (f2py_quadef_d0)
        quaref : input rank-1 array('d') with bounds (size(quadef,1))

        Returns
        -------
        anaobs : rank-2 array('d') with bounds (f2py_anaobs_d0,f2py_anaobs_d1)
        """
        api_obsens = _iarraydouble(obsens)
        api_obs = _iarraydouble(obs)
        api_obserror = _iarraydouble(obserror)
        api_quadef = _iarraydouble(quadef)
        api_quaref = _iarraydouble(quaref)
        anaobs = anaobs.ensdam_anaobs.ana_obs(api_obsens,api_obs,api_obserror,api_quadef,api_quaref)
        return anaobs

    @staticmethod
    def ana_obs_sym(obs,obserror,obsqua,quaref):
        """
        anaobs = ana_obs_sym(obs,obserror,obsqua,quaref)

        Transform observation probability distribution
        (simplified for symmetric observation error distribution)

        Parameters
        ----------
        obs : input rank-1 array('d') with bounds (size(anaobs,1))
        obserror : input rank-1 array('d') with bounds (size(anaobs,1))
        obsqua : input rank-2 array('d') with bounds (size(anaobs,1),f2py_obsqua_d1)
        quaref : input rank-1 array('d') with bounds (size(obsqua,2))

        Returns
        -------
        anaobs : rank-2 array('d') with bounds (f2py_anaobs_d0,f2py_anaobs_d1)
        """
        api_obs = _iarraydouble(obs)
        api_obserror = _iarraydouble(obserror)
        api_obsqua = _iarraydouble(obsqua)
        api_quaref = _iarraydouble(quaref)
        anaobs = anaobs.ensdam_anaobs.ana_obs_sym(api_obs,api_obserror,api_obsqua,api_quaref)
        return anaobs

    @staticmethod
    def ana_util_quaref(quadef):
        """
        quaref = ana_util_quaref(quadef)

        Compute quantiles of reference distribution

        Parameters
        ----------
        quadef : input rank-1 array('d') : definition of quantiles

        Returns
        -------
        quaref : rank-1 array('d') : quantiles of reference distribution
        """
        # Update module public variables
        anautil.ensdam_anautil.anautil_reference_cdf = anautil_reference_cdf
        anautil.ensdam_anautil.anautil_a = anautil_a
        anautil.ensdam_anautil.anautil_b = anautil_b
        # Call to Fortran function
        api_quadef = _iarraydouble(quadef)
        quaref = anautil.ensdam_anautil.ana_util_quaref(api_quadef)
        return quaref

class ensaugm:
    """
    The purpose of EnsAugm is to provide tools
    to augment an ensemble simulation with new members.
    The objective is to mitigate the effects of undersampling
    by an artifical increase of the ensemble size.

    Available variables:
    ensaugm_chain_index,ensaugm_with_renormalization,mpi_comm_ensaugm

    Available functions:
    sample_augmented_ensemble(),newproduct(),getproduct()
    """

    @staticmethod
    def sample_augmented_ensemble(maxchain,augens,ens,multiplicity):
        """
        sample_augmented_ensemble(maxchain,augens,ens,multiplicity)

        Ensemble augmentation by Schur product with large scale patterns

        Parameters
        ----------
        maxchain : input int
        augens : in/output rank-2 array('d') with bounds (f2py_augens_d0,f2py_augens_d1)
        ens : input rank-3 array('d') with bounds (size(augens,1),f2py_ens_d1,f2py_ens_d2)
        multiplicity : input rank-1 array('i') with bounds (f2py_multiplicity_d0)
        """
        # Update module public variables
        ensaugm.ensdam_ensaugm.ensaugm_chain_index = ensaugm_chain_index
        ensaugm.ensdam_ensaugm.ensaugm_with_renormalization = ensaugm_with_renormalization
        if mpi:
            ensaugm.ensdam_ensaugm.mpi_comm_ensaugm = mpi_comm_ensaugm
        # Call to Fortran function
        api_augens = _iarraydouble(augens)
        api_ens = _iarraydouble(ens)
        api_multiplicity = _iarrayint(multiplicity)
        ensaugm.ensdam_ensaugm.sample_augmented_ensemble(maxchain,api_augens,api_ens,api_multiplicity)
        augens = api_augens

    @staticmethod
    def newproduct(ens,multiplicity):
        """
        new,sample = newproduct(ens,multiplicity)

        Compute new multiple Schur product

        Parameters
        ----------
        ens : input rank-3 array('d') with bounds (size(new_bn,1),f2py_ens_d1,f2py_ens_d2)
        multiplicity : input rank-1 array('i') with bounds (f2py_multiplicity_d0)

        Returns
        -------
        new : rank-1 array('d') with bounds (f2py_new_bn_d0)
        sample : rank-1 array('i') with bounds (f2py_sample_d0)
        """
        # Update module public variables
        ensaugm.ensdam_ensaugm.ensaugm_chain_index = ensaugm_chain_index
        ensaugm.ensdam_ensaugm.ensaugm_with_renormalization = ensaugm_with_renormalization
        if mpi:
            ensaugm.ensdam_ensaugm.mpi_comm_ensaugm = mpi_comm_ensaugm
        # Call to Fortran function
        api_ens = _iarraydouble(ens)
        api_multiplicity = _iarrayint(multiplicity)
        new,sample = ensaugm.ensdam_ensaugm.newproduct(api_ens,api_multiplicity)
        return new,sample

    @staticmethod
    def getproduct(ens,multiplicity,sample):
        """
        new = newproduct(ens,multiplicity,sample)

        Compute specified multiple Schur product

        Parameters
        ----------
        ens : input rank-3 array('d') with bounds (size(new_bn,1),f2py_ens_d1,f2py_ens_d2)
        multiplicity : input rank-1 array('i') with bounds (f2py_multiplicity_d0)
        sample : rank-1 array('i') with bounds (f2py_sample_d0)

        Returns
        -------
        new : rank-1 array('d') with bounds (f2py_new_bn_d0)
        """
        # Update module public variables
        ensaugm.ensdam_ensaugm.ensaugm_chain_index = ensaugm_chain_index
        ensaugm.ensdam_ensaugm.ensaugm_with_renormalization = ensaugm_with_renormalization
        if mpi:
            ensaugm.ensdam_ensaugm.mpi_comm_ensaugm = mpi_comm_ensaugm
        # Call to Fortran function
        api_ens = _iarraydouble(ens)
        api_multiplicity = _iarrayint(multiplicity)
        api_sample = _iarrayint(sample)
        new = ensaugm.ensdam_ensaugm.newproduct(api_ens,api_multiplicity,api_sample)
        return new,sample

class ensscores:
    """
    The purpose of EnsScores is to provide tools
    to compute probabilistic scores of ensemble simulations.
    The scores evaluates the reliability, resolution and optimality
    of the simulation by comparison to verification data.

    Available variables:
    crps_missing_value,mpi_comm_score_crps,mpi_comm_score_rcrv,
    rcrv_number_of_quantiles,rcrv_with_anamorphosis,rcrv_missing_value,
    mpi_comm_score_entropy,score_entropy_base

    Available functions:
    crps_score(), crps_cumul(),crps_final(), rcrv_score(),rcrv_cumul(),
    optimality_score(),optimality_cumul(),events_score(),events_relative_entropy(),
    events_cross_entropy(),events_entropy(),events_probability()
    """

    @staticmethod
    def crps_score(ens,verif,partition=()):
        """
        crps,reliability,resolution = crps_score(ens,verif,partition)

        Compute CRPS score (with option to partition the data)

        Parameters
        ----------
        ens : input rank-2 array('d') with bounds (f2py_ens_d0,f2py_ens_d1)
        verif : input rank-1 array('d') with bounds (size(ens,1))
        partition : input rank-1 array('i') with bounds (size(ens,1))

        Returns
        -------
        crps : rank-1 array('d') with bounds (f2py_crps_d0)
        reliability : rank-1 array('d') with bounds (size(crps,1))
        resolution : rank-1 array('d') with bounds (size(crps,1))
        """
        # Update module public variables
        score_crps.ensdam_score_crps.crps_missing_value = crps_missing_value
        if mpi:
            score_crps.ensdam_score_crps.mpi_comm_score_crps = mpi_comm_score_crps

        # Apply requested function
        if partition == ():
            api_ens = _iarraydouble(ens)
            api_verif = _iarraydouble(verif)
            crps,reliability,resolution =  score_crps.ensdam_score_crps.crps_score_global(api_ens,api_verif)
        else:
            api_ens = _iarraydouble(ens)
            api_verif = _iarraydouble(verif)
            api_parttion = _iarraydouble(parttion)
            crps,reliability,resolution =  score_crps.ensdam_score_crps.crps_score_partition(api_ens,api_verif,api_parttion)

        return crps,reliability,resolution

    @staticmethod
    def crps_cumul(ens,a,aa,bb):
        """
        crps_cumul(ens,a,aa,bb)

        Accumulate data to prepare the final computation of the score

        Parameters
        ----------
        ens : input rank-1 array('d') with bounds (f2py_ens_d0)
        a : input float
        aa : in/output rank-1 array('d') with bounds (1)
        bb : in/output rank-1 array('d') with bounds (1)
        """
        # Update module public variables
        score_crps.ensdam_score_crps.crps_missing_value = crps_missing_value
        if mpi:
            score_crps.ensdam_score_crps.mpi_comm_score_crps = mpi_comm_score_crps

    @staticmethod
    def crps_final(aa,bb):
        """
        reli,resol,crps = crps_final(aa,bb)

        Compute final score from accumulated data

        Parameters
        ----------
        aa : in/output rank-1 array('d') with bounds (1)
        bb : in/output rank-1 array('d') with bounds (1)

        Returns
        -------
        reli : float
        resol : float
        crps : float
        """
        # Update module public variables
        score_crps.ensdam_score_crps.crps_missing_value = crps_missing_value
        if mpi:
            score_crps.ensdam_score_crps.mpi_comm_score_crps = mpi_comm_score_crps

    @staticmethod
    def rcrv_score(ens,verif,partition=()):
        """
        ens_bias,ens_spread = rcrv_score_partition(ens,verif,partition)

        Compute RCRV score (with option to partition the data)

        Parameters
        ----------
        ens : input rank-2 array('d') with bounds (f2py_ens_d0,f2py_ens_d1)
        verif : input rank-1 array('d') with bounds (size(ens,1))
        partition : input rank-1 array('i') with bounds (size(ens,1))

        Returns
        -------
        ens_bias : rank-1 array('d') with bounds (f2py_ens_bias_d0)
        ens_spread : rank-1 array('d') with bounds (size(ens_bias,1))
        """
        # Update module public variables
        score_rcrv.ensdam_score_rcrv.rcrv_number_of_quantiles = rcrv_number_of_quantiles
        score_rcrv.ensdam_score_rcrv.rcrv_with_anamorphosis = rcrv_with_anamorphosis
        score_rcrv.ensdam_score_rcrv.rcrv_missing_value = rcrv_missing_value
        if mpi:
            score_rcrv.ensdam_score_rcrv.mpi_comm_score_rcrv = mpi_comm_score_rcrv

        # Apply requested function
        if partition == ():
            api_ens = _iarraydouble(ens)
            api_verif = _iarraydouble(verif)
            ens_bias,ens_spread =  score_rcrv.ensdam_score_rcrv.rcrv_score_global(api_ens,api_verif)
        else:
            api_ens = _iarraydouble(ens)
            api_verif = _iarraydouble(verif)
            api_parttion = _iarraydouble(parttion)
            ens_bias,ens_spread =  score_rcrv.ensdam_score_rcrv.rcrv_score_partition(api_ens,api_verif,api_parttion)

        return ens_bias,ens_spread

    @staticmethod
    def rcrv_cumul(e,a,idx,mean,sqrs):
        """
        rcrv_cumul(e,a,idx,mean,sqrs)

        Accumulate data to prepare the final computation of the score

        Parameters
        ----------
        e : input rank-1 array('d') with bounds (f2py_e_d0)
        a : input float
        idx : input int
        mean : in/output rank-0 array(float,'d')
        sqrs : in/output rank-0 array(float,'d')
        """
        # Update module public variables
        score_rcrv.ensdam_score_rcrv.rcrv_number_of_quantiles = rcrv_number_of_quantiles
        score_rcrv.ensdam_score_rcrv.rcrv_with_anamorphosis = rcrv_with_anamorphosis
        score_rcrv.ensdam_score_rcrv.rcrv_missing_value = rcrv_missing_value
        if mpi:
            score_rcrv.ensdam_score_rcrv.mpi_comm_score_rcrv = mpi_comm_score_rcrv

    @staticmethod
    def optimality_score(ens,obs,cdf_obs,cdf_obs_extra_args=(),partition=()):
        """
        ens_optimality,ens_optimality_bias,ens_optimality_spread = optimality_score(ens,obs,cdf_obs,[cdf_obs_extra_args,partition])

        Compute optimality score (with option to partition the data)

        Parameters
        ----------
        ens : input rank-2 array('d') with bounds (f2py_ens_d0,f2py_ens_d1)
        obs : input rank-1 array('d') with bounds (f2py_obs_d0)
        cdf_obs : call-back function => callback_cdf_obs

        Other Parameters
        ----------------
        cdf_obs_extra_args : input tuple, optional
            Default: ()
        partition : input rank-1 array('i') with bounds (f2py_partition_d0)

        Returns
        -------
        ens_optimality : rank-1 array('d') with bounds (f2py_ens_optimality_d0)
        ens_optimality_bias : rank-1 array('d') with bounds (f2py_ens_optimality_bias_d0)
        ens_optimality_spread : rank-1 array('d') with bounds (f2py_ens_optimality_spread_d0)

        Notes
        -----
        Call-back functions::

          def callback_cdf_obs(o,y,obs_idx): return callback_cdf_obs
          Required arguments:
            o : input float
            y : input float
            obs_idx : input int
          Return objects:
            callback_cdf_obs : float
        """

    @staticmethod
    def optimality_cumul(e,o,idx,mean,sqrs,cdf_obs,cdf_obs_extra_args=()):
        """
        optimality_cumul(e,o,idx,mean,sqrs,cdf_obs,[cdf_obs_extra_args])

        Accumulate data to prepare the final computation of the score

        Parameters
        ----------
        e : input rank-1 array('d') with bounds (f2py_e_d0)
        o : input float
        idx : input int
        mean : in/output rank-0 array(float,'d')
        sqrs : in/output rank-0 array(float,'d')
        cdf_obs : call-back function => callback_cdf_obs

        Other Parameters
        ----------------
        cdf_obs_extra_args : input tuple, optional
            Default: ()

        Notes
        -----
        Call-back functions::

          def callback_cdf_obs(o,y,obs_idx): return callback_cdf_obs
          Required arguments:
            o : input float
            y : input float
            obs_idx : input int
          Return objects:
            callback_cdf_obs : float
        """

    @staticmethod
    def events_score(ens,pref,events_outcome,events_outcome_extra_args=()):
        """
        score = events_score(ens,pref,events_outcome,[events_outcome_extra_args])

        Compute ensemble score for the required events

        Parameters
        ----------
        ens : input rank-2 array('d') with bounds (f2py_ens_d0,f2py_ens_d1)
        pref : input rank-2 array('d') with bounds (f2py_pref_d0,f2py_pref_d1)
        events_outcome : call-back function

        Other Parameters
        ----------------
        events_outcome_extra_args : input tuple, optional
            Default: ()

        Returns
        -------
        score : rank-1 array('d') with bounds (f2py_score_d0)

        Notes
        -----
        Call-back functions::

          def events_outcome(member): return outcome
          Required arguments:
            member : input rank-1 array('d') with bounds (:)
          Return objects:
            outcome : rank-1 array('i') with bounds (:)
        """
        # Update module public variables
        score_entropy.ensdam_score_entropy.score_entropy_base = score_entropy_base
        if mpi:
            score_entropy.ensdam_score_entropy.mpi_comm_score_entropy = mpi_comm_score_entropy

    @staticmethod
    def events_relative_entropy(ens,pref,events_outcome,events_outcome_extra_args=()):
        """
        relative_entropy = events_relative_entropy(ens,pref,events_outcome,[events_outcome_extra_args])

        Compute relative entropy between ensemble distribution and reference distribution

        Parameters
        ----------
        ens : input rank-2 array('d') with bounds (f2py_ens_d0,f2py_ens_d1)
        pref : input rank-2 array('d') with bounds (f2py_pref_d0,f2py_pref_d1)
        events_outcome : call-back function

        Other Parameters
        ----------------
        events_outcome_extra_args : input tuple, optional
            Default: ()

        Returns
        -------
        relative_entropy : rank-1 array('d') with bounds (f2py_relative_entropy_d0)

        Notes
        -----
        Call-back functions::

          def events_outcome(member): return outcome
          Required arguments:
            member : input rank-1 array('d') with bounds (:)
          Return objects:
            outcome : rank-1 array('i') with bounds (:)
        """
        # Update module public variables
        score_entropy.ensdam_score_entropy.score_entropy_base = score_entropy_base
        if mpi:
            score_entropy.ensdam_score_entropy.mpi_comm_score_entropy = mpi_comm_score_entropy

    @staticmethod
    def events_cross_entropy(ens,pref,events_outcome,events_outcome_extra_args=()):
        """
        cross_entropy,entropy = events_cross_entropy(ens,pref,events_outcome,[events_outcome_extra_args])

        Compute cross entropy between ensemble distribution and reference distribution

        Parameters
        ----------
        ens : input rank-2 array('d') with bounds (f2py_ens_d0,f2py_ens_d1)
        pref : input rank-2 array('d') with bounds (f2py_pref_d0,f2py_pref_d1)
        events_outcome : call-back function

        Other Parameters
        ----------------
        events_outcome_extra_args : input tuple, optional
            Default: ()

        Returns
        -------
        cross_entropy : rank-1 array('d') with bounds (f2py_cross_entropy_d0)
        entropy : rank-1 array('d') with bounds (f2py_entropy_d0)

        Notes
        -----
        Call-back functions::

          def events_outcome(member): return outcome
          Required arguments:
            member : input rank-1 array('d') with bounds (:)
          Return objects:
            outcome : rank-1 array('i') with bounds (:)
        """
        # Update module public variables
        score_entropy.ensdam_score_entropy.score_entropy_base = score_entropy_base
        if mpi:
            score_entropy.ensdam_score_entropy.mpi_comm_score_entropy = mpi_comm_score_entropy

    @staticmethod
    def events_entropy(ens,events_outcome,number_outcome,events_outcome_extra_args=()):
        """
        entropy = events_entropy(ens,events_outcome,number_outcome,[events_outcome_extra_args])

        Compute entropy of ensemble distribution for the required events

        Parameters
        ----------
        ens : input rank-2 array('d') with bounds (f2py_ens_d0,f2py_ens_d1)
        events_outcome : call-back function
        number_outcome : input int

        Other Parameters
        ----------------
        events_outcome_extra_args : input tuple, optional
            Default: ()

        Returns
        -------
        entropy : rank-1 array('d') with bounds (f2py_entropy_d0)

        Notes
        -----
        Call-back functions::

          def events_outcome(member): return outcome
          Required arguments:
            member : input rank-1 array('d') with bounds (:)
          Return objects:
            outcome : rank-1 array('i') with bounds (:)
        """
        # Update module public variables
        score_entropy.ensdam_score_entropy.score_entropy_base = score_entropy_base
        if mpi:
            score_entropy.ensdam_score_entropy.mpi_comm_score_entropy = mpi_comm_score_entropy

    @staticmethod
    def events_probability(ens,events_outcome,events_outcome_extra_args=()):
        """
        pens = events_probability(ens,events_outcome,[events_outcome_extra_args])

        Compute events marginal probability distributions from the ensemble

        Parameters
        ----------
        ens : input rank-2 array('d') with bounds (f2py_ens_d0,f2py_ens_d1)
        events_outcome : call-back function

        Other Parameters
        ----------------
        events_outcome_extra_args : input tuple, optional
            Default: ()

        Returns
        -------
        pens : rank-2 array('d') with bounds (f2py_pens_d0,f2py_pens_d1)

        Notes
        -----
        Call-back functions::

          def events_outcome(member): return outcome
          Required arguments:
            member : input rank-1 array('d') with bounds (:)
          Return objects:
            outcome : rank-1 array('i') with bounds (:)
        """
        # Update module public variables
        score_entropy.ensdam_score_entropy.score_entropy_base = score_entropy_base
        if mpi:
            score_entropy.ensdam_score_entropy.mpi_comm_score_entropy = mpi_comm_score_entropy

class ensstat:
    """
    The purpose of EnsStat is to provide tools
    to compute basic ensemble statistics.
    This includes the computation of ensemble mean and standard deviation,*
    and the computation of ensemble covariance or correlation structure.
    """

    @staticmethod
    def ensemble_meanstd(ens,weight=()):
        """
        mean,std = ensemble_meanstd_vector(ens,[weight])

        Compute mean and standard deviation from input ensemble

        Parameters
        ----------
        ens : input rank-1 or rank-2 array('d')

        Other Parameters
        ----------------
        weight : input rank-1 or rank-2 array('d')

        Returns
        -------
        mean : float or rank-1 array('d')
        std : float or rank-1 array('d')

        """
        # Update module public variables
        if mpi:
            meanstd.ensdam_meanstd.mpi_comm_meanstd = mpi_comm_meanstd
        # Check optional argument
        if weight == ():
            weight=numpy.ones(ens.shape)
        # Apply requested function
        api_ens = _iarraydouble(ens)
        api_weight = _iarraydouble(weight)
        if api_ens.ndim == 1:
            mean, std = meanstd.ensdam_meanstd.ensemble_meanstd_variable(api_ens,api_weight)
        elif ens.ndim == 2:
            mean, std = meanstd.ensdam_meanstd.ensemble_meanstd_vector(api_ens,api_weight)
        else:
            raise ValueError("Invalid array dimension: ",ens.ndim)
        return mean,std

    @staticmethod
    def update_meanstd(vct,mean,idx=(),msqra=(),weight=(),weightsum=()):
        """
        update_meanstd_vector(vct,mean,[dx,msqra],weight,weightsum)

        Update mean and mean squared anomalies

        Parameters
        ----------
        vct : input float or rank-1 array('d')
        mean : in/output float or rank-1 array('d')

        Other Parameters
        ----------------
        idx : input int
        weight : input float or rank-1 array('d')
        weightsum : in/output float or rank-1 array(float,'d')
        msqra : in/output float or rank-1 array('d')
        """
        print "Warning: python interface only partially tested !!"
        # Update module public variables
        if mpi:
            meanstd.ensdam_meanstd.mpi_comm_meanstd = mpi_comm_meanstd 
        # Check the presence of weights
        if weight == ():
            # Check type of input
            if numpy.isscalar(vct):
                if msqra == ():
                    meanstd.ensdam_meanstd.update_meanstd_variable(vct,idx,mean)
                else:
                    meanstd.ensdam_meanstd.update_meanstd_variable(vct,idx,mean,msqra)
            else:
                api_vct = _iarraydouble(vct)
                api_mean = _iarraydouble(mean)
                if api_vct.ndim != 1:
                    raise ValueError("Invalid array dimension: ",api_vct.ndim)
                else:
                    if msqra == ():
                        msqra=numpy.zeros(api_mean.shape)
                    api_msqra = _iarraydouble(msqra)
                    meanstd.ensdam_meanstd.update_meanstd_vector(api_vct,idx,api_mean,api_msqra)
        else:
            # Check type of input
            if numpy.isscalar(vct):
                if msqra == ():
                    meanstd.ensdam_meanstd.update_meanstd_variable_weight(vct,weight,weightsum,mean)
                else:
                    meanstd.ensdam_meanstd.update_meanstd_variable_weight(vct,weight,weightsum,mean,msqra)
            else:
                api_vct = _iarraydouble(vct)
                api_mean = _iarraydouble(mean)
                api_weight = _iarraydouble(weight)
                api_weightsum = _iarraydouble(weightsum)
                if api_vct.ndim != 1:
                    raise ValueError("Invalid array dimension: ",api_vct.ndim)
                else:
                    if msqra == ():
                        msqra=numpy.zeros(api_mean.shape)
                    api_msqra = _iarraydouble(msqra)
                    meanstd.ensdam_meanstd.update_meanstd_vector_weight(api_vct,api_weight,api_weightsum,mean,api_msqra)

    @staticmethod
    def ensemble_correlation(ens,ensref,correl,weight=()):
        """
        correl = ensemble_correlation(ens,ensref,[weight])

        Compute correlation from input ensemble

        Parameters
        ----------
        ens : input rank-2 array('d') with bounds (f2py_ens_d0,f2py_ens_d1)
        ensref : input rank-1 array('d') with bounds (size(ens,2))

        Other Parameters
        ----------------
        weight : input rank-1 array('d') with bounds (size(ens,2))

        Returns
        -------
        correl : rank-1 array('d') with bounds (size(ens,1))
        """
        print "Warning: python interface untested !!"
        # Check optional argument
        if weight == ():
            weight=numpy.ones(ensref.shape)
        # Apply requested function
        api_ens = _iarraydouble(ens)
        api_ensref = _iarraydouble(ensref)
        api_weight = _iarraydouble(weight)
        correl = covariance.ensdam_covariance.ensemble_correlation(api_ens,api_ensref,api_weight)
        return correl

    @staticmethod
    def ensemble_covariance(ens,ensref,correl,weight=()):
        """
        covariance = ensemble_correlation(ens,ensref,[weight])

        Compute covariance from input ensemble

        Parameters
        ----------
        ens : input rank-2 array('d') with bounds (f2py_ens_d0,f2py_ens_d1)
        ensref : input rank-1 array('d') with bounds (size(ens,2))

        Other Parameters
        ----------------
        weight : input rank-1 array('d') with bounds (size(ens,2))

        Returns
        -------
        covariance : rank-1 array('d') with bounds (size(ens,1))
        """
        print "Warning: python interface untested !!"
        # Check optional argument
        if weight == ():
            weight=numpy.ones(ensref.shape)
        # Apply requested function
        api_ens = _iarraydouble(ens)
        api_ensref = _iarraydouble(ensref)
        api_weight = _iarraydouble(weight)
        covariance = covariance.ensdam_covariance.ensemble_correlation(api_ens,api_ensref,api_weight)
        return covariance

    @staticmethod
    def ensemble_representer(ens,ensref,correl,weight=()):
        """
        representer = ensemble_representer(ens,ensref,[weight])

        Compute representer from input ensemble

        Parameters
        ----------
        ens : input rank-2 array('d') with bounds (f2py_ens_d0,f2py_ens_d1)
        ensref : input rank-1 array('d') with bounds (size(ens,2))

        Other Parameters
        ----------------
        weight : input rank-1 array('d') with bounds (size(ens,2))

        Returns
        -------
        representer : rank-1 array('d') with bounds (size(ens,1))
        """
        print "Warning: python interface untested !!"
        # Check optional argument
        if weight == ():
            weight=numpy.ones(ensref.shape)
        # Apply requested function
        api_ens = _iarraydouble(ens)
        api_ensref = _iarraydouble(ensref)
        api_weight = _iarraydouble(weight)
        representer = covariance.ensdam_covariance.ensemble_representer(api_ens,api_ensref,api_weight)
        return representer

    @staticmethod
    def update_meancov(vct,varref,mean,meanref,mproda,memidx=(),weight=(),weightsum=()):
        """
        update_meancov(vct,varref,mean,meanref,mproda,[memidx,weight,weightsum])

        update mean and mean product of anomalies
        -using Welford's incremental algorithm (equal weights)
        -using West's incremental algorithm (unequal weights)

        Parameters
        ----------
        vct : input rank-1 array('d') with bounds (f2py_vct_d0)
        varref : input float
        mean : in/output rank-1 array('d') with bounds (size(vct,1))
        meanref : in/output float
        mproda : in/output rank-1 array('d') with bounds (size(vct,1))
        memidx : input int
        weight : input float
        weightsum : input float
        """
        print "Warning: python interface untested !!"
        if weight == ():
            # Check type of input
            if numpy.isscalar(vct):
                meanstd.ensdam_meanstd.update_meancov_variable(vct,varref,idx,mean,meanref,mproda)
            else:
                api_vct = _iarraydouble(vct)
                api_mean = _iarraydouble(mean)
                api_mproda = _iarraydouble(mproda)
                if api_vct.ndim != 1:
                    raise ValueError("Invalid array dimension: ",api_vct.ndim)
                else:
                    meanstd.ensdam_meanstd.update_meancov_vector(api_vct,varref,idx,api_mean,meanref,api_mproda)
        else:
            # Check type of input
            if numpy.isscalar(vct):
                meanstd.ensdam_meanstd.update_meancov_variable_weight(vct,varref,weight,weightsum,mean,meanref,mproda)
            else:
                api_vct = _iarraydouble(vct)
                api_mean = _iarraydouble(mean)
                api_mproda = _iarraydouble(mproda)
                if api_vct.ndim != 1:
                    raise ValueError("Invalid array dimension: ",api_vct.ndim)
                else:
                    meanstd.ensdam_meanstd.update_meancov_vector_weight(api_vct,varref,weight,weightsum,api_mean,meanref,api_mproda)

class interptools:
    """
    The purpose of InterpTools is to provide tools
    to localize data in grids and compute interpolation coefficients.

    Available variables:
    None

    Available functions:
    grid1d_locate(),grid1d_interp(),grid2d_init(),grid2d_locate(),grid2d_interp()
    """

    @staticmethod
    def grid1d_locate(kgrid,kx):
        """
        located,ki = grid1d_locate(kgrid,kx)

        Locate data point in 1D grid

        Parameters
        ----------
        kgrid : input rank-1 array('d') with bounds (f2py_kgrid_d0)
        kx : input float

        Returns
        -------
        located : int
        ki : int
        """
        api_grid = _iarraydouble(kgrid)
        located,ki = interp.ensdam_interp.grid1d_locate(api_grid,kx)
        return located,ki

    @staticmethod
    def grid1d_interp(kgrid,kx,ki):
        """
        w = grid1d_interp(kgrid,kx,ki)

        Compute interpolation weight in 1D grid

        Parameters
        ----------
        kgrid : input rank-1 array('d') with bounds (f2py_kgrid_d0)
        kx : input float
        ki : input int

        Returns
        -------
        w : float
        """
        api_grid = _iarraydouble(kgrid)
        w = interp.ensdam_interp.grid1d_interp(api_grid,kx,ki)
        return w

    @staticmethod
    def grid2d_init(kxg,kyg,gtype):
        """
        grid2d_init(kxg,kyg,gtype)

        Initialize 2D grid

        Parameters
        ----------
        kxg : input rank-2 array('d') with bounds (f2py_kxg_d0,f2py_kxg_d1)
        kyg : input rank-2 array('d') with bounds (size(kxg,1),size(kxg,2))
        gtype : input string(len=-1)
        """
        api_xg = _iarraydouble(kxg)
        api_yg = _iarraydouble(kyg)
        interp.ensdam_interp.grid2d_init(api_xg,api_yg,gtype)

    @staticmethod
    def grid2d_locate(kx,ky):
        """
        located,ki,kj = grid2d_locate(kx,ky)

        Locate data point in 2D grid

        Parameters
        ----------
        kx : input float
        ky : input float

        Returns
        -------
        grid2d_locate : int
        ki : int
        kj : int
        """
        located,ki,kj = interp.ensdam_interp.grid2d_locate(kx,ky)
        return located,ki,kj
      

    @staticmethod
    def grid2d_interp(kx,ky,ki,kj):
        """
        w = grid2d_interp(kx,ky,ki,kj)

        Compute interpolation weight in 2D grid

        Parameters
        ----------
        kx : input float
        ky : input float
        ki : input int
        kj : input int

        Returns
        -------
        kw : rank-2 array('d') with bounds (2,2)
        """
        w = interp.ensdam_interp.grid2d_interp(kx,ky,ki,kj)
        return w

class obserror:
    """
    The purpose of ObsError is to provide tools to deal with observation errors.

    Available variables:
    None

    Available functions:
    obserror_logpdf(),obserror_cdf(),obserror_sample()
    """

    @staticmethod
    def obserror_logpdf(y,x,sigma):
        """
        obserror_logpdf = obserror_logpdf(y,x,sigma)

        Compute logarithm of observation error probability density function

        Parameters
        ----------
        y : input float or rank-1 array('d') 
        x : input float or rank-1 array('d')
        sigma : input float or rank-1 array('d')

        Returns
        -------
        obserror_logpdf : float
        """
        if numpy.isscalar(y):
            logpdf = obserror.endam_obserror.obserror_logpdf_variable(y,x,sigma)
        elif y.ndim == 1:
            api_y = _iarraydouble(y)
            api_x = _iarraydouble(x)
            if numpy.isscalar(sigma):
                logpdf = obserror.endam_obserror.obserror_logpdf_vector_homogeneous(api_y,api_x,sigma)
            elif sigma.ndim == 1:
                api_sigma = _iarraydouble(sigma)
                logpdf = obserror.endam_obserror.obserror_logpdf_vector(api_y,api_x,api_sigma)
            else:
                raise ValueError("Invalid array dimension")
        else:
            raise ValueError("Invalid array dimension")
        return logpdf

    @staticmethod
    def obserror_cdf(y,x,sigma):
        """
        obserror_cdf = obserror_cdf(y,x,sigma)

        Compute cumulate distribution function for observation errors

        Parameters
        ----------
        y : input float or rank-1 array('d')
        x : input float or rank-1 array('d')
        sigma : input float or rank-1 array('d')

        Returns
        -------
        obserror_cdf : float or rank-1 array('d')
        """
        if numpy.isscalar(y):
            cdf = obserror.endam_obserror.obserror_cdf_variable(y,x,sigma)
        elif y.ndim == 1:
            api_y = _iarraydouble(y)
            api_x = _iarraydouble(x)
            if numpy.isscalar(sigma):
                cdf = obserror.endam_obserror.obserror_cdf_vector_homogeneous(api_y,api_x,sigma)
            elif sigma.ndim == 1:
                api_sigma = _iarraydouble(sigma)
                cdf = obserror.endam_obserror.obserror_cdf_vector(api_y,api_x,api_sigma)
            else:
                raise ValueError("Invalid array dimension")
        else:
            raise ValueError("Invalid array dimension")
        return cdf

    @staticmethod
    def obserror_sample(x,sigma,rank,reuse_last_rank=False):
        """
        obserror_sample = obserror_sample(x,sigma,rank,[reuse_last_rank])

        Sample probability distribution of observation errors

        Parameters
        ----------
        x : input float or rank-1 array('d')
        sigma : input float or rank-1 array('d')
        rank : input float or rank-1 array('i')

        Other Parameters
        ----------------
        reuse_last_rank : input int

        Returns
        -------
        obserror_sample : float
        """
        if numpy.isscalar(x):
            if reuse_last_rank==():
                reuse_last_rank = False
            sample = obserror.endam_obserror.obserror_sample_variable(x,sigma,rank,reuse_last_rank)
        elif y.ndim == 1:
            api_x = _iarraydouble(x)
            if numpy.isscalar(sigma):
                sample = obserror.endam_obserror.obserror_sample_vector_homogeneous(api_x,sigma,rank,reuse_last_rank)
            elif sigma.ndim == 1:
                api_sigma = _iarraydouble(sigma)
                api_rank = _iarraydouble(rank)
                sample = obserror.endam_obserror.obserror_sample_vector(api_x,api_sigma,api_rank,reuse_last_rank)
            else:
                raise ValueError("Invalid array dimension")
        else:
            raise ValueError("Invalid array dimension")
        return sample

class stochtools:
    """
    The purpose of StochTools is to provide tools
    to generate random numbers, random fields,
    or stochastic processes of various types,
    to be used in stochastic modelling or data assimilation systems.

    Available variables:
    nominal_accuracy,accuracy,maxiter,
    storfg_ylm_resolution,mpi_comm_storfg

    Available functions:
    kiss(),kiss_seed(),kiss_reset(),kiss_check(),kiss_save(),kiss_load(),
    kiss_uniform(),kiss_gaussian(),kiss_gamma(),kiss_beta(),kiss_sample()
    ranv_tg(),ran_tg(),ran_te(),
    cdf_gaussian(),pdf_gaussian(),logpdf_gaussian(),invcdf_gaussian(),
    cdf_gamma(),pdf_gamma(),logpdf_gamma(),invcdf_gamma(),
    cdf_beta(),pdf_beta(),logpdf_beta(),invcdf_beta().
    gprod_to_gau(),gau_to_gam(),gau_to_beta(),gen_field_2s()
    """

    @staticmethod
    def kiss():
        """
        kiss = kiss()

        64-bit KISS random number generator (period ~ 2^250)

        Returns
        -------
        kiss : long
        """
        v = storng.ensdam_storng.kiss()
        return v

    @staticmethod
    def kiss_seed(ix,iy,iz,iw):
        """
        kiss_seed(ix,iy,iz,iw)

        Define seeds for KISS random number generator

        Parameters
        ----------
        ix : input long
        iy : input long
        iz : input long
        iw : input long
        """
        storng.ensdam_storng.kiss_seed(ix,iy,iz,iw)

    @staticmethod
    def kiss_reset():
        """
        kiss_reset()

        Reset the default seeds
        """
        storng.ensdam_storng.kiss_reset()

    @staticmethod
    def kiss_check(check_type):
        """
        kiss_check(check_type)

        Check the KISS pseudo-random sequence

        check_type : input string(len=-1)
        """
        storng.ensdam_storng.kiss_check(check_type)

    @staticmethod
    def kiss_save():
        """
        kiss_save()

        Save current state of KISS (for future restart)
        """
        storng.ensdam_storng.kiss_save()

    @staticmethod
    def kiss_load():
        """
        kiss_load()

        Load the saved state of KISS
        """
        storng.ensdam_storng.kiss_load()

    @staticmethod
    def kiss_uniform():
        """
        kiss_uniform = kiss_uniform()

        Real random numbers with uniform distribution in [0,1]

        Returns
        -------
        kiss_uniform : float
        """
        v = storng.ensdam_storng.kiss_uniform()
        return v

    @staticmethod
    def kiss_gaussian():
        """
        kiss_gaussian = kiss_gaussian()

        Real random numbers with Gaussian distribution N(0,1)

        Returns
        -------
        kiss_gaussian : float
        """
        v = storng.ensdam_storng.kiss_gaussian()
        return v

    @staticmethod
    def kiss_gamma(k):
        """
        kiss_gamma = kiss_gamma(k)

        Real random numbers with Gamma distribution Gamma(k,1)

        Parameters
        ----------
        k : input float

        Returns
        -------
        kiss_gamma : float
        """
        v = storng.ensdam_storng.kiss_gamma(k)
        return v

    @staticmethod
    def kiss_beta(a,b):
        """
        kiss_beta = kiss_beta(a,b)

        Real random numbers with Beta distribution Beta(a,b)

        Parameters
        ----------
        a : input float
        b : input float

        Returns
        -------
        kiss_beta : float
        """
        v = storng.ensdam_storng.kiss_beta(a,b)
        return v

    @staticmethod
    def kiss_sample(a,k):
        """
        kiss_sample=kiss_sample(a,k)

        Select a random sample from a set of integers

        Parameters
        ----------
        a : input rank-1 array('i'): integers from which to sample
        k : input int: size of the output sample

        Returns
        -------
        kiss_sample : rank-1 array('i') with bounds (k)

        """
        api_a = _iarrayint(a)
        n = api_a.size
        storng.ensdam_storng.kiss_sample(api_a,n,k)
        return api_a[0:k]

    @staticmethod
    def ranv_tg(n,matarm,vecbm):
        """
        ranv_tg = ranv_tg(tgvsmpl,matarm,vecbm)

        Multidimensional truncated Gaussian random vector generator.
        Sample p(x)  ~  exp (- xT x / 2 ) , A x <= b.
        Gibbs sampler, using one-dimensional truncated Gaussian random number generator.

        Parameters
        ----------
        n : size of the output sample
        matarm : input rank-2 array('d') with bounds (r,m)
        vecbm : input rank-1 array('d') with bounds (m))

        Returns
        -------
        ran_tg : rank-2 array('d') with bounds (n,r)
        """
        print "Warning: python interface untested !!"
        api_matarm = _iarraydouble(matarm)
        api_vecbm = _iarraydouble(vecbm)
        r = matarm.shape[0] # number of dimensions
        api_smp = _iarraydouble(numpy.zeros([n,r]))
        stotge.ensdam_stotge.ran_tg(api_smp,api_matarm,api_vecbm)
        return api_smp

    @staticmethod
    def ran_tg(n,aa,bb):
        """
        ran_tg = ran_tg(aa,bb)

        One-dimensional truncated Gaussian random number generator.
        Sample p(y)  ~  exp(-y*y) , aa < y < bb.
        The algorithm is taken from Geweke, 1991.

        Parameters
        ----------
        n : size of the output sample
        aa : input float
        bb : input float

        Returns
        -------
        ran_tg : rank-1 array('d') with bounds (n)
        """
        api_smp = _iarraydouble(numpy.zeros([n]))
        stotge.ensdam_stotge.ran_tg(api_smp,aa,bb)
        return api_smp

    @staticmethod
    def ran_te(a):
        """
        ran_te = ran_te(a)

        Truncated exponential random number generator.
        Sample: p(y) = a * exp(a*a) * exp(-ay) , y > a

        Parameters
        ----------
        a : input float

        Returns
        -------
        ran_te : float
        """
        v = stotge.ensdam_stotge.ran_te(a)
        return v

    @staticmethod
    def cdf_gaussian(x):
        """
        cdf_gaussian = cdf_gaussian(x)

        Evaluate Gaussian cumulative distribution function (cdf)

        Parameters
        ----------
        x : input float

        Returns
        -------
        cdf_gaussian : float
        """
        stoutil.ensdam_stoutil.nominal_accuracy = nominal_accuracy
        stoutil.ensdam_stoutil.accuracy = accuracy
        stoutil.ensdam_stoutil.maxiter = maxiter
        v = stoutil.ensdam_stoutil.cdf_gaussian(x)
        return v

    @staticmethod
    def pdf_gaussian(x):
        """
        pdf_gaussian = pdf_gaussian(x)

        Evaluate Gaussian probability density function (pdf)

        Parameters
        ----------
        x : input float

        Returns
        -------
        pdf_gaussian : float
        """
        stoutil.ensdam_stoutil.nominal_accuracy = nominal_accuracy
        stoutil.ensdam_stoutil.accuracy = accuracy
        stoutil.ensdam_stoutil.maxiter = maxiter
        v = stoutil.ensdam_stoutil.pdf_gaussian(x)
        return v

    @staticmethod
    def logpdf_gaussian(x):
        """
        logpdf_gaussian = logpdf_gaussian(x)

        Evaluate log of Gaussian probability density function (pdf)

        Parameters
        ----------
        x : input float

        Returns
        -------
        logpdf_gaussian : float
        """
        stoutil.ensdam_stoutil.nominal_accuracy = nominal_accuracy
        stoutil.ensdam_stoutil.accuracy = accuracy
        stoutil.ensdam_stoutil.maxiter = maxiter
        v = stoutil.ensdam_stoutil.logpdf_gaussian(x)
        return v

    @staticmethod
    def invcdf_gaussian(x):
        """
        invcdf_gaussian = invcdf_gaussian(x)

        Evaluate inverse Gaussian cumulative distribution function (cdf)

        Parameters
        ----------
        x : input float

        Returns
        -------
        invcdf_gaussian : float
        """
        stoutil.ensdam_stoutil.nominal_accuracy = nominal_accuracy
        stoutil.ensdam_stoutil.accuracy = accuracy
        stoutil.ensdam_stoutil.maxiter = maxiter
        v = stoutil.ensdam_stoutil.invcdf_gaussian(x)
        return v

    @staticmethod
    def cdf_gamma(a,x):
        """
        cdf_gamma = cdf_gamma(a,x)

        Evaluate Gamma cumulative distribution function (cdf)

        Parameters
        ----------
        a : input float
        x : input float

        Returns
        -------
        cdf_gamma : float
        """
        stoutil.ensdam_stoutil.nominal_accuracy = nominal_accuracy
        stoutil.ensdam_stoutil.accuracy = accuracy
        stoutil.ensdam_stoutil.maxiter = maxiter
        v = stoutil.ensdam_stoutil.cdf_gamma(a,x)
        return v

    @staticmethod
    def pdf_gamma(a,x):
        """
        pdf_gamma = pdf_gamma(a,x)

        Evaluate Gamma probability density function (pdf)

        Parameters
        ----------
        a : input float
        x : input float

        Returns
        -------
        pdf_gamma : float
        """
        stoutil.ensdam_stoutil.nominal_accuracy = nominal_accuracy
        stoutil.ensdam_stoutil.accuracy = accuracy
        stoutil.ensdam_stoutil.maxiter = maxiter
        v = stoutil.ensdam_stoutil.pdf_gamma(a,x)
        return v

    @staticmethod
    def logpdf_gamma(a,x):
        """
        logpdf_gamma = logpdf_gamma(a,x)

        Evaluate log of Gamma probability density function (pdf)

        Parameters
        ----------
        a : input float
        x : input float

        Returns
        -------
        logpdf_gamma : float
        """
        stoutil.ensdam_stoutil.nominal_accuracy = nominal_accuracy
        stoutil.ensdam_stoutil.accuracy = accuracy
        stoutil.ensdam_stoutil.maxiter = maxiter
        v = stoutil.ensdam_stoutil.logpdf_gamma(a,x)
        return v

    @staticmethod
    def invcdf_gamma(a,x):
        """
        invcdf_gamma = invcdf_gamma(a,x)

        Evaluate inverse Gamma cumulative distribution function (cdf)

        Parameters
        ----------
        a : input float
        x : input float

        Returns
        -------
        invcdf_gamma : float
        """
        stoutil.ensdam_stoutil.nominal_accuracy = nominal_accuracy
        stoutil.ensdam_stoutil.accuracy = accuracy
        stoutil.ensdam_stoutil.maxiter = maxiter
        v = stoutil.ensdam_stoutil.invcdf_gamma(a,x)
        return v

    @staticmethod
    def cdf_beta(a,b,x):
        """
        cdf_beta = cdf_beta(a,b,x)

        Evaluate Beta cumulative distribution function (cdf)

        Parameters
        ----------
        a : input float
        b : input float
        x : input float

        Returns
        -------
        cdf_beta : float
        """
        stoutil.ensdam_stoutil.nominal_accuracy = nominal_accuracy
        stoutil.ensdam_stoutil.accuracy = accuracy
        stoutil.ensdam_stoutil.maxiter = maxiter
        v = stoutil.ensdam_stoutil.cdf_beta(a,b,x)
        return v

    @staticmethod
    def pdf_beta(a,b,x):
        """
        pdf_beta = pdf_beta(a,b,x)

        Evaluate Beta probability density function (pdf)

        Parameters
        ----------
        a : input float
        b : input float
        x : input float

        Returns
        -------
        pdf_beta : float
        """
        stoutil.ensdam_stoutil.nominal_accuracy = nominal_accuracy
        stoutil.ensdam_stoutil.accuracy = accuracy
        stoutil.ensdam_stoutil.maxiter = maxiter
        v = stoutil.ensdam_stoutil.pdf_beta(a,b,x)
        return v

    @staticmethod
    def logpdf_beta(a,b,x):
        """
        logpdf_beta = logpdf_beta(a,b,x)

        Evaluate log of Beta probability density function (pdf)

        Parameters
        ----------
        a : input float
        b : input float
        x : input float

        Returns
        -------
        logpdf_beta : float
        """
        stoutil.ensdam_stoutil.nominal_accuracy = nominal_accuracy
        stoutil.ensdam_stoutil.accuracy = accuracy
        stoutil.ensdam_stoutil.maxiter = maxiter
        v = stoutil.ensdam_stoutil.logpdf_beta(a,b,x)
        return v

    @staticmethod
    def invcdf_beta(a,b,x):
        """
        invcdf_beta = invcdf_beta(a,b,x)

        Evaluate inverse Beta cumulative distribution function (cdf)

        Parameters
        ----------
        a : input float
        b : input float
        x : input float

        Returns
        -------
        invcdf_beta : float
        """
        stoutil.ensdam_stoutil.nominal_accuracy = nominal_accuracy
        stoutil.ensdam_stoutil.accuracy = accuracy
        stoutil.ensdam_stoutil.maxiter = maxiter
        v = stoutil.ensdam_stoutil.invcdf_beta(a,b,x)
        return v

    @staticmethod
    def gprod_to_gau(x):
        """
        y = gprod_to_gau(x)

        Transform a Gaussian product random number
        (product of two N(0,1) numbers)
        into a Gaussian random number, with N(0,1) distribution

        Parameters
        ----------
        x : input float

        Returns
        -------
        y : float
        """
        v = stoanam.ensdam_stoanam.gprod_to_gau(x)
        return v

    @staticmethod
    def gau_to_gam(x,k):
        """
        y = gau_to_gam(x,k)

        Transform a Gaussian random number (with N(0,1) pdf)
        into a Gamma random number, with GAMMA(k,1)

        The shape parameter k is assumed larger than 1.
        The scale parameter theta is set to 1
        (it can be subsequently scaled to any value).

        Parameters
        ----------
        x : input float
        k : input float

        Returns
        -------
        y : float
        """
        v = stoanam.ensdam_stoanam.gau_to_gam(x,k)
        return v

    @staticmethod
    def gau_to_beta(x,mu,stdmax):
        """
        y = gau_to_beta(x,mu,stdmax)

        Transform a Gaussian random number (with N(0,1) pdf)
        into a Beta random number, with BETA(a,b)
        with a = mu * nu, and b = (1-mu) * nu

        The mean mu must be between 0 and 1 (strictly).
        The sample size nu must be positive (strictly).
        The variance is equal  to mu * (1-mu) / (1+nu)

        The sample size nu is computed from the maximum std
        (stdmax, occuring if mu=1/2), from the formula
        nu = 1 / ( 4 * stdmax **2 )  -  1
        This maximum standard deviation must be positive
        and smaller than 1/2 (strictly).

        Parameters
        ----------
        x : input float
        mu : input float
        stdmax : input float

        Returns
        -------
        y : float
        """
        v = stoanam.ensdam_stoanam.gau_to_beta(x,mu,stdmax)
        return v

    @staticmethod
    def gen_field_2s(lon,lat,pow_spect,lmin,lmax,pow_spect_extra_args=()):
        """
        ranfield = gen_field_2s_new(lon,lat,pow_spect,lmin,lmax,[pow_spect_extra_args])

        Parameters
        ----------
        lon : input rank-1 or rank-2 array('d')
        lat : input rank-1 or rank-2 array('d')
        pow_spect : call-back function => fun_pow_spect_sph
        lmin : input int
        lmax : input int

        Other Parameters
        ----------------
        pow_spect_extra_args : input tuple, optional
            Default: ()

        Returns
        -------
        ranfield : rank-1 or rank-2  array('d')

        Notes
        -----
        Call-back functions::

          def fun_pow_spect_sph(l,m): return fun_pow_spect_sph
          Required arguments:
            l : input int
            m : input int
          Return objects:
            fun_pow_spect_sph : float
        """
        print "Warning: python interface untested !!"
        api_lon = _iarraydouble(lon)
        api_lat = _iarraydouble(lat)
        if lon.ndim == 1:
            if lat.ndim != 1:
                raise ValueError("Inconsistent array dimnension: ",lat.ndim)
            api_ranfield = gen_field_2s_new1(api_lon,api_lat,pow_spect,lmin,lmax,[pow_spect_extra_args])
        elif lon.ndim == 2:
            if lat.ndim != 2:
                raise ValueError("Inconsistent array dimnension: ",lat.ndim)
            api_ranfield = gen_field_2s_new2(api_lon,api_lat,pow_spect,lmin,lmax,[pow_spect_extra_args])
        else:
            raise ValueError("Invalid array dimnesion: ",lon.ndim)
        return api_ranfield

class transpho:
    """
    The purpose of TranSpHO is to provide tools
    to transform the data assimilation problem
    by projection on the spherical harmonics.
    The objective is to separate scales and
    to make the assimilation system behave differently for different scales
    (e.g.\ different localization, different error parameterization,..).

    Available variables:
    regr_epsilon,regr_overlap,regr_maxiter,regr_type,regr_rho,regr_maxbloc,
    mpi_comm_sphylm,external_vector_decomposition

    Available functions:
    init_ylm(),init_regr_ylm(),proj_ylm(),back_ylm(),back_ylm_loc(),
    regr_ylm(),disp_ylm(),ylm(),mesh_area()
    """

    @staticmethod
    def init_ylm(kjpl,kjlmin,latmin,latmax,dlatmax):
        """
        init_ylm(kjpl,kjlmin,latmin,latmax,dlatmax)

        Initialize computation of spherical harmonics

        Parameters
        ----------
        kjpl : input int
        kjlmin : input int
        latmin : input float
        latmax : input float
        dlatmax : input float
        """
        ensdam.transpho.init_ylm(kjpl,kjlmin,latmin,latmax,dlatmax)

    @staticmethod
    def init_regr_ylm(ktype,kmaxiter,kmaxbloc,koverlap,kepsilon,krho):
        """
        init_regr_ylm(ktype,kmaxiter,kmaxbloc,koverlap,kepsilon,krho)

        Initialize regression of observations

        Parameters
        ----------
        ktype : input string(len=-1)
        kmaxiter : input int
        kmaxbloc : input int
        koverlap : input int
        kepsilon : input float
        krho : input float
        """

    @staticmethod
    def proj_ylm(ktab,klon,klat):
        """
        kproj = proj_ylm(ktab,klon,klat)

        Project on spherical harmonics

        Parameters
        ----------
        ktab : input rank-1 array('d') with bounds (f2py_ktab_d0)
        klon : input rank-1 array('d') with bounds (f2py_klon_d0)
        klat : input rank-1 array('d') with bounds (f2py_klat_d0)

        Returns
        -------
        kproj : rank-2 array('d') with bounds (f2py_kproj_d0,f2py_kproj_d1)
        """
        # Update module public variables
        sphylm.ensdam_sphylm.external_vector_decomposition = external_vector_decomposition

    @staticmethod
    def back_ylm(kproj,klon,klat):
        """
        ktab = back_ylm(kproj,klon,klat)

        Transform back on the sphere

        Parameters
        ----------
        kproj : input rank-2 array('d') with bounds (f2py_kproj_d0,f2py_kproj_d1)
        klon : input rank-1 array('d') with bounds (f2py_klon_d0)
        klat : input rank-1 array('d') with bounds (f2py_klat_d0)

        Returns
        -------
        ktab : rank-1 array('d') with bounds (f2py_ktab_d0)
        """
        # Update module public variables
        sphylm.ensdam_sphylm.external_vector_decomposition = external_vector_decomposition

    @staticmethod
    def back_ylm_loc(kproj,klon,klat,kl0,kl1):
        """
        ktab = back_ylm_loc(kproj,klon,klat,kl0,kl1)

        Transform back on the sphere (for a range of degrees)

        Parameters
        ----------
        kproj : input rank-2 array('d') with bounds (f2py_kproj_d0,f2py_kproj_d1)
        klon : input rank-1 array('d') with bounds (f2py_klon_d0)
        klat : input rank-1 array('d') with bounds (f2py_klat_d0)
        kl0 : input int
        kl1 : input int

        Returns
        -------
        ktab : rank-1 array('d') with bounds (f2py_ktab_d0)
        """
        # Update module public variables
        sphylm.ensdam_sphylm.external_vector_decomposition = external_vector_decomposition

    @staticmethod
    def regr_ylm(kwei,kobs,klon,klat,kobswei):
        """
        kregr = regr_ylm(kwei,kobs,klon,klat,kobswei)

        Regression of observations on spherical harmonics

        Parameters
        ----------
        kwei : input rank-2 array('d') with bounds (f2py_kwei_d0,f2py_kwei_d1)
        kobs : in/output rank-1 array('d') with bounds (f2py_kobs_d0)
        klon : input rank-1 array('d') with bounds (f2py_klon_d0)
        klat : input rank-1 array('d') with bounds (f2py_klat_d0)
        kobswei : input rank-1 array('d') with bounds (f2py_kobswei_d0)

        Returns
        -------
        kregr : rank-2 array('d') with bounds (f2py_kregr_d0,f2py_kregr_d1)
        """
        # Update module public variables
        sphylm.ensdam_sphylm.external_vector_decomposition = external_vector_decomposition
        sphylm.ensdam_sphylm.regr_epsilon = regr_epsilon
        sphylm.ensdam_sphylm.regr_overlap = regr_overlap
        sphylm.ensdam_sphylm.regr_maxiter = regr_maxiter
        sphylm.ensdam_sphylm.regr_type = regr_type
        sphylm.ensdam_sphylm.regr_rho = regr_rho
        sphylm.ensdam_sphylm.regr_maxbloc = regr_maxbloc

    @staticmethod
    def disp_ylm(klon,klat,kl,km):
        """
        ktab = disp_ylm(klon,klat,kl,km)

        Output one single spherical harmonics

        Parameters
        ----------
        klon : input rank-1 array('d') with bounds (f2py_klon_d0)
        klat : input rank-1 array('d') with bounds (f2py_klat_d0)
        kl : input int
        km : input int

        Returns
        -------
        ktab : rank-1 array('d') with bounds (f2py_ktab_d0)
        """
        # Update module public variables
        sphylm.ensdam_sphylm.external_vector_decomposition = external_vector_decomposition


    @staticmethod
    def ylm(kl,km,klon,klat):
        """
        ylm = ylm(kl,km,klon,klat)

        Evaluate spherical harmonics

        Parameters
        ----------
        kl : input int
        km : input int
        klon : input float
        klat : input float

        Returns
        -------
        ylm : float
        """
        v = sphylm.ensdam_sphylm.ylm(kl,km,klon,klat)
        return v

    @staticmethod
    def mesh_area(lon,lat):
        """
        area = mesh_area(lon,lat)

        Parameters
        ----------
        lon : input rank-2 array('d') with bounds (f2py_lon_d0,f2py_lon_d1)
        lat : input rank-2 array('d') with bounds (size(lon,1),size(lon,2))

        Returns
        -------
        area : rank-2 array('d') with bounds (size(lon,1)-1,size(lon,2)-2)
        """

