"""
pyensdam: Ensemble data assimilation modules
============================================

Avalaible modules:
 - pyensdam.anamorphosis : ensemble anamorphosis transformation
 - pyensdam.augment : ensemble augmentation with MCMC sampler
 - pyensdam.interpolation : interpolation in 1D and 2D grids
 - pyensdam.obserror : observation error (normal, lognormal, gamma, beta)
 - pyensdam.probability : probability distribution (normal, lognormal, gamma, beta)
 - pyensdam.random : random field generator
 - pyensdam.scores : ensemble probabilistic scores
 - pyensdam.statistics : ensemble statistics (mean, std, correlation, representer)
 - pyensdam.transpho : transformation in the basis of the spherical harmonics
 - pyensdam.update : ensemble observational update, with an MCMC sampler

"""

__version__ = '0.1.2'

from . import anamorphosis
from . import augment
from . import scores
from . import statistics
from . import update
from . import interpolation
from . import obserror
from . import random
from . import probability
from . import transpho

