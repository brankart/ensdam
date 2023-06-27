# cython: language_level=3
# cython: profile=True
"""
pyensdam.transpho: transformation in the basis of the spherical harmonics
=========================================================================

Available functions:
 - transpho.forward : forward transformation (compute spectrum)
 - transpho.backward : backward transformation
 - transpho.mesh_area : compute area of grid cells on the sphere

Module parameters:

"""

cimport cython
import tempfile
import atexit
import os

cimport numpy
import numpy

# Definition of external C-callable routines
cdef extern void c_init_ylm(int nl,int lmin,double* latmin,double* latmax,double* dlatmax) nogil
cdef extern void c_proj_ylm(int nvar,int nl,double* proj,double* tab,double* lon,double* lat) nogil
cdef extern void c_back_ylm(int nvar,int nl,double* proj,double* tab,double* lon,double* lat) nogil
cdef extern void c_back_ylm_loc(int nvar,int nl,double* proj,double* tab,double* lon,double* lat,int* l0,int* l1) nogil
cdef extern void c_mesh_area(int ni,int nj,double* lon,double* lat,double* area) nogil

# Set default values of module attributes (to precompute the Legendre polynomials)
reinitialize=True  # Reinitilize module with new attributes
lmin=0      # Minimum degree of the spherical harmonics
lmax=500    # Maximum degree of the spherical harmonics
latmin=-90. # Minimum latitude
latmax=90.  # Maximum latitude
latres=0.05 # Resolution of Legendre polynomials (in degrees)

# Public function to perform forward transformation (compute spectrum)
def forward(double[::1] field, double[::1] lon, double[::1] lat):
    """spectrum = forward(field,lon,lat)

       Forward transformation (compute spectrum)

       Inputs
       ------
       field [rank-1 double array] : field to transform (nvar)
       lon [rank-1 double array] : longitude of grid points (nvar)
       lat [rank-1 double array] : latitude of grid points (nvar)

       Returns
       -------
       spectrum [rank-2 double array] : spectrum in the basis of the spherical harmonics

    """
    # Reinitialize the precomputation of Legendre polynomials if required
    cdef int lmin_=lmin, lmax_=lmax
    cdef double latmin_=latmin, latmax_=latmax, latres_=latres
    global reinitialize
    if (reinitialize):
      c_init_ylm(lmax_,lmin_,&latmin_,&latmax_,&latres_)
      reinitialize=False

    # Perform forward transformation
    spectrum = numpy.zeros((2*lmax+1,lmax+1), dtype=numpy.double)
    cdef double[:,::1] spectrum_ = spectrum
    c_proj_ylm(<int>lon.shape[0],lmax_,&spectrum_[0,0],&field[0],&lon[0],&lat[0])

    return spectrum

# Public function to perform forward transformation (compute spectrum)
def backward(double[:,::1] spectrum, double[::1] lon, double[::1] lat, l0=None, l1=None):
    """field = forward(spectrum,lon,lat,[l0,l1])

       Backward transformation (from spectrum)

       Inputs
       ------
       spectrum [rank-2 double array] : spectrum in the basis of the spherical harmonics
       lon [rank-1 double array] : longitude of grid points (nvar)
       lat [rank-1 double array] : latitude of grid points (nvar)
       l0 [integer] : minimum degree to use in the backward transformation
       l1 [integer] : maximum degree to use in the backward transformation

       Returns
       -------
       field [rank-1 double array] : reconstructed field (nvar)

    """
    # Reinitialize the precomputation of Legendre polynomials if required
    cdef int lmin_=lmin, lmax_=lmax
    cdef double latmin_=latmin, latmax_=latmax, latres_=latres
    global reinitialize
    if (reinitialize):
      c_init_ylm(lmax_,lmin_,&latmin_,&latmax_,&latres_)
      reinitialize=False

    # Test size of input spectrum array
    if (spectrum.shape[0] != 2*lmax+1):
      raise ValueError("Invalid spectrum array dimensions in transpho")
    if (spectrum.shape[1] != lmax+1):
      raise ValueError("Invalid spectrum array dimensions in transpho")

    # Perform backward transformation
    cdef int l0_, l1_
    field = numpy.zeros_like(lon)
    cdef double[::1] field_ = field

    if ((l0==None) & (l1==None)):
      c_back_ylm(<int>lon.shape[0],lmax_,&spectrum[0,0],&field_[0],&lon[0],&lat[0])
    else:
      if (l0==None):
        l0_ = 0 ; l1_ = l1
      elif (l1==None):
        l0_ = l0 ; l1_ = lmax
      else:
        l0_ = l0 ; l1_ = l1

      c_back_ylm_loc(<int>lon.shape[0],lmax_,&spectrum[0,0],&field_[0],&lon[0],&lat[0],&l0_,&l1_)

    return field

# Public function to compute the area of grid cells on the sphere
def mesh_area(double[:,::1] lon,double[:,::1] lat):
    """area = mesh_area(lon,lat)

       Compute the area of grid cells on the sphere

       Inputs
       ------
       lon [rank-2 double array] : longitude of grid points (nx,ny)
       lat [rank-2 double array] : latitude of grid points (nx,ny)

       Returns
       -------
       area [rank-2 double array] : area of grid cells (nx,ny)

    """
    area = numpy.zeros_like(lon)
    cdef double[:,::1] area_ = area

    c_mesh_area(<int>lon.shape[1],<int>lon.shape[0],&lon[0,0],&lat[0,0],&area_[0,0])

    return area
