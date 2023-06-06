# cython: language_level=3
# cython: profile=True
"""
pyensdam.interpolation: interpolation tools
===========================================

Available functions:


Module parameters:


"""

cimport cython
import tempfile
import atexit
import os

cimport numpy
import numpy

# Definition of external C-callable routines
cdef extern void c_grid2D_init(int ni,int nj,double* xgrid, double* ygrid,int* gtype) nogil
cdef extern void c_grid1D_locate(int ngrid,double* grid,double* x,int* ix,int* located) nogil
cdef extern void c_grid2D_locate(double* x,double* y,int* ix,int* iy,int* located) nogil
cdef extern void c_grid1D_interp(int ngrid,double* grid,double* x,int* ix,double* w) nogil
cdef extern void c_grid2D_interp(double* x,double* y,int* ix,int* iy,double* w) nogil

# Public function to locate position in 1D grid and compute interpolation weights
def locate1D(double[::1] grid not None,x):
    """location, weight = locate1D(grid,x)

       Locate positions in 1D grid and compute interpolation weights

       Inputs
       ------
       grid  [rank-1 double array] : definition of the grid locations (in ascending order (ngrid)
       x [double array] : array of positions to locate in the grid

       Returns
       -------
       location [integer array] : index of the grid cell where postions are located (same shape as x)
       weight [double array] : interpolation weight to use (same shape as x)

    """
    location = numpy.zeros(x.shape), dtype=numpy.intc)
    weight = numpy.zeros(x.shape), dtype=numpy.double)
    cdef int location_, located_
    cdef double x_, weight_
    for i in range(numpy.prod(x.shape)):
      indices = numpy.unravel_index(i,x.shape)

      x_ = x[indices]
      c_grid1D_locate(<int>grid.shape[0],&grid[0],&x_,&location_,&located_)
      location[indices] = location_
      if located_ == 1:
        c_grid1D_interp(<int>grid.shape[0],&grid[0],&x_,&location_,&weight_)
        weight[indices] = weight_
      else:
        location[indices] = -1
        weight[indices] = None

    return location, weight

# Public function to apply interpolation on 1D input field
def interp1D(field,location,weight):
    """field_interpolated = interp1D(field,location,weight)

       Apply interpolation on input field

       Inputs
       ------
       field [rank-1 double array] : field from which interpolating (same dimensions are as grid)
       location [integer array] : index of the grid cell where interpolated postions are located
       weight [double array] : interpolation weight to use (same shape as location)

       Returns
       -------
       field_interpolated [double array] : interpolation weight to use (same shape as location)

    """
    for i in range(numpy.prod(location.shape)):
      indices = numpy.unravel_index(i,location.shape)

      loc0 = location[indices]-1
      loc1 = location[indices]
      w1 = weight[indices]
      w2 = 1 - weight1

      if location[indices] == -1 :
        field_interpolated[indices] = None
      else:
        field_interpolated[indices] = w1 * field[loc0] + w2 * field[loc1]

    return field_interpolated

# Public function to locate position in 1D grid and compute interpolation weights
def locate2D(x,y):
    """location, weight = locate2D(x,y)

       Locate positions in 2D grid and compute interpolation weights

       Inputs
       ------
       x [double array] : array of positions to locate in the grid (x-coordinate)
       y [double array] : array of positions to locate in the grid (y-coordinate)

       Returns
       -------
       location [integer array] : index of the grid cell where postions are located
                                  (same shape as x and y, for a 2-dimension vector)
       weight [double array] : interpolation weight to use
                               (same shape as x and y, for a '2 by 2' matrix)

    """
    location = numpy.zeros((x.shape,2)), dtype=numpy.intc)
    weight = numpy.zeros((x.shape,2,2)), dtype=numpy.double)
    cdef int locx_, locy_, located_
    cdef double[2,2] weight_
    for i in range(numpy.prod(x.shape)):
      indices = numpy.unravel_index(i,x.shape)

      x_ = x[indices]
      y_ = y[indices]
      c_grid2D_locate(&x_,&y_,&locx_,&locy_,&located_)
      location[indices,0] = locx_
      location[indices,1] = locy_
      if located_ == 1:
        c_grid2D_interp(&x_,&y_,&locx_,&locy_,&weight_)
        weight[indices,:,:] = weight_[:,:]
      else:
        location[indices,0] = -1
        location[indices,1] = -1
        weight[indices,:,:] = None

    return location, weight

# Public function to apply interpolation on 2D input field
def interp2D(field,location,weight):
    """field_interpolated = interp2D(field,location,weight)

       Apply interpolation on input field

       Inputs
       ------
       field [rank-2 double array] : field from which interpolating (same dimensions are as grid)
       location [integer array] : index of the grid cell where interpolated postions are located
       weight [double array] : interpolation weight to use (same shape as location)

       Returns
       -------
       field_interpolated [double array] : interpolation weight to use (same shape as location)

    """
    for i in range(numpy.prod(location.shape)):
      indices = numpy.unravel_index(i,location.shape)

      locx0 = location[indices,0]-1
      locx1 = location[indices,0]
      locy0 = location[indices,1]-1
      locy1 = location[indices,1]
      w1 = weight[indices,0,0]
      w2 = weight[indices,0,1]
      w3 = weight[indices,1,0]
      w4 = weight[indices,1,1]

      if location[indices,0] == -1 :
        field_interpolated[indices] = None
      else:
        field_interpolated[indices] = w1 * field[locx0,locy0] + \
                                      w2 * field[locx0,locy1] + \
                                      w3 * field[locx1,locy0] + \
                                      w4 * field[locx1,locy1]

    return field_interpolated

