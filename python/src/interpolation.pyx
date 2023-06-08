# cython: language_level=3
# cython: profile=True
"""
pyensdam.interpolation: interpolation tools
===========================================

Available functions:
 -  interpolation.locate1D : locate positions in 1D grid and compute interpolation weights
 -  interpolation.interp1D : apply interpolation on 1D input field
 -  interpolation.define2D : define 2D grid
 -  interpolation.locate2D : locate positions in 2D grid and compute interpolation weights
 -  interpolation.interp2D : apply interpolation on 2D input field

"""

cimport cython
import tempfile
import atexit
import os

cimport numpy
import numpy

# Definition of external C-callable routines
cdef extern void c_grid2d_init(int ni,int nj,double* xgrid, double* ygrid,int* gtype) nogil
cdef extern void c_grid1d_locate(int ngrid,double* grid,double* x,int* ix,int* located) nogil
cdef extern void c_grid2d_locate(double* x,double* y,int* ix,int* iy,int* located) nogil
cdef extern void c_grid1d_interp(int ngrid,double* grid,double* x,int* ix,double* w) nogil
cdef extern void c_grid2d_interp(double* x,double* y,int* ix,int* iy,double* w) nogil

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
    location = numpy.zeros(x.shape, dtype=numpy.intc)
    weight = numpy.zeros(x.shape, dtype=numpy.double)
    cdef int location_, located_
    cdef double x_, weight_
    for i in range(numpy.prod(x.shape)):
      indices = numpy.unravel_index(i,x.shape)

      x_ = x[indices]
      c_grid1d_locate(<int>grid.shape[0],&grid[0],&x_,&location_,&located_)
      location[indices] = location_
      if located_ == 1:
        c_grid1d_interp(<int>grid.shape[0],&grid[0],&x_,&location_,&weight_)
        weight[indices] = weight_
      else:
        location[indices] = -1
        weight[indices] = None

    return location, weight

# Public function to apply interpolation on 1D input field
def interp1D(double[::1] field,location,weight):
    """field_interpolated = interp1D(field,location,weight)

       Apply interpolation on 1D input field

       Inputs
       ------
       field [rank-1 double array] : field from which interpolating (same dimensions are as grid)
       location [integer array] : index of the grid cell where interpolated postions are located
       weight [double array] : interpolation weight to use (same shape as location)

       Returns
       -------
       field_interpolated [double array] : interpolation weight to use (same shape as location)

    """
    cdef int loc0,loc1
    cdef double w1, w2, result
    field_interpolated = numpy.zeros_like(weight)

    for i in range(numpy.prod(location.shape)):
      indices = numpy.unravel_index(i,location.shape)

      loc0 = location[indices]-1
      loc1 = location[indices]
      w1 = weight[indices]
      w2 = 1 - w1

      if location[indices] == -1 :
        field_interpolated[indices] = None
      else:
        result = w1 * field[loc0] + w2 * field[loc1]
        field_interpolated[indices] = result

    return field_interpolated

# Public function to define 2D grid
def define2D(double[:,::1] xgrid,double[:,::1] ygrid,grid_type='cartesian'):
    """define2D(xgrid,ygrid,[grid_type])

       Define 2D grid

       Inputs
       ------
       xgrid  [rank-2 double array] : definition of the x-coordinates
       ygrid  [rank-2 double array] : definition of the x-coordinates
       grid_type: typr of coordinates (default='cartesian', spherical)

    """
    cdef int gtype = 0
    if grid_type == 'spherical' :
      gtype = 1

    c_grid2d_init(<int>xgrid.shape[1],<int>xgrid.shape[0],&xgrid[0,0],&ygrid[0,0],&gtype)

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
    locx = numpy.zeros(x.shape, dtype=numpy.intc)
    locy = numpy.zeros(x.shape, dtype=numpy.intc)
    w00 = numpy.zeros_like(x)
    w01 = numpy.zeros_like(x)
    w10 = numpy.zeros_like(x)
    w11 = numpy.zeros_like(x)
    cdef int locx_, locy_, located_
    cdef double x_, y_
    cell_weight = numpy.zeros((2,2), dtype=numpy.double)
    cdef double[:,::1] cell_weight_ = cell_weight
    for i in range(numpy.prod(x.shape)):
      indices = numpy.unravel_index(i,x.shape)

      x_ = x[indices]
      y_ = y[indices]
      c_grid2d_locate(&x_,&y_,&locx_,&locy_,&located_)
      locx[indices] = locx_
      locy[indices] = locy_
      if located_ == 1:
        c_grid2d_interp(&x_,&y_,&locx_,&locy_,&cell_weight_[0,0])
        w00[indices] = cell_weight[0,0]
        w01[indices] = cell_weight[1,0]
        w10[indices] = cell_weight[0,1]
        w11[indices] = cell_weight[1,1]
      else:
        locx[indices] = -1
        locy[indices] = -1

    return [locx,locy], [w00,w01,w10,w11]

# Public function to apply interpolation on 2D input field
def interp2D(double[:,::1] field,location,weight):
    """field_interpolated = interp2D(field,location,weight)

       Apply interpolation on 2D input field

       Inputs
       ------
       field [rank-2 double array] : field from which interpolating (same dimensions are as grid)
       location [integer array] : index of the grid cell where interpolated postions are located
       weight [double array] : interpolation weight to use (same shape as location)

       Returns
       -------
       field_interpolated [double array] : interpolation weight to use (same shape as location)

    """
    cdef int locx0,locx1,locy0,locy1
    cdef double w1, w2, w3, w4, result
    field_interpolated = numpy.zeros_like(weight[0])

    for i in range(numpy.prod(location[0].shape)):
      indices = numpy.unravel_index(i,location[0].shape)

      if location[0][indices] == -1 :
        field_interpolated[indices] = None
      else:
        locx0 = location[0][indices]-1
        locx1 = location[0][indices]
        locy0 = location[1][indices]-1
        locy1 = location[1][indices]
        w00 = weight[0][indices]
        w01 = weight[1][indices]
        w10 = weight[2][indices]
        w11 = weight[3][indices]

        result = w00 * field[locx0,locy0] + w01 * field[locx0,locy1] + \
                 w10 * field[locx1,locy0] + w11 * field[locx1,locy1]
        field_interpolated[indices] = result

    return field_interpolated

