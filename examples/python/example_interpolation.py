import pyensdam as edam
import numpy as np

# Example illustrating the module: interpolation
# ==============================================

# The module contains the functions:
# - locate1D : locate positions in 1D grid and compute interpolation weights
# - interp1D : apply interpolation on 1D input field
# - define2D : define 2D grid
# - locate2D : locate positions in 2D grid and compute interpolation weights
# - interp2D : apply interpolation on 2D input field

grid_size = 3

print('----------------------')
print('1. Interpolation in 1D')
print('----------------------')

# Definition of the original 1D grid
x0 = 0. ; dx = 1. ; nx = grid_size
x = np.arange(x0, x0 + nx * dx, dx, dtype=np.float64)

print('Original grid: x=',x)

# Definition of the 1D field on this grid
f = np.zeros_like(x)
f[::2] = 1  # set f=1 at even indices

print('Original field: f=',f)

# Definition of the target grid
x0 = 0.2 ; dx = 1 ; nx = grid_size - 1
xnew = np.arange(x0, x0 + nx * dx, dx, dtype=np.float64)

print('Target grid: xnew=',xnew)

# Locate target points in original grid
location, weight = edam.interpolation.locate1D(x,xnew)

print('Location indices: location=',location)
print('Interpolation weight: weight=',weight)

# Interpolation original field at new locations
fnew = edam.interpolation.interp1D(f,location,weight)

print('Interpolated field: fnew=',fnew)

print('----------------------')
print('2. Interpolation in 2D')
print('----------------------')

# Definition of the original 2D grid
xx, yy = np.meshgrid(x, x)

print('Original grid: x=\n',xx)
print('Original grid: y=\n',yy)

# Definition of the 1D field on this grid
ff = np.zeros_like(xx)
ff[::2,::2] = 1  # set ff=1 at even indices

print('Original field: f=\n',ff)

#print('shape:',ff.shape)
#newshape=np.append(ff.shape,2)
#print('newshape:',newshape)

# Definition of the target grid
xxnew, yynew = np.meshgrid(xnew, xnew)

print('Target grid: xnew=\n',xxnew)
print('Target grid: ynew=\n',yynew)

# Define 2D grid in pyensdam
edam.interpolation.define2D(xx,yy)

# Locate target points in original grid
location, weight = edam.interpolation.locate2D(xxnew,yynew)

print('Location indices: location=\n',location)
print('Interpolation weight: weight=\n',weight)

# Interpolation original field at new locations
ffnew = edam.interpolation.interp2D(ff,location,weight)
print('Interpolated field: fnew=\n',ffnew)
