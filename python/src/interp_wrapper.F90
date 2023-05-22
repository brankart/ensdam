module interptools

use iso_c_binding, only: c_int, c_double

use ensdam_interp

implicit none

contains

subroutine c_grid2D_init(ni,nj,xgrid,ygrid,gtype) bind(c)
  implicit none
  integer(c_int), intent(in), value :: ni, nj
  real(c_double), intent(in) :: xgrid(ni,nj)
  real(c_double), intent(in) :: ygrid(ni,nj)
  integer(c_int), intent(in) :: gtype
  character(len=20) :: f_gtype

  select case(gtype)
  case(0)
    f_gtype = 'cartesian'
  case(1)
    f_gtype = 'spherical'
  end select

  call grid2D_init(xgrid,ygrid,f_gtype)
end subroutine

subroutine c_grid1D_locate(ngrid,grid,x,ix,located) bind(c)
  implicit none
  integer(c_int), intent(in), value :: ngrid
  real(c_double), intent(in) :: grid(ngrid)
  real(c_double), intent(in) :: x
  integer(c_int), intent(out) :: ix
  integer(c_int), intent(out) :: located
  logical :: f_located

  f_located = grid1D_locate(grid,x,ix)

  located = 0
  if (f_located) located = 1
end subroutine

subroutine c_grid2D_locate(x,y,ix,iy,located) bind(c)
  implicit none
  real(c_double), intent(in) :: x, y
  integer(c_int), intent(out) :: ix, ix
  integer(c_int), intent(out) :: located
  logical :: f_located

  f_located = grid2D_locate(x,iy,ix,iy)

  located = 0
  if (f_located) located = 1
end subroutine

subroutine c_grid1D_interp(ngrid,grid,x,ix,w) bind(c)
  implicit none
  integer(c_int), intent(in), value :: ngrid
  real(c_double), intent(in) :: grid(ngrid)
  real(c_double), intent(in) :: x
  integer(c_int), intent(out) :: ix
  real(c_double), intent(out) :: w

  call c_grid1D_interp(grid,x,ix,w)
end subroutine

subroutine c_grid2D_interp(x,y,ix,iy,w) bind(c)
  real(c_double), intent(in) :: x, y
  integer(c_int), intent(in) :: ix, ix
  real(c_double), intent(out) :: w(2,2)

  call c_grid2D_interp(x,y,ix,iy,w)
end subroutine

end module
