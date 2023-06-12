module interptools

use iso_c_binding, only: c_int, c_double

use ensdam_interp
use ensdam_unmask

implicit none

contains

! Routines to exchange global variables from Fortran modules
subroutine c_get_unmask_spval(var) bind(c)
   real(c_double), intent(out) :: var
   var = unmask_spval
end subroutine

subroutine c_set_unmask_spval(var) bind(c)
   real(c_double), intent(in) :: var
   unmask_spval = var
end subroutine

subroutine c_get_unmask_damping(var) bind(c)
   real(c_double), intent(out) :: var
   var = unmask_damping
end subroutine

subroutine c_set_unmask_damping(var) bind(c)
   real(c_double), intent(in) :: var
   unmask_damping = var
end subroutine

subroutine c_get_unmask_max(var) bind(c)
   integer(c_int), intent(out) :: var
   var = unmask_max
end subroutine

subroutine c_set_unmask_max(var) bind(c)
   integer(c_int), intent(in) :: var
   unmask_max = var
end subroutine

subroutine c_get_unmask_window(var) bind(c)
   integer(c_int), intent(out) :: var
   var = unmask_window
end subroutine

subroutine c_set_unmask_window(var) bind(c)
   integer(c_int), intent(in) :: var
   unmask_window = var
end subroutine

subroutine c_get_unmask_k_ew(var) bind(c)
   integer(c_int), intent(out) :: var
   var = k_ew
end subroutine

subroutine c_set_unmask_k_ew(var) bind(c)
   integer(c_int), intent(in) :: var
   k_ew = var
end subroutine

! Routines from ensdam_interp
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
  integer(c_int), intent(out) :: ix, iy
  integer(c_int), intent(out) :: located
  logical :: f_located

  f_located = grid2D_locate(x,y,ix,iy)

  located = 0
  if (f_located) located = 1
end subroutine

subroutine c_grid1D_interp(ngrid,grid,x,ix,w) bind(c)
  implicit none
  integer(c_int), intent(in), value :: ngrid
  real(c_double), intent(in) :: grid(ngrid)
  real(c_double), intent(in) :: x
  integer(c_int), intent(in) :: ix
  real(c_double), intent(out) :: w

  call grid1D_interp(grid,x,ix,w)
end subroutine

subroutine c_grid2D_interp(x,y,ix,iy,w) bind(c)
  real(c_double), intent(in) :: x, y
  integer(c_int), intent(in) :: ix, iy
  real(c_double), intent(out) :: w(2,2)

  call grid2D_interp(x,y,ix,iy,w)
end subroutine

! Routines from ensdam_unmask
subroutine c_unmask(ni,nj,phi) bind(c)
  integer(c_int), intent(in), value :: ni, nj
  real(c_double), intent(inout) :: phi(ni,nj)

  call unmask(phi)
end subroutine

end module
