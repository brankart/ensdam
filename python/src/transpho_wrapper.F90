module transpho

use iso_c_binding, only: c_int, c_double

use ensdam_sphylm
use ensdam_spharea

implicit none

contains

! From sphylm

subroutine c_init_ylm(nl,lmin,latmin,latmax,dlatmax) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nl
  integer(c_int), intent(in), value :: lmin
  real(c_double), intent(in) :: latmin, latmax, dlatmax

  call init_ylm(nl,lmin,latmin,latmax,dlatmax)
end subroutine

subroutine c_proj_ylm(nvar,nl,proj,tab,lon,lat) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nvar
  integer(c_int), intent(in), value :: nl
  !real(c_double), intent(out) :: proj(nl+1,2*nl+1)
  real(c_double), intent(out) :: proj(0:nl,-nl:nl)
  real(c_double), intent(in) :: tab(nvar)
  real(c_double), intent(in) :: lon(nvar)
  real(c_double), intent(in) :: lat(nvar)

  call proj_ylm(proj,tab,lon,lat)
end subroutine

subroutine c_back_ylm(nvar,nl,proj,tab,lon,lat) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nvar
  integer(c_int), intent(in), value :: nl
  !real(c_double), intent(in) :: proj(nl+1,2*nl+1)
  real(c_double), intent(out) :: proj(0:nl,-nl:nl)
  real(c_double), intent(out) :: tab(nvar)
  real(c_double), intent(in) :: lon(nvar)
  real(c_double), intent(in) :: lat(nvar)

  call back_ylm(proj,tab,lon,lat)
end subroutine

subroutine c_back_ylm_loc(nvar,nl,proj,tab,lon,lat,l0,l1) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nvar
  integer(c_int), intent(in), value :: nl
  real(c_double), intent(in) :: proj(nl+1,2*nl+1)
  real(c_double), intent(out) :: tab(nvar)
  real(c_double), intent(in) :: lon(nvar)
  real(c_double), intent(in) :: lat(nvar)
  integer(c_int), intent(in) :: l0, l1

  call back_ylm_loc(proj,tab,lon,lat,l0,l1)
end subroutine

! From spharea

subroutine c_mesh_area(ni,nj,lon,lat,area) bind(c)
  implicit none
  integer(c_int), intent(in), value :: ni, nj
  real(c_double), intent(in) :: lon(ni,nj)
  real(c_double), intent(in) :: lat(ni,nj)
  real(c_double), intent(out) :: area(ni,nj)

  call mesh_area (lon,lat,area)
end subroutine

end module
