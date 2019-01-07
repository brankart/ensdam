program random_field_on_the_sphere
  use netcdf
  use ensdam_storfg
  implicit none

  integer, parameter :: nlon=360 , nlat=180
  real(kind=8), dimension(nlon,nlat) :: x, y, ranfield

  integer, parameter :: lmax=45   ! maximum degree of the spherical harmonics
  real(kind=8), parameter :: lm=0   ! scale of power spectrum maximum
  real(kind=8), parameter :: lc=6.4 ! power spectrum characteristic scale

  integer :: i, j, n
  integer :: is, ncid, idx, idy, idv, idlon, idlat

  ! Define the output grid
  do j=1,nlat
    x(:,j) = (/ (i, i=0,nlon-1) /)
  enddo
  do i=1,nlon
    y(i,:) = (/ (j, j=0,nlat-1) /)
  enddo

  x(:,:) = x(:,:) * 360. / nlon + 360. / nlon / 2
  y(:,:) = y(:,:) * 180. / nlat - 90. + 180. / nlat / 2

  ! Generate random field
  n = nlon*nlat ! size of 1d vectors
  call gen_field_2s(ranfield,x,y,pow_spectrum,0,lmax)

  ! Write output file (in NetCDF)
  is = NF90_CREATE("random_field_on_the_sphere.nc",NF90_CLOBBER,ncid)
  is = NF90_DEF_DIM(ncid,'lon',nlon+1,idx)
  is = NF90_DEF_DIM(ncid,'lat',nlat+1,idy)
  is = NF90_DEF_VAR(ncid,'lon',NF90_FLOAT,(/idx,idy/),idlon)
  is = NF90_DEF_VAR(ncid,'lat',NF90_FLOAT,(/idx,idy/),idlat)
  is = NF90_DEF_VAR(ncid,'ranfield',NF90_FLOAT,(/idx,idy/),idv)
  is = NF90_ENDDEF(ncid)
  is = NF90_PUT_VAR(ncid,idlon,x)
  is = NF90_PUT_VAR(ncid,idlat,y)
  is = NF90_PUT_VAR(ncid,idv,ranfield)
  is = NF90_CLOSE(ncid)

contains

  ! Power spectrum in the basis of the spherical harmonics
  ! (as a function of the degree l)
  function pow_spectrum(l,m)
  implicit none
  integer, intent(in) :: l,m
  real(kind=8) :: pow_spectrum

  integer :: ll
  real(kind=8) :: norm

  ! Power spectrum
  pow_spectrum = 1. / ( 1. + (l-lm)*(l-lm)/(lc*lc) )

  ! Normalize the spectrum
  norm = 0.
  DO ll=0,lmax
    norm = norm + 1. / ( 1. + (ll-lm)*(ll-lm)/(lc*lc) )
  ENDDO
  pow_spectrum = pow_spectrum / norm

  ! Scale to account for the multiplicity of each degree
  pow_spectrum = pow_spectrum / ( 1. + 2. * l )

  end function pow_spectrum

end
