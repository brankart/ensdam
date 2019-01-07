program observation_regression_on_the_sphere
  USE netcdf
  USE ensdam_sphylm
  USE ensdam_storng
  IMPLICIT none

  INTEGER, PARAMETER :: lmin=0, lmax=10
  INTEGER, PARAMETER :: nlon=360 , nlat=180
  INTEGER, PARAMETER :: obsratio=50, nobs=nlon*nlat/obsratio ! every 50 grid point is observed
  REAL(KIND=8) :: obserror_std=0.05
  REAL(KIND=8) :: freq0=15
  REAL(KIND=8), DIMENSION(nlon,nlat) :: x, y, ranfield, area
  REAL(KIND=8), DIMENSION(nlon+1,nlat+1) :: xc, yc
  REAL(KIND=8), DIMENSION(0:lmax,-lmax:lmax) :: spectrum, std_spectrum
  REAL(KIND=8), DIMENSION(nlon*nlat) :: x_1d, y_1d, ranfield_1d
  REAL(KIND=8), DIMENSION(nobs) :: x_obs, y_obs, obs, invstd_obs

  INTEGER :: i, j, l, iobs

  INTEGER :: is, ncid, idx, idy, idv, idlon, idlat
  REAL(KIND=8) :: latmin, latmax, dlatmax, gran, norm

  ! Read input file (in NetCDF)
  is = NF90_OPEN("low_pass_filter_on_the_sphere.nc",NF90_NOWRITE,ncid)
  is = NF90_INQ_VARID(ncid,'lon',idlon)
  is = NF90_INQ_VARID(ncid,'lat',idlat)
  is = NF90_INQ_VARID(ncid,'ranfield',idv)
  is = NF90_GET_VAR(ncid,idlon,x)
  is = NF90_GET_VAR(ncid,idlat,y)
  is = NF90_GET_VAR(ncid,idv,ranfield)
  is = NF90_CLOSE(ncid)

  ! Define the corners of the grid cells
  do j=1,nlat+1
    xc(:,j) = (/ (i, i=0,nlon) /)
  enddo
  do i=1,nlon+1
    yc(i,:) = (/ (j, j=0,nlat) /)
  enddo

  xc(:,:) = xc(:,:) * 360. / nlon
  yc(:,:) = yc(:,:) * 180. / nlat - 90.

  ! Initialize computation of spherical harmonics (up to degree lmax)
  latmin = -90. ; latmax = 90. ; dlatmax = 0.1
  CALL init_ylm( lmax, lmin, latmin, latmax, dlatmax )

  ! Sample observations from the input field
  x_1d=RESHAPE(x,(/nlon*nlat/))
  y_1d=RESHAPE(y,(/nlon*nlat/))
  ranfield_1d=RESHAPE(ranfield,(/nlon*nlat/))
  DO iobs=1,nobs
    x_obs(iobs) = x_1d(obsratio*iobs)
    y_obs(iobs) = y_1d(obsratio*iobs)
    obs(iobs) = ranfield_1d(obsratio*iobs)
  ENDDO

  ! Add obervation error
  DO iobs=1,nobs
    CALL kiss_gaussian(gran)
    obs(iobs) = obs(iobs) + gran * obserror_std
  ENDDO

  ! Observation error std
  invstd_obs(:) = 1.0 / obserror_std

  ! Signal std in spectral space
  std_spectrum(:,:) = 0.0
  DO l=0,lmax
    std_spectrum(l,-l:l) = 1.0 / ( 1.0 + l*l/(freq0*freq0) )
  ENDDO
  norm = SQRT(SUM(std_spectrum(:,:)*std_spectrum(:,:)))
  std_spectrum(:,:) = std_spectrum(:,:) / norm

  ! Regression of observations on the spherical harmonics
  ! call init_regr_ylm( 'local', 50, 1, 1, 0.01_8, 1.0_8) ! default parameters are used
  call regr_ylm(spectrum,std_spectrum,obs,x_obs,y_obs,invstd_obs)

  ! Transform back into original space (up to degree lmax)
  call back_ylm(spectrum,ranfield_1d,x_1d,y_1d)
  x=RESHAPE(x_1d,(/nlon,nlat/))
  y=RESHAPE(y_1d,(/nlon,nlat/))
  ranfield=RESHAPE(ranfield_1d,(/nlon,nlat/))

  ! Write output file (in NetCDF)
  is = NF90_CREATE("observation_regression_on_the_sphere.nc",NF90_CLOBBER,ncid)
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

end
