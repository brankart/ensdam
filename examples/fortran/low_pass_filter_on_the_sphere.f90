program low_pass_filter_on_the_sphere
  USE netcdf
  USE ensdam_spharea
  USE ensdam_sphylm
  IMPLICIT none

  INTEGER, PARAMETER :: lmin=0, lmax=10
  INTEGER, PARAMETER :: nlon=360 , nlat=180
  REAL(KIND=8), DIMENSION(nlon,nlat) :: x, y, ranfield, area
  REAL(KIND=8), DIMENSION(nlon+1,nlat+1) :: xc, yc
  REAL(KIND=8), DIMENSION(0:lmax,-lmax:lmax) :: spectrum
  REAL(KIND=8), DIMENSION(nlon*nlat) :: x_1d, y_1d, ranfield_1d

  INTEGER :: i, j, l

  INTEGER :: is, ncid, idx, idy, idv, idlon, idlat
  REAL(KIND=8) :: latmin, latmax, dlatmax

  ! Read input file (in NetCDF)
  is = NF90_OPEN("random_field_on_the_sphere.nc",NF90_NOWRITE,ncid)
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

  ! Compute area of each grid cell
  call mesh_area(xc,yc,area)
  
  ! Weight input field with area of each grid cell
  ranfield = ranfield * area
  
  ! Initialize computation of spherical harmonics (up to degree lmax)
  latmin = -90. ; latmax = 90. ; dlatmax = 0.1
  CALL init_ylm( lmax, lmin, latmin, latmax, dlatmax )

  ! Compute spectrum of the input field up to degree lmax
  x_1d=RESHAPE(x,(/nlon*nlat/))
  y_1d=RESHAPE(y,(/nlon*nlat/))
  ranfield_1d=RESHAPE(ranfield,(/nlon*nlat/))
  call proj_ylm(spectrum,ranfield_1d,x_1d,y_1d)

  ! Transform back into original space (up to degree lmax)
  call back_ylm(spectrum,ranfield_1d,x_1d,y_1d)
  x=RESHAPE(x_1d,(/nlon,nlat/))
  y=RESHAPE(y_1d,(/nlon,nlat/))
  ranfield=RESHAPE(ranfield_1d,(/nlon,nlat/))

  ! Write output file (in NetCDF)
  is = NF90_CREATE("low_pass_filter_on_the_sphere.nc",NF90_CLOBBER,ncid)
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
