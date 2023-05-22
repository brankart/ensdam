module stochtools

use iso_c_binding, only: c_int, c_long, c_double, c_int32_t, c_funptr, c_char, c_null_char

use ensdam_storng
use ensdam_stotge
use ensdam_stoutil
use ensdam_storfg

implicit none

interface
  function pow_spectrum_callback_interface(l,m) bind(c)
    import c_double, c_int
    integer(c_int), intent(in) :: l, m
    real(c_double) :: pow_spectrum_callback_interface
  end function
end interface
procedure (pow_spectrum_callback_interface), pointer, save :: pow_spectrum => null()

contains

! Routine to associate callback function

subroutine c_associate_pow_spectrum_callback(pow_spectrum_in) bind(c)
   type(c_funptr), intent(in), value :: pow_spectrum_in
   call c_f_procpointer(pow_spectrum_in, pow_spectrum)
end subroutine

! Fortran side of callback function

function f_pow_spectrum(l,m)
   integer, intent(in) :: l,m
   real(kind=8) :: f_pow_spectrum

   f_pow_spectrum = pow_spectrum(l,m)
end function

! From storng

subroutine c_kiss(iran) bind(c)
  implicit none
  integer(c_long), intent(out) :: iran

  iran = kiss()
end subroutine

subroutine c_kiss_seed(ix, iy, iz, iw) bind(c)
  implicit none
  integer(c_long), intent(in), value :: ix, iy, iz, iw

  call kiss_seed(ix, iy, iz, iw)
end subroutine

subroutine c_kiss_save() bind(c)
  call kiss_save()
end subroutine

subroutine c_kiss_load() bind(c)
  call kiss_load()
end subroutine

subroutine c_kiss_reset() bind(c)
  call kiss_reset()
end subroutine

subroutine c_kiss_check(len_check_type,check_type) bind(c)
  implicit none
  integer(c_int32_t), intent(in) :: len_check_type
  character(c_char), intent(in) :: check_type(len_check_type+1)
  character(len=10) :: f_check_type

  call copy_string_ctof(check_type, f_check_type)
  call kiss_check(f_check_type)
end subroutine

subroutine c_kiss_uniform(uran) bind(c)
  implicit none
  real(c_double), intent(out) :: uran

  call kiss_uniform(uran)
end subroutine

subroutine c_kiss_gaussian(gran) bind(c)
  implicit none
  real(c_double), intent(out) :: gran

  call kiss_gaussian(gran)
end subroutine

subroutine c_kiss_gamma(gamr,k) bind(c)
  implicit none
  real(c_double), intent(out) :: gamr
  real(c_double), intent(in) :: k

  call kiss_gamma(gamr,k)
end subroutine

subroutine c_kiss_beta(betar,a,b) bind(c)
  implicit none
  real(c_double), intent(out) :: betar
  real(c_double), intent(in) :: a, b

  call kiss_beta(betar,a,b)
end subroutine

subroutine c_kiss_sample(a,n,k) bind(c)
  implicit none
  integer(c_int), intent(inout) :: a(n)
  integer(c_int), intent(in), value :: n
  integer(c_int), intent(in), value :: k

  call kiss_sample(a,n,k)
end subroutine

! From stotge

subroutine c_ran_te(teran,a) bind(c)
  implicit none
  real(c_double), intent(out) :: teran
  real(c_double), intent(in) :: a

  call ran_te(teran,a)
end subroutine

subroutine c_ran_tg(nsmpl,tgsmpl,aa,bb) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nsmpl
  real(c_double), intent(out) :: tgsmpl(nsmpl)
  real(c_double), intent(in) :: aa, bb

  call ran_tg(tgsmpl,aa,bb)
end subroutine

subroutine c_ranv_tg(nvar,ncstr,nsmpl,tgvsmpl,matArm,vecbm) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nvar
  integer(c_int), intent(in), value :: ncstr
  integer(c_int), intent(in), value :: nsmpl
  real(c_double), intent(out) :: tgvsmpl(nsmpl,nvar)
  real(c_double), intent(in) :: matArm(nvar,ncstr)
  real(c_double), intent(in) :: vecbm(ncstr)

  call ranv_tg(tgvsmpl,matArm,vecbm)
end subroutine

! From storfg

subroutine c_def_spect_init(nfreq,nspct1d,nspct2d,nspct2s) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nfreq
  integer(c_int), intent(in), value :: nspct1d
  integer(c_int), intent(in), value :: nspct2d
  integer(c_int), intent(in), value :: nspct2s

  call def_spect_init(nfreq,nspct1d,nspct2d,nspct2s)
end subroutine

subroutine c_def_spect_power(spct_type,spct_idx,nspct,spct_freq,spct_power) bind(c)
  implicit none
  integer(c_int), intent(in), value :: spct_type
  integer(c_int), intent(in), value :: spct_idx
  integer(c_int), intent(in), value :: nspct
  real(c_double), intent(in) :: spct_freq(nspct)
  real(c_double), intent(in) :: spct_power(nspct)

  call def_spect_power(spct_type,spct_idx,spct_freq,spct_power)
end subroutine

subroutine c_def_sample_size(nsmp1d,nsmp2d,nsmp2s,nseed) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nsmp1d
  integer(c_int), intent(in), value :: nsmp2d
  integer(c_int), intent(in), value :: nsmp2s
  integer(c_int), intent(in), value :: nseed

  call def_sample_size(nsmp1d,nsmp2d,nsmp2s,nseed)
end subroutine

subroutine c_sample_freq_1d(spct_idx)
  implicit none
  integer(c_int), intent(in), value :: spct_idx

  call sample_freq_1d(spct_idx)
end subroutine

subroutine c_sample_freq_2d(spct_idx)
  implicit none
  integer(c_int), intent(in), value :: spct_idx

  call sample_freq_2d(spct_idx)
end subroutine

subroutine c_gen_field_1d(spct_idx,nx,ranfield,x)
  implicit none
  integer(c_int), intent(in), value :: spct_idx
  integer(c_int), intent(in), value :: nx
  real(c_double), intent(out) :: ranfield(nx)
  real(c_double), intent(in) :: x(nx)

  call gen_field_1d(spct_idx,ranfield,x)
end subroutine

subroutine c_gen_field_2d(spct_idx,nx,ny,ranfield,x,y)
  implicit none
  integer(c_int), intent(in), value :: spct_idx
  integer(c_int), intent(in), value :: nx, ny
  real(c_double), intent(out) :: ranfield(nx,ny)
  real(c_double), intent(in) :: x(nx,ny), y(nx,ny)

  call gen_field_2d(spct_idx,ranfield,x,y)
end subroutine

subroutine c_gen_field_2s(ngrid,ranfield,lon,lat,lmin,lmax)
  implicit none
  integer(c_int), intent(in), value :: ngrid
  real(c_double), intent(out) :: ranfield(ngrid)
  real(c_double), intent(in) :: lon(ngrid)
  real(c_double), intent(in) :: lat(ngrid)
  integer(c_int), intent(in), value :: lmin, lmax

  call gen_field_2s_new(ranfield,lon,lat,f_pow_spectrum,lmin,lmax)
end subroutine

! From stoutil

subroutine c_cdf_gaussian(x,fun) bind(c)
  implicit none
  real(c_double), intent(in) :: x
  real(c_double), intent(out) :: fun

  fun = cdf_gaussian(x)
end subroutine

subroutine c_pdf_gaussian(x,fun) bind(c)
  implicit none
  real(c_double), intent(in) :: x
  real(c_double), intent(out) :: fun

  fun = pdf_gaussian(x)
end subroutine

subroutine c_logpdf_gaussian(x,fun) bind(c)
  implicit none
  real(c_double), intent(in) :: x
  real(c_double), intent(out) :: fun

  fun = logpdf_gaussian(x)
end subroutine

subroutine c_invcdf_gaussian(rank,fun) bind(c)
  implicit none
  real(c_double), intent(in) :: rank
  real(c_double), intent(out) :: fun

  fun = cdf_gaussian(rank)
end subroutine

subroutine c_cdf_gamma(a,x,fun) bind(c)
  implicit none
  real(c_double), intent(in) :: a,x
  real(c_double), intent(out) :: fun

  fun = cdf_gamma(a,x)
end subroutine

subroutine c_pdf_gamma(a,x,fun) bind(c)
  implicit none
  real(c_double), intent(in) :: a,x
  real(c_double), intent(out) :: fun

  fun = pdf_gamma(a,x)
end subroutine

subroutine c_logpdf_gamma(k,theta,x,fun) bind(c)
  implicit none
  real(c_double), intent(in) :: k,theta,x
  real(c_double), intent(out) :: fun

  fun = logpdf_gamma(k,theta,x)
end subroutine

subroutine c_invcdf_gamma(a,rank,fun) bind(c)
  implicit none
  real(c_double), intent(in) :: a,rank
  real(c_double), intent(out) :: fun

  fun = invcdf_gamma(a,rank)
end subroutine

subroutine c_cdf_beta(a,b,x,fun) bind(c)
  implicit none
  real(c_double), intent(in) :: a,b,x
  real(c_double), intent(out) :: fun

  fun = cdf_beta(a,b,x)
end subroutine

subroutine c_pdf_beta(a,b,x,fun) bind(c)
  implicit none
  real(c_double), intent(in) :: a,b,x
  real(c_double), intent(out) :: fun

  fun = pdf_beta(a,b,x)
end subroutine

subroutine c_logpdf_beta(a,b,x,fun) bind(c)
  implicit none
  real(c_double), intent(in) :: a,b,x
  real(c_double), intent(out) :: fun

  fun = logpdf_beta(a,b,x)
end subroutine

subroutine c_invcdf_beta(a,b,rank,fun) bind(c)
  implicit none
  real(c_double), intent(in) :: a,b,rank
  real(c_double), intent(out) :: fun

  fun = invcdf_beta(a,b,rank)
end subroutine

!!!!!!!!!!!!!!!!!!
!!! Utilities  !!!
!!!!!!!!!!!!!!!!!!

subroutine copy_string_ctof(stringc,stringf)
  ! utility function to convert c string to fortran string
  character(len=*), intent(out) :: stringf
  character(c_char), intent(in) :: stringc(:)
  integer j
  stringf = ''
  char_loop: do j=1,min(size(stringc),len(stringf))
     if (stringc(j)==c_null_char) exit char_loop
     stringf(j:j) = stringc(j)
  end do char_loop
end subroutine copy_string_ctof

end module
