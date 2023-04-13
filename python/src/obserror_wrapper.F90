module obserror

use iso_c_binding, only: c_ptr, c_int, c_int32_t, c_double, c_char, c_loc, c_f_pointer, C_NULL_CHAR, C_NULL_PTR

use ensdam_obserror

implicit none

contains

subroutine c_get_obstype(var) bind(c)
  character(kind=c_char), intent(out) :: var(20+1)
  call copy_string_ftoc(obserror_type, var)
end subroutine

subroutine c_set_obstype(len1, var) bind(c)
  integer(c_int32_t), intent(in) :: len1
  character(kind=c_char), intent(in) :: var(len1+1)
  call copy_string_ctof(var, obserror_type)
end subroutine


function c_obserror_logpdf_vector(nobs,y,x,sigma) bind(c)
  integer(c_int), intent(in), value :: nobs
  real(c_double), intent(in) :: y(nobs)
  real(c_double), intent(in) :: x(nobs)
  real(c_double), intent(in) :: sigma(nobs)
  real(c_double) :: c_obserror_logpdf_vector
  c_obserror_logpdf_vector = obserror_logpdf_vector(y, x, sigma)
end function

function c_obserror_logpdf_vector_homogeneous(nobs,y,x,sigma) bind(c)
  integer(c_int), intent(in), value :: nobs
  real(c_double), intent(in) :: y(nobs)
  real(c_double), intent(in) :: x(nobs)
  real(c_double), intent(in) :: sigma
  real(c_double) :: c_obserror_logpdf_vector_homogeneous
  c_obserror_logpdf_vector_homogeneous = obserror_logpdf_vector_homogeneous(y, x, sigma)
end function

function c_obserror_logpdf_variable(y,x,sigma) bind(c)
  real(c_double), intent(in) :: y
  real(c_double), intent(in) :: x
  real(c_double), intent(in) :: sigma
  real(c_double) :: c_obserror_logpdf_variable
  c_obserror_logpdf_variable = obserror_logpdf_variable(y, x, sigma)
end function


subroutine c_obserror_cdf_vector(nobs,y,x,sigma,cdf) bind(c)
  integer(c_int), intent(in), value :: nobs
  real(c_double), intent(in) :: y(nobs)
  real(c_double), intent(in) :: x(nobs)
  real(c_double), intent(in) :: sigma(nobs)
  real(c_double) :: cdf(nobs)
  call obserror_cdf_vector(y, x, sigma, cdf)
end subroutine

subroutine c_obserror_cdf_vector_homogeneous(nobs,y,x,sigma,cdf) bind(c)
  integer(c_int), intent(in), value :: nobs
  real(c_double), intent(in) :: y(nobs)
  real(c_double), intent(in) :: x(nobs)
  real(c_double), intent(in) :: sigma
  real(c_double) :: cdf(nobs)
  call obserror_cdf_vector_homogeneous(y, x, sigma, cdf)
end subroutine

subroutine c_obserror_cdf_variable(y,x,sigma,cdf) bind(c)
  real(c_double), intent(in) :: y
  real(c_double), intent(in) :: x
  real(c_double), intent(in) :: sigma
  real(c_double) :: cdf
  call obserror_cdf_variable(y, x, sigma, cdf)
end subroutine


subroutine c_obserror_sample_vector(nobs,x,sigma,sample) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nobs
  real(c_double), intent(in) :: x(nobs)
  real(c_double), intent(in) :: sigma(nobs)
  real(c_double) :: sample(nobs)
  call obserror_sample_vector(x, sigma, sample)
end subroutine

subroutine c_obserror_sample_vector_homogeneous(nobs,x,sigma,sample,rank,uniform_rank,argcase) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nobs
  real(c_double), intent(in) :: x(nobs)
  real(c_double), intent(in) :: sigma
  real(c_double), intent(in) :: rank
  integer(c_int), intent(in) :: uniform_rank
  integer(c_int), intent(in) :: argcase
  real(c_double) :: sample(nobs)

  logical :: f_uniform_rank

  select case(argcase)
  case(0)
    call obserror_sample_vector_homogeneous(x, sigma, sample)
  case(1)
    call obserror_sample_vector_homogeneous(x, sigma, sample, rank=rank)
  case(2)
    f_uniform_rank = uniform_rank.ne.0
    call obserror_sample_vector_homogeneous(x, sigma, sample, uniform_rank=f_uniform_rank)
  end select

end subroutine

subroutine c_obserror_sample_variable(x,sigma,sample,rank,reuse_last_rank,argcase) bind(c)
  implicit none
  real(c_double), intent(in) :: x
  real(c_double), intent(in) :: sigma
  real(c_double), intent(in) :: rank
  integer(c_int), intent(in) :: reuse_last_rank
  integer(c_int), intent(in) :: argcase
  real(c_double) :: sample

  logical :: f_reuse_last_rank

  select case(argcase)
  case(0)
    call obserror_sample_variable(x, sigma, sample)
  case(1)
    call obserror_sample_variable(x, sigma, sample, rank=rank)
  case(2)
    f_reuse_last_rank = reuse_last_rank.ne.0
    call obserror_sample_variable(x, sigma, sample, reuse_last_rank=f_reuse_last_rank)
  end select

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

subroutine copy_string_ftoc(stringf,stringc)
  ! utility function to convert c string to fortran string
  character(len=*), intent(in) :: stringf
  character(c_char), intent(out) :: stringc(:)
  integer j,n
  n = len_trim(stringf)   
  do j=1,n    
    stringc(j) = stringf(j:j)   
  end do
  stringc(n+1) = c_null_char
end subroutine copy_string_ftoc

end module
