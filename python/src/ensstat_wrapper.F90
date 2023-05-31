module obserror

use iso_c_binding, only: c_int, c_double

use ensdam_meanstd
use ensdam_covariance

implicit none

contains

subroutine c_ensemble_meanstd_vector(nstate, nens, ens, mean, std, weight, argcase) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nstate
  integer(c_int), intent(in), value :: nens
  real(c_double), intent(in) :: ens(nstate,nens)
  real(c_double), intent(out) :: mean(nstate)
  real(c_double), intent(out) :: std(nstate)
  real(c_double), intent(in) :: weight(nstate,nens)
  integer(c_int), intent(in) :: argcase

  select case(argcase)
  case(0)
    call ensemble_meanstd_vector(ens,mean)
  case(1)
    call ensemble_meanstd_vector(ens,mean,std=std)
  case(2)
    call ensemble_meanstd_vector(ens,mean,weight=weight)
  case(3)
    call ensemble_meanstd_vector(ens,mean,std=std,weight=weight)
  end select

end subroutine

subroutine c_ensemble_meanstd_variable(nens, ens, mean, std, weight, argcase) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nens
  real(c_double), intent(in) :: ens(nens)
  real(c_double), intent(out) :: mean
  real(c_double), intent(out) :: std
  real(c_double), intent(in) :: weight(nens)
  integer(c_int), intent(in) :: argcase

  select case(argcase)
  case(0)
    call ensemble_meanstd_variable(ens,mean)
  case(1)
    call ensemble_meanstd_variable(ens,mean,std=std)
  case(2)
    call ensemble_meanstd_variable(ens,mean,weight=weight)
  case(3)
    call ensemble_meanstd_variable(ens,mean,std=std,weight=weight)
  end select

end subroutine

subroutine c_ensemble_correlation(nstate, nens, ens, ensref, correl, weight, argcase) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nstate
  integer(c_int), intent(in), value :: nens
  real(c_double), intent(in) :: ens(nstate,nens)
  real(c_double), intent(in) :: ensref(nens)
  real(c_double), intent(out) :: correl(nstate)
  real(c_double), intent(in) :: weight(nens)
  integer(c_int), intent(in) :: argcase

  select case(argcase)
  case(0)
    call ensemble_correlation(ens, ensref, correl)
  case(1)
    call ensemble_correlation(ens, ensref, correl, weight=weight)
  end select

end subroutine

subroutine c_ensemble_representer(nstate, nens, ens, ensref, representer, weight, argcase) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nstate
  integer(c_int), intent(in), value :: nens
  real(c_double), intent(in) :: ens(nstate,nens)
  real(c_double), intent(in) :: ensref(nens)
  real(c_double), intent(out) :: representer(nstate)
  real(c_double), intent(in) :: weight(nens)
  integer(c_int), intent(in) :: argcase

  select case(argcase)
  case(0)
    call ensemble_representer(ens, ensref, representer)
  case(1)
    call ensemble_representer(ens, ensref, representer, weight=weight)
  end select

end subroutine

subroutine c_ensemble_covariance(nstate, nens, ens, ensref, cov, weight, argcase) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nstate
  integer(c_int), intent(in), value :: nens
  real(c_double), intent(in) :: ens(nstate,nens)
  real(c_double), intent(in) :: ensref(nens)
  real(c_double), intent(out) :: cov(nstate)
  real(c_double), intent(in) :: weight(nens)
  integer(c_int), intent(in) :: argcase

  select case(argcase)
  case(0)
    call ensemble_covariance(ens, ensref, cov)
  case(1)
    call ensemble_covariance(ens, ensref, cov, weight=weight)
  end select

end subroutine

end module

