module enssanam

use iso_c_binding, only: c_int, c_double

use ensdam_anaqua
use ensdam_anatra
use ensdam_anaobs

implicit none

contains

! From anaqua

subroutine c_ens_quantiles_vector(nvar,nens,nqua,qua,ens,quadef,enswei,ensweiloc,argcase) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nvar
  integer(c_int), intent(in), value :: nens
  integer(c_int), intent(in), value :: nqua
  real(c_double), intent(out) :: qua(nvar,nqua)
  real(c_double), intent(in) :: ens(nvar,nens)
  real(c_double), intent(in) :: quadef(nqua)
  real(c_double), intent(in) :: enswei(nens)
  real(c_double), intent(in) :: ensweiloc(nvar,nens)
  integer(c_int), intent(in) :: argcase

  select case(argcase)
  case(0)
    call ens_quantiles_vector(qua,ens,quadef)
  case(1)
    call ens_quantiles_vector(qua,ens,quadef,enswei=enswei)
  case(2)
    call ens_quantiles_vector(qua,ens,quadef,ensweiloc=ensweiloc)
  end select

end subroutine

subroutine c_ens_quantiles_variable(nens,nqua,qua,ens,quadef,enswei,argcase) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nens
  integer(c_int), intent(in), value :: nqua
  real(c_double), intent(out) :: qua(nqua)
  real(c_double), intent(in) :: ens(nens)
  real(c_double), intent(in) :: quadef(nqua)
  real(c_double), intent(in) :: enswei(nens)
  integer(c_int), intent(in) :: argcase

  select case(argcase)
  case(0)
    call ens_quantiles_variable(qua,ens,quadef)
  case(1)
    call ens_quantiles_variable(qua,ens,quadef,enswei=enswei)
  end select

end subroutine

! From anatra

subroutine c_ana_forward_ensemble(nvar,nens,nqua,ens,qua,quaref,rank,argcase) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nvar
  integer(c_int), intent(in), value :: nens
  integer(c_int), intent(in), value :: nqua
  real(c_double), intent(inout) :: ens(nvar,nens)
  real(c_double), intent(in) :: qua(nvar,nqua)
  real(c_double), intent(in) :: quaref(nqua)
  real(c_double), intent(in) :: rank(nens)
  integer(c_int), intent(in) :: argcase

  select case(argcase)
  case(0)
    call ana_forward_ensemble(ens,qua,quaref)
  case(1)
    call ana_forward_ensemble(ens,qua,quaref,rank=rank)
  end select

end subroutine

subroutine c_ana_forward_vector(nvar,nqua,vct,qua,quaref,rank,argcase) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nvar
  integer(c_int), intent(in), value :: nqua
  real(c_double), intent(inout) :: vct(nvar)
  real(c_double), intent(in) :: qua(nvar,nqua)
  real(c_double), intent(in) :: quaref(nqua)
  real(c_double), intent(in) :: rank
  integer(c_int), intent(in) :: argcase

  select case(argcase)
  case(0)
    call ana_forward_vector(vct,qua,quaref)
  case(1)
    call ana_forward_vector(vct,qua,quaref,rank=rank)
  end select

end subroutine

subroutine c_ana_forward_variable(nqua,var,qua,quaref,rank,argcase) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nqua
  real(c_double), intent(inout) :: var
  real(c_double), intent(in) :: qua(nqua)
  real(c_double), intent(in) :: quaref(nqua)
  real(c_double), intent(in) :: rank
  integer(c_int), intent(in) :: argcase

  select case(argcase)
  case(0)
    call ana_forward_variable(var,qua,quaref)
  case(1)
    call ana_forward_variable(var,qua,quaref,rank=rank)
  end select

end subroutine

subroutine c_ana_backward_ensemble(nvar,nens,nqua,ens,qua,quaref) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nvar
  integer(c_int), intent(in), value :: nens
  integer(c_int), intent(in), value :: nqua
  real(c_double), intent(inout) :: ens(nvar,nens)
  real(c_double), intent(in) :: qua(nvar,nqua)
  real(c_double), intent(in) :: quaref(nqua)

  call ana_backward_ensemble(ens,qua,quaref)
end subroutine

subroutine c_ana_backward_vector(nvar,nqua,vct,qua,quaref) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nvar
  integer(c_int), intent(in), value :: nqua
  real(c_double), intent(inout) :: vct(nvar)
  real(c_double), intent(in) :: qua(nvar,nqua)
  real(c_double), intent(in) :: quaref(nqua)

  call ana_backward_vector(vct,qua,quaref)
end subroutine

subroutine c_ana_backward_variable(nqua,var,qua,quaref) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nqua
  real(c_double), intent(inout) :: var
  real(c_double), intent(in) :: qua(nqua)
  real(c_double), intent(in) :: quaref(nqua)

  call ana_backward_variable(var,qua,quaref)
end subroutine

! From anaobs

subroutine c_ana_obs(nobs,nens,nqua,nsmp,anaobs,obsens,obs,obserror,quadef,quaref) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nobs
  integer(c_int), intent(in), value :: nens
  integer(c_int), intent(in), value :: nqua
  integer(c_int), intent(in), value :: nsmp
  real(c_double), intent(out) :: anaobs(nobs,nsmp)
  real(c_double), intent(in) :: obsens(nobs,nens)
  real(c_double), intent(in) :: obs(nobs)
  real(c_double), intent(in) :: obserror(nobs)
  real(c_double), intent(in) :: quadef(nqua)
  real(c_double), intent(in) :: quaref(nqua)

  call ana_obs(anaobs,obsens,obs,obserror,quadef,quaref)
end subroutine

subroutine c_ana_obs_sym(nobs,nqua,nsmp,anaobs,obs,obserror,obsqua,quaref) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nobs
  integer(c_int), intent(in), value :: nqua
  integer(c_int), intent(in), value :: nsmp
  real(c_double), intent(out) :: anaobs(nobs,nsmp)
  real(c_double), intent(in) :: obs(nobs)
  real(c_double), intent(in) :: obserror(nobs)
  real(c_double), intent(in) :: obsqua(nobs,nqua)
  real(c_double), intent(in) :: quaref(nqua)

  call ana_obs_sym(anaobs,obs,obserror,obsqua,quaref)
end subroutine

end module
