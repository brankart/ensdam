module ensaugm

use iso_c_binding, only: c_int, c_double

use ensdam_ensaugm

implicit none

contains

! Routines to exchange global variables from Fortran module
subroutine c_get_augment_chain_index(var) bind(c)
   integer(c_int), intent(out) :: var
   var = ensaugm_chain_index
end subroutine

subroutine c_set_augment_chain_index(var) bind(c)
   integer(c_int), intent(out) :: var
   ensaugm_chain_index = var
end subroutine

! Routine to sample augmented ensemble
subroutine c_sample_augmented_ensemble(nvar,nens,nres,naug,maxchain,augens,ens,multiplicity) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nvar
  integer(c_int), intent(in), value :: nens
  integer(c_int), intent(in), value :: nres
  integer(c_int), intent(in), value :: naug
  integer(c_int), intent(in) :: maxchain
  real(c_double), intent(out) :: augens(nvar,naug)
  real(c_double), intent(in) :: ens(nvar,nens,nres)
  integer(c_int), intent(in) :: multiplicity(nres)

  call sample_augmented_ensemble(maxchain,augens,ens,multiplicity)
end subroutine

end module
