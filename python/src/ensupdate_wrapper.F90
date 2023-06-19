module ensupdate

use iso_c_binding, only: c_int, c_double, c_funptr, c_f_procpointer

use ensdam_mcmc_update

implicit none

! Interfaces for callback routines
interface
   function my_jo_callback_interface(nvar,v) bind(c)
      import c_double, c_int
      integer(c_int), intent(in), value :: nvar
      real(c_double), intent(in) :: v(nvar)
      real(c_double) :: my_jo_callback_interface
   end function
end interface
procedure (my_jo_callback_interface), pointer, save :: my_jo => null()

interface
   function my_test_callback_interface(nvar,nens,nextra,upens,upxens) bind(c)
      import c_double, c_int
      integer(c_int), intent(in), value :: nvar
      integer(c_int), intent(in), value :: nens
      integer(c_int), intent(in), value :: nextra
      real(c_double), intent(in) :: upens(nvar,nens)
      real(c_double), intent(in) :: upxens(nvar,nextra)
      integer(c_int) :: my_test_callback_interface
   end function
end interface
procedure (my_test_callback_interface), pointer, save :: my_test => null()

contains

! Routines to exchange global variables from Fortran module
subroutine c_get_mcmc_index(var) bind(c)
   integer(c_int), intent(out) :: var
   var = mcmc_index
end subroutine

subroutine c_set_mcmc_index(var) bind(c)
   integer(c_int), intent(in) :: var
   mcmc_index = var
end subroutine

subroutine c_get_mcmc_zero_start(var) bind(c)
   integer(c_int), intent(out) :: var
   var = 0
   if (mcmc_zero_start) var = 1
end subroutine

subroutine c_set_mcmc_zero_start(var) bind(c)
   integer(c_int), intent(in) :: var
   mcmc_zero_start = (var/=0)
end subroutine

subroutine c_get_mcmc_control_print(var) bind(c)
   integer(c_int), intent(out) :: var
   var = mcmc_control_print
end subroutine

subroutine c_set_mcmc_control_print(var) bind(c)
   integer(c_int), intent(in) :: var
   mcmc_control_print = var
end subroutine

subroutine c_get_mcmc_convergence_check(var) bind(c)
   integer(c_int), intent(out) :: var
   var = mcmc_convergence_check
end subroutine

subroutine c_set_mcmc_convergence_check(var) bind(c)
   integer(c_int), intent(in) :: var
   mcmc_convergence_check = var
end subroutine

subroutine c_get_mcmc_convergence_stop(var) bind(c)
   integer(c_int), intent(out) :: var
   var = 0
   if (mcmc_convergence_stop) var = 1
end subroutine

subroutine c_set_mcmc_convergence_stop(var) bind(c)
   integer(c_int), intent(in) :: var
   mcmc_convergence_stop = (var/=0)
end subroutine

! Routines to associate callback functions

subroutine c_associate_my_jo_callback(my_jo_in) bind(c)
   type(c_funptr), intent(in), value :: my_jo_in
   call c_f_procpointer(my_jo_in, my_jo)
end subroutine

subroutine c_associate_my_test_callback(my_test_in) bind(c)
   type(c_funptr), intent(in), value :: my_test_in
   call c_f_procpointer(my_test_in, my_test)
end subroutine

! Fortran side of callback functions

function f_my_jo(v)
   implicit none
   real(kind=8), dimension(:), intent(in) :: v
   real(kind=8) :: f_my_jo
   integer :: nvar

   nvar = size(v,1)
   f_my_jo = my_jo(nvar,v)
end function

function f_my_test(upens,upxens)
   implicit none
   real(kind=8), dimension(:,:), intent(in) :: upens
   real(kind=8), dimension(:,:), intent(in), optional :: upxens
   logical :: f_my_test
   integer :: nvar, nens, nextra

   nvar = size(upens,1)
   nens = size(upens,2)
   nextra = 0

   if (present(upxens)) then
     nextra = size(upxens,1)
     f_my_test = my_test(nvar, nens, nextra, upens, upxens) /= 0
   else
     f_my_test = my_test(nvar, nens, nextra, upens, upens)  /= 0
   endif

end function

! Routine to sample posterior ensemble
subroutine c_mcmc_iteration(nvar,nens,nres,nup,nextra,maxchain,upens,ens,multiplicity,upxens,xens,argcase) bind(c)
  implicit none
  integer(c_int), intent(in), value :: nvar
  integer(c_int), intent(in), value :: nens
  integer(c_int), intent(in), value :: nres
  integer(c_int), intent(in), value :: nup
  integer(c_int), intent(in), value :: nextra
  integer(c_int), intent(in) :: maxchain
  real(c_double), intent(out) :: upens(nvar,nup)
  real(c_double), intent(in) :: ens(nvar,nens,nres)
  integer(c_int), intent(in) :: multiplicity(nres)
  real(c_double), intent(out) :: upxens(nextra,nup)
  real(c_double), intent(in) :: xens(nextra,nens,nres)
  integer(c_int), intent(in) :: argcase

  if (.not.associated(my_jo)) then
    stop 'callback function my_jo not associated in ensupdate'
  endif

  select case(argcase)
  case(0)
    if (associated(my_test)) then
      call mcmc_iteration(maxchain,upens,ens,multiplicity,f_my_jo,my_test=f_my_test)
    else
      call mcmc_iteration(maxchain,upens,ens,multiplicity,f_my_jo)
    endif
  case(1)
    if (associated(my_test)) then
      call mcmc_iteration(maxchain,upens,ens,multiplicity,f_my_jo,my_test=f_my_test,upxens=upxens,xens=xens)
    else
      call mcmc_iteration(maxchain,upens,ens,multiplicity,f_my_jo,upxens=upxens,xens=xens)
    endif
  end select
end subroutine

end module
