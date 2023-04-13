module ensscores

use iso_c_binding, only: c_int, c_double, c_funptr, c_f_procpointer

use ensdam_score_crps
use ensdam_score_rcrv
use ensdam_score_ranks
use ensdam_score_optimality
use ensdam_score_entropy

implicit none

interface
   function cdf_callback_interface(o,y,obs_idx) bind(c)
      import c_double, c_int
      real(c_double), intent(in) :: o
      real(c_double), intent(in) :: y
      integer(c_int), intent(in) :: obs_idx
      real(c_double) :: cdf_callback_interface
   end function
end interface
procedure (cdf_callback_interface), pointer, save :: cdf_obs => null()

interface
   subroutine events_callback_interface(nstate,nevents,member,outcome) bind(c)
      import c_double, c_int
      integer(c_int), intent(in), value :: nstate
      integer(c_int), intent(in), value :: nevents
      real(c_double), dimension(nstate), intent(in) :: member
      integer(c_int), dimension(nevents), intent(out) :: outcome
   end subroutine
end interface
procedure (events_callback_interface), pointer, save :: events_outcome => null()

contains

! Routines to exchange global variables from Fortran modules
subroutine c_get_crps_missing_value(var) bind(c)
   real(c_double), intent(out) :: var
   var = crps_missing_value
end subroutine

subroutine c_set_crps_missing_value(var) bind(c)
   real(c_double), intent(out) :: var
   crps_missing_value = var
end subroutine

subroutine c_get_rcrv_missing_value(var) bind(c)
   real(c_double), intent(out) :: var
   var = rcrv_missing_value
end subroutine

subroutine c_set_rcrv_missing_value(var) bind(c)
   real(c_double), intent(out) :: var
   rcrv_missing_value = var
end subroutine

subroutine c_get_rcrv_with_anamorphosis(var) bind(c)
   integer(c_int), intent(out) :: var
   var = 0
   if (rcrv_with_anamorphosis) var=1
end subroutine

subroutine c_set_rcrv_with_anamorphosis(var) bind(c)
   integer(c_int), intent(out) :: var
   rcrv_with_anamorphosis = (var/=0)
end subroutine

subroutine c_get_rcrv_number_of_quantiles(var) bind(c)
   integer(c_int), intent(out) :: var
   var = rcrv_number_of_quantiles
end subroutine

subroutine c_set_rcrv_number_of_quantiles(var) bind(c)
   integer(c_int), intent(out) :: var
   rcrv_number_of_quantiles = var
end subroutine

subroutine c_get_optimality_missing_value(var) bind(c)
   real(c_double), intent(out) :: var
   var = optimality_missing_value
end subroutine

subroutine c_set_optimality_missing_value(var) bind(c)
   real(c_double), intent(out) :: var
   optimality_missing_value = var
end subroutine

subroutine c_get_score_entropy_base(var) bind(c)
   real(c_double), intent(out) :: var
   var = score_entropy_base
end subroutine

subroutine c_set_score_entropy_base(var) bind(c)
   real(c_double), intent(out) :: var
   score_entropy_base = var
end subroutine

! Routines to associate callback functions

subroutine c_associate_cdf_callback(cdf_in) bind(c)
   type(c_funptr), intent(in), value :: cdf_in
   call c_f_procpointer(cdf_in, cdf_obs)
end subroutine

subroutine c_associate_events_callback(events_outcome_in) bind(c)
   type(c_funptr), intent(in), value :: events_outcome_in
   call c_f_procpointer(events_outcome_in, events_outcome)
end subroutine

! Fortran side of callback functions

function f_cdf_obs(o,y,obs_idx)
   real(kind=8), intent(in) :: o
   real(kind=8), intent(in) :: y
   integer, intent(in) :: obs_idx
   real(kind=8) :: f_cdf_obs

   f_cdf_obs = cdf_obs(o,y,obs_idx)
end function

subroutine f_events_outcome(member,outcome)
   real(kind=8), dimension(:), intent(in) :: member
   integer, dimension(:), intent(out) :: outcome
   integer :: nstate, nevents

   nstate = size(member,1)
   nevents = size(outcome,1)
   call events_outcome(nstate,nevents,member,outcome)
end subroutine

! Routines to compute ranks and rank histograms
subroutine c_compute_ranks(nstate, nens, ens, verif, ranks) bind(c)
   integer(c_int), intent(in), value :: nstate
   integer(c_int), intent(in), value :: nens
   real(c_double), intent(in) :: ens(nstate,nens)
   real(c_double), intent(in) :: verif(nstate)
   integer(c_int), intent(out) :: ranks(nstate)

   call compute_ranks(ens, verif, ranks)
end subroutine

subroutine c_compute_ranks_histogram(nstate, nens, ens, verif, ranks, rank_histogram) bind(c)
   integer(c_int), intent(in), value :: nstate
   integer(c_int), intent(in), value :: nens
   real(c_double), intent(in) :: ens(nstate,nens)
   real(c_double), intent(in) :: verif(nstate)
   integer(c_int), intent(out) :: ranks(nstate)
   integer(c_int), intent(out) :: rank_histogram(nens+1)

   call compute_ranks(ens, verif, ranks, rank_histogram)
end subroutine

! Routines to compute CRPS score
subroutine c_crps_score_global(nstate, nens, crps, reliability, resolution, ens, verif) bind(c)
   integer(c_int), intent(in), value :: nstate
   integer(c_int), intent(in), value :: nens
   real(c_double), intent(out) :: crps
   real(c_double), intent(out) :: reliability
   real(c_double), intent(out) :: resolution
   real(c_double), intent(in) :: ens(nstate,nens)
   real(c_double), intent(in) :: verif(nstate)

   call crps_score_global(crps, reliability, resolution, ens, verif)
end subroutine

subroutine c_crps_score_partition(nstate, nens, npartition, crps, reliability, resolution, ens, verif, partition) bind(c)
   integer(c_int), intent(in), value :: nstate
   integer(c_int), intent(in), value :: nens
   integer(c_int), intent(in), value :: npartition
   real(c_double), intent(out) :: crps(npartition)
   real(c_double), intent(out) :: reliability(npartition)
   real(c_double), intent(out) :: resolution(npartition)
   real(c_double), intent(in) :: ens(nstate,nens)
   real(c_double), intent(in) :: verif(nstate)
   integer(c_int), intent(in) :: partition(nstate)

   call crps_score_partition(crps, reliability, resolution, ens, verif, partition)
end subroutine

! Routines to compute RCRV score
subroutine c_rcrv_score_global(nstate, nens, ens_bias, ens_spread, ens, verif) bind(c)
   integer(c_int), intent(in), value :: nstate
   integer(c_int), intent(in), value :: nens
   real(c_double), intent(out) :: ens_bias
   real(c_double), intent(out) :: ens_spread
   real(c_double), intent(in) :: ens(nstate,nens)
   real(c_double), intent(in) :: verif(nstate)

   call rcrv_score_global(ens_bias, ens_spread, ens, verif)
end subroutine

subroutine c_rcrv_score_partition(nstate, nens, npartition, ens_bias, ens_spread, ens, verif, partition) bind(c)
   integer(c_int), intent(in), value :: nstate
   integer(c_int), intent(in), value :: nens
   integer(c_int), intent(in), value :: npartition
   real(c_double), intent(out) :: ens_bias(npartition)
   real(c_double), intent(out) :: ens_spread(npartition)
   real(c_double), intent(in) :: ens(nstate,nens)
   real(c_double), intent(in) :: verif(nstate)
   integer(c_int), intent(in) :: partition(nstate)

   call rcrv_score_partition(ens_bias, ens_spread, ens, verif, partition)
end subroutine

! Routines to compute OPTIMALITY score
subroutine c_optimality_score_global(nstate, nens, ens_optimality, ens, obs) bind(c)
   integer(c_int), intent(in), value :: nstate
   integer(c_int), intent(in), value :: nens
   real(c_double), intent(out) :: ens_optimality
   real(c_double), intent(in) :: ens(nstate,nens)
   real(c_double), intent(in) :: obs(nstate)

   if (associated(cdf_obs)) then
      call optimality_score_global(ens_optimality, ens, obs, f_cdf_obs)
   else
      stop 'callback function cdf_obs not associated in ensscores'
   endif
end subroutine

subroutine c_optimality_score_global_gaussian(nstate, nens, ens_optimality, ens, obs, std_obs) bind(c)
   integer(c_int), intent(in), value :: nstate
   integer(c_int), intent(in), value :: nens
   real(c_double), intent(out) :: ens_optimality
   real(c_double), intent(in) :: ens(nstate,nens)
   real(c_double), intent(in) :: obs(nstate)
   real(c_double), intent(in) :: std_obs(nstate)

   call optimality_score_global_gaussian(ens_optimality, ens, obs, std_obs)
end subroutine

subroutine c_optimality_score_partition(nstate, nens, npartition, ens_optimality, ens, obs, partition) bind(c)
   integer(c_int), intent(in), value :: nstate
   integer(c_int), intent(in), value :: nens
   integer(c_int), intent(in), value :: npartition
   real(c_double), intent(out) :: ens_optimality(npartition)
   real(c_double), intent(in) :: ens(nstate,nens)
   real(c_double), intent(in) :: obs(nstate)
   integer(c_int), intent(in) :: partition(nstate)

   if (associated(cdf_obs)) then
      call optimality_score_partition(ens_optimality, ens, obs, partition, f_cdf_obs)
   else
      stop 'callback function cdf_obs not associated in ensscores'
   endif
end subroutine

subroutine c_optimality_score_partition_gaussian(nstate, nens, npartition, ens_optimality, ens, obs, std_obs, partition) bind(c)
   integer(c_int), intent(in), value :: nstate
   integer(c_int), intent(in), value :: nens
   integer(c_int), intent(in), value :: npartition
   real(c_double), intent(out) :: ens_optimality(npartition)
   real(c_double), intent(in) :: ens(nstate,nens)
   real(c_double), intent(in) :: obs(nstate)
   real(c_double), intent(in) :: std_obs(nstate)
   integer(c_int), intent(in) :: partition(nstate)

   call optimality_score_partition_gaussian(ens_optimality, ens, obs, std_obs, partition)
end subroutine

! Routines to compute ENTROPY score
subroutine c_events_score(nstate, nens, nevents, noutcomes, score, ens, pref) bind(c)
   integer(c_int), intent(in), value :: nstate
   integer(c_int), intent(in), value :: nens
   integer(c_int), intent(in), value :: nevents
   integer(c_int), intent(in), value :: noutcomes
   real(c_double), intent(out) :: score(nevents)
   real(c_double), intent(in) :: ens(nstate,nens)
   real(c_double), intent(in) :: pref(nevents,noutcomes)
 
   if (associated(events_outcome)) then
      call events_score(score, ens, pref, f_events_outcome)
   else
      stop 'callback function events_outcome not associated in ensscores'
   endif
end subroutine

subroutine c_events_relative_entropy(nstate, nens, nevents, noutcomes, relative_entropy, ens, pref) bind(c)
   integer(c_int), intent(in), value :: nstate
   integer(c_int), intent(in), value :: nens
   integer(c_int), intent(in), value :: nevents
   integer(c_int), intent(in), value :: noutcomes
   real(c_double), intent(out) :: relative_entropy(nevents)
   real(c_double), intent(in) :: ens(nstate,nens)
   real(c_double), intent(in) :: pref(nevents,noutcomes)

   if (associated(events_outcome)) then
      call events_relative_entropy(relative_entropy, ens, pref, f_events_outcome)
   else
      stop 'callback function events_outcome not associated in ensscores'
   endif
end subroutine

subroutine c_events_cross_entropy(nstate, nens, nevents, noutcomes, cross_entropy, entropy, ens, pref) bind(c)
   integer(c_int), intent(in), value :: nstate
   integer(c_int), intent(in), value :: nens
   integer(c_int), intent(in), value :: nevents
   integer(c_int), intent(in), value :: noutcomes
   real(c_double), intent(out) :: cross_entropy(nevents)
   real(c_double), intent(out) :: entropy(nevents)
   real(c_double), intent(in) :: ens(nstate,nens)
   real(c_double), intent(in) :: pref(nevents,noutcomes)

   if (associated(events_outcome)) then
      call events_cross_entropy(cross_entropy, entropy, ens, pref, f_events_outcome)
   else
      stop 'callback function events_outcome not associated in ensscores'
   endif
end subroutine

subroutine c_events_entropy(nstate, nens, nevents, noutcomes, entropy, ens) bind(c)
   integer(c_int), intent(in), value :: nstate
   integer(c_int), intent(in), value :: nens
   integer(c_int), intent(in), value :: nevents
   integer(c_int), intent(in), value :: noutcomes
   real(c_double), intent(out) :: entropy(nevents)
   real(c_double), intent(in) :: ens(nstate,nens)

   if (associated(events_outcome)) then
      call events_entropy(entropy, ens, f_events_outcome, noutcomes)
   else
      stop 'callback function events_outcome not associated in ensscores'
   endif
end subroutine

subroutine c_events_probability(nstate, nens, nevents, noutcomes, pens, ens) bind(c)
   integer(c_int), intent(in), value :: nstate
   integer(c_int), intent(in), value :: nens
   integer(c_int), intent(in), value :: nevents
   integer(c_int), intent(in), value :: noutcomes
   real(c_double), intent(out) :: pens(nevents,noutcomes)
   real(c_double), intent(in) :: ens(nstate,nens)

   if (associated(events_outcome)) then
      call events_probability(pens, ens, f_events_outcome)
   else
      stop 'callback function events_outcome not associated in ensscores'
   endif
end subroutine

end module
