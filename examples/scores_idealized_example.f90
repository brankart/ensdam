program scores_idealized_example
  use ensdam_score_crps
  use ensdam_score_rcrv
  use ensdam_score_entropy
  use ensdam_score_optimality
  use ensdam_obserror
  implicit none

  integer, parameter :: m=100           ! Size of the ensemble
  integer, parameter :: n=1000          ! Size of the state vector
  real(kind=8), parameter :: sigma=0.1  ! Observation error standard deviation

  real(kind=8), dimension(n,m) :: prior_ensemble
  real(kind=8), dimension(n,m) :: posterior_ensemble
  real(kind=8), dimension(n) :: reference_truth
  real(kind=8), dimension(n) :: observations

  real(kind=8) :: crps,crps_reliability,crps_resolution  ! CRPS score
  real(kind=8) :: rcrv_bias, rcrv_spread                 ! RCRV score
  real(kind=8), dimension(2) :: entropy_score            ! entropy score
  real(kind=8) :: optimality                             ! optimality score

  real(kind=8), dimension(2,2) :: binary_pref ! Reference probability distribution for two binary events

  integer :: i, j    ! array indices

  ! Sample prior ensemble from N(0,I) distribution
  do i=1,n
  do j=1,n
    call kiss_gaussian(prior_ensemble(i,j))
  enddo
  enddo

  ! Sample reference truth from the same distribution
  do i=1,n
    call kiss_gaussian(reference_truth(i))
  enddo

  ! Generate observations by adding N(0,sigma) perturbations to the reference truth
  do i=1,n
    call kiss_gaussian(observations(i))
  enddo

  observations(:) = reference_truth(:) + sigma * observations(:)

  ! Compute the posterior ensemble by conditioning the prior ensemble on observations
  posterior_ensemble = prior_ensemble

  ! Compute CRPS score, using reference truth as verification data
  call crps_score(crps,crps_reliability,crps_resolution,prior_ensemble,reference_truth)
  print *, 'Prior CRPS reliability:',crps_reliability
  print *, 'Prior CRPS resolution:',crps_resolution

  call crps_score(crps,crps_reliability,crps_resolution,posterior_ensemble,reference_truth)
  print *, 'Posterior CRPS reliability:',crps_reliability
  print *, 'Posterior CRPS resolution:',crps_resolution

  ! Compute RCRV score, using reference truth as verification data
  call rcrv_score(rcrv_bias,rcrv_spread,prior_ensemble,reference_truth)
  print *, 'Prior RCRV bias:',rcrv_bias
  print *, 'Prior RCRV spread:',rcrv_spread

  call rcrv_score(rcrv_bias,rcrv_spread,posterior_ensemble,reference_truth)
  print *, 'Posterior RCRV bias:',rcrv_bias
  print *, 'Posterior RCRV spread:',rcrv_spread

  ! Compute entropy score
  call events_score(entropy_score,prior_ensemble,binary_pref,binary_event_outcomes)
  print *, 'Prior entropy scores:',entropy_score

  call events_score(entropy_score,posterior_ensemble,binary_pref,binary_event_outcomes)
  print *, 'Posterior entropy scores:',entropy_score

  ! Compute optimality score
  call optimality_score(optimality,prior_ensemble,observations,cdf_obs)
  print *, 'Prior optimality score:',optimality

  call optimality_score(optimality,posterior_ensemble,observations,cdf_obs)
  print *, 'Posterior optimality score:',optimality

contains


  ! Callback routine to the outcome of the user defined events
  ! for a given ensemble member
  subroutine binary_event_outcomes(x,outcome)
  implicit none
  integer, dimension(:), intent(in) :: x    ! ensemble member
  integer(kind=8), dimension(2) :: outcome  ! events outcome

  if (maxval(abs(x))<1.) then
    outcome(1) = 1
  else
    outcome(1) = 2
  endif

  if (maxval(abs(x))<1.) then
    outcome(2) = 1
  else
    outcome(2) = 2
  endif

  end subroutine binary_event_outcomes


  ! Callback routine to compute observation error cdf:
  ! compute the rank of observation yo in p(yo|Hx), given y=Hx
  function cdf_obs(o,y,obs_idx)
  implicit none
  real(kind=8), intent(in) :: o   ! observation
  real(kind=8), intent(in) :: y   ! model equivalent to observation
  integer , intent(in) :: obs_idx ! index of observation in observation vector
  real(kind=8) :: cdf_obs

  obserror_type='gaussian'
  cdf_obs = obserror_cdf( o, y, sigma )

  end function cdf_obs

end
