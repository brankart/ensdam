program scores_idealized_example
  use ensdam_score_crps
  use ensdam_score_rcrv
  use ensdam_score_entropy
  use ensdam_score_optimality
  use ensdam_storng
  use ensdam_obserror
  use ensdam_meanstd
  implicit none

  integer, parameter :: m=100            ! Size of the ensemble
  integer, parameter :: n=1000           ! Size of the state vector
  real(kind=8), parameter :: sigma=0.3   ! Observation error standard deviation

  real(kind=8), dimension(n,m) :: prior_ensemble
  real(kind=8), dimension(n,m) :: posterior_ensemble
  real(kind=8), dimension(n) :: reference_truth
  real(kind=8), dimension(n) :: observations

  real(kind=8) :: crps,crps_reliability,crps_resolution  ! CRPS score
  real(kind=8) :: rcrv_bias, rcrv_spread                 ! RCRV score
  real(kind=8), dimension(2) :: entropy_score            ! entropy score
  real(kind=8) :: optimality                             ! optimality score

  real(kind=8), dimension(n) :: ens_mean      ! mean of prior ensemble
  real(kind=8), dimension(n) :: ens_var       ! variance of prior ensemble
  real(kind=8), dimension(n) :: std_obs       ! observation error standard deviation
  real(kind=8), dimension(2,2) :: binary_pref ! Reference probability distribution for two binary events

  integer :: i, j    ! array indices

  integer :: nproc=1  ! Number of processors
  integer :: iproc=0  ! Current processor index

#if defined MPI
  ! Definition to make the example compatible with MPI version of the library
  include "mpif.h"
  integer, save :: mpi_code

  ! Initialize parallel computation
  call mpi_init(mpi_code)
  call mpi_comm_size(mpi_comm_world,nproc,mpi_code)
  call mpi_comm_rank(mpi_comm_world,iproc,mpi_code)
#endif

  ! Sample prior ensemble from N(0,I) distribution
  do j=1,m
  do i=1,n
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
  ! i) Generate perturbations to observations for each ensemble member
  do j=1,m
    do i=1,n
      call kiss_gaussian(posterior_ensemble(i,j))
    enddo
    posterior_ensemble(:,j) = sigma * posterior_ensemble(:,j) + observations(:)
  enddo
  ! ii) Compute innovation for each ensemble member
  posterior_ensemble(:,:) = posterior_ensemble(:,:) - prior_ensemble(:,:)
  ! iii) Compute mean and variance of prior ensemble
  call ensemble_meanstd(prior_ensemble,ens_mean,ens_var)
  ens_var = ens_var * ens_var
  ! iii) Multiply by Kalman gain
  do j=1,m
    posterior_ensemble(:,j) = posterior_ensemble(:,j) * ens_var(:) / ( ens_var(:) + sigma*sigma )
  enddo
  ! iv) Add prior ensemble
  posterior_ensemble(:,:) = posterior_ensemble(:,:) + prior_ensemble(:,:)

  ! Compute CRPS score, using reference truth as verification data
  call crps_score(crps,crps_reliability,crps_resolution,prior_ensemble,reference_truth)
  print '(a,2f8.5)', 'Prior CRPS reliability and resolution:    ',crps_reliability,crps_resolution

  call crps_score(crps,crps_reliability,crps_resolution,posterior_ensemble,reference_truth)
  print '(a,2f8.5)', 'Posterior CRPS reliability and resolution:',crps_reliability,crps_resolution

  ! Compute RCRV score, using reference truth as verification data
  call rcrv_score(rcrv_bias,rcrv_spread,prior_ensemble,reference_truth)
  print '(a,2f9.5)', 'Prior RCRV bias and spread:    ',rcrv_bias,rcrv_spread

  call rcrv_score(rcrv_bias,rcrv_spread,posterior_ensemble,reference_truth)
  print '(a,2f9.5)', 'Posterior RCRV bias and spread:',rcrv_bias,rcrv_spread

  ! Compute entropy score
  call events_probability(binary_pref,prior_ensemble,binary_event_outcomes)
  print '(a,2f6.3)', 'Prior probability distribution (event 1):    ',binary_pref(1,:)
  print '(a,2f6.3)', 'Prior probability distribution (event 2):    ',binary_pref(2,:)

  call events_probability(binary_pref,posterior_ensemble,binary_event_outcomes)
  print '(a,2f6.3)', 'Posterior probability distribution (event 1):',binary_pref(1,:)
  print '(a,2f6.3)', 'Posterior probability distribution (event 2):',binary_pref(2,:)

  call events_probability(binary_pref,prior_ensemble,binary_event_outcomes)
  call events_score(entropy_score,posterior_ensemble,binary_pref,binary_event_outcomes)
  print '(a,f6.3)', 'Entropy score (posterior vs prior, event 1):',entropy_score(1)
  print '(a,f6.3)', 'Entropy score (posterior vs prior, event 2):',entropy_score(2)

  ! Compute optimality score
  std_obs(:) = sigma

  call optimality_score(optimality,prior_ensemble,observations,std_obs)
  print '(a,f8.5)', 'Prior optimality score:    ',optimality

  call optimality_score(optimality,posterior_ensemble,observations,std_obs)
  print '(a,f8.5)', 'Posterior optimality score:',optimality

contains


  ! Callback routine to the outcome of the user defined events
  ! for a given ensemble member
  subroutine binary_event_outcomes(x,outcome)
  implicit none
  real(kind=8), dimension(:), intent(in) :: x    ! ensemble member
  integer, dimension(:), intent(out) :: outcome  ! events outcome

  if (sum(x*x)/(size(x,1)-1)<1.) then
    outcome(1) = 1
  else
    outcome(1) = 2
  endif

  if (maxval(abs(x))<3.4) then
    outcome(2) = 1
  else
    outcome(2) = 2
  endif

  end subroutine binary_event_outcomes


end
