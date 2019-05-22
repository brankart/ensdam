program mcmc_ensemble_update
  use ensdam_mcmc_update
  use ensdam_anaqua
  use ensdam_anatra
  use ensdam_anaobs
  use ensdam_obserror
  use ensdam_storng
  use ensdam_storfg
  use ensdam_stoanam
  use ensdam_stoutil
  use ensdam_sphylm
  use ensdam_ensaugm
  use ensdam_covariance
  use ensdam_score_crps
  use ensdam_score_optimality
  implicit none

  ! Include program parameters
  include 'mcmc_ensemble_update_parameters.h90'

  real(kind=8), parameter :: pi=3.14159265358979

  ! Size of vectors
  integer :: n ! size of state vector (for local subdomain if parallel computing)
  integer :: p ! size of observation vector (for local subdomain if parallel computing)
  integer :: ptot ! size of observation vector (for the whole problem)
  integer :: smax ! maximum degree of the spherical harmonics

  ! Main program arrays: grid, ensemble, observations
  real(kind=8), dimension(:), allocatable :: x, y ! coordinates of state variable
  real(kind=8), dimension(:), allocatable :: a ! area associated to each grid point
  integer, dimension(:,:), allocatable :: neigh ! neighbour indices
  real(kind=8), dimension(:,:,:), allocatable :: prior_ens ! prior ensemble
  real(kind=8), dimension(:), allocatable :: true_state ! synthetic true state
  real(kind=8), dimension(:), allocatable :: obs ! synthetic observations
  integer, dimension(:), allocatable :: hobs ! observation operator
  real(kind=8), dimension(:), allocatable :: obs_error_std ! obs error std after anamorphosis transformation
  real(kind=8), dimension(:,:), allocatable :: upobsens ! updated ensemble at observation locations
  real(kind=8), dimension(:,:), allocatable :: updated_ens ! updated ensemble
  real(kind=8), dimension(:,:), allocatable :: augmented_ens ! updated ensemble
  real(kind=8) :: obszero, obsmax, obsxmax, obsymax
  real(kind=8) :: tmpzero, tmpmax, tmpxmax, tmpymax
  logical :: nonlocal_observations, local_observations

  ! Arrays for anamorphosis transformation: quantiles, observation sample
  real(kind=8), dimension(q) :: quantiles_def ! quantiles definition
  real(kind=8), dimension(q) :: quantiles_ref ! quantiles of N(0,1)
  real(kind=8), dimension(:,:), allocatable :: quantiles_ens ! ensemble quantiles
  real(kind=8), dimension(:,:), allocatable :: obs_anam_sample ! sample of transformed observations

  ! Arrays for scale separation
  real(kind=8), dimension(:,:), allocatable :: spct ! spectrum in the basis of the spherical harmonics

  ! Temporary storage
  real(kind=8) :: uran, gran, amax, tmp ! temporary storage
  real(kind=8) :: stdmin, mean, std, misfit, misfit2 ! temporary storage
  real(kind=8) :: gamma_k, gamma_theta ! parameters of gamma distribution
  real(kind=8) :: obs_gamma_k, obs_gamma_theta ! parameters of observation error gamma distribution
  real(kind=8), dimension(:), allocatable :: tmp_state ! temporary storage
  real(kind=8), dimension(:), allocatable :: tmp_correl ! temporary storage
  integer, dimension(:), allocatable :: tmp_hobs ! temporary storage
  integer, dimension(:), allocatable :: tmp_sample ! temporary storage
  integer(KIND=8) :: zseed1, zseed2, zseed3, zseed4 ! seeds for random number generator
  character(len=256) :: outfile
  integer :: myunit

  ! Vector indices
  integer :: jm ! index of ensemble member
  integer :: jn ! index of state variable
  integer :: jp ! index of observation
  integer :: jq ! index of quantile
  integer :: js ! index of scale in multiple scale ensemble

  ! Parallel computation variables
  integer :: nproc=1  ! Number of processors
  integer :: iproc=0  ! Current processor index
  integer :: jiproc   ! Processor index

  ! Scores
  real(kind=8) :: score, score_reli, score_resol, score_prev
  integer :: iteration_index
  character(len=5) :: first_call='true'
  logical :: conv_test

  ! Diagnostic
  integer :: njo=0  ! number of evaluations of the cost function
  real(kind=8), dimension(:,:), allocatable :: ensref ! ensemble at reference point (to compute correlations)

#if defined MPI
  include "mpif.h"
  integer, save :: mpi_code
#endif

  nonlocal_observations=obs_max_location.or.obs_zero_surface

  ! Initialize parallel computation
#if defined MPI
  call mpi_init(mpi_code)
  call mpi_comm_size(mpi_comm_world,nproc,mpi_code)
  call mpi_comm_rank(mpi_comm_world,iproc,mpi_code)
#endif

  ! Seed random number generator
  call kiss_reset()
  do jiproc = 0, iproc + seedloop
    zseed1 = kiss() ; zseed2 = kiss() ; zseed3 = kiss() ; zseed4 = kiss()
  enddo
  call kiss_seed( zseed1, zseed2, zseed3, zseed4 )

  ! Problem definition
  ! ==================
  ! Definition of the grid locations (x,y) and
  ! domain decomposition for parallel computing
  if (iproc.eq.0) print *, 'Domain decomposition'
  call domain_decomposition()

  ! Generate the prior ensemble
  ! by sampling Gaussian random field with the specified spectrum 
  if (iproc.eq.0) print *, 'Generate prior ensemble'
  allocate(prior_ens(n,m,s))
  storfg_ylm_resolution=0.1  ! resolution for the precomputation of Legendre polynomials
  do jm = 1, m
    call gen_field_2s(prior_ens(:,jm,1),x,y,pow_spectrum,0,lmax)
  enddo

  ! Use additional random field as synthetic truth
  if (iproc.eq.0) print *, 'Generate synthetic truth'
  allocate(true_state(n))
  call gen_field_2s(true_state,x,y,pow_spectrum,0,lmax)

  ! Transform the marginal distributions to gamma distributions
  gamma_theta = gamma_nu * gamma_nu
  gamma_k = 1.0 / gamma_theta

  do jm = 1, m
  do jn = 1, n
    call gau_to_gam(tmp,prior_ens(jn,jm,1),gamma_k)
    tmp = tmp * gamma_theta
    tmp = tmp - gamma_shift
    if (tmp.lt.0.) tmp = 0.
    prior_ens(jn,jm,1) = tmp
  enddo
  enddo

  do jn = 1, n
    call gau_to_gam(tmp,true_state(jn),gamma_k)
    tmp = tmp * gamma_theta
    tmp = tmp - gamma_shift
    if (tmp.lt.0.) tmp = 0.
    true_state(jn) = tmp
  enddo

  ! Simulate local observations
  ! ===========================
  ! Generate observation operator by random sampling of state variables
  if (iproc.eq.0) print *, 'Generate observations'
  allocate(tmp_hobs(n))
  if (weighted_observation_sampling) then
    ! Compute maximum grid cell area
    amax = maxval(a)
#if defined MPI
    CALL MPI_ALLREDUCE (MPI_IN_PLACE, amax, 1, MPI_DOUBLE_PRECISION,  &
                  &     MPI_MAX,mpi_comm_world,mpi_code)
#endif
    ! Browse grid points and decide which are observed
    p=0
    do jn=1,n
      call kiss_uniform(uran)
      if (uran.lt.obs_coverage*a(jn)/amax) then
        ! this grid point is observed
        p = p + 1
        tmp_hobs(p) = jn
      endif
    enddo
    if (p.gt.0) then
      allocate(hobs(p))
      hobs(:) = tmp_hobs(1:p)
    endif
  else
    p = nint( obs_coverage * n )   ! number of observations
    if (p.gt.0) then
      p = max(1,p) ; p = min(p,n)
      allocate(hobs(p))
      if (p.eq.n) then
        hobs(:) = (/ (jp, jp=1,p) /)
      else
        tmp_hobs(:) = (/ (jn, jn=1,n) /)
        call kiss_sample(tmp_hobs,n,p)  ! observe random sample of state variables
        hobs(:) = tmp_hobs(1:p)
      endif
    endif
  endif
  deallocate(tmp_hobs)

  ! Compute total number of observations
  ptot = p
#if defined MPI
  CALL MPI_ALLREDUCE (MPI_IN_PLACE, ptot, 1, MPI_INTEGER,  &
                &     MPI_SUM,mpi_comm_world,mpi_code)
#endif
  if (iproc.eq.0) print *, 'Number of observations',ptot
  local_observations=ptot.gt.0
  if (local_observations) then
    if (p.eq.0) stop 'impossible option: each prcessor must contain at least one observation'
  else
    if (.not.nonlocal_observations) stop 'impossible option: there must be observations'
  endif

  if (local_observations) then
    ! Simulate synthetic observations from true state
    allocate(obs(p))
    obs(:) = true_state(hobs(:))

    ! Simulate observation errors
    if (iproc.eq.0) print *, 'Simulate observation errors'
    obserror_type='gamma'
    min_expected_gamma=obs_error_epsilon
    obs(:) = obserror_sample( obs(:), obs_error_ratio )
  endif

  ! Write problem in output files (in NetCDF)
  outfile=trim(expt)//'_prior_ensemble.nc'
  call write_ensemble(outfile,x,y,prior_ens(:,:,1))
  outfile=trim(expt)//'_true_state.nc'
  call write_state_vector(outfile,x,y,true_state(:))

  if (p.eq.n) then
    outfile=trim(expt)//'_observation.nc'
    call write_state_vector(outfile,x,y,obs(:))
    outfile=trim(expt)//'_obs_error.nc'
    call write_state_vector(outfile,x,y,obs(:)/true_state(:))
  else
    allocate(tmp_state(n))
    tmp_state = 0.
    if (local_observations) tmp_state(hobs(:)) = 1.
    outfile=trim(expt)//'_obs_operator.nc'
    call write_state_vector(outfile,x,y,tmp_state(:))
    deallocate(tmp_state)
  endif

  ! Simulate global observations
  ! ============================
  ! Observe maximum field location
  obsmax = get_max_location(tmpxmax,tmpymax,true_state)
  if (iproc.eq.0) print *, 'True maximum location',tmpxmax,tmpymax

  if (obs_max_location) then
    ! sample Gaussian angular perturbation
    call kiss_gaussian(gran) ; gran = obs_location_error_std * gran
    ! sample uniform azimuth
    call kiss_uniform(uran) ; uran = 2 * pi * uran
    ! compute perturbed coordinates
    tmpxmax = tmpxmax * pi / 180. ; tmpymax = tmpymax * pi / 180.
    obsymax = asin( sin(tmpymax)*cos(gran) + cos(tmpymax)*sin(gran)*cos(uran) )
    if (cos(tmpymax).eq.0.) then
      ! reference is at one pole
      obsxmax = uran
    elseif (cos(obsymax).eq.0.) then
      ! observation is at one pole
      obsxmax = tmpxmax
    else
      ! both are elsewhere
      obsxmax = tmpxmax + acos( (cos(gran) - sin(tmpymax)*sin(obsymax)) / ( cos(tmpymax)*cos(obsymax) ) )
    endif
    tmpxmax = tmpxmax / pi * 180. ; tmpymax = tmpymax / pi * 180.
    obsxmax = obsxmax / pi * 180. ; obsymax = obsymax / pi * 180.
#if defined MPI
    CALL mpi_bcast(obsxmax,1,mpi_double_precision,0,mpi_comm_world,mpi_code)
    CALL mpi_bcast(obsymax,1,mpi_double_precision,0,mpi_comm_world,mpi_code)
#endif
    if (iproc.eq.0) print *, 'Obs. maximum location',obsxmax,obsymax
    obsxmax = obsxmax * pi / 180. ; obsymax = obsymax * pi / 180.
  endif

  ! Observe surface where the field is equal to zero
  obszero = get_zero_surface(true_state)
  if (iproc.eq.0) print *, 'True zero surface',obszero

  if (obs_zero_surface) then
    obserror_type='beta'
    min_expected_beta=obs_error_epsilon
    obszero = obserror_sample( obszero, obs_surface_error_std )
#if defined MPI
    CALL mpi_bcast(obszero,1,mpi_double_precision,0,mpi_comm_world,mpi_code)
#endif
    if (iproc.eq.0) print *, 'Obs. zero surface',obszero
  endif

  ! Write ensemble equivalent to global observations
  outfile=trim(expt)//'_prior_ensemble_max_location.txt'
  open(newunit=myunit,file=outfile)
  do jm = 1, m
    tmpmax = get_max_location(tmpxmax,tmpymax,prior_ens(:,jm,1))
    write(myunit,'(i,f,f)') jm,tmpxmax,tmpymax
  enddo
  close(myunit)

  outfile=trim(expt)//'_prior_ensemble_zero_surface.txt'
  open(newunit=myunit,file=outfile)
  do jm = 1, m
    tmpzero = get_zero_surface(prior_ens(:,jm,1))
    tmp = get_true_zero_surface(prior_ens(:,jm,1))
    write(myunit,'(i,2f)') jm,tmpzero,tmp
  enddo
  close(myunit)

  ! Forward anamorphosis transformation
  ! ===================================
  ! Define ensemble quantiles and compute
  ! quantiles of reference distribution [ N(0,1) ]
  if (iproc.eq.0) print *, 'Computation of the quantiles of the prior ensemble'
  do jq=1,q
    quantiles_def(jq) = real(jq-1,8)/(q-1)
    tmp = ( 1. + quantiles_def(jq) * (m-1) ) / (m+1)
    quantiles_ref(jq) = invcdf_gaussian(tmp)
  enddo

  ! Compute ensemble quantiles
  allocate(quantiles_ens(n,q))
  call ens_quantiles(quantiles_ens, prior_ens(:,:,1), quantiles_def )

  ! Option to transform observations (!! not done by default !!)
  if ((observation_anamorphosis).and.(p.gt.0)) then

    ! Transform observations
    if (iproc.eq.0) print *, 'Forward anamorphosis of observations'
    allocate(obs_anam_sample(p,obs_sample_size))
    allocate(obs_error_std(p))
    obs_error_std(:) = obs_error_ratio
    obserror_type='gamma'
    call ana_obs( obs_anam_sample, prior_ens(hobs(:),:,1), obs, obs_error_std, quantiles_def, quantiles_ref )

    ! Write transformed observation sample
    if (p.eq.n) then
      outfile=trim(expt)//'_observation_sample_ana.nc'
      call write_ensemble(outfile,x,y,obs_anam_sample(:,:))
    endif

    ! Compute observation and observation error as mean and std of the sample
    do jp=1,p
      std = 0. ; mean = 0. ! temporary storage to store ensemble mean
      do js=1,obs_sample_size
        misfit = obs_anam_sample(jp,js) - mean
        mean = mean + misfit / js
        misfit2 = obs_anam_sample(jp,js) - mean
        std = std + misfit * misfit2
      enddo
      std = SQRT( std / (m-1) )
      obs(jp) = mean
      obs_error_std(jp) = std
    enddo

    deallocate(obs_anam_sample)

    ! Regularize observation error standard deviation
    stdmin = maxval(obs_error_std) * obs_error_min
    where(obs_error_std.lt.stdmin) obs_error_std = stdmin

  endif

  ! Transform prior ensemble
  if (iproc.eq.0) print *, 'Forward anamorphosis of prior ensemble'
  call ana_forward( prior_ens(:,:,1), quantiles_ens, quantiles_ref )

  ! Transform and diagnose true state
  allocate(tmp_state(n))
  tmp_state = true_state
  call ana_forward( tmp_state, quantiles_ens, quantiles_ref )

  outfile=trim(expt)//'_true_state_ana.nc'
  call write_state_vector(outfile,x,y,true_state(:))

  if (nonlocal_observations) then
    tmp = cost_jo_state(tmp_state(:))
  else
    tmp = cost_jo(tmp_state(hobs(:)))
  endif

  if (iproc.eq.0) print *, 'Jo for true state:',tmp

  ! Construct multiple scale ensemble
  ! =================================
  if (iproc.eq.0) print *, 'Construct multiple scale ensemble'
  smax = maxval(scales(:))  ! maximum degree of the spherical harmonics
  allocate(spct(0:smax,-smax:smax))

  external_vector_decomposition=.TRUE. ! parallelization is ruled from outside sphylm
  tmp = min( 360./nlon, 180./nlat ) / plegres  ! resolution required for Legendre polynomials
  call init_ylm( smax, 0, -90._8, 90._8, tmp )

  do jm=1,m  ! loop on ensemble members

    tmp_state(:) = prior_ens(:,jm,1) * a(:)
    call proj_ylm( spct, tmp_state, x, y )
#if defined MPI
    call mpi_allreduce (MPI_IN_PLACE, spct, (smax+1)*(2*smax+1),  &
     &                  MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,mpi_code)
#endif

    do js=2,s ! loop on additional scales to include in the prior ensemble
      prior_ens(:,jm,js) = 0.
      call back_ylm_loc( spct, prior_ens(:,jm,js), x, y, 0, scales(js) )
    enddo

  enddo

  ! Renormalize std of all additional scales to 1
  do js=2,s
    do jn=1,n
      std = 0. ; mean = 0. ! temporary storage to store ensemble mean
      do jm=1,m
        misfit = prior_ens(jn,jm,js) - mean
        mean = mean + misfit / jm
        misfit2 = prior_ens(jn,jm,js) - mean
        std = std + misfit * misfit2
      enddo
      std = SQRT( std / (m-1) )
      prior_ens(jn,:,js) = prior_ens(jn,:,js) / std
    enddo
  enddo

  ! Write files (in NetCDF)
  outfile=trim(expt)//'_prior_ensemble_ana.nc'
  call write_ensemble(outfile,x,y,prior_ens(:,:,1))
  outfile=trim(expt)//'_prior_ensemble_ana_js=2.nc'
  call write_ensemble(outfile,x,y,prior_ens(:,:,2))
  outfile=trim(expt)//'_ensemble_quantiles.nc'
  call write_ensemble(outfile,x,y,quantiles_ens(:,:))


  if ((observation_anamorphosis).and.(p.eq.n)) then
    outfile=trim(expt)//'_observation_ana.nc'
    call write_state_vector(outfile,x,y,obs(:))
    outfile=trim(expt)//'_obs_error_std_ana.nc'
    call write_state_vector(outfile,x,y,obs_error_std(:))
  endif

  ! Diagnostic of correlation structures
  ! ====================================

  if (sample_aug_ensemble) then

    ! Check shape of multiple product
    allocate( tmp_sample(SUM(multiplicity(:))) )
    if (iproc.eq.0) print *, 'multiplicity:',multiplicity
    CALL newproduct( tmp_state, prior_ens, multiplicity, tmp_sample )
    deallocate( tmp_sample )

    ! Write files (in NetCDF)
    outfile=trim(expt)//'_multiple_product.nc'
    call write_state_vector(outfile,x,y,tmp_state(:))

    ! Sample augmented ensemble
    allocate(augmented_ens(n,maug))
    if (iproc.eq.0) print *, 'Sample augmented ensemble'
    call sample_augmented_ensemble( augm_maxiter, augmented_ens, prior_ens(:,:,:), multiplicity )
    if (iproc.eq.0) print *, 'Write augmented ensemble'
    outfile=trim(expt)//'_augmented_ensemble_ana.nc'
    call write_ensemble(outfile,x,y,augmented_ens(:,:))

  endif

  if (diagnose_correlation) then

    ! Prepare to compute correlation structure
    allocate(tmp_correl(n))
    tmp_correl = 1.

    ! Define ensemble at reference point
    allocate(ensref(m,s))
    ensref =  prior_ens(n/2,:,:)

#if defined MPI
    CALL mpi_bcast(ensref,m*s,mpi_double_precision,0,mpi_comm_world,mpi_code)
#endif

    ! Diagnose correlation structure for each scale
    do js=2,s

      CALL ensemble_correlation( prior_ens(:,:,js), ensref(:,js), tmp_state )
  
      tmp_correl = tmp_correl * ( tmp_state ** multiplicity(js) )

      if (js.EQ.2) then
        outfile=trim(expt)//'_prior_correlation_js=2.nc'
        call write_state_vector(outfile,x,y,tmp_state(:))
      endif

    enddo

    outfile=trim(expt)//'_localizing_correlation.nc'
    call write_state_vector(outfile,x,y,tmp_correl(:))

    CALL ensemble_correlation( prior_ens(:,:,1), ensref(:,1), tmp_state )
    outfile=trim(expt)//'_prior_correlation.nc'
    call write_state_vector(outfile,x,y,tmp_state(:))

    tmp_correl = tmp_correl * tmp_state

    outfile=trim(expt)//'_localized_correlation.nc'
    call write_state_vector(outfile,x,y,tmp_correl(:))

    deallocate(ensref)

  endif

  if (sample_aug_ensemble.and.diagnose_correlation) then

    ! Define augmented ensemble at reference point
    allocate(ensref(maug,s))
    ensref(:,1) =  augmented_ens(n/2,:)

#if defined MPI
    CALL mpi_bcast(ensref,maug*s,mpi_double_precision,0,mpi_comm_world,mpi_code)
#endif

    ! Diagnose correlation structure of augmented ensemble
    CALL ensemble_correlation( augmented_ens(:,:), ensref(:,1), tmp_correl )

    outfile=trim(expt)//'_augmented_correlation.nc'
    call write_state_vector(outfile,x,y,tmp_correl(:))

  endif

  if (sample_aug_ensemble) then
    call ana_backward( augmented_ens(:,:), quantiles_ens, quantiles_ref )

    outfile=trim(expt)//'_augmented_ensemble.nc'
    call write_ensemble(outfile,x,y,augmented_ens(:,:))

    deallocate(augmented_ens)
  endif

  if (diagnose_correlation) deallocate(tmp_correl)

  ! Update ensemble using MCMC algorithm
  ! ====================================
  if (iproc.eq.0) print *, 'Update ensemble using MCMC algorithm'

  if (perform_observational_update.or.compute_scores) then

    ! Initialize scores and convergence test
    if (nonlocal_observations) then
      conv_test = convergence_test(prior_ens(:,:,1),prior_ens(:,:,1))
    else
      conv_test = convergence_test(prior_ens(hobs(:),:,1),prior_ens(:,:,1))
    endif

  endif

  if (perform_observational_update) then

    ! Set number of iterations between convergence checks
    mcmc_convergence_check = convergence_check

    ! Call to MCMC iteration to perform observational update
    allocate(updated_ens(n,mup))
    if (nonlocal_observations) then
      call mcmc_iteration( mcmc_maxiter, updated_ens, prior_ens, multiplicity, cost_jo_state, convergence_test )
    else
      allocate(upobsens(p,mup))
      call mcmc_iteration( mcmc_maxiter, upobsens, prior_ens(hobs(:),:,:), multiplicity, &
                        &  cost_jo, my_test=convergence_test, upxens=updated_ens, xens=prior_ens(:,:,:) )
      deallocate(upobsens)
    endif

    ! Diagnostics of the update
    if (iproc.eq.0) print *, 'Number of evaluations of the cost function:', njo

    if (laplacian_penalty) then
      tmp = lapl_penalty(updated_ens(:,1))
      if (iproc.eq.0) print *, 'Laplacian penalty for updated state:',tmp
    endif

    ! Write files (in NetCDF)
    outfile=trim(expt)//'_updated_ensemble_ana.nc'
    call write_ensemble(outfile,x,y,updated_ens(:,:))

    ! Backward anamorphosis transformation of updated ensemble
    if (iproc.eq.0) print *, 'Backward anamorphosis of updated ensemble'
    call ana_backward( updated_ens(:,:), quantiles_ens, quantiles_ref )

    ! Write files (in NetCDF)
    outfile=trim(expt)//'_updated_ensemble.nc'
    call write_ensemble(outfile,x,y,updated_ens(:,:))

    ! Write ensemble equivalent to global observations
    outfile=trim(expt)//'_updated_ensemble_max_location.txt'
    open(newunit=myunit,file=outfile)
    do jm = 1, mup
      tmpmax = get_max_location(tmpxmax,tmpymax,updated_ens(:,jm))
      write(myunit,'(i,f,f)') jm,tmpxmax,tmpymax
    enddo
    close(myunit)

    outfile=trim(expt)//'_updated_ensemble_zero_surface.txt'
    open(newunit=myunit,file=outfile)
    do jm = 1, mup
      tmpzero = get_zero_surface(updated_ens(:,jm))
      tmp = get_true_zero_surface(updated_ens(:,jm))
      write(myunit,'(i,2f)') jm,tmpzero,tmp
    enddo
    close(myunit)

  endif

contains

  ! ======================================================================

  ! Callback routine to check the convergence of iterations
  ! Intermediate diagnostics can also be computed here
  function convergence_test(upens,upxens)
  use ensdam_score_crps
  implicit none
  real(kind=8), dimension(:,:), intent(in) :: upens
  real(kind=8), dimension(:,:), intent(in), optional :: upxens
  logical :: convergence_test

  real(kind=8), dimension(:,:), allocatable :: tmpens
  real(kind=8) :: tmpdist,tmpzero, cost_zero, zsqr_zero, zsqr_max

  convergence_test = .FALSE.

  ! Initialize previous score at first call
  if (first_call/='false') then
    iteration_index = 0
    score_prev=tiny(score_prev)
    first_call='false'
  else
    iteration_index = iteration_index + convergence_check
  endif

  ! Compute optimality score
  if ((check_convergence).and.(p.ne.0)) then
    if (nonlocal_observations) then
      call optimality_score(score, upens(hobs(:),:), obs, cdf_obs)
    else
      call optimality_score(score, upens, obs, cdf_obs)
    endif
    if (iproc.eq.0) print *, 'OPTIMALITY:',iteration_index,score
    convergence_test = abs( score - score_prev ) / score_prev < convergence_epsilon
    score_prev=score
    if (convergence_test) then
      if (iproc.eq.0) print *, 'Convergence reached !'
    endif
  endif

  ! Compute other scores
  if (compute_scores) then
    allocate(tmpens(n,mup))

    if (nonlocal_observations) then
      tmpens(:,:) = upens(:,:)
    else
      tmpens(:,:) = upxens(:,:)
    endif

    call ana_backward( tmpens(:,:), quantiles_ens, quantiles_ref )
    call crps_score(score, score_reli, score_resol, tmpens, true_state)

    if (iproc.eq.0) print *, 'CRPS:',iteration_index,score_reli,score_resol

    ! DIagnose cost function associated to the nonlocal observations
    if (nonlocal_observations) then

      cost_zero = 0. ; zsqr_zero = 0. ; zsqr_max = 0.
      do jm = 1, m

        if (obs_max_location) then
          tmpmax = get_max_location(tmpxmax,tmpymax,tmpens(:,jm))
          tmpxmax = tmpxmax * pi / 180. ; tmpymax = tmpymax * pi / 180.
          tmpdist = ACOS ( MIN( SIN(tmpymax)*SIN(obsymax)*COS(tmpxmax-obsxmax) + COS(tmpymax)*COS(obsymax) , 1.) )
          tmpdist = tmpdist / obs_location_error_std
          zsqr_max = zsqr_max + tmpdist * tmpdist / 2.
        endif
        if (obs_zero_surface) then
          tmpzero = get_zero_surface(tmpens(:,jm))
          obserror_type='beta'
          min_expected_beta=obs_error_epsilon
          cost_zero = cost_zero + obserror_logpdf( obszero, tmpzero, obs_surface_error_std )
          tmpzero = obserror_cdf( obszero, tmpzero, obs_surface_error_std )
          tmpzero = invcdf_gaussian( tmpzero )
          zsqr_zero = zsqr_zero + tmpzero * tmpzero
        endif

      enddo

      cost_zero = cost_zero / m ; zsqr_zero = zsqr_zero / m ; zsqr_max = zsqr_max / m
      if (obs_zero_surface) then
        if (iproc.eq.0) print *, 'Cost zero surface:',cost_zero
        if (iproc.eq.0) print *, 'Zsqr zero surface:',zsqr_zero
      endif
      if (obs_max_location) then
        if (iproc.eq.0) print *, 'Zsqr max location:',zsqr_max
      endif

    endif

    deallocate(tmpens)
  endif

  end function

  ! ======================================================================

  ! Callback routine to compute observation error cdf:
  ! compute the rank of observation yo in p(yo|Hx), given y=Hx
  function cdf_obs(o,y,obs_idx)
  implicit none
  real(kind=8), intent(in) :: o
  real(kind=8), intent(in) :: y
  integer , intent(in) :: obs_idx
  real(kind=8) :: cdf_obs

  real(kind=8) :: ytmp, uran

  if (observation_anamorphosis) then
    obserror_type='gaussian'
    cdf_obs = obserror_cdf( o, y, obs_error_std(obs_idx) )
  else
    ytmp = y
    call ana_backward( ytmp, quantiles_ens(hobs(obs_idx),:), quantiles_ref )
    obserror_type='gamma'
    cdf_obs = obserror_cdf( o, ytmp, obs_error_ratio )
  endif

  end function cdf_obs
  ! ======================================================================

  ! Callnack routine to compute observation cost function
  ! Jo = -log p(yo|hx), computed using observation equivalent as argument
  ! This version is used if there are only local observations
  function cost_jo(hx)
  implicit none
  real(kind=8), dimension(:), intent(in) :: hx
  real(kind=8) :: cost_jo

  real(kind=8), dimension(:), allocatable :: tmpstate

  if (observation_anamorphosis) then
    obserror_type='gaussian'
    cost_jo = obserror_logpdf( obs(1:p), hx(1:p), obs_error_std(1:p) )
  else
    allocate(tmpstate(p))
    tmpstate = hx(:)
    call ana_backward( tmpstate(:), quantiles_ens(hobs(:),:), quantiles_ref )
    obserror_type='gamma'
    cost_jo = obserror_logpdf( obs(1:p), tmpstate(1:p), obs_error_ratio )
    deallocate(tmpstate)
  endif

#if defined MPI
  call MPI_ALLREDUCE (MPI_IN_PLACE, cost_jo, 1, MPI_DOUBLE_PRECISION,  &
                        &     MPI_SUM,mpi_comm_world,mpi_code)
#endif

  njo = njo + 1

  end function cost_jo

  ! ======================================================================

  ! Callnack routine to compute observation cost function
  ! Jo = -log p(yo|hx), computed using state vector as argument
  ! This version is used if there are global observations
  function cost_jo_state(state)
  implicit none
  real(kind=8), dimension(:), intent(in) :: state
  real(kind=8) :: cost_jo_state
  
  real(kind=8), dimension(:), allocatable :: tmpstate
  real(kind=8) :: tmpzero, tmpdist

  cost_jo_state = 0.

  if (local_observations) then
    if (observation_anamorphosis) then
      obserror_type='gaussian'
      cost_jo_state = obserror_logpdf( obs(1:p), state(hobs(1:p)), obs_error_std(1:p) )
    else
      allocate(tmpstate(p))
      tmpstate = state(hobs(:))
      call ana_backward( tmpstate(:), quantiles_ens(hobs(:),:), quantiles_ref )
      obserror_type='gamma'
      cost_jo_state = obserror_logpdf( obs(1:p), tmpstate(1:p), obs_error_ratio )
      deallocate(tmpstate)
    endif
  endif

#if defined MPI
  call MPI_ALLREDUCE (MPI_IN_PLACE, cost_jo_state, 1, MPI_DOUBLE_PRECISION,  &
                        &     MPI_SUM,mpi_comm_world,mpi_code)
#endif

  if (nonlocal_observations) then
    allocate(tmpstate(n))
    tmpstate = state(:)
    call ana_backward( tmpstate(:), quantiles_ens(:,:), quantiles_ref )
    if (obs_max_location) then
      tmpmax = get_max_location(tmpxmax,tmpymax,tmpstate)
      tmpxmax = tmpxmax * pi / 180. ; tmpymax = tmpymax * pi / 180.
      tmpdist = acos ( MIN( cos(tmpymax)*cos(obsymax)*cos(tmpxmax-obsxmax) + sin(tmpymax)*sin(obsymax) , 1.) )
      tmpdist = tmpdist / obs_location_error_std
      cost_jo_state = cost_jo_state + tmpdist * tmpdist / 2.
    endif
    if (obs_zero_surface) then
      obserror_type='beta'
      min_expected_beta=obs_error_epsilon
      tmpzero = get_zero_surface(tmpstate)
      cost_jo_state = cost_jo_state + obserror_logpdf( obszero, tmpzero, obs_surface_error_std )
    endif
    deallocate(tmpstate)
  endif

  if (laplacian_penalty) cost_jo_state = cost_jo_state + lapl_penalty(state)

  njo = njo + 1

  end function cost_jo_state

  ! ======================================================================

  ! Get maximum and maximum location of input field
  function get_max_location(xmax,ymax,state)
  implicit none
  real(kind=8), dimension(:), intent(in) :: state
  real(kind=8), intent(out) :: xmax,ymax
  real(kind=8) :: get_max_location

  real(kind=8) :: vmax
  real(kind=8), dimension(2) :: tmp_max
  integer :: int_type

  vmax = state(transfer(maxloc(state),int_type))
  xmax = x(transfer(maxloc(state),int_type))
  ymax = y(transfer(maxloc(state),int_type))

#if defined MPI
  tmp_max(1) = vmax
  tmp_max(2) = xmax
  CALL MPI_ALLREDUCE (MPI_IN_PLACE, tmp_max, 1, MPI_2DOUBLE_PRECISION,  &
                &     MPI_MAXLOC,mpi_comm_world,mpi_code)
  xmax = tmp_max(2)

  tmp_max(1) = vmax
  tmp_max(2) = ymax
  CALL MPI_ALLREDUCE (MPI_IN_PLACE, tmp_max, 1, MPI_2DOUBLE_PRECISION,  &
                &     MPI_MAXLOC,mpi_comm_world,mpi_code)
  ymax = tmp_max(2)

  vmax = tmp_max(1)
#endif

  get_max_location = vmax

  end function get_max_location

  ! ======================================================================

  ! Get surface where the input field is equal to zero
  function get_zero_surface(state)
  implicit none
  real(kind=8), dimension(:), intent(in) :: state
  real(kind=8) :: get_zero_surface

  real(kind=8) :: szero

  szero = sum(a,mask=(state.eq.0.))

#if defined MPI
  CALL MPI_ALLREDUCE (MPI_IN_PLACE, szero, 1, MPI_DOUBLE_PRECISION,  &
                &     MPI_SUM,mpi_comm_world,mpi_code)
#endif

  get_zero_surface = min(max(szero,0.),1.)

  end function get_zero_surface

  ! ======================================================================

  ! Get surface where both the input field and the true state are equal to zero
  function get_true_zero_surface(state)
  implicit none
  real(kind=8), dimension(:), intent(in) :: state
  real(kind=8) :: get_true_zero_surface

  real(kind=8) :: szero

  szero = sum(a,mask=((state.eq.0.).and.(true_state.eq.0.)))

#if defined MPI
  CALL MPI_ALLREDUCE (MPI_IN_PLACE, szero, 1, MPI_DOUBLE_PRECISION,  &
                &     MPI_SUM,mpi_comm_world,mpi_code)
#endif

  get_true_zero_surface = min(max(szero,0.),1.)

  end function get_true_zero_surface

  ! ======================================================================

  ! Laplacian penalty on state vector
  function lapl_penalty(state)
  implicit none
  real(kind=8), dimension(:), intent(in) :: state
  real(kind=8) :: lapl_penalty

  real(kind=8), dimension(:), allocatable :: fullstate, lat
  real(kind=8), dimension(:,:), allocatable :: nei
  real(kind=8) :: lapl, laplon, laplat, factor
  integer :: i, k

  ! allocate full state
  allocate(fullstate(nlon*nlat),lat(nlon*nlat),nei(nlon*nlat,4))
  fullstate = 0.

  ! merge processors contributions
  k=0
  do i=1,nlon*nlat
    if (mod(k-iproc,nproc).eq.0) then
      lat(i) = y(1+k/nproc) * pi / 180.
      nei(i,:) = neigh(1+k/nproc,:)
      fullstate(i) = state(1+k/nproc)
    endif
    k = k + 1
  enddo

#if defined MPI
  call MPI_ALLREDUCE (MPI_IN_PLACE,fullstate, nlon*nlat,MPI_DOUBLE_PRECISION, &
     &                MPI_SUM,mpi_comm_world,mpi_code)
#endif

  ! compute Laplacian penalty
  lapl_penalty = 0.
  do i=1+iproc,nlon*nlat,nproc
    if (all(nei(i,:).ne.0)) then
      ! longitude contribution
      factor = cos(lat(i)) ; factor = 1. / factor
      laplon = fullstate(nei(i,1)) + fullstate(nei(i,2)) - 2 * fullstate(i)
      laplon = laplon * factor * factor
      ! latitude contribution
      laplat = fullstate(nei(i,4)) * ( cos(lat(i)) + cos(lat(nei(i,4))) ) &
           & + fullstate(nei(i,3)) * ( cos(lat(i)) + cos(lat(nei(i,3))) ) &
           & - fullstate(i) * ( 2*cos(lat(i)) + cos(lat(nei(i,4))) + cos(lat(nei(i,3))) )
      laplat = 0.5 * laplat * factor
      ! total local Laplacian
      lapl = laplon + laplat
      ! add to penalty
      lapl_penalty = lapl_penalty + lapl * lapl / factor
    endif
  enddo

#if defined MPI
  call MPI_ALLREDUCE (MPI_IN_PLACE,lapl_penalty,1,MPI_DOUBLE_PRECISION, &
     &                MPI_SUM,mpi_comm_world,mpi_code)
#endif

  deallocate(fullstate,lat,nei)

  end function lapl_penalty

  ! ======================================================================

  ! Power spectrum in the basis of the spherical harmonics
  ! (as a function of the degree l)
  function pow_spectrum(l,m)
  implicit none
  integer, intent(in) :: l,m
  real(kind=8) :: pow_spectrum

  integer :: ll, mm
  real(kind=8) :: norm, rm

  ! Power spectrum
  pow_spectrum = 1. / ( 1. + (l-lm)*(l-lm)/(lc*lc) )

  ! Normalize the spectrum
  norm = 0.
  do ll=0,lmax
    norm = norm + 1. / ( 1. + (ll-lm)*(ll-lm)/(lc*lc) )
  enddo
  pow_spectrum = pow_spectrum / norm

  if (anisotropy.eq.0.) then

    ! Scale to account for the multiplicity of each degree
    pow_spectrum = pow_spectrum / ( 1. + 2. * l )

  else

    ! Modify power spectrum
    rm = real(abs(m),8)/max(l,1)
    pow_spectrum = pow_spectrum * ( 1. - rm ) ** anisotropy
    !pow_spectrum = pow_spectrum * exp ( - anisotropy * rm )

    ! Normalize the spectrum
    norm = 0.
    do mm=-l,l
      rm = real(abs(mm),8)/max(l,1)
      norm = norm + ( 1. - rm ) ** anisotropy
      !norm = norm + exp ( - anisotropy * rm )
    enddo
    pow_spectrum = pow_spectrum / norm

  endif

  end function pow_spectrum

  ! ======================================================================

  ! Domain decomposition for parallel computing
  subroutine domain_decomposition()
  use ensdam_spharea
  implicit none

  real(kind=8), dimension(:,:), allocatable :: lon, lat, area
  integer, dimension(:,:,:), allocatable :: neighbours
  real(kind=8) :: dlon, dlat
  integer :: i, j, k

  ! allocate grid on the full domain
  allocate(lon(0:nlon,0:nlat),lat(0:nlon,0:nlat),area(nlon,nlat))
  if (laplacian_penalty) then
    allocate(neighbours(0:nlon+1,0:nlat+1,4))
    neighbours = 0
  endif

  ! define the output grid
  dlon = 360. / nlon
  dlat = 180. / nlat

  ! location of mesh vertices
  do i=0,nlon
    lon(i,0) = i * dlon
  enddo
  do j=0,nlat
    lat(0,j) = j * dlat - 90.
  enddo
  do i=0,nlon
    lon(i,1:nlat) = lon(i,0)
  enddo
  do j=0,nlat
    lat(1:nlon,j) = lat(0,j)
  enddo

  ! compute area associated to each grid point
  call mesh_area (lon,lat,area)

  ! location of state variables (0 indices are not used anymore)
  lon(:,:) = lon(:,:) - dlon / 2
  lat(:,:) = lat(:,:) - dlat / 2

  ! determine size of subdomain for current processor
  n = nlon * nlat / nproc
  if (iproc.lt.mod(nlon*nlat,nproc)) n = n + 1

  ! allocate local grid arrays
  allocate(x(n),y(n),a(n))

  ! provide a subdomain to each processor
  k=0
  do j=1,nlat
  do i=1,nlon
    if (mod(k-iproc,nproc).eq.0) then
      x(1+k/nproc) = lon(i,j)
      y(1+k/nproc) = lat(i,j)
      a(1+k/nproc) = area(i,j)
      if (laplacian_penalty) then
        ! define global neigbour array
        neighbours(i+1,j,1) = k+1
        neighbours(i-1,j,2) = k+1
        neighbours(i,j+1,3) = k+1
        neighbours(i,j-1,4) = k+1
      endif
    endif
    k = k + 1
  enddo
  enddo

  ! deallocate global grid
  deallocate(lon,lat)

  ! initialize local neighbour arrays
  if (laplacian_penalty) then

#if defined MPI
    call MPI_ALLREDUCE (MPI_IN_PLACE,neighbours, (nlon+2)*(nlat+2)*4,MPI_INTEGER, &
       &                MPI_SUM,mpi_comm_world,mpi_code)
#endif

    allocate(neigh(n,4))
    k=0
    do j=1,nlat
    do i=1,nlon
      if (mod(k-iproc,nproc).eq.0) then
        neigh(1+k/nproc,:) = neighbours(i,j,:)
      endif
      k = k + 1
    enddo
    enddo

    deallocate(neighbours)

  endif

  end subroutine domain_decomposition

  ! ======================================================================

  ! Write a few members of the ensemble
  subroutine write_state_vector(filename,x,y,state)
  use netcdf
  implicit none

  character(len=256), intent(in) :: filename
  real(kind=8), dimension(:), intent(in) :: x, y, state

  real(kind=8), dimension(:,:), allocatable :: lon, lat, fullstate
  integer :: i, j, k

  integer :: is, ncid, idx, idy, idv, idlon, idlat

  ! allocate full grid and full state
  allocate(lon(nlon,nlat),lat(nlon,nlat),fullstate(nlon,nlat))
  lon = 0. ; lat = 0. ; fullstate = 0.

  ! merge processors contributions
  k=0
  do j=1,nlat
  do i=1,nlon
    if (mod(k-iproc,nproc).eq.0) then
      lon(i,j) = x(1+k/nproc)
      lat(i,j) = y(1+k/nproc)
      fullstate(i,j) = state(1+k/nproc)
    endif
    k = k + 1
  enddo
  enddo

#if defined MPI
  call MPI_ALLREDUCE (MPI_IN_PLACE,lon, nlon*nlat,MPI_DOUBLE_PRECISION, &
     &                MPI_SUM,mpi_comm_world,mpi_code)
  call MPI_ALLREDUCE (MPI_IN_PLACE,lat, nlon*nlat,MPI_DOUBLE_PRECISION, &
     &                MPI_SUM,mpi_comm_world,mpi_code)
  call MPI_ALLREDUCE (MPI_IN_PLACE,fullstate, nlon*nlat,MPI_DOUBLE_PRECISION, &
     &                MPI_SUM,mpi_comm_world,mpi_code)
#endif

  ! write full variable in output file
  if (iproc.eq.0) then
    is = NF90_CREATE(filename,NF90_CLOBBER,ncid)
    is = NF90_DEF_DIM(ncid,'lon',nlon,idx)
    is = NF90_DEF_DIM(ncid,'lat',nlat,idy)
    is = NF90_DEF_VAR(ncid,'lon',NF90_FLOAT,(/idx,idy/),idlon)
    is = NF90_DEF_VAR(ncid,'lat',NF90_FLOAT,(/idx,idy/),idlat)
    is = NF90_DEF_VAR(ncid,'state',NF90_FLOAT,(/idx,idy/),idv)
    is = NF90_ENDDEF(ncid)
    is = NF90_PUT_VAR(ncid,idlon,lon)
    is = NF90_PUT_VAR(ncid,idlat,lat)
    is = NF90_PUT_VAR(ncid,idv,fullstate)
    is = NF90_CLOSE(ncid)
  endif

  deallocate(lon,lat,fullstate)

  end subroutine write_state_vector

  ! ======================================================================

  ! Write a few members of the ensemble
  subroutine write_ensemble(filename,x,y,ensemble)
  use netcdf
  implicit none

  character(len=256), intent(in) :: filename
  real(kind=8), dimension(:), intent(in) :: x, y
  real(kind=8), dimension(:,:), intent(in) :: ensemble

  real(kind=8), dimension(:,:), allocatable :: lon, lat, fullstate
  integer :: i, j, k, imem, nmem

  integer :: is, ncid, idx, idy, idm, idv, idlon, idlat
  integer, dimension(3) :: vstart

  nmem = SIZE(ensemble,2)

  if (iproc.eq.0) then
    is = NF90_CREATE(filename,NF90_CLOBBER,ncid)
    is = NF90_DEF_DIM(ncid,'lon',nlon,idx)
    is = NF90_DEF_DIM(ncid,'lat',nlat,idy)
    is = NF90_DEF_DIM(ncid,'member',nmem,idm)
    is = NF90_DEF_VAR(ncid,'lon',NF90_FLOAT,(/idx,idy/),idlon)
    is = NF90_DEF_VAR(ncid,'lat',NF90_FLOAT,(/idx,idy/),idlat)
    is = NF90_DEF_VAR(ncid,'ensemble',NF90_FLOAT,(/idx,idy,idm/),idv)
    is = NF90_ENDDEF(ncid)
  endif

  ! allocate full grid
  allocate(lon(nlon,nlat),lat(nlon,nlat))
  lon = 0. ; lat = 0.

  ! merge location arrays
  k=0
  do j=1,nlat
  do i=1,nlon
    if (mod(k-iproc,nproc).eq.0) then
      lon(i,j) = x(1+k/nproc)
      lat(i,j) = y(1+k/nproc)
    endif
    k = k + 1
  enddo
  enddo

#if defined MPI
    call MPI_ALLREDUCE (MPI_IN_PLACE,lon, nlon*nlat,MPI_DOUBLE_PRECISION, &
       &                MPI_SUM,mpi_comm_world,mpi_code)
    call MPI_ALLREDUCE (MPI_IN_PLACE,lat, nlon*nlat,MPI_DOUBLE_PRECISION, &
       &                MPI_SUM,mpi_comm_world,mpi_code)
#endif

  if (iproc.eq.0) then
    is = NF90_PUT_VAR(ncid,idlon,lon)
    is = NF90_PUT_VAR(ncid,idlat,lat)
  endif

  ! allocate full state
  deallocate(lon,lat)
  allocate(fullstate(nlon,nlat))

  DO imem=1,nmem

    fullstate = 0.

    ! merge processors contributions
    k=0
    do j=1,nlat
    do i=1,nlon
      if (mod(k-iproc,nproc).eq.0) then
        fullstate(i,j) = ensemble(1+k/nproc,imem)
      endif
      k = k + 1
    enddo
    enddo

#if defined MPI
    call MPI_ALLREDUCE (MPI_IN_PLACE,fullstate, nlon*nlat,MPI_DOUBLE_PRECISION, &
       &                MPI_SUM,mpi_comm_world,mpi_code)
#endif

    ! write full variable in output file
    if (iproc.eq.0) then
      vstart=1 ; vstart(3)=imem
      is = NF90_PUT_VAR(ncid,idv,fullstate,start=vstart)
    endif

  enddo

  if (iproc.eq.0) then
    is = NF90_CLOSE(ncid)
  endif

  deallocate(fullstate)

  end subroutine write_ensemble

  ! ======================================================================

end program mcmc_ensemble_update
