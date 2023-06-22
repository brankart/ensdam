!---------------------------------------------------------------------
! Copyright: CNRS - Université de Grenoble Alpes
!
! Contributors : Jean-Michel Brankart
!
! Jean-Michel.Brankart@univ-grenoble-alpes.fr
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software.  You can  use,
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".
!
!---------------------------------------------------------------------
!
!                        MODULE MCMC_UPDATE
!
!---------------------------------------------------------------------
! Ensemble observational update using Monte Carlo Markov Chains
! by Jean-Michel Brankart, October 2018
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! mcmc_iteration : iterate to sample the posterior probability distribution
! ----------------------------------------------------------------------
MODULE ensdam_mcmc_update
      use ensdam_ensaugm
      use ensdam_storng
      IMPLICIT NONE
      PRIVATE

      PUBLIC mcmc_iteration

      ! Public variables parameterizing the behaviour of the MCMC chain
      INTEGER, PUBLIC, SAVE :: mcmc_index=1  ! Current iteration index in MCMC chain (intialized to 1)
      LOGICAL, PUBLIC, SAVE :: mcmc_zero_start=.TRUE. ! Start MCMC chains from zero
      INTEGER, PUBLIC, SAVE :: mcmc_control_print=1000  ! Number of iterations between control prints
      INTEGER, PUBLIC, SAVE :: mcmc_convergence_check=1000   ! Number of iterations between convergence checks
      LOGICAL, PUBLIC, SAVE :: mcmc_convergence_stop=.FALSE. ! Stop iterating at convergence
      INTEGER, PUBLIC, SAVE :: mcmc_member_test=1  ! Number of test to perform with the same Schur product
      LOGICAL, PUBLIC, SAVE :: mcmc_proposal=.FALSE. ! Input is a proposal distribution, not a prior ensemble
      REAL(KIND=8), PUBLIC, SAVE :: mcmc_proposal_std=1._8 ! Standard deviation of proposal distribution
      REAL(KIND=8), PUBLIC, SAVE :: mcmc_schedule=0._8 ! MCMC schedule, can be changed in cost_jo, as a function of Jo

#if defined MPI
      ! Public definitions for MPI
      include "mpif.h"
      INTEGER, PUBLIC, SAVE  :: mpi_comm_mcmc_update=mpi_comm_world   ! definition of module global communicator
      INTEGER, save :: mpi_code
#endif
      INTEGER, save :: iproc=0

      ! Module private definitions for observation cost function
      INTERFACE
        FUNCTION callback_jo(v)            ! callback function for jo = -log[p(yo|Hx)]
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), intent(in) :: v
        REAL(KIND=8) :: callback_jo
        END FUNCTION
      END INTERFACE

      PROCEDURE(callback_jo), POINTER, SAVE :: cost_jo
      REAL(KIND=8), DIMENSION(:), SAVE, allocatable :: cost_jo_saved  ! saved value of the observation cost function

      ! Module private definitions for convergence test
      INTERFACE
        FUNCTION callback_test(upens,upxens)   ! callback function for convergence test
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), intent(in) :: upens
        REAL(KIND=8), DIMENSION(:,:), intent(in), optional :: upxens
        LOGICAL :: callback_test
        END FUNCTION
      END INTERFACE

      PROCEDURE(callback_test), POINTER, SAVE :: convergence_test
  
   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE mcmc_iteration( maxchain, upens, ens, multiplicity, my_jo, my_test, upxens, xens)
!----------------------------------------------------------------------
! ** Purpose :   iterate to sample the posterior probability distribution
! 
! ** Arguments :
!         maxchain      : maximum number of successful iterations in the Markov chain
!         upens         : current version of updated ensemble (all variables needed to apply observation operator),
!                         initialized to output of last call
!                         output as last iterate of the current call
!         ens           : input ensemble to be updated (all variables needed to apply observation operator)),
!                         assumed available at several resolutions
!                         js=1 -> full resolution
!                         js>1 -> lower resolution versions of the same ensemble
!                         all marginal distributions must be N(0,1)
!         multiplicity  : multiplicity of each resolution to produce new members
!         my_jo         : callback routine to observation cost function: jo = -log[p(yo|Hx)]
!         my_test       : callback routine for convergence test
!         upxens        : current version of updated ensemble (extra variables, not needed to apply observation operator),
!                         initialized to output of last call
!                         output as last iterate of the current call
!         xens          : input ensemble to be updated (extra variables, not needed to apply observation operator)),
!                         assumed available at several resolutions
!                         js=1 -> full resolution
!                         js>1 -> lower resolution versions of the same ensemble
!                         all marginal distributions must be N(0,1)
!----------------------------------------------------------------------
        IMPLICIT NONE
        INTEGER, INTENT( in ) :: maxchain
        REAL(KIND=8), DIMENSION(:,:), INTENT( inout ) :: upens
        REAL(KIND=8), DIMENSION(:,:,:), INTENT( in ) :: ens
        INTEGER, DIMENSION(:), INTENT( in ) :: multiplicity
        PROCEDURE(callback_jo) :: my_jo
        PROCEDURE(callback_test), OPTIONAL :: my_test
        REAL(KIND=8), DIMENSION(:,:), INTENT( inout ), OPTIONAL :: upxens
        REAL(KIND=8), DIMENSION(:,:,:), INTENT( in ), OPTIONAL :: xens

        INTEGER :: jpi,jps,jpm,jpup,jpextra,jpfactor,js,jchain,allocstat
        LOGICAL :: extra_variables, convergence_test

#if defined MPI
        call mpi_comm_rank(mpi_comm_world,iproc,mpi_code)
#endif

        ! Definition of observation cost function: jo = -log[p(yo|Hx)]
        cost_jo => my_jo

        ! Check size of input vectors
        jpi = SIZE(ens,1)  ! Size of observation vector
        jpm = SIZE(ens,2)  ! Size of ensemble
        jps = SIZE(ens,3)  ! Number of resolutions (or filtering length scales)
        jpup = SIZE(upens,2)  ! Size of updated ensemble

        IF (SIZE(upens,1).NE.jpi) STOP 'Inconsistent size in mcmc_update'
        IF (SIZE(multiplicity,1).NE.jps) STOP 'Inconsistent size in mcmc_update'

        IF (multiplicity(1).NE.1) STOP 'Inconsistent multiplicity in mcmc_update'
        DO js = 1,jps
          IF (multiplicity(js).LT.0) STOP 'Inconsistent multiplicity in mcmc_update'
          IF (js.EQ.1) THEN
            IF (MOD(multiplicity(js),2).NE.1) STOP 'Inconsistent multiplicity in mcmc_update'
          ELSE
            IF (MOD(multiplicity(js),2).NE.0) STOP 'Inconsistent multiplicity in mcmc_update'
          ENDIF
        ENDDO

        ! Check number of factors use to sample a new member
        jpfactor = SUM(multiplicity(:))
        IF (jpfactor.GT.jpm) STOP 'More product factors than members in mcmc_update'

        ! Are there extra variables to update ?
        extra_variables = PRESENT(upxens)
        IF (extra_variables) THEN
          IF (.NOT.PRESENT(xens)) STOP 'Inconsistent optional arguments in mcmc_update'
          jpextra = SIZE(xens,1)  ! Number of extra variables
          IF (SIZE(xens,2).NE.jpm) STOP 'Inconsistent size in mcmc_update'
          IF (SIZE(xens,3).NE.jps) STOP 'Inconsistent size in mcmc_update'
          IF (SIZE(upxens,2).NE.jpup) STOP 'Inconsistent size in mcmc_update'
        ELSE
          IF (PRESENT(xens)) STOP 'Inconsistent optional arguments in mcmc_update'
        ENDIF

        ! Perform convergence test ?
        convergence_test = PRESENT(my_test)

        ! Initialize Markov chain
        CALL mcmc_init( upens )

        ! Iterate the MCMC chain
        DO jchain = 1, maxchain

          ! Get next accepted draw of the Markov Chain
          IF (extra_variables) THEN
            CALL next_accepted_draw( upens, ens, multiplicity, upxens, xens )
          ELSE
            CALL next_accepted_draw( upens, ens, multiplicity )
          ENDIF

          IF (convergence_test) THEN
            IF (MOD(jchain,mcmc_convergence_check).EQ.0) THEN
              ! Check convergence of the MCMC chain
              IF (extra_variables) THEN
                IF (my_test(upens,upxens).AND.mcmc_convergence_stop) EXIT
              ELSE
                IF (my_test(upens).AND.mcmc_convergence_stop) EXIT
              ENDIF
            ENDIF
          ENDIF

        ENDDO

        END SUBROUTINE mcmc_iteration
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE mcmc_init( upens )
!----------------------------------------------------------------------
! ** Purpose :   initialize Markov chain
! 
! ** Arguments :
!         upens : current version of the updated ensemble (all variables needed to apply observation operator),
!                 initialized to last accepted draw, output as next accepted draw
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT( inout ) :: upens

        INTEGER :: jup, jpup, allocstat

        jpup = SIZE(upens,2)  ! Size of updated ensemble

        IF (mcmc_index.EQ.1) THEN

          ! Initialize the ensemble of Markov chains
          IF (mcmc_zero_start) upens = 0.

          ! Allocate cost function (one for each chain)
          IF (allocated(cost_jo_saved)) deallocate(cost_jo_saved)
          allocate( cost_jo_saved(jpup), stat=allocstat )
          IF (allocstat.NE.0) STOP 'Allocation error in mcmc_update'

          ! Initialize cost function
          IF (mcmc_zero_start) THEN
            cost_jo_saved(:) = cost_jo( upens(:,1) )
            IF (iproc.eq.0) PRINT *, 'Initial Jo:',cost_jo_saved(1)
          ELSE
            DO jup=1,jpup
              cost_jo_saved(jup) = cost_jo( upens(:,jup) )
              IF (iproc.eq.0) PRINT *, 'Initial Jo:',cost_jo_saved(jup)
            ENDDO
          ENDIF

        ELSE

          ! Allocate cost function (one for each chain)
          IF (.NOT.allocated(cost_jo_saved)) THEN
            allocate( cost_jo_saved(jpup), stat=allocstat )
            IF (allocstat.NE.0) STOP 'Allocation error in mcmc_update'
          ENDIF

          ! Initialize cost function
          DO jup=1,jpup
            cost_jo_saved(jup) = cost_jo( upens(:,jup) )
            IF (iproc.eq.0) PRINT *, 'Restart Jo:',cost_jo_saved(jup)
          ENDDO

        ENDIF

        END SUBROUTINE mcmc_init
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE next_accepted_draw( upens, ens, multiplicity, upxens, xens )
!----------------------------------------------------------------------
! ** Purpose :   get next accepted draw of the Markov Chain
! 
! ** Arguments :
!         upens         : current version of the updated ensemble (all variables needed to apply observation operator),
!                         initialized to last accepted draw
!                         output as next accepted draw
!         ens           : input ensemble to be updated (all variables needed to apply observation operator),
!                         assumed available at several resolutions
!                         js=1 -> full resolution
!                         js>1 -> lower resolution versions of the same ensemble
!                         all marginal distributions must be N(0,1)
!         multiplicity  : multiplicity of each resolution to produce new members
!         upxens        : current version of updated ensemble (extra variables, not needed to apply observation operator),
!                         initialized to output of last call
!                         output as last iterate of the current call
!         xens          : input ensemble to be updated (extra variables, not needed to apply observation operator)),
!                         assumed available at several resolutions
!                         js=1 -> full resolution
!                         js>1 -> lower resolution versions of the same ensemble
!                         all marginal distributions must be N(0,1)
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT( inout ) :: upens
        REAL(KIND=8), DIMENSION(:,:,:), INTENT( in ) :: ens
        INTEGER, DIMENSION(:), INTENT( in ) :: multiplicity
        REAL(KIND=8), DIMENSION(:,:), INTENT( inout ), OPTIONAL :: upxens
        REAL(KIND=8), DIMENSION(:,:,:), INTENT( in ), OPTIONAL :: xens

        INTEGER, DIMENSION(:), allocatable :: sample
        REAL(KIND=8), DIMENSION(:), allocatable :: vtest, vextra
        REAL(KIND=8) :: coefficient, alpha, beta
        INTEGER :: jpi,jps,jpm,jpup,jpextra,jpfactor,js,jup,jtest,allocstat
        LOGICAL :: extra_variables

        ! Are there extra variables to update ?
        extra_variables = PRESENT(upxens)

        ! Size of arrays
        jpi = SIZE(ens,1)  ! Size of observation vector
        jpm = SIZE(ens,2)  ! Size of ensemble
        jps = SIZE(ens,3)  ! Number of resolutions (or filtering length scales)
        jpup = SIZE(upens,2)  ! Size of updated ensemble
        jpfactor = SUM(multiplicity(:))
        IF (extra_variables) THEN
          jpextra = SIZE(xens,1)  ! Number of extra variables
        ENDIF

        ! Allocate sample (member indices) used to build the new Schur product
        allocate( sample(jpfactor), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in mcmc_update'

        ! Allocate test vector to be used as tentative new draw
        allocate( vtest(jpi), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in mcmc_update'

        ! Allocate vector to store new ensemble member for extra variables
        IF (extra_variables) THEN
          allocate( vextra(jpextra), stat=allocstat )
          IF (allocstat.NE.0) STOP 'Allocation error in mcmc_update'
        ENDIF

        ! Loop on Markov chains (size of updated ensemble)
        next_chain : DO jup=1,jpup

          ! Iterate until new draw is accepted
          next_draw : DO

            ! select proposal distribution, and prepare perturbation deterministic coefficients
            IF (mcmc_proposal) THEN
              ! assuming that the input is the proposal distribution
              alpha = 1._8
              beta = mcmc_proposal_std
            ELSE
              ! assuming that the input is a prior ensemble to update
              ! 1. Default beta ** 2
              beta = 1._8 / REAL(mcmc_index,8)
              ! 2. Upgrade with user-defined schedule (default: mcmc_schedule=0.)
              beta = MAX ( beta , MIN (  0.5_8 , mcmc_schedule ) )
              ! 3. Compute alpha ** 2
              alpha = 1._8 - beta
              ! 4. Compute alpha and beta
              alpha = SQRT( alpha ) ; beta = SQRT( beta )
            ENDIF

            ! get new Schur product (using the ensemble augmentation tool)
            CALL newproduct( vtest, ens, multiplicity, sample )

            ! use it to test perturbation(s) to the current vector
            next_test : DO jtest=1,mcmc_member_test

              CALL kiss_gaussian(coefficient)   ! coefficient ~ N(0,1)
#if defined MPI
              CALL mpi_bcast(coefficient,1,mpi_double_precision,0,mpi_comm_mcmc_update,mpi_code)
#endif
              ! get draw from proposal distribution
              vtest = alpha * upens(:,jup) + beta * coefficient * vtest

              IF ( accept_new_draw( vtest, jup ) ) THEN
                ! update ensemble member
                upens(:,jup) = vtest
                ! update extra variables of the same member
                IF (extra_variables) THEN
                  ! get Schur product of extra variables corresponding to the one obtained by newproduct
                  ! (i.e. using the same sample to compute the product)
                  CALL getproduct( vextra, xens, multiplicity, sample )
                  ! perform the same linear combination, with the same random coefficient
                  upxens(:,jup) = alpha * upxens(:,jup) + beta * coefficient * vextra
                ENDIF
                ! move to next chain
                EXIT next_draw
              ENDIF

            ENDDO next_test

          ENDDO next_draw

        ENDDO next_chain

        ! Update chain index
        mcmc_index = mcmc_index + 1

        deallocate(vtest,sample)

        END SUBROUTINE next_accepted_draw
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        FUNCTION accept_new_draw( vtest, jup )
!----------------------------------------------------------------------
! ** Purpose :   acceptance test for new draw
! 
! ** Arguments :
!         vtest   : new draw to test
!         jup     : index of MCMC chain
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: vtest
        INTEGER, INTENT( in ) :: jup
        LOGICAL :: accept_new_draw

        REAL(KIND=8) :: cost_jo_test, uran, aprob
        INTEGER :: jpi, allocstat

        jpi = SIZE(vtest,1)  ! Size of observation vector

        ! Evaluate observation cost function
        cost_jo_test = cost_jo( vtest )

        ! Compute acceptance probability
        aprob = EXP ( cost_jo_saved(jup) - cost_jo_test )

        ! Accept or reject test draw
        CALL kiss_uniform( uran )  ! uran ~ U(0,1)
#if defined MPI
        CALL mpi_bcast(uran,1,mpi_double_precision,0,mpi_comm_mcmc_update,mpi_code)
#endif
        accept_new_draw = uran .LE. aprob

        ! Update cost function if accepted
        IF (accept_new_draw) THEN
          cost_jo_saved(jup) = cost_jo_test
          IF ( (iproc.EQ.0) .AND. (MOD(mcmc_index,mcmc_control_print).EQ.0) ) PRINT *, 'Jo(',mcmc_index,') =',cost_jo_test
        ENDIF

        END FUNCTION accept_new_draw
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE ensdam_mcmc_update
