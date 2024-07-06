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
!                        MODULE ENSAUGM
!
!---------------------------------------------------------------------
! Ensemble augmentation by Schur product with large scale patterns
! by Jean-Michel Brankart, October 2018
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! sample_augmented_ensemble : iterate the Markov chains to sample the augmented ensemble
! newproduct : compute new multiple Schur product
! getproduct : compute specified multiple Schur product
! ----------------------------------------------------------------------
MODULE ensdam_ensaugm
#ifdef MPI_MODULE
      use mpi
#endif
#ifdef OPENACC
      use openacc
#endif
      use ensdam_schurprod
      use ensdam_storng
      IMPLICIT NONE
      PRIVATE
#ifdef MPI_INCLUDE
      include "mpif.h"
#endif

      PUBLIC sample_augmented_ensemble, newproduct, getproduct

      ! Public variables parameterizing the behaviour of the MCMC chain
      INTEGER, PUBLIC, SAVE :: ensaugm_chain_index=1  ! Current iteration index in MCMC chain (intialized to 1)
      LOGICAL, PUBLIC, SAVE :: ensaugm_with_renormalization=.FALSE.  ! Restore N(0,1) marginal distributions

      ! Definition for MPI
#if defined MPI
      INTEGER, PUBLIC, SAVE  :: mpi_comm_ensaugm=mpi_comm_world   ! definition of module global communicator
      INTEGER, save :: mpi_code
#endif

   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE sample_augmented_ensemble( maxchain, augens, ens, multiplicity )
!----------------------------------------------------------------------
! ** Purpose :   iterate to sample the augmented ensemble
! 
! ** Arguments :
!         maxchain      : number of iterations in the Markov chain
!         augens        : current version of augmented ensemble
!                         initialized to output of last call
!                         output as last iterate of the current call
!         ens           : input ensemble to be augmented,
!                         assumed available at several resolutions
!                         js=1 -> full resolution
!                         js>1 -> lower resolution versions of the same ensemble
!                         all marginal distributions must be N(0,1)
!         multiplicity  : multiplicity of each resolution to produce new members
!----------------------------------------------------------------------
        IMPLICIT NONE
        INTEGER, INTENT( in ) :: maxchain
        REAL(KIND=8), DIMENSION(:,:), INTENT( inout ) :: augens
        REAL(KIND=8), DIMENSION(:,:,:), INTENT( in ) :: ens
        INTEGER, DIMENSION(:), INTENT( in ) :: multiplicity

        INTEGER, DIMENSION(:), allocatable :: sample
        REAL(KIND=8), DIMENSION(:), allocatable :: schur_product
        REAL(KIND=8) :: coefficient, alpha, beta
        INTEGER :: jpi,jps,jpm,jpaug,jpfactor,js,jaug,jchain,allocstat

        ! Check size of input vectors
        jpi = SIZE(ens,1)  ! Size of observation vector
        jpm = SIZE(ens,2)  ! Size of ensemble
        jps = SIZE(ens,3)  ! Number of resolutions (or filtering length scales)
        jpaug = SIZE(augens,2)  ! Size of updated ensemble

        IF (SIZE(augens,1).NE.jpi) STOP 'Inconsistent size in ensaugm'
        IF (SIZE(multiplicity,1).NE.jps) STOP 'Inconsistent size in ensaugm'

        IF (multiplicity(1).NE.1) STOP 'Inconsistent multiplicity in ensaugm'
        DO js = 1,jps
          IF (multiplicity(js).LT.0) STOP 'Inconsistent multiplicity in ensaugm'
          IF (js.EQ.1) THEN
            IF (MOD(multiplicity(js),2).NE.1) STOP 'Inconsistent multiplicity in ensaugm'
          ELSE
            IF (MOD(multiplicity(js),2).NE.0) STOP 'Inconsistent multiplicity in ensaugm'
          ENDIF
        ENDDO

        ! Check number of factors use to sample a new member
        jpfactor = SUM(multiplicity(:))
        IF (jpfactor.GT.jpm) STOP 'More product factors than members in ensaugm'

        ! Allocate sample (member indices) used to build the new Schur product
        allocate( sample(jpfactor), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in ensaugm'

        ! Allocate vector to be used to store the Schur product
        allocate( schur_product(jpi), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in ensaugm'

        IF (ensaugm_chain_index.EQ.1) THEN
          ! Initialize the ensemble of Markov chains
          augens = 0.
        ENDIF

        ! Iterate the MCMC chains
        DO jchain = 1, maxchain

          ! Loop on Markov chains (size of augmented ensemble)
          DO jaug=1,jpaug

            ! get new Schur product
            CALL newproduct( schur_product, ens, multiplicity, sample )

            ! sample random coefficient, with N(0,1) distribution
            CALL kiss_gaussian(coefficient)

#if defined MPI
            CALL mpi_bcast(coefficient,1,mpi_double_precision,0,mpi_comm_ensaugm,mpi_code)
#endif

            ! update augmented member with new Schur product
            beta = 1. / SQRT(REAL(ensaugm_chain_index,8))
            alpha = SQRT( 1. - beta*beta )
            augens(:,jaug) = alpha * augens(:,jaug) + beta * coefficient * schur_product

          ENDDO

          ! Update chain index
          ensaugm_chain_index = ensaugm_chain_index + 1

        ENDDO

        deallocate(sample,schur_product)

        END SUBROUTINE sample_augmented_ensemble
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE newproduct( new, ens, multiplicity, sample )
!----------------------------------------------------------------------
! ** Purpose :   compute new multiple Schur product
! 
! ** Arguments :
!         new           : new multiple Schur product
!         ens           : input ensemble to be augmented, available at several resolution
!                         js=1 -> full resolution
!                         js>1 -> lower resolution versions of the same ensemble
!                         all marginal distributions must be N(0,1)
!         multiplicity  : multiplicity of each resolution
!         sample        : member indices used to build the new member
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: new
        REAL(KIND=8), DIMENSION(:,:,:), INTENT( in ) :: ens
        INTEGER, DIMENSION(:), INTENT( in ) :: multiplicity
        INTEGER, DIMENSION(:), INTENT( out ) :: sample

        INTEGER :: jpi,jps,jpm,jpfactor,jm,js,jmul,jfactor,allocstat,ji
        INTEGER, DIMENSION(:), allocatable :: tmp_sample

        jpi = SIZE(ens,1)  ! Size of state vector
        jpm = SIZE(ens,2)  ! Size of ensemble
        jps = SIZE(ens,3)  ! Number of resolutions (or filtering lngth scales)

        IF (SIZE(new,1).NE.jpi) STOP 'Inconsistent size in ensaugm'
        IF (SIZE(multiplicity,1).NE.jps) STOP 'Inconsistent size in ensaugm'

        IF (multiplicity(1).NE.1) STOP 'Inconsistent multiplicity in ensaugm'
        DO js = 1,jps
          IF (multiplicity(js).LT.0) STOP 'Inconsistent multiplicity in ensaugm'
          IF (js.EQ.1) THEN
            IF (MOD(multiplicity(js),2).NE.1) STOP 'Inconsistent multiplicity in ensaugm'
          ELSE
            IF (MOD(multiplicity(js),2).NE.0) STOP 'Inconsistent multiplicity in ensaugm'
          ENDIF
        ENDDO

        ! Compute total number of vectors to be multiplied
        jpfactor = SUM(multiplicity(:))
        IF (jpfactor.GT.jpm) STOP 'More product factors than members in ensaugm'
        IF (SIZE(sample,1).NE.jpfactor) STOP 'Inconsistent size in ensaugm'

        ! Sample ensemble member indices to generate new multiple Schur product
        allocate(tmp_sample(jpm),stat=allocstat)
        IF (allocstat.NE.0) STOP 'Allocation error in ensaugm'
        tmp_sample(:) = (/ (jm, jm=1,jpm) /)
        CALL kiss_sample(tmp_sample,jpm,jpfactor)
        sample = tmp_sample(1:jpfactor)
        deallocate(tmp_sample)
#if defined MPI
        CALL mpi_bcast(sample,jpfactor,mpi_integer,0,mpi_comm_ensaugm,mpi_code)
#endif
        
        ! Perform the Schur product with randomly selected members
        jfactor = 1
#if defined OPENACC
        !$acc data copyin(sample) present(ens, new)
        !$acc parallel loop
        DO ji=1,jpi
           new(ji) = ens(ji,sample(1),1)
        ENDDO
        !$acc end parallel loop
        !$acc end data
#else
        new(:) = ens(:,sample(1),1)
#endif
        DO js = 2,jps
          DO jmul = 1,multiplicity(js)
            jfactor = jfactor+1
            IF (ensaugm_with_renormalization) THEN
              CALL schurprod( new(:), ens(:,sample(jfactor),js) )
            ELSE
#if defined OPENACC
              !$acc data copyin(sample) present(ens, new)
              !$acc parallel loop
              DO ji=1,jpi
                new(ji) = new(ji) * ens(ji,sample(jfactor),js)
              ENDDO
              !$acc end parallel loop
              !$acc end data
#else
              new(:) = new(:) * ens(:,sample(jfactor),js)
#endif
            ENDIF
          ENDDO
        ENDDO

        END SUBROUTINE newproduct
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE getproduct( new, ens, multiplicity, sample )
!----------------------------------------------------------------------
! ** Purpose :   compute required multiple Schur product
! 
! ** Arguments :
!         new           : required multiple Schur product
!         ens           : input ensemble to be augmented, available at several resolution
!                         js=1 -> full resolution
!                         js>1 -> lower resolution versions of the same ensemble
!                         all marginal distributions must be N(0,1)
!         multiplicity  : multiplicity of each resolution
!         sample        : member indices used to build the new member
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: new
        REAL(KIND=8), DIMENSION(:,:,:), INTENT( in ) :: ens
        INTEGER, DIMENSION(:), INTENT( in ) :: multiplicity
        INTEGER, DIMENSION(:), INTENT( in ) :: sample

        INTEGER :: jpi,jps,jpm,jpfactor,jm,js,jmul,jfactor,allocstat,ji

        jpi = SIZE(ens,1)  ! Size of state vector
        jpm = SIZE(ens,2)  ! Size of ensemble
        jps = SIZE(ens,3)  ! Number of resolutions (or filtering lngth scales)

        IF (SIZE(new,1).NE.jpi) STOP 'Inconsistent size in ensaugm'
        IF (SIZE(multiplicity,1).NE.jps) STOP 'Inconsistent size in ensaugm'

        IF (multiplicity(1).NE.1) STOP 'Inconsistent multiplicity in ensaugm'
        DO js = 1,jps
          IF (multiplicity(js).LT.0) STOP 'Inconsistent multiplicity in ensaugm'
          IF (js.EQ.1) THEN
            IF (MOD(multiplicity(js),2).NE.1) STOP 'Inconsistent multiplicity in ensaugm'
          ELSE
            IF (MOD(multiplicity(js),2).NE.0) STOP 'Inconsistent multiplicity in ensaugm'
          ENDIF
        ENDDO

        ! Compute total number of vectors to be multiplied
        jpfactor = SUM(multiplicity(:))
        IF (jpfactor.GT.jpm) STOP 'More product factors than members in ensaugm'
        IF (SIZE(sample,1).NE.jpfactor) STOP 'Inconsistent size in ensaugm'

        ! Perform the Schur product with required selected members
        jfactor = 1
#if defined OPENACC
        !$acc data copyin(sample) present(ens, new)
        !$acc parallel loop
        DO ji=1,jpi
           new(ji) = ens(ji,sample(1),1)
        ENDDO
        !$acc end parallel loop
        !$acc end data
#else
        new(:) = ens(:,sample(1),1)
#endif
        DO js = 2,jps
          DO jmul = 1,multiplicity(js)
            jfactor = jfactor+1
            IF (ensaugm_with_renormalization) THEN
              CALL schurprod( new(:), ens(:,sample(jfactor),js) )
            ELSE
#if defined OPENACC
              !$acc data copyin(sample) present(ens, new)
              !$acc parallel loop
              DO ji=1,jpi
                new(ji) = new(ji) * ens(ji,sample(jfactor),js)
              ENDDO
              !$acc end parallel loop
              !$acc end data
#else
              new(:) = new(:) * ens(:,sample(jfactor),js)
#endif
            ENDIF
          ENDDO
        ENDDO

        END SUBROUTINE getproduct
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE ensdam_ensaugm
