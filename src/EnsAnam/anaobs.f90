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
!                        MODULE ANAOBS
!
!---------------------------------------------------------------------
! Perform anamorphosis transformation:
! -Transformation of observation probability distribution
! by Jean-Michel Brankart, September 2018
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! ana_obs : transformation of observation probability distribution
! ana_obs_sym : transformation of observation probability distribution
!               (simplified for symmetric observation error distribution)
! ----------------------------------------------------------------------
MODULE ensdam_anaobs
      use ensdam_anaqua
      use ensdam_anatra
      use ensdam_obserror
      use ensdam_storng
      IMPLICIT NONE
      PRIVATE

      PUBLIC ana_obs, ana_obs_sym

   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE ana_obs( anaobs, obsens, obs, obserror, quadef, quaref )
!----------------------------------------------------------------------
! ** Purpose : transformation of observation probability distribution
!
! ** Arguments :
!         anaobs   : sample of transformed observation probability distribution
!         obsens   : ensemble equivalent to observations
!         obs      : observation vector to be transformed
!         obserror : spread of observation error
!                    (precise meaning depending on the shape of the distribution)
!         quadef   : definition of the quantiles used for anamorphosis
!         quaref   : quantiles of the reference/target distribution
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT( out ) :: anaobs
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: obsens
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: obs
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: obserror
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: quadef
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: quaref

        INTEGER :: jpi,jpm,jpq,jps,ji,jm,jq,js,allocstat
        REAL(KIND=8), DIMENSION(:), allocatable :: ens_pert, ens_qua
        REAL(KIND=8) :: rank

        jpi = SIZE(obsens,1)   ! size of observation vector
        jpm = SIZE(obsens,2)   ! ensemble size
        jpq = SIZE(quadef,1)   ! number of quantiles used for anamorphosis
        jps = SIZE(anaobs,2)   ! output sample size
        IF (SIZE(anaobs,1).NE.jpi) STOP 'Inconsistent size in ana_obs'
        IF (SIZE(obs,1).NE.jpi) STOP 'Inconsistent size in ana_obs'
        IF (SIZE(obserror,1).NE.jpi) STOP 'Inconsistent size in ana_obs'
        IF (SIZE(quaref,1).NE.jpq) STOP 'Inconsistent size in ana_obs'

        IF (MINVAL(quadef).LT.0.0) STOP 'Bad quantile definition in ana_obs'
        IF (MAXVAL(quadef).GT.1.0) STOP 'Bad quantile definition in ana_obs'

        allocate ( ens_pert(1:jpm), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in ana_obs'
        allocate ( ens_qua(1:jpq), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in ana_obs'

        DO js = 1,jps
        DO ji = 1,jpi

          ! Use random rank to transform observation in probability concentration
          call kiss_uniform(rank)

          ! Perturb input ensemble with observation error
          call obserror_perturbation( ens_pert(:), obsens(ji,:), obserror(ji) )
          ! Compute quantiles of perturbed ensemble
          call ens_quantiles( ens_qua, ens_pert, quadef )
          ! Forward anamorphosis of observation with perturbed transformation
          anaobs(ji,js) = obs(ji)
          call ana_forward( anaobs(ji,js), ens_qua, quaref, rank=rank )

        ENDDO
        ENDDO

        deallocate(ens_pert,ens_qua)

        END SUBROUTINE ana_obs
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE ana_obs_sym( anaobs, obs, obserror, obsqua, quaref )
!----------------------------------------------------------------------
! ** Purpose : transformation of observation probability distribution
!              (simplified for symmetric observation error distribution)
! 
! ** Arguments :
!         anaobs   : sample of transformed observation probability distribution
!         obs      : observation vector to be transformed
!         obserror : spread of observation error
!                    (precise meaning depending on the shape of the distribution)
!         obsqua   : quantiles of ensemble equivalents to observations
!         quaref   : quantiles of the reference/target distribution
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT( out ) :: anaobs
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: obs
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: obserror
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: obsqua
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: quaref

        INTEGER :: jpi,jpq,jps,ji,jq,js
        REAL(KIND=8) :: obs_pert, rank

        jpi = SIZE(obs,1)   ! size of observation vector
        jpq = SIZE(obsqua,2)   ! number of quantiles used for anamorphosis
        jps = SIZE(anaobs,2)   ! output sample size
        IF (SIZE(anaobs,1).NE.jpi) STOP 'Inconsistent size in ana_obs_sym'
        IF (SIZE(obserror,1).NE.jpi) STOP 'Inconsistent size in ana_obs_sym'
        IF (SIZE(obsqua,1).NE.jpi) STOP 'Inconsistent size in ana_obs_sym'
        IF (SIZE(quaref,1).NE.jpq) STOP 'Inconsistent size in ana_obs_sym'

        DO js = 1,jps
        DO ji = 1,jpi

          ! Use random rank to transform observation in probability concentration
          call kiss_uniform(rank)

          ! Perturb observation with observation error
          call obserror_perturbation_sym( obs_pert, obs(ji), obserror(ji) )
          ! Forward anamorphosis of perturbed observation
          anaobs(ji,js) = obs_pert
          call ana_forward( anaobs(ji,js), obsqua(ji,:), quaref, rank=rank )

        ENDDO
        ENDDO

        END SUBROUTINE ana_obs_sym
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE obserror_perturbation( ens_pert, ens, obserror )
!----------------------------------------------------------------------
! ** Purpose : perturb input ensemble with observation error
! 
! ** Arguments :
!         ens_pert : perturbed ensemble
!         ens      : input ensemble
!         obserror : spread of observation error
!                    (precise meaning depending on the shape of the distribution)
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: ens_pert
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: ens
        REAL(KIND=8), INTENT( in ) :: obserror

        INTEGER :: jpm,jm

        jpm = SIZE(ens,1)   ! ensemble size
        IF (SIZE(ens_pert,1).NE.jpm) STOP 'Inconsistent size in obserror_perturbation'

        ! Check observation error type
        SELECT CASE(obserror_type)
        CASE('normal','gaussian')
           print *, 'Warning: this algorithm is non-optimal for symmetric observation error distribution'
        CASE('lognormal','gamma','beta')
        CASE DEFAULT
           STOP 'Bad observation error cdf in anaobs'
        END SELECT

        ! Apply observation error with the same rank
        ! to all members of the ensemble
        CALL obserror_sample( ens(:), obserror, ens_pert(:), uniform_rank=.TRUE. )

        END SUBROUTINE obserror_perturbation
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE obserror_perturbation_sym( obs_pert, obs, obserror )
!----------------------------------------------------------------------
! ** Purpose : perturb input observation with observation error
! 
! ** Arguments :
!         obs_pert : perturbed observation
!         obs      : input observation
!         obserror : spread of observation error
!                    (precise meaning depending on the shape of the distribution)
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), INTENT( out ) :: obs_pert
        REAL(KIND=8), INTENT( in ) :: obs
        REAL(KIND=8), INTENT( in ) :: obserror

        ! Check observation error type
        SELECT CASE(obserror_type)
        CASE('normal','gaussian')
        CASE('lognormal','gamma','beta')
           STOP 'This algorithm is inappropriate for non-symmetric observation error distribution'
        CASE DEFAULT
           STOP 'Bad observation error cdf in anaobs'
        END SELECT

        ! Sample perturbed observation
        obs_pert = obserror_sample( obs, obserror )

        END SUBROUTINE obserror_perturbation_sym
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE ensdam_anaobs
