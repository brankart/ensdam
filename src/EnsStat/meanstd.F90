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
!                        MODULE MEANSTD
!
!---------------------------------------------------------------------
! Compute mean and standard devisation of input ensemble
! by Jean-Michel Brankart, October 2018
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! ensemble_meanstd : compute mean and standard deviation from input ensemble
! update_meanstd : update mean and standard deviation with one additional input member
! ----------------------------------------------------------------------
MODULE ensdam_meanstd
#ifdef MPI_MODULE
      use mpi
#endif
      IMPLICIT NONE
      PRIVATE
#ifdef MPI_INCLUDE
      include "mpif.h"
#endif

      PUBLIC ensemble_meanstd, update_meanstd

      INTERFACE ensemble_meanstd
        MODULE PROCEDURE ensemble_meanstd_vector, ensemble_meanstd_variable
      END INTERFACE

      INTERFACE update_meanstd
        MODULE PROCEDURE update_meanstd_vector, update_meanstd_variable, &
                       & update_meanstd_vector_weight, update_meanstd_variable_weight
      END INTERFACE

      ! Public definitions needed by python/julia APIs
      PUBLIC ensemble_meanstd_vector, ensemble_meanstd_variable
      PUBLIC update_meanstd_vector, update_meanstd_variable, &
           & update_meanstd_vector_weight, update_meanstd_variable_weight

      ! Definition for MPI
#if defined MPI
      INTEGER, public, save  :: mpi_comm_meanstd=mpi_comm_world   ! definition of module global communicator
      INTEGER, private, save :: mpi_code
#endif

   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE ensemble_meanstd_vector ( ens, mean, std, weight )
!----------------------------------------------------------------------
!                  ***  ensemble_meanstd  ***
! 
! ** Purpose :   compute mean and standard deviation from input ensemble
! 
! ** Arguments :
!         ens    : input ensemble
!         mean   : ensemble mean
!         std    : ensemble standard deviation
!         weight : weight associated to eah ensemble member
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: ens
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: mean
        REAL(KIND=8), DIMENSION(:), INTENT( out ), OPTIONAL :: std
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ), OPTIONAL :: weight

        INTEGER :: jpi,jpm,jm,allocstat
        LOGICAL :: compute_std, weighted_stat
        REAL(KIND=8), DIMENSION(:), allocatable :: weightsum

        compute_std   = PRESENT(std)
        weighted_stat = PRESENT(weight)

        jpi = SIZE(ens,1)  ! Size of state vector
        jpm = SIZE(ens,2)  ! Size of ensemble

        IF (jpm.LT.2) STOP 'Ensemble size smaller than 2 in meanstd'

        IF (SIZE(mean,1).NE.jpi) STOP 'Inconsistent size in meanstd'
        mean = 0.

        IF (compute_std) THEN
          IF (SIZE(std,1).NE.jpi) STOP 'Inconsistent size in meanstd'
          std = 0.
        ENDIF

        IF (weighted_stat) THEN
          IF (SIZE(weight,1).NE.jpi) STOP 'Inconsistent size in meanstd'
          IF (SIZE(weight,2).NE.jpm) STOP 'Inconsistent size in meanstd'
          allocate( weightsum(jpi), stat=allocstat )
          IF (allocstat.NE.0) STOP 'Allocation error in meanstd'
          weightsum(:) = 0.
        ENDIF

        IF (.NOT.weighted_stat) THEN
          IF (compute_std) THEN
            DO jm=1,jpm
              CALL update_meanstd_vector( ens(:,jm), jm, mean, msqra=std )
            ENDDO
            std = SQRT( std / ( jpm -1 ) )
          ELSE
            DO jm=1,jpm
              CALL update_meanstd_vector( ens(:,jm), jm, mean )
            ENDDO
          ENDIF
        ELSE
          IF (compute_std) THEN
            DO jm=1,jpm
              CALL update_meanstd_vector_weight( ens(:,jm), weight(:,jm), weightsum, mean, msqra=std )
            ENDDO
            std = SQRT( std / weightsum )
          ELSE
            DO jm=1,jpm
              CALL update_meanstd_vector_weight( ens(:,jm), weight(:,jm), weightsum, mean )
            ENDDO
          ENDIF
        ENDIF
        
        IF (weighted_stat) deallocate(weightsum)

        END SUBROUTINE ensemble_meanstd_vector
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE ensemble_meanstd_variable ( ens, mean, std, weight )
!----------------------------------------------------------------------
!                  ***  ensemble_meanstd  ***
! 
! ** Purpose :   compute mean and standard deviation from input ensemble
! 
! ** Arguments :
!         ens    : input ensemble
!         mean   : ensemble mean
!         std    : ensemble standard deviation
!         weight : weight associated to eah ensemble member
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: ens
        REAL(KIND=8), INTENT( out ) :: mean
        REAL(KIND=8), INTENT( out ), OPTIONAL :: std
        REAL(KIND=8), DIMENSION(:), INTENT( in ), OPTIONAL :: weight

        INTEGER :: jpi,jpm,jm,allocstat
        LOGICAL :: compute_std, weighted_stat
        REAL(KIND=8) :: weightsum

        compute_std   = PRESENT(std)
        weighted_stat = PRESENT(weight)

        jpm = SIZE(ens,1)  ! Size of ensemble

        IF (jpm.LT.2) STOP 'Ensemble size smaller than 2 in meanstd'

        mean = 0.

        IF (compute_std) THEN
          std = 0.
        ENDIF

        IF (weighted_stat) THEN
          IF (SIZE(weight,1).NE.jpm) STOP 'Inconsistent size in meanstd'
          weightsum = 0.
        ENDIF

        IF (.NOT.weighted_stat) THEN
          IF (compute_std) THEN
            DO jm=1,jpm
              CALL update_meanstd_variable( ens(jm), jm, mean, msqra=std )
            ENDDO
            std = SQRT( std / ( jpm -1 ) )
          ELSE
            DO jm=1,jpm
              CALL update_meanstd_variable( ens(jm), jm, mean )
            ENDDO
          ENDIF
        ELSE
          IF (compute_std) THEN
            DO jm=1,jpm
              CALL update_meanstd_variable_weight( ens(jm), weight(jm), weightsum, mean, msqra=std )
            ENDDO
            std = SQRT( std / weightsum )
          ELSE
            DO jm=1,jpm
              CALL update_meanstd_variable_weight( ens(jm), weight(jm), weightsum, mean )
            ENDDO
          ENDIF
        ENDIF

        END SUBROUTINE ensemble_meanstd_variable
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE update_meanstd_vector ( vct, idx, mean, msqra )
!----------------------------------------------------------------------
!                  ***  update_meanstd_vector  ***
! 
! ** Purpose :   update mean and mean squared anomalies
! 
! ** Arguments :
!         vct   : additional vector
!         idx   : index of new member
!         mean  : ensemble mean
!         msqra : ensemble mean squared anomalies (with respect to the mean)
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: vct
        INTEGER, INTENT( in ) :: idx
        REAL(KIND=8), DIMENSION(:), INTENT( inout ) :: mean
        REAL(KIND=8), DIMENSION(:), INTENT( inout ), OPTIONAL :: msqra

        INTEGER :: jpi,ji
        LOGICAL :: compute_msqra

        compute_msqra  = PRESENT(msqra)

        jpi = SIZE(vct,1)  ! Size of state vector

        IF (SIZE(mean,1).NE.jpi) STOP 'Inconsistent size in meanstd'

        IF (compute_msqra) THEN
          IF (SIZE(msqra,1).NE.jpi) STOP 'Inconsistent size in meanstd'
        ENDIF

        IF (compute_msqra) THEN
          DO ji=1,jpi
            CALL update_meanstd_variable( vct(ji), idx, mean(ji), msqra=msqra(ji) )
          ENDDO
        ELSE
          DO ji=1,jpi
            CALL update_meanstd_variable( vct(ji), idx, mean(ji) )
          ENDDO
        ENDIF

        END SUBROUTINE update_meanstd_vector
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE update_meanstd_vector_weight ( vct, weight, weightsum, mean, msqra )
!----------------------------------------------------------------------
!                  ***  update_meanstd_vector  ***
! 
! ** Purpose :   update mean and mean squared anomalies (with unequal sample weight)
! 
! ** Arguments :
!         vct   : additional vector
!         weight    : weight for new value
!         weightsum : current sum of weights
!         mean  : ensemble mean
!         msqra : ensemble mean squared anomalies (with respect to the mean)
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: vct
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: weight
        REAL(KIND=8), DIMENSION(:), INTENT( inout ) :: weightsum
        REAL(KIND=8), DIMENSION(:), INTENT( inout ) :: mean
        REAL(KIND=8), DIMENSION(:), INTENT( inout ), OPTIONAL :: msqra

        INTEGER :: jpi,ji
        LOGICAL :: compute_msqra

        compute_msqra  = PRESENT(msqra)

        jpi = SIZE(vct,1)  ! Size of state vector

        IF (SIZE(mean,1).NE.jpi) STOP 'Inconsistent size in meanstd'
        IF (SIZE(weightsum,1).NE.jpi) STOP 'Inconsistent size in meanstd'

        IF (compute_msqra) THEN
          IF (SIZE(msqra,1).NE.jpi) STOP 'Inconsistent size in meanstd'
        ENDIF

        IF (compute_msqra) THEN
          DO ji=1,jpi
            CALL update_meanstd_variable_weight( vct(ji), weight(ji), weightsum(ji), mean(ji), msqra=msqra(ji) )
          ENDDO
        ELSE
          DO ji=1,jpi
            CALL update_meanstd_variable_weight( vct(ji), weight(ji), weightsum(ji), mean(ji) )
          ENDDO
        ENDIF

        END SUBROUTINE update_meanstd_vector_weight
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE update_meanstd_variable ( var, idx, mean, msqra )
!----------------------------------------------------------------------
!                  ***  update_meanstd_variable  ***
! 
! ** Purpose :   update mean and mean squared anomalies
!                (using Welford's incremental algorithm)
! 
! ** Arguments :
!         vct   : additional variable value
!         idx   : index of new member
!         mean  : ensemble mean
!         msqra : ensemble mean squared anomalies (with respect to the mean)
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), INTENT( in ) :: var
        INTEGER, INTENT( in ) :: idx
        REAL(KIND=8), INTENT( inout ) :: mean
        REAL(KIND=8), INTENT( inout ), OPTIONAL :: msqra

        REAL(KIND=8) :: misfit1, misfit2

        misfit1 = var - mean
        mean = mean + misfit1 / idx

        IF (PRESENT(msqra)) THEN
          misfit2 = var - mean
          msqra = msqra + misfit1 * misfit2
        ENDIF

        END SUBROUTINE update_meanstd_variable
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE update_meanstd_variable_weight ( var, weight, weightsum, mean, msqra )
!----------------------------------------------------------------------
!                  ***  update_meanstd_variable  ***
! 
! ** Purpose :   update mean and mean squared anomalies (handling unequal sample weights)
!                (using West's incremental algorithm)
! 
! ** Arguments :
!         var       : new variable value
!         weight    : weight for new value
!         weightsum : current sum of weights
!         mean      : ensemble mean
!         msqra     : ensemble mean squared anomalies (with respect to the mean)
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), INTENT( in ) :: var
        REAL(KIND=8), INTENT( in ) :: weight
        REAL(KIND=8), INTENT( inout ) :: weightsum
        REAL(KIND=8), INTENT( inout ) :: mean
        REAL(KIND=8), INTENT( inout ), OPTIONAL :: msqra

        REAL(KIND=8) :: misfit1, misfit2

        IF (weight.GT.0.) THEN

          weightsum = weightsum + weight

          misfit1 = var - mean
          mean = mean + misfit1 * weight / weightsum

          IF (PRESENT(msqra)) THEN
            misfit2 = var - mean
            msqra = msqra + misfit1 * misfit2 * weight
          ENDIF

        ENDIF

        END SUBROUTINE update_meanstd_variable_weight
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE ensdam_meanstd
