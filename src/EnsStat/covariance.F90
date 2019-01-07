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
!                        MODULE COVARIANCE
!
!---------------------------------------------------------------------
! Compute covariance from input ensemble
! by Jean-Michel Brankart, October 2018
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! ensemble_correlation : compute covariance from input ensemble
! ensemble_representer : compute representer from input ensemble
! ensemble_covariance : compute covariance from input ensemble
! update_meancov : update mean and covariance with one additional input member
! ----------------------------------------------------------------------
MODULE ensdam_covariance
#if defined MPI
      use mpi
#endif
      use ensdam_meanstd
      IMPLICIT NONE
      PRIVATE

      PUBLIC ensemble_correlation, ensemble_representer, ensemble_covariance, update_meancov

      INTERFACE update_meancov
        MODULE PROCEDURE update_meancov_vector, update_meancov_variable, &
                       & update_meancov_vector_weight, update_meancov_variable_weight
      END INTERFACE

      ! Public definitions needed by python/julia APIs
      PUBLIC update_meancov_vector, update_meancov_variable, &
           & update_meancov_vector_weight, update_meancov_variable_weight

      REAL(KIND=8), public, save :: correlation_missing_value = -9999.

      ! Definition for MPI
#if defined MPI
      INTEGER, public, save :: mpi_comm_covariance=mpi_comm_world   ! definition of module global communicator
      INTEGER, private, save :: mpi_code
#endif

   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE ensemble_correlation ( ens, ensref, correl, weight )
!----------------------------------------------------------------------
!                  ***  ensemble_correlation  ***
! 
! ** Purpose :   compute correlation from input ensemble
! 
! ** Arguments :
!         ens    : input ensemble
!         ensref : scalar reference ensemble (with respect to which computing covariances)
!         correl : ensemble covariance with respect to reference variable
!         weight : weight associated to eah ensemble member
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: ens
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: ensref
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: correl
        REAL(KIND=8), DIMENSION(:), INTENT( in ), OPTIONAL :: weight

        REAL(KIND=8), DIMENSION(:), allocatable :: mean, std
        REAL(KIND=8) :: meanref, stdref

        INTEGER :: jpi,jpm,jm,ji,allocstat
        LOGICAL :: weighted_stat

        weighted_stat = PRESENT(weight)

        jpi = SIZE(ens,1)  ! Size of state vector
        jpm = SIZE(ens,2)  ! Size of ensemble

        IF (jpm.LT.2) STOP 'Ensemble size smaller than 2 in covariance'
        IF (SIZE(correl,1).NE.jpi) STOP 'Inconsistent size in covariance'
        IF (SIZE(ensref,1).NE.jpm) STOP 'Inconsistent size in covariance'

        IF (weighted_stat) THEN
          IF (SIZE(weight,1).NE.jpm) STOP 'Inconsistent size in covariance'
        ENDIF

        allocate( mean(jpi), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in covariance'
        allocate( std(jpi), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in covariance'

        IF (.NOT.weighted_stat) THEN
          CALL ensemble_meanstd ( ensref, meanref, stdref )
          IF (stdref.GT.0.) THEN
            CALL ensemble_covariance ( ens, ensref, correl )
            CALL ensemble_meanstd ( ens, mean, std )
            WHERE (std.GT.0.)
              correl = correl / std / stdref
            ELSEWHERE
              correl = correlation_missing_value
            ENDWHERE
          ELSE
            correl = correlation_missing_value
          ENDIF
        ELSE
          CALL ensemble_meanstd ( ensref, meanref, stdref, weight )
          IF (stdref.GT.0.) THEN
            CALL ensemble_covariance ( ens, ensref, correl, weight )
            DO ji=1,jpi
              CALL ensemble_meanstd ( ens(ji,:), mean(ji), std(ji), weight(:) )
            ENDDO
            WHERE (std.GT.0.)
              correl = correl / std / stdref
            ELSEWHERE
              correl = correlation_missing_value
            ENDWHERE
          ELSE
            correl = correlation_missing_value
          ENDIF
        ENDIF

        deallocate(mean,std)

        END SUBROUTINE ensemble_correlation
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE ensemble_representer ( ens, ensref, representer, weight )
!----------------------------------------------------------------------
!                  ***  ensemble_representer  ***
! 
! ** Purpose :   compute representer from input ensemble
! 
! ** Arguments :
!         ens    : input ensemble
!         ensref : scalar reference ensemble (with respect to which computiing covariances)
!         representer : ensemble covariance with respect to reference variable
!         weight : weight associated to eah ensemble member
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: ens
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: ensref
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: representer
        REAL(KIND=8), DIMENSION(:), INTENT( in ), OPTIONAL :: weight

        REAL(KIND=8) :: meanref, stdref

        INTEGER :: jpi,jpm,jm,allocstat
        LOGICAL :: weighted_stat

        weighted_stat = PRESENT(weight)

        jpi = SIZE(ens,1)  ! Size of state vector
        jpm = SIZE(ens,2)  ! Size of ensemble

        IF (jpm.LT.2) STOP 'Ensemble size smaller than 2 in covariance'
        IF (SIZE(representer,1).NE.jpi) STOP 'Inconsistent size in covariance'
        IF (SIZE(ensref,1).NE.jpm) STOP 'Inconsistent size in covariance'

        IF (.NOT.weighted_stat) THEN
          CALL ensemble_meanstd ( ensref, meanref, stdref )
          IF (stdref.GT.0.) THEN
            CALL ensemble_covariance ( ens, ensref, representer )
            representer = representer / stdref / stdref
          ELSE
            representer = correlation_missing_value
          ENDIF
        ELSE
          CALL ensemble_meanstd ( ensref, meanref, stdref, weight )
          IF (stdref.GT.0.) THEN
            CALL ensemble_covariance ( ens, ensref, representer, weight )
            representer = representer / stdref / stdref
          ELSE
            representer = correlation_missing_value
          ENDIF
        ENDIF

        END SUBROUTINE ensemble_representer
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE ensemble_covariance ( ens, ensref, cov, weight )
!----------------------------------------------------------------------
!                  ***  ensemble_meancov  ***
! 
! ** Purpose :   compute covariance from input ensemble
! 
! ** Arguments :
!         ens    : input ensemble
!         ensref : scalar reference ensemble (with respect to which computiing covariances)
!         cov    : ensemble covariance with respect to reference variable
!         weight : weight associated to eah ensemble member
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: ens
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: ensref
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: cov
        REAL(KIND=8), DIMENSION(:), INTENT( in ), OPTIONAL :: weight

        INTEGER :: jpi,jpm,jm,allocstat
        LOGICAL :: weighted_stat
        REAL(KIND=8), DIMENSION(:), allocatable :: mean
        REAL(KIND=8) :: weightsum
        REAL(KIND=8) :: meanref

        weighted_stat = PRESENT(weight)

        jpi = SIZE(ens,1)  ! Size of state vector
        jpm = SIZE(ens,2)  ! Size of ensemble

        IF (jpm.LT.2) STOP 'Ensemble size smaller than 2 in covariance'

        IF (SIZE(ensref,1).NE.jpm) STOP 'Inconsistent size in covariance'
        IF (SIZE(cov,1).NE.jpi) STOP 'Inconsistent size in covariance'
        cov = 0.

        allocate( mean(jpi), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in covariance'
        mean = 0.

        IF (weighted_stat) THEN
          IF (SIZE(weight,1).NE.jpm) STOP 'Inconsistent size in covariance'
          weightsum = 0.
        ENDIF

        IF (.NOT.weighted_stat) THEN
          DO jm=1,jpm
            CALL update_meancov_vector( ens(:,jm), ensref(jm), jm, mean, meanref, cov )
          ENDDO
          cov = cov / ( jpm -1 )
        ELSE
          DO jm=1,jpm
            CALL update_meancov_vector_weight( ens(:,jm), ensref(jm), weight(jm), weightsum, mean, meanref, cov )
          ENDDO
          cov = cov / weightsum
        ENDIF

        deallocate(mean)
        
        END SUBROUTINE ensemble_covariance
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE update_meancov_vector ( vct, varref, memidx, mean, meanref, mproda )
!----------------------------------------------------------------------
!                  ***  update_meancov_vector  ***
! 
! ** Purpose :   update mean and mean squared anomalies
! 
! ** Arguments :
!         vct    : new vector
!         varref : new reference variable (with respect to which computiing covariances)
!         memidx : index of new member
!         mean   : ensemble mean
!         meanref: reference ensemble mean
!         mproda : ensemble mean product of anomalies (with respect to the mean)
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: vct
        REAL(KIND=8), INTENT( in ) :: varref
        INTEGER, INTENT( in ) :: memidx
        REAL(KIND=8), DIMENSION(:), INTENT( inout ) :: mean
        REAL(KIND=8), INTENT( inout ) :: meanref
        REAL(KIND=8), DIMENSION(:), INTENT( inout ) :: mproda

        INTEGER :: jpi,ji
        REAL(KIND=8) :: meansav, meantmp

        jpi = SIZE(vct,1)  ! Size of state vector

        IF (SIZE(mean,1).NE.jpi) STOP 'Inconsistent size in covariance'
        IF (SIZE(mproda,1).NE.jpi) STOP 'Inconsistent size in covariance'

        meansav = meanref

        DO ji=1,jpi
          meantmp = meansav
          CALL update_meancov_variable( vct(ji), varref, memidx, mean(ji), meantmp, mproda(ji) )
        ENDDO

        meanref = meantmp

        END SUBROUTINE update_meancov_vector
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE update_meancov_vector_weight ( vct, varref, weight, weightsum, mean, meanref, mproda )
!----------------------------------------------------------------------
!                  ***  update_meancov_vector_weight  ***
! 
! ** Purpose :   update mean and mean squared anomalies (with unequal sample weight)
! 
! ** Arguments :
!         vct   : additional vector
!         varref : new reference variable (with respect to which computiing covariances)
!         weight    : weight for new value
!         weightsum : current sum of weights
!         mean  : ensemble mean
!         meanref: reference ensemble mean
!         mproda : ensemble mean squared anomalies (with respect to the mean)
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: vct
        REAL(KIND=8), INTENT( in ) :: varref
        REAL(KIND=8), INTENT( in ) :: weight
        REAL(KIND=8), INTENT( inout ) :: weightsum
        REAL(KIND=8), DIMENSION(:), INTENT( inout ) :: mean
        REAL(KIND=8), INTENT( inout ) :: meanref
        REAL(KIND=8), DIMENSION(:), INTENT( inout ) :: mproda

        INTEGER :: jpi,ji
        REAL(KIND=8) :: meansav, meantmp

        jpi = SIZE(vct,1)  ! Size of state vector

        IF (SIZE(mean,1).NE.jpi) STOP 'Inconsistent size in covariance'
        IF (SIZE(mproda,1).NE.jpi) STOP 'Inconsistent size in covariance'

        meansav = meanref

        DO ji=1,jpi
          meantmp = meansav
          CALL update_meancov_variable_weight( vct(ji), varref, weight, weightsum, &
                                             & mean(ji), meantmp, mproda(ji) )
        ENDDO

        meanref = meantmp

        END SUBROUTINE update_meancov_vector_weight
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE update_meancov_variable ( var1, var2, idx, mean1, mean2, mproda )
!----------------------------------------------------------------------
!                  ***  update_meancov_variable  ***
! 
! ** Purpose :   update mean and mean product of  anomalies
!                (using Welford's incremental algorithm)
! 
! ** Arguments :
!         var1   : additional variable value
!         var2   : additional variable value
!         idx    : index of new member
!         mean1  : ensemble mean
!         mean2  : ensemble mean
!         mproda : ensemble mean squared anomalies (with respect to the mean)
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), INTENT( in ) :: var1, var2
        INTEGER, INTENT( in ) :: idx
        REAL(KIND=8), INTENT( inout ) :: mean1, mean2
        REAL(KIND=8), INTENT( inout ) :: mproda

        REAL(KIND=8) :: misfit1, misfit2, misfit3

        misfit1 = var1 - mean1
        mean1 = mean1 + misfit1 / idx

        misfit2 = var2 - mean2
        mean2 = mean2 + misfit2 / idx

        misfit3 = var2 - mean2
        mproda = mproda + misfit1 * misfit3

        END SUBROUTINE update_meancov_variable
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE update_meancov_variable_weight ( var1, var2, weight, weightsum, mean1, mean2, mproda )
!----------------------------------------------------------------------
!                  ***  update_meancov_variable  ***
! 
! ** Purpose :   update mean and mean squared anomalies (handling unequal sample weights)
!                (using West's incremental algorithm)
! 
! ** Arguments :
!         var1      : additional variable value
!         var2      : additional variable value
!         weight    : weight for new value
!         weightsum : current sum of weights
!         mean1     : ensemble mean
!         mean2     : ensemble mean
!         mproda    : ensemble mean squared anomalies (with respect to the mean)
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), INTENT( in ) :: var1, var2
        REAL(KIND=8), INTENT( in ) :: weight
        REAL(KIND=8), INTENT( inout ) :: weightsum
        REAL(KIND=8), INTENT( inout ) :: mean1, mean2
        REAL(KIND=8), INTENT( inout ) :: mproda

        REAL(KIND=8) :: misfit1, misfit2, misfit3

        IF (weight.GT.0.) THEN

          weightsum = weightsum + weight

          misfit1 = var1 - mean1
          mean1 = mean1 + misfit1 * weight / weightsum

          misfit2 = var2 - mean2
          mean2 = mean2 + misfit2 * weight / weightsum

          misfit3 = var2 - mean2
          mproda = mproda + misfit1 * misfit3 * weight

        ENDIF

        END SUBROUTINE update_meancov_variable_weight
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE ensdam_covariance
