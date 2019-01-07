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
!                        MODULE SCORE_OPTIMALITY
!
!---------------------------------------------------------------------
! Computation of optimality score by comparing ensemble simulation to observations
! by Jean-Michel Brankart, November 2018
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! optimality_score : compute optimality score (with option to partition the data)
! optimality_cumul : accumulate data to prepare the final computation of the score
! ----------------------------------------------------------------------
MODULE ensdam_score_optimality
      IMPLICIT NONE
      PRIVATE

      PUBLIC optimality_score, optimality_cumul

      INTERFACE optimality_score
        MODULE PROCEDURE optimality_score_global, optimality_score_partition
      END INTERFACE

      REAL(KIND=8), PUBLIC, SAVE  :: optimality_missing_value = -9999.

      ! Public definitions needed by python/julia APIs
      PUBLIC optimality_score_global, optimality_score_partition

      ! Definition for MPI
#if defined MPI
      include "mpif.h"
      INTEGER, PUBLIC, SAVE  :: mpi_comm_score_optimality=mpi_comm_world   ! definition of module global communicator
      INTEGER, SAVE :: mpi_code
#endif

      ! Module private definitions for observation error cdf
      INTERFACE
        FUNCTION callback_cdf_obs(o,y,obs_idx)   ! callback function for observation error cdf
                                                 ! compute rank of observation yo in p(yo|Hx), given y=Hx
        IMPLICIT NONE
        REAL(KIND=8), intent(in) :: o  ! observation yo
        REAL(KIND=8), intent(in) :: y  ! model equivalent y=Hx
        INTEGER, intent(in) :: obs_idx ! observation index
        REAL(KIND=8) :: callback_cdf_obs
        END FUNCTION
      END INTERFACE

      ! Public definitions needed by python/julia APIs
      PUBLIC callback_cdf_obs

   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE optimality_score_global( ens_optimality, ens_optimality_bias, ens_optimality_spread, ens, obs, cdf_obs )
!----------------------------------------------------------------------
! ** Purpose :   compute optimality score
! 
! ** Arguments :
!         ens    : ensemble to evaluate (equivalent to observation data)
!         obs  : observation data
!         ens_optimality : optimality score (should be )
!         ens_optimality_bias : bias component of optimality score (should be )
!         ens_optimality_spread : spread component of optimality score (should be 1)
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), INTENT( out ) :: ens_optimality
      REAL(KIND=8), INTENT( out ) :: ens_optimality_bias
      REAL(KIND=8), INTENT( out ) :: ens_optimality_spread
      REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: ens
      REAL(KIND=8), DIMENSION(:), INTENT( in ) :: obs
      PROCEDURE(callback_cdf_obs) :: cdf_obs

      INTEGER :: jpi,jpm,nbr,ji

      jpi = SIZE(ens,1)  ! Size of state vector
      jpm = SIZE(ens,2)  ! Size of ensemble

      IF (SIZE(obs,1).NE.jpi) STOP 'Inconsistent size in score_optimality'

      ens_optimality_bias = 0. ; ens_optimality_spread = 0. ; nbr = 0

      DO ji=1,jpi
        nbr = nbr + 1
        CALL optimality_cumul(ens(ji,:),obs(ji),nbr,ens_optimality_bias,ens_optimality_spread,cdf_obs)
      ENDDO

#if defined MPI
      ens_optimality_bias = ens_optimality_bias * nbr
      ens_optimality_spread = ens_optimality_spread * nbr
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, ens_optimality_bias, 1, MPI_DOUBLE_PRECISION,  &
                      &     MPI_SUM,mpi_comm_score_optimality,mpi_code)
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, ens_optimality_spread, 1, MPI_DOUBLE_PRECISION,  &
                      &     MPI_SUM,mpi_comm_score_optimality,mpi_code)
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, nbr, 1, MPI_INTEGER,  &
                      &     MPI_SUM,mpi_comm_score_optimality,mpi_code)
      IF (nbr.GT.0) ens_optimality_bias = ens_optimality_bias / nbr
      IF (nbr.GT.0) ens_optimality_spread = ens_optimality_spread / nbr
#endif

      IF (nbr.GT.0) ens_optimality_spread = SQRT( ens_optimality_spread )

      IF (nbr.GT.0) ens_optimality = SQRT ( ens_optimality_spread**2 + ens_optimality_bias **2 )

      END SUBROUTINE optimality_score_global
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE optimality_score_partition( ens_optimality, ens_optimality_bias, ens_optimality_spread, ens, obs, partition, cdf_obs )
!----------------------------------------------------------------------
! ** Purpose :   compute optimality score with partition of the data)
! 
! ** Arguments :
!         ens    : ensemble to evaluate (equivalent to observation data)
!         obs  : observation data
!         ens_optimality : optimality score (should be )
!         ens_optimality_bias : bias component of optimality score (should be )
!         ens_optimality_spread : spread component of optimality score (should be 1)
!         partition   : partition of observation data
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:), INTENT( out ) :: ens_optimality
      REAL(KIND=8), DIMENSION(:), INTENT( out ) :: ens_optimality_bias
      REAL(KIND=8), DIMENSION(:), INTENT( out ) :: ens_optimality_spread
      REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: ens
      REAL(KIND=8), DIMENSION(:), INTENT( in ) :: obs
      INTEGER, DIMENSION(:), INTENT( in ) :: partition
      PROCEDURE(callback_cdf_obs) :: cdf_obs

      INTEGER :: jpi,jpm,submin,submax,jsub,ji,allocstat
      REAL(KIND=8), DIMENSION(:,:), allocatable :: aa,bb
      INTEGER, DIMENSION(:), allocatable :: nbr

      jpi = SIZE(ens,1)  ! Size of state vector
      jpm = SIZE(ens,2)  ! Size of ensemble

      IF (SIZE(obs,1).NE.jpi) STOP 'Inconsistent size in score_optimality'
      IF (SIZE(partition,1).NE.jpi) STOP 'Inconsistent size in score_optimality'

      submin = MINVAL(partition)  ! Maximum index of subdomains in observation data
      submax = MAXVAL(partition)  ! Maximum index of subdomains in observation data

      IF (LBOUND(ens_optimality_bias,1).NE.submin) STOP 'Inconsistent array bounds in score_optimality'
      IF (UBOUND(ens_optimality_bias,1).NE.submax) STOP 'Inconsistent array bounds in score_optimality'
      IF (LBOUND(ens_optimality_spread,1).NE.submin) STOP 'Inconsistent array bounds in score_optimality'
      IF (UBOUND(ens_optimality_spread,1).NE.submax) STOP 'Inconsistent array bounds in score_optimality'

      allocate( nbr(submin:submax), stat=allocstat )
      IF (allocstat.NE.0) STOP 'Allocation error in score_optimality'

      ens_optimality_bias = 0. ; ens_optimality_spread = 0.

      DO ji=1,jpi
        nbr(partition(ji)) = nbr(partition(ji)) + 1
        CALL optimality_cumul(ens(ji,:),obs(ji),nbr(partition(ji)), &
             & ens_optimality_bias(partition(ji)),ens_optimality_spread(partition(ji)),cdf_obs)
      ENDDO

#if defined MPI
      ens_optimality_bias = ens_optimality_bias * nbr
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, ens_optimality_bias, submax-submin+1, MPI_DOUBLE_PRECISION,  &
                      &     MPI_SUM,mpi_comm_score_optimality,mpi_code)
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, ens_optimality_spread, submax-submin+1, MPI_DOUBLE_PRECISION,  &
                      &     MPI_SUM,mpi_comm_score_optimality,mpi_code)
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, nbr, submax-submin+1, MPI_INTEGER,  &
                      &     MPI_SUM,mpi_comm_score_optimality,mpi_code)
      WHERE (nbr.GT.0)
        ens_optimality_bias = ens_optimality_bias / nbr
      ELSEWHERE
        ens_optimality_bias = optimality_missing_value
      ENDWHERE
#endif

      WHERE (nbr.GT.0)
        ens_optimality_spread = SQRT( ens_optimality_spread / nbr )
      ELSEWHERE
        ens_optimality_spread = optimality_missing_value
      ENDWHERE

      WHERE (nbr.GT.0)
        ens_optimality = SQRT ( ens_optimality_spread**2 + ens_optimality_bias **2 )
      ELSEWHERE
        ens_optimality = optimality_missing_value
      ENDWHERE

      deallocate(nbr)

      END SUBROUTINE optimality_score_partition
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE optimality_cumul(e,o,idx,mean,sqrs,cdf_obs)
!----------------------------------------------------------------------
! ** Purpose :   Accumulate one more piece of information to compute optimality score
! 
! ** Arguments :
!         e    : input ensemble
!         o    : observation
!         idx  : index of accumulated piece of information
!         mean : mean to update
!         sqrs : square sum to update
!----------------------------------------------------------------------
      use ensdam_meanstd
      use ensdam_stoutil
      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: e
      REAL(KIND=8), INTENT(IN) :: o
      INTEGER, INTENT(IN) :: idx
      REAL(KIND=8), INTENT(INOUT) :: mean, sqrs
      PROCEDURE(callback_cdf_obs) :: cdf_obs

      INTEGER :: n, i
      REAL(KIND=8) :: rank, z, ens_mean, ens_msqr

      n = size(e)  ! ensemble size

      ens_mean = 0. ; ens_msqr = 0
      DO i=1,n
        ! compute rank of observation yo in p(yo|Hx), given ensemble member e(i)
        rank = cdf_obs(o,e(i),idx)
        ! transform this rank into N(0,1) distribution
        z = invcdf_gaussian(rank)
        ! update ensemble mean and square sum of z
        CALL update_meanstd( z, i, ens_mean, ens_msqr )
      ENDDO
      ! Compute ensemble mean square from ensemble square sum
      ens_msqr = ens_msqr / n + ens_mean * ens_mean

      ! Update global averages with contribution from this observation
      CALL update_meanstd( ens_mean, idx, mean )
      CALL update_meanstd( ens_msqr, idx, sqrs )

      END SUBROUTINE optimality_cumul
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE ensdam_score_optimality