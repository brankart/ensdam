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
!                        MODULE SCORE_RCRV
!
!---------------------------------------------------------------------
! Computation of RCRV score by comparing ensemble simulation to verification data
! by Jean-Michel Brankart, November 2018
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! rcrv_score : compute RCRV score (with option to partition the data)
! rcrv_cumul : accumulate data to prepare the final computation of the score
! ----------------------------------------------------------------------
MODULE ensdam_score_rcrv
      IMPLICIT NONE
      PRIVATE

      PUBLIC rcrv_score, rcrv_cumul

      INTERFACE rcrv_score
        MODULE PROCEDURE rcrv_score_global, rcrv_score_partition
      END INTERFACE

      LOGICAL, PUBLIC, SAVE  :: rcrv_with_anamorphosis = .FALSE. ! use anamorphosis to compute reduced variable
      INTEGER, PUBLIC, SAVE  :: rcrv_number_of_quantiles = 11    ! number of quantiles to perform anamorphosis
      REAL(KIND=8), PUBLIC, SAVE  :: rcrv_missing_value = -9999.

      ! module variable defining quantiles for anamorphosis transformation
      REAL(KIND=8), DIMENSION(:), SAVE, allocatable :: quadef, quaref, qua

      ! Public definitions needed by python/julia APIs
      PUBLIC rcrv_score_global, rcrv_score_partition

      ! Definition for MPI
#if defined MPI
      include "mpif.h"
      INTEGER, PUBLIC, SAVE  :: mpi_comm_score_rcrv=mpi_comm_world   ! definition of module global communicator
      INTEGER, SAVE :: mpi_code
#endif

   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE rcrv_score_global( ens_bias, ens_spread, ens, verif )
!----------------------------------------------------------------------
! ** Purpose :   compute RCRV score
! 
! ** Arguments :
!         ens    : ensemble to evaluate (equivalent to verification data)
!         verif  : verification data
!         ens_bias : bias component of RCRV (should be 0)
!         ens_spread : spread component of RCRV (should be 1)
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), INTENT( out ) :: ens_bias
      REAL(KIND=8), INTENT( out ) :: ens_spread
      REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: ens
      REAL(KIND=8), DIMENSION(:), INTENT( in ) :: verif

      INTEGER :: jpi,jpm,nbr,ji

      jpi = SIZE(ens,1)  ! Size of state vector
      jpm = SIZE(ens,2)  ! Size of ensemble

      IF (SIZE(verif,1).NE.jpi) STOP 'Inconsistent size in score_rcrv'

      ens_bias = 0. ; ens_spread = 0. ; nbr = 0

      DO ji=1,jpi
        nbr = nbr + 1
        CALL rcrv_cumul(ens(ji,:),verif(ji),nbr,ens_bias,ens_spread)
      ENDDO

#if defined MPI
      ens_bias = ens_bias * nbr
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, ens_bias, 1, MPI_DOUBLE_PRECISION,  &
                      &     MPI_SUM,mpi_comm_score_rcrv,mpi_code)
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, ens_spread, 1, MPI_DOUBLE_PRECISION,  &
                      &     MPI_SUM,mpi_comm_score_rcrv,mpi_code)
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, nbr, 1, MPI_INTEGER,  &
                      &     MPI_SUM,mpi_comm_score_rcrv,mpi_code)
      IF (nbr.GT.0) ens_bias = ens_bias / nbr
#endif

      IF (nbr.GT.1) ens_spread = SQRT( ens_spread / (nbr-1) )

      IF (allocated(quadef)) deallocate(quadef)
      IF (allocated(quaref)) deallocate(quaref)
      IF (allocated(qua)) deallocate(qua)

      END SUBROUTINE rcrv_score_global
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE rcrv_score_partition( ens_bias, ens_spread, ens, verif, partition )
!----------------------------------------------------------------------
! ** Purpose :   compute RCRV score with partition of the data)
! 
! ** Arguments :
!         ens    : ensemble to evaluate (equivalent to verification data)
!         verif  : verification data
!         ens_bias : bias component of RCRV (should be 0)
!         ens_spread : spread component of RCRV (should be 1)
!         partition   : partition of verification data
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:), INTENT( out ) :: ens_bias
      REAL(KIND=8), DIMENSION(:), INTENT( out ) :: ens_spread
      REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: ens
      REAL(KIND=8), DIMENSION(:), INTENT( in ) :: verif
      INTEGER, DIMENSION(:), INTENT( in ) :: partition

      INTEGER :: jpi,jpm,submin,submax,jsub,ji,allocstat
      REAL(KIND=8), DIMENSION(:,:), allocatable :: aa,bb
      INTEGER, DIMENSION(:), allocatable :: nbr

      jpi = SIZE(ens,1)  ! Size of state vector
      jpm = SIZE(ens,2)  ! Size of ensemble

      IF (SIZE(verif,1).NE.jpi) STOP 'Inconsistent size in score_rcrv'
      IF (SIZE(partition,1).NE.jpi) STOP 'Inconsistent size in score_rcrv'

      submin = MINVAL(partition)  ! Maximum index of subdomains in verification data
      submax = MAXVAL(partition)  ! Maximum index of subdomains in verification data

      IF (LBOUND(ens_bias,1).NE.submin) STOP 'Inconsistent array bounds in score_rcrv'
      IF (UBOUND(ens_bias,1).NE.submax) STOP 'Inconsistent array bounds in score_rcrv'
      IF (LBOUND(ens_spread,1).NE.submin) STOP 'Inconsistent array bounds in score_rcrv'
      IF (UBOUND(ens_spread,1).NE.submax) STOP 'Inconsistent array bounds in score_rcrv'

      allocate( nbr(submin:submax), stat=allocstat )
      IF (allocstat.NE.0) STOP 'Allocation error in score_rcrv'

      ens_bias = 0. ; ens_spread = 0. ; nbr = 0

      DO ji=1,jpi
        nbr(partition(ji)) = nbr(partition(ji)) + 1
        CALL rcrv_cumul(ens(ji,:),verif(ji),nbr(partition(ji)),ens_bias(partition(ji)),ens_spread(partition(ji)))
      ENDDO

#if defined MPI
      ens_bias = ens_bias * nbr
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, ens_bias, submax-submin+1, MPI_DOUBLE_PRECISION,  &
                      &     MPI_SUM,mpi_comm_score_rcrv,mpi_code)
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, ens_spread, submax-submin+1, MPI_DOUBLE_PRECISION,  &
                      &     MPI_SUM,mpi_comm_score_rcrv,mpi_code)
      CALL MPI_ALLREDUCE (MPI_IN_PLACE, nbr, submax-submin+1, MPI_INTEGER,  &
                      &     MPI_SUM,mpi_comm_score_rcrv,mpi_code)
      WHERE (nbr.GT.0)
        ens_bias = ens_bias / nbr
      ELSEWHERE
        ens_bias = rcrv_missing_value
      ENDWHERE
#endif

      WHERE (nbr.GT.1)
        ens_spread = SQRT( ens_spread / (nbr-1) )
      ELSEWHERE
        ens_spread = rcrv_missing_value
      ENDWHERE

      deallocate(nbr)

      IF (allocated(quadef)) deallocate(quadef)
      IF (allocated(quaref)) deallocate(quaref)
      IF (allocated(qua)) deallocate(qua)

      END SUBROUTINE rcrv_score_partition
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE rcrv_cumul(e,a,idx,mean,sqrs)
!----------------------------------------------------------------------
! ** Purpose :   Accumulate one more piece of information to compute RCRV score
! 
! ** Arguments :
!         e    : input ensemble
!         a    : verification
!         idx  : index of accumulated piece of information
!         mean : mean to update
!         sqrs : square sum to update
!----------------------------------------------------------------------
      use ensdam_anaqua
      use ensdam_anatra
      use ensdam_meanstd
      use ensdam_stoutil
      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: e
      REAL(KIND=8), INTENT(IN) :: a
      INTEGER, INTENT(IN) :: idx
      REAL(KIND=8), INTENT(INOUT) :: mean, sqrs

      INTEGER :: n, i, jq, jpq, allocstat
      REAL(KIND=8) :: xmk, xsk, y

      REAL(KIND=8) :: eps=0.001

      n = size(e)  ! ensemble size

      IF (rcrv_with_anamorphosis) THEN
        jpq = rcrv_number_of_quantiles
        IF (.NOT.allocated(quadef)) THEN
          ! Allocate quantiles arrays if needed
          allocate(quadef(jpq),stat=allocstat)
          IF (allocstat.NE.0) STOP 'Allocation error in score_rcrv'
          allocate(quaref(jpq),stat=allocstat)
          IF (allocstat.NE.0) STOP 'Allocation error in score_rcrv'
          allocate(qua(jpq),stat=allocstat)
          IF (allocstat.NE.0) STOP 'Allocation error in score_rcrv'
          ! Define quantiles arrays if needed
          DO jq=1,jpq
            quadef(jq) = real(jq-1,8)/(jpq-1)
            quaref(jq) = eps + quadef(jq) * ( 1. - 2. * eps )
            !quaref(jq) = ( 1. + quadef(jq) * (n-1) ) / (n+1)
            quaref(jq) = invcdf_gaussian(quaref(jq))
          ENDDO
        ENDIF
      ENDIF

      IF (rcrv_with_anamorphosis) THEN
        y = a
        ! compute ensemble quantiles
        CALL ens_quantiles( qua, e, quadef )
        ! transform verfication data
        CALL ana_forward( y, qua, quaref )
        ! deallocate(quadef,quaref,qua)
      ELSE
        ! compute ensemble mean (xmk) and ensemble std (xsk)
        xmk=SUM(e(:))/n
        xsk=SQRT( SUM((e(:)-xmk)*(e(:)-xmk))/(n-1) )
        ! compute reduced anomaly
        y = a - xmk
        IF (xsk.NE.0.) y = y / xsk
      ENDIF

      ! Update mean and square sum of transformed/reduced anomaly
      CALL update_meanstd( y, idx, mean, sqrs )

      END SUBROUTINE rcrv_cumul
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE ensdam_score_rcrv
