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
!                        MODULE SCORE_RANKS
!
!---------------------------------------------------------------------
! Compute ranks of verification data in ensemble simulation,
! and produce a histogram of these ranks
! by Jean-Michel Brankart, September 2021
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! compute_ranks : compute ranks of verification data in ensemble simulation
! ----------------------------------------------------------------------
MODULE ensdam_score_ranks
      use ensdam_anaqua
      use ensdam_storng
      IMPLICIT NONE
      PRIVATE

      PUBLIC compute_ranks

      ! Definition for MPI
#if defined MPI
      include "mpif.h"
      INTEGER, PUBLIC, SAVE  :: mpi_comm_score_ranks=mpi_comm_world  ! definition of module global communicator
      INTEGER, SAVE :: mpi_code
#endif

   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE compute_ranks( ens, verif, ranks, rank_histogram )
!----------------------------------------------------------------------
!                  ***  compute_ranks  ***
!
! ** Purpose :   compute ranks of verification data in ensemble simulation
!
! ** Arguments :
!         ens      : input ensemble (equivalent to verification data)
!         verif    : verification data
!         ranks    : ranks of verification data in the ensemble
!         rank_histogram  : histogram of the ranks
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: ens
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: verif
        INTEGER, DIMENSION(:), INTENT( out ) :: ranks
        INTEGER, DIMENSION(0:), INTENT( out ), OPTIONAL :: rank_histogram

        REAL(KIND=8), DIMENSION(:), allocatable :: ens_sort
        INTEGER, DIMENSION(:), allocatable :: ens_ranks
        INTEGER :: jpi,jpm,ji,jm,allocstat
        INTEGER :: rankmin,rankmax,jrank
        LOGICAL :: do_histogram

        do_histogram = PRESENT(rank_histogram)

        jpi = SIZE(ens,1)
        jpm = SIZE(ens,2)
        IF (SIZE(verif,1).NE.jpi) STOP 'Inconsistent size in compute_ranks'
        IF (SIZE(ranks,1).NE.jpi) STOP 'Inconsistent size in compute_ranks'
        IF (do_histogram) THEN
          IF (SIZE(rank_histogram,1).NE.jpm+1) STOP 'Inconsistent size in compute_ranks'
        ENDIF

        allocate ( ens_sort(1:jpm), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in compute_ranks'
        allocate ( ens_ranks(0:jpm), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in compute_ranks'

        DO ji = 1,jpi

          ! Sort input ensemble
          ens_sort(:) = ens(ji,:)
          CALL heapsort(ens_sort)

          ! Get interval of possible ranks for verification data
          rankmax=0
          DO jm=1,jpm
            IF (verif(ji).GE.ens_sort(jm)) THEN
              rankmax=jm
            ELSE
              EXIT
            ENDIF
          ENDDO

          rankmin=rankmax
          DO jm=rankmax,1,-1
            IF (ens_sort(jm).EQ.ens_sort(rankmax)) THEN
              rankmin=jm
            ELSE
              EXIT
            ENDIF
          ENDDO

          IF (rankmin.GT.rankmax) STOP 'Inconsistent order in compute_ranks'

          ! Compute rank by random sampling among possibilities
          IF (rankmin==rankmax) THEN
            ranks(ji)=rankmin
          ELSEIF (rankmin==rankmax-1) THEN
            ranks(ji)=rankmin
          ELSE
            DO jrank=rankmin,rankmax-1
              ens_ranks(jrank) = jrank
            ENDDO
            CALL kiss_sample(ens_ranks(rankmin:rankmax-1),rankmax-rankmin,1)
            ranks(ji)=ens_ranks(rankmin)
          ENDIF

        ENDDO

        deallocate(ens_sort)

        IF (MINVAL(ranks).LT.0)   STOP 'Inconsistent rank in compute_ranks'
        IF (MAXVAL(ranks).GT.jpm) STOP 'Inconsistent rank in compute_ranks'

        IF (do_histogram) THEN

          rank_histogram(0:jpm) = 0
          DO ji=1,jpi
            rank_histogram(ranks(ji)) = rank_histogram(ranks(ji)) + 1
          ENDDO

#if defined MPI
          CALL MPI_ALLREDUCE (MPI_IN_PLACE, rank_histogram(0:jpm), jpm+1, MPI_INTEGER,  &
                        &     MPI_SUM,mpi_comm_score_ranks,mpi_code)
#endif

        ENDIF

      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE ensdam_score_ranks
