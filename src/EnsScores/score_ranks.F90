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
! as a prerequisite to drawing rank histograms
! by Jean-Michel Brankart, September 2021
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! compute_ranks : compute ranks of verification data in ensemble simulation
! ----------------------------------------------------------------------
MODULE ensdam_score_ranks
      use ensdam_anaqua
      use ensdam_anatra
      use ensdam_storng
      IMPLICIT NONE
      PRIVATE

      PUBLIC compute_ranks

   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE compute_ranks( ens, verif, ranks )
!----------------------------------------------------------------------
!                  ***  compute_ranks  ***
!
! ** Purpose :   compute ranks of verification data in ensemble simulation
!
! ** Arguments :
!         ens      : input ensemble (equivalent to verification data)
!         verif    : verification data
!         ranks    : ranks of verification data in th ensemble
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: ens
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: verif
        INTEGER, DIMENSION(:), INTENT( out ) :: ranks

        REAL(KIND=8), DIMENSION(:), allocatable :: ens_sort
        INTEGER, DIMENSION(:), allocatable :: ens_ranks
        INTEGER :: jpi,jpm,ji,jm,allocstat
        INTEGER :: rankmin,rankmax,jrank

        jpi = SIZE(ens,1)
        jpm = SIZE(ens,2)
        IF (SIZE(verif,1).NE.jpi) STOP 'Inconsistent size in compute_ranks'
        IF (SIZE(ranks,1).NE.jpi) STOP 'Inconsistent size in compute_ranks'

        allocate ( ens_sort(1:jpm), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in compute_ranks'
        allocate ( ens_ranks(1:jpm), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in compute_ranks'

        DO ji = 1,jpi

          ! Sort input ensemble
          ens_sort(:) = ens(ji,:)
          CALL heapsort(ens_sort)

          ! Get interval of possible ranks for verification data
          CALL qualoc_fwd(rankmin,rankmax,verif(ji),ens_sort(:))

          ! Compute rank by random sampling among possibilities
          IF (rankmin==rankmax) THEN
            ranks(ji)=rankmin
          ELSE
            DO jrank=rankmin,rankmax
              ens_ranks(jrank) = jrank
            ENDDO
            CALL kiss_sample(ens_ranks(rankmin:rankmax),rankmax-rankmin+1,1)
            ranks(ji)=ens_ranks(rankmin)
          ENDIF

        ENDDO

        deallocate(ens_sort)

      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE ensdam_score_ranks
