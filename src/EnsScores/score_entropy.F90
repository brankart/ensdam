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
!                        MODULE SCORE_ENTROPY
!
!---------------------------------------------------------------------
! Ensemble resolution score based on information theory
! by Jean-Michel Brankart, November 2018
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! events_score : compute ensemble score for the required events
! events_relative_entropy : compute relative entropy between ensemble distribution and reference distribution
! events_cross_entropy : compute cross entropy between ensemble distribution and reference distribution
! events_entropy : compute entropy of ensemble distribution for the required events
! events_probability : compute events marginal probability distributions from the ensemble
! ----------------------------------------------------------------------
MODULE ensdam_score_entropy
      IMPLICIT NONE
      PRIVATE

      PUBLIC events_score,events_relative_entropy,events_cross_entropy,events_entropy,events_probability

      ! Public variables parameterizing entropy unit
      REAL(KIND=8), PUBLIC, SAVE :: score_entropy_base=2.  ! Default = base 2 logarithm

      ! Module private definition of callback routine providing the outcome of events
      INTERFACE
        SUBROUTINE callback_events_outcome(member,outcome)
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), intent(in) :: member
        INTEGER, DIMENSION(:), intent(out) :: outcome
        END SUBROUTINE
      END INTERFACE

   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE events_score( score, ens, pref, events_outcome )
!----------------------------------------------------------------------
! ** Purpose :   compute ensemble score for the required events
! 
! ** Arguments :
!         ens    : ensemble to evaluate
!         pref   : reference probability distribution
!         events_outcome : callback routine providing the outcome of events for a given member
!         score  : ensemble score for each event
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: score
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: ens
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: pref
        PROCEDURE(callback_events_outcome) :: events_outcome

        INTEGER :: jpi,jpm,jpe,jpo,allocstat
        REAL(KIND=8), DIMENSION(:), allocatable :: entropy
        REAL(KIND=8), DIMENSION(:), allocatable :: cross_entropy

        jpi = SIZE(ens,1)  ! Size of state vector
        jpm = SIZE(ens,2)  ! Size of ensemble
        jpe = SIZE(pref,1) ! Number of events
        jpo = SIZE(pref,2) ! Maximum number of outcomes for each event

        IF (SIZE(score,1).NE.jpe) STOP 'Inconsistent size in score_entropy'

        ! Allocate entropy and cross-entropy arrays
        allocate( entropy(jpe), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in score_entropy'
        allocate( cross_entropy(jpe), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in score_entropy'

        call events_cross_entropy( cross_entropy, entropy, ens, pref, events_outcome )

        score = entropy / cross_entropy

        deallocate(entropy,cross_entropy)

        END SUBROUTINE events_score
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE events_relative_entropy( relative_entropy, ens, pref, events_outcome )
!----------------------------------------------------------------------
! ** Purpose :   compute relative entropy between ensemble distribution and reference distribution
! 
! ** Arguments :
!         ens    : ensemble to evaluate
!         pref   : reference probability distribution
!         events_outcome : callback routine providing the outcome of events for a given member
!         relative_entropy : relative entropy for each event
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: relative_entropy
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: ens
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: pref
        PROCEDURE(callback_events_outcome) :: events_outcome

        INTEGER :: jpi,jpm,jpe,jpo,allocstat
        REAL(KIND=8), DIMENSION(:), allocatable :: entropy
        REAL(KIND=8), DIMENSION(:), allocatable :: cross_entropy

        jpi = SIZE(ens,1)  ! Size of state vector
        jpm = SIZE(ens,2)  ! Size of ensemble
        jpe = SIZE(pref,1) ! Number of events
        jpo = SIZE(pref,2) ! Maximum number of outcomes for each event

        IF (SIZE(relative_entropy,1).NE.jpe) STOP 'Inconsistent size in score_entropy'

        ! Allocate entropy and cross-entropy arrays
        allocate( entropy(jpe), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in score_entropy'
        allocate( cross_entropy(jpe), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in score_entropy'

        call events_cross_entropy( cross_entropy, entropy, ens, pref, events_outcome )

        relative_entropy = cross_entropy - entropy

        deallocate(entropy,cross_entropy)

        END SUBROUTINE events_relative_entropy
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE events_cross_entropy( cross_entropy, entropy, ens, pref, events_outcome )
!----------------------------------------------------------------------
! ** Purpose :   compute cross entropy between ensemble distribution and reference distribution
! 
! ** Arguments :
!         ens    : ensemble to evaluate
!         pref   : reference probability distribution
!         events_outcome : callback routine providing the outcome of events for a given member
!         cross_entropy : cross entropy for each event
!         entropy : entropy for each event
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: cross_entropy
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: entropy
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: ens
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: pref
        PROCEDURE(callback_events_outcome) :: events_outcome

        INTEGER :: jpi,jpm,jpe,jpo,jo,allocstat
        REAL(KIND=8), DIMENSION(:,:), allocatable :: pens

        jpi = SIZE(ens,1)  ! Size of state vector
        jpm = SIZE(ens,2)  ! Size of ensemble
        jpe = SIZE(pref,1) ! Number of events
        jpo = SIZE(pref,2) ! Maximum number of outcomes for each event

        IF (SIZE(cross_entropy,1).NE.jpe) STOP 'Inconsistent size in score_entropy'
        IF (MINVAL(pref).LE.0.) STOP 'Inconsistent reference probability distribution in score_entropy'
        IF (MAXVAL(pref).GE.1.) STOP 'Inconsistent reference probability distribution in score_entropy'

        ! Allocate ensemble probability distribution
        allocate( pens(jpe,jpo), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in score_entropy'

        call events_probability( pens, ens, events_outcome )

        cross_entropy = 0. ; entropy = 0.
        DO jo=1,jpo
          WHERE(pens(:,jo).GT.0.) entropy(:) = entropy(:) - pens(:,jo) * LOG(pens(:,jo))
          cross_entropy(:) = cross_entropy(:) - pens(:,jo) * LOG(pref(:,jo))
        ENDDO

        entropy(:)       = entropy(:)       / LOG(score_entropy_base)
        cross_entropy(:) = cross_entropy(:) / LOG(score_entropy_base)

        deallocate(pens)

        END SUBROUTINE events_cross_entropy
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE events_entropy( entropy, ens, events_outcome, number_outcome )
!----------------------------------------------------------------------
! ** Purpose :   compute entropy between ensemble distribution and reference distribution
! 
! ** Arguments :
!         ens    : ensemble to evaluate
!         number_outcome : number of possible events outcome
!         events_outcome : callback routine providing the outcome of events for a given member
!         entropy : entropy for each event
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: entropy
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: ens
        PROCEDURE(callback_events_outcome) :: events_outcome
        INTEGER, INTENT(in) :: number_outcome

        INTEGER :: jpi,jpm,jpe,jpo,jo,allocstat
        REAL(KIND=8), DIMENSION(:,:), allocatable :: pens

        jpi = SIZE(ens,1)  ! Size of state vector
        jpm = SIZE(ens,2)  ! Size of ensemble
        jpe = SIZE(entropy,1) ! Number of events
        jpo = number_outcome

        ! Allocate ensemble probability distribution
        allocate( pens(jpe,jpo), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in score_entropy'

        call events_probability( pens, ens, events_outcome )

        entropy = 0.
        DO jo=1,jpo
          WHERE(pens(:,jo).GT.0.) entropy(:) = entropy(:) - pens(:,jo) * LOG(pens(:,jo))
        ENDDO

        entropy(:) = entropy(:) / LOG(score_entropy_base)

        deallocate(pens)

        END SUBROUTINE events_entropy
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE events_probability( pens, ens, events_outcome )
!----------------------------------------------------------------------
! ** Purpose :   compute events marginal probability distributions from the ensemble
! 
! ** Arguments :
!         ens    : ensemble to evaluate
!         events_outcome : callback routine providing the outcome of events for a given member
!         pens : ensemble probability distribution for each event
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT( out ) :: pens
        REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: ens
        PROCEDURE(callback_events_outcome) :: events_outcome

        INTEGER :: jpi,jpm,jpe,jpo,je,jm,allocstat
        INTEGER, DIMENSION(:), allocatable :: outcome

        jpi = SIZE(ens,1)  ! Size of state vector
        jpm = SIZE(ens,2)  ! Size of ensemble
        jpe = SIZE(pens,1) ! Number of events
        jpo = SIZE(pens,2) ! Maximum number of outcomes for each event

        ! Allocate events outcome
        allocate( outcome(jpe), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in score_entropy'

        pens =0.
        DO jm=1,jpm
          CALL events_outcome(ens(:,jm),outcome)
          IF (MINVAL(outcome).LT.1)   STOP 'Bad event outcome in score_entropy'
          IF (MAXVAL(outcome).GT.jpo) STOP 'Bad event outcome in score_entropy'
          DO je=1,jpe
            pens(je,outcome(je)) = pens(je,outcome(je)) + 1
          ENDDO
        ENDDO

        pens = pens / jpm

        deallocate(outcome)

        END SUBROUTINE events_probability
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE ensdam_score_entropy
