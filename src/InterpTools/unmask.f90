! Copyright: CNRS - Université de Grenoble

! Contributors : Jean-Michel Brankart, Charles-Emmanuel Testut, Laurent Parent,
!                Emmanuel Cosme, Claire Lauvernet, Frédéric Castruccio

! Jean-Michel.Brankart@hmg.inpg.fr

! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software.  You can  use,
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".

!-----------------------------------------------------------------------

!                        MODULE UNMASK

!-----------------------------------------------------------------------
! Unmask 2D field, by filling masked grid point
! with a 'reasonable' average of close values

! by Jean-Michel Brankart, June 2023
!
! The code and the algorithm have been derived from the module 'mod_drown'
! in the sosie package github.com/brodeau/sosie, by Laurent Brodeau
! ----------------------------------------------------------------------
      MODULE ensdam_unmask
        IMPLICIT NONE
        PRIVATE

        PUBLIC :: unmask

        INTEGER, public :: k_ew = -1   ! east-west periodicity of the grid
                ! k_ew = -1  --> no periodicity
                ! k_ew >= 0  --> periodicity with overlap of k_ew points
        INTEGER, public :: unmask_max = 400 ! how far in terms of number of grid-point we extrapolate into mask,
                ! will normally stop before 400 iterations, when all land points have been treated !
        REAL(KIND=8), public :: unmask_spval= -9999. ! Special value defining mask

      CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE unmask(phi)
!----------------------------------------------------------------------
!                  ***  unmask  ***
!
! ** Purpose :   unmask input array
!
! ** Method : weighted box average (gaussian weight) of surounding points
!
! ** Arguments:
!       phi  :  array to unmask                           (2D array)
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT(inout) :: phi

        ! Local allocatable arrays:
        INTEGER(1), ALLOCATABLE, DIMENSION(:,:) :: mask, mask_coast
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: phitmp

        INTEGER :: allocok, jext, ni, nj, ns, ji, jj, ji0
        INTEGER :: jim, jjm, jic, jjc, jip1, jim1, jjp1, jjm1
        REAL(KIND=8) :: datmsk, summsk, zweight, argexp
        LOGICAL :: ingrid

        ! Get size of input array
        ni = SIZE(phi,1)
        nj = SIZE(phi,2)

        ! Allocate temporary arrays
        allocate(mask(ni,nj), stat=allocok)
        IF (allocok.NE.0) PRINT *, 'Allocation error in ensdam_unmask'
        allocate(mask_coast(ni,nj), stat=allocok)
        IF (allocok.NE.0) PRINT *, 'Allocation error in ensdam_unmask'
        allocate(phitmp(ni,nj), stat=allocok)
        IF (allocok.NE.0) PRINT *, 'Allocation error in ensdam_unmask'

        ! Initialize the mask
        mask = 1
        WHERE (phi == unmask_spval) mask = 0

        ! Set field to zero on mask
        phi = phi * mask

        ! Iterate to reduce the mask step by step
        ! By one row of grid points at each iteration
        DO jext = 1, unmask_max

          ! Build mask of the coast-line (belonging to land points)
          mask_coast(:,:) = 0
          DO jj = 1, nj
            DO ji = 1, ni

              ! define nearest grid point in j
              jjp1 = jj+1 ; jjm1 = jj-1
              IF ( jjp1 > nj ) jjp1 = nj
              IF ( jjm1 <  1 ) jjm1 = 1

              ! define nearest grid point in i
              ji0 = ji ; jip1 = ji+1 ; jim1 = ji-1
              IF  ( k_ew >= 0 ) THEN
                IF ( ji0  > ni - k_ew ) ji0  = ji0  - ni + 2 * k_ew
                IF ( ji0  <  1 + k_ew ) ji0  = ji0  + ni - 2 * k_ew
                IF ( jip1 > ni - k_ew ) jip1 = jip1 - ni + 2 * k_ew ! W boundary
                IF ( jip1 <  1 + k_ew ) jip1 = jip1 + ni - 2 * k_ew ! W boundary
                IF ( jim1 > ni - k_ew ) jim1 = jim1 - ni + 2 * k_ew ! W boundary
                IF ( jim1 <  1 + k_ew ) jim1 = jim1 + ni - 2 * k_ew ! E boundary
              ELSE
                IF ( jip1 > ni ) jip1 = ni
                IF ( jim1 <  1 ) jim1 = 1
              END IF

              ! define a mask for the coast line
              mask_coast(ji,jj) = MIN( mask(jim1,jj) + mask(jip1,jj) + &
                                & mask(ji0,jjm1) + mask(ji0,jjp1), 1 ) &
                                & * ( 1 - mask(ji0,jj) )
            ENDDO
          ENDDO

          ! Fill the coastline points with values from the nearest unmasked points
          ns=MIN(jext,50)  ! limit of the area of research for unmasked points
          phitmp = phi     ! phitmp is going to be computed from phi
          DO jj = 1, nj
            DO ji = 1, ni

              ! update coastal points only
              IF ( mask_coast(ji,jj) == 1 ) THEN
                summsk = 0.0
                datmsk = 0.0

                ! compute weighted average in a box center on ji,jj 
                DO jjm=-ns,ns
                  DO jim=-ns,ns
                    ! box index definition
                    jic = ji + jim ; jjc = jj + jjm ;
                    IF  ( k_ew >= 0 ) THEN
                      IF ( jic > ni - k_ew ) jic = jic - ni + 2 * k_ew ! W boundary
                      IF ( jic <  1 + k_ew ) jic = jic + ni - 2 * k_ew ! E boundary
                    ENDIF
                    ! compute weighted average
                    ingrid = (jic >= 1 .AND. jic <= ni .AND. jjc >= 1 .AND. jjc <= nj)
                    IF (ingrid) THEN
                      ! compute gaussian weight
                      argexp = - REAL( jjm**2+jim**2 , 8 ) &
                             & / REAL( jext**2 , 8 )
                      zweight = EXP(argexp) * mask(jic,jjc)
                      ! compute sum of weight, and sum of weighted field
                      summsk = summsk + zweight
                      datmsk = datmsk + zweight * phi(jic,jjc)
                    ENDIF
                  ENDDO
                ENDDO

                ! compute mean over the defined box with gaussian weight (only where data valid)
                phitmp(ji,jj) = datmsk / summsk
              ENDIF
            ENDDO
          END DO

          ! Loosing land for the next iteration:
          mask = mask + mask_coast
          ! update phi (with one more coastline unmasked)
          phi = phitmp

        ENDDO

        DEALLOCATE ( mask, mask_coast, phitmp )

        END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE ensdam_unmask

