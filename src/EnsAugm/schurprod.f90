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
!                        MODULE SCHURPROD
!
!---------------------------------------------------------------------
! Perform the Schur product of two vectors with N(0,1) marginal distribution
! and restore a marginal N(0,1) distribution
! by Jean-Michel Brankart, October 2018
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! schurprod : perform Schur product
! ----------------------------------------------------------------------
MODULE ensdam_schurprod
      use ensdam_stoanam
      IMPLICIT NONE
      PRIVATE

      PUBLIC schurprod

      ! Module public variables
      LOGICAL, PUBLIC, SAVE :: schurprod_precompute = .TRUE. ! precompute transformation
      REAL(KIND=8), PUBLIC, SAVE :: schurprod_gpmin = -12.   ! min Gaussian product in precomputed table
      REAL(KIND=8), PUBLIC, SAVE :: schurprod_gpmax = 12.    ! max Gaussian product in precomputed table
      INTEGER, PUBLIC, SAVE :: schurprod_tablesize = 10000   ! number of precomputed values

      ! Module saved private variables
      LOGICAL, SAVE :: table_precomputed = .FALSE.
      REAL(KIND=8), DIMENSION(:), SAVE, allocatable :: gptable

      INTERFACE schurprod
        MODULE PROCEDURE schurprod_vector, schurprod_variable
      END INTERFACE

      ! Public definitions needed by python/julia APIs
      PUBLIC schurprod_vector, schurprod_variable

   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE schurprod_vector( prod, vct )
!----------------------------------------------------------------------
! ** Purpose :   perform Schur product of two vectors with N(0,1) marginal distribution
! 
! ** Arguments :
!         prod     : Schur product (overwritten on 1st input vector)
!         vct      : second input vector with which performing the Schur product
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( inout ) :: prod
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: vct

        INTEGER :: jpi,ji

        jpi = SIZE(prod,1)
        IF (SIZE(vct,1).NE.jpi) STOP 'Inconsistent size in schurprod'

        DO ji = 1,jpi
          CALL schurprod_variable(prod(ji),vct(ji))
        ENDDO

        END SUBROUTINE schurprod_vector
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE schurprod_variable( prod, var )
!----------------------------------------------------------------------
! ** Purpose :   perform Schur product of two variable with N(0,1) marginal distribution
!                and restore a marginal N(0,1) distribution
! ** Arguments :
!         prod     : Schur product (overwritten on 1st input vector)
!         var      : second input variable with which performing the Schur product
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), INTENT( out ) :: prod
        REAL(KIND=8), INTENT( in ) :: var

        REAL(KIND=8) :: gprod

        gprod = prod * var

        IF (schurprod_precompute) THEN
          IF (.NOT.table_precomputed) THEN
            CALL precompute_gptable()
          ENDIF
          CALL gprod_to_gau_interp(prod,gprod)
        ELSE
          CALL gprod_to_gau(prod,gprod)
        ENDIF

        END SUBROUTINE schurprod_variable
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE gprod_to_gau_interp(prod,gprod)
!----------------------------------------------------------------------
! ** Purpose :   evaluate transformation function
!                by interpolating in precomputed table
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), INTENT( out ) :: prod
        REAL(KIND=8), INTENT( in ) :: gprod

        REAL(KIND=8) :: rr, gpmin, gpmax
        INTEGER :: jr, jpr

        jpr = schurprod_tablesize
        gpmin = schurprod_gpmin
        gpmax = schurprod_gpmax

        rr = jpr*(gprod-gpmin)/(gpmax-gpmin)

        IF (rr.LT.0.) THEN
          jr = 0 ; rr = 0.
        ELSEIF (rr.GE.jpr) THEN
          jr = jpr-1 ; rr = 1.
        ELSE
          jr = INT(rr) ; rr = rr - jr
        ENDIF
 
        prod = gptable(jr) + rr * ( gptable(jr+1) - gptable(jr) )

        END SUBROUTINE gprod_to_gau_interp
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE precompute_gptable()
!----------------------------------------------------------------------
! ** Purpose :   precompute the required transformation
!                to restore a marginal N(0,1) distribution
!----------------------------------------------------------------------
        IMPLICIT NONE

        REAL(KIND=8) :: gpval, gpmin, gpmax
        INTEGER :: jtable, jptable, allocstat

        jptable = schurprod_tablesize
        gpmin = schurprod_gpmin
        gpmax = schurprod_gpmax

        allocate ( gptable(0:jptable), stat=allocstat )
        IF (allocstat.NE.0) STOP 'Allocation error in schurprod'

        DO jtable=0,jptable
           gpval=gpmin+jtable*(gpmax-gpmin)/jptable
           CALL gprod_to_gau(gptable(jtable),gpval)
        ENDDO

        table_precomputed = .TRUE.

        END SUBROUTINE precompute_gptable
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE ensdam_schurprod
