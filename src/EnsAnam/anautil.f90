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
!                        MODULE ANAUTIL
!
!---------------------------------------------------------------------
! Utilities for anamorphic transfromation
! -Compute quantiles of target distribution
! by Jean-Michel Brankart, September 2018
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! ana_util_init : initialize module parameters
! ana_util_quaref : compute quantiles of target distribution
! ----------------------------------------------------------------------
MODULE ensdam_anautil
      use ensdam_stoutil
      IMPLICIT NONE
      PRIVATE

      PUBLIC ana_util_quaref, ana_util_init

      CHARACTER(len=80), PUBLIC, SAVE :: anautil_reference_cdf='gaussian' ! default: gaussian anamorphosis
      REAL(kind=8), PUBLIC, SAVE :: anautil_a, anautil_b ! distribution parameters

   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE ana_util_init( kref_cdf, ka, kb )
!----------------------------------------------------------------------
! ** Purpose :   initialize module parameters
!----------------------------------------------------------------------
        IMPLICIT NONE
        CHARACTER, INTENT( in ) :: kref_cdf
        REAL(kind=8) :: ka, kb

        anautil_reference_cdf = kref_cdf
        anautil_a = ka  ;  anautil_b = kb

        SELECT CASE(anautil_reference_cdf)
        CASE('gaussian','gamma','beta')
        CASE DEFAULT
           STOP 'Bad reference cdf in ana_util_init'
        END SELECT

        END SUBROUTINE ana_util_init
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        SUBROUTINE ana_util_quaref( quaref, quadef )
!----------------------------------------------------------------------
! ** Purpose : compute quantiles of reference cdf
! 
! ** Arguments :
!         quaref : quantiles of reference cdf
!         quadef : definition of quantiles (in [0,1])
!----------------------------------------------------------------------
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT( out ) :: quaref
        REAL(KIND=8), DIMENSION(:), INTENT( in ) :: quadef

        INTEGER :: jpq,jq
        REAL(KIND=8) :: rank

        jpq = SIZE(quadef,1)
        IF (SIZE(quaref,1).NE.jpq) STOP 'Inconsistent size in ana_util_quaref'

        IF (MINVAL(quadef).LT.0.0) STOP 'Bad quantile definition in ana_util_quaref'
        IF (MAXVAL(quadef).GT.1.0) STOP 'Bad quantile definition in ana_util_quaref'

        DO jq=1,jpq

          SELECT CASE(anautil_reference_cdf)
          CASE('gaussian')
            ! Gaussian distribution N(0,1)
            quaref(jq) = invcdf_gaussian(quadef(jq))
          CASE('gamma')
            ! Gamma distribution
            quaref(jq) = invcdf_gamma(anautil_a,quadef(jq))
          CASE('beta')
            ! Beta distribution
            quaref(jq) = invcdf_beta(anautil_a,anautil_b,quadef(jq))
          CASE DEFAULT
             STOP 'Bad reference cdf in ana_util_quaref'
          END SELECT

        ENDDO

        END SUBROUTINE ana_util_quaref
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE ensdam_anautil
