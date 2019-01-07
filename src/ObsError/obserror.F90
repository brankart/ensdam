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
!                        MODULE OBSERROR
!
!---------------------------------------------------------------------
! Obervation error probability distribution
! by Jean-Michel Brankart, November 2018
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! obserror_logpdf : compute logarithm of observation error probability density function
! obserror_cdf : compute cumulate distribution function for observation errors
! obserror_sample : sample probability distribution of observation errors
! ----------------------------------------------------------------------
MODULE ensdam_obserror
      use ensdam_storng
      use ensdam_stoutil
      IMPLICIT NONE
      PRIVATE

      PUBLIC obserror_logpdf, obserror_cdf, obserror_sample

      ! public variable for the probability distribution of observation errors (default=gaussian)
      CHARACTER(len=80), PUBLIC, SAVE :: obserror_type='gaussian'
      ! public variable for the minimum expected value of gamma distribution (default=0)
      REAL(KIND=8), PUBLIC, SAVE :: min_expected_gamma=0.
      REAL(KIND=8), PUBLIC, SAVE :: min_expected_beta=0.

      INTERFACE obserror_logpdf
        MODULE PROCEDURE obserror_logpdf_variable, obserror_logpdf_vector, obserror_logpdf_vector_homogeneous
      END INTERFACE

      INTERFACE obserror_cdf
        MODULE PROCEDURE obserror_cdf_variable, obserror_cdf_vector, obserror_cdf_vector_homogeneous
      END INTERFACE

      INTERFACE obserror_sample
        MODULE PROCEDURE obserror_sample_variable, obserror_sample_vector, obserror_sample_vector_homogeneous
      END INTERFACE

      ! Public definitions needed by python/julia APIs
      PUBLIC obserror_logpdf_variable, obserror_logpdf_vector, obserror_logpdf_vector_homogeneous
      PUBLIC obserror_cdf_variable, obserror_cdf_vector, obserror_cdf_vector_homogeneous
      PUBLIC obserror_sample_variable, obserror_sample_vector, obserror_sample_vector_homogeneous

      ! variable to save last random number obtained
      REAL(KIND=8), SAVE :: gran_saved

   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION obserror_logpdf_vector( y, x, sigma )
!----------------------------------------------------------------------
! ** Purpose : compute logarithm of observation error probability density function
!              function = - log p(y|x) + constant
!              assuming independent observation error for different vector component
!              (to make the sum of the logarithms across the components)
! 
! ** Arguments :
!         y  : observation vector
!         x  : reference vector
!         sigma : spread of observation error (different meaning according to the distribution)
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:), INTENT( in ) :: y, x, sigma
      REAL(KIND=8) :: obserror_logpdf_vector

      REAL(KIND=8) :: ologpdf
      INTEGER :: jpi, ji

      jpi = SIZE(y,1)
      IF (SIZE(x,1).NE.jpi) STOP 'Inconsistent size in obserror'
      IF (SIZE(sigma,1).NE.jpi) STOP 'Inconsistent size in obserror'

      ologpdf = 0.
      DO ji=1,jpi
        ologpdf = ologpdf + obserror_logpdf_variable( y(ji), x(ji), sigma(ji) )
      ENDDO

      obserror_logpdf_vector = ologpdf

      END FUNCTION obserror_logpdf_vector
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION obserror_logpdf_vector_homogeneous( y, x, sigma )
!----------------------------------------------------------------------
! ** Purpose : compute logarithm of observation error probability density function
!              function = - log p(y|x) + constant
!              assuming independent observation error for different vector component
!              (to make the sum of the logarithms across the components)
! 
! ** Arguments :
!         y  : observation vector
!         x  : reference vector
!         sigma : spread of observation error (different meaning according to the distribution)
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:), INTENT( in ) :: y, x
      REAL(KIND=8), INTENT( in ) :: sigma
      REAL(KIND=8) :: obserror_logpdf_vector_homogeneous

      REAL(KIND=8) :: ologpdf
      INTEGER :: jpi, ji

      jpi = SIZE(y,1)
      IF (SIZE(x,1).NE.jpi) STOP 'Inconsistent size in obserror'

      ologpdf = 0.
      DO ji=1,jpi
        ologpdf = ologpdf + obserror_logpdf_variable( y(ji), x(ji), sigma )
      ENDDO

      obserror_logpdf_vector_homogeneous = ologpdf

      END FUNCTION obserror_logpdf_vector_homogeneous
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION obserror_logpdf_variable( y, x, sigma )
!----------------------------------------------------------------------
! ** Purpose : compute logarithm of observation error probability density function
!              function = - log p(y|x) + constant
! 
! ** Arguments :
!         y  : observation value
!         x  : reference value
!         sigma : spread of observation error (different meaning according to the distribution)
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), INTENT( in ) :: y, x, sigma
      REAL(KIND=8) :: obserror_logpdf_variable

      REAL(KIND=8) :: ologpdf
      REAL(KIND=8) :: logmean, logvar
      REAL(KIND=8) :: gamma_k, gamma_theta
      REAL(KIND=8) :: mu, nu, beta_a, beta_b

      SELECT CASE(obserror_type)
        CASE('normal','gaussian')
          ! Gaussian distribution (mean->obs, std->obserror)
          ologpdf = logpdf_gaussian( (y-x)/sigma )
        CASE('lognormal')
          ! Lognormal distribution (mean->x, std->sigma)
          ! Conditions: 0 < x ; 0 < sigma
          logvar = log ( 1. + exp(2.) * sigma * sigma / x )
          logmean = log( x ) - logvar * logvar / 2
          ologpdf = logpdf_gaussian( ( log(y) - logmean ) / sqrt( logvar) )
          ologpdf = ologpdf - log(x)
        CASE('gamma')
          ! Gamma distribution (mean->x, std->sigma*x)
          ! Conditions: 0 < x ; 0 < obserror
          gamma_k = 1.0 / (sigma * sigma)
          gamma_theta = MAX(x,min_expected_gamma) * (sigma * sigma)
          ologpdf = logpdf_gamma( gamma_k , gamma_theta, y )
        CASE('beta')
          ! Beta distribution (mean->y, stdmax->obserror)
          ! std = stdmax if mean = 1/2, std < stdmax otherwise
          ! Conditions: 0 < y < 1 ; 0 < obserror < 1/2
          mu = MIN(MAX(x,min_expected_beta),1.-min_expected_beta)
          nu = 1. / ( 4. * sigma * sigma )  -  1.
          beta_a = mu * nu ; beta_b = (1-mu) * nu
          ologpdf = logpdf_beta( beta_a, beta_b, y )
        CASE DEFAULT
           STOP 'Bad observation error distribution in obserror'
      END SELECT

      obserror_logpdf_variable = - ologpdf

      END FUNCTION obserror_logpdf_variable
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION obserror_cdf_vector( y, x, sigma )
!----------------------------------------------------------------------
! ** Purpose : compute cumulate distribution function for observation errors
! 
! ** Arguments :
!         y  : observation value
!         x  : reference value
!         sigma : spread of observation error (different meaning according to the distribution)
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:), INTENT( in ) :: y, x, sigma
      REAL(KIND=8), DIMENSION(:), POINTER :: obserror_cdf_vector

      INTEGER :: jpi, ji

      jpi = SIZE(y,1)
      IF (SIZE(x,1).NE.jpi) STOP 'Inconsistent size in obserror'
      IF (SIZE(sigma,1).NE.jpi) STOP 'Inconsistent size in obserror'

      DO ji=1,jpi
        obserror_cdf_vector(ji) = obserror_cdf_variable( y(ji), x(ji), sigma(ji) )
      ENDDO

      END FUNCTION obserror_cdf_vector
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION obserror_cdf_vector_homogeneous( y, x, sigma )
!----------------------------------------------------------------------
! ** Purpose : compute cumulate distribution function for observation errors
! 
! ** Arguments :
!         y  : observation value
!         x  : reference value
!         sigma : spread of observation error (different meaning according to the distribution)
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:), INTENT( in ) :: y, x
      REAL(KIND=8), INTENT( in ) :: sigma
      REAL(KIND=8), DIMENSION(:), POINTER :: obserror_cdf_vector_homogeneous

      INTEGER :: jpi, ji

      jpi = SIZE(y,1)
      IF (SIZE(x,1).NE.jpi) STOP 'Inconsistent size in obserror'

      DO ji=1,jpi
        obserror_cdf_vector_homogeneous(ji) = obserror_cdf_variable( y(ji), x(ji), sigma )
      ENDDO

      END FUNCTION obserror_cdf_vector_homogeneous
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION obserror_cdf_variable( y, x, sigma )
!----------------------------------------------------------------------
! ** Purpose : compute cumulate distribution function for observation errors
! 
! ** Arguments :
!         y  : observation value
!         x  : reference value
!         sigma : spread of observation error (different meaning according to the distribution)
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), INTENT( in ) :: y, x, sigma
      REAL(KIND=8) :: obserror_cdf_variable

      REAL(KIND=8) :: ocdf
      REAL(KIND=8) :: logmean, logvar
      REAL(KIND=8) :: gamma_k, gamma_theta
      REAL(KIND=8) :: mu, nu, beta_a, beta_b

      SELECT CASE(obserror_type)
        CASE('normal','gaussian')
          ! Gaussian distribution (mean->obs, std->obserror)
          ocdf = cdf_gaussian( (y-x)/sigma )
        CASE('lognormal')
          ! Lognormal distribution (mean->x, std->sigma)
          ! Conditions: 0 < x ; 0 < sigma
          logvar = log ( 1. + exp(2.) * sigma * sigma / x )
          logmean = log(x) - logvar * logvar / 2
          ocdf = cdf_gaussian( ( log(y) - logmean ) / sqrt( logvar) )
        CASE('gamma')
          ! Gamma distribution (mean->x, std->sigma*x)
          ! Conditions: 0 < x ; 0 < obserror
          gamma_k = 1.0 / (sigma * sigma)
          gamma_theta = MAX(x,min_expected_gamma) * (sigma * sigma)
          ocdf = cdf_gamma( gamma_k , y / gamma_theta )
        CASE('beta')
          ! Beta distribution (mean->y, stdmax->obserror)
          ! std = stdmax if mean = 1/2, std < stdmax otherwise
          ! Conditions: 0 < y < 1 ; 0 < obserror < 1/2
          mu = MIN(MAX(x,min_expected_beta),1.-min_expected_beta)
          nu = 1. / ( 4. * sigma * sigma )  -  1.
          beta_a = mu * nu ; beta_b = (1-mu) * nu
          ocdf = cdf_beta( beta_a, beta_b, y )
        CASE DEFAULT
           STOP 'Bad observation error distribution in obserror'
      END SELECT

      obserror_cdf_variable = ocdf

      END FUNCTION obserror_cdf_variable
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION obserror_sample_vector( x, sigma, rank )
!----------------------------------------------------------------------
! ** Purpose : sample probability distribution of observation errors
! 
! ** Arguments :
!         x  : reference value
!         sigma : spread of observation error (different meaning according to the distribution)
!         rank : value of the rank to sample (default=random)
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:), INTENT( in ) :: x, sigma
      REAL(KIND=8), DIMENSION(:), INTENT( in ), optional :: rank
      REAL(KIND=8), DIMENSION(:), POINTER :: obserror_sample_vector

      INTEGER :: jpi, ji
      LOGICAL :: random_sampling

      random_sampling = .NOT. present(rank)

      jpi = SIZE(x,1)
      IF (SIZE(sigma,1).NE.jpi) STOP 'Inconsistent size in obserror'
      IF (.NOT.random_sampling) THEN
        IF (SIZE(rank,1).NE.jpi) STOP 'Inconsistent size in obserror'
      ENDIF

      IF (random_sampling) THEN
        DO ji=1,jpi
          obserror_sample_vector(ji) = obserror_sample_variable( x(ji), sigma(ji) )
        ENDDO
      ELSE
        DO ji=1,jpi
          obserror_sample_vector(ji) = obserror_sample_variable( x(ji), sigma(ji), rank(ji) )
        ENDDO
      ENDIF

      END FUNCTION obserror_sample_vector
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION obserror_sample_vector_homogeneous( x, sigma, rank, uniform_rank )
!----------------------------------------------------------------------
! ** Purpose : sample probability distribution of observation errors
! 
! ** Arguments :
!         x  : reference value
!         sigma : spread of observation error (different meaning according to the distribution)
!         rank : value of the rank to sample (default=random)
!         uniform_rank : use same rank for all component of the input vector
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:), INTENT( in ) :: x
      REAL(KIND=8), INTENT( in ) :: sigma
      REAL(KIND=8), INTENT( in ), optional :: rank
      LOGICAL, INTENT( in ), optional :: uniform_rank
      REAL(KIND=8), DIMENSION(SIZE(x,1)) :: obserror_sample_vector_homogeneous

      INTEGER :: jpi, ji
      LOGICAL :: random_sampling, use_saved_rank

      random_sampling = .NOT. present(rank)
      use_saved_rank = .FALSE.
      IF (present(uniform_rank)) THEN
        use_saved_rank = uniform_rank
      ENDIF

      jpi = SIZE(x,1)

      IF (random_sampling) THEN
        IF (use_saved_rank) THEN
          obserror_sample_vector_homogeneous(1) = obserror_sample_variable( x(1), sigma )
          DO ji=2,jpi
            obserror_sample_vector_homogeneous(ji) = obserror_sample_variable( x(ji), sigma, reuse_last_rank=.TRUE. )
          ENDDO
        ELSE
          DO ji=1,jpi
            obserror_sample_vector_homogeneous(ji) = obserror_sample_variable( x(ji), sigma )
          ENDDO
        ENDIF
      ELSE
        obserror_sample_vector_homogeneous(1) = obserror_sample_variable( x(1), sigma, rank )
        DO ji=1,jpi
          obserror_sample_vector_homogeneous(ji) = obserror_sample_variable( x(ji), sigma, reuse_last_rank=.TRUE. )
        ENDDO
      ENDIF

      END FUNCTION obserror_sample_vector_homogeneous
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION obserror_sample_variable( x, sigma, rank, reuse_last_rank )
!----------------------------------------------------------------------
! ** Purpose : sample probability distribution of observation errors
! 
! ** Arguments :
!         x  : reference value
!         sigma : spread of observation error (different meaning according to the distribution)
!         rank : value of the rank to sample (default=random)
!         reuse_last_rank : use same rank as in the last call to the function
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), INTENT( in ) :: x, sigma
      REAL(KIND=8), INTENT( in ), optional :: rank
      LOGICAL, INTENT( in ), optional :: reuse_last_rank
      REAL(KIND=8) :: obserror_sample_variable

      REAL(KIND=8) :: uran, gran, bran, xpert
      REAL(KIND=8) :: logmean, logvar
      REAL(KIND=8) :: gamma_k, gamma_theta
      REAL(KIND=8) :: mu, nu, beta_a, beta_b
      LOGICAL :: random_sampling, use_saved_rank

      random_sampling = .NOT. present(rank)
      use_saved_rank = .FALSE.
      IF (present(reuse_last_rank)) THEN
        use_saved_rank = reuse_last_rank
      ENDIF

      SELECT CASE(obserror_type)
        CASE('normal','gaussian')
          ! Gaussian distribution (mean->x, std->sigma)
          IF (use_saved_rank) THEN
            gran = gran_saved
          ELSE
            IF (random_sampling) THEN
              call kiss_gaussian(gran)
            ELSE
              gran = invcdf_gaussian(rank)
            ENDIF
            gran_saved = gran
          ENDIF
          xpert = x + sigma * gran
        CASE('lognormal')
          ! Lognormal distribution (mean->x, std->sigma)
          ! Conditions: 0 < x ; 0 < sigma
          IF (use_saved_rank) THEN
            gran = gran_saved
          ELSE
            IF (random_sampling) THEN
              call kiss_gaussian(gran)
            ELSE
              gran = invcdf_gaussian(rank)
            ENDIF
            gran_saved = gran
          ENDIF
          logvar = log ( 1. + exp(2.) * sigma * sigma / x )
          logmean = log( x ) - logvar * logvar / 2
          xpert = exp( logmean + sqrt( logvar) * gran )
        CASE('gamma')
          ! Gamma distribution (mean->x, std->sigma*x)
          ! Conditions: 0 < x ; 0 < sigma
          gamma_k = 1.0 / (sigma * sigma)
          gamma_theta = MAX(x,min_expected_gamma) * (sigma * sigma)
          IF (use_saved_rank) THEN
            gran = gran_saved
          ELSE
            IF (random_sampling) THEN
              call kiss_gamma(gran, gamma_k)
            ELSE
              gran = invcdf_gamma(gamma_k, rank)
            ENDIF
            gran_saved = gran
          ENDIF
          xpert = gamma_theta * gran
        CASE('beta')
          ! Beta distribution (mean->x, stdmax->sigma)
          ! std = stdmax if mean = 1/2, std < stdmax otherwise
          ! Conditions: 0 < x < 1 ; 0 < sigma < 1/2
          mu = MIN(MAX(x,min_expected_beta),1.-min_expected_beta)
          nu = 1. / ( 4. * sigma * sigma )  -  1.
          beta_a = mu * nu ; beta_b = (1-mu) * nu
          IF (random_sampling) THEN
            call kiss_beta(bran, beta_a, beta_b)
          ELSE
            bran = invcdf_beta(beta_a, beta_b, rank)
          ENDIF
          xpert = bran
        CASE DEFAULT
           STOP 'Bad observation error distribution in obserror'
      END SELECT

      obserror_sample_variable = xpert

      END FUNCTION obserror_sample_variable
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE ensdam_obserror
