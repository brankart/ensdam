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
!                        MODULE STOANAM
!
!---------------------------------------------------------------------
! Anamorphosis of random numbers
! by Jean-Michel Brankart, September 2018
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! ----------------------------------------------------------------------
MODULE ensdam_stoanam
   use ensdam_stogprod
   use ensdam_stoutil
   IMPLICIT NONE
   PRIVATE

   PUBLIC gprod_to_gau, gau_to_gam, gau_to_beta

   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE gprod_to_gau(y,x)
!
! Transform a Gaussian product random number (product of two N(0,1) numbers)
! into a Gaussian random number, with N(0,1) distribution
!
      IMPLICIT NONE
      REAL(KIND=8), intent(out) :: y
      REAL(KIND=8), intent(in) :: x
!
      REAL(KIND=8), parameter :: pi=3.1415926535897932384526
      INTEGER, parameter :: maxiter = 50
      REAL(KIND=8), parameter :: yacc = 1.e-6
!
      INTEGER :: iter, ierr
      REAL(KIND=8) :: z, g, dg, dy, zero

! Compute Gaussian product cdf corresponding to x
      zero = 0.0
      CALL fnprod(zero,zero,zero,x,z,ierr)
      IF (ierr.NE.0) STOP 'Error in gprod_to_gau'
!
! Find the value for which the Gaussian cdf
! is equal to z using the Newton-Raphson method.
      y = 0.0
      DO iter=1,maxiter
        g = cdf_gaussian(y) - z
        dg = exp(-y*y/2.0) / sqrt(pi*2.0)
        dy = g/dg
        y = y-dy
        IF (ABS(dy).LT.yacc) EXIT
        IF (iter.EQ.maxiter) THEN
          print *, 'Warning: No convergence in gprod_to_gau'
        ENDIF
      ENDDO

      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE gau_to_gam(y,x,k)
!
! Transform a Gaussian random number (with N(0,1) pdf)
! into a Gamma random number, with GAMMA(k,1)
!
! The shape parameter k is assumed larger than 1.
! The scale parameter theta is set to 1
! (it can be subsequently scaled to any value).
!
! --- Module declaration
      IMPLICIT NONE
! --- Variable declaration 
      REAL(KIND=8), intent(out) :: y
      REAL(KIND=8), intent(in) :: x, k
!
      INTEGER, parameter :: maxiter = 50
      REAL(KIND=8), parameter :: yacc = 1.e-6
!
      INTEGER :: iter
      REAL(KIND=8) :: z, g, dg, dy
!
      IF (k.LE.1.0) STOP 'Invalid Gamma parameter'

! Compute cumulative N(0,1) pdf corresponding to x
      z = cdf_gaussian(x)

! Find the value for which the cumulative Gamma distribution is equal to z
      y = invcdf_gamma(k,z)

      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE gau_to_beta(y,x,mu,stdmax)

! Transform a Gaussian random number (with N(0,1) pdf)
! into a Beta random number, with BETA(a,b)
! with a = mu * nu, and b = (1-mu) * nu
!
! The mean mu must be between 0 and 1 (strictly).
! The sample size nu must be positive (strictly).
! The variance is equal  to mu * (1-mu) / (1+nu)
!
! The sample size nu is computed from the maximum std
! (stdmax, occuring if mu=1/2), from the formula
! nu = 1 / ( 4 * stdmax **2 )  -  1
! This maximum standard deviation must be positive
! and smaller than 1/2 (strictly).
!
      IMPLICIT NONE
! --- Variable declaration 
      REAL(KIND=8), intent(out) :: y
      REAL(KIND=8), intent(in) :: x, mu, stdmax
!
      INTEGER, parameter :: maxiter = 50
      REAL(KIND=8), parameter :: yacc = 1.e-6
!
      INTEGER :: iter
      REAL(KIND=8) :: z, g, dg, dy, a, b, nu
!
      IF (mu.LE.0.0) STOP 'Invalid Beta parameter'
      IF (mu.GE.1.0) STOP 'Invalid Beta parameter'
      IF (stdmax.LE.0.0) STOP 'Invalid Beta parameter'
      IF (stdmax.GE.0.5) STOP 'Invalid Beta parameter'
!
! Compute standard parameters of the Beta distribution
      nu = 1. / ( 4. * stdmax * stdmax )  -  1.
      a = mu * nu ; b = (1-mu) * nu
!
! Compute cumulative N(0,1) pdf corresponding to x
      z = cdf_gaussian(x)
!
! Find the value for which the cumulative Beta distribution is equal to z
      y = invcdf_beta(a,b,z)

      END SUBROUTINE
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE ensdam_stoanam
