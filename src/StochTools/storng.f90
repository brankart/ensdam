! ----------------------------------------------------------------------
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
!                        MODULE STORNG
!
!---------------------------------------------------------------------
! Module to generate random numbers
! with Uniform, Gaussian, Gamma or Beta probability distribution
! by Jean-Michel Brankart, June 2011
! ----------------------------------------------------------------------
! The module is based on (and includes) the
! 64-bit KISS (Keep It Simple Stupid) random number generator
! distributed by George Marsaglia :
!
! http://groups.google.com/group/comp.lang.fortran/
!        browse_thread/thread/a85bf5f2a97f5a55
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! kiss          : 64-bit KISS random number generator (period ~ 2^250)
! kiss_seed     : Define seeds for KISS random number generator
! kiss_save     : Save current state of KISS (for future restart)
! kiss_load     : Load the saved state of KISS
! kiss_reset    : Reset the default seeds
! kiss_check    : Check the KISS pseudo-random sequence
! kiss_uniform  : Real random numbers with uniform distribution in [0,1]
! kiss_gaussian : Real random numbers with Gaussian distribution N(0,1)
! kiss_gamma    : Real random numbers with Gamma distribution Gamma(k,1)
! kiss_beta     : Real random numbers with Beta distribution Beta(a,b)
! kiss_sample   : Select a random sample from a set of integers
! ----------------------------------------------------------------------
MODULE ensdam_storng
  IMPLICIT NONE
  PRIVATE
! Public functions/subroutines
  PUBLIC :: kiss, kiss_seed, kiss_save, kiss_load, kiss_reset, kiss_check
  PUBLIC :: kiss_uniform, kiss_gaussian, kiss_gamma, kiss_beta, kiss_sample
! Default/initial seeds
  INTEGER(KIND=8) :: x=1234567890987654321_8
  INTEGER(KIND=8) :: y=362436362436362436_8
  INTEGER(KIND=8) :: z=1066149217761810_8
  INTEGER(KIND=8) :: w=123456123456123456_8
! Parameters to generate real random variates
  REAL(KIND=8), PARAMETER :: huge64=9223372036854775808.0_8  ! +1
  REAL(KIND=8), PARAMETER :: zero=0.0_8, half=0.5_8, one=1.0_8, two=2.0_8
! Variables to store 2 Gaussian random numbers with current index (ig)
  INTEGER(KIND=8), SAVE :: ig=1
  REAL(KIND=8), SAVE :: gran1, gran2
CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  FUNCTION kiss()
! --------------------------------------------------------------------
! The 64-bit KISS (Keep It Simple Stupid) random number generator.
! Components:
!  (1) Xorshift (XSH), period 2^64-1,
!  (2) Multiply-with-carry (MWC), period (2^121+2^63-1)
!  (3) Congruential generator (CNG), period 2^64.
! Overall period:
!  (2^250+2^192+2^64-2^186-2^129)/6 ~= 2^(247.42) or 10^(74.48)
! Set your own seeds with statement "CALL kiss_seed(ix,iy,iz,iw)".
! --------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER(KIND=8) :: kiss, t

  t = ishft(x,58) + w
  IF (s(x).eq.s(t)) THEN
    w = ishft(x,-6) + s(x)
  ELSE
    w = ishft(x,-6) + 1 - s(x+t)
  ENDIF
  x = t + x
  y = m( m( m(y,13_8), -17_8 ), 43_8 )
  z = 6906969069_8 * z + 1234567_8

  kiss = x + y + z

  CONTAINS

    FUNCTION s(k)
      INTEGER(KIND=8) :: s, k
      s = ishft(k,-63)
    END FUNCTION s

    FUNCTION m(k, n)
      INTEGER(KIND=8) :: m, k, n
      m = ieor (k, ishft(k, n) )
    END FUNCTION m

  END FUNCTION kiss
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  SUBROUTINE kiss_seed(ix, iy, iz, iw)
! --------------------------------------------------------------------
! Define seeds for KISS random number generator
! --------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER(KIND=8), INTENT(in) :: ix, iy, iz, iw

  x = ix
  y = iy
  z = iz
  w = iw

  END SUBROUTINE kiss_seed
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  SUBROUTINE kiss_reset()
! --------------------------------------------------------------------
! Reset the default seeds for KISS random number generator
! --------------------------------------------------------------------
  IMPLICIT NONE

  x=1234567890987654321_8
  y=362436362436362436_8
  z=1066149217761810_8
  w=123456123456123456_8

  END SUBROUTINE kiss_reset
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  SUBROUTINE kiss_check(check_type)
! --------------------------------------------------------------------
! Check the KISS pseudo-random sequence
! (i.e. that it reproduces the correct sequence from the default seed)
! --------------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(in) :: check_type
  INTEGER(KIND=8) :: iter, niter, correct, iran
  LOGICAL :: print_success

! Save current state of KISS
  CALL kiss_save()
! Reset the default seed
  CALL kiss_reset()

! Select check type
  SELECT CASE(check_type)
  CASE('short')
    niter = 5_8
    correct = 542381058189297533_8
    print_success = .FALSE.
  CASE('long')
    niter = 100000000_8
    correct = 1666297717051644203_8 ! Check provided by G. Marsaglia
    print_success = .TRUE.
  CASE('default')
  CASE DEFAULT
    STOP 'Bad check type in kiss_check'
  END SELECT

! Run kiss for the required number of iterations (niter)
  DO iter=1,niter
    iran = kiss()
  ENDDO

! Check that last iterate is correct
  IF (iran.NE.correct) THEN
    STOP 'Check failed: KISS internal error !!'
  ELSE
    IF (print_success) PRINT *, 'Check successful: 100 million calls to KISS OK'
  ENDIF

! Reload the previous state of KISS
  CALL kiss_load()

  END SUBROUTINE kiss_check
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  SUBROUTINE kiss_save
! --------------------------------------------------------------------
! Save current state of KISS random number generator
! --------------------------------------------------------------------
  IMPLICIT NONE

  OPEN(UNIT=30,FILE='.kiss_restart')
  WRITE(30,*) x
  WRITE(30,*) y
  WRITE(30,*) z
  WRITE(30,*) w
  CLOSE(30)

  END SUBROUTINE kiss_save
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  SUBROUTINE kiss_load
! --------------------------------------------------------------------
! Load the saved state of KISS random number generator
! --------------------------------------------------------------------
  IMPLICIT NONE
  LOGICAL :: filexists

  INQUIRE(FILE='.kiss_restart',EXIST=filexists)
  IF (filexists) THEN
    OPEN(UNIT=30,FILE='.kiss_restart')
    READ(30,*) x
    READ(30,*) y
    READ(30,*) z
    READ(30,*) w
    CLOSE(30)
  ENDIF

  END SUBROUTINE kiss_load
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  SUBROUTINE kiss_uniform(uran)
! --------------------------------------------------------------------
! Real random numbers with uniform distribution in [0,1]
! --------------------------------------------------------------------
  IMPLICIT NONE
  REAL(KIND=8), INTENT(OUT) :: uran

  uran = half * ( one + REAL(kiss(),8) / huge64 )
  
  END SUBROUTINE kiss_uniform
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  SUBROUTINE kiss_gaussian(gran)
!---------------------------------------------------------------------
! Real random numbers with Gaussian distribution N(0,1)
!---------------------------------------------------------------------
  IMPLICIT NONE
  REAL(KIND=8), INTENT(OUT) :: gran
  REAL(KIND=8) :: u1, u2, rsq, fac

! Generate 2 new Gaussian draws (gran1 and gran2)
! from 2 uniform draws on [-1,1] (u1 and u2),
! using the Marsaglia polar method
! (see Devroye, Non-Uniform Random Variate Generation, p. 235-236)
  IF (ig.EQ.1) THEN
    rsq = two
    DO WHILE ( (rsq.GE.one).OR. (rsq.EQ.zero) )
      u1 = REAL(kiss(),8) / huge64
      u2 = REAL(kiss(),8) / huge64
      rsq = u1*u1 + u2*u2
    ENDDO
    fac = SQRT(-two*LOG(rsq)/rsq)
    gran1 = u1 * fac
    gran2 = u2 * fac
  ENDIF

! Output one of the 2 draws
  IF (ig.EQ.1) THEN
    gran = gran1 ; ig = 2
  ELSE
    gran = gran2 ; ig = 1
  ENDIF

  END SUBROUTINE kiss_gaussian
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  SUBROUTINE kiss_gamma(gamr,k)
!---------------------------------------------------------------------
! Real random numbers with Gamma distribution Gamma(k,1)
!---------------------------------------------------------------------
  IMPLICIT NONE
  REAL(KIND=8), INTENT(out) :: gamr
  REAL(KIND=8), INTENT(in) :: k
  REAL(KIND=8), PARAMETER :: p1 = 4.5_8
  REAL(KIND=8), PARAMETER :: p2 = 2.50407739677627_8  ! 1+LOG(9/2)
  REAL(KIND=8), PARAMETER :: p3 = 1.38629436111989_8  ! LOG(4)
  REAL(KIND=8) :: u1, u2, b, c, d, xx, yy, zz, rr, ee
  LOGICAL :: accepted

  IF (k.GT.one) THEN
! Cheng's rejection algorithm
! (see Devroye, Non-Uniform Random Variate Generation, p. 413)
    b = k - p3 ; d = SQRT(two*k-one) ; c = k + d

    accepted=.FALSE.
    DO WHILE (.NOT.accepted)
      CALL kiss_uniform(u1)
      yy = LOG(u1/(one-u1)) / d  ! Mistake in Devroye: "* k" instead of "/ d"
      xx = k * EXP(yy)
      rr = b + c * yy - xx
      CALL kiss_uniform(u2)
      zz = u1 * u1 * u2
    
      accepted = rr .GE. (zz*p1-p2)
      IF (.NOT.accepted) accepted =  rr .GE. LOG(zz)
    ENDDO

    gamr = xx

  ELSEIF (k.LT.one) THEN
! Rejection from the Weibull density
! (see Devroye, Non-Uniform Random Variate Generation, p. 415)
    c = one/k ; d = (one-k) * EXP( (k/(one-k)) * LOG(k) )
    
    accepted=.FALSE.
    DO WHILE (.NOT.accepted)
      CALL kiss_uniform(u1)
      zz = -LOG(u1)
      xx = EXP( c * LOG(zz) )
      CALL kiss_uniform(u2)
      ee = -LOG(u2)

      accepted = (zz+ee) .GE. (d+xx)  ! Mistake in Devroye: "LE" instead of "GE"
    ENDDO

    gamr = xx

  ELSE
! Exponential distribution

    CALL kiss_uniform(u1)
    gamr = -LOG(u1)

  ENDIF

  END SUBROUTINE kiss_gamma
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  SUBROUTINE kiss_beta(betar,a,b)
!---------------------------------------------------------------------
! Real random numbers with Beta distribution Beta(a,b)
!---------------------------------------------------------------------
  IMPLICIT NONE
  REAL(KIND=8), INTENT(in) :: a,b
  REAL(KIND=8), INTENT(out) :: betar
  REAL(KIND=8) :: x,y

  CALL kiss_gamma(x,a)
  CALL kiss_gamma(y,b)

  betar = x / ( x + y )

  END SUBROUTINE kiss_beta
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  SUBROUTINE kiss_sample(a,n,k)
!---------------------------------------------------------------------
! Select a random sample of size k from a set of n integers
! The sample is output in the first k elements of a
! Set k equal to n to obtain a random permutation
!   of the whole set of integers
!---------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, DIMENSION(:), INTENT(inout) :: a
  INTEGER, INTENT(in) :: n, k
  INTEGER :: i, j, atmp
  REAL(KIND=8) :: uran

! Select the sample using the swapping method
! (see Devroye, Non-Uniform Random Variate Generation, p. 612)
  DO i=1,k
! Randomly select the swapping element between i and n (inclusive)
    CALL kiss_uniform(uran)
!   j = i + INT( REAL(n-i+1,8) * uran )
    j = i - 1 + CEILING( REAL(n-i+1,8) * uran )
! Swap elements i and j
    atmp = a(i) ; a(i) = a(j) ; a(j) = atmp
  ENDDO

  END SUBROUTINE kiss_sample
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE ensdam_storng
