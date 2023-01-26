!---------------------------------------------------------------------
! Copyright: CNRS - Universite de Grenoble Alpes
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
!                        MODULE STOUTIL
!
!---------------------------------------------------------------------
! Special functions required by other tools
! by Jean-Michel Brankart, September 2018
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! ----------------------------------------------------------------------
MODULE ensdam_stoutil
   IMPLICIT NONE
   PRIVATE

! Public functions/subroutines
   PUBLIC :: cdf_gaussian, pdf_gaussian, logpdf_gaussian, invcdf_gaussian
   PUBLIC :: cdf_gamma, pdf_gamma, logpdf_gamma, invcdf_gamma
   PUBLIC :: cdf_beta, pdf_beta, logpdf_beta, invcdf_beta

! Accuracy parameters
   LOGICAL, PUBLIC, save :: nominal_accuracy = .false.  ! use nominal accuracy of fundamental routines
   REAL(KIND=8), PUBLIC, save :: accuracy = 1.d-3  ! required accuracy (if not nominal accuracy)
   INTEGER, PUBLIC, save :: maxiter = 50 ! maximum number of Newton-Raphson iterations

   REAL(KIND=8), parameter :: pi=3.1415926535897932384526

   INTEGER :: iter
   REAL(KIND=8) :: g, dg, y, dy

   CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION cdf_gaussian(x)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(in) :: x
      REAL(KIND=8) :: cdf_gaussian

      REAL(KIND=8) :: gerf, a, xx

      a = 0.5
      xx = 0.5 * x * x

      if(x.lt.0.)then
        gerf=1.0-cdf_gamma(a,xx)
      elseif(x.gt.0.)then
        gerf=1.0+cdf_gamma(a,xx)
      else
        gerf=1.0
      endif

      cdf_gaussian=0.5*gerf

      END FUNCTION
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION pdf_gaussian(x)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(in) :: x
      REAL(KIND=8) :: pdf_gaussian

      REAL(KIND=8) :: xx

      xx = 0.5 * x * x
      pdf_gaussian=exp(-xx) / sqrt(pi*2.0)

      END FUNCTION
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION logpdf_gaussian(x)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(in) :: x
      REAL(KIND=8) :: logpdf_gaussian

      ! log of pdf (minus a constant)
      logpdf_gaussian = - 0.5 * x * x

      END FUNCTION
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION invcdf_gaussian(rank)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(in) :: rank
      REAL(KIND=8) :: invcdf_gaussian

      y=0.0

      IF (rank.NE.0.5) THEN
        DO iter=1,maxiter
          g = cdf_gaussian(y) - rank
          dg = pdf_gaussian(y)
          dy = g/dg
          y = y-dy
          IF (ABS(dy)/y.LT.accuracy) EXIT
          IF (iter.EQ.maxiter) print *, 'Warning: no convergence in invcdf_gaussian',rank,y
        ENDDO
      ENDIF

      invcdf_gaussian=y

      END FUNCTION
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION cdf_gamma(a,x)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(in) :: a, x
      REAL(KIND=8) :: cdf_gamma

      INTEGER :: stat

      if(x.lt.0..or.a.le.0.) stop 'bad arguments in cdf_gamma'

      if (x.eq.0.) then
        cdf_gamma=0.0
      else
        cdf_gamma=gamain(x,a,stat)
      endif

      if (stat.ne.0) then
        if (stat.ne.3) stop 'error in cdf_gamma'
      endif

      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION pdf_gamma(a,x)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(in) :: a, x
      REAL(KIND=8) :: pdf_gamma

      REAL(KIND=8) :: gln, dln
      INTEGER :: stat

      if(x.le.0..or.a.le.0.) stop 'bad arguments in pdf_gamma'

      gln=alogam(a,stat)

      if(stat.ne.0) stop 'error in pdf_gamma'

      dln=(a-1)*log(x)-x-gln

      pdf_gamma=EXP(dln)

      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION logpdf_gamma(k,theta,x)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(in) :: k,theta,x
      REAL(KIND=8) :: logpdf_gamma

      if(x.le.0..or.k.le.0..or.theta.le.0.) stop 'bad arguments in logpdf_gamma'

      ! log of pdf (minus a fucntion of k)
      logpdf_gamma=(k-1)*log(x)-x/theta-k*log(theta)

      END FUNCTION
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION invcdf_gamma(a,rank)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(in) :: a,rank
      REAL(KIND=8) :: invcdf_gamma

      if(a.le.1.) stop 'bad arguments in invpdf_gamma'

      y=a
      DO iter=1,maxiter
        g = cdf_gamma(a,y) - rank
        dg = pdf_gamma(a,y)
        dy = g/dg
        y = y-dy
        IF (ABS(dy)/y.LT.accuracy) EXIT 
        IF (iter.EQ.maxiter) print *, 'Warning: no convergence in invcdf_gamma',rank,y,dg
      ENDDO

      invcdf_gamma=y

      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION cdf_beta(a,b,x)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(in) :: a,b,x
      REAL(KIND=8) :: cdf_beta

      REAL(KIND=8) :: glna,glnb,glnab,bln
      INTEGER :: stat

      if(a.le.0.) stop 'bad argument a in cdf_beta'
      if(b.le.0.) stop 'bad argument b in cdf_beta'
      if(x.lt.0..or.x.gt.1.) stop 'bad argument in x cdf_beta'

      glna=alogam(a,stat)
      if(stat.ne.0) stop 'error in cdf_beta'
      glnb=alogam(b,stat)
      if(stat.ne.0) stop 'error in cdf_beta'
      glnab=alogam(a+b,stat)
      if(stat.ne.0) stop 'error in cdf_beta'

      bln = glna + glnb - glnab
      cdf_beta=betain(x,a,b,bln,stat)
      
      if(stat.ne.0) stop 'error in cdf_beta'

      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION pdf_beta(a,b,x)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(in) :: a,b,x
      REAL(KIND=8) :: pdf_beta

      pdf_beta=exp(logpdf_beta(a,b,x))

      END FUNCTION
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION logpdf_beta(a,b,x)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(in) :: a,b,x
      REAL(KIND=8) :: logpdf_beta

      REAL(KIND=8) :: glna,glnb,glnab,bln
      INTEGER :: stat

      if(a.le.0.) stop 'bad argument a in logpdf_beta'
      if(b.le.0.) stop 'bad argument b in logpdf_beta'
      if(x.lt.0..or.x.gt.1.) stop 'bad argument x in logpdf_beta'

      if(x.eq.0..and.a.lt.1.) stop 'bad argument in logpdf_beta'
      if(x.eq.1..and.b.lt.1.) stop 'bad argument in logpdf_beta'

      glna=alogam(a,stat)
      if(stat.ne.0) stop 'error in logpdf_beta'
      glnb=alogam(b,stat)
      if(stat.ne.0) stop 'error in logpdf_beta'
      glnab=alogam(a+b,stat)
      if(stat.ne.0) stop 'error in logpdf_beta'

      bln = glna + glnb - glnab

      if(x.eq.0.)then
        if(a.eq.1.)then
          logpdf_beta=(b-1)*log(1.-x)-bln
        else
          logpdf_beta=0.
        endif
      elseif(x.eq.1.)then
        if(b.eq.1.)then
          logpdf_beta=(a-1.)*log(x)-bln
        else
          logpdf_beta=0.
        endif
      else
        logpdf_beta=(a-1.)*log(x)+(b-1)*log(1.-x)-bln
      endif

      END FUNCTION
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION invcdf_beta(a,b,rank)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(in) :: a,b,rank
      REAL(KIND=8) :: invcdf_beta

      REAL(KIND=8) :: glna,glnb,glnab,bln
      INTEGER :: stat

      IF (a.eq.0.) THEN
        invcdf_beta=0.
      ELSEIF (b.eq.0.) THEN
        invcdf_beta=1.
      ELSE
        glna=alogam(a,stat)
        if(stat.ne.0) stop 'error in invcdf_beta'
        glnb=alogam(b,stat)
        if(stat.ne.0) stop 'error in invcdf_beta'
        glnab=alogam(a+b,stat)
        if(stat.ne.0) stop 'error in invcdf_beta'

        bln = glna + glnb - glnab

        invcdf_beta=xinbta ( a, b, bln, rank, stat )
        if(stat.ne.0) stop 'error in invcdf_beta'
      ENDIF

      END FUNCTION
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! Original code by John Burkardt ditributed on page:
! http://people.sc.fsu.edu/~jburkardt/f_src
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
function alogam ( x, ifault )

!*****************************************************************************80
!
!! ALOGAM computes the logarithm of the Gamma function.
!
!  Modified:
!
!    28 March 1999
!
!  Author:
!
!    Malcolm Pike, David Hill.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Malcolm Pike, David Hill,
!    Algorithm 291:
!    Logarithm of Gamma Function,
!    Communications of the ACM,
!    Volume 9, Number 9, September 1966, page 684.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Gamma function.
!    X should be greater than 0.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error.
!    1, X <= 0.
!
!    Output, real ( kind = 8 ) ALOGAM, the logarithm of the Gamma
!    function of X.
!
  implicit none

  real ( kind = 8 ) alogam
  real ( kind = 8 ) f
  integer ( kind = 4 ) ifault
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  if ( x <= 0.0D+00 ) then
    ifault = 1
    alogam = 0.0D+00
    return
  end if

  ifault = 0
  y = x

  if ( x < 7.0D+00 ) then

    f = 1.0D+00
    z = y

    do while ( z < 7.0D+00 )
      f = f * z
      z = z + 1.0D+00
    end do

    y = z
    f = - log ( f )

  else

    f = 0.0D+00

  end if

  z = 1.0D+00 / y / y

  alogam = f + ( y - 0.5D+00 ) * log ( y ) - y &
    + 0.918938533204673D+00 + &
    ((( &
    - 0.000595238095238D+00   * z &
    + 0.000793650793651D+00 ) * z &
    - 0.002777777777778D+00 ) * z &
    + 0.083333333333333D+00 ) / y

  return
end function

function alngam ( xvalue, ifault )

!*****************************************************************************80
!
!! ALNGAM computes the logarithm of the gamma function.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Allan Macleod
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Allan Macleod,
!    Algorithm AS 245,
!    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
!    Applied Statistics,
!    Volume 38, Number 2, 1989, pages 397-402.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the Gamma function.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error occurred.
!    1, XVALUE is less than or equal to 0.
!    2, XVALUE is too big.
!
!    Output, real ( kind = 8 ) ALNGAM, the logarithm of the gamma function of X.
!
  implicit none

  real ( kind = 8 ) alngam
  real ( kind = 8 ), parameter :: alr2pi = 0.918938533204673D+00
  integer ( kind = 4 ) ifault
  real ( kind = 8 ), dimension ( 9 ) :: r1 = (/ &
    -2.66685511495D+00, &
    -24.4387534237D+00, &
    -21.9698958928D+00, &
     11.1667541262D+00, &
     3.13060547623D+00, &
     0.607771387771D+00, &
     11.9400905721D+00, &
     31.4690115749D+00, &
     15.2346874070D+00 /)
  real ( kind = 8 ), dimension ( 9 ) :: r2 = (/ &
    -78.3359299449D+00, &
    -142.046296688D+00, &
     137.519416416D+00, &
     78.6994924154D+00, &
     4.16438922228D+00, &
     47.0668766060D+00, &
     313.399215894D+00, &
     263.505074721D+00, &
     43.3400022514D+00 /)
  real ( kind = 8 ), dimension ( 9 ) :: r3 = (/ &
    -2.12159572323D+05, &
     2.30661510616D+05, &
     2.74647644705D+04, &
    -4.02621119975D+04, &
    -2.29660729780D+03, &
    -1.16328495004D+05, &
    -1.46025937511D+05, &
    -2.42357409629D+04, &
    -5.70691009324D+02 /)
  real ( kind = 8 ), dimension ( 5 ) :: r4 = (/ &
     0.279195317918525D+00, &
     0.4917317610505968D+00, &
     0.0692910599291889D+00, &
     3.350343815022304D+00, &
     6.012459259764103D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ), parameter :: xlge = 5.10D+05
  real ( kind = 8 ), parameter :: xlgst = 1.0D+30
  real ( kind = 8 ) xvalue
  real ( kind = 8 ) y

  x = xvalue
  alngam = 0.0D+00
!
!  Check the input.
!
  if ( xlgst <= x ) then
    ifault = 2
    return
  end if

  if ( x <= 0.0D+00 ) then
    ifault = 1
    return
  end if

  ifault = 0
!
!  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
!
  if ( x < 1.5D+00 ) then

    if ( x < 0.5D+00 ) then

      alngam = - log ( x )
      y = x + 1.0D+00
!
!  Test whether X < machine epsilon.
!
      if ( y == 1.0D+00 ) then
        return
      end if

    else

      alngam = 0.0D+00
      y = x
      x = ( x - 0.5D+00 ) - 0.5D+00

    end if

    alngam = alngam + x * (((( &
        r1(5)   * y &
      + r1(4) ) * y &
      + r1(3) ) * y &
      + r1(2) ) * y &
      + r1(1) ) / (((( &
                  y &
      + r1(9) ) * y &
      + r1(8) ) * y &
      + r1(7) ) * y &
      + r1(6) )

    return

  end if
!
!  Calculation for 1.5 <= X < 4.0.
!
  if ( x < 4.0D+00 ) then

    y = ( x - 1.0D+00 ) - 1.0D+00

    alngam = y * (((( &
        r2(5)   * x &
      + r2(4) ) * x &
      + r2(3) ) * x &
      + r2(2) ) * x &
      + r2(1) ) / (((( &
                  x &
      + r2(9) ) * x &
      + r2(8) ) * x &
      + r2(7) ) * x &
      + r2(6) )
!
!  Calculation for 4.0 <= X < 12.0.
!
  else if ( x < 12.0D+00 ) then

    alngam = (((( &
        r3(5)   * x &
      + r3(4) ) * x &
      + r3(3) ) * x &
      + r3(2) ) * x &
      + r3(1) ) / (((( &
                  x &
      + r3(9) ) * x &
      + r3(8) ) * x &
      + r3(7) ) * x &
      + r3(6) )
!
!  Calculation for 12.0 <= X.
!
  else

    y = log ( x )
    alngam = x * ( y - 1.0D+00 ) - 0.5D+00 * y + alr2pi

    if ( x <= xlge ) then

      x1 = 1.0D+00 / x
      x2 = x1 * x1

      alngam = alngam + x1 * ( ( &
             r4(3)   * &
        x2 + r4(2) ) * &
        x2 + r4(1) ) / ( ( &
        x2 + r4(5) ) * &
        x2 + r4(4) )

    end if

  end if

  return
end function

function gamain ( x, p, ifault )

!*****************************************************************************80
!
!! GAMAIN computes the incomplete gamma ratio.
!
!  Discussion:
!
!    A series expansion is used if P > X or X <= 1.  Otherwise, a
!    continued fraction approximation is used.
!
!  Modified:
!
!    17 January 2008
!
!  Author:
!
!    G Bhattacharjee
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    G Bhattacharjee,
!    Algorithm AS 32:
!    The Incomplete Gamma Integral,
!    Applied Statistics,
!    Volume 19, Number 3, 1970, pages 285-287.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, P, the parameters of the incomplete 
!    gamma ratio.  0 <= X, and 0 < P.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no errors.
!    1, P <= 0.
!    2, X < 0.
!    3, underflow.
!    4, error return from the Log Gamma routine.
!
!    Output, real ( kind = 8 ) GAMAIN, the value of the incomplete
!    gamma ratio.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) acu
  real ( kind = 8 ) an
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ) dif
  real ( kind = 8 ) factor
  real ( kind = 8 ) g
  real ( kind = 8 ) gamain
  real ( kind = 8 ) gin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifault
  real ( kind = 8 ), parameter :: oflo = 1.0D+37
  real ( kind = 8 ) p
  real ( kind = 8 ) pn(6)
  real ( kind = 8 ) rn
  real ( kind = 8 ) term
  real ( kind = 8 ), parameter :: uflo = 1.0D-37
  real ( kind = 8 ) x
!
! Define accuracy
!
  if (nominal_accuracy) then
    acu = 1.d-8
  else
    acu = accuracy
  endif
!
!  Check the input.
!
  if ( p <= 0.0D+00 ) then
    ifault = 1
    gamain = 0.0D+00
    return
  end if

  if ( x < 0.0D+00 ) then
    ifault = 2
    gamain = 0.0D+00
    return
  end if

  if ( x == 0.0D+00 ) then
    ifault = 0
    gamain = 0.0D+00
    return
  end if

  g = alngam ( p, ifault )

  if ( ifault /= 0 ) then
    ifault = 4
    gamain = 0.0D+00
    return
  end if

  arg = p * log ( x ) - x - g

  if ( arg < log ( uflo ) ) then
    ifault = 3
    gamain = 0.0D+00
    return
  end if

  ifault = 0
  factor = exp ( arg )
!
!  Calculation by series expansion.
!
  if ( x <= 1.0D+00 .or. x < p ) then

    gin = 1.0D+00
    term = 1.0D+00
    rn = p

    do

      rn = rn + 1.0D+00
      term = term * x / rn
      gin = gin + term

      if ( term <= acu ) then
        exit
      end if

    end do

    gamain = gin * factor / p
    return

  end if
!
!  Calculation by continued fraction.
!
  a = 1.0D+00 - p
  b = a + x + 1.0D+00
  term = 0.0D+00

  pn(1) = 1.0D+00
  pn(2) = x
  pn(3) = x + 1.0D+00
  pn(4) = x * b

  gin = pn(3) / pn(4)

  do

    a = a + 1.0D+00
    b = b + 2.0D+00
    term = term + 1.0D+00
    an = a * term
    do i = 1, 2
      pn(i+4) = b * pn(i+2) - an * pn(i)
    end do

    if ( pn(6) /= 0.0D+00 ) then

      rn = pn(5) / pn(6)
      dif = abs ( gin - rn )
!
!  Absolute error tolerance satisfied?
!
      if ( dif <= acu ) then
!
!  Relative error tolerance satisfied?
!
        if ( dif <= acu * rn ) then
          gamain = 1.0D+00 - factor * gin
          exit
        end if

      end if

      gin = rn

    end if

    do i = 1, 4
      pn(i) = pn(i+2)
    end do

    if ( oflo <= abs ( pn(5) ) ) then

      do i = 1, 4
        pn(i) = pn(i) / oflo
      end do

    end if

  end do

  return
end function

function betain ( x, p, q, beta, ifault )

!*****************************************************************************80
!
!! BETAIN computes the incomplete Beta function ratio.
!
!  Modified:
!
!    12 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    KL Majumder, GP Bhattacharjee,
!    Algorithm AS 63:
!    The incomplete Beta Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 409-411.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument, between 0 and 1.
!
!    Input, real ( kind = 8 ) P, Q, the parameters, which
!    must be positive.
!
!    Input, real ( kind = 8 ) BETA, the logarithm of the complete
!    beta function.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error.
!    nonzero, an error occurred.
!
!    Output, real ( kind = 8 ) BETAIN, the value of the incomplete
!    Beta function ratio.
!
  implicit none

  real ( kind = 8 ) ai
  real ( kind = 8 ) beta
  real ( kind = 8 ) betain
  real ( kind = 8 ) acu
  real ( kind = 8 ) cx
  integer ( kind = 4 ) ifault
  logical indx
  integer ( kind = 4 ) ns
  real ( kind = 8 ) p
  real ( kind = 8 ) pp
  real ( kind = 8 ) psq
  real ( kind = 8 ) q
  real ( kind = 8 ) qq
  real ( kind = 8 ) rx
  real ( kind = 8 ) temp
  real ( kind = 8 ) term
  real ( kind = 8 ) x
  real ( kind = 8 ) xx

!
! Define accuracy
!
  if (nominal_accuracy) then
    acu = 0.1d-14
  else
    acu = accuracy
  endif
!
  betain = x
  ifault = 0
!
!  Check the input arguments.
!
  if ( p <= 0.0D+00 .or. q <= 0.0D+00 ) then
    ifault = 1
    return
  end if

  if ( x < 0.0D+00 .or. 1.0D+00 < x ) then
    ifault = 2
    return
  end if
!
!  Special cases.
!
  if ( x == 0.0D+00 .or. x == 1.0D+00 ) then
    return
  end if
!
!  Change tail if necessary and determine S.
!
  psq = p + q
  cx = 1.0D+00 - x

  if ( p < psq * x ) then
    xx = cx
    cx = x
    pp = q
    qq = p
    indx = .true.
  else
    xx = x
    pp = p
    qq = q
    indx = .false.
  end if

  term = 1.0D+00
  ai = 1.0D+00
  betain = 1.0D+00
  ns = int ( qq + cx * psq )
!
!  Use Soper's reduction formula.
!
  rx = xx / cx
  temp = qq - ai
  if ( ns == 0 ) then
    rx = xx
  end if

  do

    term = term * temp * rx / ( pp + ai )
    betain = betain + term
    temp = abs ( term )

    if ( temp <= acu .and. temp <= acu * betain ) then

      betain = betain * exp ( pp * log ( xx ) &
        + ( qq - 1.0D+00 ) * log ( cx ) - beta ) / pp

      if ( indx ) then
        betain = 1.0D+00 - betain
      end if

      exit

    end if

    ai = ai + 1.0D+00
    ns = ns - 1

    if ( 0 <= ns ) then
      temp = qq - ai
      if ( ns == 0 ) then
        rx = xx
      end if
    else
      temp = psq
      psq = psq + 1.0D+00
    end if

  end do

  return
end function

function xinbta ( p, q, beta, alpha, ifault )

!*****************************************************************************80
!
!! XINBTA computes the inverse of the incomplete Beta function.
!
!  Discussion:
!
!    The accuracy exponent SAE was loosened from -37 to -30, because
!    the code would not otherwise accept the results of an iteration
!    with p = 0.3, q = 3.0, alpha = 0.2.
!
!  Modified:
!
!    25 September 2014
!
!  Author:
!
!    Original FORTRAN77 version by GW Cran, KJ Martin, GE Thomas.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    GW Cran, KJ Martin, GE Thomas,
!    Remark AS R19 and Algorithm AS 109:
!    A Remark on Algorithms AS 63: The Incomplete Beta Integral
!    and AS 64: Inverse of the Incomplete Beta Integeral,
!    Applied Statistics,
!    Volume 26, Number 1, 1977, pages 111-114.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, Q, the parameters of the incomplete
!    Beta function.
!
!    Input, real ( kind = 8 ) BETA, the logarithm of the value of
!    the complete Beta function.
!
!    Input, real ( kind = 8 ) ALPHA, the value of the incomplete Beta
!    function.  0 <= ALPHA <= 1.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Output, real ( kind = 8 ) XINBTA, the argument of the incomplete
!    Beta function which produces the value ALPHA.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) SAE, requests an accuracy of about 10^SAE.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) acu
  real ( kind = 8 ) adj
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) fpu
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) iex
  integer ( kind = 4 ) ifault
  logical indx
  real ( kind = 8 ) p
  real ( kind = 8 ) pp
  real ( kind = 8 ) prev
  real ( kind = 8 ) q
  real ( kind = 8 ) qq
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ), parameter :: sae = -30.0D+00
  real ( kind = 8 ) sq
  real ( kind = 8 ) t
  real ( kind = 8 ) tx
  real ( kind = 8 ) w
  real ( kind = 8 ) value
  real ( kind = 8 ) xin
  real ( kind = 8 ) xinbta
  real ( kind = 8 ) y
  real ( kind = 8 ) yprev

  fpu = 10.0D+00 ** sae

  ifault = 0
  value = alpha
!
!  Test for admissibility of parameters.
!
  if ( p <= 0.0D+00 ) then
    ifault = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XINBTA - Fatal error!'
    write ( *, '(a)' ) '  P <= 0.0'
    stop 1
  end if

  if ( q <= 0.0D+00 ) then
    ifault = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XINBTA - Fatal error!'
    write ( *, '(a)' ) '  Q <= 0.0'
    stop 1
  end if

  if ( alpha < 0.0D+00 .or. 1.0D+00 < alpha ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XINBTA - Fatal error!'
    write ( *, '(a)' ) '  ALPHA not between 0 and 1.'
    ifault = 2
    stop 1
  end if
!
!  Return immediately if the answer is easy to determine.  
!
  if ( alpha == 0.0D+00 ) then
    value = 0.0D+00
    xinbta = value
    return
  end if

  if ( alpha == 1.0D+00 ) then
    value = 1.0D+00
    xinbta = value
    return
  end if
!
!  Change tail if necessary.
!
  if ( 0.5D+00 < alpha ) then
    a = 1.0D+00 - alpha
    pp = q
    qq = p
    indx = .true.
  else
    a = alpha
    pp = p
    qq = q
    indx = .false.
  end if
!
!  Calculate the initial approximation.
!
  r = sqrt ( - log ( a * a ) )

  y = r - ( 2.30753D+00 + 0.27061D+00 * r ) &
    / ( 1.0D+00 + ( 0.99229D+00 + 0.04481D+00 * r ) * r )

  if ( 1.0D+00 < pp .and. 1.0D+00 < qq ) then

    r = ( y * y - 3.0D+00 ) / 6.0D+00
    s = 1.0D+00 / ( pp + pp - 1.0D+00 )
    t = 1.0D+00 / ( qq + qq - 1.0D+00 )
    h = 2.0D+00 / ( s + t )
    w = y * sqrt ( h + r ) / h - ( t - s ) &
      * ( r + 5.0D+00 / 6.0D+00 - 2.0D+00 / ( 3.0D+00 * h ) )
    value = pp / ( pp + qq * exp ( w + w ) )

  else

    r = qq + qq
    t = 1.0D+00 / ( 9.0D+00 * qq )
    t = r * ( 1.0D+00 - t + y * sqrt ( t ) ) ** 3

    if ( t <= 0.0D+00 ) then
      value = 1.0D+00 - exp ( ( log ( ( 1.0D+00 - a ) * qq ) + beta ) / qq )
    else

      t = ( 4.0D+00 * pp + r - 2.0D+00 ) / t

      if ( t <= 1.0D+00 ) then
        value = exp ( ( log ( a * pp ) + beta ) / pp )
      else
        value = 1.0D+00 - 2.0D+00 / ( t + 1.0D+00 )
      end if

    end if

  end if
!
!  Solve for X by a modified Newton-Raphson method,
!  using the function BETAIN.
!
  r = 1.0D+00 - pp
  t = 1.0D+00 - qq
  yprev = 0.0D+00
  sq = 1.0D+00
  prev = 1.0D+00

  if ( value < 0.0001D+00 ) then
    value = 0.0001D+00
  end if

  if ( 0.9999D+00 < value ) then
    value = 0.9999D+00
  end if

  if (nominal_accuracy) then
    iex = max ( - 5.0D+00 / pp / pp - 1.0D+00 / a ** 0.2D+00 - 13.0D+00, sae )
    acu = 10.0D+00 ** iex
  else
    acu = accuracy
  endif
!
!  Iteration loop.
!
  do

    y = betain ( value, pp, qq, beta, ifault )

    if ( ifault /= 0 ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'XBINTA - Fatal error!'
      write ( *, '(a,i6)' ) '  BETAIN returned IFAULT = ', ifault
      stop 1
    end if

    xin = value
    y = ( y - a ) * exp ( beta + r * log ( xin ) + t * log ( 1.0D+00 - xin ) )

    if ( y * yprev <= 0.0D+00 ) then
      prev = max ( sq, fpu )
    end if

    g = 1.0D+00

    do
!
!  Choose damping factor.
!
      do

        adj = g * y
        sq = adj * adj

        if ( sq < prev ) then

          tx = value - adj

          if ( 0.0D+00 <= tx .and. tx <= 1.0D+00 ) then
            exit
          end if

        end if

        g = g / 3.0D+00

      end do
!
!  Check whether current estimate is acceptable.
!  The change "VALUE = TX" was suggested by Ivan Ukhov.
!
      if ( prev <= acu .or. y * y <= acu ) then
        value = tx
        if ( indx ) then
          value= 1.0D+00 - value
        end if
        xinbta = value
        return
      end if

      if ( tx /= 0.0D+00 .and. tx /= 1.0D+00 ) then
        exit
      end if

      g = g / 3.0D+00

    end do

    if ( tx == value ) then
      exit
    end if

    value = tx
    yprev = y

  end do

  if ( indx ) then
    value = 1.0D+00 - value
  end if

  xinbta = value

  return
end function

!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! End of code by John Burkardt ditributed on page:
! http://people.sc.fsu.edu/~jburkardt/f_src
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE ensdam_stoutil
