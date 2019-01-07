! =====================================================
! Module adapted from "Alan Miller's Fortran Software"
! "https://jblevins.org/mirror/amiller/"
! This is basically file "fnprod.f90" from the website
! transformed into one single module
! =====================================================

! An algorithm to compute the CDF of the product of two normal random
! variables (which may be correlated).
! Requires the use of the module constants_NSWC from file: constant.f90,
! and the module adapt_quad from file: qxgs.f90

! This version assembled by Alan Miller
! http://users.bigpond.net.au/amiller
! e-mail:  amiller@bigpond.net.au

! Latest revision - 10 December 1997
! N.B. Function FUNCGE has been modified for the case X = 0.


MODULE ensdam_stogprod
   USE ensdam_constants_NSWC
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: fnprod

   REAL (dp), PUBLIC, SAVE :: xmuyp, xmuxp, zp, rhop, rootr

CONTAINS

!----------------------------------------------------------------
SUBROUTINE fnprod (xmux, xmuy, rho, z, answer, ier, abserr, last)
!----------------------------------------------------------------

!         SUBROUTINE TO COMPUTE THE PROBABILITY PR(XY < Z) WHERE
!         WHERE X AND Y ARE BIVARIATE NORMAL
!         WITH MEANS XMUX AND XMUY RESPECTIVELY,
!         UNIT VARIANCES, AND CORRELATION RHO

! INPUTS:

!       XMUX         MEAN OF X
!
!       XMUY         MEAN OF Y

!       RHO          CORRELATION COEFFICIENT

!       Z            POINT BELOW WHICH PROBABILITY IS TO BE COMPUTED

! OUTPUTS:

!       ANSWER       COMPUTED PROBABILITY PR(XY < Z)

!       IER          RETURN CONDITION INDICATOR
!                   -1 IF ABS(RHO) = 1 AND THE ALGORITHM FOR THE
!                         DEGENERATE DISTRIBUTION WAS USED.
!                    0 IF ABS(RHO) < 1 AND IF NO ERRORS WERE
!                         DETECTED IN QXGS (DQAGS).
!                         SEE QXGS (DQAGS) DOCUMENTATION FOR MEANING OF
!                         VALUES OF IER > O

!       ABSERR (Optional) Estimated absolute error in ANSWER.

!       LAST   (Optional) The number of sub-intervals of the range of
!                         integration used by QXGS.

! REFERENCE:
! Meeker, W.Q. & Escobar, L.A. (1994) `An algorithm to compute the CDF
!       of the product of two normal random variables', Commun. in
!       Statist.-Simula., vol.23(1), 271-280.

! AUXILIARY ROUTINES:
! DERFC & QXGS (N.B. QXGS is used instead of DQAGS)

!----------------------------------------------------------------
USE ensdam_constants_NSWC
USE ensdam_adapt_quad
IMPLICIT NONE
REAL (dp), INTENT(IN)            :: xmux, xmuy, rho, z
REAL (dp), INTENT(OUT)           :: answer
INTEGER, INTENT(OUT)             :: ier
REAL (dp), INTENT(OUT), OPTIONAL :: abserr
INTEGER, INTENT(OUT), OPTIONAL   :: last

! Local variables
REAL (dp)            :: errabs, dmrho
REAL (dp), PARAMETER :: zero = 0.0D0, one = 1.0D0
INTEGER              :: final

!     CONSTANTS THAT ONE MIGHT WANT TO CHANGE TO
!     ACHIEVE A HIGHER DEGREE OF ACCURACY FROM THE ALGORITHM

REAL (dp), PARAMETER :: xm = 8.0_dp, eps = 1.0E-10_dp
INTEGER, PARAMETER   :: limit = 100

ier = -1

!        CHECK TO SEE IF RHO IS CLOSE TO -1 OR 1

dmrho = SQRT(one - rho*rho)
IF (dmrho > zero) GO TO 40

!        CALL SPECIAL ROUTINE WHEN ABS(RHO) IS EQUAL TO ONE

CALL fprode(xmux, xmuy, rho, z, answer)
RETURN
40 rootr = one/dmrho

!       DEFINE OTHER CONSTANTS NEEDED TO COMPUTE THE INTEGRAND

xmuxp = xmux
xmuyp = xmuy
rhop = rho
zp = z

!        DO THE INTEGRATION

CALL qxgs(funcge, xmux-xm, xmux+xm, eps, eps, answer, errabs, ier, limit,  &
          final)
IF (PRESENT( abserr )) abserr = errabs
IF (PRESENT( last )) last = final

RETURN
END SUBROUTINE fnprod


!----------------------------------------------------------------
FUNCTION funcge(x) RESULT(fn_val)
!----------------------------------------------------------------

!        FUNCTION TO COMPUTE INTEGRAND FOR COMPUTING PROBABILITY
!        OVER THE NONRECTANGULAR REGION

!         THE  FOLLOWING VARIABLES ARE COMMUNICATED THROUGH MODULE

!              XMUXP MEAN OF VARIABLE ON THE Y AXIS
!              XMUYP MEAN OF VARIABLE ON THE X AXIS
!              ROOTR = ONE/DSQRT(ONE-RHOR*RHOR)
!              RHOP   CORRELATION BETWEEN VARIABLES ON X AND Y AXES
!              ZP     SPECIFIED VALUE OF Z FOR THE PRODUCT

USE ensdam_constants_NSWC
IMPLICIT NONE
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp)            :: prob, xdiff, xmucy
REAL (dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp

IF (x == zero) GO TO 99
xdiff = x - xmuxp
xmucy = xmuyp + rhop*xdiff
prob = fcdfn(SIGN(one,x)*rootr*(zp/x-xmucy))
fn_val = prob*phin(xdiff)
RETURN

! 99 fn_val = phin(-xmuxp)
99 fn_val = zero
RETURN
END FUNCTION funcge


!----------------------------------------------------------------
SUBROUTINE fprode(xmux, xmuy, rho, z, answer)
!----------------------------------------------------------------

!          SUBROUTINE TO COMPUTE THE PROBABILITY P(Y<Z) WHERE
!          Y = X(1)*X(2) AND X(1) AND X(2) ARE BIVARIATE NORMAL
!          WITH MEANS XMUX AND XMUY RESPECTIVELY, SIGMA = 1, AND
!          CORRELATION RHO =  -1 OR 1

USE ensdam_constants_NSWC
IMPLICIT NONE
REAL (dp), INTENT(IN)  :: xmux, xmuy, rho, z
REAL (dp), INTENT(OUT) :: answer

! Local variables
REAL (dp)            :: discr, dsqdis, rho1, root1, root2, xmdiff
REAL (dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp, &
                        four = 4.0_dp

rho1 = SIGN(one,rho)
xmdiff = xmuy - rho1*xmux
discr = xmdiff**2 + four*rho1*z
IF(discr < zero) GO TO 95
dsqdis = SQRT(discr)
root1 = (-xmdiff-dsqdis) / (two*rho1)
root2 = (-xmdiff+dsqdis) / (two*rho1)
IF(rho1 > zero) THEN
  answer = fcdfn(root2-xmux) - fcdfn(root1-xmux)
ELSE
  answer = fcdfn(xmux-root1) + fcdfn(root2-xmux)
END IF
RETURN

95 IF(rho1 > zero) THEN
  answer = zero
ELSE
  answer = one
END IF

RETURN
END SUBROUTINE fprode


!----------------------------------------------------------------
FUNCTION phin(x) RESULT(fn_val)
!----------------------------------------------------------------

!        STANDARD NORMAL DENSITY

USE ensdam_constants_NSWC
IMPLICIT NONE
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER  :: cval = .39894228040143_dp, half = 0.5_dp

fn_val = cval*EXP(-half*x*x)
RETURN
END FUNCTION phin


!----------------------------------------------------------------
FUNCTION fcdfn(z) RESULT(fn_val)
!----------------------------------------------------------------

!        STANDARD NORMAL CDF

USE ensdam_constants_NSWC
IMPLICIT NONE
REAL (dp), INTENT(IN) :: z
REAL (dp)             :: fn_val

! Local variables
REAL (dp)             :: zroot
REAL (dp), PARAMETER  :: half = 0.5_dp, root = .7071067811865475_dp

zroot = z*root
fn_val = half*derfc(-zroot)
RETURN
END FUNCTION fcdfn



FUNCTION derfc(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION
!-----------------------------------------------------------------------
USE ensdam_constants_NSWC
IMPLICIT NONE
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: ax, t, w
INTEGER   :: i, k
REAL (dp), PARAMETER :: a(21) = (/ .1283791670955125738961589031215D+00,  &
                                  -.3761263890318375246320529677070D+00,  &
                                   .1128379167095512573896158902931D+00,  &
                                  -.2686617064513125175943235372542D-01,  &
                                   .5223977625442187842111812447877D-02,  &
                                  -.8548327023450852832540164081187D-03,  &
                                   .1205533298178966425020717182498D-03,  &
                                  -.1492565035840625090430728526820D-04,  &
                                   .1646211436588924261080723578109D-05,  &
                                  -.1636584469123468757408968429674D-06,  &
                                   .1480719281587021715400818627811D-07,  &
                                  -.1229055530145120140800510155331D-08,  &
                                   .9422759058437197017313055084212D-10,  &
                                  -.6711366740969385085896257227159D-11,  &
                                   .4463222608295664017461758843550D-12,  &
                                  -.2783497395542995487275065856998D-13,  &
                                   .1634095572365337143933023780777D-14,  &
                                  -.9052845786901123985710019387938D-16,  &
                                   .4708274559689744439341671426731D-17,  &
                                  -.2187159356685015949749948252160D-18,  &
                                   .7043407712019701609635599701333D-20 /)
!-------------------------------

!                     ABS(X) <= 1

ax = ABS(x)
IF (ax <= 1.d0) THEN
  t = x * x
  w = a(21)
  DO i = 1, 20
    k = 21 - i
    w = t * w + a(k)
  END DO
  fn_val = 0.5D0 + (0.5D0-x*(1.d0+w))
  RETURN
END IF

!                       X < -1

IF (x <= 0.d0) THEN
  fn_val = 2.d0
  IF (x < -8.3D0) RETURN
  t = x * x
  fn_val = 2.d0 - EXP(-t) * derfc0(ax)
  RETURN
END IF

!                       X > 1

fn_val = 0.d0
IF (x > 100.d0) RETURN
t = x * x
IF (t > -dxparg(1)) RETURN
fn_val = EXP(-t) * derfc0(x)
RETURN

CONTAINS


FUNCTION derfc0(x) RESULT(fn_val)
!-----------------------------------------------------------------------
USE ensdam_constants_NSWC
IMPLICIT NONE
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

!           EVALUATION OF EXP(X**2)*ERFC(X) FOR X >= 1

!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN, VIRGINIA
!        APRIL 1992
!-------------------------------
REAL (dp)            :: t, u, v, z
REAL (dp), PARAMETER :: rpinv = .56418958354775628694807945156077259D0
REAL (dp), PARAMETER :: p0 = .16506148041280876191828601D-03,  &
                        p1 =  .15471455377139313353998665D-03,  &
                        p2 =  .44852548090298868465196794D-04,  &
                        p3 = -.49177280017226285450486205D-05,  &
                        p4 = -.69353602078656412367801676D-05,  &
                        p5 = -.20508667787746282746857743D-05,  &
                        p6 = -.28982842617824971177267380D-06,  &
                        p7 = -.17272433544836633301127174D-07,  &
                        q1 =  .16272656776533322859856317D+01,  &
                        q2 =  .12040996037066026106794322D+01,  &
                        q3 =  .52400246352158386907601472D+00,  &
                        q4 =  .14497345252798672362384241D+00,  &
                        q5 =  .25592517111042546492590736D-01,  &
                        q6 =  .26869088293991371028123158D-02,  &
                        q7 =  .13133767840925681614496481D-03
REAL (dp), PARAMETER :: r0 =  .145589721275038539045668824025D+00,  &
                        r1 = -.273421931495426482902320421863D+00,  &
                        r2 =  .226008066916621506788789064272D+00,  &
                        r3 = -.163571895523923805648814425592D+00,  &
                        r4 =  .102604312032193978662297299832D+00,  &
                        r5 = -.548023266949835519254211506880D-01,  &
                        r6 =  .241432239725390106956523668160D-01,  &
                        r7 = -.822062115403915116036874169600D-02,  &
                        r8 =  .180296241564687154310619200000D-02
REAL (dp), PARAMETER :: a0 = -.45894433406309678202825375D-03,   &
                        a1 = -.12281298722544724287816236D-01,  &
                        a2 = -.91144359512342900801764781D-01,  &
                        a3 = -.28412489223839285652511367D-01,  &
                        a4 =  .14083827189977123530129812D+01,  &
                        a5 =  .11532175281537044570477189D+01,  &
                        a6 = -.72170903389442152112483632D+01,  &
                        a7 = -.19685597805218214001309225D+01,  &
                        a8 =  .93846891504541841150916038D+01,  &
                        b1 =  .25136329960926527692263725D+02,  &
                        b2 =  .15349442087145759184067981D+03,  &
                        b3 = -.29971215958498680905476402D+03,  &
                        b4 = -.33876477506888115226730368D+04,  &
                        b5 =  .28301829314924804988873701D+04,  &
                        b6 =  .22979620942196507068034887D+05,  &
                        b7 = -.24280681522998071562462041D+05,  &
                        b8 = -.36680620673264731899504580D+05,  &
                        b9 =  .42278731622295627627042436D+05,  &
                        b10=  .28834257644413614344549790D+03,  &
                        b11=  .70226293775648358646587341D+03
REAL (dp), PARAMETER :: c0 = -.7040906288250128001000086D-04,   &
                        c1 = -.3858822461760510359506941D-02,  &
                        c2 = -.7708202127512212359395078D-01,  &
                        c3 = -.6713655014557429480440263D+00,  &
                        c4 = -.2081992124162995545731882D+01,  &
                        c5 =  .2898831421475282558867888D+01,  &
                        c6 =  .2199509380600429331650192D+02,  &
                        c7 =  .2907064664404115316722996D+01,  &
                        c8 = -.4766208741588182425380950D+02,  &
                        d1 =  .5238852785508439144747174D+02,  &
                        d2 =  .9646843357714742409535148D+03,  &
                        d3 =  .7007152775135939601804416D+04,  &
                        d4 =  .8515386792259821780601162D+04,  &
                        d5 = -.1002360095177164564992134D+06,  &
                        d6 = -.2065250031331232815791912D+06,  &
                        d7 =  .5695324805290370358175984D+06,  &
                        d8 =  .6589752493461331195697873D+06,  &
                        d9 = -.1192930193156561957631462D+07
REAL (dp), PARAMETER :: e0 = .540464821348814822409610122136D+00,  &
                        e1 = -.261515522487415653487049835220D-01, &
                        e2 = -.288573438386338758794591212600D-02, &
                        e3 = -.529353396945788057720258856000D-03
REAL (dp), PARAMETER :: s1 = .75000000000000000000D+00,   &
        s2  = -.18750000000000000000D+01, s3  = .65625000000000000000D+01,  &
        s4  = -.29531250000000000000D+02, s5  = .16242187500000000000D+03,  &
        s6  = -.10557421875000000000D+04, s7  = .79180664062500000000D+04,  &
        s8  = -.67303564453125000000D+05, s9  = .63938386230468750000D+06,  &
        s10 = -.67135305541992187500D+07, s11 = .77205601373291015625D+08
!-------------------------------
!     RPINV = 1/SQRT(PI)
!-------------------------------

!                     1 <= X <= 2

IF (x <= 2.d0) THEN
  u = ((((((p7*x + p6)*x + p5)*x + p4)*x + p3)*x + p2)*x + p1) * x + p0
  v = ((((((q7*x + q6)*x + q5)*x + q4)*x + q3)*x + q2)*x + q1) * x + 1.d0

  t = (x-3.75D0) / (x+3.75D0)
  fn_val = (((((((((u/v)*t + r8)*t + r7)*t + r6)*t + r5)*t + r4)*t + r3)*t + &
           r2)*t + r1) * t + r0
  RETURN
END IF

!                     2 < X <= 4

IF (x <= 4.d0) THEN
  z = 1.d0 / (2.5D0 + x*x)
  u = (((((((a8*z + a7)*z + a6)*z + a5)*z + a4)*z + a3)*z + a2)*z + a1) * z + a0
  v = ((((((((((b11*z + b10)*z + b9)*z + b8)*z + b7)*z + b6)*z + b5)*z +  &
      b4)*z + b3)*z + b2)*z + b1) * z + 1.d0

  t = 13.d0 * z - 1.d0
  fn_val = ((((u/v)*t + e2)*t + e1)*t + e0) / x
  RETURN
END IF

!                     4 < X < 50

IF (x < 50.d0) THEN
  z = 1.d0 / (2.5D0 + x*x)
  u = (((((((c8*z + c7)*z + c6)*z + c5)*z + c4)*z + c3)*z + c2)*z + c1) * z + &
      c0
  v = ((((((((d9*z + d8)*z + d7)*z + d6)*z + d5)*z + d4)*z + d3)*z + d2)*z +  &
      d1)*z + 1.d0

  t = 13.d0 * z - 1.d0
  fn_val = (((((u/v)*t + e3)*t + e2)*t + e1)*t + e0) / x
  RETURN
END IF

!                        X >= 50

t = (1.d0/x) ** 2
z = (((((((((((s11*t + s10)*t + s9)*t + s8)*t + s7)*t + s6)*t + s5)*t +  &
    s4)*t + s3)*t + s2)*t + s1)*t - 0.5D0) * t + 1.d0
fn_val = rpinv * (z/x)
RETURN
END FUNCTION derfc0

END FUNCTION derfc

END MODULE ensdam_stogprod
