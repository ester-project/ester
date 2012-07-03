      SUBROUTINE CONST(METHOD,VECSTG,REQSTG,LINTPL)
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
*************************************************
C
C  Purpose:   Set formula definitions and formula characteristics for
C             selected method. Return storage requirements for the
C             selected method.
C
C  Input:     METHOD
C  Output:    VECSTG, REQSTG, LINTPL
C
C  Common:    Initializes:    /RKCOM4/ A(*,*), B(*), C(*), BHAT(*), R(*),
C                                      E(*), PTR(*), NSTAGE, METHD, INTP, MINTP
C                             /RKCOM5/ TOOSML, COST, SAFETY, EXPON, STBRAD,
C                                      TANANG, RS, RS1, RS2, RS3, RS4, ORDER,
C                                      LSTSTG, MAXTRY, NSEC, FSAL
C             Reads:          /RKCOM7/ RNDOFF
C             Alters:         none
C
C  Comments:
C  =========
C  Runge-Kutta formula pairs are described by a set of coefficients
C  and by the setting of a number of parameters that describe the
C  characteristics of the pair.  The input variable METHD indicates
C  which of the three pairs available is to be set. In the case of
C  METHD = 2 additional coefficients are defined that make interpolation
C  of the results possible and provide an additional error estimator.
C  VECSTG is the number of columns of workspace required to compute the
C  stages of a METHD. For interpolation purposes the routine returns via
C  COMMON the logical variable INTP indicating whether interpolation is
C  possible and via the call list:
C  REQSTG - whether the stages are required to form the
C           interpolant (set .FALSE. if INTP=.FALSE.)
C  LINTPL - the number of extra columns of storage required for use
C           with UT (set 0 if INTP=.FALSE.)
C
C  Quantities set in common blocks:
C  METHD - copy of METHOD
C  A, B, C, BHAT - coefficients of the selected method
C  R      - extra coefficents for interpolation with METHD = 2
C  E      - extra coefficients for additional local error estimate
C           with METHD = 2
C  PTR    - vector of pointers indicating how individual stages are to
C           be stored.  With it zero coefficients of the formulas can
C           be exploited to reduce the storage required
C  NSTAGE - number of stages for the specified METHD
C  INTP   - indicates whether there is an associated interpolant
C           (depending on the method, the user may have to supply
C           extra workspace)
C  MINTP  - the degree of the interpolating polynomial, if one exists
C  FSAL   - indicates whether the last stage of a step can be used as
C           the first stage of the following step
C  LSTSTG - pointer to location of last stage for use with FSAL=.TRUE.
C  ORDER  - the lower order of the pair of Runge-Kutta formulas that
C           constitute a METHD
C  TANANG, 
C  STBRAD - the stability region of the formula used to advance
C           the integration is approximated by a sector in the left half
C           complex plane.  TANANG is the tangent of the interior angle
C           of the sector and STBRAD is the radius of the sector.
C  COST   - cost of a successful step in function evaluations
C  MAXTRY - limit on the number of iterations in the stiffness check. As
C           set, no more than 24 function evaluations are made in the check.
C  NSEC   - each step of size H in the primary integration corresponds to
C           NSEC steps of size H/NSEC in the secondary integration when
C           global error assessment is done.
C  EXPON  - used to adjust the step size; this code implements an error
C           per step control for which EXPON = 1/(ORDER + 1).
C  SAFETY - quantity used in selecting the step size
C  TOOSML - quantity used to determine when a step size is too small for
C           the precision available
C  RS, RS1,
C  RS2, RS3,
C  RS4    - quantities used in determining the maximum and minimum change
C           change in step size (set independently of METHD)
C
C  Further comments on SAFETY:
C ============================
C  The code estimates the largest step size that will yield the specified
C  accuracy.  About half the time this step size would result in a local
C  error that is a little too large, and the step would be rejected.  For
C  this reason a SAFETY factor is used to reduce the "optimal" value to one
C  more likely to succeed.  Unfortunately, there is some art in choosing this
C  value. The more expensive a failed step, the smaller one is inclined to 
C  take this factor. However, a small factor means that if the prediction were
C  good, more accuracy than desired would be obtained and the behavior of the
C  error would then be irregular.  The more stringent the tolerance, the better
C  the prediction, except near limiting precision. Thus the general range of 
C  tolerances expected influences the choice of SAFETY.
C
C     .. Scalar Arguments ..
      INTEGER           LINTPL, METHOD, VECSTG
      LOGICAL           REQSTG
C     .. Common Block to hold Formula Definitions ..
      DOUBLE PRECISION  A(13,13), B(13), C(13), BHAT(13), R(11,6),
     &                  E(7)
      INTEGER           PTR(13), NSTAGE, METHD, MINTP
      LOGICAL           INTP
      COMMON /RKCOM4/   A, B, C, BHAT, R, E, PTR, NSTAGE, METHD,
     &                  MINTP, INTP
      SAVE   /RKCOM4/
C     .. Common Block to hold Formula Characterisitcs ..
      DOUBLE PRECISION  TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
     &                  RS, RS1, RS2, RS3, RS4
      INTEGER           ORDER, LSTSTG, MAXTRY, NSEC
      LOGICAL           FSAL
      COMMON /RKCOM5/   TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
     &                  RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG, MAXTRY,
     &                  NSEC, FSAL
      SAVE   /RKCOM5/
C     .. Common Block for Environment Parameters ..
      DOUBLE PRECISION  MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY
      INTEGER           OUTCH
      COMMON /RKCOM7/   MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY,
     &                  OUTCH
      SAVE   /RKCOM7/
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO, TWO, FIFTY, FIVEPC
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0,TWO=2.0D+0,FIFTY=50.D+0,
     &                  FIVEPC=0.05D+0)
C     .. Local Scalars ..
      DOUBLE PRECISION  CDIFF, DIFF
      INTEGER           I, J
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, INT, MAX, MIN
C     .. Executable Statements ..
C
      METHD = METHOD
C
      GO TO (20,40,100) METHD
C
C  METHD = 1.
C    This pair is from "A 3(2) Pair of Runge-Kutta Formulas" by P. Bogacki
C    and L.F. Shampine, Appl. Math. Lett., 2, pp. 321-325, 1989.  The authors
C    are grateful to P. Bogacki for his assistance in implementing the pair.
C
   20 CONTINUE
      NSTAGE = 4
      FSAL = .TRUE.
      ORDER = 2
      TANANG = 8.9D0
      STBRAD = 2.3D0
      SAFETY = 0.8D0
      INTP = .TRUE.
      MINTP = 3
      REQSTG = .FALSE.
      LINTPL = 2
      NSEC = 3
C
      PTR(1) = 0
      PTR(2) = 1
      PTR(3) = 2
      PTR(4) = 3
C
      A(2,1) = 1.D0/2.D0
      A(3,1) = 0.D0
      A(3,2) = 3.D0/4.D0
      A(4,1) = 2.D0/9.D0
      A(4,2) = 1.D0/3.D0
      A(4,3) = 4.D0/9.D0
C
C  The coefficients BHAT(*) refer to the formula used to advance the
C  integration, here the one of order 3.  The coefficients B(*) refer
C  to the other formula, here the one of order 2. For this pair, BHAT(*)
C  is not needed since FSAL = .TRUE.
C
      B(1) = 7.D0/24.D0
      B(2) = 1.D0/4.D0
      B(3) = 1.D0/3.D0
      B(4) = 1.D0/8.D0
C
      C(1) = 0.D0
      C(2) = 1.D0/2.D0
      C(3) = 3.D0/4.D0
      C(4) = 1.D0
C
      GO TO 120
C
C  METHD = 2
C    This pair is from "An Efficient Runge-Kutta (4,5) Pair" by P. Bogacki
C    and L.F. Shampine, Rept. 89-20, Math. Dept., Southern Methodist
C    University, Dallas, Texas, USA, 1989.  The authors are grateful to
C    P. Bogacki for his assistance in implementing the pair.  Shampine and
C    Bogacki subsequently modified the formula to enhance the reliability of
C    the pair.  The original fourth order formula is used in an estimate of
C    the local error.  If the step fails, the computation is broken off.  If
C    the step is acceptable, the first evaluation of the next step is done,
C    i.e., the pair is implemented as FSAL and the local error of the step
C    is again estimated with a fourth order formula using the additional data.
C    The step must succeed with both estimators to be accepted.  When the
C    second estimate is formed, it is used for the subsequent adjustment of
C    the step size because it is of higher quality.  The two fourth order
C    formulas are well matched to leading order, and only exceptionally do
C    the estimators disagree -- problems with discontinuous coefficients are
C    handled more reliably by using two estimators as is global error
C    estimation.
C
   40 CONTINUE
      NSTAGE = 8
      FSAL = .TRUE.
      ORDER = 4
      TANANG = 5.2D0
      STBRAD = 3.9D0
      SAFETY = 0.8D0
      INTP = .TRUE.
      REQSTG = .TRUE.
      MINTP = 6
      LINTPL = 6
      NSEC = 2
C
      PTR(1) = 0
      PTR(2) = 1
      PTR(3) = 2
      PTR(4) = 3
      PTR(5) = 4
      PTR(6) = 5
      PTR(7) = 6
      PTR(8) = 7
C
      A(2,1) = 1.D0/6.D0
      A(3,1) = 2.D0/27.D0
      A(3,2) = 4.D0/27.D0
      A(4,1) = 183.D0/1372.D0
      A(4,2) = -162.D0/343.D0
      A(4,3) = 1053.D0/1372.D0
      A(5,1) = 68.D0/297.D0
      A(5,2) = -4.D0/11.D0
      A(5,3) = 42.D0/143.D0
      A(5,4) = 1960.D0/3861.D0
      A(6,1) = 597.D0/22528.D0
      A(6,2) = 81.D0/352.D0
      A(6,3) = 63099.D0/585728.D0
      A(6,4) = 58653.D0/366080.D0
      A(6,5) = 4617.D0/20480.D0
      A(7,1) = 174197.D0/959244.D0
      A(7,2) = -30942.D0/79937.D0
      A(7,3) = 8152137.D0/19744439.D0
      A(7,4) = 666106.D0/1039181.D0
      A(7,5) = -29421.D0/29068.D0
      A(7,6) = 482048.D0/414219.D0
      A(8,1) = 587.D0/8064.D0
      A(8,2) = 0.D0
      A(8,3) = 4440339.D0/15491840.D0
      A(8,4) = 24353.D0/124800.D0
      A(8,5) = 387.D0/44800.D0
      A(8,6) = 2152.D0/5985.D0
      A(8,7) = 7267.D0/94080.D0
C
C  The coefficients B(*) refer to the formula of order 4.
C
      B(1) = 2479.D0/34992.D0
      B(2) = 0.D0
      B(3) = 123.D0/416.D0
      B(4) = 612941.D0/3411720.D0
      B(5) = 43.D0/1440.D0
      B(6) = 2272.D0/6561.D0
      B(7) = 79937.D0/1113912.D0
      B(8) = 3293.D0/556956.D0
C
C  The coefficients E(*) refer to an estimate of the local error based on
C  the first formula of order 4.  It is the difference of the fifth order
C  result, here located in A(8,*), and the fourth order result.  By
C  construction both E(2) and E(7) are zero.
C
      E(1) = -3.D0/1280.D0
      E(2) = 0.D0
      E(3) = 6561.D0/632320.D0
      E(4) = -343.D0/20800.D0
      E(5) = 243.D0/12800.D0
      E(6) = -1.D0/95.D0
      E(7) = 0.D0
C
      C(1) = 0.D0
      C(2) = 1.D0/6.D0
      C(3) = 2.D0/9.D0
      C(4) = 3.D0/7.D0
      C(5) = 2.D0/3.D0
      C(6) = 3.D0/4.D0
      C(7) = 1.D0
      C(8) = 1.D0
C
C  To do interpolation with this pair, some extra stages have to be computed.
C  The following additional A(*,*) and C(*) coefficients are for this purpose.
C  In addition there is an array R(*,*) that plays a role for interpolation
C  analogous to that of BHAT(*) for the basic step.
C
      C(9) = 1.D0/2.D0
      C(10) = 5.D0/6.D0
      C(11) = 1.D0/9.D0
C
      A(9,1) = 455.D0/6144.D0
      A(10,1) = -837888343715.D0/13176988637184.D0
      A(11,1) = 98719073263.D0/1551965184000.D0
      A(9,2) = 0.D0
      A(10,2) = 30409415.D0/52955362.D0
      A(11,2) = 1307.D0/123552.D0
      A(9,3) = 10256301.D0/35409920.D0
      A(10,3) = -48321525963.D0/759168069632.D0
      A(11,3) = 4632066559387.D0/70181753241600.D0
      A(9,4) = 2307361.D0/17971200.D0
      A(10,4) = 8530738453321.D0/197654829557760.D0
      A(11,4) = 7828594302389.D0/382182512025600.D0
      A(9,5) = -387.D0/102400.D0
      A(10,5) = 1361640523001.D0/1626788720640.D0
      A(11,5) = 40763687.D0/11070259200.D0
      A(9,6) = 73.D0/5130.D0
      A(10,6) = -13143060689.D0/38604458898.D0
      A(11,6) = 34872732407.D0/224610586200.D0
      A(9,7) = -7267.D0/215040.D0
      A(10,7) = 18700221969.D0/379584034816.D0
      A(11,7) = -2561897.D0/30105600.D0
      A(9,8) = 1.D0/32.D0
      A(10,8) = -5831595.D0/847285792.D0
      A(11,8) = 1.D0/10.D0
      A(10,9) = -5183640.D0/26477681.D0
      A(11,9) = -1.D0/10.D0
      A(11,10) = -1403317093.D0/11371610250.D0
C
      DO 60 I = 1, 11
         R(I,1) = 0.D0
   60 CONTINUE
      DO 80 I = 1, 6
         R(2,I) = 0.D0
   80 CONTINUE
      R(1,6) = -12134338393.D0/1050809760.D0
      R(1,5) = -1620741229.D0/50038560.D0
      R(1,4) = -2048058893.D0/59875200.D0
      R(1,3) = -87098480009.D0/5254048800.D0
      R(1,2) = -11513270273.D0/3502699200.D0
C
      R(3,6) = -33197340367.D0/1218433216.D0
      R(3,5) = -539868024987.D0/6092166080.D0
      R(3,4) = -39991188681.D0/374902528.D0
      R(3,3) = -69509738227.D0/1218433216.D0
      R(3,2) = -29327744613.D0/2436866432.D0
C
      R(4,6) = -284800997201.D0/19905339168.D0
      R(4,5) = -7896875450471.D0/165877826400.D0
      R(4,4) = -333945812879.D0/5671036800.D0
      R(4,3) = -16209923456237.D0/497633479200.D0
      R(4,2) = -2382590741699.D0/331755652800.D0
C
      R(5,6) = -540919.D0/741312.D0
      R(5,5) = -103626067.D0/43243200.D0
      R(5,4) = -633779.D0/211200.D0
      R(5,3) = -32406787.D0/18532800.D0
      R(5,2) = -36591193.D0/86486400.D0
C
      R(6,6) = 7157998304.D0/374350977.D0
      R(6,5) = 30405842464.D0/623918295.D0
      R(6,4) = 183022264.D0/5332635.D0
      R(6,3) = -3357024032.D0/1871754885.D0
      R(6,2) = -611586736.D0/89131185.D0
C
      R(7,6) = -138073.D0/9408.D0
      R(7,5) = -719433.D0/15680.D0
      R(7,4) = -1620541.D0/31360.D0
      R(7,3) = -385151.D0/15680.D0
      R(7,2) = -65403.D0/15680.D0
C
      R(8,6) = 1245.D0/64.D0
      R(8,5) = 3991.D0/64.D0
      R(8,4) = 4715.D0/64.D0
      R(8,3) = 2501.D0/64.D0
      R(8,2) = 149.D0/16.D0
      R(8,1) = 1.D0
C
      R(9,6) = 55.D0/3.D0
      R(9,5) = 71.D0
      R(9,4) = 103.D0
      R(9,3) = 199.D0/3.D0
      R(9,2) = 16.0D0
C
      R(10,6) = -1774004627.D0/75810735.D0
      R(10,5) = -1774004627.D0/25270245.D0
      R(10,4) = -26477681.D0/359975.D0
      R(10,3) = -11411880511.D0/379053675.D0
      R(10,2) = -423642896.D0/126351225.D0
C
      R(11,6) = 35.D0
      R(11,5) = 105.D0
      R(11,4) = 117.D0
      R(11,3) = 59.D0
      R(11,2) = 12.D0
C
      GO TO 120
C
C  METHD = 3
C    This pair is from "High Order Embedded Runge-Kutta Formulae" by P.J.
C    Prince and J.R. Dormand, J. Comp. Appl. Math.,7, pp. 67-75, 1981.  The
C    authors are grateful to P. Prince and J. Dormand for their assistance in
C    implementing the pair.
C
  100 CONTINUE
      NSTAGE = 13
      FSAL = .FALSE.
      ORDER = 7
      TANANG = 11.0D0
      STBRAD = 5.2D0
      SAFETY = 0.8D0
      INTP = .FALSE.
      REQSTG = .FALSE.
      MINTP = 0
      LINTPL = 0
      NSEC = 2
C
      PTR(1) = 0
      PTR(2) = 1
      PTR(3) = 2
      PTR(4) = 1
      PTR(5) = 3
      PTR(6) = 2
      PTR(7) = 4
      PTR(8) = 5
      PTR(9) = 6
      PTR(10) = 7
      PTR(11) = 8
      PTR(12) = 9
      PTR(13) = 1
C
      A(2,1) = 5.55555555555555555555555555556D-2
      A(3,1) = 2.08333333333333333333333333333D-2
      A(3,2) = 6.25D-2
      A(4,1) = 3.125D-2
      A(4,2) = 0.D0
      A(4,3) = 9.375D-2
      A(5,1) = 3.125D-1
      A(5,2) = 0.D0
      A(5,3) = -1.171875D0
      A(5,4) = 1.171875D0
      A(6,1) = 3.75D-2
      A(6,2) = 0.D0
      A(6,3) = 0.D0
      A(6,4) = 1.875D-1
      A(6,5) = 1.5D-1
      A(7,1) = 4.79101371111111111111111111111D-2
      A(7,2) = 0.D0
      A(7,3) = 0.0D0
      A(7,4) = 1.12248712777777777777777777778D-1
      A(7,5) = -2.55056737777777777777777777778D-2
      A(7,6) = 1.28468238888888888888888888889D-2
      A(8,1) = 1.6917989787292281181431107136D-2
      A(8,2) = 0.D0
      A(8,3) = 0.D0
      A(8,4) = 3.87848278486043169526545744159D-1
      A(8,5) = 3.59773698515003278967008896348D-2
      A(8,6) = 1.96970214215666060156715256072D-1
      A(8,7) = -1.72713852340501838761392997002D-1
      A(9,1) = 6.90957533591923006485645489846D-2
      A(9,2) = 0.D0
      A(9,3) = 0.D0
      A(9,4) = -6.34247976728854151882807874972D-1
      A(9,5) = -1.61197575224604080366876923982D-1
      A(9,6) = 1.38650309458825255419866950133D-1
      A(9,7) = 9.4092861403575626972423968413D-1
      A(9,8) = 2.11636326481943981855372117132D-1
      A(10,1) = 1.83556996839045385489806023537D-1
      A(10,2) = 0.D0
      A(10,3) = 0.D0
      A(10,4) = -2.46876808431559245274431575997D0
      A(10,5) = -2.91286887816300456388002572804D-1
      A(10,6) = -2.6473020233117375688439799466D-2
      A(10,7) = 2.84783876419280044916451825422D0
      A(10,8) = 2.81387331469849792539403641827D-1
      A(10,9) = 1.23744899863314657627030212664D-1
      A(11,1) = -1.21542481739588805916051052503D0
      A(11,2) = 0.D0
      A(11,3) = 0.D0
      A(11,4) = 1.66726086659457724322804132886D1
      A(11,5) = 9.15741828416817960595718650451D-1
      A(11,6) = -6.05660580435747094755450554309D0
      A(11,7) = -1.60035735941561781118417064101D1
      A(11,8) = 1.4849303086297662557545391898D1
      A(11,9) = -1.33715757352898493182930413962D1
      A(11,10) = 5.13418264817963793317325361166D0
      A(12,1) = 2.58860916438264283815730932232D-1
      A(12,2) = 0.D0
      A(12,3) = 0.D0
      A(12,4) = -4.77448578548920511231011750971D0
      A(12,5) = -4.3509301377703250944070041181D-1
      A(12,6) = -3.04948333207224150956051286631D0
      A(12,7) = 5.57792003993609911742367663447D0
      A(12,8) = 6.15583158986104009733868912669D0
      A(12,9) = -5.06210458673693837007740643391D0
      A(12,10) = 2.19392617318067906127491429047D0
      A(12,11) = 1.34627998659334941535726237887D-1
      A(13,1) = 8.22427599626507477963168204773D-1
      A(13,2) = 0.D0
      A(13,3) = 0.D0
      A(13,4) = -1.16586732572776642839765530355D1
      A(13,5) = -7.57622116690936195881116154088D-1
      A(13,6) = 7.13973588159581527978269282765D-1
      A(13,7) = 1.20757749868900567395661704486D1
      A(13,8) = -2.12765911392040265639082085897D0
      A(13,9) = 1.99016620704895541832807169835D0
      A(13,10) = -2.34286471544040292660294691857D-1
      A(13,11) = 1.7589857770794226507310510589D-1
      A(13,12) = 0.D0
C
C  The coefficients BHAT(*) refer to the formula used to advance the
C  integration, here the one of order 8.  The coefficients B(*) refer
C  to the other formula, here the one of order 7.
C
      BHAT(1) = 4.17474911415302462220859284685D-2
      BHAT(2) = 0.D0
      BHAT(3) = 0.D0
      BHAT(4) = 0.D0
      BHAT(5) = 0.D0
      BHAT(6) = -5.54523286112393089615218946547D-2
      BHAT(7) = 2.39312807201180097046747354249D-1
      BHAT(8) = 7.0351066940344302305804641089D-1
      BHAT(9) = -7.59759613814460929884487677085D-1
      BHAT(10) = 6.60563030922286341461378594838D-1
      BHAT(11) = 1.58187482510123335529614838601D-1
      BHAT(12) = -2.38109538752862804471863555306D-1
      BHAT(13) = 2.5D-1
C
      B(1) = 2.9553213676353496981964883112D-2
      B(2) = 0.D0
      B(3) = 0.D0
      B(4) = 0.D0
      B(5) = 0.D0
      B(6) = -8.28606276487797039766805612689D-1
      B(7) = 3.11240900051118327929913751627D-1
      B(8) = 2.46734519059988698196468570407D0
      B(9) = -2.54694165184190873912738007542D0
      B(10) = 1.44354858367677524030187495069D0
      B(11) = 7.94155958811272872713019541622D-2
      B(12) = 4.44444444444444444444444444445D-2
      B(13) = 0.D0
C
      C(1) = 0.D0
      C(2) = 5.55555555555555555555555555556D-2
      C(3) = 8.33333333333333333333333333334D-2
      C(4) = 1.25D-1
      C(5) = 3.125D-1
      C(6) = 3.75D-1
      C(7) = 1.475D-1
      C(8) = 4.65D-1
      C(9) = 5.64865451382259575398358501426D-1
      C(10) = 6.5D-1
      C(11) = 9.24656277640504446745013574318D-1
      C(12) = 1.D0
      C(13) = C(12)
C
      GO TO 120
C
C  The definitions of all pairs come here for the calculation of
C  LSTSTG, RS1, RS2, RS3, RS4, COST, MAXTRY, EXPON, TOOSML, and VECSTG.
C
  120 CONTINUE
      LSTSTG = PTR(NSTAGE)
      IF (FSAL) THEN
         COST = DBLE(NSTAGE-1)
      ELSE
         COST = DBLE(NSTAGE)
      END IF
C
C  MAXTRY - limit on the number of iterations of a computation made in
C  diagnosing stiffness.  There are at most Q = 3 function calls per
C  iteration. MAXTRY is determined so that  Q*MAXTRY <= 5% of the cost of
C  50 steps and 1 <= MAXTRY <= 8. This limits the number of calls to FCN
C  in each diagnosis of stiffness to 24 calls.
C
      MAXTRY = MIN(8,MAX(1,INT(FIVEPC*COST*FIFTY)))
C
      EXPON = ONE/(ORDER+ONE)
C
C     In calculating CDIFF it is assumed that there will be a non-zero
C     difference |C(I) - C(J)| less than one. If C(I) = C(J) for any I not
C     equal to J, they should be made precisely equal by assignment.
C
      CDIFF = ONE
      DO 160 I = 1, NSTAGE - 1
         DO 140 J = I + 1, NSTAGE
            DIFF = ABS(C(I)-C(J))
            IF (DIFF.NE.ZERO) CDIFF = MIN(CDIFF,DIFF)
  140    CONTINUE
  160 CONTINUE
      TOOSML = RNDOFF/CDIFF
C
C  Determine the number of columns needed in STAGES(1:NEQ,*) (this will be
C  at most NSTAGE-1 since the first stage is held in a separate array).
C  The PTR array contains the column positions of the stages.
C
      VECSTG = 0
      DO 180 I = 2, NSTAGE
         VECSTG = MAX(PTR(I),VECSTG)
  180 CONTINUE
C
      RS = TWO
      RS1 = ONE/RS
      RS2 = RS**2
      RS3 = RS**3
      RS4 = ONE/RS3
C
      RETURN
      END
