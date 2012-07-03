      SUBROUTINE STIFFA(F,X,Y,HNOW,HAVG,XEND,MAXFCN,WT,FXY,V0,UNSURE,
     &                  STIF,V1,V2,V3,VTEMP)
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
C************************************************
C
C  External:     F
C  Input:        X, Y(*), HNOW, HAVG, XEND, MAXFCN, WT(*), FXY(*)
C  Input/Output  V0(*)
C  Output:       UNSURE, STIF
C  Workspace:    V1(*), V2(*), V3(*), VTEMP(*)
C
C  Common:       Initializes:    none
C                Reads:          /RKCOM1/ TND, NEQN
C                                /RKCOM5/ COST, STBRAD, TANANG
C                                /RKCOM7/ SQRRMC, CUBRMC
C                Alters:         none
C
C  STIFFA diagnoses stiffness for an explicit Runge-Kutta code.  When it
C  is called, either many step failures have been observed, or a lot of
C  work has been done.
C
C  The NEQ equations of the problem are defined by the subroutine F(X,Y,YP).
C  When STIFFA is called, the integration has reached X where the approximate
C  solution is Y(*).  The vector FXY(*) is defined by a call of F(X,Y,FXY).
C  It is an input argument because it is usually available from the integrator.
C
C  The last successful step was of size HNOW, and an average step size is
C  HAVG.  A weighted norm is used to measure the local error with the error
C  in solution component L divided by the positive weight WT(L) provided in 
C  the vector WT(*).
C
C  Explicit Runge - Kutta codes estimate the local error of Y(*) by
C  forming the difference of two approximate solutions.  This difference
C  must be provided in the vector V0(*).  When this difference is too
C  small to be significant, STIFFA will replace it with a "random" vector.
C
C  STIF is set .TRUE. when the average step size appears to be restricted
C  on grounds of stability.  In certain cases the variable UNSURE is set 
C  .TRUE.; the value of STIF is then not defined.
C
C  The stability region of the explicit Runge-Kutta formula is described
C  by quantities TANANG and STBRAD that are communicated by the setup routine
C  via COMMON.  Stability regions often change sharply near the imaginary
C  axis so that it is difficult to classify the stiffness of a problem with
C  eigenvalues of a local Jacobian that are "near" the imaginary axis.  For
C  this reason,  we consider only points Z in the upper left half complex
C  plane for which TAN( IMAG(Z)/( - RE(Z))) <= TANANG. Eigenvalues outside
C  this region are one reason for the code being UNSURE.  The stability
C  region is approximated by the intersection of a disk with this sector.
C  The radius of this disk is called STBRAD.
C
C  Working storage must be provided via the four vectors V1(*),V2(*),
C  V3(*),VTEMP(*).  These vectors must be of length at least NEQ.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HAVG, HNOW, X, XEND
      INTEGER           MAXFCN
      LOGICAL           STIF, UNSURE
C     .. Array Arguments ..
      DOUBLE PRECISION  FXY(*), V0(*), V1(*), V2(*), V3(*), VTEMP(*),
     &                  WT(*), Y(*)
C     .. Subroutine Arguments ..
      EXTERNAL          F
C     .. Common Block for Problem Definition ..
      DOUBLE PRECISION  TSTRT, TND, DIR, HSTRT, TOLR
      INTEGER           NEQN
      COMMON /RKCOM1/   TSTRT, TND, DIR, HSTRT, TOLR, NEQN
      SAVE   /RKCOM1/
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
      DOUBLE PRECISION  LARGE
      PARAMETER         (LARGE=1.0D+10)
      DOUBLE PRECISION  ZERO, P001, P9, ONE, TWO, FIVE, FIFTH
      PARAMETER         (ZERO=0.0D+0,P001=0.001D+0,P9=0.9D+0,ONE=1.0D+0,
     &                  TWO=2.0D+0,FIVE=5.0D+0,FIFTH=0.2D+0)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA1, ALPHA2, BETA1, BETA2, D1, D2, DET1,
     &                  DET2, DIST, RES2, RHO, RHO2, ROLD, SCALE, V0NRM,
     &                  V0V0, V0V1, V0V2, V1V1, V1V2, V1V3, V2V2, V2V3,
     &                  V3NRM, V3V3, XTRFCN, YNRM
      INTEGER           L, NTRY
      LOGICAL           ROOTRE
C     .. Local Arrays ..
      DOUBLE PRECISION  R1(2), R2(2), ROOT1(2), ROOT2(2)
C     .. External Functions ..
      DOUBLE PRECISION  DOTPRD
      EXTERNAL          DOTPRD
C     .. External Subroutines ..
      EXTERNAL          STIFFB, STIFFC, STIFFD
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN, SQRT
C     .. Executable Statements ..
C
C  If the current step size differs substantially from the average,
C  the problem is not stiff.
C
      IF (ABS(HNOW/HAVG).GT.FIVE .OR. ABS(HNOW/HAVG).LT.FIFTH) THEN
         STIF = .FALSE.
         UNSURE = .FALSE.
         RETURN
      ELSE
         UNSURE = .TRUE.
      END IF
C
C  The average step size is used to predict the cost in function evaluations
C  of finishing the integration to XEND.  If this cost is no more than MAXFCN,
C  the problem is declared not stiff: If the step size is being restricted on
C  grounds of stability, it will stay close to HAVG.  The prediction will
C  then be good, but the cost is too low to consider the problem stiff.  If
C  the step size is not close to HAVG, the problem is not stiff.  Either way
C  there is no point to testing for a step size restriction due to stability.
C
      XTRFCN = COST*ABS((XEND-X)/HAVG)
      IF (XTRFCN.LE.MAXFCN) THEN
         STIF = .FALSE.
         UNSURE = .FALSE.
         RETURN
      ELSE
         UNSURE = .TRUE.
      END IF
C
C  There have been many step failures or a lot of work has been done.  Now 
C  we must determine if this is due to the stability characteristics of the
C  formula.  This is done by calculating the dominant eigenvalues of the
C  local Jacobian and then testing whether HAVG corresponds to being on the
C  boundary of the stability region.
C
C  The size of Y(*) provides scale information needed to approximate
C  the Jacobian by differences.
C
      YNRM = SQRT(DOTPRD(Y,Y,WT,NEQN))
      SCALE = YNRM*SQRRMC
      IF (SCALE.EQ.ZERO) THEN
C
C  Degenerate case.  Y(*) is (almost) the zero vector so the scale is not 
C  defined.  The input vector V0(*) is the difference between Y(*) and a 
C  lower order approximation to the solution that is within the error 
C  tolerance.  When Y(*) vanishes, V0(*) is itself an acceptable approximate
C  solution, so we take SCALE from it, if this is possible.
C
         YNRM = SQRT(DOTPRD(V0,V0,WT,NEQN))
         SCALE = YNRM*SQRRMC
         IF (SCALE.EQ.ZERO) THEN
            UNSURE = .TRUE.
            RETURN
         END IF
      END IF
C
      V0V0 = DOTPRD(V0,V0,WT,NEQN)
      IF (V0V0.EQ.ZERO) THEN
C
C  Degenerate case.  V0(*) is (almost) the zero vector so cannot
C  be used to define a direction for an increment to Y(*).  Try a
C  "random" direction.
C
         DO 20 L = 1, NEQN
            V0(L) = ONE
   20    CONTINUE
         V0V0 = DOTPRD(V0,V0,WT,NEQN)
      END IF
      V0NRM = SQRT(V0V0)
      DO 40 L = 1, NEQN
         V0(L) = V0(L)/V0NRM
   40 CONTINUE
      V0V0 = ONE
C
C  Use a nonlinear power method to estimate the two dominant eigenvalues.
C  V0(*) is often very rich in the two associated eigenvectors.  For this 
C  reason the computation is organized with the expectation that a minimal 
C  number of iterations will suffice.  Indeed, it is necessary to recognize 
C  a kind of degeneracy when there is a dominant real eigenvalue.  The
C  subroutine STIFFB does this.  In the first try, NTRY = 1, a Rayleigh 
C  quotient for such an eigenvalue is initialized as ROLD.  After each 
C  iteration, REROOT computes a new Rayleigh quotient and tests whether the
C  two approximations agree to one tenth of one per cent and the eigenvalue,
C  eigenvector pair satisfy a stringent test on the residual.  ROOTRE = .TRUE.
C  signals that a single dominant real root has been found.
C
      NTRY = 1
   60 CONTINUE
C
      CALL STIFFD(V0,HAVG,X,Y,F,FXY,WT,SCALE,V0V0,V1,V1V1,VTEMP)
C
C  The quantity SQRT(V1V1/V0V0) is a lower bound for the product of HAVG
C  and a Lipschitz constant.  If it should be LARGE, stiffness is not
C  restricting the step size to the stability region.  The principle is
C  clear enough, but the real reason for this test is to recognize an
C  extremely inaccurate computation of V1V1 due to finite precision
C  arithmetic in certain degenerate circumstances.
C
      IF (SQRT(V1V1).GT.LARGE*SQRT(V0V0)) THEN
         UNSURE = .TRUE.
         RETURN
      END IF
C
      V0V1 = DOTPRD(V0,V1,WT,NEQN)
      IF (NTRY.EQ.1) THEN
         ROLD = V0V1/V0V0
C
C  This is the first Rayleigh quotient approximating the product of HAVG
C  and a dominant real eigenvalue.  If it should be very small, the
C  problem is not stiff.  It is important to test for this possibility so
C  as to prevent underflow and degeneracies in the subsequent iteration.
C
         IF (ABS(ROLD).LT.CUBRMC) THEN
            UNSURE = .FALSE.
            STIF = .FALSE.
            RETURN
         END IF
      ELSE
         CALL STIFFB(V1V1,V0V1,V0V0,ROLD,RHO,ROOT1,ROOT2,ROOTRE)
         IF (ROOTRE) GO TO 100
      END IF
      CALL STIFFD(V1,HAVG,X,Y,F,FXY,WT,SCALE,V1V1,V2,V2V2,VTEMP)
      V0V2 = DOTPRD(V0,V2,WT,NEQN)
      V1V2 = DOTPRD(V1,V2,WT,NEQN)
      CALL STIFFB(V2V2,V1V2,V1V1,ROLD,RHO,ROOT1,ROOT2,ROOTRE)
      IF (ROOTRE) GO TO 100
C
C  Fit a quadratic in the eigenvalue to the three successive iterates
C  V0(*),V1(*),V2(*) of the power method to get a first approximation to
C  a pair of eigenvalues.  A test made earlier in STIFFB implies that
C  the quantity DET1 here will not be too small.
C
      DET1 = V0V0*V1V1 - V0V1**2
      ALPHA1 = (-V0V0*V1V2+V0V1*V0V2)/DET1
      BETA1 = (V0V1*V1V2-V1V1*V0V2)/DET1
C
C  Iterate again to get V3, test again for degeneracy, and then fit a
C  quadratic to V1(*),V2(*),V3(*) to get a second approximation to a pair
C  of eigenvalues.
C
      CALL STIFFD(V2,HAVG,X,Y,F,FXY,WT,SCALE,V2V2,V3,V3V3,VTEMP)
      V1V3 = DOTPRD(V1,V3,WT,NEQN)
      V2V3 = DOTPRD(V2,V3,WT,NEQN)
      CALL STIFFB(V3V3,V2V3,V2V2,ROLD,RHO,ROOT1,ROOT2,ROOTRE)
      IF (ROOTRE) GO TO 100
      DET2 = V1V1*V2V2 - V1V2**2
      ALPHA2 = (-V1V1*V2V3+V1V2*V1V3)/DET2
      BETA2 = (V1V2*V2V3-V2V2*V1V3)/DET2
C
C  First test the residual of the quadratic fit to see if we might
C  have determined a pair of eigenvalues.
C
      RES2 = ABS(V3V3+V2V2*ALPHA2**2+V1V1*BETA2**2+TWO*V2V3*ALPHA2+
     &       TWO*V1V3*BETA2+TWO*V1V2*ALPHA2*BETA2)
      IF (RES2.LE.V3V3*P001**2) THEN
C
C  Calculate the two approximate pairs of eigenvalues.
C
         CALL STIFFC(ALPHA1,BETA1,R1,R2)
         CALL STIFFC(ALPHA2,BETA2,ROOT1,ROOT2)
C
C  The test for convergence is done on the larger root of the second
C  approximation.  It is complicated by the fact that one pair of roots 
C  might be real and the other complex.  First calculate the spectral 
C  radius RHO of HAVG*J as the magnitude of ROOT1.  Then see if one of 
C  the roots R1,R2 is within one per cent of ROOT1.  A subdominant root 
C  may be very poorly approximated if its magnitude is much smaller than 
C  RHO -- this does not matter in our use of these eigenvalues.
C
         RHO = SQRT(ROOT1(1)**2+ROOT1(2)**2)
         D1 = (ROOT1(1)-R1(1))**2 + (ROOT1(2)-R1(2))**2
         D2 = (ROOT1(1)-R2(1))**2 + (ROOT1(2)-R2(2))**2
         DIST = SQRT(MIN(D1,D2))
         IF (DIST.LE.P001*RHO) GO TO 100
      END IF
C
C  Do not have convergence yet.  Because the iterations are cheap, and
C  because the convergence criterion is stringent, we are willing to try
C  a few iterations.
C
      IF (NTRY.LT.MAXTRY) THEN
         NTRY = NTRY + 1
         V3NRM = SQRT(V3V3)
         DO 80 L = 1, NEQN
            V0(L) = V3(L)/V3NRM
   80    CONTINUE
         V0V0 = ONE
         GO TO 60
      ELSE
         UNSURE = .TRUE.
         RETURN
      END IF
C
C                        **************
C
C  We now have the dominant eigenvalues.  Decide if the average step
C  size is being restricted on grounds of stability.  Check the real
C  parts of the eigenvalues.  First see if the dominant eigenvalue is
C  in the left half plane -- there won't be a stability restriction
C  unless it is. If there is another eigenvalue of comparable magnitude
C  with a positive real part, the problem is not stiff. If the dominant
C  eigenvalue is too close to the imaginary axis, we cannot diagnose
C  stiffness.
C
  100 CONTINUE
      IF (ROOT1(1).GT.ZERO) THEN
         STIF = .FALSE.
         UNSURE = .FALSE.
         RETURN
      END IF
      RHO2 = SQRT(ROOT2(1)**2+ROOT2(2)**2)
      IF (RHO2.GE.P9*RHO .AND. ROOT2(1).GT.ZERO) THEN
         STIF = .FALSE.
         UNSURE = .FALSE.
         RETURN
      END IF
      IF (ABS(ROOT1(2)).GT.ABS(ROOT1(1))*TANANG) THEN
         UNSURE = .TRUE.
         RETURN
      END IF
C
C  If the average step size corresponds to being well within the
C  stability region, the step size is not being restricted because
C  of stability.
C
      STIF = RHO .GE. P9*STBRAD
      UNSURE = .FALSE.
      RETURN
      END
