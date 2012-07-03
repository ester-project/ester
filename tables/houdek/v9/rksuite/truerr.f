      SUBROUTINE TRUERR(F,NEQ,Y,TOL,WEIGHT,ZY,ZYP,ZERROR,ZYNEW,ZERRES,
     &                  ZSTAGE,IER)
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
C************************************************
C
C  Purpose:      Compute a running RMS measure of the true (global) error
C                for a general Runge-Kutta pair.
C
C
C  Input:        NEQ, Y(*), TOL, WEIGHT(*),
C  Input/output: ZY(*), ZYP(*), ZERROR(*)
C  Workspace:    ZYNEW(*), ZERRES(*), ZSTAGE(NEQ,*)
C  Output:       IER
C  External:     F
C
C  Common:       Initializes:    none
C                Reads:          /RKCOM2/ T, HOLD
C                                /RKCOM5/ TOOSML, ORDER, NSEC
C                                /RKCOM7/ TINY
C                Alters:         /RKCOM6/ MAXERR, LOCMAX, GNFCN
C
C  Comments:
C  =========
C  A secondary integration is performed using a fraction of the step size 
C  of the primary integration. ZY(*) and ZYP(*) are the approximate solution
C  and first derivative of this secondary integration. ZERRES(*) contains the 
C  error estimates for the secondary integration. ZYNEW(*) and ZSTAGE(*,*) are
C  workspace for taking a step. The error assessment is computed using the
C  difference of the primary and secondary solutions at the primary
C  integration points as an estimate of the true error there.  The weights 
C  used are those of the error test of the primary integration. This error 
C  assessment is maintained in the vector ZERROR(*).  MAXERR and LOCMAX 
C  contain the maximum contribution to the assessment and its location,
C  respectively.  The number of calls to F is counted by GNFCN.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           IER, NEQ
C     .. Array Arguments ..
      DOUBLE PRECISION  WEIGHT(*), Y(*), ZERRES(*), ZERROR(*),
     &                  ZSTAGE(NEQ,*), ZY(*), ZYNEW(*), ZYP(*)
C     .. Subroutine Arguments ..
      EXTERNAL          F
C     .. Common Block to hold Problem Status ..
      DOUBLE PRECISION  T, H, TOLD, HOLD
      INTEGER           NFCN, SVNFCN, OKSTP, FLSTP
      LOGICAL           FIRST, LAST
      COMMON /RKCOM2/   T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP, FLSTP,
     &                  FIRST, LAST
      SAVE   /RKCOM2/
C     .. Common Block to hold Formula Characterisitcs ..
      DOUBLE PRECISION  TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
     &                  RS, RS1, RS2, RS3, RS4
      INTEGER           ORDER, LSTSTG, MAXTRY, NSEC
      LOGICAL           FSAL
      COMMON /RKCOM5/   TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
     &                  RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG, MAXTRY,
     &                  NSEC, FSAL
      SAVE   /RKCOM5/
C     .. Common Block for Global Error Assessment ..
      DOUBLE PRECISION  MAXERR, LOCMAX
      INTEGER           GNFCN, PRZSTG, PRZY, PRZYP, PRZERS, PRZERR,
     &                  PRZYNU
      LOGICAL           ERASON, ERASFL
      COMMON /RKCOM6/   MAXERR, LOCMAX, GNFCN, PRZSTG, PRZY, PRZYP,
     &                  PRZERS, PRZERR, PRZYNU, ERASON, ERASFL
      SAVE   /RKCOM6/
C     .. Common Block for Environment Parameters ..
      DOUBLE PRECISION  MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY
      INTEGER           OUTCH
      COMMON /RKCOM7/   MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY,
     &                  OUTCH
      SAVE   /RKCOM7/
C     .. Parameters ..
      DOUBLE PRECISION  PT1, TEN, DUMMY
      PARAMETER         (PT1=0.1D0,TEN=10.0D0,DUMMY=1.0D0)
C     .. Local Scalars ..
      DOUBLE PRECISION  DIFF, ERRMAX, HMIN, HSEC, MXERLC, TSEC, ZLERR,
     &                  ZTEST1, ZTEST2
      INTEGER           ISTEP, L, LEVEL
      LOGICAL           LDUMMY, MAIN
C     .. Local Arrays ..
      DOUBLE PRECISION  DUMARR(1)
C     .. External Subroutines ..
      EXTERNAL          STEP
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX
C     .. Executable Statements ..
      TSEC = T - HOLD
      HSEC = HOLD/DBLE(NSEC)
      HMIN = MAX(TINY,TOOSML*MAX(ABS(TSEC),ABS(T)))
      IF (ABS(HSEC).LT.HMIN) THEN
         IER = 6
         GO TO 120
      END IF
      ZTEST1 = TOL/DBLE(NSEC)
      ZTEST2 = TOL/TEN
      LEVEL = 0
C
C  The subroutine STEP is used to take a step.  In its use in the primary
C  integration provision is made for getting on scale in the first step.
C  In this situation only the subroutine might reduce the step size.  By
C  setting MAIN = .FALSE., the subroutine will take a step of the size input.
C  In this use of the subroutine, all items of the call list appearing after
C  MAIN are dummy variables.
C
C  Perform secondary integration.
      MAIN = .FALSE.
      LDUMMY = .FALSE.
      DO 60 ISTEP = 1, NSEC
C
C  Take a step.
         CALL STEP(F,NEQ,TSEC,ZY,ZYP,ZSTAGE,ZTEST1,HSEC,WEIGHT,ZYNEW,
     &             ZERRES,ZLERR,MAIN,DUMMY,DUMARR,LDUMMY)
C
C  The primary integration is using a step size of HUSED and the secondary
C  integration is using the smaller step size HSEC = HUSED/NSEC.  If steps
C  of this size were taken from the same starting point and the asymptotic
C  behavior were evident, the smaller step size would result in a local error
C  that is considerably smaller, namely by a factor of 1/(NSEC**(ORDER+1)).
C  If the two approximate solutions are close and TOLR is neither too large nor
C  too small, this should be approximately true.  The step size is chosen in
C  the primary integration so that the local error ERR is no larger than TOLR.
C  The local error, ZLERR, of the secondary integration is compared to TOLR in
C  an attempt to diagnose a secondary integration that is not rather more
C  accurate than the primary integration.
C
         IF (ZLERR.GE.ZTEST1) THEN
            LEVEL = 2
         ELSE IF (ZLERR.GT.ZTEST2) THEN
            LEVEL = LEVEL + 1
         END IF
         IF (LEVEL.GE.2) THEN
            IER = 6
            GO TO 120
         END IF
C
C  Advance TSEC and the dependent variables ZY(*) and ZYP(*).
         TSEC = T - DBLE(NSEC-ISTEP)*HSEC
         DO 20 L = 1, NEQ
            ZY(L) = ZYNEW(L)
   20    CONTINUE
C
         IF (FSAL) THEN
C
C  When FSAL = .TRUE., the derivative ZYP(*) is the last stage of the step.
            DO 40 L = 1, NEQ
               ZYP(L) = ZSTAGE(L,LSTSTG)
   40       CONTINUE
         ELSE
C
C  Call F to evaluate ZYP(*).
            CALL F(TSEC,ZY,ZYP)
            GNFCN = GNFCN + 1
         END IF
C
   60 CONTINUE
C
C  Update the maximum error seen, MAXERR, and its location, LOCMAX.
C  Use local variables ERRMAX and MXERLC.
C
      ERRMAX = MAXERR
      MXERLC = LOCMAX
      DO 80 L = 1, NEQ
         DIFF = ABS(ZY(L)-Y(L))/WEIGHT(L)
         IF (DIFF.GT.ERRMAX) THEN
            ERRMAX = DIFF
            MXERLC = T
         END IF
   80 CONTINUE
C
C  If the global error is greater than 0.1D0, the solutions have diverged so
C  far that comparing them may not provide a reliable estimate of the global
C  error. The test is made before ZERROR(*) and MAXERR, LCMXER are updated so
C  that on a failure, they refer to the last reliable results.
C
      IF (ERRMAX.GT.PT1) THEN
         IER = 6
         GO TO 120
      ELSE
         MAXERR = ERRMAX
         LOCMAX = MXERLC
         DO 100 L = 1, NEQ
            DIFF = ABS(ZY(L)-Y(L))/WEIGHT(L)
            ZERROR(L) = ZERROR(L) + DIFF**2
  100    CONTINUE
         IER = 1
      END IF
C
C  Exit point for TRUERR
  120 CONTINUE
C
      RETURN
      END
