      SUBROUTINE STEP(F,NEQ,TNOW,Y,YP,STAGES,TOL,HTRY,WEIGHT,YNEW,
     &                ERREST,ERR,MAIN,HMIN,THRES,PHASE2)
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
C************************************************
C
C  Purpose:      To compute a step of an explicit Runge-Kutta
C                method and estimate the local error of the step.
C
C  Input:        NEQ, TNOW, Y(*), YP(*), TOL, MAIN, HMIN, THRES(*)
C  Input/output: HTRY, PHASE2, LAST, WEIGHT(*)
C  Output:       STAGES(NEQ,*), YNEW(*), ERREST(*), ERR
C
C  Common:       Initializes:    none
C                Reads:          /RKCOM1/ TND
C                                /RKCOM2/ LAST
C                                /RKCOM4/ A, B, C, BHAT, PTR, NSTAGE, METHD
C                                /RKCOM5/ FSAL
C                Alters:         /RKCOM2/ NFCN, LAST
C                                /RKCOM6/ GNFCN
C
C  Comments:
C  =========
C  From an approximate solution Y(*) at TNOW and first derivative there,
C  YP(*) = F(TNOW,Y,YP), a step is taken to get an approximation YNEW(*)
C  at TNOW + HTRY. The Runge-Kutta method and how it is used are defined
C  by A, B, C, BHAT, PTR, NSTAGE, METHD and FSAL. Intermediate stages
C  of the method are stored in the array STAGES(NEQ,*). The error in
C  each solution component is estimated and returned in ERREST(*). A
C  weighted maximum norm of the local error, ERR, is formed. For some
C  methods an intermediate error estimate can be computed before completion
C  of the step (see routine STEPB); if the estimate is greater than the
C  specified tolerance TOL, the computation of the step is terminated.
C
C  When global error estimation is desired, two integrations are done.
C  The usual integration is referred to as the "primary", or "main",
C  integration (MAIN=.TRUE.).  For global error estimation another,
C  "secondary" integration (MAIN=.FALSE.) is carried out with a smaller
C  step size.  The weight vector WEIGHT(*) used in computing ERR is
C  determined by the main integration.  Thus this argument is output when
C  MAIN = .TRUE. and input when MAIN = .FALSE..
C
C  When taking the first step in an integration, the logical variable
C  PHASE2 may be input as .TRUE. and if the first step is the whole of
C  the range of integration, then LAST will be .TRUE.. When PHASE2=.TRUE.,
C  the first three stages are monitored to help assure that the step
C  size H is small enough for the integration to be stable and for the
C  estimate of the error of the step to be credible. Calls are made to
C  the subroutine STEPA for this purpose. If necessary, H will be
C  reduced in STEPA (and LAST altered accordingly) and the step retried
C  in STEP until an acceptable value is found.
C
C  In the primary integration the number of calls to F is counted by
C  NFCN, and in the secondary integration, by GNFCN.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ERR, HMIN, HTRY, TNOW, TOL
      INTEGER           NEQ
      LOGICAL           MAIN, PHASE2
C     .. Array Arguments ..
      DOUBLE PRECISION  ERREST(*), STAGES(NEQ,*), THRES(*), WEIGHT(*),
     &                  Y(*), YNEW(*), YP(*)
C     .. Subroutine Arguments ..
      EXTERNAL          F
C     .. Common Block for Problem Definition ..
      DOUBLE PRECISION  TSTRT, TND, DIR, HSTRT, TOLR
      INTEGER           NEQN
      COMMON /RKCOM1/   TSTRT, TND, DIR, HSTRT, TOLR, NEQN
      SAVE   /RKCOM1/
C     .. Common Block to hold Problem Status ..
      DOUBLE PRECISION  T, H, TOLD, HOLD
      INTEGER           NFCN, SVNFCN, OKSTP, FLSTP
      LOGICAL           FIRST, LAST
      COMMON /RKCOM2/   T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP, FLSTP,
     &                  FIRST, LAST
      SAVE   /RKCOM2/
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
C     .. Common Block for Global Error Assessment ..
      DOUBLE PRECISION  MAXERR, LOCMAX
      INTEGER           GNFCN, PRZSTG, PRZY, PRZYP, PRZERS, PRZERR,
     &                  PRZYNU
      LOGICAL           ERASON, ERASFL
      COMMON /RKCOM6/   MAXERR, LOCMAX, GNFCN, PRZSTG, PRZY, PRZYP,
     &                  PRZERS, PRZERR, PRZYNU, ERASON, ERASFL
      SAVE   /RKCOM6/
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
C     .. Local Scalars ..
      DOUBLE PRECISION  AVGY, TSTG
      INTEGER           I, J, L
      LOGICAL           CUTBAK
C     .. External Subroutines ..
      EXTERNAL          STEPA, STEPB
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SIGN
C     .. Executable Statements ..
C
C  Many of the following loops over L = 1, NEQ have constant array values
C  inside. The code is written with clarity in mind.  Any optimizing
C  compiler will identify these occurrences and take appropriate action.
C  A check for zero multipliers has been included so as to prevent
C  needless computation resulting from the storing of zero coefficients
C  in the arrays for the sake of clarity.  The array ERREST(*) is used
C  for working storage in this computation.
C
   20 CONTINUE
      IF (MAIN) THEN
         IF (PHASE2) THEN
C
C  Initialize weights for measuring the local error.
            DO 40 L = 1, NEQ
               WEIGHT(L) = MAX(THRES(L),ABS(Y(L)))
   40       CONTINUE
         END IF
      END IF
C
      DO 140 I = 2, NSTAGE
         DO 100 J = 1, I - 1
            IF (J.EQ.1) THEN
               DO 60 L = 1, NEQ
                  ERREST(L) = A(I,1)*YP(L)
   60          CONTINUE
            ELSE
               IF (A(I,J).NE.ZERO) THEN
                  DO 80 L = 1, NEQ
                     ERREST(L) = ERREST(L) + A(I,J)*STAGES(L,PTR(J))
   80             CONTINUE
               END IF
            END IF
  100    CONTINUE
         DO 120 L = 1, NEQ
            YNEW(L) = Y(L) + HTRY*ERREST(L)
  120    CONTINUE
C
C  METHD = 2 is special in that an estimate of the local error can be
C  formed before the step is completed.  If the step is a failure,
C  return immediately.  Otherwise, complete the step and compute a more
C  accurate error estimate.
         IF (METHD.EQ.2 .AND. I.EQ.7) THEN
            CALL STEPB(NEQ,Y,YP,HTRY,YNEW,STAGES,THRES,ERR,MAIN,WEIGHT)
            IF (ERR.GT.TOL) RETURN
         END IF
C
         TSTG = TNOW + C(I)*HTRY
         IF (MAIN .AND. LAST .AND. C(I).EQ.ONE) TSTG = TND
         CALL F(TSTG,YNEW,STAGES(1,PTR(I)))
C
C  Increment the counter for the number of function evaluations
C  depending on whether the primary or secondary integration is taking
C  place.
         IF (MAIN) THEN
            NFCN = NFCN + 1
         ELSE
            GNFCN = GNFCN + 1
         END IF
C
C----------------------------------------------------------------------
C  When PHASE2 is .TRUE. we are in the second phase of the automatic
C  selection of the initial step size.  The results of the first three
C  stages are monitored in the subroutine STEPA for evidence that H is
C  too large -- instability and/or an unreliable estimate of the error
C  of the step is then possible.  When the subroutine believes H to be
C  too large, it returns CUTBAK = .TRUE. and a suitably reduced H for
C  another try.
C
         IF (MAIN) THEN
            IF (PHASE2) THEN
               IF (I.LE.3 .AND. ABS(HTRY).GT.HMIN) THEN
                  CALL STEPA(TNOW,Y,YP,TSTG,YNEW,STAGES(1,PTR(I)),
     &                       HTRY,WEIGHT,CUTBAK)
                  IF (CUTBAK) THEN
                     LAST = .FALSE.
C
C  Make sure that STEPA does not reduce the step size below the
C  minimum. If it does, reset H to HMIN and deactivate PHASE2.
                     IF (ABS(HTRY).LE.HMIN) THEN
                        HTRY = SIGN(HMIN,HTRY)
                        PHASE2 = .FALSE.
                     END IF
                     GO TO 20
                  END IF
               END IF
            END IF
         END IF
C----------------------------------------------------------------------
C
  140 CONTINUE
C
C  Some formulas are constructed so that the last stage represents
C  the result of the step (FSAL=.TRUE.), hence if the step is acceptable,
C  it will be the first stage for the next step. When FSAL=.FALSE., we
C  have to complete the computation of the step.
C
      IF (.NOT.FSAL) THEN
         DO 200 I = 1, NSTAGE
            IF (I.EQ.1) THEN
               DO 160 L = 1, NEQ
                  ERREST(L) = BHAT(1)*YP(L)
  160          CONTINUE
            ELSE
               IF (BHAT(I).NE.ZERO) THEN
                  DO 180 L = 1, NEQ
                     ERREST(L) = ERREST(L) + BHAT(I)*STAGES(L,PTR(I))
  180             CONTINUE
               END IF
            END IF
  200    CONTINUE
         DO 220 L = 1, NEQ
            YNEW(L) = Y(L) + HTRY*ERREST(L)
  220    CONTINUE
      END IF
C
C  Form an estimate of the error in the lower order formula by comparing
C  it to the higher order formula of the pair. ERREST(*) has been used
C  as working storage above.  The higher order approximation has been
C  formed as YNEW(*) = Y(*) + HTRY*ERREST(*) where ERREST(*) is a linear
C  combination of the stages of the formula. The lower order result also
C  has the form Y(*) plus HTRY times a different linear combination of
C  the stages. Hence, this different linear combination of stages for
C  the lower order formula can just be subtracted from the combination
C  stored in ERREST(*) to produce the errors. The result is then
C  multiplied by HTRY to obtain the error estimate.
C
      DO 280 I = 1, NSTAGE
         IF (I.EQ.1 .AND. B(1).NE.ZERO) THEN
            DO 240 L = 1, NEQ
               ERREST(L) = ERREST(L) - B(1)*YP(L)
  240       CONTINUE
         ELSE
            IF (B(I).NE.ZERO) THEN
               DO 260 L = 1, NEQ
                  ERREST(L) = ERREST(L) - B(I)*STAGES(L,PTR(I))
  260          CONTINUE
            END IF
         END IF
  280 CONTINUE
      DO 300 L = 1, NEQ
         ERREST(L) = HTRY*ERREST(L)
  300 CONTINUE
C
C  The error in a solution component is measured relative to a weight
C  that is the larger of a threshold and the size of the solution over
C  the step.  Using the magnitude of a solution component at both ends
C  of the step in the definition of "size" increases the robustness of
C  the test. When global error estimation is specified, the weight
C  vector WEIGHT(*) is defined by the primary integration and is then
C  used in the secondary integration.
C
      IF (MAIN) THEN
         DO 320 L = 1, NEQ
            AVGY = HALF*(ABS(Y(L))+ABS(YNEW(L)))
            WEIGHT(L) = MAX(AVGY,THRES(L))
  320    CONTINUE
      END IF
C
      ERR = ZERO
      DO 340 L = 1, NEQ
         ERR = MAX(ERR,ABS(ERREST(L)/WEIGHT(L)))
  340 CONTINUE   
C
      RETURN
      END
