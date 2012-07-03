      SUBROUTINE CT(F,TNOW,YNOW,YPNOW,WORK,CFLAG)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C  If you are not familiar with the code CT and how it is used in
C  conjunction with SETUP to solve initial value problems, you should study
C  the document file rksuite.doc carefully before attempting to use the code.
C  The following "Brief Reminder" is intended only to remind you of the
C  meaning, type, and size requirements of the arguments.
C
C  NAME DECLARED IN AN EXTERNAL STATEMENT IN THE CALLING PROGRAM:
C
C     F         - name of the subroutine for evaluating the differential
C                 equations.
C
C  The subroutine F must have the form
C
C  SUBROUTINE F(T,Y,YP)
C  DOUBLE PRECISION T,Y(*),YP(*)
C     Using the input values of the independent variable T and the solution
C     components Y(*), for each L = 1,2,...,NEQ evaluate the differential
C     equation for the derivative of the Lth solution component and place the
C     value in YP(L).  Do not alter the input values of T and Y(*).
C  RETURN
C  END
C
C  OUTPUT VARIABLES
C
C     TNOW      - DOUBLE PRECISION
C                 Current value of the independent variable.
C     YNOW(*)   - DOUBLE PRECISION array of length NEQ
C                 Approximation to the true solution at TNOW.
C     YPNOW(*)  - DOUBLE PRECISION array of length NEQ
C                 Approximation to the first derivative of the
C                 true solution at TNOW.
C
C  WORKSPACE
C
C     WORK(*)   - DOUBLE PRECISION array as used in SETUP
C                 Do not alter the contents of this array.
C
C  OUTPUT VARIABLE
C
C     CFLAG     - INTEGER
C
C                       SUCCESS.  A STEP WAS TAKEN TO TNOW.
C                 = 1 - Complete success.
C
C                       "SOFT" FAILURES
C                 = 2 - Warning:  You have obtained an answer by integrating
C                       to TEND (TNOW = TEND).  You have done this at least
C                       100 times, and monitoring of the computation reveals
C                       that this way of getting output has degraded the
C                       efficiency of the code. If you really need answers at
C                       so many specific points, it would be more efficient to
C                       get them with INTPO.  (If METHOD = 3, you would need
C                       to change METHOD and restart from TNOW, YNOW(*) by a
C                       call to SETUP.)  If you wish to continue as you are,
C                       you may.
C                 = 3 - Warning:  A considerable amount of work has been
C                       expended. To continue the integration, just call
C                       CT again.
C                 = 4 - Warning:  It appears that this problem is "stiff".
C                       You really should change to another code that is
C                       intended for such problems, but if you insist, you 
C                       can continue with CT by calling it again.
C
C                       "HARD" FAILURES
C                 = 5 - You are asking for too much accuracy. You cannot
C                       continue integrating this problem.
C                 = 6 - The global error assessment may not be reliable beyond
C                       the current point in the integration.  You cannot
C                       continue integrating this problem.
C
C                       "CATASTROPHIC" FAILURES
C                 = 911 - The nature of the catastrophe is reported on
C                         the standard output channel. Unless special
C                         provision was made in advance (see rksuite.doc),
C                         the computation then comes to a STOP.
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TNOW
      INTEGER           CFLAG
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(*), YNOW(*), YPNOW(*)
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
C     .. Common Block for General Workspace Pointers ..
      INTEGER           PRTHRS, PRERST, PRWT, PRYOLD, PRSCR, PRY, PRYP,
     &                  PRSTGS, PRINTP, LNINTP
      COMMON /RKCOM3/   PRTHRS, PRERST, PRWT, PRYOLD, PRSCR, PRY, PRYP,
     &                  PRSTGS, PRINTP, LNINTP
      SAVE   /RKCOM3/
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
C     .. Common Block for Integrator Options ..
      LOGICAL           MSG, UTASK
      COMMON /RKCOM8/   MSG, UTASK
      SAVE   /RKCOM8/
C     .. Common Block for Error Message ..
      CHARACTER*80      REC(10)
      COMMON /RKCOM9/   REC
      SAVE   /RKCOM9/
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='CT')
      LOGICAL           ASK, TELL
      PARAMETER         (ASK=.TRUE.,TELL=.FALSE.)
      INTEGER           MINUS1, MINUS2
      PARAMETER         (MINUS1=-1,MINUS2=-2)
      INTEGER           MAXFCN
      PARAMETER         (MAXFCN=5000)
      DOUBLE PRECISION  ZERO, PT1, PT9, ONE, TWO, HUNDRD
      PARAMETER         (ZERO=0.0D+0,PT1=0.1D+0,PT9=0.9D+0,ONE=1.0D+0,
     &                  TWO=2.0D+0,HUNDRD=100.0D+0)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, BETA, ERR, ERROLD, HAVG, HMIN, HTRY, TAU,
     &                  TEMP1, TEMP2, YPNORM
      INTEGER           IER, JFLSTP, L, NREC, NTEND, POINT, STATE, YNEW,
     &                  YPOLD
      LOGICAL           CHKEFF, FAILED, MAIN, PHASE1, PHASE2, PHASE3,
     &                  TOOMCH
C     .. External Subroutines ..
      EXTERNAL          RKMSG, RKSIT, STEP, STIFF, TRUERR
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SIGN
C     .. Save statement ..
      SAVE              JFLSTP, NTEND, ERROLD, HAVG, PHASE2, YNEW,
     &                  YPOLD, CHKEFF
C     .. Executable Statements ..
C
      IER = 1
      NREC = 0
C
C  Is it permissible to call CT?
C
      CALL RKSIT(ASK,'SETUP',STATE)
      IF (STATE.EQ.911) THEN
         IER = 912
         NREC = 1
         WRITE (REC,'(A)')
     &' ** A catastrophic error has already been detected elsewhere.'
         GO TO 180
      END IF
      IF (STATE.EQ.MINUS1) THEN
         IER = 911
         NREC = 1
         WRITE (REC,'(A)')
     &' ** You have not called SETUP, so you cannot use CT.'
         GO TO 180
      END IF
      IF (UTASK) THEN
         CALL RKSIT(ASK,'UT',STATE)
         IF (STATE.NE.MINUS2) THEN
            IER = 911
            NREC = 2
            WRITE (REC,'(A/A)')
     &' ** You have called CT after you specified in SETUP that ',
     &' ** you were going to use UT. This is not permitted.'
            UTASK = .FALSE.
            GO TO 180
         END IF
      END IF
      CALL RKSIT(ASK,SRNAME,STATE)
      IF (STATE.EQ.5 .OR. STATE.EQ.6) THEN
         IER = 911
         NREC = 3
         WRITE (REC,'(A/A)')
     &' ** CT has already returned with a flag value of 5 or 6.',
     &' ** You cannot continue integrating this problem. You must ',
     &' ** call SETUP to start another problem.'
         GO TO 180
      END IF
C
      IF (FIRST) THEN
C
C  First call in an integration -- initialize everything.
C
         CHKEFF = .FALSE.
         NTEND = 0
         JFLSTP = 0
C
C  A scratch area of WORK(*) starting at PRSCR is used to hold two
C  arrays in this subroutine: the higher order approximate solution at
C  the end of a step and the approximate derivative of the solution at
C  the end of the last step. To make this clear, local pointers YNEW and
C  YPOLD are used.
         YNEW = PRSCR
         YPOLD = PRSCR
C
C  For this first step T was initialized to TSTRT in SETUP and the
C  starting values YSTART(*) were loaded into the area of WORK(*) reserved
C  for the current solution approximation starting at location PRY. The
C  derivative is now computed and stored in WORK(*) starting at PRYP.
C  Subsequently these arrays are copied to the output vectors YNOW(*)
C  and YPNOW(*).
         CALL F(T,WORK(PRY),WORK(PRYP))
         NFCN = NFCN + 1
         DO 20 L = 1, NEQN
            YNOW(L) = WORK(PRY-1+L)
            YPNOW(L) = WORK(PRYP-1+L)
   20    CONTINUE
C
C  Set dependent variables for error assessment.
         IF (ERASON) THEN
            DO 40 L = 1, NEQN
               WORK(PRZY-1+L) = YNOW(L)
               WORK(PRZYP-1+L) = YPNOW(L)
   40       CONTINUE
         END IF
C
C  The weights for the control of the error depend on the size of the
C  solution at the beginning and at the end of the step. On the first
C  step we do not have all this information. Whilst determining the
C  initial step size we initialize the weight vector to the larger of
C  abs(Y(L)) and the threshold for this component.
         DO 60 L = 1, NEQN
            WORK(PRWT-1+L) = MAX(ABS(YNOW(L)),WORK(PRTHRS-1+L))
   60    CONTINUE
C
C  If HSTRT is equal to zero, the code is to find an on-scale initial step
C  size H.  CT has an elaborate scheme of three phases for finding such an H,
C  and some preparations are made earlier.  In SETUP an upper bound is placed
C  on H that reflects the scale of the independent variable. When UTASK is
C  .TRUE., UT refines this bound using the first output point.  Here in CT
C  PHASE1 applies a rule of thumb based on the error control, the order of the
C  the formula, and the size of the initial slope to get a crude approximation
C  to an on-scale H.  PHASE2 may reduce H in the course of taking the first
C  step.  PHASE3 repeatedly adjusts H and retakes the first step until H is
C  on scale.
C
C  A guess for the magnitude of the first step size H can be provided to SETUP
C  as HSTART.  If it is too big or too small, it is ignored and the automatic
C  determination of an on-scale initial step size is activated.  If it is
C  acceptable, H is set to HSTART in SETUP.  Even when H is supplied to CT,
C  PHASE3 of the scheme for finding an on-scale initial step size is made
C  active so that the code can deal with a bad guess.
C
         PHASE1 = HSTRT .EQ. ZERO
         PHASE2 = PHASE1
         PHASE3 = .TRUE.
         IF (PHASE1) THEN
            H = ABS(H)
            YPNORM = ZERO
            DO 80 L = 1, NEQN
               IF (ABS(YNOW(L)).NE.ZERO) THEN
                  YPNORM = MAX(YPNORM,ABS(YPNOW(L))/WORK(PRWT-1+L))
               END IF
   80       CONTINUE
            TAU = TOLR**EXPON
            IF (H*YPNORM.GT.TAU) H = TAU/YPNORM
            HMIN = MAX(TINY,TOOSML*MAX(ABS(TSTRT),ABS(TND)))
            H = DIR*MAX(H,HMIN)
            PHASE1 = .FALSE.
         END IF
C
      ELSE
C
C Continuation call
C
         IF (LAST) THEN
            IER = 911
            NREC = 3
            WRITE (REC,'(A,D13.5,A/A/A)')
     &' ** You have already reached TEND ( = ', TND, ').',
     &' ** To integrate further with the same problem you must ',
     &' ** call the routine RESET with a new value of TEND.'
            GO TO 180
         END IF
      END IF
C
C  Begin computation of a step here.
C
      FAILED = .FALSE.
C
  100 CONTINUE
      H = SIGN(ABS(H),DIR)
C
C  Reduce the step size if necessary so that the code will not step
C  past TND.  "Look ahead" to prevent unnecessarily small step sizes.
      LAST = DIR*((T+H)-TND) .GE. ZERO
      IF (LAST) THEN
         H = TND - T
      ELSE IF (DIR*((T+TWO*H)-TND).GE.ZERO) THEN
         H = (TND-T)/TWO
      END IF
C
C  When the integrator is at T and attempts a step of H, the function
C  defining the differential equations will be evaluated at a number of
C  arguments between T and T+H.  If H is too small, these arguments cannot
C  be clearly distinguished in the precision available.
C
      HMIN = MAX(TINY,TOOSML*MAX(ABS(T),ABS(T+H)))
      IF (ABS(H).LT.HMIN) THEN
         IER = 5
         NREC = 3
         WRITE (REC,'(A/A,D13.5,A,D13.5/A)')
     &' ** In order to satisfy your error requirements CT would ',
     &' ** have to use a step size of ',H,' at TNOW = ',T,
     &' ** This is too small for the machine precision.'
         GO TO 180
      END IF
C
C  Monitor the impact of output on the efficiency of the integration.
C
      IF (CHKEFF) THEN
         NTEND = NTEND + 1
         IF (NTEND.GE.100 .AND. NTEND.GE.OKSTP/3) THEN
            IER = 2
            NREC = 6
            WRITE (REC,'(A/A/A/A/A/A)')
     &' ** More than 100 output points have been obtained by ',
     &' ** integrating to TEND.  They have been sufficiently close ',
     &' ** to one another that the efficiency of the integration has ',
     &' ** been degraded. It would probably be (much) more efficient ',
     &' ** to obtain output by interpolating with INTPO (after ',
     &' ** changing to METHOD=2 if you are using METHOD = 3).'
            NTEND = 0
            GO TO 180
         END IF
      END IF
C
C  Check for stiffness and for too much work.  Stiffness can be
C  checked only after a successful step.
C
      IF (.NOT.FAILED) THEN
C
C  Check for too much work.
         TOOMCH = NFCN .GT. MAXFCN
         IF (TOOMCH) THEN
            IER = 3
            NREC = 3
           WRITE (REC,'(A,I6,A/A/A)')
     &' ** Approximately ',MAXFCN,' function evaluations have been ',
     &' ** used to compute the solution since the integration ',
     &' ** started or since this message was last printed.'
C
C  After this warning message, NFCN is reset to permit the integration
C  to continue.  The total number of function evaluations in the primary
C  integration is SVNFCN + NFCN.
            SVNFCN = SVNFCN + NFCN
            NFCN = 0
         END IF
C
C  Check for stiffness.  NREC is passed on to STIFF because when
C  TOOMCH = .TRUE. and stiffness is diagnosed, the message about too
C  much work is augmented inside STIFF to explain that it is due to
C  stiffness.
         CALL STIFF(F,HAVG,JFLSTP,TOOMCH,MAXFCN,WORK,IER,NREC)
C
         IF (IER.NE.1) GO TO 180
      END IF
C
C  Take a step.  Whilst finding an on-scale H (PHASE2 = .TRUE.), the input
C  value of H might be reduced (repeatedly), but it will not be reduced
C  below HMIN.  The local error is estimated, a weight vector is formed,
C  and a weighted maximum norm, ERR, of the local error is returned.
C  The variable MAIN is input as .TRUE. to tell STEP that this is the
C  primary, or "main", integration.
C
C  H resides in the common block /RKCOM2/ which is used by both CT and STEP;
C  since it may be changed inside STEP, a local copy is made to ensure
C  portability of the code.
C
      MAIN = .TRUE.
      HTRY = H
      CALL STEP(F,NEQN,T,WORK(PRY),WORK(PRYP),WORK(PRSTGS),TOLR,HTRY,
     &          WORK(PRWT),WORK(YNEW),WORK(PRERST),ERR,MAIN,HMIN,
     &          WORK(PRTHRS),PHASE2)
      H = HTRY
C
C  Compare the norm of the local error to the tolerance.
C
      IF (ERR.GT.TOLR) THEN
C
C  Failed step.  Reduce the step size and try again.
C
C  First step:  Terminate PHASE3 of the search for an on-scale step size.
C               The step size is not on scale, so ERR may not be accurate;
C               reduce H by a fixed factor.  Failed attempts to take the
C               first step are not counted.
C  Later step:  Use ERR to compute an "optimal" reduction of H.  More than
C               one failure indicates a difficulty with the problem and an
C               ERR that may not be accurate, so reduce H by a fixed factor.
C
         IF (FIRST) THEN
            PHASE3 = .FALSE.
            ALPHA = RS1
         ELSE
            FLSTP = FLSTP + 1
            JFLSTP = JFLSTP + 1
            IF (FAILED) THEN
               ALPHA = RS1
            ELSE
               ALPHA = SAFETY*(TOLR/ERR)**EXPON
               ALPHA = MAX(ALPHA,RS1)
            END IF
         END IF
         H = ALPHA*H
         FAILED = .TRUE.
         GO TO 100
      END IF
C
C  Successful step.
C
C  Predict a step size appropriate for the next step.  After the first
C  step the prediction can be refined using an idea of H.A. Watts that
C  takes account of how well the prediction worked on the previous step.
      BETA = (ERR/TOLR)**EXPON
      IF (.NOT.FIRST) THEN
         TEMP1 = (ERR**EXPON)/H
         TEMP2 = (ERROLD**EXPON)/HOLD
         IF (TEMP1.LT.TEMP2*HUNDRD .AND. TEMP2.LT.TEMP1*HUNDRD) THEN
            BETA = BETA*(TEMP1/TEMP2)
         END IF
      END IF
      ALPHA = RS3
      IF (SAFETY.LT.BETA*ALPHA) ALPHA = SAFETY/BETA
C
C  On the first step a search is made for an on-scale step size.  PHASE2
C  of the scheme comes to an end here because a step size has been found
C  that is both successful and has a credible local error estimate. Except
C  in the special case that the first step is also the last, the step is
C  repeated in PHASE3 as long as an increase greater than RS2 appears
C  possible.  An increase as big as RS3 is permitted.  A step failure
C  terminates PHASE3.
C
      IF (FIRST) THEN
         PHASE2 = .FALSE.
         PHASE3 = PHASE3 .AND. .NOT. LAST .AND. (ALPHA.GT.RS2)
         IF (PHASE3) THEN
            H = ALPHA*H
            GO TO 100
         END IF
      END IF
C
C  After getting on scale, step size changes are more restricted.
      ALPHA = MIN(ALPHA,RS)
      IF (FAILED) ALPHA = MIN(ALPHA,ONE)
      ALPHA = MAX(ALPHA,RS1)
      HOLD = H
      H = ALPHA*H
C
C  For the diagnosis of stiffness, an average accepted step size, HAVG,
C  must be computed and SAVEd.
      IF (FIRST) THEN
         HAVG = HOLD
      ELSE
         HAVG = PT9*HAVG + PT1*HOLD
      END IF
C
      FIRST = .FALSE.
      ERROLD = ERR
      TOLD = T
C
C  Take care that T is set to precisely TND when the end of the
C  integration is reached.
      IF (LAST) THEN
         T = TND
      ELSE
         T = T + HOLD
      END IF
C
C  Increment counter on accepted steps.  Note that successful steps
C  that are repeated whilst getting on scale are not counted.
      OKSTP = OKSTP + 1
C
C  Advance the current solution and its derivative.  (Stored in WORK(*)
C  with the first location being PRY and PRYP, respectively.)  Update the
C  previous solution and its derivative.  (Stored in WORK(*) with the first
C  location being PRYOLD and YPOLD, respectively.)  Note that the previous
C  derivative will overwrite YNEW(*).
C
      DO 120 L = 1, NEQN
         WORK(PRYOLD-1+L) = WORK(PRY-1+L)
         WORK(PRY-1+L) = WORK(YNEW-1+L)
         WORK(YPOLD-1+L) = WORK(PRYP-1+L)
  120 CONTINUE
C
      IF (FSAL) THEN
C
C  When FSAL = .TRUE., YP(*) is the last stage of the step.
         POINT = PRSTGS + (LSTSTG-1)*NEQN
         DO 140 L = 1, NEQN
            WORK(PRYP-1+L) = WORK(POINT-1+L)
  140    CONTINUE
      ELSE
C
C  Call F to evaluate YP(*).
         CALL F(T,WORK(PRY),WORK(PRYP))
         NFCN = NFCN + 1
      END IF
C
C  If global error assessment is desired, advance the secondary
C  integration from TOLD to T.
C
      IF (ERASON) THEN
         CALL TRUERR(F,NEQN,WORK(PRY),TOLR,WORK(PRWT),WORK(PRZY),
     &               WORK(PRZYP),WORK(PRZERR),WORK(PRZYNU),WORK(PRZERS),
     &               WORK(PRZSTG),IER)
         IF (IER.EQ.6) THEN
C
C  The global error estimating procedure has broken down. Treat it as a
C  failed step. The solution and derivative are reset to their values at
C  the beginning of the step since the last valid error assessment refers
C  to them.
            OKSTP = OKSTP - 1
            ERASFL = .TRUE.
            LAST = .FALSE.
            T = TOLD
            H = HOLD
            DO 160 L = 1, NEQN
               WORK(PRY-1+L) = WORK(PRYOLD-1+L)
               WORK(PRYP-1+L) = WORK(YPOLD-1+L)
  160       CONTINUE
            IF (OKSTP.GT.1) THEN
               NREC = 2
               WRITE (REC,'(A/A,D13.5,A)')
     &' ** The global error assessment may not be reliable for T past ',
     &' ** TNOW = ',T,'.  The integration is being terminated.'
            ELSE
               NREC = 2
               WRITE (REC,'(A/A)')
     &' ** The global error assessment algorithm failed at the start',
     &' ** the integration.  The integration is being terminated.'
            END IF
            GO TO 180
         END IF
      END IF
C
C
C  Exit point for CT
C
  180 CONTINUE
C
C  Set the output variables and flag that interpolation is permitted
C
      IF (IER.LT.911) THEN
         TNOW = T
         LAST = TNOW .EQ. TND
         CHKEFF = LAST
         DO 200 L = 1, NEQN
            YNOW(L) = WORK(PRY-1+L)
            YPNOW(L) = WORK(PRYP-1+L)
  200    CONTINUE
         IF (IER.EQ.1) THEN
            STATE = MINUS2
            CALL RKSIT(TELL,'INTPO',STATE)
         END IF
      END IF
C
C  Call RKMSG to report what happened and set CFLAG.
C
      CALL RKMSG(IER,SRNAME,NREC,CFLAG)
C
      RETURN
      END
