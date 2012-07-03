      SUBROUTINE UT(F,TWANT,TGOT,YGOT,YPGOT,YMAX,WORK,UFLAG)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C  If you are not familiar with the code UT and how it is used in
C  conjunction with SETUP to solve initial value problems, you should study
C  the document file rksuite.doc carefully before proceeding further.  The
C  following "Brief Reminder" is intended only to remind you of the meaning,
C  type, and size requirements of the arguments.
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
C     Given input values of the independent variable T and the solution
C     components Y(*), for each L = 1,2,...,NEQ evaluate the differential
C     equation for the derivative of the Ith solution component and place the
C     value in YP(L).  Do not alter the input values of T and Y(*).
C  RETURN
C  END
C
C  INPUT VARIABLE
C
C     TWANT     - DOUBLE PRECISION
C                 The next value of the independent variable where a
C                 solution is desired.
C
C                 Constraints: TWANT must lie between the previous value
C                 of TGOT (TSTART on the first call) and TEND. TWANT can be
C                 equal to TEND, but it must be clearly distinguishable from
C                 the previous value of TGOT (TSTART on the first call) in 
C                 the precision available.
C
C  OUTPUT VARIABLES
C
C     TGOT      - DOUBLE PRECISION
C                 A solution has been computed at this value of the
C                 independent variable.
C     YGOT(*)   - DOUBLE PRECISION array of length NEQ
C                 Approximation to the true solution at TGOT. Do not alter
C                 the contents of this array
C     YPGOT(*)  - DOUBLE PRECISION array of length NEQ
C                 Approximation to the first derivative of the true
C                 solution at TGOT.
C     YMAX(*)   - DOUBLE PRECISION array of length NEQ
C                 YMAX(L) is the largest magnitude of YGOT(L) computed at any
C                 time in the integration from TSTART to TGOT. Do not alter
C                 the contents of this array.
C
C  WORKSPACE
C
C     WORK(*)   - DOUBLE PRECISION array as used in SETUP
C                 Do not alter the contents of this array.
C
C  OUTPUT VARIABLE
C
C     UFLAG     - INTEGER
C
C                       SUCCESS.  TGOT = TWANT.
C                 = 1 - Complete success.
C
C                       "SOFT" FAILURES
C                 = 2 - Warning:  You are using METHOD = 3 inefficiently
C                       by computing answers at many values of TWANT.  If
C                       you really need answers at so many specific points,
C                       it would be more efficient to compute them with
C                       METHOD = 2.  To do this you would need to restart
C                       from TGOT, YGOT(*) by a call to SETUP.  If you wish
C                       to continue as you are, you may.
C                 = 3 - Warning:  A considerable amount of work has been
C                       expended.  If you wish to continue on to TWANT, just
C                       call UT again.
C                 = 4 - Warning:  It appears that this problem is "stiff".
C                       You really should change to another code that is
C                       intended for such problems, but if you insist, you can
C                       continue with UT by calling it again.
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
      DOUBLE PRECISION  TGOT, TWANT
      INTEGER           UFLAG
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(*), YGOT(*), YMAX(*), YPGOT(*)
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
      PARAMETER         (SRNAME='UT')
      LOGICAL           ASK, TELL
      PARAMETER         (ASK=.TRUE.,TELL=.FALSE.)
      INTEGER           MINUS1, MINUS2
      PARAMETER         (MINUS1=-1,MINUS2=-2)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Local Scalars ..
      DOUBLE PRECISION  HMIN, TLAST, TNOW, UTEND
      INTEGER           CFLAG, IER, L, NREC, STATE
      LOGICAL           BADERR, GOBACK
C     .. External Subroutines ..
      EXTERNAL          CHKFL, CT, INTPO, RESET, RKMSG, RKSIT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Save statement ..
      SAVE              UTEND, TLAST
C     .. Executable Statements ..
      IER = 1
      NREC = 0
      GOBACK = .FALSE.
      BADERR = .FALSE.
C
C  Is it permissible to call UT?
C
      CALL RKSIT(ASK,'SETUP',STATE)
      IF (STATE.EQ.911) THEN
         IER = 912
         NREC = 1
         WRITE (REC,'(A)')
     &' ** A catastrophic error has already been detected elsewhere.'
         GO TO 100
      END IF
      IF (STATE.EQ.MINUS1) THEN
         IER = 911
         NREC = 1
         WRITE (REC,'(A)')
     &' ** You have not called SETUP, so you cannot use UT.'
         GO TO 100
      END IF
      IF (.NOT.UTASK) THEN
         IER = 911
         NREC = 2
         WRITE (REC,'(A/A)')
     &' ** You have called UT after you specified in SETUP that ',
     &' ** you were going to use CT. This is not permitted.'
         GO TO 100
      END IF
      CALL RKSIT(ASK,SRNAME,STATE)
      IF (STATE.EQ.5 .OR. STATE.EQ.6) THEN
         IER = 911
         NREC = 1
         WRITE (REC,'(A/A)')
     &' ** This routine has already returned with a hard failure.',
     &' ** You must call SETUP to start another problem.'
         GO TO 100
      END IF
      STATE = MINUS2
      CALL RKSIT(TELL,SRNAME,STATE)
C
      IF (FIRST) THEN
C
C  First call.
C
C  A value of TND is specified in SETUP. When INTP = .FALSE., as with
C  METHD = 3, output is obtained at the specified TWANT by resetting TND
C  to TWANT.  At this point, before the integration gets started, this can
C  be done with a simple assignment.  Later it is done with a call to RESET.
C  The original TND is SAVEd as a local variable UTEND.
C
         UTEND = TND
         IF (.NOT.INTP) TND = TWANT
C
C  The last TGOT returned is SAVEd in the variable TLAST.  T (a variable
C  passed through the common block RKCOM2) records how far the integration
C  has advanced towards the specified TND.  When output is obtained by
C  interpolation, the integration goes past the TGOT returned (T is closer
C  to the specified TND than TGOT).  Initialize these variables and YMAX(*).
         TLAST = TSTRT
         TGOT = TSTRT
         DO 20 L = 1, NEQN
            YMAX(L) = ABS(WORK(PRY-1+L))
   20    CONTINUE
C
C  If the code is to find an on-scale initial step size H, a bound was placed
C  on H in SETUP.  Here the first output point is used to refine this bound.
         IF (HSTRT.EQ.ZERO) THEN
            H = MIN(ABS(H),ABS(TWANT-TSTRT))
            HMIN = MAX(TINY,TOOSML*MAX(ABS(TSTRT),ABS(TND)))
            H = MAX(H,HMIN)
         END IF
C
      ELSE
C
C  Subsequent call.
C
         IF (TLAST.EQ.UTEND) THEN
            IER = 911
            NREC = 3
            WRITE (REC,'(A/A/A)')
     &' ** You have called UT after reaching TEND. (Your last    ',
     &' ** call to UT resulted in TGOT = TEND.)  To start a new ',
     &' ** problem, you will need to call SETUP.'
            GO TO 100
         END IF
C
      END IF
C
C  Check for valid TWANT.
C
      IF (DIR*(TWANT-TLAST).LE.ZERO) THEN
         IER = 911
         NREC = 4
         WRITE (REC,'(A/A/A/A)')
     &' ** You have made a call to UT with a TWANT that does   ',
     &' ** not lie between the previous value of TGOT (TSTART  ',
     &' ** on the first call) and TEND. This is not permitted. ',
     &' ** Check your program carefully.'
         GO TO 100
      END IF
      IF (DIR*(TWANT-UTEND).GT.ZERO) THEN
         HMIN = MAX(TINY,TOOSML*MAX(ABS(TWANT),ABS(UTEND)))
         IF (ABS(TWANT-UTEND).LT.HMIN) THEN
            IER = 911
            NREC = 5
            WRITE (REC,'(A/A/A/A)')
     &' ** You have made a call to UT with a TWANT that does      ',
     &' ** not lie between the previous value of TGOT (TSTART on  ',
     &' ** the first call) and TEND. This is not permitted. TWANT ',
     &' ** is very close to TEND, so you may have meant to set    ',
     &' ** it to be TEND exactly.  Check your program carefully.  '
         ELSE
            IER = 911
            NREC = 4
            WRITE (REC,'(A/A/A/A)')
     &' ** You have made a call to UT with a TWANT that does   ',
     &' ** not lie between the previous value of TGOT (TSTART  ',
     &' ** on the first call) and TEND. This is not permitted. ',
     &' ** Check your program carefully.'
         END IF
         GO TO 100
      END IF
      IF (.NOT.INTP) THEN
         HMIN = MAX(TINY,TOOSML*MAX(ABS(TLAST),ABS(TWANT)))
         IF (ABS(TWANT-TLAST).LT.HMIN) THEN
            IER = 911
            NREC = 4
            WRITE (REC,'(A/A/A/A,D13.5,A)')
     &' ** You have made a call to UT with a TWANT that is not ',
     &' ** sufficiently different from the last value of TGOT  ',
     &' ** (TSTART on the first call).  When using METHOD = 3, ',
     &' ** it must differ by at least ',HMIN,'.'
            GO TO 100
         END IF
C
C  We have a valid TWANT. There is no interpolation with this METHD and
C  therefore we step to TWANT exactly by resetting TND with a call to RESET.
C  On the first step this matter is handled differently as explained above.
C
         IF (.NOT.FIRST) THEN
            CALL RESET(TWANT)
            CALL CHKFL(ASK,BADERR)
            IF (BADERR) GO TO 100
         END IF
      END IF
C
C  Process output, decide whether to take another step.
C
   40 CONTINUE
C
      IF (INTP) THEN
C
C  Interpolation is possible with this METHD.  The integration has
C  already reached T. If this is past TWANT, GOBACK is set .TRUE. and
C  the answers are obtained by interpolation.
C
         GOBACK = DIR*(T-TWANT) .GE. ZERO
         IF (GOBACK) THEN
            CALL INTPO(TWANT,'Both solution and derivative',NEQN,YGOT,
     &                 YPGOT,F,WORK,WORK(PRINTP),LNINTP)
            CALL CHKFL(ASK,BADERR)
            IF (BADERR) GO TO 100
            TGOT = TWANT
         END IF
      ELSE
C
C  Interpolation is not possible with this METHD, so output is obtained
C  by integrating to TWANT = TND.  Both YGOT(*) and YPGOT(*) are then 
C  already loaded with the solution at TWANT by CT.
C
         GOBACK = T .EQ. TWANT
         IF (GOBACK) TGOT = TWANT
      END IF
C
C  Updating of YMAX(*) is done here to account for the fact that when
C  interpolation is done, the integration goes past TGOT.  Note that YGOT(*)
C  is not defined until CT is called.  YMAX(*) was initialized at TSTRT
C  from values stored in WORK(*), so only needs to be updated for T
C  different from TSTRT.
      IF (T.NE.TSTRT) THEN
         DO 60 L = 1, NEQN
            YMAX(L) = MAX(YMAX(L),ABS(YGOT(L)))
   60    CONTINUE
      END IF
C
C  If done, go to the exit point.
      IF (GOBACK) GO TO 100
C
C  Take a step with CT in the direction of TND.  On exit, the solution is
C  advanced to TNOW.  The way CT is written, the approximate solution at
C  TNOW is available in both YGOT(*) and in WORK(*).  If output is obtained by
C  stepping to the end (TNOW = TWANT = TND), YGOT(*) can be returned directly.
C  If output is obtained by interpolation, the subroutine INTPO that does this
C  uses the values in WORK(*) for its computations and places the approximate
C  solution at TWANT in the array YGOT(*) for return to the calling program.
C  The approximate derivative is handled in the same way. TNOW is output from
C  CT and is actually a copy of T declared above in a common block.
C
      CALL CT(F,TNOW,YGOT,YPGOT,WORK,CFLAG)
      IER = CFLAG
C
C  A successful step by CT is indicated by CFLAG = 1 or = 2.
      IF (CFLAG.EQ.1) THEN
         GO TO 40
      ELSE IF (CFLAG.EQ.2) THEN
C
C  Supplement the warning message written in CT.
         NREC = 3
         WRITE (REC,'(A/A/A)')
     &' ** The last message was produced on a call to CT from UT.  ',
     &' ** In UT the appropriate action is to change to METHOD = 2,',
     &' ** or, if insufficient memory is available, to METHOD = 1. '
      ELSE IF (CFLAG.LE.6) THEN
         NREC = 1
         WRITE (REC,'(A)')
     &' ** The last message was produced on a call to CT from UT.'
      ELSE
         BADERR = .TRUE.
      END IF
      TGOT = T
C
C  Update YMAX(*) before the return.
      DO 80 L = 1, NEQN
         YMAX(L) = MAX(YMAX(L),ABS(YGOT(L)))
   80 CONTINUE
C
C  Exit point for UT.
C
  100 CONTINUE
C
      IF (BADERR) THEN
         IER = 911
         NREC = 4
         WRITE (REC,'(A/A/A/A)')
     &' ** An internal call by UT to a subroutine resulted in an  ',
     &' ** error that should not happen.  Check your program      ',
     &' ** carefully for array sizes, correct number of arguments,',
     &' ** type mismatches ... .'
      END IF
C
      TLAST = TGOT
C
C  All exits are done here after a call to RKMSG to report
C  what happened and set UFLAG.
C
      CALL RKMSG(IER,SRNAME,NREC,UFLAG)
C
      RETURN
      END
