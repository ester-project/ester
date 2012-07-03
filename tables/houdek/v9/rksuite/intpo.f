      SUBROUTINE INTPO(TWANT,REQEST,NWANT,YWANT,YPWANT,F,WORK,WRKINT,
     &                 LENINT)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C  If you are not familiar with the code INTPO and how it is used in
C  conjunction with CT to solve initial value problems, you should study the
C  document file rksuite.doc carefully before attempting to use the code. The
C  following "Brief Reminder" is intended only to remind you of the meaning,
C  type, and size requirements of the arguments.
C
C  When integrating with METHOD = 1 or 2, answers may be obtained inexpensively
C  between steps by interpolation. INTPO is called after a step by CT from a
C  previous value of T, called TOLD below, to the current value of T to get
C  an answer at TWANT. You can specify any value of TWANT you wish, but
C  specifying a value outside the interval [TOLD,T] is likely to yield
C  answers with unsatisfactory accuracy.
C
C  INPUT VARIABLE
C
C     TWANT     - DOUBLE PRECISION
C                 The value of the independent variable where a solution
C                 is desired.
C
C  The interpolant is to be evaluated at TWANT to approximate the solution
C  and/or its first derivative there.  There are three cases:
C
C  INPUT VARIABLE
C
C     REQEST    - CHARACTER*(*)
C                 Only the first character of REQEST is significant.
C                 REQEST(1:1) = `S' or `s'- compute approximate `S'olution 
C                                           only.
C                             = `D' or `d'- compute approximate first 
C                                           `D'erivative of the solution only.
C                             = `B' or `b'- compute `B'oth approximate solution
C                                           and first derivative.
C                 Constraint: REQEST(1:1) must be `S',`s',`D',`d',`B' or `b'.
C
C  If you intend to interpolate at many points, you should arrange for the
C  the interesting components to be the first ones because the code
C  approximates only the first NWANT components.
C
C  INPUT VARIABLE
C
C     NWANT     - INTEGER
C                 Only the first NWANT components of the answer are to be
C                 computed.
C                 Constraint:  NEQ >= NWANT >= 1
C
C  OUTPUT VARIABLES
C
C     YWANT(*)  - DOUBLE PRECISION array of length NWANT
C                 Approximation to the first NWANT components of the true
C                 solution at TWANT when REQESTed.
C     YPWANT(*) - DOUBLE PRECISION array of length NWANT
C                 Approximation to the first NWANT components of the first
C                 derivative of the true solution at TWANT when REQESTed.
C
C  NAME DECLARED IN AN EXTERNAL STATEMENT IN THE PROGRAM CALLING INTPO:
C
C     F         - name of the subroutine for evaluating the differential
C                 equations as provided to CT.
C
C  WORKSPACE
C
C     WORK(*)   - DOUBLE PRECISION array as used in SETUP and CT
C                 Do not alter the contents of this array.
C
C     WRKINT(*) - DOUBLE PRECISION array of length LENINT
C                 Do not alter the contents of this array.
C
C     LENINT    - INTEGER
C                 Length of WRKINT. If
C                 METHOD = 1, LENINT must be at least 1
C                        = 2, LENINT must be at least NEQ+MAX(NEQ,5*NWANT)
C                        = 3--CANNOT BE USED WITH THIS SUBROUTINE
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TWANT
      INTEGER           LENINT, NWANT
      CHARACTER*(*)     REQEST
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(*), WRKINT(*), YPWANT(*), YWANT(*)
C     .. Subroutine Arguments ..
      EXTERNAL          F
C     .. Common Block for Problem Definition ..
      DOUBLE PRECISION  TSTRT, TND, DIR, HSTRT, TOLR
      INTEGER           NEQN
      COMMON /RKCOM1/   TSTRT, TND, DIR, HSTRT, TOLR, NEQN
      SAVE   /RKCOM1/
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
      PARAMETER         (SRNAME='INTPO')
      LOGICAL           ASK
      INTEGER           PLUS1, MINUS1, MINUS2
      PARAMETER         (ASK=.TRUE.,PLUS1=1,MINUS1=-1,MINUS2=-2)
C     .. Local Scalars ..
      INTEGER           FLAG, ICHK, IER, NREC, NWNTSV, STARTP, STATE,
     &                  STATE1
      LOGICAL           ININTP, LEGALR
      CHARACTER         REQST1
C     .. External Subroutines ..
      EXTERNAL          EVALI, FORMI, RKMSG, RKSIT
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Save statement ..
      SAVE              NWNTSV, ININTP, STARTP
C     .. Data statements ..
      DATA              NWNTSV/MINUS1/
C     .. Executable Statements ..
C
      IER = 1
      NREC = 0
C
C  Is it permissible to call INTPO?
C
      CALL RKSIT(ASK,'CT',STATE)
      IF (STATE.EQ.911) THEN
         IER = 912
         NREC = 1
         WRITE (REC,'(A)')
     &' ** A catastrophic error has already been detected elsewhere.'
         GO TO 20
      END IF
      IF (UTASK) THEN
         CALL RKSIT(ASK,'UT',STATE1)
         IF (STATE1.NE.MINUS2) THEN
            IER = 911
            NREC = 2
            WRITE (REC,'(A/A)')
     &' ** You have called INTPO after you specified to SETUP ',
     &' ** that you were going to use UT. This is not permitted.'
            GO TO 20
         END IF
      END IF
      IF (STATE.EQ.MINUS1) THEN
         IER = 911
         NREC = 1
         WRITE (REC,'(A)')
     &' ** You have not called CT, so you cannot use INTPO.'
         GO TO 20
      END IF
      IF (STATE.GT.PLUS1) THEN
         IER = 911
         NREC = 2
         WRITE (REC,'(A/A)')
     &' ** CT has returned with a flag value greater than 1.',
     &' ** You cannot call INTPO in this circumstance.'
         GO TO 20
      END IF
C
C  Check input
C
      REQST1 = REQEST(1:1)
      LEGALR = REQST1 .EQ. 'S' .OR. REQST1 .EQ. 's' .OR.
     &         REQST1 .EQ. 'D' .OR. REQST1 .EQ. 'd' .OR.
     &         REQST1 .EQ. 'B' .OR. REQST1 .EQ. 'b'
      IF (.NOT.LEGALR) THEN
         IER = 911
         NREC = 3
         WRITE (REC,'(A/A,A,A/A)')
     &' ** You have set the first character of ',
     &' ** REQEST to be ''',REQST1,'''. It must be one of ',
     &' ** ''S'',''s'',''D'',''d'',''B'' or ''b''.'
         GO TO 20
      END IF
C
      IF (NWANT.GT.NEQN) THEN
         IER = 911
         NREC = 3
         WRITE (REC,'(A,I6,A/A,I6,A/A)')
     &' ** You have specified the value of NWANT to be ',NWANT,'. This',
     &' ** is greater than ',NEQN,', which is the number of equations ',
     &' ** in the system being integrated.'
         GO TO 20
      ELSE IF (NWANT.LT.1) THEN
         IER = 911
         NREC = 3
         WRITE (REC,'(A,I6,A/A/A)')
     &' ** You have specified the value of NWANT to be ',NWANT,', but ',
     &' ** this is less than 1. You cannot interpolate a zero or ',
     &' ** negative number of components.'
         GO TO 20
      END IF
C
      IF (METHD.EQ.1) THEN
         IF (LENINT.LT.1) THEN
            IER = 911
            NREC = 2
            WRITE (REC,'(A,I6,A/A)')
     &' ** You have specified LENINT to be ',LENINT,'.',
     &' ** This is too small. LENINT must be at least 1.'
            GO TO 20
         END IF
         STARTP = 1
      ELSE IF (METHD.EQ.2) THEN
         ICHK = NEQN + MAX(NEQN,5*NWANT)
         IF (LENINT.LT.ICHK) THEN
            IER = 911
            NREC = 3
            WRITE (REC,'(A,I6,A/A/A,I6,A)')
     &' ** You have specified LENINT to be ',LENINT,'. This is too',
     &' ** small. NINT must be at least NEQ + MAX(NEQ, 5*NWANT) ',
     &' ** which is ', ICHK,'.'
            GO TO 20
         END IF
         STARTP = NEQN + 1
      ELSE IF (METHD.EQ.3) THEN
         IER = 911
         NREC = 5
         WRITE (REC,'(A/A/A/A/A)')
     &' ** You have been using CT with METHOD = 3 to integrate your  ',
     &' ** equations. You have just called INTPO, but interpolation  ',
     &' ** is not available for this METHOD. Either use METHOD = 2,  ',
     &' ** for which interpolation is available, or use RESET to make',
     &' ** CT step exactly to the points where you want output.'
         GO TO 20
      END IF
C
C  Has the interpolant been initialised for this step?
C
      CALL RKSIT(ASK,SRNAME,STATE)
      ININTP = STATE .NE. MINUS2
C
C  Some initialization must be done before interpolation is possible.
C  To reduce the overhead, the interpolating polynomial is formed for
C  the first NWANT components.  In the unusual circumstance that NWANT
C  is changed while still interpolating within the span of the current
C  step, the scheme must be reinitialized to accomodate the additional
C  components.
C
      IF (.NOT.ININTP .OR. NWANT.NE.NWNTSV) THEN
C
C  At present the derivative of the solution at the previous step, YPOLD(*),
C  is stored in the scratch array area starting at PRSCR. In the case of
C  METHD = 1 we can overwrite the stages.
C
         IF (METHD.EQ.1) THEN
            CALL FORMI(F,NEQN,NWANT,WORK(PRY),WORK(PRYP),WORK(PRYOLD),
     &                 WORK(PRSCR),WORK(PRSTGS),.NOT.ININTP,
     &                 WORK(PRSTGS),WORK(PRSTGS))
         ELSE
            CALL FORMI(F,NEQN,NWANT,WORK(PRY),WORK(PRYP),WORK(PRYOLD),
     &                 WORK(PRSCR),WORK(PRSTGS),.NOT.ININTP,WRKINT,
     &                 WRKINT(STARTP))
         END IF
C
C  Set markers to show that interpolation has been initialized for
C  NWANT components.
         NWNTSV = NWANT
         ININTP = .TRUE.
      END IF
C
C  The actual evaluation of the interpolating polynomial and/or its first
C  derivative is done in EVALI.
C
      IF (METHD.EQ.1) THEN
         CALL EVALI(WORK(PRY),WORK(PRYP),WORK(PRSTGS),TWANT,REQEST,
     &              NWANT,YWANT,YPWANT)
      ELSE
         CALL EVALI(WORK(PRY),WORK(PRYP),WRKINT(STARTP),TWANT,REQEST,
     &              NWANT,YWANT,YPWANT)
      END IF
C
   20 CONTINUE
C
      CALL RKMSG(IER,SRNAME,NREC,FLAG)
C
      RETURN
      END
