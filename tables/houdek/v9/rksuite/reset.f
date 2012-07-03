      SUBROUTINE RESET(TENDNU)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C  If you are not familiar with the code RESET and how it is used in
C  conjunction with CT to solve initial value problems, you should study the
C  document file rksuite.doc carefully before attempting to use the code. The
C  following "Brief Reminder" is intended only to remind you of the meaning,
C  type, and size requirements of the arguments.
C
C  The integration using CT proceeds from TSTART in the direction of TEND, and
C  is now at TNOW.  To reset TEND to a new value TENDNU, just call RESET with
C  TENDNU as the argument.  You must continue integrating in the same
C  direction, so the sign of (TENDNU - TNOW) must be the same as the sign of
C  (TEND - TSTART). To change direction you must restart by a call to SETUP.
C
C  INPUT VARIABLE
C
C     TENDNU    - DOUBLE PRECISION
C                 The new value of TEND.
C                 Constraint: TENDNU and TNOW must be clearly distinguishable
C                 in the precision used.  The sign of (TENDNU - TNOW) must be
C                 the same as the sign of (TEND - TSTART).
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TENDNU
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
      PARAMETER         (SRNAME='RESET')
      LOGICAL           ASK
      INTEGER           MINUS1, MINUS2
      PARAMETER         (ASK=.TRUE.,MINUS1=-1,MINUS2=-2)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Local Scalars ..
      DOUBLE PRECISION  HMIN, TDIFF
      INTEGER           FLAG, IER, NREC, STATE, STATE1
C     .. External Subroutines ..
      EXTERNAL          RKMSG, RKSIT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
      IER = 1
      NREC = 0
C
C  Is it permissible to call RESET?
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
     &' ** You have called RESET after you specified to SETUP that ',
     &' ** you were going to use UT. This is not permitted.'
            GO TO 20
         END IF
      END IF
      IF (STATE.EQ.MINUS1) THEN
         IER = 911
         NREC = 1
         WRITE (REC,'(A)')
     &' ** You have not called CT, so you cannot use RESET.'
         GO TO 20
      END IF
      IF (STATE.EQ.5 .OR. STATE.EQ.6) THEN
         IER = 911
         NREC = 2
         WRITE (REC,'(A,I1,A/A)')
     &' ** CT has returned with CFLAG =  ',STATE,'.',
     &' ** You cannot call RESET in this circumstance.'
         GO TO 20
      END IF
C
C  Check value of TENDNU
C
      IF (DIR.GT.ZERO .AND. TENDNU.LE.T) THEN
         IER = 911
         NREC = 4
         WRITE (REC,'(A/A,D13.5/A,D13.5,A/A)')
     &' ** Integration is proceeding in the positive direction. The ',
     &' ** current value for the independent variable is ',T,
     &' ** and you have set TENDNU = ',TENDNU,'.  TENDNU must be ',
     &' ** greater than T.'
      ELSE IF (DIR.LT.ZERO .AND. TENDNU.GE.T) THEN
         IER = 911
         NREC = 4
         WRITE (REC,'(A/A,D13.5/A,D13.5,A/A)')
     &' ** Integration is proceeding in the negative direction. The ',
     &' ** current value for the independent variable is ',T,
     &' ** and you have set TENDNU = ',TENDNU,'.  TENDNU must be ',
     &' ** less than T.'
      ELSE
         HMIN = MAX(TINY,TOOSML*MAX(ABS(T),ABS(TENDNU)))
         TDIFF = ABS(TENDNU-T)
         IF (TDIFF.LT.HMIN) THEN
            IER = 911
            NREC = 4
            WRITE (REC,'(A,D13.5,A/A,D13.5,A/A/A,D13.5,A)')
     &' ** The current value of the independent variable T is ',T,'.',
     &' ** The TENDNU you supplied has ABS(TENDNU-T) = ',TDIFF,'.',
     &' ** For the METHOD and the precision of the computer being ',
     &' ** used, this difference must be at least ',HMIN,'.'
         END IF
      END IF
      IF (IER.EQ.911) GO TO 20
C
      TND = TENDNU
      LAST = .FALSE.
C
   20 CONTINUE
C
      CALL RKMSG(IER,SRNAME,NREC,FLAG)
C
      RETURN
      END
