      SUBROUTINE STAT(TOTFCN,STPCST,WASTE,STPSOK,HNEXT)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C  If you are not familiar with the code STAT and how it is used in
C  conjunction with the integrators CT and UT, you should study the 
C  document file rksuite.doc carefully before attempting to use the code. 
C  The following "Brief Reminder" is intended only to remind you of the 
C  meaning, type, and size requirements of the arguments.
C
C  STAT is called to obtain some details about the integration.
C
C  OUTPUT VARIABLES
C
C     TOTFCN    - INTEGER
C                 Total number of calls to F in the integration so far --
C                 a measure of the cost of the integration.
C     STPCST    - INTEGER
C                 Cost of a typical step with this METHOD measured in
C                 calls to F.
C     WASTE     - DOUBLE PRECISION
C                 The number of attempted steps that failed to meet the 
C                 local error requirement divided by the total number of 
C                 steps attempted so far in the integration.
C     STPSOK    - INTEGER
C                 The number of accepted steps.
C     HNEXT     - DOUBLE PRECISION
C                 The step size the integrator plans to use for the next step.
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HNEXT, WASTE
      INTEGER           STPCST, STPSOK, TOTFCN
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
      PARAMETER         (SRNAME='STAT')
      LOGICAL           ASK
      INTEGER           MINUS1, MINUS2
      PARAMETER         (ASK=.TRUE.,MINUS1=-1,MINUS2=-2)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Local Scalars ..
      INTEGER           FLAG, IER, NREC, STATE
C     .. External Subroutines ..
      EXTERNAL          RKMSG, RKSIT
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
      IER = 1
      NREC = 0
C
C  Is it permissible to call STAT?
C
      CALL RKSIT(ASK,SRNAME,STATE)
      IF (STATE.EQ.911) THEN
         IER = 912
         NREC = 1
         WRITE (REC,'(A)')
     &' ** A catastrophic error has already been detected elsewhere.'
         GO TO 20
      END IF
      IF (STATE.EQ.MINUS2) THEN
         IER = 911
         NREC = 3
         WRITE (REC,'(A/A/A)')
     &' ** You have already made a call to STAT after a hard   ',
     &' ** failure was reported from the integrator. You cannot',
     &' ** call STAT again.'
         GO TO 20
      END IF
      CALL RKSIT(ASK,'CT',STATE)
      IF (STATE.EQ.MINUS1) THEN
         IER = 911
         NREC = 1
         IF (UTASK) THEN
            WRITE (REC,'(A)')
     &' ** You have not called UT, so you cannot use STAT.'
         ELSE
            WRITE (REC,'(A)')
     &' ** You have not called CT, so you cannot use STAT.'
         END IF
         GO TO 20
      END IF
C
C  Set flag so that the routine can only be called once after a hard 
C  failure from the integrator.
      IF (STATE.EQ.5 .OR. STATE.EQ.6) IER = MINUS2
C
      TOTFCN = SVNFCN + NFCN
      IF (ERASON) TOTFCN = TOTFCN + GNFCN
      STPCST = COST
      STPSOK = OKSTP
      IF (OKSTP.LE.1) THEN
         WASTE = ZERO
      ELSE
         WASTE = DBLE(FLSTP)/DBLE(FLSTP+OKSTP)
      END IF
      HNEXT = H
C
   20 CONTINUE
C
      CALL RKMSG(IER,SRNAME,NREC,FLAG)
C
      RETURN
      END
