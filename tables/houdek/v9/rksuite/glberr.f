      SUBROUTINE GLBERR(RMSERR,ERRMAX,TERRMX,WORK)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C  If you are not familiar with the code GLBERR and how it is used in
C  conjunction with UT and CT to solve initial value problems, you should
C  study the document file rksuite.doc carefully before attempting to use 
C  the code.  The following "Brief Reminder" is intended only to remind you 
C  of the meaning, type, and size requirements of the arguments.
C
C  If ERRASS was set .TRUE. in the call to SETUP, then after any call to UT
C  or CT to advance the integration to TNOW or TWANT, the subroutine GLBERR
C  may be called to obtain an assessment of the true error of the integration.
C  At each step and for each solution component Y(L), a more accurate "true"
C  solution YT(L), an average magnitude "size(L)" of its size, and its error
C                abs(Y(L) - YT(L))/max("size(L)",THRES(L))
C  are computed.  The assessment returned in RMSERR(L) is the RMS (root-mean-
C  square) average of the error in the Lth solution component over all steps
C  of the integration from TSTART through TNOW.
C
C  OUTPUT VARIABLES
C
C     RMSERR(*) - DOUBLE PRECISION array of length NEQ
C                 RMSERR(L) approximates the RMS average of the true error 
C                 of the numerical solution for the Ith solution component,
C                 L = 1,2,...,NEQ.  The average is taken over all steps from
C                 TSTART to TNOW.
C     ERRMAX    - DOUBLE PRECISION
C                 The maximum (approximate) true error taken over all
C                 solution components and all steps from TSTART to TNOW.
C     TERRMX    - DOUBLE PRECISION
C                 First value of the independent variable where the
C                 (approximate) true error attains the maximum value ERRMAX.
C
C  WORKSPACE
C
C     WORK(*)   - DOUBLE PRECISION array as used in SETUP and UT or CT
C                 Do not alter the contents of this array.
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ERRMAX, TERRMX
C     .. Array Arguments ..
      DOUBLE PRECISION  RMSERR(*), WORK(*)
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
      PARAMETER         (SRNAME='GLBERR')
      LOGICAL           ASK
      PARAMETER         (ASK=.TRUE.)
      INTEGER           MINUS1, MINUS2
      PARAMETER         (MINUS1=-1,MINUS2=-2)
C     .. Local Scalars ..
      INTEGER           FLAG, IER, L, NREC, STATE
C     .. External Subroutines ..
      EXTERNAL          RKMSG, RKSIT
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
C
      IER = 1
      NREC = 0
C
C  Is it permissible to call GLBERR?
C
      CALL RKSIT(ASK,SRNAME,STATE)
      IF (STATE.EQ.911) THEN
         IER = 912
         NREC = 1
         WRITE (REC,'(A)')
     &' ** A catastrophic error has already been detected elsewhere.'
         GO TO 40
      END IF
      IF (STATE.EQ.MINUS2) THEN
         IER = 911
         NREC = 3
         WRITE (REC,'(A/A/A)')
     &' ** You have already made a call to GLBERR after a hard ',
     &' ** failure was reported from the integrator. You cannot',
     &' ** call GLBERR again.'
         GO TO 40
      END IF
      CALL RKSIT(ASK,'CT',STATE)
      IF (STATE.EQ.MINUS1) THEN
         IER = 911
         NREC = 1
         IF (UTASK) THEN
            WRITE (REC,'(A)')
     &' ** You have not yet called UT, so you cannot call GLBERR.'
         ELSE
            WRITE (REC,'(A)')
     &' ** You have not yet called CT, so you cannot call GLBERR.'
         END IF
         GO TO 40
      END IF
C
C  Set flag so that the routine can only be called once after a hard 
C  failure from the integrator.
      IF (STATE.EQ.5 .OR. STATE.EQ.6) IER = MINUS2
C
C  Check that ERRASS was set properly for error assessment in SETUP.
C
      IF (.NOT.ERASON) THEN
         IER = 911
         NREC = 3
         WRITE (REC,'(A/A/A)')
     &' ** No error assessment is available since you did not ',
     &' ** ask for it in your call to the routine SETUP.',
     &' ** Check your program carefully.'
         GO TO 40
      END IF
C
C  Check to see if the integrator has not actually taken a step.
C
      IF (OKSTP.EQ.0) THEN
         IER = 911
         NREC = 3
         WRITE (REC,'(A/A/A)')
     &' ** The integrator has not actually taken any successful ',
     &' ** steps.  You cannot call GLBERR in this circumstance. ',
     &' ** Check your program carefully.'
         GO TO 40
      END IF
C
C  Compute RMS error and set output variables.
C
      ERRMAX = MAXERR
      TERRMX = LOCMAX
      DO 20 L = 1, NEQN
         RMSERR(L) = SQRT(WORK(PRZERR-1+L)/OKSTP)
   20 CONTINUE
C
   40 CONTINUE
C
      CALL RKMSG(IER,SRNAME,NREC,FLAG)
C
      RETURN
      END
