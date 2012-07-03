      SUBROUTINE RKMSG(IER,SRNAME,NREC,FLAG)
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
C************************************************
C
C  Purpose:      To process error messages and terminate the program
C                in the event of a "catastrophic" failure.
C
C  Input:        IER, SRNAME, NREC
C  Output:       FLAG
C
C  Common:       Initializes:    none
C                Reads:          /RKCOM7/ OUTCH
C                                /RKCOM8/ MSG, UTASK
C                                /RKCOM9/ REC
C                Alters:         none
C
C  Comments:
C  =========
C  The output variable FLAG is assigned the value of the input variable IER.
C
C  IER = -2  reports a successful call of the subroutine SRNAME and
C            indicates that special action is to be taken elsewhere
C            in the suite.  FLAG is set and a return is effected.
C
C  IER = 1   reports a successful call of the subroutine SRNAME.  FLAG
C            is set and a return is effected.
C
C  1 < IER < 911 and MSG = .TRUE.: a message of NREC records contained in
C            the array REC(*) is written to the standard output channel, 
C            OUTCH.  FLAG is set and a return is effected.
C
C  IER = 911 reports a "catastrophic" error was detected in SRNAME.  A
C            message is written to OUTCH regardless of the value of MSG and
C            normally the execution of the program is terminated.  The
C            execution is not terminated if the error is the result of an
C            indirect call to CT, RESET, or INTPO through UT (UTASK = .TRUE.).
C            Termination can be prevented by using the subroutine SOFTFL.
C
C  IER = 912 reports that a "catastrophic" error was detected earlier and
C            termination was prevented, but the user has failed to take
C            appropriate remedial action.  Execution is terminated.
C
C     .. Scalar Arguments ..
      INTEGER           FLAG, IER, NREC
      CHARACTER*(*)     SRNAME
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
      INTEGER           PLUS1
      LOGICAL           ASK, TELL
      PARAMETER         (PLUS1=1,ASK=.TRUE.,TELL=.FALSE.)
C     .. Local Scalars ..
      INTEGER           I
      LOGICAL           BADERR, OK, ON, UTCALL
C     .. External Subroutines ..
      EXTERNAL          CHKFL, RKSIT, SOFTFL
C     .. Executable Statements ..
C
C  Check where the call came from - if it is an indirect call from UT,
C  the run is not STOPped.
      UTCALL = (SRNAME.EQ.'RESET' .OR. SRNAME.EQ.'CT' .OR.
     &          SRNAME.EQ.'INTPO') .AND. UTASK
C
C  Check if can continue with integrator.
      OK = (SRNAME.EQ.'CT' .OR. SRNAME.EQ.'UT') .AND.
     &     (IER.EQ.2 .OR. IER.EQ.3 .OR. IER.EQ.4)
C
C  Check if program termination has been overridden.
      CALL SOFTFL(ASK,ON)
C
      IF ((MSG.AND.IER.GT.PLUS1) .OR. IER.GE.911) THEN
         WRITE (OUTCH,'(/A)') ' **'
         WRITE (OUTCH,'(A)') (REC(I),I=1,NREC)
         IF (IER.GE.911) THEN
            WRITE (OUTCH,'(A/A,A,A/A/)') 
     &' **',
     &' ** Catastrophic error detected in ', SRNAME, '.', 
     &' **'
            IF ((.NOT.(UTCALL.OR.ON).AND.IER.EQ.911) .OR.
     &          IER.EQ.912) THEN
               WRITE (OUTCH,'(A/A/A)') 
     &' **',
     &' ** Execution of your program is being terminated.',
     &' **'
               STOP
            END IF
         ELSE IF (OK) THEN
            WRITE (OUTCH,'(A/A,A,A,I2,A/A/A)') 
     &' **',
     &' ** Warning from routine ', SRNAME, ' with flag set ',IER, '.',
     &' ** You can continue integrating this problem.',
     &' **'
         ELSE
            WRITE (OUTCH,'(A/A,A,A,I2,A/A/A)') 
     &' **',
     &' ** Warning from routine ', SRNAME, ' with flag set ',IER, '.',
     &' ** You cannot continue integrating this problem.', 
     &' **'
         END IF
      END IF
      DO 20 I = NREC + 1, 10
         REC(I) = ' '
   20 CONTINUE
      FLAG = IER
C
C  TELL RKSIT the status of the routine associated with SRNAME
      CALL RKSIT(TELL,SRNAME,FLAG)
C
C  Indicate that a catastrophic error has been detected
      BADERR = FLAG .GE. 911
      CALL CHKFL(TELL,BADERR)
C
      RETURN
C
      END
