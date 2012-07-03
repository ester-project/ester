      SUBROUTINE CHKFL(ASK,ERROR)
C
C  Purpose:      Enquiry routine used in conjunction with SOFTFL.
C                Reports whether a "catastrophic" error was detected.
C
C  Input:        ASK
C  Input/output: ERROR
C
C  Comments:
C  =========
C  When a "catastrophic" failure is detected, the default action of
C  RKSUITE is to write an explanation to the standard output channel,
C  OUTCH, and STOP.  SOFTFL can be used to prevent the STOP and so
C  allow the main program to continue.  It is then necessary to call
C  CHKFL with ASK = .TRUE. after every call to a user-callable routine 
C  in RKSUITE to check whether a catastrophic error occurred and take 
C  appropriate action if it did.  If there was a catastrophic error, 
C  ERROR is returned .TRUE.  Of course, you may call SETUP at any time 
C  to start a new problem, but calling any other user-callable routine 
C  in RKSUITE after a catastrophic error will lead to a STOP (even when
C  "soft failure" has been set "on").
C
C  When a catastrophic failure (IER = 911) is detected in one of
C  the routines in RKSUITE, it calls CHKFL with ASK = .FALSE. and
C  ERROR = .TRUE.  This value of ERROR is SAVEd.
C
C     .. Scalar Arguments ..
      LOGICAL           ASK, ERROR
C     .. Local Scalars ..
      LOGICAL           SAVERR
C     .. Save statement ..
      SAVE              SAVERR
C     .. Data statements ..
      DATA              SAVERR/.FALSE./
C     .. Executable Statements ..
C
      IF (ASK) THEN
         ERROR = SAVERR
      ELSE
         SAVERR = ERROR
      END IF
C
      RETURN
      END
