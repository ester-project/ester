      SUBROUTINE SOFTFL(ASK,ON)
C
C  Purpose:      To prevent a program STOP after a "catastrophic"
C                failure when using a routine from RKSUITE.
C
C  Input:        ASK
C  Input/output: ON
C
C  Comments:
C  =========
C  When a "catastrophic" failure is detected, the default action of
C  RKSUITE is to write an explanation to the standard output channel,
C  OUTCH, and STOP.  This subroutine can be used to prevent the STOP and
C  so allow the main program to continue.  To do this, you call SOFTFL with
C  ASK = .FALSE. and ON = .TRUE.  You must then call the subroutine CHKFL
C  after every call to a user-callable routine in RKSUITE to check whether
C  a catastrophic error occurred and take appropriate action if it did.  Of
C  course, you may call SETUP at any time to start a new problem, but calling
C  any other user-callable routine in RKSUITE after a catastrophic error will
C  lead to a STOP (even when "soft failure" has been set "on").
C
C  When ON is set by a call to SOFTFL with ASK = .FALSE., the value of ON
C  is SAVEd.  The subroutine RKMSG in RKSUITE calls SOFTFL with ASK = .TRUE.
C  to find out the SAVEd value of ON.
C
C     .. Scalar Arguments ..
      LOGICAL           ASK, ON
C     .. Local Scalars ..
      LOGICAL           SOFT
C     .. Save statement ..
      SAVE              SOFT
C     .. Data statements ..
      DATA              SOFT/.FALSE./
C     .. Executable Statements ..
C
      IF (ASK) THEN
         ON = SOFT
      ELSE
         SOFT = ON
      END IF
C
      RETURN
      END
