      SUBROUTINE STIFF(F,HAVG,JFLSTP,TOOMCH,MAXFCN,WORK,IER,NREC)
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
C************************************************
C
C  Purpose:      Diagnose stiffness.  This depends on two things: whether
C                the step size is being restricted on grounds of stability
C                and whether the integration to TND can be completed in no
C                more than MAXFCN function evaluations.
C
C  Input:        HAVG, TOOMCH, MAXFCN, WORK(*)
C  Input/output: JFLSTP
C  Output:       IER, NREC
C  Workspace:    WORK(*)
C  External:     F
C
C  Common:       Initializes:    /RKCOM9/ REC
C                Reads:          /RKCOM1/ TND, NEQN
C                                /RKCOM2/ T, H, NFCN, SVNFCN, OKSTP
C                                /RKCOM3/ PRY, PRYP, PRTHRS, PRWT, PRSCR,
C                                         PRSTGS, PRYOLD
C                                /RKCOM5/ COST
C                Alters:         /RKCOM2/ NFCN
C                                /RKCOM9/ REC
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HAVG
      INTEGER           IER, JFLSTP, MAXFCN, NREC
      LOGICAL           TOOMCH
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(*)
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
C     .. Common Block for Error Message ..
      CHARACTER*80      REC(10)
      COMMON /RKCOM9/   REC
      SAVE   /RKCOM9/
C     .. Parameters ..
      DOUBLE PRECISION  HALF
      PARAMETER         (HALF=0.5D0)
C     .. Local Scalars ..
      DOUBLE PRECISION  AVGY, XTRAWK
      INTEGER           L
      LOGICAL           LOTSFL, STIF, UNSURE
C     .. External Subroutines ..
      EXTERNAL          STIFFA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, MOD
C     .. Executable Statements ..
C
      IF (MOD(OKSTP-10,40).EQ.0) THEN
         LOTSFL = JFLSTP .GE. 10
         JFLSTP = 0
      ELSE
         LOTSFL = .FALSE.
      END IF
C
C  If either too much work has been done or there are lots of failed steps,
C  test for stiffness.
C
      IF (TOOMCH .OR. LOTSFL) THEN
C
C  Regenerate weight vector
         DO 20 L = 1, NEQN
            AVGY = HALF*(ABS(WORK(PRY-1+L))+ABS(WORK(PRYOLD-1+L)))
            WORK(PRWT-1+L) = MAX(AVGY,WORK(PRTHRS-1+L))
   20    CONTINUE
C
C  STIFFA determines whether the problem is STIFF. In some circumstances it
C  is UNSURE.  The decision depends on two things: whether the step size is
C  being restricted on grounds of stability and whether the integration to
C  TND can be completed in no more than MAXFCN function evaluations.  The
C  last four arguments of STIFFA are vectors of length NEQN used for working
C  storage.  Some storage in WORK(*) reserved for the stages (there are a
C  minimum of three such vectors reserved for the METHDs implemented) and
C  the scratch vector starting at PRSCR are used for this purpose.
C
         CALL STIFFA(F,T,WORK(PRY),H,HAVG,TND,MAXFCN,WORK(PRWT),
     &               WORK(PRYP),WORK(PRERST),UNSURE,STIF,WORK(PRSTGS),
     &               WORK(PRSTGS+NEQN),WORK(PRSTGS+2*NEQN),WORK(PRSCR))
         IF (.NOT.UNSURE) THEN
            IF (STIF) THEN
C
C  Predict how much eXTRA WorK will be needed to reach TND.
               XTRAWK = (COST*ABS((TND-T)/HAVG))/DBLE(SVNFCN+NFCN)
               IER = 4
               WRITE (REC(NREC+1),'(A)')
     &' ** Your problem has been diagnosed as stiff.  If the '
               WRITE (REC(NREC+2),'(A,D13.5)')
     &' ** situation persists, it will cost roughly ', XTRAWK
               WRITE (REC(NREC+3),'(A)')
     &' ** times as much to reach TEND as it has cost to reach TNOW.'
               WRITE (REC(NREC+4),'(A)')
     &' ** You should probably change to a code intended for '
               WRITE (REC(NREC+5),'(A)')
     &' ** stiff problems. '
               NREC = NREC + 5
            END IF
         END IF
      END IF
C
      RETURN
      END
