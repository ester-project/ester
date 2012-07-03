      SUBROUTINE STEPA(TNOW,Y,YP,TSTG,YSTG,YPSTG,HTRY,WEIGHT,CUTBAK)
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
C************************************************
C
C  Purpose:      To calculate an "on-scale" step size for phase 2 of
C                the initial step size computation.
C
C  Input:        TNOW, Y(*), YP(*), TSTG, YSTG(*), YPSTG(*)
C  Input/output: HTRY, WEIGHT
C  Output:       CUTBAK
C
C  Common:       Initializes:    none
C                Reads:          /RKCOM1/ TND, NEQ
C                                /RKCOM5/ STBRAD, RS1, RS4
C                                /RKCOM7/ RNDOFF
C                Alters:         none
C
C  Comments:
C  =========
C  This subroutine is used during the first three stages of the first step.
C  A Lipschitz constant L for the differential equation in autonomous form
C  is approximated, and the product abs(HTRY)*L is compared to an approximate
C  radius, STBRAD, of the stability region of the method. The step size is 
C  reduced as necessary, within a range specified by the step size control 
C  parameters RS1 and RS4, to assure stability and give some confidence in 
C  the error estimator.  If HTRY is reduced, CUTBAK is set .TRUE..
C
C  Y(*) and YP(*) contain the solution and its derivative at TNOW and
C  similarly YSTG(*) and YPSTG(*) contain approximations at TSTG.
C
C  Normally the weights used in the control of the error depend on the
C  size of the solution at the beginning and at the end of the step, but
C  at this time we do not have a solution at the end of the step.  Each
C  stage YSTG(*) of the Runge - Kutta process represents a low order
C  approximation to the solution at TSTG.  Because the initial value of
C  WEIGHT(*) provided in the first phase of the scheme is based only on
C  the solution at T and THRES(*), it is continually updated in STEPA to
C  account for the size of the solution throughout the step as revealed
C  by the intermediate stages YSTG(*). Inside this subroutine only, the
C  differential equation is converted to autonomous form. After the
C  conversion, the end of the interval of integration, TND, is used
C  to define a suitable weight for the independent variable.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HTRY, TNOW, TSTG
      LOGICAL           CUTBAK
C     .. Array Arguments ..
      DOUBLE PRECISION  WEIGHT(*), Y(*), YP(*), YPSTG(*), YSTG(*)
C     .. Common Block for Problem Definition ..
      DOUBLE PRECISION  TSTRT, TND, DIR, HSTRT, TOLR
      INTEGER           NEQN
      COMMON /RKCOM1/   TSTRT, TND, DIR, HSTRT, TOLR, NEQN
      SAVE   /RKCOM1/
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
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Local Scalars ..
      DOUBLE PRECISION  ARGDIF, FDIFF, SCL, TDIFF, TWT, WT, YNRM, YSTGNM
      INTEGER           L
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Executable Statements ..
C
C  Update the weights to account for the current intermediate solution
C  approximation YSTG(*).  Compute the sizes of Y(*) and YSTG(*) in the
C  new norm.  The size of the Lipschitz constant is assessed by a difference
C  in the arguments Y(*), YSTG(*) and a difference in the function evaluated
C  at these arguments.
C
      YNRM = ZERO
      YSTGNM = ZERO
      ARGDIF = ZERO
      FDIFF = ZERO
      DO 20 L = 1, NEQN
         WT = MAX(WEIGHT(L),ABS(YSTG(L)))
         WEIGHT(L) = WT
         YNRM = MAX(YNRM,ABS(Y(L))/WT)
         YSTGNM = MAX(YSTGNM,ABS(YSTG(L))/WT)
         ARGDIF = MAX(ARGDIF,ABS(YSTG(L)-Y(L))/WT)
         FDIFF = MAX(FDIFF,ABS(YPSTG(L)-YP(L))/WT)
   20 CONTINUE
C
C  The transformation of the equation to autonomous form is done
C  implicitly.  The difference of the arguments must take into account
C  the difference between the values of the independent variable T and
C  TSTG. The difference of the corresponding component of the function
C  is zero because of the way the standard transformation is done.
C
      TDIFF = TSTG - TNOW
      TWT = ABS(TND-TNOW)
      YNRM = MAX(YNRM,ABS(TNOW)/TWT)
      YSTGNM = MAX(YSTGNM,ABS(TSTG)/TWT)
      ARGDIF = MAX(ARGDIF,ABS(TDIFF)/TWT)
C
C  The ratio FDIFF/ARGDIF is a lower bound for, and an approximation to, a
C  Lipschitz constant L for the differential equation written in autonomous
C  form.  First we must ask if the difference ARGDIF is significant in the 
C  precision available.  If it appears to be, we insist that abs(HTRY)*L be 
C  less than an approximate radius, STBRAD, of the stability region of the
C  method.  This is more stringent than necessary for stability, possibly a
C  lot more stringent, but the aim is to get an HTRY small enough that the
C  error estimate for the step is credible.  The reduction is required to be
C  at least as much as the step control parameter RS1. It is necessary to 
C  limit the reduction of HTRY at any one time because we may be misled in 
C  the size of the reduction that is appropriate due to nonlinearity of the 
C  differential equation and to inaccurate weights caused by HTRY much too 
C  large.  The reduction is not permitted to be more than the step control 
C  parameter RS4.
C
      CUTBAK = .FALSE.
      IF (ARGDIF.GT.RNDOFF*MAX(YNRM,YSTGNM)) THEN
         IF ((ABS(HTRY)*FDIFF).GT.(STBRAD*ARGDIF)) THEN
            SCL = (STBRAD*ARGDIF)/(ABS(HTRY)*FDIFF)
            SCL = MIN(SCL,RS1)
            SCL = MAX(SCL,RS4)
            HTRY = SCL*HTRY
            CUTBAK = .TRUE.
         END IF
      END IF
C
      RETURN
      END
