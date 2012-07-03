      SUBROUTINE MCONST(METHOD)
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
C************************************************
C
C  Purpose:   Sets machine-dependent global quantities
C
C  Common:    Initializes:    /RKCOM7/ OUTCH, MCHEPS, DWARF, RNDOFF,
C                                      SQRRMC, CUBRMC, TINY
C             Reads:          none
C             Alters:         none
C
C  Comments:
C  =========
C  OUTCH, MCHEPS, DWARF are pure environmental parameters with values
C  obtained from a call to ENVIRN. The other quantities set depend on
C  the environmental parameters, the implementation, and, possibly,
C  METHOD. At present the METHODs implemented in the RK suite do not
C  influence the values of these quantities.
C  OUTCH  - Standard output channel
C  MCHEPS - Largest positive number such that 1.0D0 + MCHEPS = 1.0D0
C  DWARF  - Smallest positive number
C  RNDOFF - 10 times MCHEPS
C  SQRRMC - Square root of MCHEPS
C  CUBRMC - Cube root of MCHEPS
C  TINY   - Square root of DWARF
C
C     .. Scalar Arguments ..
      INTEGER           METHOD
C     .. Common Block for Environment Parameters ..
      DOUBLE PRECISION  MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY
      INTEGER           OUTCH
      COMMON /RKCOM7/   MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY,
     &                  OUTCH
      SAVE   /RKCOM7/
C     .. Parameters ..
      DOUBLE PRECISION  TEN, THIRD
      PARAMETER         (TEN=10.0D+0,THIRD=1.0D+0/3.0D+0)
C     .. External Subroutines ..
      EXTERNAL          ENVIRN
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
C
      CALL ENVIRN(OUTCH,MCHEPS,DWARF)
C
      RNDOFF = TEN*MCHEPS
      SQRRMC = SQRT(MCHEPS)
      CUBRMC = MCHEPS**THIRD
      TINY = SQRT(DWARF)
C
      RETURN
      END
