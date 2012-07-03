      SUBROUTINE STIFFC(ALPHA,BETA,R1,R2)
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
C************************************************
C
C  Input:  ALPHA, BETA
C  Output: R1(*), R2(*)
C
C  This subroutine computes the two complex roots R1 and R2 of
C  the quadratic equation X**2 + ALPHA*X + BETA = 0.  The magnitude
C  of R1 is greater than or equal to the magnitude of R2. R1 and R2 are
C  returned as vectors of two components with the first being the real
C  part of the complex number and the second being the imaginary part.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA, BETA
C     .. Array Arguments ..
      DOUBLE PRECISION  R1(2), R2(2)
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, TWO
      PARAMETER         (ZERO=0.0D+0,TWO=2.0D+0)
C     .. Local Scalars ..
      DOUBLE PRECISION  DISC, SQDISC, TEMP
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      TEMP = ALPHA/TWO
      DISC = TEMP**2 - BETA
      IF (DISC.EQ.ZERO) THEN
C
C  Double root.
C
         R1(1) = -TEMP
         R1(2) = ZERO
         R2(1) = R1(1)
         R2(2) = R1(2)
         RETURN
      END IF
C
      SQDISC = SQRT(ABS(DISC))
      IF (DISC.LT.ZERO) THEN
C
C  Complex conjugate roots.
C
         R1(1) = -TEMP
         R1(2) = SQDISC
         R2(1) = R1(1)
         R2(2) = -R1(2)
      ELSE
C
C  Real pair of roots.  Calculate the bigger one in R1(1).
C
         IF (TEMP.GT.ZERO) THEN
            R1(1) = -TEMP - SQDISC
         ELSE
            R1(1) = -TEMP + SQDISC
         END IF
         R1(2) = ZERO
         R2(1) = BETA/R1(1)
         R2(2) = ZERO
      END IF
C
      RETURN
      END
