      SUBROUTINE STIFFB(V1V1,V0V1,V0V0,ROLD,RHO,ROOT1,ROOT2,ROOTRE)
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
C************************************************
C
C  Input:        V1V1, V0V1, V0V0
C  Input/output: ROLD
C  Output:       RHO, ROOT1(*),ROOT2(*),ROOTRE
C
C  Decide if the iteration has degenerated because of a strongly
C  dominant real eigenvalue.  Have just computed the latest iterate.
C  V1V1 is its dot product with itself, V0V1 is the dot product
C  of the previous iterate with the current one, and V0V0 is the
C  dot product of the previous iterate with itself.  ROLD is a
C  previous Rayleigh quotient approximating a dominant real
C  eigenvalue.  It must be computed directly the first time the
C  subroutine is called.  It is updated each call to STIFFB, hence
C  is available for subsequent calls.
C
C  If there is a strongly dominant real eigenvalue, ROOTRE is set
C  .TRUE., ROOT1(*) returns the eigenvalue, RHO returns the magnitude
C  of the eigenvalue, and ROOT2(*) is set to zero.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RHO, ROLD, V0V0, V0V1, V1V1
      LOGICAL           ROOTRE
C     .. Array Arguments ..
      DOUBLE PRECISION  ROOT1(2), ROOT2(2)
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, P001
      PARAMETER         (ZERO=0.0D+0,P001=0.001D+0)
C     .. Local Scalars ..
      DOUBLE PRECISION  DET, R, RES
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
C
      R = V0V1/V0V0
      RHO = ABS(R)
      DET = V0V0*V1V1 - V0V1**2
      RES = ABS(DET/V0V0)
      ROOTRE = DET .EQ. ZERO .OR. (RES.LE.V1V1*P001**2 .AND.
     &         ABS(R-ROLD).LE.P001*RHO)
      IF (ROOTRE) THEN
         ROOT1(1) = R
         ROOT1(2) = ZERO
         ROOT2(1) = ZERO
         ROOT2(2) = ZERO
      END IF
      ROLD = R
C
      RETURN
      END
