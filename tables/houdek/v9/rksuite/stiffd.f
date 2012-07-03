      SUBROUTINE STIFFD(V,HAVG,X,Y,F,FXY,WT,SCALE,VDOTV,Z,ZDOTZ,VTEMP)
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
C************************************************
C
C  External:     F
C  Input:        V(*), HAVG, X, Y(*), FXY(*), WT(*), SCALE, VDOTV,
C  Output:       Z(*), ZDOTZ
C  Workspace:    VTEMP(*)
C
C  For an input vector V(*) of length NEQ, this subroutine computes a vector
C  Z(*) that approximates the product HAVG*J*V where HAVG is an input scalar 
C  and J is the Jacobian matrix of a function F evaluated at the input 
C  arguments (X,Y(*)).  This function is defined by a subroutine of the form
C  F(T,U,F) that when given T and U(*), returns the value of the function in 
C  F(*).  The input vector FXY(*) is defined by F(X,Y,FXY).  Scaling is a 
C  delicate matter.  A weighted Euclidean norm is used with the (positive) 
C  weights provided in WT(*).  The input scalar SCALE is the square root of 
C  the unit roundoff times the norm of Y(*).  The square of the norm of the
C  input vector V(*) is input as VDOTV.  The routine outputs the square of
C  the norm of the output vector Z(*) as ZDOTZ.  The subroutine calls the
C  DOUBLE PRECISION FUNCTION DOTPRD(U,V,WT,NEQ) to compute the dot (inner)
C  product.  The vector VTEMP(*) is used for working storage.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HAVG, SCALE, VDOTV, X, ZDOTZ
C     .. Array Arguments ..
      DOUBLE PRECISION  FXY(*), V(*), VTEMP(*), WT(*), Y(*), Z(*)
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
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP1, TEMP2
      INTEGER           L
C     .. External Functions ..
      DOUBLE PRECISION  DOTPRD
      EXTERNAL          DOTPRD
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
C
C  Scale V(*) so that it can be used as an increment to Y(*)
C  for an accurate difference approximation to the Jacobian.
C
      TEMP1 = SCALE/SQRT(VDOTV)
      DO 20 L = 1, NEQN
         VTEMP(L) = Y(L) + TEMP1*V(L)
   20 CONTINUE
C
      CALL F(X,VTEMP,Z)
      NFCN = NFCN + 1
C
C  Form the difference approximation.  At the same time undo
C  the scaling of V(*) and introduce the factor of HAVG.
C
      TEMP2 = HAVG/TEMP1
      DO 40 L = 1, NEQN
         Z(L) = TEMP2*(Z(L)-FXY(L))
   40 CONTINUE
C
      ZDOTZ = DOTPRD(Z,Z,WT,NEQN)
C
      RETURN
      END
