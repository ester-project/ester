      SUBROUTINE EVALI(Y,YP,P,TWANT,REQEST,NWANT,YWANT,YPWANT)
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
*************************************************
C
C  Purpose:    Evaluation of an interpolating polynomial and/or its
C              first derivative.
C
C  Input:      Y(*), YP(*), P(NWANT,*), TWANT, REQEST, NWANT
C  Output:     YWANT(*), YPWANT(*)
C
C  Common:     Initializes:    none
C              Reads:          /RKCOM2/ HOLD, T
C                              /RKCOM4/ MINTP
C              Alters:         none
C
C  Comments:
C  =========
C  The interpolant is evaluated at TWANT to approximate the solution,
C  YWANT, and/or its first derivative there, YPWANT. Only the first
C  NWANT components of the answer are computed. There are three cases
C  that are indicated by the first character of REQEST:
C    REQEST(1:1) = `S' or `s'- compute approximate `S'olution only.
C                = `D' or `d'- compute approximate first `D'erivative
C                              of the solution only.
C                = `B' or `b'- compute `B'oth approximate solution and
C                              first derivative.
C  The coefficents of the polynomial are contained in Y(*), YP(*) and
C  P(NWANT,*).
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TWANT
      INTEGER           NWANT
      CHARACTER*(*)     REQEST
C     .. Array Arguments ..
      DOUBLE PRECISION  P(NWANT,*), Y(*), YP(*), YPWANT(*), YWANT(*)
C     .. Common Block to hold Problem Status ..
      DOUBLE PRECISION  T, H, TOLD, HOLD
      INTEGER           NFCN, SVNFCN, OKSTP, FLSTP
      LOGICAL           FIRST, LAST
      COMMON /RKCOM2/   T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP, FLSTP,
     &                  FIRST, LAST
      SAVE   /RKCOM2/
C     .. Common Block to hold Formula Definitions ..
      DOUBLE PRECISION  A(13,13), B(13), C(13), BHAT(13), R(11,6),
     &                  E(7)
      INTEGER           PTR(13), NSTAGE, METHD, MINTP
      LOGICAL           INTP
      COMMON /RKCOM4/   A, B, C, BHAT, R, E, PTR, NSTAGE, METHD,
     &                  MINTP, INTP
      SAVE   /RKCOM4/
C     .. Local Scalars ..
      DOUBLE PRECISION  SIGMA
      INTEGER           K, L
      CHARACTER         REQST1
C     .. Executable Statements ..
C
C  Evaluate the interpolating polynomial of degree MINTP in terms of the
C  shifted and scaled independent variable SIGMA.
C
      SIGMA = (TWANT-T)/HOLD
C
      REQST1 = REQEST(1:1)
      IF (REQST1.EQ.'S' .OR. REQST1.EQ.'s' .OR. 
     &    REQST1.EQ.'B' .OR. REQST1.EQ.'b') THEN
C
         DO 20 L = 1, NWANT
            YWANT(L) = P(L,MINTP-1)*SIGMA
   20    CONTINUE
         DO 60 K = MINTP - 2, 1, -1
            DO 40 L = 1, NWANT
               YWANT(L) = (YWANT(L)+P(L,K))*SIGMA
   40       CONTINUE
   60    CONTINUE
         DO 80 L = 1, NWANT
            YWANT(L) = (YWANT(L)+HOLD*YP(L))*SIGMA + Y(L)
   80    CONTINUE
      END IF
C
C  Evaluate the derivative of the interpolating polynomial.
C
      IF (REQST1.EQ.'D' .OR. REQST1.EQ.'d' .OR. 
     &    REQST1.EQ.'B' .OR. REQST1.EQ.'b') THEN
C
C  The derivative of the interpolating polynomial with respect to TWANT 
C  is the derivative with respect to S divided by HOLD.
C
         DO 100 L = 1, NWANT
            YPWANT(L) = MINTP*P(L,MINTP-1)*SIGMA
  100    CONTINUE
         DO 140 K = MINTP - 1, 2, -1
            DO 120 L = 1, NWANT
               YPWANT(L) = (YPWANT(L)+K*P(L,K-1))*SIGMA
  120       CONTINUE
  140    CONTINUE
         DO 160 L = 1, NWANT
            YPWANT(L) = (YPWANT(L)+HOLD*YP(L))/HOLD
  160    CONTINUE
      END IF
C
      RETURN
      END
