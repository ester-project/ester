
c********************************************************************

	SUBROUTINE horner(k,der,z,a,x,p)

c	routine private du module mod_numerique

c	algorithme de horner

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	version: 03 06 91
c	F95 29 08 00
c	CESAM2k

c entrées
c	k : degré du polynôme
c	der : ordre de dérivation
c	z : point du calcul

c entrées/sorties
c	a : coefficients / nouveaux coefficients
c	x : points de centrage / nouveaux points

c sorties
c	p : polynôme et dérivées

c---------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: z
	INTEGER, INTENT(in) :: k, der		
	REAL (kind=dp), INTENT(inout), DIMENSION(0:) :: a, x
	REAL (kind=dp), INTENT(out), DIMENSION(0:) :: p	

	REAL (kind=dp) :: fac
	INTEGER :: i, j

c--------------------------------------------------------------------

2000	FORMAT(8es10.3)

c	PRINT*,'k,der,z/a/x/p',k,der,z ; WRITE(*,2000)a(0:k)
c	WRITE(*,2000)x(0:k)

	fac=1.d0
	DO j=0,der
	 IF(j /= 0)THEN
	  DO i=k-1,0,-1
	   x(i+1)=x(i)
	  ENDDO
	  x(0)=z ; fac=fac*j
	 ENDIF
	 DO i=k-1,0,-1
	  a(i)=a(i+1)*(z-x(i))+a(i)
	 ENDDO
	 p(j)=fac*a(j)
	ENDDO
c	WRITE(*,2000)p(0:k)

	RETURN

	END SUBROUTINE horner
