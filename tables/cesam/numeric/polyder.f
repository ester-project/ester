	
c************************************************************************
			
	SUBROUTINE polyder(a,n,der,x,p)
	
c	subroutine public du module mod_numerique

c	calcul, au point x, des dérivées jusqu'à l'ordre der du polynôme
c	a_0 + a_1 x + ... +a_n x**n algorithme de horner

c	le tableau a est conservé

c	dérivées dans p(0), p(1), ... , p(der)

c	Auteur: P. Morel, Département J.D. Cassini, O.C.A.
c	version: 03 06 91
c	CESAM2k

c---------------------------------------------------------------------

	USE mod_kind
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION(0:) :: a
	REAL (kind=dp), INTENT(in) :: x	
	INTEGER, INTENT(in) :: der, n
		
	REAL (kind=dp), INTENT(out), DIMENSION(0:) :: p 

	REAL (kind=dp) :: fac	

	INTEGER :: i, j
	
c-------------------------------------------------------------

	p(0:n)=a(0:n)

	DO j=0,der
	 DO i=n-1,j,-1
	  p(i)=p(i+1)*x+p(i)
	 ENDDO	!i
	ENDDO	!j

	IF(der > 1)THEN
	 fac=1.d0
	 DO i=2,der
	  fac=fac*i ; p(i)=p(i)*fac
	 ENDDO	!i
	ENDIF

	RETURN

	END SUBROUTINE polyder
	
