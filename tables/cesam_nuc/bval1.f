
c******************************************************************

	SUBROUTINE bval1(x,y,m,l,q,d)
	
c	subroutine public de mod_numerique

c	calcul les B-splines normalisees d'ordre m > 1
c	et leurs derivees premieres au point x tel que y(l) <= x < y(l+1)
c	si x est l'extremite droite : y(l) < x = y(l+1) on obtient
c	la limite a droite
c	les valeurs des m b-splines non nulles : l-m+1, ..., l sont 
c	dans q(1), ..., q(m), les derivees dans d(1), ... , d(m)
c	!! ATTENTION : q doit etre dimensionne au moins a m+1 !!

c	d'apres l'algorithme 5-5 p.192 de schumaker 
c	et addition d'une instruction pour les derivees

c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A.

c	derive de bvald.f

c entrees:
c	x: point de calcul
c	y(knot): points de table
c	m: ordre des splines
c	l: indice

c sorties:
c	q(m+1): table des valeurs des m B-splines
c	d(m): table des valeurs des derivees premieres

c------------------------------------------------------------------------

	USE mod_kind
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) , DIMENSION(:) :: y
	REAL (kind=dp), INTENT(in) :: x
	INTEGER, INTENT(in) :: l, m
	REAL (kind=dp), INTENT(out) , DIMENSION(:) :: d, q

	REAL (kind=dp) :: a1, a2, denom
	
	INTEGER :: i, j

c-------------------------------------------------------------------------

	DO j=1,m-1
	 q(j)=0.d0
	ENDDO
	
	q(m)=1.d0/(y(l+1)-y(l)) ; q(m+1)=0.d0

	DO j=2,m-1
	 DO i=m-j+1,m
	  denom=y(i+l-m+j)-y(i+l-m) ; a1=(x-y(i+l-m))/denom
	  a2=1.d0-a1 ; q(i)=a1*q(i)+a2*q(i+1)
	 ENDDO
	ENDDO

	DO i=1,m
	 d(i)=(q(i)-q(i+1))*(m-1.d0) ; q(i)=(x-y(i+l-m))*q(i)+(y(i+l)-x)*q(i+1)
	ENDDO

	RETURN

	END SUBROUTINE bval1
