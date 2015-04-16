
c**********************************************************************

	SUBROUTINE bval0(x,y,m,l,q)

c	subroutine private du module mod_numerique

c	calcul les B-splines normalisees d'ordre m > 1
c	au point x tel que y(l) .le. x. lt. y(l+1)
c	si x est l'extremite droite : y(l) .lt. x .eq. y(l+1) on obtient
c	la limite a droite
c	les valeurs des m b-splines non nulles : l-m+1, ..., l sont 
c	dans q(1), ..., q(m)
c	!! ATTENTION : q doit etre dimensionné au moins a m+1 !!

c	d'apres l'algorithme 5-5 p.192 de schumaker

c	derive sbval.f

c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice

c entrees:
c	x: point de calcul
c	y(knot): points de table
c	m: ordre des splines
c	l: indice

c sorties:
c	q(m+1): table des valeurs des m B-splines

c-----------------------------------------------------------------	

	USE mod_kind
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION (:) :: y
	REAL (kind=dp), INTENT(in) :: x	
	INTEGER, INTENT(in) :: l, m	
	REAL (kind=dp), INTENT(out), DIMENSION (:) :: q	

	REAL (kind=dp) :: a1, a2, denom
	
	INTEGER :: i, j
		
c-----------------------------------------------------------------		
	
	q=0.d0 ; q(m)=1.d0/(y(l+1)-y(l))

	DO j=2,m-1
	 DO i=m-j+1,m
	  denom=y(i+l-m+j)-y(i+l-m) ; a1=(x-y(i+l-m))/denom
	  a2=1.d0-a1 ; q(i)=a1*q(i)+a2*q(i+1)
	 ENDDO
	ENDDO

	DO i=1,m
	 q(i)=(x-y(i+l-m))*q(i)+(y(i+l)-x)*q(i+1)
	ENDDO

	RETURN

	END SUBROUTINE bval0
