
c**********************************************************

	SUBROUTINE schu58_n(n,f,xx,xt,m,l,val)
	
c	subroutine private du module mod_numerique

c	interpolation pour n fonctions  par l'algorithme 5-8 p.194 de schumaker
c	au point xx, xt(l) .le. xx. lt. xt(l+1) par spline d'ordre m
c	dont les coefficients sont dans le tableau f

c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A.
c	25 11 00: version F95

c entrees
c	n, f: nombre de fonctions f en spline
c	xx, xt: abscisse d'interpolation et vecteur nodal
c	m: ordre des splines
c	l: localisation

c sortie
c	val: vecteur des interpolees

c---------------------------------------------------------------

	USE mod_kind
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION(:,:) :: f
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xt	 	
	REAL (kind=dp), INTENT(in) :: xx
	INTEGER, INTENT(in) :: l, m, n	
	REAL (kind=dp), INTENT(out), DIMENSION(n) :: val			
	 
	REAL (kind=dp), DIMENSION(n,m) :: cx		
	REAL (kind=dp) :: a1, a2	
	INTEGER :: i, j, k
	
c----------------------------------------------------------------

2000	FORMAT(8es10.3)

	cx=0.d0
	DO k=1,n
	 DO j=1,m
	  IF(j+l-m .gt. 0)cx(k,j)=f(k,j+l-m)
	 ENDDO
	ENDDO

	DO j=2,m
	 DO i=m,j,-1
	  a1=(xx-xt(i+l-m))/(xt(i+l-j+1)-xt(i+l-m)) ; a2=1.d0-a1
	  DO k=1,n
	   cx(k,i)=a1*cx(k,i)+a2*cx(k,i-1)
	  ENDDO
	 ENDDO
	ENDDO

	val=cx(:,m)

	RETURN

	END SUBROUTINE schu58_n
