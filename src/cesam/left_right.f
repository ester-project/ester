
c**********************************************************

	SUBROUTINE left_right(n,f,xt,m,knot,xx,l,fx,dfx)

c routine PUBLIC du module mod_numerique

c calcul des n fonction f et des dérivées premières de part et d'autre de xx
c qui diffèrent quand xx est un point du vecteur nodal
c en entrée prendre l=m pour le premier appel, et réutiliser la valeur de sortie

c entrées
c	n: nombre de fonctions
c	xt : points nodaux
c	m: ordre des splines
c	xx: point d'interpolation
c      knot-m: dimension du vecteur nodal

c entrées/sorties
c	l: localisation entrée prendre l=m pour le premier appel,
c	et réutiliser la valeur de sortie

c sorties
c	fx(0,:)/dfx(0,:) valeurs/dérivées à gauche
c	fx(1,:)/dfx(1,:) valeurs/dérivées à droite

c Auteur: P. Morel, Département J.D. Cassini, O.C.A.

c-----------------------------------------------------------------------

	USE mod_kind
	        
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:,:) :: f
	REAL (kind=dp), INTENT(in), DIMENSION (:) :: xt	
	REAL (kind=dp), INTENT(in) :: xx
	
	INTEGER, INTENT(in) :: knot,m, n	

	INTEGER, INTENT(inout) :: l                   
	REAL (kind=dp), INTENT(out), DIMENSION(0:1,n) :: dfx, fx
      
	REAL (kind=dp), DIMENSION(m) :: d, xg
	REAL (kind=dp), DIMENSION(m+1) :: q
	REAL (kind=dp) :: xx1

	INTEGER :: i, j, l1

c----------------------------------------------------------------

2000	FORMAT(8es10.3)

c localisation de xx
	xx1=xx
	no_croiss=xx1 < xt(1) .OR. xx1 > xt(knot)
	IF(no_croiss)THEN
	 PRINT*,'dans left_right le point sort de la grille'
	 PRINT*,'xt(1)=',xt(1),' xx=',xx1,' xt(knot)=',xt(knot)
	 PRINT*,'n=',n,', m=',m,', knot=',knot
	 xx1=MIN(xt(knot),MAX(xx1,xt(1)))
	ENDIF
	
c les valeurs à droite, dans linf no_croiss sera .FALSE.	
	CALL linf(xx1,xt,knot,l) ; no_croiss=xx1 /= xx

c mises à 0
	fx=0.d0 ; dfx=0.d0

c fx(1,:)/dfx(1,:) valeurs/dérivées à droite
	CALL bval1(xx1,xt,m,l,q,d)
	DO j=1,n
	 DO i=1,m
	  fx(1,j)=fx(1,j)  +q(i)*f(j,l-m+i)
	  dfx(1,j)=dfx(1,j)+d(i)*f(j,l-m+i)
	 ENDDO
	ENDDO

c fx(0,:)/dfx(0,:) valeurs/dérivées à gauche
	xg(1:m)=xt(l-m+1:l)
	l1=l-COUNT(xg == xx)
	
c xx est un point du vecteur nodal	
	IF(l1 > 0)THEN
	 CALL bval1(xx1,xt,m,l1,q,d)
	 DO j=1,n
	  DO i=1,m
	   fx(0,j)=fx(0,j)  +q(i)*f(j,l1-m+i)
	   dfx(0,j)=dfx(0,j)+d(i)*f(j,l1-m+i)
	  ENDDO
	 ENDDO
	 
c xx est un point quelconque 	  
	ELSE
	 fx(0,:)=fx(1,:) ; dfx(0,:)=dfx(1,:)
	ENDIF	 

	RETURN

	END SUBROUTINE left_right
