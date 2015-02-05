
c**********************************************************

	SUBROUTINE bsp1dn(nf,f,x,xt,n,m,knot,init,xx,l,fx,dfxdx,contour)

c     routine public du module mod_numerique

c     interpolation dans les tableaux f(nf,n) avec une spline polynomiale
c     d'ordre m > 1, au point xx : x(1) <= xx <= x(n)
c     et  xt(l) <= xx <= xt(l+1)

c     en entrée prendre l quelconque : 1 <= l <= n+2m
c     au premier appel, init=.false. il y a initialisation :
c     calcul de knot=2m+n-2, formation du tableau xt(knot) des points de table
c     et de f(nf,n) coefficients des B-splines
c     !!le tableau des données f est donc modifié!!

c     version F95 de sbsp1dn.f

c     Auteur: P. Morel, Département J.D. Cassini, O.C.A.

c entrées
c     nf: nombre de fonctions
c     x: abscisses
c     n: nombre  de points
c     m: ordre des splines
c     xx: point d'interpolation
c     contour(présent) et .TRUE. création du vecteur nodal pour lissage par contour

c entrées/sorties
c     f(nf,n): fonctions a interpoler / coefficients des splines
c     xt(2m+n-2): points de raccord
c     l: localisation
c     knot: nb de noeuds

c sorties
c     fx, dfdx: fonctions, dérivées 1-ières

c-----------------------------------------------------------------------

	USE mod_kind
          
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) :: xx
	REAL (kind=dp), INTENT(in), DIMENSION (:) :: x
	INTEGER, INTENT(in) :: m, n, nf
	LOGICAL, INTENT(in) :: init
	LOGICAL, INTENT(in), OPTIONAL :: contour	
	REAL (kind=dp), INTENT(inout), DIMENSION(:,:) :: f
	REAL (kind=dp), INTENT(inout), DIMENSION (:) :: xt
	INTEGER, INTENT(inout) :: knot, l                   
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: dfxdx, fx

	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: a        
	REAL (kind=dp), DIMENSION(m) :: d
	REAL (kind=dp), DIMENSION(m+1) :: q
	REAL (kind=dp) :: xx1

	INTEGER, DIMENSION(:), ALLOCATABLE :: indpc
	INTEGER :: i, j

	LOGICAL :: cont, inversible

c----------------------------------------------------------------

2000	FORMAT(8es10.3)

c	si contour est présent et est .TRUE. lissage par contour
c	il n'y a que création du vecteur nodal	

	IF(PRESENT(contour))THEN
	 cont=contour
	ELSE
	 cont=.FALSE.
	ENDIF

c     initialisation des B-splines

	IF(.NOT.init)THEN
	 CALL noein(x,xt,n,m,knot) ; IF(no_croiss .OR. cont)RETURN
	 ALLOCATE(a(n,m),indpc(n))       
	 l=m
	 DO i=1,n
	  CALL linf(x(i),xt,knot,l) ; CALL bval0(x(i),xt,m,l,q)
	  DO j=1,m
	   a(i,j)=q(j)
	  ENDDO !j
	  indpc(i)=l-m+1
	 ENDDO  !i
	 CALL gauss_band(a,f,indpc,n,n,m,nf,inversible)
	 IF(.NOT.inversible)THEN
	  WRITE(*,1)
1	  FORMAT('dans bsp1dn la matrice a est singulière')
	  no_croiss=.TRUE. ; RETURN
	 ENDIF
	 DEALLOCATE(a,indpc)
	ENDIF   !initialisation des B-splines

c     localisation de xx
	xx1=xx
	no_croiss=xx1 < xt(1) .OR. xx1 > xt(knot)
	IF(no_croiss)THEN
	 PRINT*,'dans bsp1dn le point sort de la grille'
	 PRINT*,'xt(1)=',xt(1),' xx=',xx1,' xt(knot)=',xt(knot)
	 PRINT*,'n=',n,', m=',m,', knot=',knot
	 xx1=MIN(xt(knot),MAX(xx1,xt(1)))
	ENDIF
	CALL linf(xx1,xt,knot,l)

c     interpolation par B-splines

	CALL bval1(xx1,xt,m,l,q,d)
	fx=0.d0 ; dfxdx=0.d0
	DO j=1,nf
	 DO i=1,m
	  fx(j)=fx(j)      +q(i)*f(j,l-m+i)
	  dfxdx(j)=dfxdx(j)+d(i)*f(j,l-m+i)
	 ENDDO
	ENDDO

	RETURN

	END SUBROUTINE bsp1dn
