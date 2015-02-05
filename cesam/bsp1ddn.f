
c**********************************************************

      SUBROUTINE bsp1ddn(nf,f,x,xt,n,m,knot,init,xx,l,fx)

c     routine public du module mod_mumerique

c     interpolation dans les tableaux f(nf,n) avec une spline polynomiale
c     d'ordre m > 1, au point xx : x(1) <= xx <= x(n)
c     et  xt(l) <= xx < xt(l+1)
c     avec calcul de toutes les dérivées

c     en entrée prendre l quelconque : 1 .le. l .le. n+2m
c     au premier appel, init=.false. il y a initialisation :
c     calcul de knot=2m+n, formation du tableau xt(knot) des points de table
c     et de f(n) coefficients des B-splines
c     !!le tableau des données f est donc modifié!!

c     Auteur: P. Morel, Departement J.D. Cassini, O.C.A.
c     version: 07 12 92 (F77), adaptée F95, nov.2002

c entrees
c     nf: nombre de fonctions
c     f(nf,n): fonctions a interpoler
c     x: abscisses
c     n: nombre  de points
c     m: ordre des splines
c     knot: nb de noeuds
c     xx: point d'interpolation

c entrees/sorties
c     f(nf,n): fonctions à interpoler/coefficients des splines
c     xt: points de raccord
c     l: localisation

c sorties
c     fx(i,der) : dérivée ordre der de la i-ieme fonction 
c     fx(3,1) : dérivée 1-iere de la 3-ieme fonction

c-----------------------------------------------------------------------

      USE mod_kind

      IMPLICIT NONE

      REAL (kind=dp), INTENT(in), DIMENSION(:) :: x
      REAL (kind=dp), INTENT(in) :: xx
      INTEGER, INTENT(in) :: n, m, nf
      LOGICAL, INTENT(in)  :: init

      REAL (kind=dp), INTENT(inout), DIMENSION(:,:) :: f
      REAL (kind=dp), INTENT(inout), DIMENSION(:) :: xt
      INTEGER, INTENT(inout) :: knot, l
      REAL (kind=dp), INTENT(out), DIMENSION(:,0:) :: fx

      REAL (kind=dp), DIMENSION(:,:), ALLOCATABLE :: a
      REAL (kind=dp), DIMENSION(0:m-1,m) :: d
      REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: q
      REAL (kind=dp) :: xx1

      INTEGER, ALLOCATABLE, DIMENSION(:) :: indpc
      INTEGER :: i, j, k

      LOGICAL :: inversible

c-----------------------------------------------------------------------

2000  FORMAT(8es10.3)

c     initialisations

      IF(.NOT.init)THEN       !les coefficient sont à calculer

       CALL noein(x,xt,n,m,knot) ; IF(no_croiss)RETURN      
       ALLOCATE(a(n,m),indpc(n),q(m+1))       
       l=m
       DO i=1,n
        CALL linf(x(i),xt,knot,l)
        CALL bval0(x(i),xt,m,l,q)
        DO j=1,m
         a(i,j)=q(j)
        ENDDO !j
        indpc(i)=l-m+1
c       PRINT*,indpc(i)
c       WRITE(*,2000)(a(i,j)
       ENDDO  !i
       CALL gauss_band(a,f,indpc,n,n,m,nf,inversible)
       IF(.NOT.inversible)THEN
        PRINT*,'dans bsp1ddn la matrice a est non inversible, arret'
        STOP
       ENDIF
       DEALLOCATE(a,indpc,q)
      ENDIF

c     localisation de xx

      IF(xx < xt(1) .OR. xx > xt(knot))THEN
       PRINT*,'dans bsp1ddn le point sort de la grille'
       PRINT*,'xt(1)=',xt(1),' xx=',xx,' xt(knot)=',xt(knot)
       PRINT*,'n=',n,' m=',m,' knot=',knot
      ENDIF
      xx1=MIN(xt(knot),MAX(xx,xt(1))) ; CALL linf(xx,xt,knot,l)

c     d(r,j):= dérivée d'ordre r,  0 <= r < m-1, de
c     la l-m+j-ieme B-splines, 1 <= j <= m, d'ordre m non nulle
c     au point x, xt(l) <=  x <  xt(l+1). 

      CALL bvald(xx1,xt,m,l,m-1,d) 

c     interpolation par B-splines

      fx=0.d0
      DO j=1,nf       !fonction
       DO i=0,m-1     !dérivée
        DO k=1,m      !spline
         fx(j,i)=fx(j,i)+f(j,l-m+k)*d(i,k)
        ENDDO
       ENDDO
      ENDDO

      RETURN

      END SUBROUTINE bsp1ddn
