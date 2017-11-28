
c**********************************************************

	SUBROUTINE noedif(x,n,m,r,xt,knot)
	
c subroutine public du module mod_numerique		

c forme le tableau xt des knotx points nodaux à partir du
c tableau x de n points de raccord pour l'integration d'une
c equation differentielle d'ordre r par collocation avec
c des B-splines d'ordre m
c on place m+r points en x(1) et en x(n), sinon 1 point

c Auteur: P. Morel, Departement J.D. Cassini, O.C.A.

c entrées:
c	x: table des points de raccord
c	n: nombre de points de raccord
c	r: ordre de l'equa. diff.

c sorties:
c	xt: vecteur nodal
c	knot: dimension du vecteur nodal

c--------------------------------------------------------------------

	USE mod_kind
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: x
	INTEGER, INTENT(in) :: m, n, r	
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: xt
	INTEGER, INTENT(out) :: knot
	
	INTEGER :: i, j
	
c--------------------------------------------------------------------

2000	FORMAT(8es10.3)

c vérification de la stricte croissance de la suite des x
	DO i=1,n-1
	 IF(x(i) >= x(i+1))THEN
	  PRINT*,'dans noedif la suite des abscisses n''est pas',
	1 ' strictement croissante en i=',i
	  WRITE(*,2000)x(1:i-1) ; PRINT*,x(i-2),x(i),x(i+1)
	  WRITE(*,2000)x(i+1:n) ; no_croiss=.TRUE. ; RETURN
	 ENDIF
	ENDDO

c formation du tableau des xt
	knot=0
	DO i=1,m+r
	 knot=knot+1 ; xt(knot)=x(1)
	ENDDO

	DO i=2,n-1
	 DO j=1,m
	  knot=knot+1 ; xt(knot)=x(i)
	 ENDDO
	ENDDO

	DO i=1,m+r
	 knot=knot+1 ; xt(knot)=x(n)
	ENDDO

	RETURN

	END SUBROUTINE noedif	
