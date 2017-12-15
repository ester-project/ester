
c******************************************************************

	SUBROUTINE noein(x,xt,n,m,knot)
	
c subroutine public du module mod_numerique	

c	détermine la séquence de noeuds de raccord xt(knot)
c	pour une interpolation "optimale"
c	par B-splines d'ordre m sur la suite strictement
c	croissante x de n points de donnée, cf. de Boor p.219 formule (10)
c	aux limites le polynôme d'interpolation s'appuie sur m points
c	de donnée

c	Auteur: P. Morel, Département J.D. Cassini, O.C.A.

c entrées:
c	x(n): table des abscisses
c	n: nombre d'abscisses
c	m: ordre des splines

c sorties:
c	xt(knot): points de table
c	knot: nombre de points de table, knot=n+m

c--------------------------------------------------------------------

	USE mod_kind
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION (:) ::	x
	INTEGER, INTENT(in) :: m, n	
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: xt	
	INTEGER, INTENT(out) :: knot
	
	REAL (kind=dp), PARAMETER :: dx=1.d-12
	REAL (kind=dp) :: eps, mm1
			
	INTEGER :: i, j
			
c-----------------------------------------------------------------

c pour l'interpolation spline il faut n >= m
	IF(n < m)THEN
	 WRITE(*,11)n,m ; STOP
11	 FORMAT('ARRET: dans noein n=',i3,' < ',i3,' = m')
	ENDIF

c vérification de la stricte croissance de la suite des x
	DO i=1,n-1
	 no_croiss=x(i) >= x(i+1)
	 IF(no_croiss)THEN
	  PRINT*,'dans noein la suite des abscisses n''est pas',
	1 ' strictement croissante en i=',i
	  PRINT*,'nombre de points: ',n ; PRINT*,'abscisse x=',x(i)
	  WRITE(*,1)x(1:i) ; WRITE(*,1)x(i) ; WRITE(*,1)x(i+1:n) ; RETURN
1	  FORMAT(3es22.15)	  
	 ENDIF
	ENDDO
	
c petit écart "normalisé"	
	eps=dx*ABS(x(1)+x(n))/2.d0	

c formation du tableau des xt
	mm1=m-1.d0 ; knot=0

c petit écart à gauche	
	DO i=1,m
	 knot=knot+1 ; xt(knot)=x(1)-eps
	ENDDO

c les points nodaux de de Boor
	DO i=1,n-m
	 knot=knot+1 ; xt(knot)=0.d0
	 DO j=i+1,i+m-1
	  xt(knot)=xt(knot)+x(j)
	 ENDDO
	 xt(knot)=xt(knot)/mm1
	ENDDO

c petit écart à droite
	DO i=1,m
	 knot=knot+1 ; xt(knot)=x(n)+eps
	ENDDO

	RETURN

	END SUBROUTINE noein
