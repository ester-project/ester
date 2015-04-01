
c*************************************************************************

	SUBROUTINE sum_n(n,f,xt,m,knot,init,c,d,sum)

c	somme de c à d des n fonctions f mise sous forme de spline
c	(par sbsp1dn par exemple) aux n points
c	x(1)<...<x(n), c avec d dans [x(1), x(n)]

c	!!! ATTENTION le tableau des splines f est modifie !!!!

c	Auteur: P. Morel, Département J.D. Cassini, O.C.A., Observatoire de Nice
c	25 11 00: version F95

c entrees
c	n: nombre de fonctions dans f
c	xt, m, knot: vecteur nodal, ordre des splines et DIMENSION
c	c, d: intervalle d'intégration
c	init=.TRUE.: le tableau des f est deja adapte a l'integration

c entrees/sorties
c	f(n,knot-m): coefficients des splines modifié si init=.FALSE.

c sorties
c	sum: vecteur des integrales

c---------------------------------------------------------------------
	
	USE mod_kind
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xt
	REAL (kind=dp), INTENT(in) :: c, d
	INTEGER, INTENT(in) :: knot, m, n			
	LOGICAL, INTENT(in) :: init	
	REAL (kind=dp), INTENT(inout), DIMENSION(:,:) :: f	
	REAL (kind=dp), INTENT(out), DIMENSION(n) :: 	sum
			
	REAL (kind=dp) :: val(n), bid
	
	INTEGER, SAVE :: l=1				
	INTEGER :: i, j

c---------------------------------------------------------------------

2000	FORMAT(8es10.3)

c l'intégrale de la spline est une spline d'ordre m+1
c calcul des ci !!!! mis dans f !!!! algorithme 5.19 de Schumaker
	IF(.NOT.init)THEN	!calcul de la spline de l'integrale
	 DO j=1,n
	  bid=0.d0
	  DO i=1,knot-m
	   bid=bid+(xt(i+m)-xt(i))*f(j,i) ; f(j,i)=bid
	  ENDDO
	 ENDDO
	ENDIF

c calcul de somme de c -> d de f(x) dx
	CALL linf(d,xt,knot,l)
	IF(no_croiss)PRINT*,'Pb en 1 dans sum_n'		
	CALL schu58_n(n,f,d,xt,m+1,l,sum)
	
	CALL linf(c,xt,knot,l)
	IF(no_croiss)PRINT*,'Pb en 2 dans sum_n'	
	CALL schu58_n(n,f,c,xt,m+1,l,val)
	
	sum=(sum-val)/REAL(m)

	RETURN

	END SUBROUTINE sum_n
