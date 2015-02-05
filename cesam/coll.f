
	SUBROUTINE coll(x,n,m,xcoll)

c	routine public du module mod_numerique
	
c	détermination des abscisses des points de collocation
c	de de Boor pour l'intégration des équations différentielles

c	il y a (n-1)*m points de collocation

c	dérivée de colloc.f

c entrées:
c	m: ordre des splines
c	n: nb. de points de raccord
c	x: abscisses des points de raccord

c sorties
c	xcoll: table des points de collocation

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	cesam2D

c------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: x
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: xcoll
	INTEGER, INTENT(in) :: m, n
		
	INTEGER :: i, l, ncoll	
	
c----------------------------------------------------------------	

2000	FORMAT(8es10.3)

c	WRITE(*,2000)x
	ncoll=(n-1)*m
	DO l=1,n-1
	 DO i=1,m	!(n-1)*(i-1)+l	!indice du point de collocation
	  xcoll((n-1)*(i-1)+l)=colpnt(i,l,m,x)
	 ENDDO		!ix
	ENDDO		!ipx
c	PRINT*,ncoll,n,m ; WRITE(*,2000)xcoll
	CALL shell(ncoll,xcoll)
	
c	verification de la stricte croissance de la suite des xcoll

	DO i=1,ncoll-1
	 IF(xcoll(i) >= xcoll(i+1))THEN
	  PRINT*,'dans coll la suite des points de collocation n''est pas',
	1 ' strictement croissante en i=',i
	  WRITE(*,2000)xcoll(1:i-1) ; PRINT*,xcoll(i-1),xcoll(i),xcoll(i+1)
	  WRITE(*,2000)xcoll(i+1:(n-1)*m) ; no_croiss=.TRUE. ; RETURN
	 ENDIF
	ENDDO
	
	RETURN
		
	END SUBROUTINE coll
