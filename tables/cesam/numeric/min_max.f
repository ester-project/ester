			
c***********************************************************************

	SUBROUTINE min_max(x,n,xmax,xmin)
	
c subroutine public du module mod_exploit
c calcul de min/max pour le plot, variables en simple précision

c entrées :
c	x : tableau
c	n : nombre de points

c sorties :
c	xmax, xmin : MAX, MIN étendus

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c---------------------------------------------------------------------

	USE mod_kind
	
	IMPLICIT NONE

	REAL (kind=sp), INTENT(in), DIMENSION(:) :: x
	INTEGER, INTENT(in) :: n			
	REAL (kind=sp), INTENT(out) :: xmax, xmin
	
	REAL (kind=sp) :: dx
	
c---------------------------------------------------------------------	


2000	FORMAT(8es10.3)

	xmax=MAXVAL(x(1:n)) ; xmin=MINVAL(x(1:n))
	
c	PRINT*,n ; WRITE(*,2000)x(1:n)
c	WRITE(*,2000)xmax ; WRITE(*,2000)xmin ; PAUSE'min_max'	

	IF(xmax == xmin)THEN
	 IF(x(1) == 0.d0)THEN
	  dx=0.01d0
	 ELSE
	  dx=0.01d0*x(1)
	 ENDIF
	ELSEIF(xmax+xmin == 0.d0)THEN
	 dx=(xmax-xmin)/100.d0
	ELSEIF((xmax-xmin)/(xmax+xmin) < 1.d-2)THEN
	 dx=0.01d0*(xmax+xmin)			
	ELSE
	 dx=(xmax-xmin)/100.d0
	ENDIF
	xmax=xmax+ABS(dx) ; xmin=xmin-ABS(dx)

	RETURN

	END SUBROUTINE min_max
