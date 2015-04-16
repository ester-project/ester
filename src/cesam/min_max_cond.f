
c-----------------------------------------------------------------------

	SUBROUTINE min_max_cond(f,cond,n,xmax,xmin)

c subroutine public du module mod_exploit
	
c déterxmination de xmin et xmax sous condition pour des plots

c entrées :
c	n : nombre de points
c	f : tableau où les xmin xmax sont recherchés
c	cond : tableau des conditions

c sorties :
c	xmin : xminimum pour plot
c	xmax : xmaximum pour plot

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c----------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE
	
	REAL (kind=sp), INTENT(in), DIMENSION(:) :: f
	INTEGER, INTENT(in) :: n
	LOGICAL, INTENT(in), DIMENSION(:) :: cond
	REAL (kind=sp), INTENT(out) :: xmax, xmin

	REAL (kind=sp) :: dx
	
	INTEGER :: i

c------------------------------------------------------------------

2000	FORMAT(8es10.3)

c initialisation
	xmin=HUGE(1.) ; xmax=-xmin
	
	DO i=1,n
	 IF(cond(i))THEN
	  xmax=MAX(xmax,f(i)) ; xmin=MIN(xmin,f(i))
	 ENDIF	
	ENDDO
			
	IF(xmax == xmin)THEN
	 IF(f(1) == 0.d0)THEN
	  dx=0.01d0
	 ELSE
	  dx=0.01d0*f(1)
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
	
	END SUBROUTINE min_max_cond
