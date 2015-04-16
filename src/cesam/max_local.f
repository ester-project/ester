
c******************************************************************
	
	SUBROUTINE max_local(f,nf,x,nx,xmin,xmax,fmax,fmin)

c routine PUBLIC du module mod_numerique	
c détection des maxima de f(nf,nx) dans l'intervalle [xmin,xmax]	

c entrées :
c	f(nf,nx) : fonctions	
c	x(nx) : abscisses
c	xmin,xmax : intervalle considéré

c sorties :
c	fmax(nf),fmin(nf) : MAX, MIN sur l'intervalle	
	
c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k	

c---------------------------------------------------

	USE mod_kind

	IMPLICIT NONE

	REAL (kind=sp), INTENT(in), DIMENSION(nf,nx) :: f
	REAL (kind=sp), INTENT(in), DIMENSION(nx) :: x
	REAL (kind=sp), INTENT(in) :: xmin, xmax
	INTEGER, INTENT(in) :: nf, nx
	REAL (kind=sp), INTENT(out), DIMENSION(nf) :: fmax
	REAL (kind=sp), INTENT(out), OPTIONAL, DIMENSION(nf) :: fmin	

	INTEGER :: i
	
c-----------------------------------------------------------------------

2000	FORMAT(8es10.3)
	
	fmax=-HUGE(1.) ; IF(PRESENT(fmin))fmin=HUGE(1.)	
	DO i=1,nx
	 IF(x(i) == MIN(xmax,MAX(x(i),xmin)))THEN
	  fmax(:)=MAX(fmax(:),f(:,i))
	  IF(PRESENT(fmin))fmin(:)=MIN(fmin(:),f(:,i))
	 ENDIF
	ENDDO
		
	RETURN
		
	END SUBROUTINE max_local
