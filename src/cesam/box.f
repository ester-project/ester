
c*********************************************************************

	SUBROUTINE box(x,dx,y,dy)

c	routine public du module mod_numerique
	
c	dessine une boite asymétrique autour de x,y

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c-------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE

	REAL (kind=sp), INTENT(in), DIMENSION(2) :: dx, dy	
	REAL (kind=sp), INTENT(in) ::	x, y
	
	REAL (kind=sp), DIMENSION(5) :: desx, desy

c---------------------------------------------------------------
 	RETURN
 	
 	END SUBROUTINE box
