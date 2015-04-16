
c*********************************************************************

	SUBROUTINE boite(x,y,dx,dy)

c routine public du module mod_numerique
	
c dessine une boite de dimensions +/- dx X dy autour de x,y

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c-------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE

	REAL (kind=sp), INTENT(in) ::	x, y, dx, dy
	
	REAL (kind=sp), DIMENSION(5) :: desx, desy

c---------------------------------------------------------------
 	RETURN
 	
 	END SUBROUTINE boite
