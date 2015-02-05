
c*********************************************************************

	SUBROUTINE boite(x,y,dx,dy)

c	routine public du module mod_numerique
	
c	dessine une boite de dimensions +/- dx X dy autour de x,y

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c-------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE

	REAL (kind=sp), INTENT(in) ::	x, y, dx, dy
	
	REAL (kind=sp), DIMENSION(5) :: desx, desy

c---------------------------------------------------------------
		
	desx(1)=x-dx
	desx(2)=x+dx	
	desx(3)=x+dx	
	desx(4)=x-dx	
	desx(5)=x-dx	
 	desy(1)=y-dy	
 	desy(2)=y-dy
 	desy(3)=y+dy
 	desy(4)=y+dy
 	desy(5)=y-dy
 	
 	CALL pgsls(2)
 	CALL pgline(5,desx,desy)
 	CALL pgsls(1)
 	CALL pgpoint(1,x,y,16)
 	
 	RETURN
 	
 	END SUBROUTINE boite
