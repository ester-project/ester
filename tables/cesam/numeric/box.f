
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
		
	desx(1)=x+dx(2)
	desx(2)=x+dx(1)	
	desx(3)=x+dx(1)	
	desx(4)=x+dx(2)	
	desx(5)=x+dx(2)	
 	desy(1)=y+dy(2)	
 	desy(2)=y+dy(2)
 	desy(3)=y+dy(1)
 	desy(4)=y+dy(1)
 	desy(5)=y+dy(2)
 	
 	! CALL pgsls(2)
 	! CALL pgline(5,desx,desy)
 	! CALL pgsls(1)
 	! CALL pgpoint(1,x,y,16)
 	
 	RETURN
 	
 	END SUBROUTINE box
