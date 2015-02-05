 
c****************************************************************

	SUBROUTINE chim_gram(xchim,dxchim)
	
c routine public du module mod_variables

c	transforme les abondances par mole ---> par gramme
	
c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c entrées/sorties
c	xchim,dxchim : comp. chim et dérivée par mole ---> par gramme

c-----------------------------------------------------------------------

	USE mod_donnees, ONLY : nchim, nucleo
	USE mod_kind	

	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(inout), DIMENSION(:) :: dxchim, xchim	

c--------------------------------------------------------------
	
2000	FORMAT(8es10.3)

	xchim(1:nchim)=xchim(1:nchim)*nucleo(1:nchim)
	dxchim(1:nchim)=dxchim(1:nchim)*nucleo(1:nchim)

	RETURN
	
	END SUBROUTINE chim_gram
