	
c******************************************************************      
      
	SUBROUTINE static(fait,cx,li,y,be,ae,compt,dt,reprend,ip)

c	subroutine générique du calcul des coefficients de l'équilibre
c	quasi-statique

c	routine private du module mod_static

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c------------------------------------------------------------------

	USE mod_donnees, ONLY : en_masse
	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:,0:) :: y
	REAL (kind=dp), INTENT(in) :: dt
	INTEGER, INTENT(in) :: cx, compt, fait, li, ip

	REAL (kind=dp), INTENT(out), DIMENSION(:,:,0:) :: ae
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: be
	LOGICAL, INTENT(out) :: reprend

c------------------------------------------------------------------

	IF(en_masse)THEN
	 CALL static_m (fait,cx,li,y,be,ae,compt,dt,reprend,ip)
	ELSE
	 CALL static_r (fait,cx,li,y,be,ae,compt,dt,reprend,ip)
	ENDIF

	RETURN
	
	END SUBROUTINE static
