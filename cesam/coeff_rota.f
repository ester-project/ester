
c*************************************************************

	SUBROUTINE coeff_rota(dt,nu,y,frl,coeff)

c routine générique private du module mod_evol

c calcul des coefficients pour la rotation

c entrées :
c	dt : pas temporel
c	nu : abscisse en m^2/3
c	y : variables de la rotation

c sorties :
c	frl : coefficients de rotation
c	deff, dh, dv : coefficients de diffusion du moment cinétique

c	Auteur: P. Morel, Département Cassiopée, O.C.A., CESAM2k

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : Krot, nchim, nrl, nrot
	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(nrot,0:1) :: y
	REAL (kind=dp), INTENT(in) :: dt, nu

	REAL (kind=dp), INTENT(out), OPTIONAL, DIMENSION(ncoeff+nchim) :: coeff	
	REAL (kind=dp), INTENT(out), DIMENSION(nrl) :: frl
		
c---------------------------------------------------------------------

	SELECT CASE(Krot)
	CASE(3)
	 CALL coeff_rota3(dt,nu,y,frl,coeff)	
	CASE(4)
	 CALL coeff_rota4(dt,nu,y,frl,coeff)	
	CASE DEFAULT
	 PRINT*,'ERREUR, pas de coefficient de rotation pour Krot=',Krot ; STOP
	END SELECT
	
	RETURN

	END SUBROUTINE coeff_rota
