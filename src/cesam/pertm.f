	
c************************************************************************

         SUBROUTINE pertm(dt)

c	routine générique pour le calcul de la perte/gain de masse  
c	il y a des appels différents suivant nom_pertm

c	routine private du module mod_static

c entrées
c	bp,q,n,qt,knot,chim,mc,nc,mct,knotc : modèle au temps t
c	dt : pas temporel
c	age : age

c entrée/sortie
c	mstar : masse totale avec perte de masse
c	en entrée au temps t, en sortie au temps t+dt 

c sorties
c	old_ptm, x_ptm, n_ptm, xt_ptm, knot_ptm : interpolation de
c	l'ancienne masse en fonction de la nouvelle (en m^2/3 ou m) 
c	old_ptm et x_ptm sont identiques ce n'est pas le cas si
c	on tient compte de la perte de masse due a E=mc**2

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c----------------------------------------------------------------

	USE mod_donnees, ONLY : nom_pertm
	USE mod_kind
	USE mod_variables, ONLY : sortie
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: dt

c---------------------------------------------------------------- 

	SELECT CASE(nom_pertm)
	CASE ('pertm_ext')
	 CALL pertm_ext(dt)
	CASE ('pertm_msol')
	 CALL pertm_msol(dt)
	CASE ('pertm_tot')
	 CALL pertm_tot(dt)
	CASE ('pertm_waldron')
	 CALL pertm_waldron(dt)	 
	CASE DEFAULT
	 PRINT*,'routine de perte/gain de masse inconnue: ',nom_pertm
	 PRINT*,'routines connues: pertm_ext, pertm_msol, pertm_tot'
	 PRINT*,'pertm_waldron'
	 PRINT*,'arrêt' ; CALL sortie
	END SELECT

	RETURN
	
	END SUBROUTINE pertm
