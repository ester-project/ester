
c************************************************************************

	SUBROUTINE conv(r,krad,gravite,delta,cp,ro,hp,taur,gradrad,
	1 gradad,der,gradconv,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
	2 dgradhp,dgradtaur,dgradgrad,dgradgad,
	3 gam,dgamkra,dgamgra,dgamdel,dgamcp,dgamro,
	4 dgamhp,dgamtaur,dgamgrad,dgamgad)

c routine générique de calcul de la convection

c routine public du module mod_conv

c entrées
c	r : rayon
c	krad : conductivité radiative
c	gravite : gravité
c	delta : delta
c	cp : c_ p
c	ro : densité
c	hp : échelle de hauteur de pression
c	taur :épaisseur optique de la bulle convective
c	gradrad : gradient radiatif 
c	gradad : gradient adiabatique
c	der=.TRUE. : calcul des dérivées

c sorties 
c 	gradconv, dgrad*~: gradient de température et dérivées
c	gam, dgam*~: efficacité de la convection et dérivées.

c	Auteur: P. Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : nom_conv
	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) :: r, krad, gravite, delta, cp, ro, hp,
	1 taur, gradrad, gradad     
	LOGICAL, INTENT(in) :: der
	REAL (kind=dp), INTENT(out) :: gradconv, dgradkra, dgradgra,
	1 dgradel, dgradcp, dgradro, dgradhp, dgradtaur, dgradgrad,
	2 dgradgad, gam, dgamkra, dgamgra, dgamdel, dgamcp, dgamro,
	3 dgamhp, dgamtaur, dgamgrad, dgamgad

c---------------------------------------------------------------------

	SELECT CASE(nom_conv)
	CASE ('conv_a0')
	 CALL conv_a0(krad,gravite,delta,cp,ro,hp,taur,gradrad,gradad,der,
	1 gradconv,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
	2 dgradhp,dgradtaur,dgradgrad,dgradgad,
	3 gam,dgamkra,dgamgra,dgamdel,dgamcp,dgamro,
	4 dgamhp,dgamtaur,dgamgrad,dgamgad)
	CASE ('conv_cm')
	 CALL conv_cm(krad,gravite,cp,ro,hp,gradrad,gradad,der,
	1 gradconv,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
	2 dgradhp,dgradtaur,dgradgrad,dgradgad,
	3 gam,dgamkra,dgamgra,dgamdel,dgamcp,dgamro,
	4 dgamhp,dgamtaur,dgamgrad,dgamgad)
	CASE('conv_cml')
	 CALL conv_cml(r,krad,gravite,cp,ro,hp,gradrad,gradad,der,
	1  gradconv,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
	2  dgradhp,dgradtaur,dgradgrad,dgradgad,
	3  gam,dgamkra,dgamgra,dgamdel,dgamcp,dgamro,
	4  dgamhp,dgamtaur,dgamgrad,dgamgad)
	CASE ('conv_cm_reza')
	 CALL conv_cm_reza(krad,gravite,delta,cp,ro,hp,gradrad,gradad,der,
	1  gradconv,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
	2  dgradhp,dgradtaur,dgradgrad,dgradgad,
	3  gam,dgamkra,dgamgra,dgamdel,dgamcp,dgamro,
	4  dgamhp,dgamtaur,dgamgrad,dgamgad)
	CASE ('conv_cgm_reza')
	 CALL conv_cgm_reza(krad,gravite,delta,cp,ro,hp,gradrad,gradad,der,
	1  gradconv,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
	2  dgradhp,dgradtaur,dgradgrad,dgradgad,
	3  gam,dgamkra,dgamgra,dgamdel,dgamcp,dgamro,
	4  dgamhp,dgamtaur,dgamgrad,dgamgad)
	CASE ('conv_jmj')
	 CALL conv_jmj(krad,gravite,delta,cp,ro,hp,taur,gradrad,gradad,der,
	1  gradconv,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
	2  dgradhp,dgradtaur,dgradgrad,dgradgad,
	3  gam,dgamkra,dgamgra,dgamdel,dgamcp,dgamro,
	4  dgamhp,dgamtaur,dgamgrad,dgamgad)
	CASE DEFAULT
	 PRINT*,'routine de convection inconnue: ',nom_conv
	 PRINT*,'routines connues: conv_a0, conv_cm, conv_cml,'
	 PRINT*,'conv_cm_reza,conv_cgm_reza, conv_jmj'
	 PRINT*,'arrêt' ; STOP
	END SELECT
       
	RETURN
       
	END SUBROUTINE conv
