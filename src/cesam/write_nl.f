
c******************************************************************

	SUBROUTINE write_nl

c écriture du fichier de données *.don
c routine public du module mod_exploit

c utilisé, en particulier, pour une calibration après la modification
c des données

c variables numériques
c	mtot: masse totale initiale
c	mdot: taux de perte de masse Msol/an
c	tau_max: épaisseur optique au fond de l'atmosphere
c	agemax: âge maximun a atteindre
c	dtlist: intervalle de temps min. entre deux listings
c	complets du modèle
c	x_stop: arrêt si Xcentre < = xstop
c	log_teff: arrêt si Log(Teff) <> +/-log_teff
c	t_stop: arrêt si Tcentre > t_stop
c	x0: abondance de X initiale
c	y0: abondance de Y initiale
c	zsx0: rapport Z/X initial zsx0=0 X0 est utilise
c	d_turb: coefficient de diffusion turbulente
c	re_nu: paramètre de diffusivite radiative
c	w_rot: vitesse angulaire initiale
c	alpha: paramètre de longueur de melange
c	cpturb: coefficient de pression turbulente
c	ovshts: coefficient d'overshoot supérieur
c	ovshti: coefficient d'overshoot inférieur

c variables logiques
c	lim_ro=.true.: condition limite externe de l'atmosphère sur ro
c	grille_fixe=.true.: utilisation d'une grille fixe en masse
c	pour interpolation de la comp.chim.
c	rot_solid=.true.: la rotation est solide
c	jpz=.true.: overshoot selon JpZ
c	ledoux=.true.: utilisation du critère de Ledoux
c	diffusion=.true.: calcul avec diffusion
c	mitler=.true.: effet d'écran selon Mitler

c variables sous forme de chaines de caractères
c	precision: niveau de précision requis
c	arret: arrêt sur zams, post, cohe, coca ou autre
c	nom_des: nom de la routine de dessin on line à utiliser
c	nom_ctes_phys: nom la routine de constantes physiques à utiliser
c	nom_perte: nom la routine de perte de masse à utiliser
c	nom_tdetau: nom la routine de loit T(tau) à utiliser
c	nom_atm: nom la routine de restitution d'atmosphère à utiliser
c	nom_conv: nom la routine de convection à utiliser
c	nom_nuc: nom la routine de réactions nucléaires à utiliser
c	nom_nuc_cpl: nom la routine de compilation de reac. nuc. à utiliser
c	nom_abon: nom la routine d'abondances initiales à utiliser
c	nom_diffm: nom la routine de diffusion microscopique à utiliser
c	nom_difft: nom la routine de  diffusion turbulente à utiliser
c	nom_etat: nom la routine d'équation d'état à utiliser
c	nom_opa: nom la routine d'opacité à utiliser

c variables sous forme de chaines de caractères de noms de fichiers
c	f_eos: noms des fichiers d'équation d'état
c	f_opa: noms des fichiers d'opacité

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c---------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE

	CHARACTER (len=20) :: nom_des_rot
	CHARACTER (len=50) :: chain

	NAMELIST/nl_cesam/nom_chemin,nom_ctes,nom_des,nom_output,
	1 n_max,precision
	NAMELIST/nl_mass/mtot,nom_pertm,mdot
	NAMELIST/nl_evol/agemax,arret,dtlist,log_teff,nb_max_modeles,
	1 he_core,t_stop,x_stop
	NAMELIST/nl_chim/grille_fixe,nom_abon,modif_chim,garde_xish,
	1 x0,y0,zsx0
	NAMELIST/nl_conv/nom_conv,alpha,ovshts,ovshti,jpz,cpturb,ledoux
	NAMELIST/nl_diff/diffusion,nom_diffm,nom_difft,d_turb,re_nu,
	1 nom_frad
	NAMELIST/nl_rot/w_rot,unit,nom_diffw,nom_thw,nom_pertw,p_pertw,
	1 lim_jpz, nom_des_rot
	NAMELIST/nl_etat/nom_etat,f_eos
	NAMELIST/nl_opa/nom_opa,f_opa
	NAMELIST/nl_nuc/nom_nuc,nom_nuc_cpl,mitler
	NAMELIST/nl_atm/nom_atm,nom_tdetau,tau_max,lim_ro

c---------------------------------------------------------------------

c identification de nom_des_rot
	SELECT CASE(Kdes_rot)
	CASE(0)
	 nom_des_rot='no_des'
	CASE(1)
	 nom_des_rot='end_evol'
	CASE(2)
	 nom_des_rot='all_mod'
	CASE(3)
	 nom_des_rot='all_iter'
	CASE(4)
	 nom_des_rot='end_mod'
	CASE DEFAULT
	 nom_des_rot='no_des'
	END SELECT
	
c grad_ov*=.TRUE. gradient adiabatique dans les parties overshootées des ZC
c grad_ov*=.FALSE. gradient radiatif dans les parties overshootées des ZC
	IF(.NOT.grad_ovi)ovshti=-ovshti	
	IF(.NOT.grad_ovs)ovshts=-ovshts	

c écriture des NAMELISTS
	chain=TRIM(nom_fich2)//'.don'
	OPEN(unit=3,form='formatted',status='unknown',delim='apostrophe',
	1 file=TRIM(chain))

	WRITE(*,1)chain ; WRITE(2,1)chain
1	FORMAT(t10,'NAMELISTS du fichier: ',a,/)

	WRITE(3,nl_cesam) ; WRITE(*,nl_cesam)
	WRITE(3,nl_mass)  ; WRITE(*,nl_mass)
	WRITE(3,nl_evol)  ; WRITE(*,nl_evol)
	WRITE(3,nl_chim)  ; WRITE(*,nl_chim)
	WRITE(3,nl_conv)  ; WRITE(*,nl_conv)
	WRITE(3,nl_diff)  ; WRITE(*,nl_diff)
	WRITE(3,nl_rot)   ; WRITE(*,nl_rot)
	WRITE(3,nl_etat)  ; WRITE(*,nl_etat)
	WRITE(3,nl_opa)   ; WRITE(*,nl_opa)
	WRITE(3,nl_nuc)   ; WRITE(*,nl_nuc)
	WRITE(3,nl_atm)   ; WRITE(*,nl_atm)

	CLOSE(unit=3) ; WRITE(*,*)

	RETURN

	END SUBROUTINE write_nl
