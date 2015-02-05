
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	MODULE mod_donnees

c	module contenant les quantités fixées au cours de l'évolution
c	- du fichier de données initialisées dans lit_nl
c	- les contantes physiques initialisées dans ini_ctes_phys
c	- les paramètres de précision initialisées dans  cesam
c	- les paramètres d'évolution initialisées dans cesam
c	- les paramètres de composition chimique initialisées dans les
c       routines de réac. nuc., etc...	

c       paramètres public:
c	dtmin : pas temporel minimum en Myrs
c	n_min : nombre minimum de couches
c	pnzc : nombre max de zones convectives
c	r_qs : ordre des eq. diff. pour l'équilibre quasi-statique
c	version : numéro de version

c       variables public:

c---------------------real dp--------------------------

c	ab_ini : abondances initiales, initialisé dans nuc 
c	ab_min : abondances minimales, initialisé dans nuc 
c	nucleo : masses atomiques des éléments, initialisé dans tabul_nuc  
c	zi : charges des éléménts, initialisé dans tabul_nuc
	
c	abe7 : masse atomique en amu du béryllium 7, initialisé dans
c	ini_ctes 
c	abe9 : masse atomique en amu du béryllium 9, initialisé dans
c	ini_ctes 
c	ac12 : masse atomique en amu du carbone 12, initialisé dans
c	ini_ctes 
c	ac13 : masse atomique en amu du carbone 13, initialisé dans
c	ini_ctes
c	afe56 : masse atomique en amu du fer 56, initialisé dans ini_ctes
c	af18 : masse atomique en amu du fluor 18, initialisé dans ini_ctes
c	af19 : masse atomique en amu du fluor 19, initialisé dans ini_ctes
c	agemax : age maximum à atteindre, initialisé dans lit_nl
c	ah : masse atomique en amu de l'hydrogène, initialisé dans ini_ctes
c	ah2 : masse atomique en amu du deutérium, initialisé dans ini_ctes
c	ahe3 : masse atomique en amu de l'hélium 3, initialisé dans
c	ini_ctes
c	ahe4 : masse atomique en amu de l'hélium 4, initialisé dans
c	ini_ctes
c	ali6 : masse atomique en amu du lithium 6 , initialisé dans
c	ini_ctes
c	ali7 :  masse atomique en amu du lithium 7, initialisé dans
c	ini_ctes
c	alpha : longueur de mélange, initialisé dans lit_nl 
c	amg23 : masse atomique en amu du magnésium 23, initialisé
c	dans ini_ctes
c	amg24 : masse atomique en amu du magnésium 24, initialisé
c	dans ini_ctes
c	amg25 : masse atomique en amu du magnésium 25, initialisé
c	dans ini_ctes
c	amg26 : masse atomique en amu du magnésium 26, initialisé
c	dans ini_ctes
c	amu : masse atom. unite, Avogadro=1/amu, initialisé dans ini_ctes
c	an : masse atomique en amu du neutron, initialisé dans ini_ctes
c	ana23 : masse atomique en amu du sodium 23, initialisé dans
c	ini_ctes
c	ane20 : masse atomique en amu du néon 20, initialisé dans ini_ctes
c	ane21 : masse atomique en amu du néon 21, initialisé dans ini_ctes
c	ane22 : masse atomique en amu du néon 22, initialisé dans ini_ctes
c	an13 : masse atomique en amu de l'azote 13, initialisé dans
c	ini_ctes
c	an14 : masse atomique en amu de l'azote 14, initialisé dans
c	ini_ctes
c	an15 : masse atomique en amu de l'azote 15, initialisé dans
c	ini_ctes
c	ap : masse atomique en amu du proton, initialisé dans ini_ctes
c	ao16 : masse atomique en amu de l'oxygène 16, initialisé dans
c	ini_ctes
c	ao17 : masse atomique en amu de l'oxygène 17, initialisé dans
c	ini_ctes
c	ao18 : masse atomique en amu de l'oxygène 18, initialisé dans
c	ini_ctes
c	aradia : cte. de la radiation, initialisé dans ini_ctes
c	clight : célérite de la lumière, initialisé dans ini_ctes
c	cpturb : paramètre de pression turbulente, initialisé dans lit_nl
c	ctem : facteur de répartition en masse, initialisé dans cesam
c	ctep : facteur de répartition en pression, initialisé dans cesam
c	cter : facteur de répartition en rayon, initialisé dans cesam
c	dpsi : variation max pour modification de n_qs, initialisé
c	dans cesam
c	dn_fixe : taux pour modif de grille fixe en comp.chim, initialisé
c	dans cesam
c	dtlist : intervalle de temps entre deux listings complets,
c	initialisé dans lit_nl
c	dtmax : pas temporel maximum, initialisé dans cesam
c	dt0 : pas temporel initial, initialisé dans cesam
c	pmw : paramètre libre de perte de moment cinétique
c	d_grav : variation maximale du TdS, initialisé dans cesam
c	d_turb : coefficient de diffusion turbulente, initialisé dans
c	lit_nl
c	echarg : charge de l'électron, initialisé dans ini_ctes
c	eve : électron volt, initialisé dans ini_ctes
c	g : gravité déduit de gmsol et de Msol, initialisé dans ini_ctes
c	gmsol : valeur observée de G Msol, initialisé dans ini_ctes
c	granr : constante des gaz parfaits, initialisé dans ini_ctes
c	hpl : cte. de Planck, initialisé dans ini_ctes
c	kbol : cte. de Boltzman, initialisé dans ini_ctes
c	lbol0 : point 0 des Mbol, initialisé dans ini_ctes
c	ln10 : ln(10), initialisé dans ini_ctes
c	loc_zc : précision de la localisation des limites ZR/ZC,
c	initialisé dans cesam
c	log_teff : limite en Teff, initialisé dans lit_nl
c	lsol : luminosité solaire, initialisé dans ini_ctes
c	mdot : taux de perte de masse, initialisé dans lit_nl
c	me : masse électron, initialisé dans ini_ctes
c	msol : masse solaire, initialisé dans ini_ctes
c	mtot : masse initiale, initialisé dans lit_nl
c	ovshti : taux d'overshooting inférieur, initialisé dans lit_nl
c	ovshts : taux d'overshooting supérieur, initialisé dans lit_nl
c	pi : pi, initialisé dans ini_ctes
c	precit : précision pour l'intégration temporelle de la composition
c	chimique, initialisé dans cesam
c	precix : précision pour intégration équilibre quasi statique,
c	initialisé dans cesam
c	re_nu : paramètre de diffusivité radiative, initialisé dans lit_nl
c	ro_test : test de variation du TdS si ro > ro_test, initialisé
c	dans cesam
c	rsol : rayon solaire, initialisé dans ini_ctes
c	secon6 : nb. de s. en 10**6 ans, initialisé dans ini_ctes
c	sigma : cte. de Stefan, initialisé dans ini_ctes
c	tau_max : épaisseur optique au fond de l'atmosphère, initialisé
c	dans lit_nl
c	t_inf, t_sup : limtes de la tabulation des réactions
c	thermonucléaires, initialisé dans tabul_nuc
c	t_stop : limite supérieure de la température centrale, initialisé
c	dans lit_nl
c	w_rot : vitesse angulaire initiale, initialisé dans lit_nl
c	x0 : abondance initiale en hydrogène, initialisé dans lit_nl
c	x_stop : valeur limite de l'abondance centrale en  hydrogène,
c	initialisé dans lit_nl
c	y0 : abondance initiale en hélium, initialisé dans lit_nl
c	zsx_sol : Z/X solaire, initialisé dans ini_ctes
c	zsx0 : valeur initiale de Z/X, initialisé dans lit_nl
c	z0 : abondance initiale en métaux, initialisé dans lit_nl

c-------------------real sp-----------------------------

c	dfesh_des : barres d'erreur pour le dessin de [Fe/H], initialisé
c	dans des
c	dl_des : barres d'erreur pour le dessin de LOG L, initialisé dans
c	dans des
c	dteff_des : barres d'erreur pour le dessin de LOG Teff, initialisé
c	dans des
c	zoom_l : extensions asymétriques en luminosité, initialisé dans des
c	zoom_t : extensions asymétriques en Teff, initialisé dans des	

c	fesh_des : valeur de [Fe/H] en surface, initialisé dans des
c	l_des : cible en L/Lsol, initialisé dans des
c	teff_des : cible en Teff, initialisé dans des	

c--------------------integer------------------------

c	ihe4 : indice de l'hélium 4, initialisé dans tabul_nuc
c	ini0 : nb. iter. N-R avec réestim. de comp.chim. lim. ZR/ZC,
c	initialisé dans cesam
c	Krot : indice de la vitesse angulaire, initialisé dans lit_nl
c	m_ch : ordre des splines pour interpolation de la comp.chim.
c	initialisé dans cesam
c	m_ptm : ordre des spl. pour inter. de la masse (perte de masse),
c	initialisé dans cesam
c	m_qs : ordre des splines pour les variables quasi-static,
c	initialisé dans cesam
c	m_tds : ordre des splines pour interpolation de TdS, initialisé
c	dans cesam
c	nchim : nombre d'éléments strictement chimiques, initialisé dans
c	nuc
c	ne : nombre d'équations pour l'équilibre quasi-statique,
c	initialisé dans cesam
c	n_atm : nombre de couches de l'atmosphère, initialisé dans cesam
c	n_max : nombre maximum de couches pour l'équilibre quasi-statique,
c	initialisé dans cesam
c	ordre : ordre d'intégration de l'évolution temporelle de la
c	composition chimique, sans diffusion, initialisé dans cesam
c	ord_qs : ordre d'intégration  pour l'équilibre quasi-statique,
c	initialisé dans cesam

c-------------logical-----------------------

c	diffusion : il y a diffusion microscopique, initialisé dans lit_nl
c	en_masse : utilisation de la variable lgrangienne, initialisé dans
c	cesam
c	grille_fixe : utilisation d'une grille fixe pour la composition
c	chimique, initialisé dans lit_nl
c	jpz : overshoot suivant JpZh, initialisé dans lit_nl
c	kipp : utilisation de l'approximation de Kippenhahn, initialisé
c	dans cesam
c	ledoux : utilisation du critère de Ledoux , initialisé dans lit_nl
c	lim_ro : limite en densité pour l'atmosphère, initialisé dans
c	lit_nl
c	mitler : utilisation de l'écrantage suivant Mitler, initialisé dans
c	lit_nl
c	mvt_dis : ajustement de comp.chim. suivant mvt. des discontinutés,
c	initialisé dans cesam
c	pturb : utilisation de la pression turbulente, initialisé dans
c	cesam
c	rep_atm : on reprend l'atmosphère en binaire, initialisé dans cesam
c	rot_solid : il y a rotation solide, initialisé dans lit_nl
	
c---------character------------

c	precision : type de calcul, initialisé dans lit_nl
c	nom_elem : nom des éléments du vecteur de composition chimique
c	généralisé, initialisé dans tabul_nuc
c	arret : type d'arrêt, initialisé dans lit_nl
c	nom_atm : désignation de la routine de restitution d'atmosphère,
c	initialisé dans lit_nl
c	nom_abon : désignation de la routine d'abondances initiales,
c	initialisé dans lit_nl
c	nom_conv : désignation de la routine de convection, initialisé
c	dans lit_nl
c	nom_ctes : désignation de la routine de constantes physiques,
c	initialisé dans lit_nl
c	nom_des : désignation de la routine de dessin, initialisé dans
c	lit_nl
c	nom_diffm : désignation de la routine de diffusion micoscopique,
c	initialisé dans lit_nl
c	nom_difft : désignation de la routine de diffusion turbulente,
c	initialisé dans lit_nl
c	nom_etat : désignation de la routine d'équation d'état, initialisé
c	dans lit_nl
c	nom_nuc : désignation de la routine de réactions thermonucléaires,
c	initialisé dans lit_nl
c	nom_nuc_cpl : désignation de la compilation de réactions
c	thermonucléaires, initialisé dans lit_nl
c	nom_output : désignation du type de fichier ASCII en sortie,
c	initialisé dans lit_nl
c	nom_pertm : désignation de la routine de perte de masse, initialisé
c	dans lit_nl
c	nom_pertw : désignation de la routine de perte de moment cinètique,
c	initialisé dans lit_nl
c	nom_tdetau : désignation de la routine de loi T(tau), initialisé
c	dans lit_nl
c	nom_fich2 : identificateur des fichiers du modèle, initialisé dans
c	cesam
c	f_eos : noms des fichiers d'opacité, initialisé dans lit_nl
c	f_opa : noms des fichiers d'équation d'état, initialisé dans lit_nl
c	nom_chemin : chemin du directory SUN_STAR_DATA des données,
c	initialisé dans lit_nl
c	nom_opa : désignation de la routine d'opacité, initialisé dans
c	lit_nl
c	source : nom de la source des constantes physiques, initialisé
c	dans ini_ctes
c	methode : description de la méthode de calcul, initialisé dans
c	cesam

c       NOTATIONS (hélas incohérentes) pour les développements sur B-splines
c	n_ch : nombre VARIABLE de points élément de mod_variables
c	nch : nombre FIXE de fonctions élément de mod_donnees
c	m_ch : ordre FIXE des splines élément de mod_donnees 
c	mch(n_ch) : abscisses VARIABLES élément de mod_variables

c       Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c       CESAM2k

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	USE mod_kind
	
	IMPLICIT NONE

c       paramètres public:

	REAL (kind=dp), PARAMETER, PUBLIC :: dtmin=5.d-7, dx_tams=1.d-4,
	1    x_tams=0.01d0
	INTEGER, PARAMETER, PUBLIC :: n_min=150, pnzc=8, r_qs=1
	
c	CHARACTER (len=7), PARAMETER, PUBLIC :: version='V1.1.4'
	INCLUDE 'journal'
	
c       variables public

	REAL (kind=dp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:) :: ab_ini,
	1    ab_min, nucleo, rot_min, xvent, zi	
	REAL (kind=dp), SAVE, PUBLIC, DIMENSION(28) :: abon_m
	REAL (kind=dp), SAVE, PUBLIC, DIMENSION(5) :: pmw
	REAL (kind=dp), SAVE, PUBLIC :: abe7, abe9, ab11, ac12, ac13,
	1    afe56, af18, af19, agemax, ah, ah2, ahe3, ahe4, ali6, ali7,
	2    alpha, amg23, amg24, amg25, amg26, amu, an, ana23, ane20, ane21,
	3    ane22, an13, an14, an15, ap, ao16, ao17, ao18, aradia, beta_cgm, clight,
	4    cpturb, ctel, ctem, ctep, cter, ctet, dpsi, dn_fixe, dtlist,
	5    dtmax, dt0, d_conv=1.d13, d_grav, d_turb, echarg, eve, f_cgm, fesh_sol,
	6    g, gmsol, granr, he_core, hhe_core, hpl, kbol, lbol0, li_ini,
	7    lnt_stop, ln_Tli, ln10, loc_zc, log_teff, lsol, mdot, me, msol, mtot,
	8    ovshti, ovshts, pi, precit, precix, p_pertw, p_vent, pw_extend,
	9    re_nu, ro_test, rsol, secon6, sigma, tau_max, t_inf, t_sup,
	1    t_stop, w_form, w_rot, x0, x_stop, y0, zeta_cgm, zsx_sol, zsx0, z0,
	2    beta_v, zeta_v, frac_vp
	
	REAL (kind=sp), SAVE, PUBLIC, DIMENSION(2) :: dfesh_des, dl_des,
	1    dteff_des, zoom_l=0., zoom_t=0.

c       pour un écran 1280 X 1024	
	REAL (kind=sp), SAVE, PUBLIC :: dh=1.5, dl=2., h=7., ld=10.
	
c       pour un écran 1280 X 1600	
c	REAL (kind=sp), SAVE, PUBLIC ::  dh=2.5, dl=2.5, h=7., ld=11.3,
	
	REAL (kind=sp), SAVE, PUBLIC :: fesh_des=1000.,
	1    l_des=-100., teff_des=-100., logteff_max=-100.,
	2    logteff_min=-100., logl_max=-100., logl_min=-100.,
	3    xleft=1.8, ybot=1.4, y_age=1.3
	
	INTEGER, SAVE, PUBLIC :: Krot, ife56=0, ihe4, ini0, Ipg, i_ex,
	1    m_ch, m_ptm, m_qs, m_rl, m_rot, m_tds, nb_max_modeles, nchim,
	2    ne, nrot, n_atm, n_max, ordre, ord_qs

	LOGICAL, SAVE, PUBLIC :: diffusion, en_masse,
	1    garde_xish, grille_fixe, He_ajuste, jpz, kipp, ledoux, lim_ro,
	2    lov_ad=.TRUE., lvent=.TRUE., mitler, modif_chim, mvt_dis, pturb,
	3    rep_atm=.FALSE., rot_solid, t_ajuste, x_ajuste

	CHARACTER (len=2), SAVE, PUBLIC :: precision	
	CHARACTER (len=4), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:) ::
	1    nom_elem, nom_rot
	CHARACTER (len=4), SAVE, PUBLIC :: arret, nom_xheavy
	CHARACTER (len=5), SAVE, PUBLIC :: unit		
	CHARACTER (len=10), SAVE, PUBLIC :: langue	
	CHARACTER (len=20), SAVE, PUBLIC :: nom_atm, nom_abon,
	1    nom_conv, nom_ctes, nom_des, nom_diffm, nom_diffw, nom_difft,
	2    nom_etat, nom_frad, nom_nuc, nom_nuc_cpl, nom_output, nom_pertm,
	3    nom_pertw, nom_tdetau
	CHARACTER (len=31), SAVE, PUBLIC :: nom_fich2
	CHARACTER (len=50), SAVE, PUBLIC, DIMENSION(8) :: f_eos, f_opa
	CHARACTER (len=50), SAVE, PUBLIC :: nom_opa, source
	CHARACTER (len=80), SAVE, PUBLIC :: methode
	CHARACTER (len=100), SAVE, PUBLIC :: device='/xw'
	CHARACTER (len=255), SAVE, PUBLIC :: nom_chemin		

	PRIVATE
	PUBLIC :: lit_nl, ini_ctes, print_ctes

	CONTAINS

c-------------------------------------------------------------------
	
	INCLUDE 'ctes_85.f'
	INCLUDE 'ctes_94.f'
	INCLUDE 'ini_ctes.f'
	INCLUDE 'lit_nl.f'
	INCLUDE 'print_ctes.f'

	END MODULE mod_donnees
