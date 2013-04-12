
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	MODULE mod_donnees

! module contenant les quantités fixées au cours de l'évolution
!	- du fichier de données initialisées dans lit_nl
!	- les contantes physiques initialisées dans ini_ctes
!	- les paramètres de précision initialisées dans  cesam
!	- les paramètres d'évolution initialisées dans cesam
!	- les paramètres de composition chimique initialisées dans les
!	  routines de réac. nuc., etc...

! NOTATIONS (hélas incohérentes) pour les développements sur B-splines
!	n_ch : nombre VARIABLE de points élément de mod_variables
!	nch : nombre FIXE de fonctions élément de mod_donnees
!	m_ch : ordre FIXE des splines élément de mod_donnees 
!	mch(n_ch) : abscisses VARIABLES élément de mod_variables

! La signification des variables est décrite au paragraphe F4 de la notice

! Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     
	USE mod_kind
      
	IMPLICIT NONE

! paramètres public : 
	REAL (kind=dp), PARAMETER, PUBLIC :: d_conv=1.d13,
	1 dx_tams=1.d-4, x_tams=0.01d0
	INTEGER, PARAMETER, PUBLIC :: nvth=10, n_min=150, pnzc=15, r_qs=1

	INCLUDE 'journal.f'
	
! variables public
	REAL (kind=dp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:) :: ab_ini,
	1 ab_min, nucleo, rot_min, xvent, zi	
	REAL (kind=dp), SAVE, PUBLIC, DIMENSION(28) :: abon_m
	REAL (kind=dp), SAVE, PUBLIC :: aal27, abe7, abe9, ab11, ac12, ac13,
	1 afe56, af18, af19, agemax, ah, ah2, ahe3, ahe4, ali6, ali7,
	2 alpha, amg23, amg24, amg25, amg26, amu, an, ana23, ane20, ane21,
	3 ane22, an13, an14, an15, ap, ap31, ao16, ao17, ao18, asi28, as32,
	4 aradia, clight, cpturb, ctel, ctem, ctep, cter, ctet, dlntc, dpsi,
	5 dn_fixe, dtlist, dtmax, dt0, dtmin=1.d-8, d_grav, d_turb, echarg, eve,
	6 fesh_sol, fmin_abon, g, gmsol, granr, he_core, hhe_core, hpl, kbol,
	7 lbol0, li_ini, lnt_stop, ln_Tli, ln10, loc_zc, log_teff, lsol, mdot,
	8 me, msol, mterre, mtot, ovshti, ovshts, pi, precit,
	9 precix, p_pertw, q0, re_nu, ro_test, rsol, secon6, sigma,
	1 tau_max, t_inf, t_sup, t_stop, ua, w_form=0.d0, w_rot, x0, x_stop, y0,
	2 zsx_sol, zsx0, z0
	
	REAL (kind=sp), SAVE, PUBLIC, DIMENSION(2) :: dfesh_des, dl_des,
	1 dteff_des, zoom_l=0., zoom_t=0.

! pour un écran 1280 X 1024	
	REAL (kind=sp), SAVE, PUBLIC :: dh=1.5, dl=2., h=7., ld=10.
	
! pour un écran 1600 X 1280	
!	REAL (kind=sp), SAVE, PUBLIC ::  dh=2.5, dl=2.5, h=7., ld=11.3,
	
	REAL (kind=sp), SAVE, PUBLIC :: fesh_des=1000.,
	1 l_des=-100., teff_des=-100., logteff_max=-100.,
	2 logteff_min=-100., logl_max=-100., logl_min=-100.,
	3 xleft=1.8, ybot=1.4, y_age=1.3
	
	INTEGER, SAVE, PUBLIC :: Kdes_rot, Krot, ife56=0, ihe4, iLi7=0,
	1 ini0, Ipg, i_ex, m_ch, m_ptm, m_qs, m_rot, m_tds, nb_max_modeles,
	2 nchim, ne, nrl, nrot, n_atm, n_max, ordre, ord_qs, ord_rot

	LOGICAL, SAVE, PUBLIC :: ajuste, all_output, all_rep, diffusion,
	1 en_masse, garde_xish, grad_ovi, grad_ovs, grille_fixe, He_ajuste,
	2 jpz, kipp, ledoux, lim_ro, lim_jpz, lisse, mitler, modif_chim,
	3 mu_saha, mvt_dis, nc_max=.FALSE., pturb, rep_atm=.FALSE., rot_solid,
	4 t_ajuste, x_ajuste

! baratine=.FALSE. permet de dérouter sur les fichiers 
! mon_modele_static, _atmos, _evol
! les informations concernant le déroulement des calculs en ce qui concerne,
! respectivement, la résolution de l'équilibre quasi-statique, la restitution de
! l'atmosphère, l'évolution de la composition chimique, de la vitesse angulaire.
!	LOGICAL, SAVE, PUBLIC :: baratine=.FALSE.
	LOGICAL, SAVE, PUBLIC :: baratine=.TRUE., l_demi, new_bv	

	CHARACTER (len=2), SAVE, PUBLIC :: precision	
	CHARACTER (len=4), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:) ::
	1 nom_elem, nom_rot
	CHARACTER (len=4), SAVE, PUBLIC :: arret, nom_xheavy
	CHARACTER (len=5), SAVE, PUBLIC :: unit		
	CHARACTER (len=10), SAVE, PUBLIC :: langue	
	CHARACTER (len=20), SAVE, PUBLIC :: nom_atm, nom_abon, nom_conv,
	1 nom_ctes, nom_des, nom_diffm, nom_diffw, nom_difft,
	2 nom_etat, nom_frad, nom_nuc, nom_nuc_cpl, nom_output, nom_pertm,
	3 nom_pertw, nom_tdetau, nom_thw
	CHARACTER (len=31), SAVE, PUBLIC :: nom_fich2
	CHARACTER (len=33), PARAMETER, PUBLIC, DIMENSION(0:5) :: thw=
	1 (/ 'sans rotation                    ',
	2    'rotation solide                  ',
	3    'rot. solide cons. glob. mnt. cin.',
	4    'diff. mnt.cin. Talon & Zahn 1997 ',
	5    'diff. mnt.cin. Matis & Zahn 2004 ',
	6    'cons. loc. mnt. cin.             ' /)	
	CHARACTER (len=50), SAVE, PUBLIC, DIMENSION(8) :: f_eos, f_opa	
	CHARACTER (len=50), SAVE, PUBLIC :: nom_opa, source
	CHARACTER (len=80), SAVE, PUBLIC :: methode
	CHARACTER (len=100), SAVE, PUBLIC :: device='/xw'
	CHARACTER (len=255), SAVE, PUBLIC :: nom_chemin		

	PRIVATE
	PUBLIC :: lit_nl, ini_ctes, print_ctes

	CONTAINS

!-------------------------------------------------------------------
 	
	INCLUDE 'ini_ctes.f'
	INCLUDE 'lit_nl.f'
	INCLUDE 'print_ctes.f'

	END MODULE mod_donnees
