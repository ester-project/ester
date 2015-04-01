
c***********************************************************************

	SUBROUTINE cesam

c routine public du module mod_cesam

c calcul de modèles de structure interne et d'évolution stellaire
c programme principal: gestion générale des calculs, des écritures

c méthode collocation de de Boor p. 277
c adapté à la structure interne, avec fonction de répartition

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c	07 10 96 : gravité effective dans vaissala
c	10 10 96 : introduction de la rotation dans thermo
c	21 03 97 : bug dans calcul de z_opa (signalé par Sacha)
c	08 05 97 : introduction de l'option kipp pour calcul du TdS
c	26 05 97 : fichiers cohe, coca, début de combustion helium, carbone
c	17 09 06 : fichiers coox, début de combustion oxygène
c	11 06 97 : simplification de la gestion de t_inf
c	26 06 97 : remplacement du moment angulaire par la vitesse angulaire
c	18 07 97 : correction d'énergie exédentaire de rotation
c	18 08 97 : modif nom elem chim
c	25 08 97 : mise en place des variables eulériennes
c	22 09 97 : abondances externes dans le fichier *.HR
c	04 11 97 : abondances centrales dans le fichier *.HR
c	18 11 97 : mstar dans le fichier *.HR
c	12 11 99 : calcul de mue avec les taux d'ionisation
c	19 11 99 : rectification de CALL time, des dates, de thermo,
c	utilisation de saha pour nh1, nhe1, nhe2, degene
c	18 01 00 : correction si e_nuc=0
c	02 02 00 : e_grav < 50% e_nuc ----> série pple (au lieu de 99%)
c	18 04 00 ; coeff_diff ---> diffm, difft
c	30 07 00 : introduction F95

c	avec pression turbulente 7 inconnues
c	sans pression turbulente 6 inconnues, Ptot=Pgaz

c	noms des routines génériques de physique
c	etat: équation d'état
c	opa: opacité
c	conv: convection
c	nuc: réactions thermonucléaires
c	atm: condition limite externe
c	tdetau: loi t(tau)
c	diffm: coefficients de diffusion microscopique
c	difft: coefficients de diffusion turbulente
c	diffw: coefficients de diffusion turbulente du moment angulaire
c	pertm: perte de masse
c	pertw: perte de moment cinétique
c	des : dessin du modèle durant l'exécution
c	ini_ctes : initialisation des constantes de physique

c	fonctions des principales routines internes numériques
c	resout : initialisation, formation, résolution du pb. aux limites
c	evol : gestion de l'évolution de la composition chimique
c	static_m/r : formation des équations de l'équilibre quasi statique
c	lim_zc : détermination des limites ZR/ZC
c	update : transfert t+dt ---> t ---> t-dt
c	diffus : calcul de la diffusion
c	coll_atm, eq_atm : restitution d'atmosphère

c	le centre est à la couche 1, la limite externe en couche n_qs
c	le nb. de couches est VARIABLE il est fixé dans le SSP lim_zc
c	une couche est placée à chaque limite ZR/ZC (à ~5% près)

c	pour avoir un algorithme unique pour la diffusion, la composition
c	chimique est toujours tabulée en fonction de mu=(m/Msol)^2/3
c	que ce soit en lagrangien ou en eulérien

c	l'énergie graviphique, TdS/dt=tds est tabulée en fonction
c	de m^2/3 en lagrangien et de m en eulérien,
c	incrément temporel : dt=delta t / 10**6 ans

c	variable indépendante : q l'indice (real) de couche
c	variables dépendantes :
c		ksi=ln P pression cgs
c		eta=ln T température K
c		dzeta=(r/rsol)**2
c		lambda=(l/lsol)**2/3
c		mu=(m/mtot)**2/3
c		psi=dQ/dq=cte avec Q(mu,t)=fonction de répartition

c	composition chimique : xchim(i),i=1,nchim : abondances / mole
c	en f(m/Mtot**2/3) + W  (w : vitesse angulaire)

c	ordre et indices d'identification des éléments definis dans
c	reac_nuc, Exemple:
c	H1 H2 He3 He4 Li7 Be7 C12 C13 N14 N15 O16 O17 W
c	1  2  3   4   5   6   7   8   9   10   11  12 13

c	avant le début d'un pas temporel:
c	mc_t,mct_t,nc_t,knotc_t,chim_t:   comp. chim. à âge-dt
c	mc,mct,nc,knotc,chim: comp. chim. à âge
c	bp,q,qt,knot,p,t,r,l,m:  var. pples. à âge
c	x_tds,xt_tds,n_tds,knot_tds,tds: TdS/dt à âge t

c	en_masse = .TRUE. variables lagrangiennes m23=m**23, r2=r**2
c	en_masse = .FALSE. variables eulériennes m23=m, r2=r

c	Fichiers utilisés

c	les fichiers de sortie ont le nom choisi par l'utilisateur
c	mon_modele suivit d'un attribut qui les distingue
c	l'extension "_B" est spécifique aux fichiers binaires

c	FOR002 : listing de sortie, créé dans cesam fermé et garde dans
c	cesam, nom : mon_modèle.lis

c	FOR003 : lecture des données, utilisé temporairement dans resout
c	nom : mon_modele.don

c	FOR004 : modèle initial ou repris (binaire),
c	utilisé temporairement dans cesam

c	FOR011 : lecture des tables d'opacité
c	utilisé temporairement dans routines d'opacité

c	FOR011-->017 : lecture des tables d'équation d'état
c	utilisé temporairement dans certains SSP d'équation d'état

c	FOR024 : écriture modèle initial d'âge 0 homogène en binaire,
c	PMS, ZAMS ouvert et fermé dans cesam
c	nom : mon_modele_B.hom/pms

c	FOR025 : dernier modèle, ouvert et fermé dans cesam
c	nom : mon_modele_B.dat en binaire

c	FOR026  point de reprise, ouvert et fermé dans cesam
c	nom : mon_modele_B.rep en binaire

c	FOR030  liste pour oscillation ouvert et fermé dans osc
c	nom : mon_modele.osc en ASCII

c	FOR031 : dernier modèle calculé de l'atmosphère restituée
c	ouvert et fermé dans lim_atm
c	nom : mon_modele_B.atm en binaire

c	FOR053 : fichier pour tracer un diag. HR et les ZC,
c	ouvert et fermé dans cesam
c	nom : mon_modele.HR en ASCII

c	FOR060 : tables d'EOS d'OPAL

c~~~~B-splines, interpolation, équations différentielles résolution par
c collocation, par éléments finis

c pour une interpolation et pour une équation équation différentielle
c résolue par éléments finis l'ordre des splines est noté m_*, exemple m_ch 
c pour une équation équation différentielle résolue par collocation
c l'ordre des splines est noté ord_*, exemple : ord_qs
c l'ordre est la somme d'un odre "primaire" noté m_*, exemple m_qs, et de
c l'ordre de l'équation différenteielle, exemple r_qs : ord_qs=m_qs+r_qs
c dans CESAM2k les équations différentielles résolues par collocation
c sont résolues sous la forme de systèmes du premier ordre ord_qs=1
c ce qui évite le calcul de dérivées d'ordre élevé des coefficients
 
c~~~~
c formation des points de collocation, qui serviront pour rota et poisson
c pour un ordre d'interpolation m_rot+r_qs=ord_rot, il y a m_rot points de
c collocation dans chaque intervalle, la multiplicité des noeuds étant mrot
c il y a ncoll_rot=(n_rot-1)*m_rot points de collocation et r_qs points limite
c r_qs ordre des équations différentielles
c ord_rot=m_rot+r_qs ordre d'interpolation
c avec nd discontinuités il convient d'ajouter nzc*r_qs points limite
c soit (nd+1)*r_qs points limite
c sans discontinuité interne dimension du vecteur nodal knotr=dim_rot+m_rot+r_qs
c avec discontinuités internes knotr=dim_rot+m_rot+r_qs(nzc+1)
c dim_rot, knotr,  éléments de mod_variables, définis dans base_rota
c ord_rot est défini dans base_rota
c m_rot, r_qs sont des éléments de mod_donnees

c voir ligne 4232 pour écrire un modèle d'initialisation en ASCII

c---------------------------------------------------------------------

	USE mod_atm, ONLY: atm, dlpp_atm, m_atm, p_atm, pt_atm,
	1 r_atm, tau, tdetau, thermo_atm, t_atm
	USE mod_donnees, ONLY: ab_min, ab_ini, agemax, ajuste, all_output,
	1 all_rep, alpha, aradia, arret, baratine, cpturb, ctel, ctem, ctep,
	2 cter, ctet, diffusion, dlntc, dn_fixe, dpsi, dtlist, dtmax, dtmin,
	3 dt0, d_grav, en_masse, fmin_abon, f_eos, f_opa, g, gmsol, granr,
	4 grille_fixe, he_core, ihe4, ini0, Ipg, Kdes_rot, Krot, kipp, langue,
	5 lim_ro, lit_nl, lisse, li_ini, lnt_stop, ln_Tli, loc_zc, LOG_teff,
	6 lsol, l_demi, methode, msol, mdot, mtot, mu_saha, mvt_dis, m_ch,
	7 m_ptm, m_qs, m_rot, m_tds, ne, nb_max_modeles, nchim, new_bv,
	8 nom_conv, nom_elem, nom_des, nom_etat, nom_fich2, nom_opa, nom_nuc,
	9 nom_nuc_cpl, nom_abon, nom_atm, nom_tdetau, nom_pertm, nom_pertw,
	1 nom_diffm, nom_difft, nom_diffw, nom_ctes, nom_output, nom_thw, nrot,
	2 nucleo, nvth, n_atm, n_max, ordre, ord_qs, ord_rot, precit, precix,
	3 print_ctes, pturb, pi, precision, q0, rep_atm, re_nu, 
	4 rsol, ro_test, rot_min, rot_solid, r_qs, secon6, sigma, thw,
	5 t_stop, unit, version, w_rot, x_tams, x_stop, x0, y0, zi, z0
	USE mod_etat, ONLY : etat, saha
	USE mod_evol, ONLY : ecrit_rota, initialise_rota, tab_vth
	USE mod_kind
	USE mod_nuc, ONLY: nuc, planetoides, vent
	USE mod_numerique, ONLY: arb_rom, bsp1ddn, bsp1dn, noedif, newspl,
	1 newton, no_croiss, sum_n, zoning
	USE mod_static, ONLY: iter_qs, resout, thermo
	USE mod_variables, ONLY: age, bp, chim, chim_gram, c_iben, dim_ch,
	1 dim_qs, dim_rot, jlim, knot, knotc, knotr, knot_tds,
	2 lconv, lhe_stop, lim, lt_stop, lx_stop, mc, mct, mc_fixe,
	3 model_num, mstar, mw_tot, mrot, mrott, m_zc, m23,
	4 nb_modeles, nc_fixe, n_ch, n_qs, n_rot, n_tds, psi0, q, qt,
	5 rota, rota_t, rstar, r_ov, r_zc, r2, tds, vth, wrot, x_tds,
	6 xt_tds

	IMPLICIT NONE

	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: compchim,
	1 compg, dcompchim, dcompg, ener, epsilon, fd, ioni, lnta
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: b, bb, jac, new_esp,
	1 var
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: alfa, alfa_atm,
	1 anube7, anub8, anun13, anuo15, anupep, anupp, beta, beta_atm,
	2 delta, delta_atm, cp, cp_atm, dcapdr, dcapdt, depsdr, depsdt, depsx,
	3 dm23dlnr, dlnpdlnr, dlntdlnr, dlpp, dxchim, dxchimg, gamma,
	4 gamma_atm, glob, grad, gradad, gradconv, gradrad, grad_mj,grad_mj_a,
	5 grad_atm, grad_mu, grada_atm, gradc_atm, gradr_atm, gtsg, hp, hp_atm,
	6 kap, k_atm, l, degene, ldcapdr_a, ldcapdt_a, lnpt_a, lnpt_at,
	7 lnt_a, m, mu, mu_atm, mue, mue_atm, p, pt, r, ro, ro_atm, t, u,
	8 vaissala, vais_atm, w, xchim, xchimg, xchim1, xchim1g, z, zbar
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: comp, dcomp, esp,
	1 ex, mt, qb, qbt
	REAL (kind=dp), DIMENSION(nvth) :: dfvth, fvth
	REAL (kind=dp), DIMENSION(8) :: dfqs, dfrot, fqs, frot
	REAL (kind=dp), DIMENSION(0:5) :: absc, fonc, poly
	REAL (kind=dp), DIMENSION(5) :: epsilo
	REAL (kind=dp), DIMENSION(4) :: dfen, fen
	REAL (kind=dp), DIMENSION(1) :: dfatm, dftds, fatm, ftds
	REAL (kind=dp), SAVE :: agep=1.d-30, cte1, cte2, cte3, cte4, cte5,
	1 cte6, cte8, cte9, dt, dtp=-100.d0, dts=0.d0, lteffp=-100.d0, n_init
	REAL (kind=dp) :: alphap, bid, be7, b8, dcpp, dcpt, dcpx,
	1 deltap, deltat, deltax, depsp, depst, df_tau,dgaml, dgamlpp,
	2 dgamm, dgamp, dgampt, dgamr, dgamrs, dgamt, dgamtau,dgamx,
	3 dgradadp, dgradadt, dgradadx, dgradl, dgradlpp, dgradm, dgradp,
	4 dgradpt, dgradr, dgradrs, dgradt, dgradtau, dgradx, dhpm, dhpp,
	5 dhppt, dhpt, dhpr, dhpx, dkapp, dkapt, dkapx, dml, dmr,
	6 dpl, dpr, dptdl, dptdr, dro_grav, drop, drot, dro_teff, drox,
	7 dtdl, dtdr, dtsdg, dtsdteff, dtsdtau, dup, dut, dux,d2f_tau,
	8 d2p, d2ro, ero, et, e_cno, e_gr, e_pp, e_nuc, e_3al,
	9 f_tau, f17, gam, grav, gtmax=-1.d0, hh, lext, lp1, lp2, LOGg,
	1 lteff, mext, mgtmax, mtotp, mtot_ascii, nel, n13, pertmm, pext,
	2 ptext, o15, rext, rgtmax, rp1, rp2, ro_ext, teff, text, tp1,
	3 tp2, w_max, w_rotp, x0_ascii, y0_ascii

	INTEGER, PARAMETER, DIMENSION(2) ::  Ktest=(/ 3, 4 /)
	INTEGER, DIMENSION(8) :: values
	INTEGER, PARAMETER :: ntot_ecrit=50
	INTEGER, DIMENSION(1) :: igtmax
	INTEGER, SAVE :: l0, n_ecrit=0
	INTEGER :: dim_qsp, i, iglob, ihe4_ascii, ipn, itot, ivar,
	1 Krotp, Krot_ascii, i_cno, i_gr, i_pp, i_3a, j, k, knotb,
	2 knote, lq, m_chp, m_qsp, m_rotp, n, nadd, nchimp,
	3 nchim_ascii, nep, nrotp, ord_qsp, un23

	LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: convec
	LOGICAL, SAVE :: coca=.FALSE., coca_e=.FALSE., coca_pr=.FALSE.,
	1 cohe=.FALSE., cohe_e=.FALSE., cohe_pr=.FALSE.,
	2 coox=.FALSE., coox_e=.FALSE., coox_pr=.FALSE.,
	3 der, llist=.FALSE., list_pms=.FALSE., list_rep=.FALSE.,
	4 list_sort=.FALSE., list_zams=.FALSE., passe=.FALSE.,
	5 post=.FALSE., post_e=.FALSE., post_pr=.FALSE.,
	6 zams=.FALSE., zams_e=.FALSE., zams_pr=.FALSE.
	LOGICAL :: diffusionp, ecrHR, ecrit, en_massep, lim_rop, ok,
	1 rad_atm, rot_solidp, sort

	CHARACTER (len=1) :: oui
	CHARACTER (len=2) :: precisionp
	CHARACTER (len=3) :: text1, text2
	CHARACTER (len=4), ALLOCATABLE, DIMENSION(:) :: nom_elemp,
	1 nom_elem_ascii
	CHARACTER (len=4), SAVE :: nom_vwrot='Wrot'
	CHARACTER (len=4) :: number
	CHARACTER (len=5) :: zone
	CHARACTER (len=8) :: date
	CHARACTER (len=9), PARAMETER, DIMENSION(12) :: mois=(/
	1 'Janvier  ','Février  ','Mars     ','Avril    ',
	2 'Mai      ','Juin     ','Juillet  ','Aout     ',
	3 'Septembre','Octobre  ','Novembre ','Décembre ' /)
	CHARACTER (len=9), PARAMETER, DIMENSION(12) :: month=(/
	1 'January  ','February ','March    ','April    ',
	2 'May      ','June     ','July     ','August   ',
	3 'September','October  ','November ','December ' /)
	CHARACTER (len=10), ALLOCATABLE, SAVE, DIMENSION(:,:) :: tx_ioni
	CHARACTER (len=10) :: time
	CHARACTER (len=20), SAVE :: nom_save
	CHARACTER (len=20) :: nom_ctesp, nom_pertmp, nom_pertwp,
	1 nom_tdetaup, nom_atmp, nom_convp, nom_nucp, nom_nuc_cplp,
	2 nom_diffmp, nom_difftp, nom_diffwp, nom_etatp
	CHARACTER (len=50), DIMENSION(8) :: f_eosp, f_opap
	CHARACTER (len=50) :: chaine, nom_opap, nom_fich1, modelb
	CHARACTER (len=80) :: chain, chain_atm, chain_don, chain_rep, titre

	NAMELIST/nl_rlg/m_qs,m_ch,m_rot,m_tds,m_ptm,ordre,precix,
	1 precit,ro_test,psi0,d_grav,loc_zc,dtmax,dt0,ini0,n_atm,kipp,
	2 en_masse,ctel,ctep,ctem,cter,ctet,mvt_dis,dn_fixe,dpsi,
	3 mu_saha,ajuste,lisse,q0,l0,new_bv,fmin_abon,dlntc,iter_qs
	NAMELIST/nl_langue/langue
	NAMELIST/nl_blabla/baratine

c--------------------------------------------------------------------------

2000	FORMAT(8es10.3)

c recherche d'un éventuel fichier langue dans l'environnement pour fixer
c une langue différente du français
	INQUIRE(file='langue',exist=ok)
	IF(ok)THEN
	 OPEN(unit=30,form='formatted',status='old',file='langue')
	 READ(30,nl_langue) ; CLOSE(UNIT=30)
	ELSE
	 WRITE(*,30)
30	 FORMAT(/,'--------------------',/,
	1 'CESAM speaks a bit of english if you include in the working',/,
	2 'directory a file named ''langue'' with the statements :'
	3 ,/,'&NL_LANGUE',/,"langue='english'",/,'/',/,
	4 'cf. aide_mem2k, chapter Personnalisation',/,
	5 '-------------------')
	 langue='francais'
	ENDIF
c	langue='english'

c menu d'entrée
	B1: DO
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,29)
29	  FORMAT(/,'To stop : enter 0 then click on RETURN',/,
	1 'To continue the evolution of a stellar model further on,',/,
	2 'enter 1 then click on RETURN',/,
	3 'To start the evolution with a homogeneous ZAMS model,',/,
	4 'enter 2 then click on RETURN',/,
	5 'To start the evolution with a homogeneous PMS model,',/,
	6 'enter 3 then click on RETURN')
	 CASE DEFAULT
	  WRITE(*,28)
28	  FORMAT(/,'Pour arrêter : taper 0 puis RETURN',/,
	1 'Pour poursuivre une évolution : taper 1 puis RETURN',/,
	2 'Pour initialiser un modèle de ZAMS : taper 2 puis RETURN',/,
	3 'Pour initialiser un modèle de PMS : taper 3 puis RETURN')
	 END SELECT

c on entre un23 = 1, 2 ou 3 comme indiqué ci dessus
c	puis, suivant les cas on donné à un23 les valeurs suivantes:

c	un23 : 1 poursuite d'une évolution, modèle repris en binaire
c	       2 modèle initial en binaire, ZAMS
c              3 modèle initial en binaire, PMS
c	      -1 poursuite d'une évolution, modèle repris en ASCII
c             -2 modèle initial en ASCII, ZAMS
c             -3 modèle initial en ASCII, PMS

c	pms : modèle de PMS
c	zams : modèle de ZAMS
c	post : modèle après la ZAMS
c	cohe : modèle avec combustion de l'helium
c	coca : modèle avec combustion du carbone
c	coox : modèle avec combustion de l'oxygène

	READ*,un23

c gestion des modèles repris
	SELECT CASE(un23)
	CASE(0)
	 SELECT CASE(langue)
	 CASE('english')
	  PRINT*,'STOP'
	 CASE DEFAULT
	  PRINT*,'ARRET'
	 END SELECT
	 STOP

c poursuite d'un modèle
c rep_atm=.TRUE.: si existe, reprise atm. restituée en binaire
	CASE(1)
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1031)
1031	  FORMAT('Is the  model which evolution is to be continued',/,
	1'stored in a binary file? enter o for yes; enter n for no')
	 CASE DEFAULT
	  WRITE(*,31)
31	  FORMAT('Le modèle à poursuivre est-il donné en binaire ? o/n')
	 END SELECT

c on doit entrer o ou n
	 B4: DO
	  READ(*,'(a)')oui
	  SELECT CASE(oui)
	  CASE('n', 'o')
	   EXIT B4
	  CASE DEFAULT
	   SELECT CASE(langue)
	   CASE('english')
	    PRINT*,'enter o or n, you have enter : ',oui
	   CASE DEFAULT
	    PRINT*,'entrer o ou n, vous avez entré :  : ',oui
	   END SELECT
	  END SELECT
	 ENDDO B4
	 IF(oui /= 'n')THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1033)
1033	   FORMAT('Enter the name of the binary file for the model',/,
	1  'which evolution is to be continued.',/,
	2  'For instance: my_model_B.rep, my_model_B.hom, my_model_B.pms')
	  CASE DEFAULT
	   WRITE(*,33)
33	   FORMAT('nom du fichier binaire à poursuivre',/,
	1  'Ex.: mon_modele_B.rep, mon_modele_B.hom, mon_modele_B.pms')
	  END SELECT
	  READ*,nom_fich1
	  INQUIRE(file=TRIM(nom_fich1),exist=ok)
	  IF(ok)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1035)TRIM(nom_fich1)
1035	    FORMAT('CESAM will use the model : ',a)
	   CASE DEFAULT
	    WRITE(*,35)TRIM(nom_fich1)
35	    FORMAT('CESAM utilisera le modèle : ',a)
	   END SELECT
	   zams=.FALSE. ; post=.FALSE. ; cohe=.FALSE. ; coca=.FALSE.
	   coox=.FALSE. ; rep_atm=.TRUE. ;  EXIT B1
	  ELSE
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1036)TRIM(nom_fich1)
1036	    FORMAT('This input file is unknown : ',a,/,'STOP')
	   CASE DEFAULT
	    WRITE(*,36)TRIM(nom_fich1)
36	    FORMAT('Fichier d''entrée inconnu : ',a/,'ARRET')
	   END SELECT
	   STOP
	  ENDIF           !ok
	 ELSE
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1081)
1081	   FORMAT('Name of the ASCII file for the model which ',/,
	1 ' evolution is to be continued; for instance: my_model.ascii')
	  CASE DEFAULT
	   WRITE(*,81)
81	   FORMAT('Nom du fichier ASCII du modèle à poursuivre',/,
	1 ' Exemple: mon_modele.ascii')
	  END SELECT
	  READ*,nom_fich1
	  INQUIRE(file=TRIM(nom_fich1),exist=ok)
	  IF(ok)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1082)TRIM(nom_fich1)
1082	    FORMAT('CESAM use the model : ',a)
	   CASE DEFAULT
	    WRITE(*,82)TRIM(nom_fich1)
82	    FORMAT('CESAM utilise le modèle : ',a)
	   END SELECT
	   rep_atm=.FALSE.
	   zams=.FALSE. ; post=.FALSE. ; coox=.FALSE.
	   cohe=.FALSE. ; coca=.FALSE. ; un23=-1 ; EXIT B1
	  ELSE
	   SELECT CASE(langue)
	   CASE('english')
	     WRITE(*,1085)TRIM(nom_fich1)
1085	     FORMAT('STOP, this ASCII file is unknown : ',a)
	    CASE DEFAULT
	     WRITE(*,85)TRIM(nom_fich1)
85	     FORMAT('ARRET, fichier ASCII inconnu : ',a)
	    END SELECT
	    STOP
	   ENDIF
	  ENDIF

c modèle initial de ZAMS
	 CASE(2)
	  zams=.TRUE. ; post=.FALSE. ; cohe=.FALSE.
	  coca=.FALSE. ; coox=.FALSE. ; rep_atm=.FALSE.
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1084)
1084	   FORMAT('Is the input ZAMS model stored in a binary file?',/,
	1  'enter o for yes; n for no')
	  CASE DEFAULT
	   WRITE(*,84)
84	   FORMAT('Le modèle initial de ZAMS est-il en binaire ? o/n')
	  END SELECT
	  READ*,oui
	  IF(oui /= 'n')THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1086)
1086	    FORMAT('Enter the name of the binary file for the initial',/,
	1   'model. For instance : my_model_B.hom')
	   CASE DEFAULT
	    WRITE(*,86)
86	    FORMAT('Entrer le nom du fichier binaire du modèle initial',/,
	1   'Exemples: mon_modele_B.hom')
	   END SELECT
	   READ*,nom_fich1 ; INQUIRE(file=TRIM(nom_fich1),exist=ok)
	   IF(ok)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1087)TRIM(nom_fich1)
1087	     FORMAT('CESAM will use the following binary file : ',a)
	    CASE DEFAULT
	     WRITE(*,87)TRIM(nom_fich1)
87	     FORMAT('CESAM utilise le modèle binaire : ',a)
	    END SELECT
	    un23=2  !fichier en binaire (intruction redondante)
	   EXIT B1
	  ELSE
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1088)TRIM(nom_fich1)
1088	    FORMAT('STOP unknown input file : ',a)
	   CASE DEFAULT
	    WRITE(*,88)TRIM(nom_fich1)
88	    FORMAT('ARRET fichier d''entrée inconnu : ',a)
	   END SELECT
	   STOP
	  ENDIF    !ok
	 ELSE      !oui
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1089)
1089	   FORMAT('Enter the name of the ASCII file for the initial ',/,
	1  'model; for instance: m010.zams, m020.zams, m050.zams, m150.zams')
	  CASE DEFAULT
	   WRITE(*,89)
89	   FORMAT('Entrer le nom du fichier ASCII du modèle initial',/,
	1  'Exemples: m010.zams, m020.zams, m050.zams, m150.zams')
	  END SELECT
	  READ*,nom_fich1
	  INQUIRE(file=TRIM(nom_fich1),exist=ok)
	  IF(ok)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1090)TRIM(nom_fich1)
1090	    FORMAT('CESAM will use the model : ',a)
	   CASE DEFAULT
	    WRITE(*,90)TRIM(nom_fich1)
90	    FORMAT('CESAM utilisera le modèle : ',a)
	   END SELECT
	   PRINT*; un23=-2       !fichier en ASCII
	   EXIT B1
	  ELSE
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1091)TRIM(nom_fich1)
1091	    FORMAT('STOP, unknown input file : ',a)
	   CASE DEFAULT
	    WRITE(*,91)TRIM(nom_fich1)
91	    FORMAT('ARRET, fichier d''entrée inconnu : ',a)
	   END SELECT
	   STOP
	  ENDIF    !ok
	 ENDIF           !oui

c PMS
	 CASE(3)
	  zams=.FALSE. ; post=.FALSE.; cohe=.FALSE.
	  coca=.FALSE. ; coox=.FALSE. ; rep_atm=.FALSE.
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1092)
1092	   FORMAT('Is the initial PMS model stored in a binary file?',/,
	1  'enter o for yes, enter n for no')
	  CASE DEFAULT
	   WRITE(*,92)
92	   FORMAT('Le modèle initial de PMS est-il donné en binaire ? o/n')
	  END SELECT
	  READ*,oui
	  IF(oui == 'o')THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1093)
1093	    FORMAT('Enter the name of the binary file for the inital',/,
	1   'model; for instance : my_model_B.pms')
	   CASE DEFAULT
	    WRITE(*,93)
93	    FORMAT('Entrer le nom du fichier binaire du modèle initial',/,
	1   'Exemples: mon_modele_B.pms')
	   END SELECT
	   READ*,nom_fich1 ; INQUIRE(file=TRIM(nom_fich1),exist=ok)
	   IF(ok)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1094)TRIM(nom_fich1)
1094	     FORMAT('CESAM uses the initial PMS model : ',a)
	    CASE DEFAULT
	     WRITE(*,94)TRIM(nom_fich1)
94	     FORMAT('CESAM utilise le modèle initial de PMS : ',a)
	    END SELECT
	    PRINT*
	    un23=3 !fichier en binaire et, pour la PMS, pas de modèle homo.
	    EXIT B1
	   ELSE
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1095)TRIM(nom_fich1)
1095	     FORMAT('STOP, unknown input file : ',a)
	    CASE DEFAULT
	     WRITE(*,95)TRIM(nom_fich1)
95	     FORMAT('ARRET, fichier d''entrée inconnu : ',a)
	    END SELECT
	    STOP
	   ENDIF    !ok
	  ELSE      !oui
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1096)
1096	    FORMAT('CESAM uses an input PMS model in an ASCII file',/,
	1   'enter the name of the ASCII file for the initial PMS model.',
	2   /,'For instance: 2d-2.pms(Tc=0.2MK), 5d-4.pms(Tc=0.5MK)',/,
	3   '8d-5.pms(Tc=1.0MK), 4d-2.pms (for M*>10Msol)')
	   CASE DEFAULT
	    WRITE(*,96)
96	    FORMAT('CESAM utilise un modèle initial de PMS en ASCII',/,
	1   'entrer le nom du fichier ASCII du modèle PMS initial',/,
	2   'Exemples: 2d-2.pms(Tc=0.2MK), 5d-4.pms(Tc=0.5MK)',/,
	3   ' 8d-5.pms(Tc=1.0MK), 4d-2.pms (pour M*>10Msol)')
	   END SELECT
	   READ*,nom_fich1; INQUIRE(file=TRIM(nom_fich1),exist=ok)
	   IF(ok)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1097)TRIM(nom_fich1)
1097	     FORMAT('CESAM uses the model in the ASCII file : ',a)
	    CASE DEFAULT
	     WRITE(*,97)TRIM(nom_fich1)
97	     FORMAT('CESAM utilise le modèle ASCII : ',a)
	    END SELECT
	    un23=-3  !fichier en ASCII et, pour la PMS, pas de modèle homo.
	    EXIT B1
	   ELSE
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1098)TRIM(nom_fich1)
1098	     FORMAT('STOP, unknown input file : ',a)
	    CASE DEFAULT
	     WRITE(*,98)TRIM(nom_fich1)
98	     FORMAT('ARRET, fichier d''entrée inconnu : ',a)
	    END SELECT
	    STOP
	   ENDIF    !ok
	  ENDIF           !oui
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1099)
1099	   FORMAT('Error: enter 0, 1, 2 or 3; it was entered : ',i4)
	  CASE DEFAULT
	   WRITE(*,99)un23
99	   FORMAT('ERREUR, entrer 0, 1, 2 ou 3, vous avez tapé : ',i4)
	  END SELECT
	  CYCLE B1
	 END SELECT             !menu d'entrée
	ENDDO B1

c identificateur des fichiers des données, du point de reprise,
c des résultats, du modèle d'âge 0, des oscillations
	SELECT CASE(langue)
	CASE('english')
	 WRITE(*,1100)
1100	 FORMAT('enter the generic name of the model : ',/,
	1 'for instance : my_model')
	CASE DEFAULT
	 WRITE(*,100)
100	 FORMAT('entrer l''identificateur du modèle : ',/,
	1 'Exemple : mon_modele')
	END SELECT
	READ*,nom_fich2
	SELECT CASE(langue)
	CASE('english')
	 WRITE(*,1103)TRIM(nom_fich2)
1103	 FORMAT('generic name of the model files : ',a)
	CASE DEFAULT
	 WRITE(*,103)TRIM(nom_fich2)
103	 FORMAT('identificateur des fichiers du modèle : ',a)
	END SELECT

c ouverture du fichier du listing mon_modele.lis
	OPEN(unit=2,form='formatted',status='unknown',
	1 file=TRIM(nom_fich2)//'.lis') !fichiers du listing
	SELECT CASE(langue)
	CASE('english')
	 WRITE(2,1001)version ; WRITE(*,1001)version
1001	 FORMAT(/,t10,
	1'***************************************************************',
	2 //,t11,
	3 'MODEL OF INTERNAL STRUCTURE computed with CESAM2k version : ',a,
	4 //,t10,
	5'***************************************************************',
	6/)
	CASE DEFAULT
	 WRITE(2,1)version ; WRITE(*,1)version
1	 FORMAT(/,t10,
	1 '***************************************************************',
	2 //,t11,
	3 'MODELE DE STRUCTURE INTERNE calculé par CESAM2k version : ',a,
	4 //,t10,
	5 '***************************************************************',
	6/)
	END SELECT
	CALL date_and_time(date,time,zone,values)
	SELECT CASE(langue)
	CASE('english')
	 titre='The calculation started on the : '//date(7:8)//' '
	1 //TRIM(month(values(2)))//' '//date(1:4)//' at '
	2 //time(1:2)//'h'//time(3:4)
	CASE DEFAULT
	 titre='Début du calcul le : '//date(7:8)//' '
	1 //TRIM(mois(values(2)))//' '//date(1:4)//' à '
	2 //time(1:2)//'h'//time(3:4)
 	END SELECT
	WRITE(*,3)TRIM(titre) ; WRITE(2,3)TRIM(titre)
3	FORMAT(t10,a,/)

c lecture du fichier des données mon_modele.don
	CALL lit_nl(wrot)
	
c recherche d'un éventuel fichier blabla dans l'environnement qui permet de
c détourner une partie des sorties sur le miniteur vers des fichiers ASCII
c cf. NOTICE § Limitation des sorties
	INQUIRE(file='blabla',exist=ok)
	IF(ok)THEN
	 OPEN(unit=30,form='formatted',status='old',file='blabla')
	 READ(30,nl_blabla) ; CLOSE(UNIT=30)
	 IF(.NOT.baratine)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1301)TRIM(nom_fich2)//'_static',TRIM(nom_fich2)//'_atm',
	1  TRIM(nom_fich2)//'_evol'
1301	   FORMAT(/,'Part of the on line outputs rerouted towards the files :',/,
	1  a,/,a,/,a,/)	 
	  CASE DEFAULT
	   WRITE(*,301)TRIM(nom_fich2)//'_static',TRIM(nom_fich2)//'_atm',
	1  TRIM(nom_fich2)//'_evol'
301	   FORMAT(/,'Détournement de sorties on line vers les fichiers :',/,
	1  a,/,a,/,a,/)
	  END SELECT
	 ENDIF
	ENDIF

c écriture des constantes initialisées dans lit_nl
	CALL print_ctes(6) ; CALL print_ctes(2)
	SELECT CASE(langue)
	CASE('english')
	 WRITE(*,1005) ; WRITE(2,1005)
1005	 FORMAT('Input physics :',/)
	CASE DEFAULT
	 WRITE(*,5) ; WRITE(2,5)
5	 FORMAT(t14,'Physique utilisée :',/)
	END SELECT
	WRITE(*,2)nom_etat, nom_opa, nom_conv, nom_nuc, nom_nuc_cpl,
	1 nom_abon, nom_atm, nom_tdetau, nom_pertm, nom_pertw, nom_diffm,
	2 nom_difft, nom_diffw, nom_ctes
	WRITE(2,2)nom_etat, nom_opa, nom_conv, nom_nuc, nom_nuc_cpl,
	1 nom_abon, nom_atm, nom_tdetau, nom_pertm, nom_pertw, nom_diffm,
	2 nom_difft, nom_diffw, nom_ctes
2	FORMAT('EOS : ',a,t30,'Opa. : ',a,t60,'conv. : ',a,/,
	1 'réac. nuc. : ',a,t30,'compil. nuc. : ',a,t60,'abon_ini. : ',a,/,
	2 'atmos.: ',a,t30,'T(tau): ',a,t60,'perte m. : ',a,/,
	3 'perte w. : ',a,t30,'diff. micro. : ',a,t60,'diff. turb. : ',a,/
	4 'diff. rot. : ',a,t30,'ctes. : ',a,/)

c quelques constantes physiques
	cte1=gmsol/rsol**2
	cte2=mtot**2*gmsol/lsol/rsol*msol*2.d0/secon6
	cte3=4.d0*pi*rsol**3/msol
	cte4=2.d0/3.d0*rsol**3/g/msol
	cte5=-4.d0*pi*g/3.d0*rsol**2
	cte6=0.5d0/g/pi
	cte8=SQRT(SQRT(lsol/(4.d0*pi*rsol**2*sigma)))	!pour Teff
	cte9=gmsol/rsol/rsol

c initialisation initiale des constantes de réglage
c correspondant à la précision par défaut : pr "réaliste"
c elles seront adaptées pour chaque type de précision
	m_qs=2		!ordre des splines pour les variables quasi-stat.
	m_ch=2		!ordre des splines pour interpolation de la comp.chim.
	m_rot=2		!ordre des splines pour la rotation
	m_tds=2		!ordre des splines pour interpolation du TdS
	m_ptm=2		!ordre des spl. pour inter. de la masse (perte de masse)
	ordre=2		!ordre du schéma d'int. des réac.nuc. avec rk_imps
	precix=1.d-3	!précision sur les itérations équilibre quasi-statique
	precit=0.15d0	!var. rela. max. pour int. temp. de la comp. chim.
	ro_test=0.1d0	!test de variation du TdS si ro > ro_test
	psi0=0.08d0	!constante de répartition à assurer
	d_grav=0.5d0	!variation maximale du TdS
	loc_zc=1.d-3	!précision de la localisation des limites ZR/ZC
	dtmax=200.d0	!pas temporel maximum
	dlntc=0.075d0	!variation relative max de Tcentrale sur un pas temporel
	IF(mtot <= 1.1d0)THEN	!pas temporel initial pour évolution 
	 dt0=10.d0		!à partir d'un modèle de ZAMS
	ELSE
	 dt0=1.d0
	ENDIF	
	ini0=4		!nb. iter. N-R avec réestim. de comp.chim. lim. ZR/ZC
	n_atm=75	!nb. couches pour l'atmosphère restituée
	kipp=.FALSE.	!approximation de Kippenhahn ou TdS=dE+PdV
	en_masse=.TRUE.	!variables quasi_stat. en lagrangien
	ctel=0.d0	!constante de répartition pour la luminosité
	ctep=-1.d0	!constante de répartition pour la pression
	ctem=15.d0	!constante de répartition pour la masse
	cter=0.d0	!constante de répartition pour le rayon
	ctet=-1.d0	!constante de répartition pour la température
	mvt_dis=.FALSE.	!ajustement de comp.chim. suivant mvt. des discon.
	dn_fixe=0.05d0	!taux pour modif de grille fixe en comp.chim
	dpsi=0.05d0	!variation max pour modification de n_qs
	mu_saha=.TRUE.	!utilisation de saha pour le calcul de mu
	ajuste=.FALSE.	!on ajustera T, He core etc...
	lisse=.FALSE.	!comp chim et rot lissés par contour
	q0=0.05d0!un point à q0 du centre dans les fichiers ASCII de sortie
	l0=0	!nombre de points ajoutés au voisinage des limites ZR/ZC
	new_bv=.TRUE.	!Vaissala utilisant phi
	fmin_abon=0.05d0	!facteur reliant ab_min et abon_ini
	iter_qs=(/ 0, 0, 0, 0, 0, 0, 7 /) !contrôle des variables q-stat
	n_max=MAX(n_max,1000)	!nb. max de couches
	l_demi=.TRUE.	!point milieu pour la comp_chim	un pas sur 2
	n_init=300	!nombre de couches du modèle initial

c initialisation des paramètres numériques suivant le type de précision
	SELECT CASE(langue)
	CASE('english')
         WRITE(*,1218) ; WRITE(2,1218)
1218	 FORMAT('Numerical parameters which are used : ')
	CASE DEFAULT
         WRITE(*,218) ; WRITE(2,218)
218      FORMAT('Paramètres numériques utilisés : ')
	END SELECT

	SELECT CASE(precision)

c défaut avec ajustement de T, He core etc...	 
	CASE('aj')	 
	 ajuste=.TRUE.
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1134) ; WRITE(2,1134)
1134	  FORMAT('model computed with a realist precision and adjustments')
	 CASE DEFAULT
	  WRITE(*,134) ; WRITE(2,134)
134	  FORMAT('modèle calculé avec précision réaliste et ajustements')
	 END SELECT

c pour stades avancés 'long run' en eulérien
	CASE('ar')	
	 en_masse=.FALSE.
	 psi0=0.2d0
	 n_init=1500 	
	 mu_saha=.FALSE.
	 ajuste=.FALSE.
	 lisse=.TRUE.
	 d_grav=5.d0	!elimination du test
	 q0=0.d0!un point à q0 du centre dans les fichiers ASCII de sortie
	 l0=0	!nombre de points ajoutés au voisinage des limites ZR/ZC
	 dtmin=1.d-10
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1135) ; WRITE(2,1135)
1135	  FORMAT('The model is computed with the eulerian long run parameters')
	 CASE DEFAULT
	  WRITE(*,135) ; WRITE(2,135)
135	  FORMAT('Modèle calculé avec les paramètres euleriens stades avancés')
	 END SELECT
	 
c pour stades avancés 'long run'
	CASE('av')	 	
	 mu_saha=.FALSE.
	 ajuste=.FALSE.
	 lisse=.TRUE.
	 d_grav=50.d0	!elimination du test
	 q0=0.d0!un point à q0 du centre dans les fichiers ASCII de sortie
	 l0=0	!nombre de points ajoutés au voisinage des limites ZR/ZC
	 dtmin=1.d-10	 
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1133) ; WRITE(2,1133)
1133	  FORMAT('The model is computed with the long run parameters')
	 CASE DEFAULT
	  WRITE(*,133) ; WRITE(2,133)
133	  FORMAT('Modèle calculé avec les paramètres pour les stades avancés')
	 END SELECT

c précision réaliste en eulérien c'est le "défaut"
	CASE('er')
	 m_qs=1
	 mrot=1	
	 en_masse=.FALSE.
	 ctem=0.d0
	 cter=15.0d0
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1316) ; WRITE(2,1316)
1316	  FORMAT('the model is computed with a realistic Eulerian precision')
	 CASE DEFAULT
	  WRITE(*,316) ; WRITE(2,316)
316	  FORMAT('Modèle calculé en précision eulérien réaliste')
	 END SELECT

c précision CoRoT
	CASE('co')
	 precix=1.d-4
	 precit=0.05d0	 
	 psi0=0.06d0		
	 d_grav=0.5d0	!variation maximale du TdS 
	 loc_zc=1.d-4
	 dtmax=MIN(50.d0,dtmax)	!pas temporel maximum	 
	 ini0=5	
	 n_atm=100
	 q0=0.01d0
	 l0=5
	 lisse=.FALSE.
	 fmin_abon=0.01d0
	 dlntc=0.05d0	!var. relative max de Tcentrale sur un pas temporel 
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1230) ; WRITE(2,1230)
1230	  FORMAT('the model is computed with a CoROT precision')
	 CASE DEFAULT
	  WRITE(*,230) ; WRITE(2,230)
230	  FORMAT('modèle calculé avec la précision CoROT')
	 END SELECT
	 
c Hyper précision
	CASE('hp')
	 m_qs=3
	 m_ch=4
	 m_rot=3
	 ordre=4
	 precix=1.d-6
	 precit=0.01d0
	 ro_test=0.05d0
	 psi0=0.05d0
	 d_grav=0.05d0
	 loc_zc=1.d-6
	 dtmax=MIN(50.d0,dtmax)
	 ini0=6
	 n_atm=150
	 n_max=MAX(n_max,2500)	 
	 q0=0.01d0
	 l0=7
	 fmin_abon=0.01d0	 
	 dlntc=0.03d0	!var. relative max de Tcentrale sur un pas temporel
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1217) ; WRITE(2,1217)
1217	  FORMAT('the model is computed with a hyper precision')
	 CASE DEFAULT
	  WRITE(*,217) ; WRITE(2,217)
217	  FORMAT('modèle calculé en hyper précision')
	 END SELECT

c pour petites masses très évoluées
	CASE('lm')
	 precit=0.2d0
	 dtmax=MIN(300.d0,dtmax)
	 ini0=5
	 kipp=.TRUE.
	 dn_fixe=0.1d0
	 mu_saha=.FALSE.
	 q0=0.d0!un point à q0 du centre dans les fichiers ASCII de sortie
	 l0=0	!nombre de points ajoutés au voisinage des limites ZR/ZC
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1306) ; WRITE(2,1306)
1306	  FORMAT('for highly evolved models with small masses')
	 CASE DEFAULT
	  WRITE(*,306) ; WRITE(2,306)
306	  FORMAT('pour petites masses très evoluées')
	 END SELECT

c nombre de couches maximal
	CASE('mx')
	 precix=1.d-5
	 precit=0.02d0
	 psi0=1.d-8	
	 d_grav=0.5d0	!variation maximale du TdS
	 loc_zc=1.d-5
	 dtmax=MIN(50.d0,dtmax)
	 ini0=5
	 n_atm=100
	 lisse=.FALSE.	!la composition chimique est lissée par contour
	 q0=0.01d0
	 l0=4
	 fmin_abon=0.01d0
	 dlntc=0.05d0	!variation relative max de Tcentrale sur un pas temporel	 
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1310)n_max ; WRITE(2,1310)n_max
1310	  FORMAT('the model is computed with a maximum shell number=',i5)
	 CASE DEFAULT
	  WRITE(*,310)n_max ; WRITE(2,310)n_max
310	  FORMAT('modèle calculé avec nb. couches max.=',i5)
	 END SELECT

c précision normale en variables lagrangiennes
	CASE('np')
	 m_qs=1
	 m_rot=1
	 ordre=1
	 precix=5.d-3
	 precit=0.3d0
	 psi0=0.1d0	 
	 loc_zc=5.d-3
	 dtmax=MIN(300.d0,dtmax)
	 ini0=2
	 n_atm=50
	 kipp=.TRUE.
	 mvt_dis=.FALSE. !ajustement comp.chim. suivant mvt. des discon	 
	 mu_saha=.FALSE.
	 ajuste=.FALSE.
	 q0=0.d0!un point à q0 du centre dans les fichiers ASCII de sortie
	 l0=0	!nombre de points ajoutés au voisinage des limites ZR/ZC
	 dlntc=0.1d0
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1216) ; WRITE(2,1216)
1216	  FORMAT('model computed with a normal Lagragian precision')
	 CASE DEFAULT
	  WRITE(*,216) ; WRITE(2,216)
216	  FORMAT('modèle Lagrangien calculé en précision normale')
	 END SELECT
	 
c précision normale en variables eulériennes
	CASE('nr')
	 m_qs=1
	 m_rot=1
	 ordre=1
	 precix=5.d-3
	 precit=0.3d0
	 psi0=0.1d0
	 en_masse=.FALSE.
	 ctem=0.d0
	 cter=15.0d0
	 loc_zc=5.d-3
	 dtmax=MIN(300.d0,dtmax)
	 ini0=3
	 n_atm=50
	 mvt_dis=.FALSE.
	 mu_saha=.FALSE.
	 ajuste=.FALSE.
	 kipp=.TRUE.
	 q0=0.d0!un point à q0 du centre dans les fichiers ASCII de sortie
	 l0=0	!nombre de points ajoutés au voisinage des limites ZR/ZC
	 dlntc=0.1d0
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1224) ; WRITE(2,1224)
1224	  FORMAT('the model is computed with a normal Eulerian precision')
	 CASE DEFAULT
	  WRITE(*,224) ; WRITE(2,224)
224	  FORMAT('modèle calculé en précision normale eulérien')
	 END SELECT	 
	 
c précision réaliste c'est le "défaut"
	CASE('pr')
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1130) ; WRITE(2,1130)
1130	  FORMAT('the model is computed with a realist precision')
	 CASE DEFAULT
	  WRITE(*,130) ; WRITE(2,130)
130	  FORMAT('modèle calculé avec une précision réaliste')
	 END SELECT
	 
c avec réglages externes
	CASE('rg')
	 chain=TRIM(nom_fich2)//'.rg'
	 INQUIRE(file=TRIM(chain),exist=ok)
	 IF(ok)THEN
	  OPEN(unit=3,form='formatted',status='old',delim='apostrophe',
	1  file=TRIM(chain))
	 ELSE
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1024)TRIM(chain) ; WRITE(2,1024)TRIM(chain)
1024	   FORMAT('The parameter setting file : ',a,/,
	1 'is unknown, research and use of the file named reglages')
	  CASE DEFAULT
	   WRITE(*,24)TRIM(chain) ; WRITE(2,24)TRIM(chain)
24	   FORMAT('le fichier des réglages : ',a,/,
	1  'non trouvé, recherche et utilisation du fichier reglages')
	  END SELECT
	  chain='reglages'
	  INQUIRE(file=TRIM(chain),exist=ok)
	  IF(ok)THEN
	   OPEN(unit=3,form='formatted',status='old',delim='apostrophe',
	1  file=TRIM(chain))
	  ELSE
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1027)chain ; WRITE(2,1027)chain ; STOP
1027	    FORMAT('STOP, The parameter setting file : ',a,' is unknown',/,
	1   'in the directory where the calculations are performed')
	   CASE DEFAULT
	    WRITE(*,27)chain ; WRITE(2,27)chain ; STOP
27	    FORMAT('ARRET, le fichier des réglages : ',a,' est inconnu',/,
	1   'dans le directory où sont effectués les calculs')
	   END SELECT
	  ENDIF
	 ENDIF
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1311)chain ; WRITE(2,1311)chain
1311	  FORMAT('The model is computed with a choice of parameters set',/,
	1 'in the file : ',a,/)
	 CASE DEFAULT
	  WRITE(*,311)chain ; WRITE(2,311)chain
311	  FORMAT('modèle utilisant les réglages externes du fichier : ',a,/)
	 END SELECT
	 READ(3,nl_rlg) ; CLOSE(UNIT=3)
	 WRITE(*,nl_rlg) ; WRITE(2,nl_rlg)	 

c Solar accuracy
	CASE('sa')
	 precix=1.d-5
	 precit=0.02d0	 
	 psi0=0.06d0
	 loc_zc=1.d-5
	 dtmax=MIN(50.d0,dtmax)
	 ini0=5
	 n_atm=100
	 lisse=.FALSE.	 
	 q0=0.01d0
	 l0=5
	 fmin_abon=0.01d0
	 new_bv=.FALSE.
	 dlntc=0.05d0	!var. relative max de Tcentrale sur un pas temporel
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1300) ; WRITE(2,1300)
1300	  FORMAT('the model is computed with the solar precision')
	 CASE DEFAULT
	  WRITE(*,300) ; WRITE(2,300)
300	  FORMAT('modèle calculé avec la précision solaire')
	 END SELECT

c super précision lagrangien
	CASE('sp')
	 precix=1.d-4
	 precit=0.05d0
	 psi0=0.06d0	
	 d_grav=0.3d0	!variation maximale du TdS	 
	 loc_zc=1.d-4
	 dtmax=MIN(50.d0,dtmax)
	 ini0=5
	 n_atm=100
	 fmin_abon=0.01d0
	 dlntc=0.05d0	!var. relative max de Tcentrale sur un pas temporel 
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1215) ; WRITE(2,1215)
1215	  FORMAT('The model is computed with a high Lagragian precision')
	 CASE DEFAULT
	  WRITE(*,215) ; WRITE(2,215)
215	  FORMAT('Modèle calculé en super précision Lagragien')
	 END SELECT
	 
c super précision eulérien
	CASE('sr')
	 m_qs=1
	 mrot=1	
	 en_masse=.FALSE.
	 ctem=0.d0
	 cter=15.0d0
	 precix=1.d-4
	 precit=0.05d0
	 psi0=0.06d0	
	 d_grav=0.5d0	!variation maximale du TdS	 
	 loc_zc=1.d-4
	 dtmax=MIN(50.d0,dtmax)
	 ini0=5
	 n_atm=100
	 fmin_abon=0.01d0
	 dlntc=0.05d0	!var. relative max de Tcentrale sur un pas temporel
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1315) ; WRITE(2,1315)
1315	  FORMAT('the model is computed with a high Eulerian precision')
	 CASE DEFAULT
	  WRITE(*,315) ; WRITE(2,315)
315	  FORMAT('Modèle calculé en super précision eulérien')
	 END SELECT	 

c problème
	CASE DEFAULT
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1312)	; WRITE(2,1312)	; STOP
1312	  FORMAT('STOP, this class of  precision is unknown : ',a,/,
	1 'defined classes of precision : sp, sa, np, co, lm, mx, rg, pr, av')
	 CASE DEFAULT
	  WRITE(*,312)precision ; WRITE(2,312)precision ; STOP
312	  FORMAT('ARRET, type de précision non prévu : ',a,/,
	1 'types de précision prévus : sp, sa, np, co, lm, mx, rg, pr, av')
	 END SELECT
	END SELECT  !type de calcul

c avec rotation, on lisse les discontinuités
	 SELECT CASE (Krot)
	 CASE(3,4)
	  IF(precision /= 'rg')THEN
	   precix=1.d-3
	   loc_zc=1.d-3	
	   mu_saha=.FALSE.
	   ajuste=.FALSE.
	   ord_rot=m_rot+r_qs
c	   lisse=.TRUE.
c	   d_grav=50.d0
c	   Kipp=.TRUE.
c	   l_demi=.FALSE.
	   grille_fixe=.TRUE.
	   ini0=3
	  ENDIF
	 CASE(5)
	  ord_rot=m_rot+1	    
	 END SELECT
c	 PRINT*,Krot,dtmax ; PAUSE'dtmax'

c ne : nombre d'inconnues du modèle quasi-statique
c Ipg : indice de ln Pgaz avec pression turbulente (1 sans Pturb)
c pturb=.TRUE. on tient compte de la pression turbulente
c avec pression turbulente 7 inconnues
c sans pression turbulente 6 inconnues, Ptot=Pgaz
	der=cpturb < 0.d0 !on tient compte de dln Pgaz/dln Ptot
	pturb=ABS(cpturb) /= 0.d0
	IF(pturb)THEN
	 ne=7 ; Ipg=7
	ELSE
	 ne=6 ; Ipg=1
	ENDIF

	SELECT CASE (Krot)
	CASE(3,4,5)
	 WRITE(*,299)m_qs,m_ch,m_rot,ordre,precix,loc_zc,precit,
	1 ro_test,psi0,d_grav,mvt_dis,ctel,ctem,ctep,cter,ctet,kipp,ajuste,
	2 lisse,new_bv,fmin_abon
	 WRITE(2,299)m_qs,m_ch,m_rot,ordre,precix,loc_zc,precit,
	1 ro_test,psi0,d_grav,mvt_dis,ctel,ctem,ctep,cter,ctet,kipp,ajuste,
	2 lisse,new_bv,fmin_abon
299	 FORMAT('m_qs=',i2,', m_ch=',i2,', m_rot=',i2,
	1 ', ordre=',i2,', precix=',es8.1,', loc_zc=',es8.1,/,
	2 'precit=',es8.1,', ro_test=',es8.1,', psi0=',es8.1,
	3 ', d_grav=',es8.1,', mvt_dis=',l2,/,'ctel=',es10.3,
	4 ', ctem=',es10.3,', ctep=',es10.3,', cter=',es10.3,
	5 ', ctet=',es10.3,/,'kipp=',l2,', ajuste=',l2,
	5  ', lisse=',l2,', new_bv=',l2,', fmin_abon=',es10.3,/)
	CASE DEFAULT
	 WRITE(*,305)m_qs,m_ch,ordre,precix,loc_zc,precit,ro_test,psi0,
	1 d_grav,mvt_dis,ctel,ctem,ctep,cter,ctet,kipp,ajuste,
	2 lisse,new_bv,fmin_abon
	 WRITE(2,305)m_qs,m_ch,ordre,precix,loc_zc,precit,ro_test,psi0,
	1 d_grav,mvt_dis,ctel,ctem,ctep,cter,ctet,kipp,ajuste,
	2 lisse,new_bv,fmin_abon
305	 FORMAT('m_qs=',i2,', m_ch=',i2,', ordre=',i2,', precix=',es8.1,
	1 ', loc_zc=',es8.1,/,'precit=',es8.1,', ro_test=',es8.1,
	2 ', psi0=',es8.1,', d_grav=',es8.1,', mvt_dis=',l2,/,
	3 'ctel=',es10.3,', ctem=',es10.3,', ctep=',es10.3,
	4 ', cter=',es10.3,', ctet=',es10.3,/,'kipp=',l2,', ajuste=',l2,
	5  ', lisse=',l2,', new_bv=',l2,', fmin_abon=',es10.3,/)
	END SELECT

	IF(en_masse)THEN
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1220) ; WRITE(2,1220)
1220	  FORMAT('use of Lagrangean variables')
	 CASE DEFAULT
	  WRITE(*,220) ; WRITE(2,220)
220	  FORMAT('utilisation des variables lagrangiennes')
	 END SELECT
	ELSE
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1221) ; WRITE(2,1221)
1221	  FORMAT('use of Eulerian variables')
	 CASE DEFAULT
	  WRITE(*,221) ; WRITE(2,221)
221	  FORMAT('utilisation des variables eulériennes')
	 END SELECT
	ENDIF
	IF(grille_fixe)THEN
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1222) ; WRITE(2,1222)
1222	  FORMAT('Rigid Lagrangean grid for the chemical composition')
	 CASE DEFAULT
	  WRITE(*,222) ; WRITE(2,222)
222	  FORMAT('grille fixe lagrangienne pour la composition chimique')
	 END SELECT
	ELSE
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1223) ; WRITE(2,1223)
1223	  FORMAT('Adaptative Lagrangean grid for the chemical composition')
	 CASE DEFAULT
	  WRITE(*,223) ; WRITE(2,223)
223	  FORMAT('grille variable lagrangienne pour la comp. chim.')
	 END SELECT
	ENDIF

c initialisations diverses et écritures
	IF(diffusion)THEN
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1240) ; WRITE(2,1240)
1240	  FORMAT('Diffusion of the chemical composition in radiative (ZR)',
	1 /,'and convective zones (ZC)')
	 CASE DEFAULT
	  WRITE(*,240) ; WRITE(2,240)
240	  FORMAT('diffusion de la comp. chim. dans ZR et ZC')
	 END SELECT

	ELSE
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1241) ; WRITE(2,1241)
1241	  FORMAT('model without diffusion')
	 CASE DEFAULT
	  WRITE(*,241) ; WRITE(2,241)
241	  FORMAT('modèle sans diffusion')
 	 END SELECT
	ENDIF
	IF(mdot == 0.d0)THEN
 	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1242) ; WRITE(2,1242)
1242	  FORMAT('model without mass loss')
	 CASE DEFAULT
	  WRITE(*,242) ; WRITE(2,242)
242	  FORMAT('modèle sans perte de masse')
	 END SELECT
	ELSEIF(mdot < 0.d0)THEN
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1212)mdot ; WRITE(2,1212)mdot
1212	  FORMAT('model with mass loss :',es10.3,'Msol/an')
	 CASE DEFAULT
	  WRITE(*,212)mdot ; WRITE(2,212)mdot
212	  FORMAT('model with mass gain :',es10.3,'Msol/an')
	 END SELECT
	ELSE
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1211)ABS(mdot) ; WRITE(2,1211)ABS(mdot)
1211	  FORMAT('modèle avec gain de masse :',es10.3,'Msol/an')
 	 CASE DEFAULT
	  WRITE(*,211)ABS(mdot) ; WRITE(2,211)ABS(mdot)
211	  FORMAT('modèle avec gain de masse :',es10.3,'Msol/an')
	 END SELECT
	ENDIF

	WRITE(*,17)thw(Krot) ; WRITE(2,17)thw(Krot)
17	FORMAT('Type de rotation utilisé : ',a)

c formation de la chaine precision
	WRITE(text1,213)m_qs ; WRITE(text2,213)m_ch
213	FORMAT(i3)
	methode='CESAM2k version '//version
	IF(en_masse)THEN
	 methode=TRIM(methode)//' lagr colloc'//text1//text2
	ELSE
	 methode=TRIM(methode)//' eulr colloc'//text1//text2
	ENDIF
	methode=TRIM(methode)//' '//precision
	IF(diffusion)THEN
	 methode=TRIM(methode)//' diffus'
	ELSE
	 methode=TRIM(methode)//' no diffus'
	ENDIF

c initialisation des modèles
	chain_atm=TRIM(nom_fich2)//'.atm'
	chain_don=TRIM(nom_fich2)//'.don'
	chain_rep=TRIM(nom_fich1)

	SELECT CASE(un23)

c poursuite d'une évolution modèle repris en binaire
	CASE(1)
	 OPEN(unit=4,form='unformatted',status='old',file=nom_fich1)
	 READ(4)nep,m_qsp,n_qs,knot,nchimp,n_ch,m_chp,knotc,Krotp,nrotp,
	1 n_rot,m_rotp,knotr,n_tds,knot_tds,mtotp,alphap,w_rotp,lim_rop,
	3 diffusionp,rot_solidp,precisionp,en_massep,f_eosp,f_opap,
	4 nom_ctesp,nom_pertmp,nom_pertwp,nom_tdetaup,nom_atmp,nom_convp,
	5 nom_nucp,nom_nuc_cplp,nom_diffmp,nom_difftp,nom_diffwp,
	6 nom_etatp,nom_opap
c	 PRINT*,nep,m_qsp,n_qs,knot,nchimp,n_ch
c	 PRINT*,m_chp,knotc,n_tds,knot_tds,mtotp ; PAUSE'case 1'

c vérification des concordances
	 IF(mtotp /= mtot)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1700)mtotp,TRIM(chain_rep),mtot,TRIM(chain_don)
	   WRITE(2,1700)mtotp,TRIM(chain_rep),mtot,TRIM(chain_don)
1700	   FORMAT('STOP, the initial mass :',es10.3,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial mass : ',es10.3,
	3  ' read the input file : ',a)
	  CASE DEFAULT
	   WRITE(*,700)mtotp,TRIM(chain_rep),mtot,TRIM(chain_don)
	   WRITE(2,700)mtotp,TRIM(chain_rep),mtot,TRIM(chain_don)
700	   FORMAT('ARRET, masse init. modèle repris : ',es10.3,
	1 ', fichier : ',a,/,'/= masse init. des données : ',es10.3,
	2' fichier : ',a)
	  END SELECT
	  STOP
	 ENDIF

	 IF(alphap /= alpha)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1701)alphap,TRIM(chain_rep),alpha,TRIM(chain_don)
	   WRITE(2,1701)alphap,TRIM(chain_rep),alpha,TRIM(chain_don)
1701	   FORMAT('STOP, the initial mixing length :',es10.3,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial mixing length : ',es10.3,
	3  ' read the input file : ',a,/,
	4  'does one keep on computing  with the continued model?',/,
	5  'enter o for yes; n for no')
	  CASE DEFAULT
	   WRITE(*,701)alphap,TRIM(chain_rep),alpha,TRIM(chain_don)
	   WRITE(2,701)alphap,TRIM(chain_rep),alpha,TRIM(chain_don)
701	   FORMAT('ARRET, long. mel. modèle repris : ',es10.3,
	1 ', fichier : ',a,/,'/= long. mel. des données : ',es10.3,
	2 ', fichier : ',a,/,
	3 'Poursuit-on avec la nouvelle valeur? entrer o/n')
	  END SELECT
	  READ*,oui
	  IF(oui == 'n')THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1702) ; WRITE(*,1702) ; STOP
1702	    FORMAT('STOP')
	   CASE DEFAULT
	    WRITE(*,702) ; WRITE(*,702) ; STOP
702	    FORMAT('ARRET')
	   END SELECT
	  ELSE
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1703)chain_don ; WRITE(2,1703)chain_don
1703	    FORMAT('The computation is carried on with the data of',/,
	1   'the file : ',a)
	   CASE DEFAULT
	    WRITE(*,703)chain_don ; WRITE(2,703)chain_don
703	    FORMAT('on poursuit avec les données du fichier : ',a)
	   END SELECT
	  ENDIF
	 ENDIF

c il faut les mêmes paramètres de rotation
	 IF(nrotp /= nrot)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1129)nrotp,TRIM(chain_rep),nrot,TRIM(chain_don)
	   WRITE(2,1129)nrotp,TRIM(chain_rep),nrot,TRIM(chain_don)
1129	   FORMAT('STOP, nb. of var. for rotation :',i2,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the nb. of var. for rotation : ',i2,
	3  ' read the input file : ',a)
	  CASE DEFAULT
	   WRITE(*,129)nrotp,TRIM(chain_rep),nrot,TRIM(chain_don)
	   WRITE(2,129)nrotp,TRIM(chain_rep),nrot,TRIM(chain_don)
129	   FORMAT('ARRET, nb. coeff. de rotation modèle repris : ',i2,
	1  ', fichier : ',a,/,'/= nb. coeff. de rotation des données : ',
	2  i2,' fichier : ',a)
	  END SELECT
	  STOP
	 ENDIF

	 IF(m_rotp /= m_rot)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1131)m_rotp,TRIM(chain_rep),m_rot,TRIM(chain_don)
	   WRITE(2,1131)m_rotp,TRIM(chain_rep),m_rot,TRIM(chain_don)
1131	   FORMAT('STOP, the spline''s order for the rotation :',i2,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from those : ',i2,' read the input file : ',a)
	  CASE DEFAULT
	   WRITE(*,131)m_rotp,TRIM(chain_rep),m_rot,TRIM(chain_don)
	   WRITE(2,131)m_rotp,TRIM(chain_rep),m_rot,TRIM(chain_don)
131	   FORMAT('ARRET, l''ordre des splines :',i2,
	1   ' pour la rotation du modèle repris dans le fichier : ',a,/,
	2   'diffère de celui des données : ',i2,' fichier : ',a)
	  END SELECT
	  STOP
	 ENDIF

	 IF(Krotp /= Krot)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1132)Krotp,TRIM(chain_rep),Krot,TRIM(chain_don)
	   WRITE(2,1132)Krotp,TRIM(chain_rep),Krot,TRIM(chain_don)
1132	   FORMAT('STOP, the flag for rotation Krot :',i2,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial flag for rotation  Krot : ',i2,
	3  ' read the input file : ',a)
	  CASE DEFAULT
	   WRITE(*,132)Krotp,TRIM(chain_rep),Krot,TRIM(chain_don)
	   WRITE(2,132)Krotp,TRIM(chain_rep),Krot,TRIM(chain_don)
132	   FORMAT('ARRET, le flag de rotation du modèle repris Krot : ',i2,
	1  ', fichier : ',a,/,'/= flag de rotation des données Krot : ',i2,
	2  ' fichier : ',a)
	  END SELECT
	  STOP
	 ENDIF

c si la rotation a été entrée en kms/s w = v / R
	 IF(w_rotp /= w_rot .AND. (unit /= 'kms/s'))THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1704)w_rotp,TRIM(chain_rep),w_rot,TRIM(chain_don)
	   WRITE(2,1704)w_rotp,TRIM(chain_rep),w_rot,TRIM(chain_don)
1704	   FORMAT('STOP, the initial angular velocity :',es10.3,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial angular velocity : ',es10.3,
	3  ' read the input file : ',a)
	  CASE DEFAULT
	   WRITE(*,704)w_rotp,TRIM(chain_rep),w_rot,TRIM(chain_don)
	   WRITE(2,704)w_rotp,TRIM(chain_rep),w_rot,TRIM(chain_don)
704	   FORMAT('ARRET, vit. ang. ini. modèle repris : ',es10.3,
	1  ', fichier : ',a,/,'/= vit. ang. init. des données : ',es10.3,
	2  ' fichier : ',a)
	  END SELECT
	  STOP
	 ENDIF

c il convient d'utiliser la même condition limite pour l'atmosphère
	 IF(lim_rop .NEQV. lim_ro)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1705)lim_rop,TRIM(chain_rep),lim_ro,TRIM(chain_don)
	   WRITE(2,1705)lim_rop,TRIM(chain_rep),lim_ro,TRIM(chain_don)
1705	   FORMAT('WARNING, the boundary conditions lim_ro : ',l1,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the boundary conditions lim_ro : ',l1,/,
	3  'read the input file : ',a,/,
	4  'The computation is carried on with the data of the input file')
	  CASE DEFAULT
	   WRITE(*,705)lim_rop,TRIM(chain_rep),lim_ro,TRIM(chain_don)
	   WRITE(2,705)lim_rop,TRIM(chain_rep),lim_ro,TRIM(chain_don)
705	   FORMAT('ATTENTION, la condition limite lim_ro : ',l1,/,
	1  'du modèle repris du fichier : ',a,/,
	2  'diffèrent des conditions limites lim_ro : ',l1,/,
	3  'du fichier des données: ',a,/,
	4  'on poursuit avec la valeur du fichier des données')
	  END SELECT
	 ENDIF

c on peut poursuivre sans diffusion un modèle avec diffusion et réciproquement
	 IF(diffusionp .NEQV. diffusion)THEN
c	 IF(.FALSE.)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1706)diffusionp,TRIM(chain_rep),diffusion,TRIM(chain_don)
	   WRITE(2,1706)diffusionp,TRIM(chain_rep),diffusion,TRIM(chain_don)
1706	   FORMAT('WARNIG, the diffusion parameter diffusion : ',l1,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial diffusion parameter diffusion : ',l1,/,
	3  'read the input file : ',a)
	  CASE DEFAULT
	   WRITE(*,706)diffusionp,TRIM(chain_rep),diffusion,TRIM(chain_don)
	   WRITE(2,706)diffusionp,TRIM(chain_rep),diffusion,TRIM(chain_don)
706	   FORMAT('ATTENTION, le paramètre de diffusion : ',l1,/,
	1  'du modèle repris du fichier : ',a,/,
	2  'diffère du paramètre de diffusion : ',l1,/,
	3  'du fichier des données : ',a)
	  END SELECT
c	  STOP
	 ENDIF

c on ne peut poursuivre en rotation solide par un modèle en rotation
c non solide et réciproquement
	 IF(rot_solidp .NEQV. rot_solid)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1707)rot_solidp,TRIM(chain_rep),rot_solid,TRIM(chain_don)
	   WRITE(2,1707)rot_solidp,TRIM(chain_rep),rot_solid,TRIM(chain_don)
1707	   FORMAT('STOP, the rotation parameter rot_solid : ',l1,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial rotation parameter rot_solid : ',l1,/,
	3  'read the input file : ',a)
	  CASE DEFAULT
	   WRITE(*,707)rot_solidp,TRIM(chain_rep),rot_solid,TRIM(chain_don)
	   WRITE(2,707)rot_solidp,TRIM(chain_rep),rot_solid,TRIM(chain_don)
707	   FORMAT('ARRET, le paramètre de rotation rot_solid : ',l1,/,
	1  'du modèle repris du fichier : ',a,/,
	2  'diffère du paramètre de rotation rot_solid : ',l1,/,
	3  'du fichier des données : ',a)
	  END SELECT
	  STOP
	 ENDIF

c on ne peut changer de masse initiale en cours d'évolution
	 IF(en_massep .NEQV. en_masse)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1708)en_massep,TRIM(chain_rep),en_masse,TRIM(chain_don)
	   WRITE(2,1708)en_massep,TRIM(chain_rep),en_masse,TRIM(chain_don)
1708	   FORMAT('STOP, the structure parameter en_masse : ',l1,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial structure parameter en_masse : ',l1,/,
	3  'read the input file : ',a)
	  CASE DEFAULT
	   WRITE(*,708)en_massep,TRIM(chain_rep),en_masse,TRIM(chain_don)
	   WRITE(2,708)en_massep,TRIM(chain_rep),en_masse,TRIM(chain_don)
708	   FORMAT('ARRET, le paramètre de structure en_masse : ',l1,/,
	1  'du modèle repris du fichier : ',a,/,
	2  'diffère du paramètre de structure en_masse : ',l1,/,
	3  'du fichier des données : ',a)
	  END SELECT
	  STOP
	 ENDIF

c il faut garder le même type de précision
	 IF(precisionp /= precision)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1709)precisionp,TRIM(chain_rep),precision,TRIM(chain_don)
	   WRITE(2,1709)precisionp,TRIM(chain_rep),precision,TRIM(chain_don)
1709	   FORMAT('STOP, the precision parameter precision : ',a,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial precision parameter precision : ',a,/,
	3  'read the input file : ',a)
	  CASE DEFAULT
	   WRITE(*,709)precisionp,TRIM(chain_rep),precision,TRIM(chain_don)
	   WRITE(2,709)precisionp,TRIM(chain_rep),precision,TRIM(chain_don)
709	   FORMAT('ARRET, le paramètre de précision : ',a,/,
	1  'du modèle repris du fichier : ',a,/,
	2  'diffère du paramètre de precision : ',a,/,
	3  'du fichier des données : ',a)
	  END SELECT

c arret sauf s'il y a réglage
	  IF(.NOT.(precision == 'rg' .OR. precisionp == 'rg'))STOP
	 ENDIF

c il convient de conserver les mêmes valeurs des constantes
	 IF(nom_ctesp /= nom_ctes)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1710)nom_ctesp,TRIM(chain_rep),nom_ctes,TRIM(chain_don)
	   WRITE(2,1710)nom_ctesp,TRIM(chain_rep),nom_ctes,TRIM(chain_don)
1710	   FORMAT('WARNING, the physical constants parameter : ',a,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial physical constants parameter : ',
	3  a,/, 'read the input file : ',a,/,
	4  'The computation is carried on with the data of the input file')
	  CASE DEFAULT
	   WRITE(*,710)nom_ctesp,TRIM(chain_rep),nom_ctes,TRIM(chain_don)
	   WRITE(2,710)nom_ctesp,TRIM(chain_rep),nom_ctes,TRIM(chain_don)
710	   FORMAT('ATTENTION, le paramètre de constantes physiques : ',a,/,
	1  'du modèle repris du fichier : ',a,/,
	2  'diffère du paramètre de constantes physiques : ',a,/,
	3  'du fichier des données : ',a,/,
	4  'on poursuit avec la valeur du fichier des données')
	  END SELECT
	 ENDIF

c il convient de conserver le même type de perte de masse
	 IF(nom_pertmp /= nom_pertm)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1711)nom_pertmp,TRIM(chain_rep),nom_pertm,TRIM(chain_don)
	   WRITE(2,1711)nom_pertmp,TRIM(chain_rep),nom_pertm,TRIM(chain_don)
1711	   FORMAT('WARNING, The mass loss parameter : ',a,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial mass loss parameter : ',a,/,
	3  'read the input file : ',a,/,
	4  'The computation is carried on with the data of the input file')
	  CASE DEFAULT
	   WRITE(*,711)nom_pertmp,TRIM(chain_rep),nom_pertm,TRIM(chain_don)
	   WRITE(2,711)nom_pertmp,TRIM(chain_rep),nom_pertm,TRIM(chain_don)
711	   FORMAT('ATTENTION, le paramètre de perte de masse : ',a,/,
	1  'du modèle repris du fichier : ',a,/,
	2  'diffère du paramètre de perte de masse : ',a,/,
	3  'du fichier des données : ',a,/,
	4  'on poursuit avec la valeur du fichier des données')
	  END SELECT
	 ENDIF

c il convient de conserver le même type de perte de moment cinétique
	 IF((Krot == 3 .OR. Krot == 4) .AND. nom_pertwp /= nom_pertw)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1699)nom_pertwp,TRIM(chain_rep),nom_pertw,TRIM(chain_don)
	   WRITE(2,1699)nom_pertwp,TRIM(chain_rep),nom_pertw,TRIM(chain_don)
1699	   FORMAT('WARNING, The ang. momentum loss parameter : ',a,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial ang. momentum loss parameter : ',a,/,
	3  'read the input file : ',a,/,
	4  'The computation is carried on with the data of the input file')
	  CASE DEFAULT
	   WRITE(*,699)nom_pertwp,TRIM(chain_rep),nom_pertw,TRIM(chain_don)
	   WRITE(2,699)nom_pertwp,TRIM(chain_rep),nom_pertw,TRIM(chain_don)
699	   FORMAT('ATTENTION, le paramètre de perte de moment cin. : ',a,/,
	1  'du modèle repris du fichier : ',a,/,
	2  'diffère du paramètre de perte de moment cin. : ',a,/,
	3  'du fichier des données : ',a,/,
	4  'on poursuit avec la valeur du fichier des données')
	  END SELECT
	 ENDIF

c il conviendrait de conserver la même loi T(tau)
	 IF(nom_tdetaup /= nom_tdetau)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1712)nom_tdetaup,TRIM(chain_rep),nom_tdetau,TRIM(chain_don)
	   WRITE(2,1712)nom_tdetaup,TRIM(chain_rep),nom_tdetau,TRIM(chain_don)
1712	   FORMAT('WRANING, the T(tau)-law parameter : ',a,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial T(tau)-law parameter : ',a,/,
	3  'read the input file : ',a,/,
	4  'The computation is carried on with the data of the input file')
	  CASE DEFAULT
	   WRITE(*,712)nom_tdetaup,TRIM(chain_rep),nom_tdetau,TRIM(chain_don)
	   WRITE(2,712)nom_tdetaup,TRIM(chain_rep),nom_tdetau,TRIM(chain_don)
712	   FORMAT('ATTENTION, le paramètre de loi T(tau) : ',a,/,
	1  'du modèle repris du fichier : ',a,/,
	2  'diffère du paramètre de loi T(tau) : ',a,/,
	3  'du fichier des données : ',a,/,
	4  'on poursuit avec la valeur du fichier des données')
	  END SELECT
	 ENDIF

c il convienfrait de conserver la même atmosphère
	 IF(nom_atmp /= nom_atm)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1713)nom_atmp,TRIM(chain_rep),nom_atm,TRIM(chain_don)
	   WRITE(2,1713)nom_atmp,TRIM(chain_rep),nom_atm,TRIM(chain_don)
1713	   FORMAT('WARNING, the atmosphere parameter : ',a,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial atmosphere parameter : ',a,/,
	3  'read the input file : ',a,/,
	4  'The computation is carried on with the data of the input file')
	  CASE DEFAULT
	   WRITE(*,713)nom_atmp,TRIM(chain_rep),nom_atm,TRIM(chain_don)
	   WRITE(2,713)nom_atmp,TRIM(chain_rep),nom_atm,TRIM(chain_don)
713	   FORMAT('ATTENTION, le paramètre d''atmosphère : ',a,/,
	1  'du modèle repris du fichier : ',a,/,
	2  'diffère du paramètre d''atmosphère : ',a,/,
	3  'du fichier des données : ',a,/,
	4  'on poursuit avec la valeur du fichier des données')
	  END SELECT
	 ENDIF

c si on poursuit avec une atmosphère simplifiée il faut n_atm=0
	 IF(nom_atm /= 'lim_atm')n_atm=0
	 IF(nom_convp /= nom_conv)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1714)nom_convp,TRIM(chain_rep),nom_conv,TRIM(chain_don)
	   WRITE(2,1714)nom_convp,TRIM(chain_rep),nom_conv,TRIM(chain_don)
1714	   FORMAT('WARNING, the convection parameter : ',a,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial convection parameter : ',a,/,
	3  'read the input file : ',a,/,
	4  'The computation is carried on with the data of the input file')
	  CASE DEFAULT
	   WRITE(*,714)nom_convp,TRIM(chain_rep),nom_conv,TRIM(chain_don)
	   WRITE(2,714)nom_convp,TRIM(chain_rep),nom_conv,TRIM(chain_don)
714	   FORMAT('ATTENTION, le paramètre de convection : ',a,/,
	1  'du modèle repris du fichier : ',a,/,
	2  'diffère du paramètre de convection : ',a,/,
	3  'du fichier des données : ',a,/,
	4  'on poursuit avec la valeur du fichier des données')
	  END SELECT
	 ENDIF

c on ne peut poursuivre avec un réseau nucléaire différent
	 IF(nom_nucp /= nom_nuc)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1715)nom_nucp,TRIM(chain_rep),nom_nuc,TRIM(chain_don)
	   WRITE(2,1715)nom_nucp,TRIM(chain_rep),nom_nuc,TRIM(chain_don)
1715	   FORMAT('The nuclear reactions parameter : ',a,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial nuclear reactions parameter : ',a,/,
	3  'read in the input file : ',a,/,
	4  'The computation is carried on with the data of the input file',
	5  /,'only if the chemical elements are identical')
	  CASE DEFAULT
	   WRITE(*,715)nom_nucp,TRIM(chain_rep),nom_nuc,TRIM(chain_don)
	   WRITE(2,715)nom_nucp,TRIM(chain_rep),nom_nuc,TRIM(chain_don)
715	   FORMAT('Le paramètre de réactions nucléaires : ',a,/,
	1  'du modèle repris du fichier : ',a,/,
	2  'diffère du paramètre de réactions nucléaires : ',a,/,
	3  'du fichier des données : ',a,/,
	4  'on poursuit avec la valeur du fichier des données que si les',/,
	5  'éléments chimiques utilisés sont identiques')
	  END SELECT
	  IF((TRIM(nom_nucp) == 'ppcno9' .AND. TRIM(nom_nuc) == 'ppcno3a9').OR. 
	1 (TRIM(nom_nucp) == 'ppcno3a9' .AND. TRIM(nom_nuc) == 'ppcno9'))THEN	
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1787) ; WRITE(2,1787)	
1787	    FORMAT('The computation is carried on.')	
	   CASE DEFAULT
	    WRITE(*,787) ; WRITE(2,787)	
787	    FORMAT('Le calcul est poursuivi.')	  	  
	   END SELECT
	  ELSE
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1788) ; WRITE(2,1788)	
1788	    FORMAT('STOP, the chemical elements are not identical.')	
	   CASE DEFAULT
	    WRITE(*,788) ; WRITE(2,788)	
788	    FORMAT('ARRET, les éléments chimiques utilisés diffèrent ')	  	  
	   END SELECT
	   STOP
	  ENDIF
	 ENDIF

c il convient de poursuivre avec la même compilation de réactions nucléaires
	 IF(nom_nuc_cplp /= nom_nuc_cpl)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1716)nom_nuc_cplp,TRIM(chain_rep),nom_nuc_cpl,TRIM(chain_don)
	   WRITE(2,1716)nom_nuc_cplp,TRIM(chain_rep),nom_nuc_cpl,TRIM(chain_don)
1716	   FORMAT('WARNING, the compilation of nuclear reactions : ',a,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial compilation of nuclear reactions : ',
	3  a,/, 'read the input file : ',a,/,
	4  'The computation is carried on with the data of the input file',
	5  /,'if the chemical elements are identical')
	  CASE DEFAULT
	   WRITE(*,716)nom_nuc_cplp,TRIM(chain_rep),nom_nuc_cpl,TRIM(chain_don)
	   WRITE(2,716)nom_nuc_cplp,TRIM(chain_rep),nom_nuc_cpl,TRIM(chain_don)
716	   FORMAT('ATTENTION, la compilation des réactions nucléaires : ',a,/,
	1  'du modèle repris du fichier : ',a,/,
	2  'diffère de la compilation des réactions nucléaires : ',a,/,
	3  'du fichier des données : ',a,/,
	4  'on poursuit avec la valeur du fichier des données si les',/,
	5  'éléments chimiques utilisés sont identiques')
	  END SELECT
	 ENDIF

c il convient de poursuivre avec la même diffusion microscopique
	 IF(nom_diffmp /= nom_diffm)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1717)nom_diffmp,TRIM(chain_rep),nom_diffm,TRIM(chain_don)
	   WRITE(2,1717)nom_diffmp,TRIM(chain_rep),nom_diffm,TRIM(chain_don)
1717	   FORMAT('WARNING, the microscopic diffusion parameter : ',a,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial microscopic diffusion parameter : ',
	3  a,/, 'read the input file : ',a,/,
	4  'The computation is carried on with the data of the input file')
	  CASE DEFAULT
	   WRITE(*,717)nom_diffmp,TRIM(chain_rep),nom_diffm,TRIM(chain_don)
	   WRITE(2,717)nom_diffmp,TRIM(chain_rep),nom_diffm,TRIM(chain_don)
717	   FORMAT('ATTENTION, le paramètre de diffusion microscopique : ',a,/,
	1  'du modèle repris du fichier : ',a,/,
	2  'diffère du paramètre de diffusion microscopique : ',a,/,
	3  'du fichier des données : ',a,/,
	4  'on poursuit avec la valeur du fichier des données')
	  END SELECT
	 ENDIF

c il convient de poursuivre avec la même diffusion du moment cinétique
	 IF(nom_diffwp /= nom_diffw)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1724)nom_diffwp,TRIM(chain_rep),nom_diffw,TRIM(chain_don)
	   WRITE(2,1724)nom_diffwp,TRIM(chain_rep),nom_diffw,TRIM(chain_don)
1724	   FORMAT('WARNING, the turbulent diffusion parameter : ',a,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial turbulent diffusion parameter : ',
	3  a,/, 'read the input file : ',a,/,
	4  'The computation is carried on with the data of the input file')
	  CASE DEFAULT
	   WRITE(*,724)nom_diffwp,TRIM(chain_rep),nom_diffw,TRIM(chain_don)
	   WRITE(2,724)nom_diffwp,TRIM(chain_rep),nom_diffw,TRIM(chain_don)
724	   FORMAT('ATTENTION, le paramètre de diffusion turbulente : ',a,/,
	1  'du modèle repris du fichier : ',a,/,
	2  'diffère du paramètre de diffusion turbulente : ',a,/,
	3  'du fichier des données : ',a,/,
	4  'on poursuit avec la valeur du fichier des données')
	  END SELECT
	 ENDIF

c il convient de poursuivre avec la même diffusion turbulente
	 IF(nom_difftp /= nom_difft)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1718)nom_difftp,TRIM(chain_rep),nom_difft,TRIM(chain_don)
	   WRITE(2,1718)nom_difftp,TRIM(chain_rep),nom_difft,TRIM(chain_don)
1718	   FORMAT('WARNING, the turbulent diffusion routine : ',a,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial turbulent diffusion routine : ',
	3  a,/, 'read the input file : ',a,/,
	4  'The computation is carried on with the data of the input file')
	  CASE DEFAULT
	   WRITE(*,718)nom_difftp,TRIM(chain_rep),nom_difft,TRIM(chain_don)
	   WRITE(2,718)nom_difftp,TRIM(chain_rep),nom_difft,TRIM(chain_don)
718	   FORMAT('ATTENTION, la routine de diffusion turbulente : ',a,/,
	1  'du modèle repris du fichier : ',a,/,
	2  'diffère de la routine de diffusion turbulente : ',a,/,
	3  'du fichier des données : ',a,/,
	4  'on poursuit avec la valeur du fichier des données')
	  END SELECT
	 ENDIF

c il convient de poursuivre avec la même EOS
	 IF(nom_etatp /= nom_etat)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1719)nom_etatp,TRIM(chain_rep),nom_etat,TRIM(chain_don)
	   WRITE(2,1719)nom_etatp,TRIM(chain_rep),nom_etat,TRIM(chain_don)
1719	   FORMAT('WARNING, the equation of state parameter : ',a,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial equation of state parameter : ',
	3  a,/, 'read the input file : ',a,/,
	4  'The computation is carried on with the data of the input file',
	5  /,'The program stops if the EOS file is not correct')
	  CASE DEFAULT
	   WRITE(*,719)nom_etatp,TRIM(chain_rep),nom_etat,TRIM(chain_don)
	   WRITE(2,719)nom_etatp,TRIM(chain_rep),nom_etat,TRIM(chain_don)
719	   FORMAT('ATTENTION, le paramètre d''équation d''état : ',a,/,
	1  'du modèle repris du fichier : ',a,/,
	2  'diffère du paramètre d''équation d''état : ',a,/,
	3  'du fichier des données : ',a,/,
	4  'on poursuit avec la valeur du fichier des données PLANTAGE ',/,
	5  'si le nom du fichier d''équation d''état est incorrect')
	  END SELECT
	 ENDIF

c il convient de poursuivre avec la même opacité
	 IF(nom_opap /= nom_opa)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1720)nom_opap,TRIM(chain_rep),nom_opa,TRIM(chain_don)
	   WRITE(2,1720)nom_opap,TRIM(chain_rep),nom_opa,TRIM(chain_don)
1720	   FORMAT('WARNING, the opacity parameter : ',a,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the initial opacity parameter : ',
	3  a,/, 'read the input file : ',a,/,
	4  'The computation is carried on with the data of the input file',
	5  /,'The program stops if the opacity file is not correct')
	  CASE DEFAULT
	   WRITE(*,720)nom_opap,TRIM(chain_rep),nom_opa,TRIM(chain_don)
	   WRITE(2,720)nom_opap,TRIM(chain_rep),nom_opa,TRIM(chain_don)
720	   FORMAT('ATTENTION, le paramètre d''opacité : ',a,/,
	1  'du modèle repris du fichier : ',a,/,
	2  'diffère du paramètre d''opacité : ',a,/,
	3  'du fichier des données : ',a,/,
	4  'on poursuit avec la valeur du fichier des données PLANTAGE ',/,
	5  'si le nom du fichier d''opacité est incorrect')
	  END SELECT
	 ENDIF

c il convient de poursuivre avec le même ordre de splines pour les
c variables de structure
	 IF(m_qsp /= m_qs)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1004)m_qsp,TRIM(chain_rep),m_qs,TRIM(chain_don)
	   WRITE(2,1004)m_qsp,TRIM(chain_rep),m_qs,TRIM(chain_don)
1004	   FORMAT('WARNING, the order of the spline for the quasi-static',/,
	1  'equilibrium m_qs=',i3,/,
	2  'read in the file of the model to be continued : ',a,/,
	3  'differs from the order of the spline m_qs=',i3,
	4  ', read the input file : ',a,/,'one keeps on computing with',/,
	5  'the value read in the file of the model to be continued')
	  CASE DEFAULT
	   WRITE(*,4)m_qsp,TRIM(chain_rep),m_qs,TRIM(chain_don)
	   WRITE(2,4)m_qsp,TRIM(chain_rep),m_qs,TRIM(chain_don)
4	   FORMAT('ATTENTION, l''ordre de la spline m_qs=',i3,/,
	1  'pour l''équilibre quasi-statique du modèle : ',a,/,
	2  'diffère de l''ordre de la spline m_qs=',i3,/,
	3  'du modèle à poursuivre : ',a,/,
	4  'on garde l''ordre de la spline du modèle repris')
	  END SELECT
	 ENDIF
         m_qs=m_qsp

c il convient de poursuivre avec le même ordre de splines pour la
c composition chimique
	 IF(m_chp /= m_ch)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1006)m_chp,TRIM(chain_rep),m_ch,TRIM(chain_don)
	   WRITE(2,1006)m_chp,TRIM(chain_rep),m_ch,TRIM(chain_don)
1006	   FORMAT('WARNING, the order of the spline for the interpolation of',/,
	1  'the chemical composition, m_ch=',i3,/,
	2  'read in the file of the model to be continued : ',a,/,
	3  'differs from the order of the spline m_ch=',i3,
	4  'read the input file : ',a,/,'one keeps on computing with',/,
	5  'the value read in the file of the model to be continued')
	  CASE DEFAULT
	   WRITE(*,6)m_chp,TRIM(chain_rep),m_ch,TRIM(chain_don)
	   WRITE(2,6)m_chp,TRIM(chain_rep),m_ch,TRIM(chain_don)
6	   FORMAT('ATTENTION, l''ordre de la spline m_ch=',i3,/,
	1  'pour la composition chimique du modèle repris : ',a,/,
	2  'diffère de l''ordre de la spline m_ch=',i3,/,
	3  'du modèle à poursuivre : ',a,/,
	4  'on garde l''ordre de la spline du modèle repris')
	  END SELECT
	 ENDIF
         m_ch=m_chp

c reprise du modèle, on ajuste les dim. à celles du modèle repris
	 ord_qs=m_qs+r_qs ; dim_qs=knot-ord_qs ; dim_ch=knotc-m_ch
	 dim_rot=knotr-ord_rot ; m_tds=knot_tds-n_tds
	 ALLOCATE(bp(nep,dim_qs),q(n_qs),qt(knot),nom_elemp(nchimp),
	1 chim(nchimp,dim_ch),mc(n_ch),mct(knotc),mrot(n_rot),
	2 mrott(knotr),tds(1,knot_tds-m_tds),x_tds(n_tds),
	3 xt_tds(knot_tds),m23(n_qs),r2(n_qs),rota(nrot,dim_rot))
	 REWIND(unit=4)
	 READ(4)ne,m_qs,n_qs,knot,nchim,n_ch,m_ch,knotc,Krot,nrot,n_rot,
	1 m_rot,knotr,n_tds,knot_tds,mtotp,alphap,w_rotp,
	2 lim_rop,diffusionp,rot_solidp,precisionp,en_massep,f_eosp,
	3 f_opap,nom_ctesp,nom_pertmp,nom_pertwp,nom_tdetaup,nom_atmp,
	4 nom_convp,nom_nucp,nom_nuc_cplp,nom_diffmp,nom_difftp,nom_diffwp,
	5 nom_etatp,nom_opap,nom_elemp,
	6 bp,q,qt,chim,mc,mct,rota,mrot,mrott,tds,x_tds,xt_tds,m23,
	7 r2,m_zc,r_zc,r_ov,age,dtp,dts,mstar,rstar,mw_tot,wrot,jlim,
	8 lconv,lim,model_num
	 CLOSE(unit=4)
	 model_num=model_num-1	!est augmenté d'une unité ligne ~4155
c	 WRITE(*,2000)dtp,dts ; PAUSE'lect'

	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1021)age ; WRITE(2,1021)age
1021	  FORMAT('Start from binary model aged',es10.3,/)
	 CASE DEFAULT
	  WRITE(*,21)age ; WRITE(2,21)age
21	  FORMAT('Reprise du modèle binaire d''âge',es10.3,/)
	 END SELECT

c si on utilise une grille fixe pour la composition chimique
	 IF(grille_fixe)THEN
	  ALLOCATE(mc_fixe(n_qs))
	  nc_fixe=n_qs ; mc_fixe=m23
	  IF(.NOT.en_masse)mc_fixe=mc_fixe**(2.d0/3.d0)
	 ENDIF

c appel d'initialisation pour tabulation des réactions nucléaires
	 ALLOCATE(comp(0),dcomp(0),jac(0,0),ex(0))
	 CALL nuc(1.5d+07,1.5d+02,comp,dcomp,jac,.FALSE.,0,
	1 epsilo,et,ero,ex,hh,be7,b8,n13,o15,f17)

c la rotation
	 WRITE(*,17)thw(Krot) ; WRITE(2,17)thw(Krot)

c vitesse angulaire minimale
	 IF(Krot >=3)rot_min(1)=wrot*1.d-2

	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1023)wrot ; WRITE(2,1023)wrot
1023	  FORMAT('initial angular velocity =',es10.3,' rd/sec',/)
	 CASE DEFAULT
	  WRITE(*,23)wrot ; WRITE(2,23)wrot
23	  FORMAT('vitesse angulaire initiale =',es10.3,' rd/sec',/)
	 END SELECT

c détermination des abondances initiales
	 DEALLOCATE(comp) ; ALLOCATE(comp(nchim))
	 CALL nuc(1.5d+07,1.5d+02,comp,dcomp,jac,.FALSE.,1,
	1 epsilo,et,ero,ex,hh,be7,b8,n13,o15,f17)
c	 PRINT*,nom_elem ; PRINT*,nom_elemp ; PAUSE'elem'

	 IF(nchimp /= nchim)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1010)nchimp,TRIM(chain_rep),nchim,TRIM(chain_don)
1010	   FORMAT('STOP, the initial number of chemical elements : ',i3,/,
	1  'read in the file of the model to be continued : ',a,/,
	2  'differs from the number of chemical elements : ',i3,/,
	3  'read the input file : ',a)
	  CASE DEFAULT
	   WRITE(*,10)nchimp,TRIM(chain_rep),nchim,TRIM(chain_don)
10	   FORMAT('ARRET, le nombre d''éléments chimimiques = ',i3,/,
	1  'du modèle repris : ',a,/,
	2  'diffère du nombre d''éléments chimiques = ',i3,/,
	3  'du modèle à poursuivre : ',a)
	  END SELECT
	  STOP
	 ENDIF
	 ok=.TRUE.
	 DO i=1,nchim
	  ok=ok .AND. nom_elem(i) == nom_elemp(i)
	 ENDDO
	 IF(.NOT.ok)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1011)nom_elemp
	   WRITE(*,1012)TRIM(chain_rep),nom_elem
	   WRITE(*,1032)TRIM(chain_don)
	   WRITE(2,1011)nom_elemp
	   WRITE(2,1012)TRIM(chain_rep),nom_elem
	   WRITE(2,1032)TRIM(chain_don)
1011	   FORMAT('STOP,the chemical elements : ',30a)
1012	   FORMAT('read the input file : ',a,/,
	1  'differ from the chemical elements : ',30a)
1032	   FORMAT('read the input file : ',a)
	  CASE DEFAULT
	   WRITE(*,11)nom_elemp
	   WRITE(*,12)TRIM(chain_rep),nom_elem
	   WRITE(*,32)TRIM(chain_don)
	   WRITE(2,11)nom_elemp
	   WRITE(2,12)TRIM(chain_rep),nom_elem
	   WRITE(2,32)TRIM(chain_don)
11	   FORMAT('ARRET, les élements chimiques : ',30a)
12	   FORMAT('du modèle repris : ',a,/,'diffèrent de ceux : ',30a)
32	   FORMAT('du modèle à poursuivre : ',a)
	  END SELECT
	  STOP
	 ENDIF

c tout semble OK, let's go

c suppression/adjonction de la pression turbulente
c si le modèle repris est avec/sans Pturb
c et que le modèle à calculer est sans/avec Pturb
c il faut enlever/ajouter la variable ln Pgaz
c ne/nep : nb. inconnues du modèle à calculer/repris
	 IF(ne /= nep)THEN
	  ALLOCATE(b(nep,dim_qs)) ; b=bp
	  DEALLOCATE(bp) ; ALLOCATE(bp(ne,dim_qs))

c changement du nb. inconnues
	  IF(ne == 6)THEN !suppression de Pgaz on passe de ne= 7 à 6 inconnues
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1721) ; WRITE(2,721)
1721	    FORMAT('Pturb included in the starting model,',/,
	1   'but not to be included in the model to be calculated')
	   CASE DEFAULT
	    WRITE(*,721) ; WRITE(2,721)
721	    FORMAT('modèle repris avec Pturb,',/,
	1   'modèle à calculer sans Pturb')
	   END SELECT
	   bp=b(1:6,:)
	  ELSE      !adjonction de Pgaz on passe de 6 à 7 inconnues
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1722) ; WRITE(2,722)
1722	    FORMAT('Pturb not included in the starting model,'/,
	1   'but to be included in the model to be calculated')
	   CASE DEFAULT
	    WRITE(*,722) ; WRITE(2,722)
722	    FORMAT('modèle repris sans Pturb,',/,
	1   'modèle à calculer avec Pturb')
	   END SELECT
	   bp(1:6,:)=b ; bp(Ipg,:)=b(1,:)
	  ENDIF     !passage de 6 <---> 7 inconnues
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1723)nep,ne ; WRITE(2,1723)nep,ne
1723	   FORMAT('we switch from ne=',i3,' unknowns to',/,
	1  'ne=',i3,' unknowns')
	  CASE DEFAULT
	   WRITE(*,723)nep,ne ; WRITE(2,723)nep,ne
723	   FORMAT('on passe de ne=',i3,' inconnues à',/,
	1  'ne=',i3,' inconnues')
	  END SELECT

c b est désormais inutile
	  DEALLOCATE(b)
	 ENDIF      !transfert et changement du nb. inconnues
c	 WRITE(*,2000)age,agemax,dtp,dts,dt ; PAUSE'dt'


c	PRINT*,dts ; PAUSE'1'

c gestion du premier pas temporel
	 IF(dts > 0.d0 .AND. agemax > 0.d0)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1111)age,agemax,dtp,dts
1111	   FORMAT('age of the continued model (in units of 10**6 y)',
	1  es10.3,', agemax=',es10.3,/,'previous time step=',es10.3,/,
	2  'estimated optimum time step=',es10.3,/,
	3  'does one choose the optimum time step?',/,
	4  'enter o for yes, enter n for no')
	  CASE DEFAULT
	   WRITE(*,111)age,agemax,dtp,dts
111	   FORMAT('âge (10**6 ans) du modèle repris=',es10.3,
	1  ', agemax=',es10.3,/,'pas temporel précédent=',es10.3,/,
	2  'pas temporel optimal estimé =',es10.3,/,
	3  'utilise-t-on le pas temporel optimal ? entrer o/n')
	  END SELECT
	  READ*,oui
	  IF(oui /= 'o')THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1113)
1113	    FORMAT('enter the new time step, in units of 10**6 years')
	   CASE DEFAULT
	    WRITE(*,113)
113	    FORMAT('entrer le nouveau pas temporel, unité: 10**6 ans')
	   END SELECT
	   READ*,dts
	  ENDIF     !oui /= 'o'

c	  PRINT*,dts,dtmax ; PAUSE'2'

	  IF(dts > dtmax)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1112)dtmax
1112	    FORMAT('dt > dtmax, the time step is set back to dtmax=',es10.3)
	   CASE DEFAULT
	    WRITE(*,112)dtmax
112	    FORMAT('dt > dtmax, le pas temporel est ramené à dtmax=',es10.3)
	   END SELECT
	  ENDIF    !dt > dtmax
	  dt=MIN(dts,dtmax) ; dts=dt !on garde le pas temporel précédent

c	  PRINT*,dt,dts,dtmax ; PAUSE'3'

	  IF(age+dt > agemax)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1053)age+dt,agemax ; WRITE(2,1053)age+dt,agemax
1053	    FORMAT('âge+dt=',es10.3,' > agemax =',es10.3,/,
	1   'tuning of dt')
	   CASE DEFAULT
	    WRITE(*,53)age+dt,agemax ; WRITE(2,53)age+dt,agemax
53	    FORMAT('âge+dt=',es10.3,' > agemax =',es10.3,/,'ajustement de dt')
	   END SELECT
	   dt=MIN(dtmax,MAX(0.d0,agemax-age)) ; dts=dt
	   IF(dt <= 0.d0)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1054) ; WRITE(2,1054)
1054	     FORMAT('STOP, agemax is reached, no further evolution',/,
	1    'before exit, write listing of the starting model?',/,
	2    'enter o for yes, n for no')
	    CASE DEFAULT
	     WRITE(*,54) ; WRITE(2,54)
54	     FORMAT('ARRET, agemax est atteint, pas de poursuite',/,
	1    'liste du modèle repris avant de sortir? entrer o/n')
	    END SELECT
	    READ*,oui
	    IF(oui == 'o')THEN
	     list_sort=.TRUE.
	    ELSE
	     STOP
	    ENDIF
	   ENDIF    !dt <= 0.d0
	  ENDIF     !age+dt > agemax

c	  PRINT*,dt,dts,dtmax ; PAUSE'4'

	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(2,1050)dt ; WRITE(*,1050)dt
1050	   FORMAT(/,'The evolution is continued with the time step',
	1  es10.3,' (million years)',/)
	  CASE DEFAULT
	   WRITE(2,50)dt ; WRITE(*,50)dt
50	   FORMAT(/,'Poursuite de l''évolution avec le pas temporel',
	1  es10.3,' (million d''années)',/)
 	  END SELECT
	  list_rep=.TRUE. !listing du modèle repris
	 ELSEIF(agemax > 0.d0)THEN    !dtp=0, on poursuit un modèle initial
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1114)dt0
1114	   FORMAT('does one use dt0=',es10.3,' for the time step?',/,
	1  'enter o for yes, enter n for no')
	  CASE DEFAULT
	   WRITE(*,114)dt0
114	   FORMAT('utilise-t-on dt0=',es10.3,' pour le pas temporel? o/n')
	  END SELECT
	  READ*,oui
	  IF(oui == 'o')THEN
	   dt=MIN(dt0,dtmax) ; dts=dt
	  ELSE
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1115)
1115	    FORMAT('enter the new time step, in units of 10**6 years')
	   CASE DEFAULT
	    WRITE(*,115)
115	    FORMAT('entrer le nouveau pas temporel, unité : 10**6 ans')
	   END SELECT
	   READ*,dt
	   IF(dt > dtmax)THEN
	    dt=dtmax ; dts=dt ; WRITE(*,112)dtmax
	   ENDIF
	  ENDIF           !oui /= 'o'
	  IF(age+dt > agemax)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1053)age+dt,agemax
	   CASE DEFAULT
	    WRITE(*,53)age+dt,agemax
	   END SELECT
	   dt=MIN(dtmax,MAX(0.d0,agemax-age))

c	  PRINT*,dt,dts,dtmax ; PAUSE'5'

	   IF(dt <= 0.d0)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1054)
	    CASE DEFAULT
	     WRITE(*,54)
	    END SELECT
	    READ*,oui
	    IF(oui == 'o')THEN
	     list_sort=.TRUE.
	    ELSE
	     STOP
	    ENDIF
	   ENDIF   !dt <= 0.d0
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(2,1051)dt ; WRITE(*,1051)dt ; list_rep=.TRUE.
1051	    FORMAT('Evolution starting at age=0 with dt=',es10.3)
	   CASE DEFAULT
	    WRITE(2,51)dt ; WRITE(*,51)dt ; list_rep=.TRUE.
51	    FORMAT('Evolution à partir de l''âge=0 avec dt=',es10.3)
	   END SELECT
	  ENDIF           !age+dt > agemax
	  list_rep=.TRUE. !listing du modèle repris
	 ELSE       !agemax=0 dtp quelconque
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1009)age,agemax ; WRITE(2,1009)age,agemax
1009	   FORMAT('age of the continued model',es10.3,'> agemax : ',es10.3,
	1  /,'one writes the listing of the continued model and exits')
	  CASE DEFAULT
	   WRITE(*,9)age,agemax ; WRITE(2,9)age,agemax
9	   FORMAT('âge du modèle repris : ',es10.3,' > agemax:',es10.3,/,
	1  'on liste le modèle repris et on sort')
	  END SELECT
	  list_sort=.TRUE.
	 ENDIF		 !dtp > 0 et agemax > 0

c initialisation des paramètres pour le vent et les chutes de planétoïdes
	 CALL vent ; CALL planetoides

c Poursuite d'une évolution modèle repris en ASCII
	CASE(-1)
	 IF(agemax <= 0.d0)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1610)agemax ; WRITE(2,1610)agemax
1610	   FORMAT('With a model to be continued stored in ASCII,',/,
	1  'agemax=',es10.3,' is not consistant, STOP')
	  CASE DEFAULT
	   WRITE(*,610)agemax ; WRITE(2,610)agemax
610	   FORMAT('Avec un modèle à poursuivre repris en ASCII, agemax=',
	1  es10.3,' est incohérent, ARRET')
	  END SELECT
	  STOP
	 ENDIF
	 OPEN(unit=4,form='formatted',status='old',file=nom_fich1,
	1 delim='apostrophe')
	 READ(4,2004)n,nchim_ascii,ihe4_ascii,Krot_ascii
2004	 FORMAT(10i5)
	 ALLOCATE(nom_elem_ascii(nchim_ascii))
	 READ(4,2006)nom_elem_ascii(1:nchim_ascii)
2006	 FORMAT(20a4)
	 READ(4,2005)age,dtp,mtot_ascii,x0_ascii,y0_ascii
	 IF(x0_ascii /= x0 .OR. y0_ascii /= y0)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1606)x0_ascii,x0,y0_ascii,y0
	   WRITE(2,1606)x0_ascii,x0,y0_ascii,y0
1606	   FORMAT('STOP, x0_ascii=',es10.3,'is different from x0=',es10.3,
	1  'or y0_ascii=',es10.3,'is different of y0',es10.3)
	  CASE DEFAULT
	   WRITE(*,606)x0_ascii,x0,y0_ascii,y0
	   WRITE(2,606)x0_ascii,x0,y0_ascii,y0
606	   FORMAT('ARRET,x0_ascii=',es10.3,' diffère de x0=',es10.3,
	1  ' ou y0_ascii=',es10.3,' diffère de y0=',es10.3)
	   END SELECT
	  STOP
	 ENDIF
	 IF(mtot /= mtot_ascii)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1607)mtot_ascii,mtot,TRIM(chain_don)
	   WRITE(2,1607)mtot_ascii,mtot,TRIM(chain_don)
1607	   FORMAT('the mass of the continued model in ascii mtot=',es10.3,
	1  'is different from mtot=',es10.3,/,'read in the file : ',a)
	  CASE DEFAULT
	   WRITE(*,607)mtot_ascii,mtot,TRIM(chain_don)
	   WRITE(2,607)mtot_ascii,mtot,TRIM(chain_don)
607	   FORMAT('STOP, la masse du modèle repris en ascii mtot=',es10.3,
	1  ' diffère de mtot=',es10.3,/,'lu dans le fichier : ',a)
	  END SELECT
	 ENDIF

c appel d'initialisation pour tabulation des réactions nucléaires
c allocations fictives

	 ALLOCATE(comp(0),dcomp(0),jac(0,0),ex(0))
	 CALL nuc(1.5d+07,1.5d+02,comp,dcomp,jac,.FALSE.,0,
	1 epsilo,et,ero,ex,hh,be7,b8,n13,o15,f17)

c il faut Krot >= 3
	 IF(Krot >= 3 .AND. Krot /= Krot_ascii)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1609)Krot_ascii,TRIM(chain_rep),Krot,TRIM(chain_don)
	   WRITE(2,1609)Krot_ascii,TRIM(chain_rep),Krot,TRIM(chain_don)
1609	   FORMAT('STOP, the option for the rotation Krot=',i3,/,
	1  'is different for the model to be continued : ',a,/,
	1  'from the  option Krot=',i3,/,'of the continued model : ',a)
	  CASE DEFAULT
	   WRITE(*,609)Krot_ascii,TRIM(chain_rep),Krot,TRIM(chain_don)
	   WRITE(2,609)Krot_ascii,TRIM(chain_rep),Krot,TRIM(chain_don)
609	   FORMAT('ARRET, le type de rotation Krot=',i3,/,
	1  'du modèle repris : ',a,/,'est différent de celui Krot=',i3,/,
	2  'du modèle à poursuivre : ',a)
	  END SELECT
	  STOP
	 ELSE
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1023)wrot ; WRITE(2,1023)wrot
	  CASE DEFAULT
	   WRITE(*,23)wrot ; WRITE(2,23)wrot
	  END SELECT
	 ENDIF

c la rotation, avec un modèle en ASCII la rotation avec diffusion du moment
c cinétique ou conservation locale du moment cinétique
c n'est pas incluse dans le modèle repris (modèle repris en ASCII)
	 WRITE(*,17)thw(Krot) ; WRITE(2,17)thw(Krot)
	 SELECT CASE(Krot)
	 CASE(0,1,2)
	  ALLOCATE(rota(nrot,n_rot),mrot(n_rot),mrott(knotr))
	 CASE(3,4,5)
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1023)wrot ; WRITE(2,1023)wrot
	  CASE DEFAULT
	   WRITE(*,23)wrot ; WRITE(2,23)wrot
	  END SELECT
	 END SELECT

c détermination des abondances initiales (modèle repris en ASCII)
	 DEALLOCATE(comp) ; ALLOCATE(comp(nchim))
	 CALL nuc(1.5d+07,1.5d+02,comp,dcomp,jac,.FALSE.,1,
	1  epsilo,et,ero,ex,hh,be7,b8,n13,o15,f17)
c	 PRINT*,nom_elem ; PRINT*,nom_elemp ; PAUSE'elem'

c test de cohérence de la composition chimique
	 IF(nchim_ascii /= nchim)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1601)nchim_ascii,TRIM(chain_rep),nchim,TRIM(chain_don)
	   WRITE(2,1601)nchim_ascii,TRIM(chain_rep),nchim,TRIM(chain_don)
1601	   FORMAT('STOP, the number of elements nchim= ,',i3,/,
	1  'of the continued model in ascii : ,'a,/,
	2  'is different nchim= ,',i3,/,
	3  'from the one built from the data of the input file : ,'a)
	  CASE DEFAULT
	   WRITE(*,601)nchim_ascii,TRIM(chain_rep),nchim,TRIM(chain_don)
	   WRITE(2,601)nchim_ascii,TRIM(chain_rep),nchim,TRIM(chain_don)
601	   FORMAT('ARRET, le nombre d''éléments nchim= ,',i3,/,
	1  'du modèle repris en ASCII : ',a,/,
	2  'diffère de celui nchim= ,',i3,/,
	3  'du modèle à calculer avec les données du fichier : ',a)
	  END SELECT
	  STOP
	 ELSEIF(nchim_ascii /= nchim)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1602)nchim_ascii,TRIM(chain_rep),nchim,TRIM(chain_don)
	   WRITE(2,1602)nchim_ascii,TRIM(chain_rep),nchim,TRIM(chain_don)
1602	   FORMAT('STOP, the number of elements nchim= ,',i3,/,
	1  'of the continued model in ascii : ,'a,/,
	2  'is different nchim= ,',i3,/,
	3  'from the one built from the data of the input file : ,'a)
	  CASE DEFAULT
	   WRITE(*,602)nchim_ascii,TRIM(chain_rep),nchim,TRIM(chain_don)
	   WRITE(2,602)nchim_ascii,TRIM(chain_rep),nchim,TRIM(chain_don)
602	   FORMAT('ARRET, le nombre d''éléments nchim= ,',i3,/,
	1  'du modèle repris en ASCII : ',a,/,
	2  'diffère de celui nchim= ,',i3,/,
	3  'du modèle à calculer avec les données du fichier : ',a)
	  END SELECT
	  STOP
	 ELSEIF(ihe4_ascii /= ihe4)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1605)ihe4_ascii,TRIM(chain_rep),ihe4,TRIM(chain_don)
	   WRITE(2,1605)ihe4_ascii,TRIM(chain_rep),ihe4,TRIM(chain_don)
1605	   FORMAT('STOP, the subscript of He4 ihe4=,',i3,/,
	1  'of the continued model in ascii : ,'a,/,
	2  'is different ih4= ,',i3,/,
	3  'from the one built from the data of the input file : ,'a)
	  CASE DEFAULT
	   WRITE(*,605)ihe4_ascii,TRIM(chain_rep),ihe4,TRIM(chain_don)
	   WRITE(2,605)ihe4_ascii,TRIM(chain_rep),ihe4,TRIM(chain_don)
605	   FORMAT('ARRET, l''indice de He4 ihe4= ,',i3,/,
	1  'du modèle repris en ASCII : ',a,/,
	2  'diffère de celui ihe4= ,',i3,/,
	3  'du modèle à calculer avec les données du fichier : ',a)
	  END SELECT
	  STOP
	 ENDIF
	 ok=.TRUE.
	 DO i=1,nchim
	  ok=ok .AND. nom_elem(i) == nom_elem_ascii(i)
	 ENDDO
	 IF(.NOT.ok)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1011)nom_elem_ascii
	   WRITE(*,1603)TRIM(chain_rep),nom_elem
	   WRITE(*,1032)TRIM(chain_don)
	   WRITE(2,1011)nom_elem_ascii
	   WRITE(2,1603)TRIM(chain_rep),nom_elem
	   WRITE(2,1032)TRIM(chain_don)
1603	   FORMAT('read the input ASCII file : ',a,/,
	1  'differ from the chemical elements : ',30a)
	  CASE DEFAULT
	   WRITE(*,11)nom_elem_ascii
	   WRITE(*,603)TRIM(chain_rep),nom_elem
	   WRITE(*,32)TRIM(chain_don)
	   WRITE(2,11)nom_elem_ascii
	   WRITE(2,603)TRIM(chain_rep),nom_elem
	   WRITE(2,32)TRIM(chain_don)
603	   FORMAT('du modèle en ASCII repris : ',a,/,
	1  'diffèrent de ceux : ',30a)
	  END SELECT
	  STOP
	 ENDIF

c suite de la lecture du modèle ASCII
	 ALLOCATE(bp(ne,n),tds(1,n),chim(nchim,n)) ; bp=0.d0
	 DO k=1,n
	  READ(4,2005)bp(5,k),bp(1,k),bp(2,k),bp(3,k),bp(4,k),tds(1,k),
	1  chim(1:nchim,k)
2005	  FORMAT(5es22.15)
	 ENDDO
	 CLOSE(unit=4)

c initialisation de la spline d'inter. de la comp. chim.
c la composition chimique est tabulée par mole (modèle repris en ASCII)
	 DO i=1,nchim
	  chim(i,:)=chim(i,:)/nucleo(i)
	 ENDDO

c pour avoir un algorithme unique pour la diffusion, la composition
c chimique est toujours tabulée en fonction de mu=(m/Msol)^2/3
c que ce soit en lagrangien ou en eulérien
	 n_ch=n
	 ALLOCATE(mc(n_ch),mct(n_ch+m_ch))
	 mc(:)=bp(5,:)**(2.d0/3.d0)
	 CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,knotc,.FALSE.,mc(1),lq,
	1 xchim,dxchim)
	 IF(no_croiss)PRINT*,'Pb. at 1 in cesam'
	 dim_ch=knotc-m_ch

c si on utilise une grille fixe pour la composition chimique
	 IF(grille_fixe)THEN
	  ALLOCATE(mc_fixe(n_ch))
	  nc_fixe=n_ch ; mc_fixe=mc
	 ENDIF

c on passe en variables d'intégration
c m--> m**2/3-1, p--> ln p, t--> ln t, r--> r**2, l-->l**2/3
c ou bien m--> m, p--> ln p, t--> ln t, r--> r, l-->l
	 bp(1:2,1:n)=LOG(bp(1:2,1:n))	!pour Ptot et T
	 IF(pturb)bp(Ipg,:)=bp(1,:) 	!Pturb
	 IF(en_masse)THEN
	  bp(3,1:n)=bp(3,1:n)**2                 !pour R
	  bp(4:5,1:n)=bp(4:5,1:n)**(2.d0/3.d0)   !pour  L et M
	 ENDIF
	 bp(6,:)=1.d0          !pour psi (fictif)
	 ord_qsp=2       !initialisation sur snoein

c initialisation de la spline d'inter. de TdS sur m ou m^2/3
	 n_tds=n
	 ALLOCATE(x_tds(n_tds),xt_tds(n_tds+m_tds)) ; x_tds(:)=bp(5,:)
	 CALL bsp1dn(1,tds,x_tds,xt_tds,n_tds,m_tds,knot_tds,.FALSE.,
	1 x_tds(1),lq,ftds,dftds)
	 IF(no_croiss)PRINT*,'Pb. at 22 in cesam'
	 ALLOCATE(q(n),qt(n+ord_qsp)) ; q=(/ (i, i=1,n) /)
	 CALL bsp1dn(ne,bp,q,qt,n,ord_qsp,knot,.FALSE.,q(1),lq,fqs,dfqs)
	 IF(no_croiss)PRINT*,'Pb. at 23 in cesam'

c recherche de la fonction d'espacement (modèle repris en ASCII)
	 ALLOCATE(esp(n))
	 DO i=1,n
	  CALL bsp1dn(ne,bp,q,qt,n,ord_qsp,knot,.TRUE.,q(i),lq,fqs,dfqs)
	  IF(no_croiss)PRINT*,'Pb. at 2 in cesam'
	  esp(i)=ctep*fqs(1)+ctet*fqs(2)+cter*fqs(3)+ctel*fqs(4)+ctem*fqs(5)
c	  WRITE(*,2000)q(i),(f(j),j=1,ne),esp(i)
	  IF(i > 1)THEN
	   DO WHILE(esp(i) < esp(i-1))
	    esp(i)=esp(i)+0.01d0
	   ENDDO
	  ENDIF
	 ENDDO

c on dispose les ncouches pour assurer une répartition approx.
c uniforme de la fonction de répartition sur n_qs=n=n_ascii couches
	 n_qs=n
	 ALLOCATE(qb(n_qs))
	 CALL zoning(n,q,esp,n_qs,qb) !choix des nouveaux q dans qb

c esp est désormais inutile (modèle repris en ASCII)
	 DEALLOCATE(esp)

c dans bb on se place aux n_qs points qb
c on calcule la nouvelle fonction de répartition
	 ALLOCATE(bb(ne,n_qs),new_esp(1,n_qs))
	 DO i=1,n_qs
	  CALL bsp1dn(ne,bp,q,qt,n,ord_qsp,knot,.TRUE.,qb(i),lq,fqs,dfqs)
	  IF(no_croiss)PRINT*,'Pb. at 3 in cesam'
	  bb(1:ne,i)=fqs(1:ne)    !; WRITE(*,2000)bb(1:ne,i)
	  new_esp(1,i)=ctep*fqs(1)+ctet*fqs(2)+cter*fqs(3)+ctel*fqs(4)
	1 +ctem*fqs(5)
c	  WRITE(*,2000)q(i),(f(j),j=1,ne),new_esp(1,i)
	  IF(i > 1)THEN
	   DO WHILE(new_esp(1,i) < new_esp(1,i-1))
	    new_esp(1,i)=new_esp(1,i)+0.01d0
	   ENDDO
	  ENDIF
	 ENDDO
c	 PAUSE'new_esp'

c spline d'ordre 2 sur les nouveaux qb sur n_qs points pour new_esp
	 qb=(/ (i, i=1,n_qs) /)
	 ALLOCATE(qbt(n_qs+ord_qsp))
	 CALL bsp1dn(1,new_esp,qb,qbt,n_qs,2,knotb,.FALSE.,qb(1),lq,fqs,dfqs)
	 IF(no_croiss)PRINT*,'Pb. at 24 in cesam'

c dérivée de la fonction de répartition dans bb(6,:) (modèle repris en ASCII)
	 DO i=1,n_qs
	  CALL bsp1dn(1,new_esp,qb,qbt,n_qs,2,knotb,.TRUE.,qb(i),lq,fqs,dfqs)
	  IF(no_croiss)PRINT*,'Pb. at 4 in cesam'
	  bb(6,i)=dfqs(1)
	 ENDDO

c on spline bb (modèle repris en ASCII)
	 CALL bsp1dn(ne,bb,qb,qbt,n_qs,2,knotb,.FALSE.,qb(1),lq,fqs,dfqs)
	 IF(no_croiss)PRINT*,'Pb. at 25 in cesam'
c	 DO i=1,n_qs
c	  CALL bsp1dn(ne,bb,qb,qbt,n_qs,2,knotb,.TRUE.,qb(i),lq,fqs,dfqs)
c	  WRITE(*,2000)qb(i),fqs
c	 ENDDO
c	 PAUSE'bb'

c dim. de la base de splines pour l'intégration equi. quasi. stat.
	 dim_qs=(n_qs-1)*m_qs+r_qs ; ord_qs=m_qs+r_qs

c vecteur nodal pour l'intégration sur n_qs points
	 DEALLOCATE(q,qt)
	 ALLOCATE(q(n_qs),qt(dim_qs+ord_qs)) ; q=(/ (i, i=1,n_qs) /)
	 CALL noedif(q,n_qs,m_qs,r_qs,qt,knot)
	 IF(no_croiss)THEN
	  PRINT*,'ARRET 6 dans cesam / STOP 6 in cesam' ; STOP
	 ENDIF
c	 PRINT*,ord_qs,knot,dim_qs ; WRITE(*,2000)qt(1:knot)
c	 PAUSE'noedif'

c on place la solution initiale dans la base de q, bb ==> bp
	 DEALLOCATE(bp) ; ALLOCATE(bp(ne,dim_qs))
	 CALL newspl(ne,qb,qbt,knotb,2,qt,knot,ord_qs,bb,bp)
c	 PAUSE'newspl'

c bb, qb, qbt sont désormais inutiles
	 DEALLOCATE(bb,qb,qbt)

c au centre R=L=M=0
	 bp(3:5,1)=0.d0

c vérification
c	 DO i=1,n_qs
c	  CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),lq,fqs,dfqs)
c	  WRITE(*,2000)q(i),fqs
c	 ENDDO
c	 PAUSE'vérification'

c initialisations et allocations diverses
	 mstar=mtot

c tabulation de r2 et m23 qui accélèrent la recherche dans le ssp.
c inter valable pour variables lagrangiennes alors m23=m^2/3, r2=r^2
c ou eulériennes alors m23=m, r2=r
	 ALLOCATE(r2(n_qs),m23(n_qs))
	 DO i=1,n_qs
c(test)   CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),lq,fqs,dfqs)
	  r2(i)=bp(3,m_qs*(i-1)+1) ; m23(i)=bp(5,m_qs*(i-1)+1)
	 ENDDO

c estimation provisoire de rstar (modèle repris en ASCII)
	 rstar=r2(n_qs) ; IF(en_masse)rstar=SQRT(rstar)

c vérification concernant la vitesse angulaire
c pour la couche externe : force centifuge < force de gravité / 2.
	 IF(wrot == 0.d0)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1083) ; WRITE(2,1083)
1083	   FORMAT('without rotation nor loss of angular momentum.')
	  CASE DEFAULT
	   WRITE(*,83) ; WRITE(2,83)
83	   FORMAT('sans rotation ni perte de moment. cinétique.')
	  END SELECT
	 ELSE
	  CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(n_qs),lq,fqs,dfqs)
	  IF(no_croiss)PRINT*,'Pb. at 5 in cesam'
	  IF(en_masse)THEN
	   fqs(5)=SQRT(fqs(5))**3 ; fqs(3)=SQRT(fqs(3))
	  ENDIF
	  fqs(5)=fqs(5)*msol ; fqs(3)=fqs(3)*rsol
	  w_max=SQRT(g*fqs(5)/fqs(3)**3)/2.d0

c si la rotation est entrée en kms/s w = v / R (modèle repris en ASCII)
	  IF(unit == 'kms/s')THEN
	   wrot=wrot*1.d5/fqs(3) ; w_rot=wrot
	  ENDIF
	  IF(wrot > w_max)THEN
	  SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1400)wrot,w_max,TRIM(chain_don)
1400	    FORMAT('STOP, the centrifugal force is too large w_rot=',
	1   es10.3,' > w_max=',es10.3,/,
	2   'decrease w_rot in the input file : ',a)
	   CASE DEFAULT
	    WRITE(*,400)wrot,w_max,TRIM(chain_don)
400	    FORMAT('ARRET, la force centrifuge est trop grande w_rot=',
	1   es10.3,' > w_max=',es10.3,/,
	2   'réduire w_rot dans le fichier de données : ',a)
	   END SELECT
	   STOP
	  ENDIF

c rot_min : valeur négligeable, initialisation de rota
	  IF(Krot >= 3)THEN
	   rot_min(1)=wrot*1.d-2 ; CALL initialise_rota
	  ENDIF
	 ENDIF

	 DEALLOCATE(nom_elem_ascii)

c gestion du premier pas temporel
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1058)dtp
1058	  FORMAT('the time step of the continued model in ASCII, dtp=',
	1 es10.3,/,'is recommanded for the first time step, OK?',/,
	2 'enter o for yes , enter n for no')
	 CASE DEFAULT
	  WRITE(*,58)dtp
58	  FORMAT('le pas temporel du modèle repris en ASCII, dtp=',es10.3,/,
	1 'est conseillé pour le premier pas temporel, ok? o/n')
	 END SELECT
	 READ*,oui
	 IF(oui == 'o')THEN
	  dt=MIN(dtp,dtmax) ; dts=dt
	 ELSE
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1059)
1059	   FORMAT('enter the new time step, in units of 10**6 years')
	  CASE DEFAULT
	   WRITE(*,59)
59	   FORMAT('entrer le nouveau pas temporel, unité : 10**6 ans')
	  END SELECT
	  READ*,dt
	  IF(dt > dtmax)THEN
	   dt=dtmax ; dts=dt
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1112)dtmax
	   CASE DEFAULT
	    WRITE(*,112)dtmax
	   END SELECT
	  ENDIF
	 ENDIF           !oui /= 'o'
	 IF(age+dt > agemax)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1053)age+dt,agemax
	  CASE DEFAULT
	   WRITE(*,53)age+dt,agemax
	  END SELECT
	  dt=MIN(dtmax,MAX(0.d0,agemax-age))
	  IF(dt <= 0.d0)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1054) ; WRITE(2,1054)
	   CASE DEFAULT
	    WRITE(*,54) ; WRITE(2,54)
	   END SELECT
	   READ*,oui
	   IF(oui == 'o')THEN
	    list_sort=.TRUE.
	   ELSE
	    STOP
	   ENDIF
	  ENDIF   !dt == 0.d0
	 ELSE
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1051)dt ; WRITE(2,1051)dt
	  CASE DEFAULT
	   WRITE(*,51)dt ; WRITE(2,51)dt
	  END SELECT
	  list_rep=.TRUE.   !liste du modèle
	 ENDIF           !age+dt > agemax

c initialisation des paramètres pour le vent et les chutes de planétoïdes
	 CALL vent ; CALL planetoides

c initialisation avec un modèle de ZAMS
	CASE(-2, 2)

c modèle d'initialisation en binaire/ASCII, modèle de ZAMS
c pour éviter de prendre un modèle d'atmosphère trop différent d'une
c atmosphère de ZAMS, suppression du modèle d'atmosphère s'il existe
	 INQUIRE(file=TRIM(chain_atm),exist=ok)
	 IF(ok)CALL system('rm '//TRIM(chain_atm))

	 IF(un23 > 0)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1008)nom_fich1 ; WRITE(2,1008)nom_fich1
1008	   FORMAT('One chooses as a binary initial model : ',a)
	  CASE DEFAULT
	   WRITE(*,8)nom_fich1 ; WRITE(2,8)nom_fich1
8	   FORMAT('On prend pour modèle initial en binaire : ',a)
	  END SELECT
	  OPEN(unit=4,form='unformatted',status='old',file=nom_fich1)
	  READ(4)nep,m_qsp,n,knot,nchimp
c	  PRINT*,nep,m_qsp,n,knot,nchimp
	  ord_qsp=m_qsp+r_qs ; dim_qsp=knot-ord_qsp
	  ALLOCATE(b(nep,dim_qsp),q(n),qt(knot),nom_elemp(nchimp))
	  REWIND(unit=4)
	  READ(4)nep,m_qsp,n,knot,nchimp,n_ch,m_chp,knotc,Krotp,nrotp,
	1 n_rot,m_rotp,knotr,n_tds,knot_tds,mtotp,alphap,w_rotp,lim_rop,
	2 diffusionp,rot_solidp,precisionp,en_massep,f_eosp,f_opap,
	3 nom_ctesp,nom_pertmp,nom_pertwp,nom_tdetaup,nom_atmp,nom_convp,
	4 nom_nucp,nom_nuc_cplp,nom_diffmp,nom_difftp,nom_diffwp,
	5 nom_etatp,nom_opap,nom_elemp,b,q,qt
	  CLOSE(unit=4)

c	  PRINT*,n,ord_qsp,knot,nep,n_tds,knot_tds
c	  WRITE(*,2000)b(1,1),q(1),qt(knot) ; PAUSE'lec'

c	  DO i=1,n
c	   CALL bsp1dn(nep,b,q,qt,n,ord_qsp,knot,.TRUE.,q(i),lq,fqs,dfqs)
c	   WRITE(*,2000)q(i),fqs
c	   PRINT*,i,fqs(5)
c	  ENDDO
c	  PAUSE'test'

c on met modèle repris dans bp en l'adaptant, modèle de ZAMS
	  ALLOCATE(bp(ne,dim_qsp))

c suppression/adjonction de la pression turbulente
c si le modèle repris est avec/sans Pturb
c et que le modèle à calculer est sans/avec Pturb
c il faut enlever/ajouter la variable ln Pgaz
c ne/nep : nb. inconnues du modèle à calculer/repris
	  IF(ne == nep)THEN
	   bp=b     	!on a le même nombre d'inconnues
          ELSE		!changement du nb. inconnues
	   IF(ne == 6)THEN  !suppression de Pgaz on passe de ne=7 à 6 inconnues
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1721) ; WRITE(2,721)
	    CASE DEFAULT
	     WRITE(*,721) ; WRITE(2,721)
	    END SELECT
	    bp=b(1:6,:)
	   ELSE      !adjonction de Pgaz on passe de 6 à 7 inconnues
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1722) ; WRITE(2,722)
	    CASE DEFAULT
	     WRITE(*,722) ; WRITE(2,722)
	    END SELECT
	    bp(1:6,:)=b ; bp(Ipg,:)=b(1,:)
	   ENDIF     !passage de 6 <---> 7 inconnues
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1723)nep ; WRITE(2,1723)ne
	   CASE DEFAULT
	    WRITE(*,723)nep ; WRITE(2,723)ne
	   END SELECT
	  ENDIF     !transfert et changement du nb. inconnues

c b est désormais inutile
	  DEALLOCATE(b)

c on adapte le type variables du modèle repris à celles
c du modèle à calculer, modèle de ZAMS
c	  PRINT*,en_massep,en_masse
	  IF(en_massep .NEQV. en_masse)THEN
	   IF(en_massep)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1225) ; WRITE(2,1225)
1225	     FORMAT('change of variables: from Lagrangean to Eulerian',/)
	    CASE DEFAULT
	     WRITE(*,225) ; WRITE(2,225)
225	     FORMAT('chgt. de variables : lagrangien --> eulérien',/)
	    END SELECT
	    bp(3,:)=SQRT(ABS(bp(3,:)))      !R**2   ==> R
	    bp(4,:)=SQRT(ABS(bp(4,:)))**3   !L**2/3 ==> L
	    bp(5,:)=SQRT(ABS(bp(5,:)))**3   !M**2/3 ==> M
	   ELSE
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1226) ; WRITE(2,1226)
1226	     FORMAT('change of variables: from Eulerian to Lagrangean',/)
	    CASE DEFAULT
	     WRITE(*,226) ; WRITE(2,226)
226	     FORMAT('chgt. de variables : eulérien --> lagrangien',/)
	    END SELECT
	    bp(3:5,:)=bp(3:5,:)**2          !R, L, M ==> (R, L, M)**2
	    bp(4,:)=bp(4,:)**(1.d0/3.d0)    !L**2 ==> L**2/3
	    bp(5,:)=bp(5,:)**(1.d0/3.d0)    !M**2 ==> M**2/3
	   ENDIF
	   
	  ELSE
	   SELECT CASE(langue)
	   CASE('english')	   
	     WRITE(*,1227)en_massep,en_masse ; WRITE(2,1227)en_massep,en_masse
1227	     FORMAT('no change of variables: en_massep: ',l1,', en_masse: ',
	1    l1,/)
	    CASE DEFAULT
	     WRITE(*,227)en_massep,en_masse ; WRITE(2,227)en_massep,en_masse
227	     FORMAT('pas de chgt. de variables, en_massep: ',l1,', en_masse: ',
	1    l1,/)
	    END SELECT
	  ENDIF

c	  DO i=1,n
c	  CALL bsp1dn(ne,bp,q,qt,n,ord_qsp,knot,.TRUE.,q(i),lq,fqs,dfqs)
c	   WRITE(*,2000)fqs
c	  ENDDO
c	  bid=bp(5,1) ; lq=ord_qsp
c	  DO i=2,n
c	   CALL bsp1dn(ne,bp,q,qt,n,ord_qsp,knot,.TRUE.,q(i),lq,fqs,dfqs)
c	   IF(bp(5,i-1) >= bp(5,i))PRINT*,i,bp(5,i-1:i)	   
c	   IF(bid >= fqs(5))PRINT*,i,bid,fqs(5),bp(5,i-1:i)
c	   bid=fqs(5)
c	  ENDDO	    	   
c	  PAUSE'reprise'

	 ELSE

c modèle d'initialisation en ASCII de ZAMS
	  OPEN(unit=4,form='formatted',status='old',file=nom_fich1)

c détermination du nombre de couches
	  n=0
	  DO
	   READ(4,*,END=15)fqs(1) ; n=n+1
	  ENDDO
15	  REWIND(unit=4)
c	  PRINT*,n ; WRITE(*,2000)mtot ; PAUSE'reprise ASCII'

c les masses sont en 1-M/Mtot
	  ALLOCATE(bp(ne,n)) ; bp=0.d0
	  DO i=1,n
	   READ(4,*)bp(5,i),bp(1:4,i) ; bp(5,i)=(1.d0-bp(5,i))*mtot
	  ENDDO
	  CLOSE(unit=4)		! PAUSE'lect'
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1016)nom_fich1,n ; WRITE(2,1016)nom_fich1,n
1016	   FORMAT('initializing model in ASCII : ',a,/,
	1  'number of shells : ',i4)
	  CASE DEFAULT
	   WRITE(*,16)nom_fich1,n ; WRITE(2,16)nom_fich1,n
16	   FORMAT('modèle d''initialisation en ASCII : ',a,/,
	1  'nombre de couches :',i4)
	  END SELECT

c retournement 1 au centre, n à la surface
	  bp(:,1:n:+1)=bp(:,n:1:-1)

c on passe en variables d'intégration
c m--> m**2/3-1, p--> ln p, t--> ln t, r--> r**2, l-->l**2/3
c ou bien m--> m, p--> ln p, t--> ln t, r--> r, l-->l
	  bp(1:2,1:n)=LOG(bp(1:2,1:n))            !pour Ptot et T
	  IF(pturb)bp(Ipg,:)=bp(1,:)	          !Pturb
	  IF(en_masse)THEN
	   bp(3,1:n)=bp(3,1:n)**2                 !pour R
	   bp(4:5,1:n)=ABS(bp(4:5,1:n))**(2.d0/3.d0)   !pour  L et M
	  ENDIF
	  bp(6,:)=1.d0	       !pour psi (fictif)
	  ord_qsp=2	 !initialisation sur snoein
	  ALLOCATE(q(n),qt(n+ord_qsp)) ; q=(/ (i, i=1,n) /)
	  CALL bsp1dn(ne,bp,q,qt,n,ord_qsp,knot,.FALSE.,q(1),lq,fqs,dfqs)
	  IF(no_croiss)PRINT*,'Pb. at 26 in cesam'	  
	 ENDIF      !initialisation binaire/ASCII

c on fixe le nb. de couches du modèle d'initialisation à ncouche=300
c 300 étant un bon ordre de grandeur pour un modèle solaire
c ce nombre sera ensuite ajusté de facon à ce que dQ/dq = psi0
c on détermine les ncouche points de masse m(ncouche) pour
c assurer une répartition approximativement uniforme en
c ctep lnP + ctet lnT + cter zeta + ctel lambda + ctem mu

c recherche de la fonction d'espacement, modèle de ZAMS
c le modèle d'initialisation est dans bp(ne,n) si on a repris
c en ASCII ou dans b(ne,dim_qsp) si on a repris en binaire, comme
c on fait qu'interpoler le fait que n puisse être différent de
c dim_qsp n'a pas d'importance
	 ALLOCATE(esp(n))
	 DO i=1,n
	  CALL bsp1dn(ne,bp,q,qt,n,ord_qsp,knot,.TRUE.,q(i),lq,fqs,dfqs)
	  IF(no_croiss)PRINT*,'Pb. at 6 in cesam'
	  esp(i)=ctep*fqs(1)+ctet*fqs(2)+cter*fqs(3)+ctel*fqs(4)+ctem*fqs(5)	  
c	  WRITE(*,2000)q(i),fqs,esp(i)
	  IF(i > 1)THEN
	   DO WHILE(esp(i) < esp(i-1))
	    esp(i)=esp(i)+0.01d0
	   ENDDO
	  ENDIF
	 ENDDO

c	 PRINT*,'ord_qsp,knot,n',ord_qsp,knot,n
c	 WRITE(*,2000)ctel,ctep,ctet,ctem,cter ; WRITE(*,2000)esp
c	 WRITE(*,2000)q ; WRITE(*,2000)qt ; PAUSE'esp'

c on dispose les ncouches pour assurer une répartition approx.
c uniforme de la fonction de répartition sur n_qs couches	 
	 n_qs=n_init ; ALLOCATE(qb(n_qs))
	 	 
	 CALL zoning(n,q,esp,n_qs,qb) !choix des nouveaux q dans qb

c	 PRINT*,'anciens q',n ; WRITE(*,2000)q(1:n) ; PRINT*,'esp'
c	 WRITE(*,2000)esp(1:n) ; PRINT*,'new q' ; WRITE(*,2000)qb(1:n_qs)
c	 PAUSE'après zoning'

c esp est désormais inutile
	 DEALLOCATE(esp)

c dans bb on se place aux n_qs points qb
c on calcule la nouvelle fonction de répartition, modèle de ZAMS
	 ALLOCATE(bb(ne,n_qs),new_esp(1,n_qs))
	 DO i=1,n_qs
	  CALL bsp1dn(ne,bp,q,qt,n,ord_qsp,knot,.TRUE.,qb(i),lq,fqs,dfqs)
	  IF(no_croiss)PRINT*,'Pb. at 7 in cesam'
	  bb(1:ne,i)=fqs(1:ne)    !; WRITE(*,2000)bb(1:ne,i)
	  new_esp(1,i)=ctep*fqs(1)+ctet*fqs(2)+cter*fqs(3)+ctel*fqs(4)
	1 +ctem*fqs(5)
	  IF(i > 1)THEN
	   DO WHILE(new_esp(1,i) < new_esp(1,i-1))
	    new_esp(1,i)=new_esp(1,i)+0.01d0
	   ENDDO
	  ENDIF
	 ENDDO
c	 PAUSE'new_esp'

c spline d'ordre 2 sur les nouveaux qb sur n_qs points pour new_esp
	 ALLOCATE(qbt(n_qs+ord_qsp)) ; qb=(/ (i, i=1,n_qs) /)
	 CALL bsp1dn(1,new_esp,qb,qbt,n_qs,2,knotb,.FALSE.,qb(1),lq,fqs,dfqs)
	 IF(no_croiss)PRINT*,'Pb. at 27 in cesam'

c dérivée de la fonction de répartition dans bb(6,:)
	 DO i=1,n_qs
	  CALL bsp1dn(1,new_esp,qb,qbt,n_qs,2,knotb,.TRUE.,qb(i),lq,fqs,dfqs)
	  IF(no_croiss)PRINT*,'Pb. at 8 in cesam'
	  bb(6,i)=dfqs(1)
	 ENDDO

c on spline bb
	 CALL bsp1dn(ne,bb,qb,qbt,n_qs,2,knotb,.FALSE.,qb(1),lq,fqs,dfqs)
	 IF(no_croiss)PRINT*,'Pb. at 28 in cesam'

c vérification
c	 DO i=1,n_qs
c	 CALL bsp1dn(ne,bb,qb,qbt,n_qs,2,knotb,.TRUE.,qb(i),lq,fqs,dfqs)
c	  WRITE(*,2000)qb(i),fqs(1:5) !; PRINT*,qb(i),fqs(5)
c	 ENDDO
c	 PAUSE'bb'

c dimension de la base de splines pour l'intégration equi. quasi. stat.,
c modèle de ZAMS
	 dim_qs=(n_qs-1)*m_qs+r_qs ; ord_qs=m_qs+r_qs

c vecteur nodal pour l'intégration sur n_qs points
	 DEALLOCATE(q,qt)
	 ALLOCATE(q(n_qs),qt(dim_qs+ord_qs)) ; q=(/ (i, i=1,n_qs) /)
	 CALL noedif(q,n_qs,m_qs,r_qs,qt,knot)
	 IF(no_croiss)THEN
	  PRINT*,'ARRET 10 dans cesam / STOP 10 in cesam' ; STOP
	 ENDIF

c	 PRINT*,ord_qs,knot,dim_qs ; WRITE(*,2000)qt(1:knot)
c	 PAUSE'noedif'

c on place la solution initiale dans la base de q, bb ==> bp, modèle de ZAMS
	 DEALLOCATE(bp) ; ALLOCATE(bp(ne,dim_qs))
	 CALL newspl(ne,qb,qbt,knotb,2,qt,knot,ord_qs,bb,bp)
c	 PAUSE'newspl'

c bb, qb, qbt sont désormais inutiles
	 DEALLOCATE(bb,qb,qbt)

c au centre R=L=M=0
	 bp(3,1)=0.d0 ; bp(4,1)=0.d0 ; bp(5,1)=0.d0

c vérification
c	 DO i=1,n_qs
c	  CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),lq,fqs,dfqs)
c	  WRITE(*,2000)q(i),f(1:ne)
c	  PRINT*,i,fqs(5)
c	 ENDDO
c	 PRINT*,'ord_qs=',ord_qs
c	 PAUSE'vérification'

c modèle de ZAMS
c initialisation du vecteur de composition chimique
c pour la couche externe : force centifuge < force de gravité / 2.
	 IF(wrot == 0.d0)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1083) ; WRITE(2,1083)
	  CASE DEFAULT
	   WRITE(*,83) ; WRITE(2,83)
	  END SELECT

	 ELSE
	  CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(n_qs),lq,fqs,dfqs)
	  IF(no_croiss)PRINT*,'Pb. at 9 in cesam'
	  IF(en_masse)THEN
	   fqs(5)=SQRT(fqs(5))**3 ; fqs(3)=SQRT(fqs(3))
	  ENDIF
	  fqs(5)=fqs(5)*msol ; fqs(3)=fqs(3)*rsol
	  w_max=SQRT(g*fqs(5)/fqs(3)**3)/2.d0

c si la rotation a été entrée en kms/s w = v / R
	  IF(unit == 'kms/s')THEN
	   wrot=wrot*1.d5/fqs(3) ; w_rot=wrot
	  ENDIF
	  IF(wrot > w_max)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1400)wrot,w_max,TRIM(chain_don)
	   CASE DEFAULT
	    WRITE(*,400)wrot,w_max,TRIM(chain_don)
	   END SELECT
	   STOP
	  ENDIF

c rot_min : valeur négligeable
	  IF(Krot >= 3)rot_min(1)=wrot*1.d-2
	 ENDIF

c appel d'initialisation pour tabulation des réactions nucléaires
c allocations fictives, modèle de ZAMS
	 ALLOCATE(comp(0),dcomp(0),jac(0,0),ex(0))
	 CALL nuc(1.5d+07,1.5d+02,comp,dcomp,jac,.FALSE.,0,
	1 epsilo,et,ero,ex,hh,be7,b8,n13,o15,f17)

c détermination des abondances initiales, modèle de ZAMS
	 DEALLOCATE(comp) ; ALLOCATE(comp(nchim))
	 CALL nuc(1.5d+07,1.5d+02,comp,dcomp,jac,.FALSE.,1,
	1 epsilo,et,ero,ex,hh,be7,b8,n13,o15,f17)

c une trop grande abondance initiale de Lithium est incompatible
c avec un modèle de ZAMS homogène, modèle de ZAMS
	 Bli: DO i=1,nchim
	  IF(nom_elem(i)(1:2) == 'Li')THEN
	   IF(comp(i) > 1.d-15)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1014)nom_elem(i),comp(i),EXP(ln_Tli)
	     WRITE(2,1014)nom_elem(i),comp(i),EXP(ln_Tli)
1014	     FORMAT('With an initial abundance of ',a4,'=',es10.3,/,
	1    'the ZAMS might not be able to converge',/,
	2    'does one try without lithium for t >',es10.3,
	3    'enter o for yes, enter n for no')
	    CASE DEFAULT
	     WRITE(*,14)nom_elem(i),comp(i),EXP(ln_Tli)
	     WRITE(2,14)nom_elem(i),comp(i),EXP(ln_Tli)
14	     FORMAT('Avec une abondance initiale de ',a4,'=',es10.3,/,
	1    'le modèle de ZAMS risque de ne pas converger',/,
	2    'tente-on en supprimant le lithium pour T >',es10.3,/,
	3    'entrer o/n.')
	    END SELECT
	    READ*,oui
	    IF(oui == 'o')THEN
	     SELECT CASE(langue)
	     CASE('english')
	      WRITE(*,1025)EXP(ln_Tli) ; WRITE(2,1025)EXP(ln_Tli)
1025	      FORMAT('one tries, without lithium for T > ',es10.3)
	     CASE DEFAULT
	      WRITE(*,25)EXP(ln_Tli) ; WRITE(2,25)EXP(ln_Tli)
25	      FORMAT('on tente, en supprimant le lithium pour T > ',es10.3)
	     END SELECT
	     EXIT Bli
	    ELSE
	     SELECT CASE(langue)
	     CASE('english')
	      WRITE(*,1026) ; WRITE(2,1026)
1026	      FORMAT('STOP, it is suggested :',/,
	1     'either to compute the PMS or to initialize Li7 < 1.d-15')
	     CASE DEFAULT
	      WRITE(*,26) ; WRITE(2,26)
26	      FORMAT('ARRET, suggestions :',/,
	1     '- calculer la PMS ou encore initialiser Li7 < 1.d-15')
	     END SELECT
	     STOP
	    ENDIF
	   ENDIF
	  ENDIF
	 ENDDO Bli

c la rotation, modèle de ZAMS
	 WRITE(*,17)thw(Krot) ; WRITE(2,17)thw(Krot)

	 SELECT CASE(Krot)
	 CASE(0,1,2)
	  nrot=0 ; n_rot=0 ; knotr=0
	  ALLOCATE(rota(nrot,n_rot),mrot(n_rot),mrott(knotr))
	 END SELECT
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1023)wrot ; WRITE(2,1023)wrot
	 CASE DEFAULT
	  WRITE(*,23)wrot ; WRITE(2,23)wrot
	 END SELECT

c pour avoir un algorithme unique pour la diffusion, la composition
c chimique est toujours tabulée en fonction de mu=(m/Msol)^2/3
c que ce soit en lagrangien ou en eulérien
	 n_ch=n_qs
	 ALLOCATE(chim(nchim,n_ch),mc(n_ch))
	 DO i=1,n_qs
	  CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),lq,fqs,dfqs)
	  IF(no_croiss)PRINT*,'Pb. at 10 in cesam'
c	  WRITE(*,2000)fqs
	  fqs(5)=MAX(fqs(5),0.d0)
	  IF(en_masse)THEN
	   mc(i)=fqs(5)                 !interpolation en m**2/3
	  ELSE
	   mc(i)=ABS(fqs(5))**(2.d0/3.d0)
	  ENDIF
	  chim(:,i)=comp(:)

c sur la ZAMS, annulation de l'abondance initiale du lithium
	  IF(fqs(2) > ln_Tli)THEN
	   DO j=1,nchim
	    IF(nom_elem(j)(1:2) == 'Li')chim(j,i)=li_ini
	   ENDDO
	  ENDIF
	 ENDDO

c initialisation de la spline d'inter. de la comp. chim., modèle de ZAMS
	 ALLOCATE(dxchim(nchim),mct(n_ch+m_ch),xchim(nchim))
	 CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,knotc,.FALSE.,mc(1),lq,
	1 xchim,dxchim)
	 IF(no_croiss)PRINT*,'Pb. at 29 in cesam'
	 dim_ch=knotc-m_ch ; DEALLOCATE(dxchim,xchim)

c initialisations et allocations diverses
	 mstar=mtot

c tabulation de r2 et m23 qui accélèrent la recherche dans le ssp.
c inter valable pour variables lagrangiennes alors m23=m^2/3, r2=r^2
c ou eulériennes alors m23=m, r2=r
	 ALLOCATE(r2(n_qs),m23(n_qs))
	 DO i=1,n_qs
c(test)   CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),lq,fqs,dfqs)
	  r2(i)=bp(3,m_qs*(i-1)+1) ; m23(i)=bp(5,m_qs*(i-1)+1)
	 ENDDO

c initialisation de rota, modèle de ZAMS
	 IF(Krot >= 3)CALL initialise_rota

c si on utilise une grille fixe pour la composition chimique
	 IF(grille_fixe)THEN
	  ALLOCATE(mc_fixe(n_qs)) ; nc_fixe=n_qs ; mc_fixe=m23
	  IF(.NOT.en_masse)mc_fixe=mc_fixe**(2.d0/3.d0)
         ENDIF

c estimation provisoire de rstar
	 rstar=r2(n_qs) ; IF(en_masse)rstar=SQRT(rstar)

c	 PAUSE'avant call resout'

c initialisation de l'âge et du pas temporel pour le modèle ZAMS
c calcul du modèle de ZAMS homogène
c pour écriture sur (26), on impose dt=0 (et non pas dt_minimum)
	 age=0.d0 ; dt=0.d0 ; CALL resout(2,dt,dts) ; dt=0.d0 ; dts=0.d0

c	 WRITE(*,2000)bp(6,1) ; PAUSE'bp(6,1)'

c on garde le modèle de ZAMS d'âge 0 sur le fichier xxxx_B.hom
	 model_num=0
	 IF(all_rep)THEN
	  WRITE(number,70)model_num
70	  FORMAT(i4.4)
	  modelb=TRIM(nom_fich2)//number//'_B.hom'
	 ELSE
	  modelb=TRIM(nom_fich2)//'_B.hom'
	 ENDIF
	 OPEN(unit=26,form='unformatted',status='unknown',file=TRIM(modelb))
	 WRITE(26)ne,m_qs,n_qs,knot,nchim,n_ch,m_ch,knotc,Krot,nrot,n_rot,
	1 m_rot,knotr,n_tds,knot_tds,mtot,alpha,w_rot,lim_ro,diffusion,
	2 rot_solid,precision,en_masse,f_eos,f_opa,
	3 nom_ctes,nom_pertm,nom_pertw,nom_tdetau,nom_atm,nom_conv,nom_nuc,
	4 nom_nuc_cpl,nom_diffm,nom_difft,nom_diffw,nom_etat,nom_opa,
	5 nom_elem,bp,q,qt,chim,mc,mct,rota,mrot,mrott,tds,x_tds,xt_tds,
	6 m23,r2,m_zc,r_zc,r_ov,age,dt,dts,mstar,rstar,mw_tot,wrot,
	7 jlim,lconv,lim,model_num
	 CLOSE(unit=26)
c	 PRINT*,all_rep,model_num,modelb ; PAUSE'.hom'

c le numéro du modèle sera augmenté d'une unité ligne ~4155
	 model_num=-1

c pas temporel initial pour évolution à partir de la ZAMS
	 dtp=0.d0 ; dt=MIN(dt0,dtmax) ; dts=dt

c pour l'entête du modèle du fichier mon_modèle.lis
	 list_zams=.TRUE.

c initialisation des paramètres pour le vent et les chutes de planétoïdes
	 CALL vent ; CALL planetoides

c initialisation d'un modèle PMS
	CASE(-3, 3)

c pour éviter de prendre un modèle d'atmosphère trop différent d'une
c atmosphère de PMS, suppression du modèle d'atmosphère s'il existe
	 INQUIRE(file=TRIM(chain_atm),exist=ok)
	 IF(ok)CALL system('rm '//TRIM(chain_atm))

c modèle d'initialisation en binaire/ASCII
	 IF(un23 > 0)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1008)nom_fich1 ; WRITE(2,1008)nom_fich1
	  CASE DEFAULT
	   WRITE(*,8)nom_fich1 ; WRITE(2,8)nom_fich1
	  END SELECT
	  OPEN(unit=4,form='unformatted',status='old',file=nom_fich1)
	  READ(4)nep,m_qsp,n,knot,nchimp
c	  PRINT*,nep,m_qsp,n,knot,nchimp
	  ord_qsp=m_qsp+r_qs ; dim_qsp=knot-ord_qsp
	  ALLOCATE(b(nep,dim_qsp),q(n),qt(knot),nom_elemp(nchimp))
	  REWIND(unit=4)
	  READ(4)nep,m_qsp,n,knot,nchimp,n_ch,m_chp,knotc,Krotp,nrotp,
	1 n_rot,m_rotp,knotr,n_tds,knot_tds,mtotp,alphap,w_rotp,lim_rop,
	2 diffusionp,rot_solidp,precisionp,en_massep,f_eosp,f_opap,
	3 nom_ctesp,nom_pertmp,nom_pertwp,nom_tdetaup,nom_atmp,nom_convp,
	4 nom_nucp,nom_nuc_cplp,nom_diffmp,nom_difftp,nom_diffwp,
	5 nom_etatp,nom_opap,nom_elemp,b,q,qt
	  CLOSE(unit=4)

c	  PRINT*,n,ord_qsp,knot,nep,knot_tds
c	  WRITE(*,2000)b(1,1),q(1),qt(knot) ; PAUSE'lec'

c	  DO i=1,n
c	   CALL bsp1dn(nep,b,q,qt,n,ord_qsp,knot,.TRUE.,q(i),lq,fqs,dfqs)
c	   WRITE(*,2000)q(i),fqs ; PAUSE'1'
c	  ENDDO
c	  PAUSE'test'

c on met modèle repris dans bp en l'adaptant
	  ALLOCATE(bp(ne,dim_qsp))

c suppression/adjonction de la pression turbulente
c si le modèle repris est avec/sans Pturb
c et que le modèle à calculer est sans/avec Pturb
c il faut enlever/ajouter la variable ln Pgaz
c ne/nep : nb. inconnues du modèle à calculer/repris

	  IF(ne == nep)THEN
	   bp=b                 !on a le même nombre d'inconnues
	  ELSE	!changement du nb. inconnues
	   IF(ne == 6)THEN  !suppression de Pgaz on passe de ne=7 à 6 inconnues
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1721) ; WRITE(2,721)
	    CASE DEFAULT
	     WRITE(*,721) ; WRITE(2,721)
	    END SELECT
	    bp=b(1:7,:)
	   ELSE      !adjonction de Pgaz on passe de 6 à 7 inconnues
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1722) ; WRITE(2,722)
	    CASE DEFAULT
	     WRITE(*,722) ; WRITE(2,722)
	    END SELECT
	    bp(1:6,:)=b ; bp(Ipg,:)=b(1,:)
	   ENDIF     !passage de 6 <---> 7 inconnues
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1723)nep ; WRITE(2,1723)ne
	   CASE DEFAULT
	    WRITE(*,723)nep ; WRITE(2,723)ne
	   END SELECT
	  ENDIF     !transfert et changement du nb. inconnues

c b est désormais inutile
	  DEALLOCATE(b)

c on adapte le type variables du modèle repris à celles du modèle à calculer
	  IF(en_massep .NEQV. en_masse)THEN
	   IF(en_massep)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1225) ; WRITE(2,1225)
	    CASE DEFAULT
	     WRITE(*,225) ; WRITE(2,225)
	    END SELECT
	    bp(3,:)=SQRT(ABS(bp(3,:)))      !R**2   ==> R
	    bp(4,:)=SQRT(ABS(bp(4,:)))**3   !L**2/3 ==> L
	    bp(5,:)=SQRT(ABS(bp(5,:)))**3   !M**2/3 ==> M
	   ELSE
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1226) ; WRITE(2,1226)
	    CASE DEFAULT
	     WRITE(*,226) ; WRITE(2,226)
	    END SELECT
	    bp(3:5,:)=bp(3:5,:)**2          !R, L, M ==> (R, L, M)**2
	    bp(4,:)=bp(4,:)**(1.d0/3.d0)    !L**2 ==> L**2/3
	    bp(5,:)=bp(5,:)**(1.d0/3.d0)    !M**2 ==> M**2/3
	   ENDIF
	  ENDIF

c	  DO i=1,n
c	   CALL bsp1dn(ne,bp,q,qt,n,ord_qsp,knot,.TRUE.,q(i),lq,fqs,dfqs)
c	   WRITE(*,2000)fqs
c	  ENDDO
c	  PAUSE'reprise'

c modèle d'initialisation de PMS en ASCII
	 ELSE

	  OPEN(unit=4,form='formatted',status='old',file=nom_fich1)

c détermination du nombre de couches
	  n=0
	  DO
	   READ(4,*,END=160)fqs(1) ; n=n+1
	  ENDDO
160	  REWIND(unit=4)
c	  PRINT*,n ; PAUSE'nb.couches'

c les masses sont en 1-M/Mtot
	  ALLOCATE(bp(ne,n)) ; bp=0.d0
	  DO i=1,n
	   READ(4,*)bp(5,i),bp(1:4,i)
	   bp(5,i)=(1.d0-bp(5,i))*mtot
c	   WRITE(*,2000)bp(1:5,i)
	  ENDDO
	  CLOSE(unit=4)
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1016)nom_fich1,n ; WRITE(2,1016)nom_fich1,n
	  CASE DEFAULT
	   WRITE(*,16)nom_fich1,n ; WRITE(2,16)nom_fich1,n
	  END SELECT

c retournement 1 au centre, n à la surface
	  bp(:,1:n:+1)=bp(:,n:1:-1)

c on passe en variables d'intégration
c m--> m**2/3-1, p--> ln p, t--> ln t, r--> r**2, l-->l**2/3
c ou bien m--> m, p--> ln p, t--> ln t, r--> r, l-->l
	  bp(1:2,1:n)=LOG(bp(1:2,1:n))            !pour Ptot et T
	  IF(pturb)bp(Ipg,:)=bp(1,:)          !Pgas avec Pturb
	  IF(en_masse)THEN
	   bp(3,1:n)=bp(3,1:n)**2                 !pour R
	   bp(4:5,1:n)=bp(4:5,1:n)**(2.d0/3.d0)   !pour  L et M
	  ENDIF
	  bp(6,:)=1.d0	       !pour psi (fictif)
	  ord_qsp=2	 !initialisation sur snoein
	  ALLOCATE(q(n),qt(n+ord_qsp)) ; q=(/ (i, i=1,n) /)
	  CALL bsp1dn(ne,bp,q,qt,n,ord_qsp,knot,.FALSE.,q(1),lq,fqs,dfqs)
	  IF(no_croiss)PRINT*,'Pb. at 30 in cesam'
	 ENDIF      !initialisation binaire/ASCII

c on fixe le nb de couches du modèle d'initialisation à ncouche=300
c 300 étant un bon ordre de grandeur pour un modèle solaire
c ce nombre sera ensuite ajusté de facon à ce que dQ/dq = psi0
c on détermine les ncouche points de masse m(ncouche) pour
c assurer une répartition approximativement uniforme en
c ctep lnP + ctet lnT  + cter zeta + ctel lambda + ctem mu

c recherche de la fonction d'espacement
	 ALLOCATE(esp(n))
	 DO i=1,n
	  CALL bsp1dn(ne,bp,q,qt,n,ord_qsp,knot,.TRUE.,q(i),lq,fqs,dfqs)
	  IF(no_croiss)PRINT*,'Pb. at 11 in cesam'
	  esp(i)=ctep*fqs(1)+ctet*fqs(2)+cter*fqs(3)+ctel*fqs(4)+ctem*fqs(5)  
c	  WRITE(*,2000)q(i),fqs,esp(i)
	  IF(i > 1)THEN
	   DO WHILE(esp(i) < esp(i-1))
	    esp(i)=esp(i)+0.01d0
	   ENDDO
	  ENDIF
	 ENDDO

c	 PRINT*,'ord_qsp,knot,n',ord_qsp,knot,n
c	 WRITE(*,2000)ctel,ctep,ctet,ctem,cter ; WRITE(*,2000)esp
c	 WRITE(*,2000)q ; WRITE(*,2000)qt ; PAUSE'esp'

c on dispose les ncouches pour assurer une répartition approx.
c uniforme de la fonction de répartition sur n_qs couches
	 n_qs=n_init ; ALLOCATE(qb(n_qs))
	  
	 CALL zoning(n,q,esp,n_qs,qb) !choix des nouveaux q dans qb

c	 PRINT*,'anciens q',n ; WRITE(*,2000)q(1:n) ; PRINT*,'esp'
c	 WRITE(*,2000)esp(1:n) ; PRINT*,'new q' ; WRITE(*,2000)qb(1:n_qs)
c	 PAUSE'après zoning'

c esp est désormais inutile
	 DEALLOCATE(esp)

c dans bb on se place aux n_qs points qb
c on calcule la nouvelle fonction de répartition
	 ALLOCATE(bb(ne,n_qs),new_esp(1,n_qs))
	 DO i=1,n_qs
	  CALL bsp1dn(ne,bp,q,qt,n,ord_qsp,knot,.TRUE.,qb(i),lq,fqs,dfqs)
	  IF(no_croiss)PRINT*,'Pb. at 12 in cesam'
	  bb(1:ne,i)=fqs(1:ne)    !; WRITE(*,2000)bb(1:ne,i)
	  new_esp(1,i)=ctep*fqs(1)+ctet*fqs(2)+cter*fqs(3)+ctel*fqs(4)
	1 +ctem*fqs(5)
c 	  WRITE(*,2000)q(i),fqs,new_esp(1,i)
	  IF(i > 1)THEN
	   DO WHILE(new_esp(1,i) < new_esp(1,i-1))
	    new_esp(1,i)=new_esp(1,i)+0.01d0
	   ENDDO
	  ENDIF
	 ENDDO
c	 PAUSE'new_esp'

c spline d'ordre 2 sur les nouveaux qb sur n_qs points pour new_esp
	 ALLOCATE(qbt(n_qs+ord_qsp)) ; qb=(/ (i, i=1,n_qs) /)
	 CALL bsp1dn(1,new_esp,qb,qbt,n_qs,2,knotb,.FALSE.,qb(1),lq,fqs,dfqs)
	 IF(no_croiss)PRINT*,'Pb. at 31 in cesam'

c dérivée de la fonction de répartition dans bb(6,:)
	 DO i=1,n_qs
	  CALL bsp1dn(1,new_esp,qb,qbt,n_qs,2,knotb,.TRUE.,qb(i),lq,fqs,dfqs)
	  IF(no_croiss)PRINT*,'Pb. at 13 in cesam'
	  bb(6,i)=dfqs(1)
	 ENDDO

c on spline bb
	 CALL bsp1dn(ne,bb,qb,qbt,n_qs,2,knotb,.FALSE.,qb(1),lq,fqs,dfqs)
	 IF(no_croiss)PRINT*,'Pb. at 32 in cesam'

c	 DO i=1,n_qs
c	 CALL bsp1dn(ne,bb,qb,qbt,n_qs,2,knotb,.TRUE.,qb(i),lq,fqs,dfqs)
c	  WRITE(*,2000)qb(i),fqs
c	 ENDDO
c	 PAUSE'bb'

c dimension de la base de splines pour l'intégration equi. quasi. stat.
	 dim_qs=(n_qs-1)*m_qs+r_qs ; ord_qs=m_qs+r_qs

c vecteur nodal pour l'intégration sur n_qs points
	 DEALLOCATE(q,qt)
	 ALLOCATE(q(n_qs),qt(dim_qs+ord_qs)) ; q=(/ (i, i=1,n_qs) /)
	 CALL noedif(q,n_qs,m_qs,r_qs,qt,knot)
	 IF(no_croiss)THEN
	  PRINT*,'ARRET 15 dans cesam / STOP 15 in cesam' ; STOP
	 ENDIF

c	 PRINT*,ord_qs,knot,dim_qs ; WRITE(*,2000)qt(1:knot)
c	 PAUSE'noedif'

c on place la solution initiale dans la base de q, bb ==> bp
	 DEALLOCATE(bp) ; ALLOCATE(bp(ne,dim_qs))
	 CALL newspl(ne,qb,qbt,knotb,2,qt,knot,ord_qs,bb,bp)
c	 PAUSE'newspl'

c bb, qb, qbt sont désormais inutiles
	 DEALLOCATE(bb,qb,qbt)

c au centre R=L=M=0
	 bp(3,1)=0.d0 ; bp(4,1)=0.d0 ; bp(5,1)=0.d0

c	 vérification
c	 DO i=1,n_qs
c	  CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),lq,fqs,dfqs)
c	  WRITE(*,2000)q(i),fqs
c	 ENDDO
c	 PAUSE'vérification'

c initialisation du vecteur de rotation
c pour la couche externe : force centifuge < force de gravité / 2.
	 IF(wrot == 0.d0)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1083) ; WRITE(2,1083)
	  CASE DEFAULT
	   WRITE(*,83) ; WRITE(2,83)
	  END SELECT
	 ELSE
	  CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(n_qs),lq,fqs,dfqs)
	  IF(no_croiss)PRINT*,'Pb. at 14 in cesam'
	  IF(en_masse)THEN
	   fqs(5)=SQRT(fqs(5))**3 ; fqs(3)=SQRT(fqs(3))
	  ENDIF
	  fqs(5)=fqs(5)*msol ; fqs(3)=fqs(3)*rsol
	  w_max=SQRT(g*fqs(5)/fqs(3)**3)/2.d0

c si la rotation est entrée en kms/s w = v / R
	  IF(unit == 'kms/s')THEN
	   wrot=wrot*1.d5/fqs(3) ; w_rot=wrot
	  ENDIF
	  IF(wrot > w_max)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1400)wrot,w_max,TRIM(chain_don)
	   CASE DEFAULT
	    WRITE(*,400)wrot,w_max,TRIM(chain_don)
	   END SELECT
	   STOP
	  ENDIF

c rot_min : valeur négligeable
	  IF(Krot == 3 .OR. Krot == 4)rot_min(1)=wrot*1.d-2
	 ENDIF

c appel d'initialisation pour tabulation des réactions nucléaires
c allocations fictives
	 ALLOCATE(comp(0),dcomp(0),jac(0,0),ex(0))
	 CALL nuc(1.5d+07,1.5d+02,comp,dcomp,jac,.FALSE.,0,
	1 epsilo,et,ero,ex,hh,be7,b8,n13,o15,f17)

c détermination des abondances initiales
	 DEALLOCATE(comp) ; ALLOCATE(comp(nchim))
	 CALL nuc(1.5d+07,1.5d+02,comp,dcomp,jac,.FALSE.,1,
	1 epsilo,et,ero,ex,hh,be7,b8,n13,o15,f17)

c la rotation
	 WRITE(*,17)thw(Krot) ; WRITE(2,17)thw(Krot)
	 SELECT CASE(Krot)
	 CASE(0,1,2)
	  nrot=0 ; n_rot=0 ; knotr=0
	  ALLOCATE(rota(nrot,n_rot),mrot(n_rot),mrott(knotr))
	 END SELECT
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1023)wrot ; WRITE(2,1023)wrot
	 CASE DEFAULT
	  WRITE(*,23)wrot ; WRITE(2,23)wrot
	 END SELECT

c pour avoir un algorithme unique pour la diffusion, la composition
c chimique est toujours tabulee en fonction de mu=(m/Msol)^2/3
c que ce soit en lagragien ou en eulérien
	 n_ch=n_qs
	 ALLOCATE(chim(nchim,n_ch),mc(n_ch))
	 DO i=1,n_qs
	  CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),lq,fqs,dfqs)
	  IF(no_croiss)PRINT*,'Pb. at 15 in cesam'
c	  WRITE(*,2000)f(1:ne)
	  fqs(5)=MAX(fqs(5),0.d0)
	  IF(en_masse)THEN
	   mc(i)=fqs(5)                 !interpolation en m**2/3
	  ELSE
	   mc(i)=fqs(5)**(2.d0/3.d0)
	  ENDIF
	  chim(1:nchim,i)=comp(:)
	 ENDDO

c initialisation de la spline d'inter. de la comp. chim.
	 ALLOCATE(dxchim(nchim),mct(n_ch+m_ch),xchim(nchim))
	 CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,knotc,.FALSE.,mc(1),lq,
	1 xchim,dxchim)
	 IF(no_croiss)PRINT*,'Pb. at 34 in cesam'
	 dim_ch=knotc-m_ch ; DEALLOCATE(dxchim,xchim)

c initialisations et allocations diverses
	 mstar=mtot

c tabulation de r2 et m23 qui accélèrent la recherche dans le ssp.
c inter valable pour variables lagrangiennes alors m23=m^2/3, r2=r^2
c ou eulériennes alors m23=m, r2=r
	 ALLOCATE(r2(n_qs),m23(n_qs))
	 DO i=1,n_qs
c(test)   CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),lq,fqs,dfqs)
	  r2(i)=bp(3,m_qs*(i-1)+1) ; m23(i)=bp(5,m_qs*(i-1)+1)
	 ENDDO

c initialisation de rota, modèle de PMS
	 IF(Krot >= 3)CALL initialise_rota

c si on utilise une grille fixe pour la composition chimique
	 IF(grille_fixe)THEN
	  ALLOCATE(mc_fixe(n_qs)) ; nc_fixe=n_qs ; mc_fixe=m23
	  IF(.NOT.en_masse)mc_fixe=mc_fixe**(2.d0/3.d0)
         ENDIF

c estimation provisoire de rstar
	 rstar=r2(n_qs) ; IF(en_masse)rstar=SQRT(rstar)

c initialisation de l'âge et du pas temporel pour le modèle initial
c  et calcul du premier modèle de PMS homogène
c le premier pas temporel en PMS étant petit on n'en modifie pas la
c si on diffuse le moment cinétique (Krot=3, 4,et 5)
	 age=0.d0 ; dt=0.d0 ; dts=dt

c on remplace la routine de réac. thermonucléaire lue dans le
c fichier de données mon_modèle.don par iben
	 nom_save=nom_nuc ; nom_nuc='iben'
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1800)
1800	  FORMAT('-----------------------------------',/,
	1 'calculation of the initial PMS model')
	 CASE DEFAULT
	  WRITE(*,800)
800	  FORMAT('-----------------------------------',/,
	1 'calcul du modèle initial de PMS')
	 END SELECT
	 B3: DO           !1-iere boucle infinie pour le choix de c_iben
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1801)
1801	   FORMAT('enter the contraction constant c_iben, examples : ',/,
	1  'for Tc=100 000 K, entrer c_iben= 2.d-2',/,
	2  'for Tc=500 000 K, enter c_iben= 5.d-4',/,
	3  'for Tc=1 000 000 K, enter c_iben= 8.d-5',/,
	4  'for M*>10 Msun, enter c_iben= 4.d-2')
	  CASE DEFAULT
	   WRITE(*,801)
801	   FORMAT('entrer la constante de contraction c_iben, exemples : ',/
	1  'pour Tc=100 000 K, entrer c_iben= 2.d-2',/,
	2  'pour Tc=500 000 K, entrer c_iben= 5.d-4',/,
	3  'pour Tc=1 000 000 K, entrer c_iben= 8.d-5'/,
	4  'pour M*>10 Msol, entrer c_iben= 4.d-2')
	  END SELECT
	  READ*,c_iben

c résolution du Pb. quasi-static avec c_iben
c	  PAUSE'avant resout'
	  CALL resout(3,dt,dts)

	  B2: DO    !2-sde boucle infinie pour le choix de c_iben
	   rp1=bp(3,m_qs*(n_qs-1)+1) ; lp1=bp(4,m_qs*(n_qs-1)+1)
	   IF(en_masse)THEN
	    rp1=SQRT(rp1) ; lp1=SQRT(lp1**3)
	   ENDIF
	   tp1=EXP(bp(2,1))
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1055)c_iben,rp1,lp1,tp1
1055	    FORMAT('Pre main sequence first model : C=',es10.3,/,
	1   'Rext=',es10.3,', Lext=',es10.3,', Tc=',es10.3,/,
	2   'ok? enter o for yes , n for no, to stop enter s')
	   CASE DEFAULT
	    WRITE(*,55)c_iben,rp1,lp1,tp1
55	    FORMAT('Pré Séquence Principale premier modèle : C=',es10.3,/,
	1   'Rext=',es10.3,', Lext=',es10.3,', Tc=',es10.3,/,
	2   'ok? entrer o/n, pour arrêter entrer s')
	   END SELECT
	   READ*,oui
c	   PRINT*,'oui ',oui
	   IF(oui == 'o')THEN
	    c_iben=c_iben*1.1d0
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1802)
1802	     FORMAT('calculation of a new model with 1.1 C')
	    CASE DEFAULT
	     WRITE(*,802)
802	     FORMAT('calcul d''un nouveau modèle avec 1.1 C')
	    END SELECT
	    CALL resout(-3,dt,dts)
	    EXIT B2
	   ELSEIF(oui == 's')THEN
	    STOP
	   ELSE
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1056)c_iben
1056	     FORMAT('c_iben requires a new value, previous value : ',es10.3)
	    CASE DEFAULT
	     WRITE(*,56)c_iben
56	     FORMAT('il faut nouvelle valeur pour c_iben,'
	1    ' valeur précédente : ',es10.3)
	    END SELECT
	    CYCLE B3
	   ENDIF
	  ENDDO B2

c le second modèle est-il satisfaisant ?
	  rp2=bp(3,m_qs*(n_qs-1)+1) ; lp2=bp(4,m_qs*(n_qs-1)+1)
	  IF(en_masse)THEN
	   rp2=SQRT(rp2) ; lp2=SQRT(lp2**3)
	  ENDIF
	  tp2=EXP(bp(2,1))

c estimation du dt correspondant
	  dt=cte2*ABS(rp1-rp2)/rp1/rp2/(lp1+lp2)
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1057)c_iben,rp2,lp2,tp2,dt
1057	   FORMAT('Pre main sequence second model, : C=',es10.3,/,
	1  'Rext=',es10.3,', Lext=',es10.3,', Tc=',es10.3,
	2  ', dt=',es10.3,/,'ok? enter o for yes, n for no,',
	3  ' to stop enter s')
	  CASE DEFAULT
	   WRITE(*,57)c_iben,rp2,lp2,tp2,dt
57	   FORMAT('Pré Séquence Principale second modèle : C=',es10.3,/,
	1  'Rext=',es10.3,', Lext=',es10.3,', Tc=',es10.3,
	2  ', dt=',es10.3,/,' ok? entrer o/n, pour arrêter entrer s')
	  END SELECT
	  READ*,oui
	  IF(oui == 'n')THEN          !nouvelle cte de contraction
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1803)
1803	    FORMAT('a new value is required for c_iben')
	   CASE DEFAULT
	    WRITE(*,803)
803	    FORMAT('il faut une nouvelle valeur pour c_iben')
	   END SELECT
	   dt=0.d0 ; CYCLE B3
	  ELSEIF(oui == 's')THEN
	   STOP
	  ELSE                        !on est content

c on remplace la routine de réactions thermonucléaires fictive iben
c par celle du fichier de données mon_modèle.don
	   nom_nuc=nom_save

c on garde le premier modèle de PMS d'âge 0 dans fichier xxxx_B.pms
	   model_num=0
	   IF(all_rep)THEN
	    WRITE(number,70)model_num
	    modelb=TRIM(nom_fich2)//number//'_B.pms'
	   ELSE
	    modelb=TRIM(nom_fich2)//'_B.pms'
	   ENDIF
c	   WRITE(*,1)ADJUSTL(modelb) ; WRITE(2,1)ADJUSTL(modelb)
	   OPEN(unit=26,form='unformatted',status='unknown',file=modelb)
	   WRITE(26)ne,m_qs,n_qs,knot,nchim,n_ch,m_ch,knotc,Krot,nrot,n_rot,
	1  m_rot,knotr,n_tds,knot_tds,mtot,alpha,w_rot,lim_ro,diffusion,
	2  rot_solid,precision,en_masse,f_eos,f_opa,
	3  nom_ctes,nom_pertm,nom_pertw,nom_tdetau,nom_atm,nom_conv,nom_nuc,
	4  nom_nuc_cpl,nom_diffm,nom_difft,nom_diffw,nom_etat,nom_opa,
	5  nom_elem,bp,q,qt,chim,mc,mct,rota,mrot,mrott,tds,x_tds,xt_tds,
	6  m23,r2,m_zc,r_zc,r_ov,age,dt,dts,mstar,rstar,mw_tot,wrot,
	7  jlim,lconv,lim,model_num
	   CLOSE(unit=26)
	   model_num=-1	!est augmenté d'une unité ligne ~4155

c dtp=0.d0 : on poursuit le modèle PMS obtenu
	   dtp=0.d0 ; dts=dt
	   EXIT B3
	  ENDIF
	 ENDDO B3

c pour l'entête du modèle du fichier mon_modèle.lis
	 list_pms=.TRUE.

c initialisation des paramètres pour le vent et les chutes de planétoïdes
	 CALL vent ; CALL planetoides

	CASE DEFAULT
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1804)un23
1804	  FORMAT('ARRET un23=',i3,/,
	1 'is different from -3, -2, -1, 1, 2, 3')
	 CASE DEFAULT
	  WRITE(*,804)un23
804	  FORMAT('ARRET, un23=',i3,/,
	1 'est différent de -3, -2, -1, 1, 2, 3')
	 END SELECT
	 STOP
	END SELECT        !initialisation des modèles

c	PAUSE'après initialisations'

c les modèles de départ, soit de PMS, soit de ZAMS, soit de reprise
c sont établis ou repris, on démarre l'évolution

c----------EVOLUTION--------------

c ouverture du fichier du diagramme HR
c pour les modèles initiaux, PMS ou ZAMS ABS(un23) = 2 ou 3 on ignore
c ce qui existe dans le fichier, pour une reprise ABS(un23) = 1 on
c écrit à la suite de ce qui existe en omettant celles du modèle repris
	SELECT CASE(ABS(un23))
	CASE(1)           !on écrira à la suite
 	 OPEN(unit=53,form='formatted',status='unknown',access='append',
	1 file=TRIM(nom_fich2)//'.HR')  !fichier du diagramme HR
	 ecrHR=.FALSE.
	CASE(2, 3)
	 OPEN(unit=53,form='formatted',status='unknown',
	1 file=TRIM(nom_fich2)//'.HR')  !fichier du diagramme HR
	 ecrHR=.TRUE.
	END SELECT

c on ignore désormais les entrées binaire/ASCII, ZAMS, PMS, ou reprise
	un23=1

c allocations indépendantes du nombre de couches du modèle quasi-statique
	ALLOCATE(alfa_atm(n_atm),beta_atm(n_atm),cp_atm(n_atm),
	1 delta_atm(n_atm),depsx(nchim),dxchim(nchim),dxchimg(nchim),
	2 fd(ne,0:ord_qs-1),hp_atm(n_atm),gamma_atm(n_atm),
	3 grad_atm(n_atm),grada_atm(n_atm),gradc_atm(n_atm),
	4 gradr_atm(n_atm),grad_mj_a(n_atm),ioni(0:NINT(MAXVAL(zi)),nchim),
	5 k_atm(n_atm),ldcapdt_a(n_atm),ldcapdr_a(n_atm),
	6 lnpt_a(n_atm),lnpt_at(n_atm+m_ch),lnta(1,n_atm),
	7 lnt_a(n_atm),mu_atm(n_atm),mue_atm(n_atm),ro_atm(n_atm),
	8 vais_atm(n_atm),xchim(nchim),xchimg(nchim),xchim1(nchim),
	9 xchim1g(nchim),zbar(nchim))

c-------évolution temporelle : boucle infinie Bev------------

	Bev: DO

c allocations générales pour les écritures, il faut changer
c les dimensions lorsque le nombre de couches a varié
	 IF(ALLOCATED(p))THEN
	  IF(SIZE(p) /= n_qs)THEN
	   DEALLOCATE(alfa,anube7,anub8,anun13,
	1  anuo15,anupep,anupp,beta,compchim,compg,convec,cp,
	2  dcapdr,dcapdt,delta,depsdr,depsdt,dcompchim,dcompg,dm23dlnr,dlnpdlnr,
	3  dlntdlnr,dlpp,ener,epsilon,
	4  gamma,grad,gradad,gradconv,gradrad,grad_mj,grad_mu,hp,kap,
	5  l,degene,m,mu,mue,p,pt,r,ro,t,tx_ioni,u,vaissala,w,z)
	   ALLOCATE(alfa(n_qs),anube7(n_qs),anub8(n_qs),anun13(n_qs),
	1  anuo15(n_qs),anupep(n_qs),anupp(n_qs),beta(n_qs),
	2  compchim(nchim,n_qs),compg(nchim,n_qs),convec(n_qs),cp(n_qs),
	3  dcapdr(n_qs),dcapdt(n_qs),delta(n_qs),
	4  depsdr(n_qs),depsdt(n_qs),dcompchim(nchim,n_qs),dcompg(nchim,n_qs),
	5  dm23dlnr(n_qs),dlnpdlnr(n_qs),dlntdlnr(n_qs),
	6  dlpp(n_qs),ener(4,n_qs),epsilon(5,n_qs),
	7  gamma(n_qs),grad(n_qs),gradad(n_qs),gradconv(n_qs),
	8  gradrad(n_qs),grad_mj(n_qs),grad_mu(n_qs),hp(n_qs),kap(n_qs),
	9  l(n_qs),degene(n_qs),m(n_qs),mu(n_qs),mue(n_qs),p(n_qs),pt(n_qs),
	1  r(n_qs),ro(n_qs),t(n_qs),tx_ioni(nchim,n_qs),
	2  u(n_qs),vaissala(n_qs),w(n_qs),z(n_qs))
	  ENDIF
	 ELSE
	  ALLOCATE(alfa(n_qs),anube7(n_qs),anub8(n_qs),anun13(n_qs),
	1  anuo15(n_qs),anupep(n_qs),anupp(n_qs),beta(n_qs),
	2  compchim(nchim,n_qs),compg(nchim,n_qs),convec(n_qs),cp(n_qs),
	3  dcapdr(n_qs),dcapdt(n_qs),delta(n_qs),
	4  depsdr(n_qs),depsdt(n_qs),dcompchim(nchim,n_qs),dcompg(nchim,n_qs),
	5  dm23dlnr(n_qs),dlnpdlnr(n_qs),dlntdlnr(n_qs),
	6  dlpp(n_qs),ener(4,n_qs),epsilon(5,n_qs),
	7  gamma(n_qs),grad(n_qs),gradad(n_qs),gradconv(n_qs),
	8  gradrad(n_qs),grad_mj(n_qs),grad_mu(n_qs),hp(n_qs),kap(n_qs),
	9  l(n_qs),degene(n_qs),m(n_qs),mu(n_qs),mue(n_qs),p(n_qs),pt(n_qs),
	1  r(n_qs),ro(n_qs),t(n_qs),tx_ioni(nchim,n_qs),
	2  u(n_qs),vaissala(n_qs),w(n_qs),z(n_qs))
	 ENDIF

c pour le test \tilde g / g avec diffusion du moment cinétique
	 IF(COUNT(Krot == Ktest) >= 1)THEN
	  IF(ALLOCATED(gtsg))DEALLOCATE(gtsg)
	  ALLOCATE(gtsg(n_qs)) ; gtsg=0.d0
	 ENDIF

c extraction des valeurs physiques des variables
c calcul des énergies PP, CNO, 3alpha, grav pour intégration des énergies
	 IF(ALLOCATED(ex))DEALLOCATE(ex) ; ALLOCATE(mt(n_qs+2),ex(nchim))
	 DO i=1,n_qs                  !extraction p,t,r,l,m
	  CALL bsp1ddn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),lq,fd)
	  IF(no_croiss)PRINT*,'Pb. at 16 in cesam'
	  pt(i)=EXP(fd(1,0))    !Ptot
	  IF(pturb)THEN         !avec pression turbulente 7 inconnues
	   p(i)=EXP(fd(Ipg,0))    !Pgas
	  ELSE                  !sans pression turbulente 6 inconnues
	   p(i)=pt(i)
	  ENDIF
	  IF(pturb .AND. der)THEN
	   dlpp(i)=fd(Ipg,1)/fd(1,1)    !dln Pgaz/dln Ptot
	  ELSE
	   dlpp(i)=1.d0
	  ENDIF
	  t(i)=EXP(fd(2,0))     !variable ln T
	  IF(en_masse)THEN
	   r(i)=SQRT(ABS(fd(3,0)))    !rayon/Rsol
	   l(i)=SQRT(ABS(fd(4,0)))**3 !l/Lsol
	   m(i)=SQRT(ABS(fd(5,0)))**3 !m/Msol
	   IF(fd(3,0) > 0.d0)THEN
	    grav=cte1*m(i)/fd(3,0)
	    dlnpdlnr(i)=2.d0*fd(3,0)*fd(1,1)/fd(3,1)	!fd(3,0)=r**2
	    dlntdlnr(i)=2.d0*fd(3,0)*fd(2,1)/fd(3,1)
	    dm23dlnr(i)=2.d0*fd(3,0)*fd(5,1)/fd(3,1)
	   ELSE
	    dlnpdlnr(i)=0.d0 ; dlntdlnr(i)=0.d0 ; dm23dlnr(i)=0.d0
	   ENDIF
	  ELSE
	   r(i)=ABS(fd(3,0))    !rayon/Rsol
	   l(i)=fd(4,0)	     !l/Lsol
	   m(i)=ABS(fd(5,0))    !m/Msol
	   fd(5,0)=m(i)**(2.d0/3.d0)
	   IF(fd(3,0) > 0.d0)THEN
	    grav=cte1*m(i)/fd(3,0)**2
	    dlnpdlnr(i)=fd(3,0)*fd(1,1)/fd(3,1)	!fd(3,0)=r
	    dlntdlnr(i)=fd(3,0)*fd(2,1)/fd(3,1)
	    dm23dlnr(i)=fd(3,0)*fd(5,1)/fd(3,1)
	   ELSE
	    dlnpdlnr(i)=0.d0 ; dlntdlnr(i)=0.d0 ; dm23dlnr(i)=0.d0
	   ENDIF
	  ENDIF
	  grad_mj(i)=fd(2,1)/fd(1,1)        !le gradient reel pour MJo
c	  WRITE(*,2000)p(i),t(i),r(i),l(i),m(i),grad_mj(i)
c	  WRITE(*,2000)r(i),r2(i),m(i),m23(i)

c la composition chimique
	  CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1  knotc,.TRUE.,MAX(mc(1),MIN(fd(5,0),mc(n_ch))),lq,xchim,dxchim)
	  IF(no_croiss)PRINT*,'Pb. at 17 in cesam'
	  IF(fd(5,0) > 0.d0)THEN
	   dxchimg=dxchim*2.d0/3.d0/SQRT(fd(5,0))
	  ELSE
	   dxchimg=0.d0
	  ENDIF
	  xchimg=xchim ; CALL chim_gram(xchimg,dxchimg)
	  compchim(:,i)=xchim(:) ; compg(:,i)=xchimg(:)
	  dcompchim(:,i)=dxchim(:) ; dcompg(:,i)=dxchimg(:)

c cas particuliers
	  IF(i == 1)THEN
	   lx_stop=lx_stop .OR. xchimg(1) <= x_stop
	  ELSEIF(i == n_qs)THEN		!pour l'atmosphère
	   xchim1=xchim ; xchim1g=xchimg
	  ENDIF
	  
c éléments lourds
	  IF(ihe4 > 1)THEN
c	   z(i)=1.d0-xchimg(1)-xchimg(ihe4)-xchimg(ihe4-1)
	   z(i)=SUM(xchimg(ihe4+1:nchim))
	  ELSE
	   z(i)=z0
	  ENDIF

c équation d'état
	  CALL etat(p(i),t(i),xchimg,.FALSE.,
	1 ro(i),drop,drot,drox,u(i),dup,dut,dux,
	2 delta(i),deltap,deltat,deltax,cp(i),dcpp,dcpt,dcpx,
	3 gradad(i),dgradadp,dgradadt,dgradadx,alfa(i),beta(i),gamma(i))

c vitesse angulaire
	  SELECT CASE(Krot)
	  CASE(0,1,2)
	   w(i)=wrot
	  CASE(3,4)
	   CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,
	1  knotr,.TRUE.,MAX(mrot(1),MIN(fd(5,0),mrot(n_rot))),lq,frot,dfrot)
	   IF(no_croiss)PRINT*,'Pb. at 18 in cesam'
	   w(i)=frot(1)

c MAX \tilde g / g
	   IF(i == 1 .AND. w(1) > 0.d0)THEN
	    gtsg(1)=cte6*w(1)**2/ro(1)
	   ELSE
	    gtsg(i)=cte4*r(i)**3*w(i)**2/m(i)
	   ENDIF

	  CASE(5)
	   CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,
	1  knotr,.TRUE.,MAX(mrot(1),MIN(fd(5,0),mrot(n_rot))),lq,frot,dfrot)
	   IF(no_croiss)PRINT*,'Pb. at 181 in cesam'
	   w(i)=frot(1)
	  END SELECT

c l'énergie graviphique tds=TdS/dt variable globale du module
c mod_static, a été calculée dans resout à l'issue
c du calcul de ZAMS ou PMS ou est celle du modèle à poursuivre
	  CALL nuc(t(i),ro(i),xchim,dxchim,jac,.FALSE.,3,
	1 epsilo,et,ero,ex,hh,be7,b8,n13,o15,f17)
	  CALL bsp1dn(1,tds,x_tds,xt_tds,n_tds,m_tds,knot_tds,.TRUE.,
	1 MIN(x_tds(i),x_tds(n_tds)),lq,ftds,dftds)
	  IF(no_croiss)PRINT*,'Pb. at 19 in cesam'

c après l'appel à nuc on a epsilo(PP+CNO+(3a+C+O), PP, CNO, 3a+C+O, ?)
c correction d'Yveline pour le signe du TdS, ABS(TdS) --> -TdS
	  epsilo(5)=ftds(1)/secon6
	  epsilon(:,i)=epsilo	!PP+CNO+(3a+C+O), PP, CNO, 3a+C+O, grav
	  epsilon(1,i)=epsilon(1,i)-epsilo(5) !PP+CNO+(3a+C+O)-(+TdS)

c les énergies locales
	  ener(1,i)=epsilo(2)	!PP
	  ener(2,i)=epsilo(3)	!CNO
	  ener(3,i)=epsilo(4)	!3alpha + C + O
	  ener(4,i)=-epsilo(5)	!grav
	 ENDDO

c tabulation et intégration des énergies
	 CALL bsp1dn(4,ener,m,mt,n_qs,2,knote,.FALSE.,m(1),lq,fen,dfen)
	 IF(no_croiss)PRINT*,'Pb. at 20 in cesam'
	 CALL sum_n(4,ener,mt,2,knote,.FALSE.,m(1),m(n_qs),fen)
	 DEALLOCATE(mt)

c les énergies globales en Lsol et proportions en %/ énergie nucléaire totale
	 fen=fen*msol/lsol ; e_nuc=SUM(fen(1:3))
	 e_pp=fen(1) ; e_cno=fen(2) ; e_3al=fen(3) ; e_gr=fen(4)	 
c	 PRINT*,'fen,e_nuc,e_gr' ; WRITE(*,2000)fen,e_nuc,e_gr
	 
c pourcentages respectifs
	 IF(e_nuc > 0.d0)THEN	  
	  i_pp=NINT(100.d0*(e_pp/e_nuc)) ; i_cno=NINT(100.d0*(e_cno/e_nuc))
	  i_3a=NINT(100.d0*(e_3al/e_nuc))

c l'énergie graviphique totale peut être négative
	  i_gr=NINT(MIN(100.d0,100.d0*ABS(e_gr/e_nuc)))	   
	 ELSE
	  i_pp=0 ; i_cno=0 ; i_3a=0 ; i_gr=100
	 ENDIF
c	 PRINT*,'énergies',i_pp,i_cno,i_3a,i_gr
c	 WRITE(*,2000)e_nuc,e_pp,e_cno,e_3al,e_gr,dt ; PAUSE'après i_gr'

c estimation de la nature du modèle
	 zams_pr=zams ; post_pr=post ; cohe_pr=cohe ; coca_pr=coca
	 coox_pr=coox
	 IF(compg(1,1) > x_tams)THEN
	  IF(e_nuc >= ABS(e_gr))THEN	  	  
	   chaine='modèle de la série principale '
	   zams=.TRUE. ; post=.FALSE. ; cohe=.FALSE. ; coca=.FALSE.
	   coox=.FALSE.
	  ELSE
	   SELECT CASE(langue)
	   CASE('english')
	    chaine='pre-main sequence model'
	   CASE DEFAULT
	    chaine='modèle de pré-série principale'
	   END SELECT
	   zams=.FALSE. ; post=.FALSE. ; cohe=.FALSE. ; coca=.FALSE.
	   coox=.FALSE.
	  ENDIF
	 ELSEIF(t(1) < 9.5d7)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   chaine='post-main sequence model'
	  CASE DEFAULT
	   chaine='modèle post série principale'
	  END SELECT
	  zams=.FALSE. ; post=.TRUE. ; cohe=.FALSE. ; coca=.FALSE.
	  coox=.FALSE.
	 ELSEIF(t(1) < 6.d8)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   chaine= 'model with helium burning'
	  CASE DEFAULT
	   chaine='modèle avec combustion hélium'
	  END SELECT
	  zams=.FALSE. ; post=.FALSE. ; cohe=.TRUE. ; coca=.FALSE.
	  coox=.FALSE.
	 ELSEIF(t(1) < 1.d9)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   chaine='model with carbon burning'
	  CASE DEFAULT
	   chaine= 'modèle avec combustion carbone'
	  END SELECT
	  zams=.FALSE. ; post=.FALSE. ; cohe=.FALSE. ; coca=.TRUE.
	  coox=.FALSE.
	 ELSE
	  SELECT CASE(langue)
	  CASE('english')
	   chaine='model with oxygen burning'
	  CASE DEFAULT
	   chaine= 'modèle avec combustion oxygène'
	  END SELECT
	  zams=.FALSE. ; post=.FALSE. ; cohe=.FALSE. ; coca=.FALSE.
	  coox=.TRUE.
	 ENDIF

c l'énergie graviphique totale peut être négative	 
	 i_gr=i_gr*(SIGN(1.d0,e_gr))	 
c	 PRINT*,i_pp,i_cno,i_3a,i_gr,chaine ; PAUSE'chaine'

c avec diffusion du moment cinétique, estimation du maximum de \tilde g / g
	 IF(COUNT(Krot == Ktest) >= 1)THEN
	  gtmax=MAXVAL(gtsg) ; igtmax=MAXLOC(gtsg)
	  rgtmax=r(igtmax(1))/rstar ; mgtmax=m(igtmax(1))/mstar
	 ENDIF

c actualisation du nombre de modèles
c nb_modeles : nombre de modèles depuis le début du calcul ou la reprise
c model_num : nombre de modèles depuis l'initialisation ZAMS ou PMS
	 nb_modeles=nb_modeles+1 ; model_num=model_num+1

c pour les "long run" vers les stades avancés
c après la séquence principale on impose l'approximation de Kippenhahn
c ce qui évite le calcul de dV/dt à travers une discontinuité de
c composition chimique en particulier à la limite d'un coeur convectif
c	 IF(precision == 'av')kipp=post .OR. cohe .OR. coca .OR. coox

c sauf pour le premier passage (alors passe=.FALSE.)
	 IF(passe)THEN
	  zams_e= .NOT.zams_pr .AND. zams   !on passe PMS--->ZAMS
	  post_e= .NOT.post_pr .AND. post   !on passe ZAMS-->POST
	  cohe_e= .NOT.cohe_pr .AND. cohe   !on passe POST-->COHE
	  coca_e= .NOT.coca_pr .AND. coca   !on passe COHE-->COCA
	  coox_e= .NOT.coox_pr .AND. coox   !on passe COCA-->coox
	 ENDIF
	 passe=.TRUE.

c écriture des fichiers mon_modèle.zams, mon_modèle.post,
c mon_modèle.cohe, mon_modèle.coca oumon_modèle.coox suivant le cas
	 IF(zams_e .OR. post_e .OR. cohe_e .OR. coca_e .OR. coox_e)THEN
	  IF(all_rep)THEN
	   WRITE(number,70)model_num
	  ELSE
	   number='    '
	  ENDIF
	  IF(zams_e)THEN
	   modelb=TRIM(TRIM(nom_fich2)//number)//'_B.zams'
	  ELSEIF(post_e)THEN
	   modelb=TRIM(TRIM(nom_fich2)//number)//'_B.post'
c augmentation du nombre de couches et diminution du pas temporel
c	   psi0=psi0*0.9d0 ;
	   dt=MIN(dt,0.1d0)
	  ELSEIF(cohe_e)THEN
	   modelb=TRIM(TRIM(nom_fich2)//number)//'_B.cohe'
c augmentation du nombre de couches (diminution de psi0)
c	   psi0=psi0*0.9d0
	  ELSEIF(coca_e)THEN
	   modelb=TRIM(TRIM(nom_fich2)//number)//'_B.coca'
c augmentation du nombre de couches (diminution de psi0)
c	   psi0=psi0*0.6d0
	  ELSEIF(coox_e)THEN
	   modelb=TRIM(TRIM(nom_fich2)//number)//'_B.coox'
c augmentation du nombre de couches (diminution de psi0)
c	   psi0=psi0*0.6d0
	  ENDIF
	  OPEN(unit=26,form='unformatted',status='unknown',file=TRIM(modelb))
	  WRITE(26)ne,m_qs,n_qs,knot,nchim,n_ch,m_ch,knotc,Krot,nrot,n_rot,
	1 m_rot,knotr,n_tds,knot_tds,mtot,alpha,w_rot,lim_ro,diffusion,
	2 rot_solid,precision,en_masse,f_eos,f_opa,
	3 nom_ctes,nom_pertm,nom_pertw,nom_tdetau,nom_atm,nom_conv,nom_nuc,
	4 nom_nuc_cpl,nom_diffm,nom_difft,nom_diffw,nom_etat,nom_opa,
	5 nom_elem,bp,q,qt,chim,mc,mct,rota,mrot,mrott,tds,x_tds,xt_tds,
	6 m23,r2,m_zc,r_zc,r_ov,age,dt,dts,mstar,rstar,mw_tot,wrot,
	7 jlim,lconv,lim,model_num
	  CLOSE(unit=26)
	 ENDIF

c llist est initialisé .FALSE.
c mettre llist=.TRUE. pour écrire les modèles sur l'écran
c et dans le fichier sortie.dat
c par exemple, pour récuperer un modèle en ASCII
c servant à l'initialisation PMS ou ZAMS
c la différence entre Ptot et Pgaz est ici ignorée

c	 llist=.TRUE.
	 IF(llist)THEN
	  OPEN(unit=4,form='formatted',status='new',file='sortie.dat')
	  WRITE(*,101)age
101	  FORMAT('age=',es10.3)
	  WRITE(*,102)
102	  FORMAT(t9,'1-m',t24,'p',t36,'t',t49,'r',t60,'l',t71,'X')
	  DO k=n_qs,1,-1
	   WRITE(*,2001)1.d0-m(k)/mstar,p(k),t(k),r(k),l(k),compg(1,k)
	   WRITE(4,2001)1.d0-m(k)/mstar,p(k),t(k),r(k),l(k),compg(1,k)
2001	   FORMAT(es17.10,4es12.5,es10.3)
	  ENDDO

c pour la fermeture du dessin TRUE pour appel à pgend
	  CLOSE(unit=4) ; CALL des(.TRUE.,dt,teff)
	  PRINT*,'arrêt après écriture du modèle en ASCII' ; STOP
	 ENDIF

c estimation de la température effective et dessins en ligne
	 IF(n_atm > 1)THEN
	  rext=r(n_qs) ; lext=l(n_qs) ; teff=cte8*SQRT(SQRT(lext)/rext)
	 ELSE
	  CALL atm(.TRUE.,l(n_qs),r(n_qs),xchim1g,ptext,dptdl,dptdr,
	1 text,dtdl,dtdr,mext,dml,dmr,pext,dpl,dpr,teff)
	 ENDIF

	 IF(age > 0.d0)THEN
c	  WRITE(*,2000)teff,r(n_qs),l(n_qs) ; PAUSE'Teff2'
	  CALL des(.FALSE.,dtp,teff)
	 ELSE
	  CALL des(.FALSE.,0.d0,teff)
	 ENDIF
c	 WRITE(*,2000)teff; PAUSE'des1' ; CALL des(.FALSE.,dt,teff)
c	 PAUSE'des2'

c---------------------------------------------------------------------

c           E C R I T U R E S

c---------------------------------------------------------------------

	 LOGg=LOG10(cte9*mstar) !pour écriture de Log g
	 lteff=LOG10(teff)    !teff approche pour le test d'arrêt

	 sort=(age >= agemax-1.d-5)
	1 .OR. (nb_modeles >= nb_max_modeles)
	2 .OR. lx_stop.OR. lt_stop .OR. lhe_stop
	3 .OR. (zams .AND. arret == 'zams')
	4 .OR. (post .AND. arret == 'post')
	5 .OR. (cohe .AND. arret == 'cohe')
	6 .OR. (coca .AND. arret == 'coca')
	7 .OR. (coox .AND. arret == 'coox')
	8 .OR. list_sort

c formation du fichier de rotation
	 IF(nb_modeles > 0)THEN
	  SELECT CASE(Kdes_rot)
	  CASE(1)
	   IF(sort)CALL ecrit_rota(dt)
	  CASE(2,4)
	   CALL ecrit_rota(dt)
	  END SELECT
	 ENDIF

c au cours de l'évolution (lteffp > 0) on sort si:
c en se déplacant vers les Teff croissantes Log Teff > LOG_teff
c en se déplacant vers les Teff décroissantes Log Teff < -LOG_teff
c LOG_teff étant positif ou négatif selon le cas
	 IF(lteffp > 0.d0 .AND. age > 100.d0)THEN
	  IF(LOG_teff > 0.d0 .AND. lteff > lteffp)THEN
	   sort=sort .OR. lteff > LOG_teff        !à gauche dans HR
	  ELSEIF(LOG_teff < 0.d0 .AND. lteff < lteffp)THEN
	   sort=sort .OR. lteff < -LOG_teff             !à droite dans HR
	  ENDIF
	 ENDIF
	 lteffp=lteff ; n_ecrit=n_ecrit+1

	 ecrit=sort .OR. zams_e .OR. post_e .OR. cohe_e .OR. coca_e .OR. coox_e
	1 .OR. age-agep >= dtlist .OR. list_rep .OR. list_zams
	2 .OR. list_pms .OR. n_ecrit >= ntot_ecrit
	 IF(ecrit)n_ecrit=0

c si écrit: écriture, ensuite si sort: sortie
c	 PRINT*,sort,ecrit,list_rep,list_sort ; PAUSE'list_sort'
	 IF(ecrit .OR. all_output)THEN
	 
c interpolation de variables thermodynamiques
	  CALL tab_vth
	 
c pour la stucture interne	 
	  DO i=1,n_qs
	   xchim=compchim(:,i) ; dxchim=dcompchim(:,i)
	   xchimg=compg(:,i) ; dxchimg=dcompg(:,i)	   
	   CALL thermo(pt(i),p(i),t(i),m(i),l(i),r(i),dlpp(i),xchim,dxchim,
	1  ro(i),drop,drot,drox,u(i),dup,dut,dux,
	2  grad(i),dgradpt,dgradp,dgradt,dgradx,dgradm,dgradl,dgradr,
	3  dgradlpp,gam,dgampt,dgamp,dgamt,dgamx,dgamm,dgaml,dgamr,
	4  dgamlpp,epsilo,depsp,depst,depsx,kap(i),dkapp,dkapt,dkapx,
	5  delta(i),deltap,deltat,deltax,cp(i),dcpp,dcpt,dcpx,
	6  gradad(i),dgradadp,dgradadt,dgradadx,
	7  hp(i),dhppt,dhpp,dhpt,dhpx,dhpr,dhpm,grad_mu(i),
	8  gradrad(i),alfa(i),beta(i),gamma(i),convec(i))

c dérivées de l'opacité et de l'énergie nucléaire
	   dcapdr(i)=dkapp/drop ; dcapdt(i)=dkapt-dcapdr(i)*drot
	   depsdr(i)=depsp/drop ; depsdt(i)=depst-depsdr(i)*drot

c thermo sort avec radiatif=convec=.TRUE. dans ZR
	   convec(i)=.NOT.convec(i)  !thermo sort avec convec=.TRUE. dans ZR
	   IF(convec(i))THEN
	    gradconv(i)=grad(i)
	   ELSE
	    gradconv(i)=0.d0
	   ENDIF

c estimation des taux d'ionisation et calcul des poids
c moléculaires moyens, si T < 4500K on suppose les éléments
c totalement recombinés (saha ne converge pas)	   
	   IF(t(i) > 4.5d3)THEN
	    IF(mu_saha)THEN	 
             CALL saha(xchim,t(i),ro(i),ioni,zbar,nel,degene(i))	  
	    ELSE
	     zbar=zi	  
	    ENDIF
            mue(i)=1.d0/DOT_PRODUCT(zbar,xchim)
	   ELSE
	    zbar=0.d0 ; mue(i)=1.d30 	 
	   ENDIF	 
	   mu(i)=1.d0/DOT_PRODUCT((1.d0+zbar),xchim)

c les taux d'ionisation	
c arb_rom : transformation nb. arabes --> nb. romains   
	   DO j=1,nchim
	    tx_ioni(j,i)=arb_rom(NINT(zbar(j))+1)
	   ENDDO
c	   PRINT*,tx_ioni(:,i) ; PAUSE'tx_ioni'

c	   PRINT*,i,pturb
c	   WRITE(*,2000)pt(i),p(i),t(i),xchim(1),m(i),l(i),r(i),dlpp(i)
c	   WRITE(*,2000)drop,dkapp,dcapdr(i)
c	   WRITE(*,2000)dcapdt(i),dkapt,dcapdr(i),drot
c	   WRITE(*,2000)depsdr(i),depsp,depsdt(i),depst,depsdr(i)

c calcul d'une des formes de Vaissala  (Kippenhahan & Weigert p. 42)
c 1/gamma1 dlnP/dln r - dln ro/dln r =
	   IF(hp(i) <= 0.d0)THEN            !au centre
	    vaissala(i)=0.d0
	   ELSEIF(new_bv)THEN
c	    vaissala(i)=r(i)*rsol/hp(i)*(alfa(i)*(dlpp(i)-1.d0)+
c	1   delta(i)*(gradad(i)-grad_mj(i)
c	2   +beta(i)/(4.d0-3.d0*beta(i))*grad_mu(i)))

c utilisation de la tabulation vth	
c fvth(1)=lnP, fvth(2)=lnT, fvth(3)=r**2, fvth(4)=l**2/3, fvth(5)=ln ro, 
c fvth(6)=cp, fvth(7)=delta, fvth(8)=gamma1, fvth(9)=ln µ, fvth(10)=ln kap
	    CALL bsp1dn(nvth,vth,mc,mct,n_ch,m_ch,knotc,.TRUE.,
	1   MAX(mc(1),MIN(m(i)**(2.d0/3.d0),mc(n_ch))),lq,fvth,dfvth)
	    IF(no_croiss)PRINT*,'Pb. at 211 in cesam'	
		
	    vaissala(i)=(dfvth(1)/fvth(8)-dfvth(5))*fvth(3)/dfvth(3)*2.d0

c ancienne formulation
	   ELSE	   
	    vaissala(i)=r(i)*rsol/hp(i)*(alfa(i)*(dlpp(i)-1.d0)
	1   +delta(i)*(gradad(i)-grad_mj(i)))-cte3*r(i)**3*drox*dxchimg(1)
	   ENDIF
c	   WRITE(*,*)'m,r,drot,t,p,grad,gradad,drox/vai,dchimg,nuc'
c	   WRITE(*,2000)m(i),r(i),drot,t(i),p(i),grad(i),gradad(i),drox
c	   WRITE(*,2000)vaissala(i),dxchimg(1)

c production de neutrinos
	   CALL nuc(t(i),ro(i),xchim,dcomp,jac,.FALSE.,4,
	1  epsilo,et,ero,ex,anupp(i),anupep(i),anub8(i),anube7(i),
	2  anun13(i),anuo15(i))
	  ENDDO

c atmosphère, température effective, rayon extérieur
c	  WRITE(*,*)'wrot,l(n_qs),r(n_qs),xchim en n'
c	  WRITE(*,2000)wrot,l(n_qs),r(n_qs),chim1g(1:nchim))
c	  PAUSE'lim ext 1'
	  CALL atm(.TRUE.,l(n_qs),r(n_qs),xchim1g,ptext,dptdl,dptdr,
	1 text,dtdl,dtdr,mext,dml,dmr,pext,dpl,dpr,teff)

c l'atmosphère n'est détaillée que si elle est reconstituée ie. n_atm > 1
	  IF(n_atm > 1)THEN
	   dxchim=0.d0
	   DO i=1,n_atm         !thermo pour l'atmosphère
	    grav=cte1*m_atm(i)/r_atm(i)**2
	    CALL tdetau(tau(i),teff,grav,bid,dtsdtau,dtsdteff,dtsdg,
	1   ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)
c	    WRITE(*,2000)tau,f_tau,df_tau,d2f_tau

c les gradients
	    CALL thermo_atm(pt_atm(i),p_atm(i),t_atm(i),xchim1g,m_atm(i),
	1   l(n_qs),r_atm(i),dlpp_atm(i),tau(i),df_tau,d2f_tau,rstar,
	2   ro_atm(i),drop,drot,k_atm(i),dkapp,dkapt,
	3   grada_atm(i),dgradadp,dgradadt,
	4   grad_atm(i),dgradpt,dgradp,dgradt,dgradr,dgradrs,dgradm,
	5   dgradtau,dgradlpp,gam,dgampt,dgamp,dgamt,dgamr,dgamrs,dgamm,
	6   dgamtau,dgamlpp,hp_atm(i),dhppt,dhpp,dhpt,dhpr,dhpm,
	7   delta_atm(i),deltap,deltat,cp_atm(i),gradr_atm(i),
	8   alfa_atm(i),beta_atm(i),gamma_atm(i),rad_atm,.TRUE.)
	    IF(rad_atm)THEN
	     gradc_atm(i)=0.d0
	    ELSE
	     gradc_atm(i)=grad_atm(i)
	    ENDIF

	    lnpt_a(i)=LOG(pt_atm(i)) ; lnt_a(i)=LOG(t_atm(i))
c	    WRITE(*,2000)p_atm(i),t_atm(i),xchim1(1),m_atm(i),l(n_qs),
c	1   r_atm(i)
	    ldcapdr_a(i)=dkapp/drop ; ldcapdt_a(i)=dkapt-ldcapdr_a(i)*drot
c	    WRITE(*,2000)m_atm(i),r_atm(i)

c dans l'atmosphère pas de gradient de comp. chim
c pour vaissala on prend ou bien le vrai gradient grad_mj
c ou celui calculé aux points de grille, X=cte dans l'atmosphère
	    vais_atm(i)=r_atm(i)*rsol/hp_atm(i)*(alfa_atm(i)*
	1   (dlpp_atm(i)-1.d0)+delta_atm(i)*(grada_atm(i)-grad_atm(i)))

c estimation des taux d'ionisation et calcul des poids
c moléculaires moyens, si T < 4500K on suppose les éléments
c totalement recombinés (saha ne converge pas)
	    IF(t_atm(i) < 4.5d3)THEN
	     zbar=0.d0 ; mue_atm(i)=1.d30
	    ELSE
	     CALL saha(xchim1,t_atm(i),ro_atm(i),ioni,zbar,nel,bid)
	     mue_atm(i)=1.d0/DOT_PRODUCT(zbar,xchim1)
	    ENDIF
	    mu_atm(i)=1.d0/DOT_PRODUCT((1.d0+zbar),xchim1)
	   ENDDO

c estimation du gradient dans l'atmosphère les abs. devant être croissantes
c on prend l'opposé pour tabuler ln T en fct. de ln P
	   IF(n_atm > m_ch)THEN
	    lnpt_a=-lnpt_a ; lnta(1,:)=lnt_a
	    CALL bsp1dn(1,lnta,lnpt_a,lnpt_at,n_atm,m_ch,
	1   k,.FALSE.,lnpt_a(1),lq,fatm,dfatm)
	    IF(no_croiss)PRINT*,'Pb. at 212 in cesam'	
	    DO i=1,n_atm
	     CALL bsp1dn(1,lnta,lnpt_a,lnpt_at,n_atm,m_ch,k,.TRUE.,
	1    lnpt_a(i),lq,fatm,dfatm)
	     IF(no_croiss)PRINT*,'Pb. at 21 in cesam'
	     grad_mj_a(i)=fatm(1)
	    ENDDO
	    lnpt_a=-lnpt_a ; grad_mj_a=-grad_mj_a
	   ELSE
	    grad_mj_a=0.4d0
	   ENDIF
	  ENDIF

c calcul de dérivées secondes au centre
	  d2p=cte5*ro(1)**2*rstar**2/p(1)
	  fonc(0)=ro(1) ; fonc(1)=ro(2) ; fonc(2)=ro(3)
	  fonc(3)=ro(4) ; fonc(4)=ro(5) ; fonc(5)=ro(6)
	  absc(0)=r(1) ; absc(1)=r(2) ; absc(2)=r(3)
	  absc(3)=r(4) ; absc(4)=r(5) ; absc(5)=r(6)

c	  PRINT*,'ro/r/deriv'
c	  WRITE(*,2000)(fonc(i),i=0,5) ; WRITE(*,2000)(absc(i),i=0,5)

	  CALL newton(0,m_qs+1,fonc,absc,0.d0,poly,2)
	  d2ro=rstar**2/ro(1)*poly(2)

c	  WRITE(*,2000)(poly(i),i=0,2)
c	  PRINT*,'rstar,d2ro,ro,poly(2)'
c	  WRITE(*,2000)rstar,d2ro,ro(n_qs),poly(2)
c	  PRINT*,'d2p,d2ro,rstar,p(n_qs),ro(n_qs)'
c	  WRITE(*,2000)d2p,d2ro,rstar,p(n_qs),ro(n_qs)

c vérification de la solution: estimation des premiers et 2d membres
c des équations d'évolution
	  IF(.FALSE.)THEN
c	  IF(.TRUE.)THEN
	   DO i=2,n_qs-1
	    WRITE(*,2000)(t(i+1)-t(i-1))/(m(i+1)-m(i-1))/(mtot*msol),
	1   -g*m(i)*mtot*msol/4.d0/pi/(rsol*r(i))**4*grad(i)*t(i)/p(i),
	2   (pt(i+1)-pt(i-1))/(m(i+1)-m(i-1))/(mtot*msol),
	3   -g*m(i)*mtot*msol/4.d0/pi/(rsol*r(i))**4,
	4   (r(i+1)-r(i-1))*rsol/(m(i+1)-m(i-1))/(mtot*msol),
	5   1.d0/4.d0/pi/(rsol*r(i))**2/ro(i),
	6   (l(i+1)-l(i-1))*lsol/(m(i+1)-m(i-1))/(mtot*msol),
	7   epsilon(1,i)
	   ENDDO
	   PAUSE'verif'
	  ENDIF
	  IF(list_sort .OR. list_rep)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(2,1805)model_num
1805	    FORMAT('Model number to be continued :'i4)
	   CASE DEFAULT
	    WRITE(2,805)model_num
805	    FORMAT('Reprise du modèle n° : ',i4)
	   END SELECT
	  ELSEIF(list_zams)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(2,1806)
1806	    FORMAT('Homogeneous ZAMS model :')
	   CASE DEFAULT
	    WRITE(2,806)
806	    FORMAT('Modèle de ZAMS homogène :')
	   END SELECT
	  ELSEIF(list_pms)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(2,1807)
1807	    FORMAT('Initial PMS model :')
	   CASE DEFAULT
	    WRITE(2,807)
807	    FORMAT('Modèle initial de PMS :')
	   END SELECT
	  ENDIF

c sortie
	  IF(sort)THEN
	   CALL list(alfa,anub8,anube7,anun13,anuo15,anupep,anupp,beta,
	1  compg,cp,delta,dcapdr,dcapdt,depsdr,depsdt,d2p,d2ro,
	2  chaine,convec,sort,epsilon,gamma,gamma_atm,gradad,grada_atm,
	3  gradconv,gradc_atm,gradrad,gradr_atm,hp,i_cno,i_gr,i_pp,i_3a,
	4  kap,l,m,mu,mue,m_atm,p,pt,pt_atm,p_atm,r,ro,ro_atm,r_atm,t,tau,
	5  teff,tx_ioni,t_atm,u,vaissala,w,z,degene)

c output ASCII
	   IF(nom_output /= 'no_output')CALL ecrit_ascii

c sorties sur écran~~~~~~~~~~~~~~

c Alerte avec diffusion du moment cinétique si \tilde g \ g > 20%
	   IF(COUNT(Krot == Ktest) >= 1)THEN
	    IF(gtmax > 0.2d0)THEN
	     SELECT CASE(langue)
	     CASE('english')
	      WRITE(*,1808)nom_thw,gtmax,igtmax,rgtmax,mgtmax
	      WRITE(2,1808)nom_thw,gtmax,igtmax,rgtmax,mgtmax
1808	      FORMAT(/,'WARNING The formalism of ',a,/,
	1     'is no longer valuable : MAX acce.centr./gravity=',es10.3,
	2     ', shell ',i4,', R/R*=',es10.3,', M/M*=',es10.3)
	     CASE DEFAULT
	      WRITE(*,808)nom_thw,gtmax,igtmax,rgtmax,mgtmax
	      WRITE(2,808)nom_thw,gtmax,igtmax,rgtmax,mgtmax
808	      FORMAT(/,'ATTENTION le formalisme de ',a,/,
	1     'est injustifié : MAX accé.centr./gravité=',es10.3,', couche ',i4,
	2     ', R/R*=',es10.3,', M/M*=',es10.3)
	     END SELECT
	    ELSE
	     SELECT CASE(langue)
	     CASE('english')
	      WRITE(*,1809)gtmax,igtmax,rgtmax,mgtmax
	      WRITE(2,1809)gtmax,igtmax,rgtmax,mgtmax
1809	      FORMAT(/,'MAX acce.centr./gravity=',es10.3,
	1     ', shell ',i4,', R/R*=',es10.3,', M/M*=',es10.3)
	     CASE DEFAULT
	      WRITE(*,809)gtmax,igtmax,rgtmax,mgtmax
	      WRITE(2,809)gtmax,igtmax,rgtmax,mgtmax
809	      FORMAT(/,'MAX accé.centr./gravité=',es10.3,', couche ',i4,
	1     ', R/R*=',es10.3,', M/M*=',es10.3)
	     END SELECT
	    ENDIF
	   ENDIF

c écriture du cartouche
c	   PRINT*,mstar,mtot
	   pertmm=mstar/mtot-1.d0
	   IF(ihe4 <= 1)THEN
	    WRITE(*,205)age,LOG10(teff),LOG10(l(n_qs)),LOG10(rstar),
	1   LOGg+LOG10(mstar/rstar**2),pt(1),t(1),ro(1),compg(1,1),
	2   i_pp,i_cno,i_3a,i_gr,pertmm,mstar,TRIM(chaine)
205	    FORMAT(/,'*********',/,'âge=',es10.3,', LogTeff=',es10.3,
	1   ', LogL/Lsol=',es10.3,', LogR/Rsol=',es10.3,/,
	2   'Log10 g=',es10.3,', Pc=',es10.3,', Tc=',es10.3,
	3   ', Roc=',es10.3,', Xc=',es10.3,/,
	4   'ePP/eNUC=',i3,'%, eCNO/eNUC=',i3,'%, e3a+C+O/eNUC=',i3,		
	5   '%, eGRAV/eNUC=',i3,'%',/,'Var. rel. de masse : ',es10.3,
	6   ', M*=',es10.3,'Msol, ',a,/,'*********',/)

	   ELSEIF(w(n_qs) > 0.d0)THEN
	    WRITE(*,207)age,LOG10(teff),LOG10(l(n_qs)),LOG10(rstar),
	1   LOGg+LOG10(mstar/rstar**2),pt(1),t(1),ro(1),compg(1,1),
	2   i_pp,i_cno,i_3a,i_gr,compg(ihe4,1),pertmm,mstar,TRIM(chaine),
	3   2.d0*pi/24.d0/3600.d0/w(n_qs),
	4   rstar*w(n_qs)*rsol*1.d-5
207	    FORMAT(/,'*********',/,'âge=',es10.3,', LogTeff=',es10.3,
	1   ', LogL/Lsol=',es10.3,', LogR/Rsol=',es10.3,/,
	2   'Log g=',es10.3,', Pc=',es10.3,', Tc=',es10.3,
	3   ', Roc=',es10.3,', Xc=',es10.3,/,
	4   'ePP/eNUC=',i3,'%, eCNO/eNUC=',i3,'%, e3a+C+O/eNUC=',i3,		
	5   '%, eGRAV/eNUC=',i3,'%',', Yc=',es10.3,/,	
	6   'Var. rel. de masse : ',es10.3,
	7   ', M*=',es10.3,'Msol, ',a,/,'Période de rota.=',es10.3,
	8   ' jours,',' vitesse de rotation à la surface=',es10.3,
	9   ' km/s',/,'*********',/)

	   ELSE
	    WRITE(*,206)age,LOG10(teff),LOG10(l(n_qs)),LOG10(rstar),
	1   LOGg+LOG10(mstar/rstar**2),pt(1),t(1),ro(1),compg(1,1),
	2   i_pp,i_cno,i_3a,i_gr,compg(ihe4,1),pertmm,mstar,TRIM(chaine)
206	    FORMAT(/,'*********',/,'âge=',es10.3,', LogTeff=',es10.3,
	1   ', LogL/Lsol=',es10.3,', LogR/Rsol=',es10.3,/,
	2   'Log g=',es10.3,', Pc=',es10.3,', Tc=',es10.3,
	3   ', Roc=',es10.3,', Xc=',es10.3,/,
	4   'ePP/eNUC=',i3,'%, eCNO/eNUC=',i3,'%, e3a+C+O/eNUC=',i3,		
	5   '%, eGRAV/eNUC=',i3,'%',', Yc=',es10.3,/,	
	6   'Var. rel. de masse : ',es10.3,
	7   ', M*=',es10.3,'Msol, ',a,/,'*********',/)
	   ENDIF

	   titre=TRIM(version)
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1501)titre ; WRITE(2,1501)titre
1501	    FORMAT('End of evolution with  CESAM2k version : ',a)
	   CASE DEFAULT
	    WRITE(*,501)titre ; WRITE(2,501)titre
501	    FORMAT('Fin d''évolution avec CESAM2k version : ',a)
	   END SELECT
	   IF(age >= agemax-1.d-5)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1409) ; WRITE(2,409)
1409	     FORMAT('Stop because agemax has been reached')
	    CASE DEFAULT
	     WRITE(*,409) ; WRITE(2,409)
409	     FORMAT('Sortie car agemax atteint')
	    END SELECT
	   ELSEIF(nb_modeles >= nb_max_modeles)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1410) ; WRITE(2,410)
1410	     FORMAT('Stop : the maximum number of models has been reached')
	    CASE DEFAULT
	     WRITE(*,410) ; WRITE(2,410)
410	     FORMAT('Sortie : le nombre maximum de modèles est atteint')
	    END SELECT
	   ELSEIF(lx_stop)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1401) ; WRITE(2,1401)
1401	     FORMAT('Stop because  H < x_stop  at the center')
	    CASE DEFAULT
	     WRITE(*,401) ; WRITE(2,401)
401	     FORMAT('Sortie car H < x_stop au centre')
	    END SELECT
	   ELSEIF(t(1) >= t_stop)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1402) ; WRITE(2,1402)
1402	     FORMAT('Stop because  T > t_stop at the center')
	    CASE DEFAULT
	     WRITE(*,402) ; WRITE(2,402)
402	     FORMAT('Sortie  car T > t_stop au centre')
	    END SELECT
	   ELSEIF(zams .AND. arret == 'zams')THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1403) ; WRITE(2,1403)
1403	     FORMAT('Stop at the ZAMS stage')
	    CASE DEFAULT
	     WRITE(*,403) ; WRITE(2,403)
403	     FORMAT('Sortie sur la ZAMS')
	    END SELECT
	   ELSEIF(post .AND. arret == 'post')THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1404) ; WRITE(2,1404)
1404	     FORMAT('stop at the end of the ZAMS stage (TAMS)')
	    CASE DEFAULT
	     WRITE(*,404) ; WRITE(2,404)
404	     FORMAT('Sortie à la fin de la ZAMS (TAMS)')
	    END SELECT
	   ELSEIF(cohe .AND. arret == 'cohe')THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1405) ; WRITE(2,1405)
1405	     FORMAT('Stop at the He burning onset')
	    CASE DEFAULT
	     WRITE(*,405) ; WRITE(2,405)
405	     FORMAT('Sortie à la combustion de He')
	    END SELECT
	   ELSEIF(coca .AND. arret == 'coca')THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1406) ; WRITE(2,1406)
1406	     FORMAT('Stop at the C combustion onset')
	    CASE DEFAULT
	     WRITE(*,406) ; WRITE(2,406)
406	     FORMAT('Sortie à la combustion de C')
	    END SELECT
	   ELSEIF(coox .AND. arret == 'coox')THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1412) ; WRITE(2,1412)
1412	     FORMAT('Stop at the O combustion onset')
	    CASE DEFAULT
	     WRITE(*,412) ; WRITE(2,412)
412	     FORMAT('Sortie à la combustion de O')
	    END SELECT
	   ELSEIF(LOG_teff > 0.d0 .AND. lteff > LOG_teff)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1407) ; WRITE(2,1407)
1407	     FORMAT('Stop because, at center, LOG(T) > LOG_teff')
	    CASE DEFAULT
	     WRITE(*,407) ; WRITE(2,407)
407	     FORMAT('Sortie car, au centre, LOG(T) > LOG_teff')
	    END SELECT
	   ELSEIF(LOG_teff < 0.d0 .AND. lteff < -LOG_teff)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1408) ; WRITE(2,1408)
1408	     FORMAT('Stop because, at center, LOG(T) < -LOG_teff')
	    CASE DEFAULT
	     WRITE(*,408) ; WRITE(2,408)
408	     FORMAT('Sortie car, au centre, LOG(T) < -LOG_teff')
	    END SELECT
	   ELSEIF(lhe_stop)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1411)he_core ; WRITE(2,1411)he_core
1411	     FORMAT('Stop because the he_core extends until ',es10.3,'Msol')
	    CASE DEFAULT
	     WRITE(*,411)he_core ; WRITE(2,411)he_core
411	     FORMAT('Sortie car le noyau d''hélium atteint ',es10.3,'Msol')
	    END SELECT
	   ENDIF
	   CALL date_and_time(date,time,zone,values)
	   SELECT CASE(langue)
	   CASE('english')
	    titre='End of the calculation on '//date(7:8)//
	1   ' '//TRIM(month(values(2)))//
	2   ' '//date(1:4)//' à '//time(1:2)//'h'//time(3:4)
	   CASE DEFAULT
	    titre='Fin du calcul le '//date(7:8)//' '//TRIM(mois(values(2)))//
	1   ' '//date(1:4)//' à '//time(1:2)//'h'//time(3:4)
	   END SELECT
	   WRITE(*,502)titre ; WRITE(2,502)titre
502	   FORMAT(a)
	   titre=TRIM(nom_fich2)
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1516)titre ; WRITE(2,1516)titre
1516	    FORMAT('name of the model : ',/,a)
	   CASE DEFAULT
	    WRITE(*,516)titre ; WRITE(2,516)titre
516	    FORMAT('nom du modèle : ',/,a)
	   END SELECT
	   titre=TRIM(nom_fich2)//'.don'
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1503)titre ; WRITE(2,1503)titre
1503	    FORMAT('filename of the data NAMELISTs : ',a)
	   CASE DEFAULT
	    WRITE(*,503)titre ; WRITE(2,503)titre
503	    FORMAT('nom du fichier des NAMELISTs des données : ',a)
	   END SELECT
	   titre=TRIM(nom_fich2)//'.lis'
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1504)titre ; WRITE(2,1504)titre
1504	    FORMAT('filename of the ASCII listing of the results : ',a)
	   CASE DEFAULT
	    WRITE(*,504)titre ; WRITE(2,504)titre
504	    FORMAT('nom du fichier du listing des resultats en ASCII : ',a)
 	   END SELECT
	   titre=TRIM(nom_fich2)//'_B.pms'
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1505)titre ; WRITE(2,1505)titre
1505	    FORMAT('filename of the PMS zero age binary model : ',a)
	   CASE DEFAULT
	    WRITE(*,505)titre ; WRITE(2,505)titre
505	    FORMAT('nom du fichier du modèle d''âge 0 PMS en binaire : ',a)
	   END SELECT
	   titre=TRIM(nom_fich2)//'_B.hom'
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1506)titre ; WRITE(2,1506)titre
1506	    FORMAT('filename of the zero age homogenous binary model : ',a)
	   CASE DEFAULT
	    WRITE(*,506)titre ; WRITE(2,506)titre
506	    FORMAT('nom du fichier du modèle d''âge 0 homogène en binaire : ',a)
	   END SELECT
	   titre=TRIM(nom_fich2)//'_B.zams'
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1507)titre ; WRITE(2,1507)titre
1507	    FORMAT('filename of the binary ZAMS model : ',a)
	   CASE DEFAULT
	    WRITE(*,507)titre ; WRITE(2,507)titre
507	    FORMAT('nom du fichier du modèle de ZAMS en binaire : ',a)
	   END SELECT
	   titre=TRIM(nom_fich2)//'_B.post'
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1508)titre ; WRITE(2,1508)titre
1508	    FORMAT('filename of the POST binary model : ',a)
	   CASE DEFAULT
	    WRITE(*,508)titre ; WRITE(2,508)titre
508	    FORMAT('nom du fichier du modèle de POST en binaire : ',a)
	   END SELECT
	   titre=TRIM(nom_fich2)//'_B.cohe'
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1509)titre ; WRITE(2,1509)titre
1509	    FORMAT('filename of the COHE binary model : ',a)
	   CASE DEFAULT
	    WRITE(*,509)titre ; WRITE(2,509)titre
509	    FORMAT('nom du fichier du modèle de COHE en binaire : ',a)
	   END SELECT
	   titre=TRIM(nom_fich2)//'_B.coca'
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1510)titre ; WRITE(2,1510)titre
1510	    FORMAT('filename of the COCA binary model : ',a)
	   CASE DEFAULT
	    WRITE(*,510)titre ; WRITE(2,510)titre
510	    FORMAT('nom du fichier du modèle de COCA en binaire : ',a)
	   END SELECT
	   titre=TRIM(nom_fich2)//'_B.coox'
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1520)titre ; WRITE(2,1520)titre
1520	    FORMAT('filename of the coox binary model : ',a)
	   CASE DEFAULT
	    WRITE(*,520)titre ; WRITE(2,520)titre
520	    FORMAT('nom du fichier du modèle de coox en binaire : ',a)
	   END SELECT
	   titre=TRIM(nom_fich2)//'_B.dat'
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1511)titre ; WRITE(2,1511)titre
1511	    FORMAT('filename of the  final binary model : ',a)
	   CASE DEFAULT
	    WRITE(*,511)titre ; WRITE(2,511)titre
511	    FORMAT('nom du fichier du modèle final en binaire : ',a)
	   END SELECT
	   titre=TRIM(nom_fich2)//'_B.rep'
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1512)titre ; WRITE(2,1512)titre
1512	    FORMAT('filename of the intermediate binary model : ',a)
	   CASE DEFAULT
	    WRITE(*,512)titre ; WRITE(2,512)titre
512	    FORMAT('nom du fichier du modèle intermédiaire en binaire : ',a)
	   END SELECT
	   titre=TRIM(nom_fich2)//'.HR'
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1513)titre ; WRITE(2,1513)titre
1513	    FORMAT('filename for he plots of the HR and ZC diagrams : ',a)
	   CASE DEFAULT
	    WRITE(*,513)titre ; WRITE(2,513)titre
513	    FORMAT('nom du fichier pour dessins du diag. HR et ZC : ',a)
	   END SELECT
	   titre=TRIM(nom_fich2)//'-ad.osc'
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1517)titre ; WRITE(2,1517)titre
1517	    FORMAT('filename for the computation of adiabatic osc. : ',a)
	   CASE DEFAULT
	    WRITE(*,517)titre ; WRITE(2,517)titre
517	    FORMAT('nom du fichier pour le calcul des osc. adia. : ',a)
	   END SELECT
	   titre=TRIM(nom_fich2)//'-nad.osc'
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1518)titre ; WRITE(2,1518)titre
1518	    FORMAT('filename for the computation of non ad. osc. : ',a)
	   CASE DEFAULT
	    WRITE(*,518)titre ; WRITE(2,518)titre
518	    FORMAT('nom du fichier pour le calcul des osc. non adia. : ',a)
	   END SELECT
	   titre=TRIM(nom_fich2)//'-inv.osc'
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1519)titre ; WRITE(2,1519)titre
1519	    FORMAT('filename for inversion computations : ',a)
	   CASE DEFAULT
	    WRITE(*,519)titre ; WRITE(2,519)titre
519	    FORMAT('nom du fichier pour le calcul des inversions : ',a)
	   END SELECT
	   titre=TRIM(nom_fich2)//'_B.atm'
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1514)titre ; WRITE(2,1514)titre
1514	    FORMAT('filename of the binary model atmosphere : ',a)
	   CASE DEFAULT
	    WRITE(*,514)titre ; WRITE(2,514)titre
514	    FORMAT('nom du fichier du modèle d''atmosphère en binaire : ',a)
	   END SELECT
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1515) ; WRITE(2,1515)
1515	    FORMAT('The stellar evolution code CESAM has been elaborated',/,
	1   'in the framework of the research group',/,
	2   'Internal Structure of Stars and Giant Planets.',/,
	3   'If you found it useful the objectives of',/,
	4   'all CESAM contributors has been reached',//,
	5   'P.Morel, ON. Décembre 1989, CESAM1',/,
	6   'P.Morel, OCA. Octobre 1991, CESAM2',/,
	7   'P.Morel, OCA. Avril 1993, CESAM3',/,
	8   'P.Morel, OCA. Décembre 1997, CESAM4',/,
	9   'P.Morel, OCA. Décembre 2002, CESAM5',/,
	1   'P.Morel, B.Pichon OCA. Septembre 2003, CESAM2k',/,
	2   'English translation : ',/,
	3   'P.Morel, OCA., Y.Lebreton, MJo.Goupil OBSPM, Février 2004',/,
	4   'P.Morel, OCA. Mars 2005, diffusion of the angular momentum',/,
	5   'P.Morel, OCA. Août 2006, towards advanced stellar evolution',//,
	6   '***********************************',/)
	   CASE DEFAULT
	    WRITE(*,515) ; WRITE(2,515)
515	    FORMAT('Le code d''évolution stellaire CESAM a été élaboré',/,
	1   'dans le cadre du Groupement de Recherche Structure Interne',/,
	2   'des Etoiles et des Planètes Géantes. Si son utilisation ',/,
	3   'vous a donné satisfaction, le but poursuivi par tous',/,
	4   'ceux qui y ont contribué aura été atteint.',//,
	5   'P.Morel, ON. Décembre 1989, CESAM1',/,
	6   'P.Morel, OCA. Octobre 1991, CESAM2',/,
	7   'P.Morel, OCA. Avril 1993, CESAM3',/,
	8   'P.Morel, OCA. Décembre 1997, CESAM4',/,
	9   'P.Morel, OCA. Décembre 2002, CESAM5',/,
	1   'P.Morel, B.Pichon OCA. Septembre 2003, CESAM2k',/,
	3   'P.Morel, OCA, Y.Lebreton, MJo Goupil OBSPM, Février 2004, Version anglaise',/,
	4   'P.Morel, OCA, A.Moya OBSPM, Mars 2005, diff. du moment cinétique',/,
	5   'P.Morel, OCA. Août 2006, vers les stades avancés',/
	6   'B.Pichon, P.Morel OCA. Décembre 2008, lois T(tau) de MARCS',/,
	7   'Y.Lebreton OBSPM, P.Morel OCA. Décembre 2008, mixtures+opacités AGS05',//,	
	8   '***********************************',/)
	   END SELECT

c on garde le modèle final, en fixant le pas temporel extrapolé
c dts au pas dt
	   dts=dt
	   IF(all_rep)THEN
	    WRITE(number,70)model_num
	    modelb=TRIM(nom_fich2)//number//'_B.dat'
	   ELSE
	    modelb=TRIM(nom_fich2)//'_B.dat'
	   ENDIF
	   OPEN(unit=26,form='unformatted',status='unknown',file=TRIM(modelb))
	   WRITE(26)ne,m_qs,n_qs,knot,nchim,n_ch,m_ch,knotc,Krot,nrot,n_rot,
	1  m_rot,knotr,n_tds,knot_tds,mtot,alpha,w_rot,lim_ro,diffusion,
	2  rot_solid,precision,en_masse,f_eos,f_opa,
	3  nom_ctes,nom_pertm,nom_pertw,nom_tdetau,nom_atm,nom_conv,nom_nuc,
	4  nom_nuc_cpl,nom_diffm,nom_difft,nom_diffw,nom_etat,nom_opa,
	5  nom_elem,bp,q,qt,chim,mc,mct,rota,mrot,mrott,tds,x_tds,xt_tds,
	6  m23,r2,m_zc,r_zc,r_ov,age,dt,dts,mstar,rstar,mw_tot,wrot,
	7  jlim,lconv,lim,model_num
	   CLOSE(unit=26)
c	   PRINT*,all_rep,model_num,modelb ; WRITE(*,2000)dt,dts ; PAUSE'.dat'

c dernière écriture sur le fichier mon_modèle.HR et fermeture
	   WRITE(53,52)age,nchim,Krot,lim,model_num,lconv(1:lim)
	   WRITE(53,2002)LOG10(teff),LOG10(l(n_qs)),LOG10(rstar),mstar,
	1  (mstar-m_zc(i),r_zc(i),r_ov(i),i=1,lim)
	   WRITE(53,2003)(nom_elem(i),compg(i,1),compg(i,n_qs),i=1,nchim)
	   IF(Krot >= 3)THEN	!pour Krot = 3, 4 et 5
	    WRITE(53,2003)nom_vwrot,rstar*rsol*1.d-5*w(n_qs),w(n_qs)
	   ENDIF
	   CLOSE(unit=53)

c fermetures du fichier .lis (2)
	   CLOSE(unit=2) ; EXIT Bev
	  ENDIF           !on sort
	 ENDIF            !ecrit

c output ASCII éventuel
	 IF(all_output)CALL ecrit_ascii

c Alerte avec diffusion du moment cinétique si \tilde g \ g > 20%
	 IF(COUNT(Krot == Ktest) >= 1)THEN
	  IF(gtmax > 0.2d0)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1808)nom_thw,gtmax,igtmax,rgtmax,mgtmax
	    WRITE(2,1808)nom_thw,gtmax,igtmax,rgtmax,mgtmax
	   CASE DEFAULT
	    WRITE(*,808)nom_thw,gtmax,igtmax,rgtmax,mgtmax
	    WRITE(2,808)nom_thw,gtmax,igtmax,rgtmax,mgtmax
	   END SELECT
	  ELSE
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1809)gtmax,igtmax,rgtmax,mgtmax
	    WRITE(2,1809)gtmax,igtmax,rgtmax,mgtmax
	   CASE DEFAULT
	    WRITE(*,809)gtmax,igtmax,rgtmax,mgtmax
	    WRITE(2,809)gtmax,igtmax,rgtmax,mgtmax
	   END SELECT
	  ENDIF
	 ENDIF

c écriture, totale si ecrit, partielle si .not.ecrit
c	 PRINT*,ecrit ; PAUSE'ecrit'
	 CALL list(alfa,anub8,anube7,anun13,anuo15,anupep,anupp,beta,
	1 compg,cp,delta,dcapdr,dcapdt,depsdr,depsdt,d2p,d2ro,
	2 chaine,convec,ecrit,epsilon,gamma,gamma_atm,gradad,grada_atm,
	3 gradconv,gradc_atm,gradrad,gradr_atm,hp,i_cno,i_gr,i_pp,i_3a,
	4 kap,l,m,mu,mue,m_atm,p,pt,pt_atm,p_atm,r,ro,ro_atm,r_atm,t,tau,
	5 teff,tx_ioni,t_atm,u,vaissala,w,z,degene)

c sorties sur écran
c	 PRINT*,mstar,mtot
	 pertmm=mstar/mtot-1.d0
	 IF(ihe4 <= 1)THEN
	  WRITE(*,205)age,LOG10(teff),LOG10(l(n_qs)),LOG10(rstar),
	1 LOGg+LOG10(mstar/rstar**2),pt(1),t(1),ro(1),compg(1,1),
	2 i_pp,i_cno,i_3a,i_gr,pertmm,mstar,chaine
	ELSEIF(w(n_qs) > 0.d0)THEN
	 WRITE(*,207)age,LOG10(teff),LOG10(l(n_qs)),LOG10(rstar),
	1 LOGg+LOG10(mstar/rstar**2),pt(1),t(1),ro(1),compg(1,1),
	2 i_pp,i_cno,i_3a,i_gr,compg(ihe4,1),pertmm,mstar,TRIM(chaine),
	3 2.d0*pi/24.d0/3600.d0/w(n_qs),rstar*w(n_qs)*rsol*1.d-5
	 ELSE
	  WRITE(*,206)age,LOG10(teff),LOG10(l(n_qs)),LOG10(rstar),
	1 LOGg+LOG10(mstar/rstar**2),pt(1),t(1),ro(1),compg(1,1),
	2 i_pp,i_cno,i_3a,i_gr,compg(ihe4,1),pertmm,mstar,chaine
	 ENDIF

c écriture sur le fichier mon_modèle.HR
	 IF(ecrHR)THEN
	  WRITE(53,52)age,nchim,Krot,lim,model_num,lconv(1:lim)
52	  FORMAT(es22.15,3i3,i4,20l2)
	  WRITE(53,2002)LOG10(teff),LOG10(l(n_qs)),LOG10(rstar),mstar,
	1 (mstar-m_zc(i),r_zc(i),r_ov(i),i=1,lim)
2002	  FORMAT(6es13.6)
	  WRITE(53,2003)(nom_elem(i),compg(i,1),compg(i,n_qs),i=1,nchim)
2003	  FORMAT(a4,2es12.5)
	  IF(Krot >= 3)THEN	!pour Krot=3, 4 et 5
	   WRITE(53,2003)nom_vwrot,rstar*rsol*1.d-5*w(n_qs),w(n_qs)
	  ENDIF
	 ELSE
	  ecrHR=.TRUE.
	 ENDIF

c avant de poursuivre l'évolution on garde le modèle intermédiaire
c	 WRITE(*,2000)age,dt,dts ; PAUSE'age,dt,dts'
	 IF(.NOT.(list_rep .OR. list_pms .OR. list_zams))THEN
	  IF(all_rep)THEN
	   WRITE(number,70)model_num
	   modelb=TRIM(nom_fich2)//number//'_B.rep'
	  ELSE
	   modelb=TRIM(nom_fich2)//'_B.rep'
	  ENDIF
	  OPEN(unit=26,form='unformatted',status='unknown',file=TRIM(modelb))
	  WRITE(26)ne,m_qs,n_qs,knot,nchim,n_ch,m_ch,knotc,Krot,nrot,n_rot,
	1 m_rot,knotr,n_tds,knot_tds,mtot,alpha,w_rot,lim_ro,diffusion,
	2 rot_solid,precision,en_masse,f_eos,f_opa,
	3 nom_ctes,nom_pertm,nom_pertw,nom_tdetau,nom_atm,nom_conv,nom_nuc,
	4 nom_nuc_cpl,nom_diffm,nom_difft,nom_diffw,nom_etat,nom_opa,
	5 nom_elem,bp,q,qt,chim,mc,mct,rota,mrot,mrott,tds,x_tds,xt_tds,
	6 m23,r2,m_zc,r_zc,r_ov,age,dt,dts,mstar,rstar,mw_tot,wrot,
	7 jlim,lconv,lim,model_num
	  CLOSE(unit=26)
c	  PRINT*,all_rep,model_num,modelb ; WRITE(*,2000)dt,dts; PAUSE'.rep'
	 ELSE
	  list_pms=.FALSE. ; list_rep=.FALSE. ; list_zams=.FALSE.
	 ENDIF

c	 PAUSE'modèle suivant'
	 
c pas temporel suivant	 
	 CALL resout(1,dt,dts) ; age=age+dt ; dtp=dt

	ENDDO Bev

c	PAUSE'fin de cesam'

	IF(age > 0.d0)THEN
	 CALL des(.TRUE.,dtp,teff)
	ELSE
	 CALL des(.TRUE.,0.d0,teff)
	ENDIF

	STOP

	CONTAINS
	 INCLUDE'ecrit_ascii.f'

	END SUBROUTINE cesam
