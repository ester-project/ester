
c******************************************************************

	SUBROUTINE lit_nl(wrot)

c routine public des modules mod_donnees et mod_exploit

c lecture des NAMELISTs du fichier de données *.don

c le fichier *.don doit se trouver dans le directory où sont
c effectués les calculs

c variables numériques
c	mtot: masse totale initiale
c	mdot: taux de perte de masse Msol/an
c	tau_max: épaisseur optique au fond de l'atmosphère
c	agemax: âge maximun a atteindre
c	dtlist: intervalle de temps min. entre deux listings
c	complets du modèle
c	x_stop: arrêt si X centre < = xstop
c	log_teff: arrêt si Log(Teff) <> +/-log_teff
c	t_stop: arrêt si T centre > t_stop
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
c	nom_pertm: nom la routine de perte de masse à utiliser
c	nom_pertw: nom la routine de perte de moment cinétique à utiliser
c	nom_tdetau: nom la routine de loit T(tau) à utiliser
c	nom_atm: nom la routine de restitution d'atmosphère à utiliser
c	nom_conv: nom la routine de convection à utiliser
c	nom_nuc: nom la routine de réactions nucléaires à utiliser
c	nom_nuc_cpl: nom la routine de compilation de reac. nuc. à utiliser
c	nom_abon: nom la routine d'abondances initiales à utiliser
c	nom_diffm: nom la routine de diffusion microscopique à utiliser
c	nom_difft: nom la routine de  diffusion turbulente à utiliser
c	nom_diffw: nom la routine de  diffusion du moment cinétique à utiliser
c	nom_thw: nom la théorie de la diffusion du moment cinétique à utiliser
c	nom_etat: nom la routine d'équation d'état à utiliser
c	nom_opa: nom la routine d'opacité à utiliser

c variables sous forme de chaines de caractères de noms de fichiers
c	f_eos: noms des fichiers d'équation d'état
c	f_opa: noms des fichiers d'opacité

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c---------------------------------------------------------------------

	USE mod_kind
	USE mod_numerique, ONLY : pause

	IMPLICIT NONE

	REAL (kind=dp), INTENT(out) :: wrot
	LOGICAL :: ok

	CHARACTER (len=20) :: nom_des_rot
	CHARACTER (len=50) :: chain

c NAMELISTs de CESAM2k avec diffusion du moment cinétique
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

2000	FORMAT(8es10.3)

	chain=TRIM(nom_fich2)//'.don'
	INQUIRE(file=TRIM(chain),exist=ok)
	IF(ok)THEN
	 OPEN(unit=3,form='formatted',status='old',delim='apostrophe',
	1 file=TRIM(chain))
	ELSE
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1010)TRIM(chain) ; WRITE(2,1010)TRIM(chain)
1010	  FORMAT('STOP, unknown input file in lit_nl : ',a)
	 CASE DEFAULT
	  WRITE(*,10)TRIM(chain) ; WRITE(2,10)TRIM(chain)
10	  FORMAT('ARRET, fichier de données inconnu dans lit_nl : ',a)
	 END SELECT
	 STOP
	ENDIF

	SELECT CASE(langue)
	CASE('english')
	 WRITE(*,1001)chain ; WRITE(2,1001)chain
1001	 FORMAT(t10,'NAMELISTS of the input file : ',a,/)
	CASE DEFAULT
	 WRITE(*,1)chain ; WRITE(2,1)chain
1	 FORMAT(t10,'NAMELISTS du fichier: ',a,/)
	END SELECT

	READ(3,nl_cesam,ERR=100,END=100)
	WRITE(*,nl_cesam) ; WRITE(2,nl_cesam)
	READ(3,nl_mass,ERR=100,END=100) ; WRITE(*,nl_mass) ; WRITE(2,nl_mass)
	READ(3,nl_evol,ERR=100,END=100) ; WRITE(*,nl_evol) ; WRITE(2,nl_evol)
	READ(3,nl_chim,ERR=100,END=100) ; WRITE(*,nl_chim) ; WRITE(2,nl_chim)
	READ(3,nl_conv,ERR=100,END=100) ; WRITE(*,nl_conv) ; WRITE(2,nl_conv)
	READ(3,nl_diff,ERR=100,END=100) ; WRITE(*,nl_diff) ; WRITE(2,nl_diff)
	READ(3,nl_rot,ERR=100,END=100)  ; WRITE(*,nl_rot)  ; WRITE(2,nl_rot)
	READ(3,nl_etat,ERR=100,END=100) ; WRITE(*,nl_etat) ; WRITE(2,nl_etat)
	READ(3,nl_opa,ERR=100,END=100)  ; WRITE(*,nl_opa)  ; WRITE(2,nl_opa)
	READ(3,nl_nuc,ERR=100,END=100)  ; WRITE(*,nl_nuc)  ; WRITE(2,nl_nuc)
	READ(3,nl_atm,ERR=100,END=100)  ; WRITE(*,nl_atm)  ; WRITE(2,nl_atm)
	CLOSE(unit=3) ; WRITE(*,*); WRITE(2,*)
	
c limite JPZ pour la rotation supprimée
	IF(lim_jpz)THEN
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1030) ; WRITE(2,1030)
1030	  FORMAT('Boundary condition JPZ, unavailable, sorry',/,
	1 'use lim_jpz=.FALSE., STOP')
	 CASE DEFAULT
	  WRITE(*,30) ; WRITE(2,30)
30	  FORMAT('Condition limite JPZ indisponible, désolé',/,
	1 'utiliser lim_jpz=.FALSE., ARRET')
	 END SELECT
	 STOP
	ENDIF	

c tams est synonyme de post
	IF(arret == 'tams')arret='post'

c all_output=.TRUE. on écrit les fichiers ASCII	de TOUS les modèles
	all_output=(nom_output == 'all_adia') .OR. (nom_output == 'all_invers')
	1 .OR. (nom_output == 'all_nadia') .OR. (nom_output == 'all_ascii')

c all_rep=TRUE il y a écriture de TOUS les fichiers mon_modèle.rep
	all_rep=nb_max_modeles < 0 ; nb_max_modeles=ABS(nb_max_modeles)

c arrêt dès que T au centre dépasse t_stop
c t_stop en ln, on élimine les cas invraisemblables
	t_ajuste=t_stop > 1.d5
	IF(t_ajuste)THEN
	 lnt_stop=LOG(t_stop)
	ELSE
	 lnt_stop=LOG(1.d20)
	ENDIF

c grad_ov*=.TRUE. gradient adiabatique dans les parties overshootées des ZC
c grad_ov*=.FALSE. gradient radiatif dans les parties overshootées des ZC
	grad_ovi=ovshti > 0.d0 ; grad_ovs=ovshts > 0.d0
	ovshti=ABS(ovshti) ; ovshts=ABS(ovshts)

	IF(ovshti > 0.d0)THEN
	 IF(grad_ovi)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1021) ; WRITE(2,1021)
1021	   FORMAT('Use of the adiabatic gradient in undershooted parts')
	  CASE DEFAULT
	   WRITE(*,21) ; WRITE(2,21)
21	   FORMAT('Utilisation du grad. adia. dans les parties underhootées')
	  END SELECT
	 ELSE
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1022) ; WRITE(2,1022)
1022	   FORMAT('Use of the radiative gradient in undershooted parts')
	  CASE DEFAULT
	   WRITE(*,22) ; WRITE(2,22)
22	   FORMAT('Utilisation du grad.radia. dans les parties underhootées')
	  END SELECT
	 ENDIF
	ENDIF

	IF(ovshts > 0.d0)THEN
	 IF(grad_ovs)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1023) ; WRITE(2,1023)
1023	   FORMAT('Use of the adiabatic gradient in overshooted area')
	  CASE DEFAULT
	   WRITE(*,23) ; WRITE(2,23)
23	   FORMAT('Utilisation du grad. adia. dans les parties overhootées')
	  END SELECT
	 ELSE
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1024) ; WRITE(2,1024)
1024	   FORMAT('Use of the radiative gradient in overshooted area')
	  CASE DEFAULT
	   WRITE(*,24) ; WRITE(2,24)
24	   FORMAT('Utilisation du grad. radia dans les parties overhootées')
	  END SELECT
	 ENDIF
	ENDIF

c la composition chimique initiale élimination d'inconsistences
	IF(y0 < 0.d0)THEN
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1015)y0 ; WRITE(2,1015)y0
1015	  FORMAT('STOP,Y0=',es10.3)
	 CASE DEFAULT
	  WRITE(*,15)y0 ; WRITE(2,15)y0
15	  FORMAT('ARRET, incohérence : Y0=',es10.3)
	 END SELECT
	 STOP
	ENDIF

	IF(zsx0 <= 0.d0 .AND. x0 < 0.d0)THEN
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1014)zsx0,x0 ; WRITE(2,1014)zsx0,x0
1014	  FORMAT('STOP, unconsistency :  zsx0=',es10.3,', x0=',es10.3)
	 CASE DEFAULT
	  WRITE(*,14)zsx0,x0 ; WRITE(2,14)zsx0,x0
14	  FORMAT('ARRET, incohérence : zsx0=',es10.3,', x0=',es10.3)
	 END SELECT
	 STOP
	ENDIF

c détermination de la composition chimique initiale
	IF(zsx0 > 0.d0)THEN
	 x0=(1.d0-y0)/(1.d0+zsx0) ; z0=1.d0-x0-y0
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1003)chain,x0,y0,z0 ; WRITE(2,1003)chain,x0,y0,z0
1003	  FORMAT('X0 et Z0 are calculated from Y0 et Z/X0, read on the',/,
	1 'file : ',a,/,'X0=',es10.3,', Y0=',es10.3,', Z0=',es10.3)
	 CASE DEFAULT
	  WRITE(*,3)chain,x0,y0,z0 ; WRITE(2,3)chain,x0,y0,z0
3	  FORMAT('X0 et Z0 sont déduits de Y0 et Z/X0, du fichier: ',a,/,
	1 'X0=',es10.3,', Y0=',es10.3,', Z0=',es10.3)
	 END SELECT
	ELSEIF(x0 < 0.d0)THEN
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1004)TRIM(chain) ; WRITE(2,1004)TRIM(chain)
1004	  FORMAT('STOP, from the file : ',a,/,'X0, Y0 or ZSX0 < 0')
	 CASE DEFAULT
	  WRITE(*,4)TRIM(chain) ; WRITE(2,4)TRIM(chain)
4	  FORMAT('ARRET, dans le fichier : ',a,/,'X0, Y0 ou ZSX0 < 0')
	 END SELECT
	 STOP
	ELSEIF(garde_xish)THEN
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1013)y0,TRIM(nom_abon) ; WRITE(2,1013)y0,TRIM(nom_abon)
1013	  FORMAT('X0 and Z0 will be deduced from Y0=',es10.3,/,
	1  'and of the ratio Xi/H of the initial mixture : ',a)
	 CASE DEFAULT
	  WRITE(*,13)y0,TRIM(nom_abon) ; WRITE(2,13)y0,TRIM(nom_abon)
13	  FORMAT('X0 et Z0 seront déduits de Y0=',es10.3,/,
	1  'et du Xi/H de la mixture initiale : ',a)
	 END SELECT
	ELSE
	 z0=1.d0-x0-y0
	 SELECT CASE(langue)
	 CASE('english')
	  IF(x0 == 0.d0)CALL pause('WARNING : x0=0')
	  WRITE(*,1009)z0 ; WRITE(2,1009)z0
1009	  FORMAT('The initial metallicity is computed from X0 and Y0, Z0=',
	1   es10.3)
	 CASE DEFAULT
	  IF(x0 == 0.d0)CALL pause('ATTENTION : x0=0')
	  WRITE(*,9)z0 ; WRITE(2,9)z0
9	  FORMAT('Métallicité initiale déduite de X0 et Y0, Z0=',es10.3)
	 END SELECT
	ENDIF

c arrêt dès que X au centre devient inférieur à x_stop
	x_ajuste= x_stop > 0.d0 .AND. x_stop < x0

c arrêt dès que X devient < X_tams=0.01 au niveau m=He_core < Mtot/2
	he_ajuste=he_core > 0.d0 .AND. he_core < mtot/2.d0
	IF(he_ajuste)THEN
	 hhe_core=he_core**(2.d0/3.d0)
	ELSE
	 hhe_core=1.d-30
	ENDIF

c la convection
	IF(.FALSE.)THEN		!suppression
c	IF(ledoux .AND. .NOT.diffusion)THEN
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1005)chain ; WRITE(2,1005)chain
1005	  FORMAT('STOP, with ledoux=.TRUE. we must have diffusion=.TRUE.',
	1 /,'Modify the data file : ',a)
	 CASE DEFAULT
	  WRITE(*,5)chain ; WRITE(2,5)chain
5	  FORMAT('ARRET, avec ledoux=.TRUE. on doit avoir diffusion=.TRUE.',
	1 /,'Modifier le fichier de données : ',a)
	 END SELECT
	 STOP
	ENDIF

c cas de pp1
	IF(nom_nuc == 'pp1' .AND. diffusion)THEN
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1008)chain ; WRITE(2,1008)chain
1008	  FORMAT('STOP, NOM_NUC = pp1 and DIFFUSION = TRUE',/,
	1 'are unconsistant. Modify the data file : ',a)
	 CASE DEFAULT
	  WRITE(*,8)chain ; WRITE(2,8)chain
8	  FORMAT('ARRET, NOM_NUC = pp1 et DIFFUSION = TRUE incompatibles',
	1 /,'Modifier le fichier de données: ',a)
	 END SELECT
	 STOP
	ENDIF

c initialisation des constantes
	CALL ini_ctes

c choix entre les formalismes utilisés pour la rotation
c initialisation de Krot : flag pour le traitement de la rotation
	SELECT CASE(nom_thw)
	CASE('rot_0')
	 IF(w_rot /= 0.d0)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1025)w_rot,TRIM(nom_thw)
1025	   FORMAT('WARNING, inconsistency, w_rot=',es10.3,' with NOM_THW=',a,/,
	1  'the model is computed wihout rotation' )
	  CASE DEFAULT
	   WRITE(*,25)w_rot,TRIM(nom_thw)
25	   FORMAT('ATTENTION, w_rot=',es10.3,' avec NOM_THW=',a,/,
	1  'le modèle sera calculé sans rotation')
	  END SELECT
	  w_rot=0.d0
	 ENDIF
	 Krot=0 ; nrot=0 ; rot_solid=.TRUE. ; nom_des_rot='no_des'

	CASE('rot_cte')
	 IF(w_rot <= 0.d0)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1025)w_rot ,TRIM(nom_thw)
	  CASE DEFAULT
	   WRITE(*,25)w_rot,TRIM(nom_thw)
	  END SELECT
	  w_rot=0.d0 ; Krot=0
	 ELSE
	  Krot=1
	 ENDIF
	 nrot=0 ; rot_solid=.TRUE. ; nom_des_rot='no_des'

	CASE('cons_glob_mnt_cin')
	 IF(w_rot <= 0.d0)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1025)w_rot ,TRIM(nom_thw)
	  CASE DEFAULT
	   WRITE(*,25)w_rot,TRIM(nom_thw)
	  END SELECT
	  w_rot=0.d0 ; Krot=0
	 ELSE
	  Krot=2
	 ENDIF
	 nrot=0 ; rot_solid=.TRUE. ; nom_des_rot='no_des'

	CASE('diff_tz97')
	 IF(w_rot <= 0.d0)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1110)w_rot ,TRIM(nom_thw)
1110	   FORMAT('STOP, inconsistency, w_rot=',es10.3,' with NOM_THW=',a)
	  CASE DEFAULT
	   WRITE(*,110)w_rot,TRIM(nom_thw)
110	   FORMAT('ARRET, inconsistance, w_rot=',es10.3,' avec NOM_THW=',a)
	  END SELECT
	  STOP
	 ENDIF
	 Krot=3 ; nrot=7 ; nrl=24 ; rot_solid=.FALSE.
	 ALLOCATE(nom_rot(nrot),rot_min(nrot))

c noms et valeurs non significatives des variables de rotation
	 nom_rot=(/ 'Omeg',' U  ','Teta','Lamb','Psi ','Flux','Upsi' /)
	 rot_min=(/ 1.d-7, 1.d-5, 1.d-3, 1.d-2, 1.d-2, 1.d-13, 1.d-13 /)

	CASE('diff_mz04')
	 IF(w_rot <= 0.d0)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1110)w_rot ,TRIM(nom_thw)
	  CASE DEFAULT
	   WRITE(*,110)w_rot,TRIM(nom_thw)
	  END SELECT
	  STOP
	 ENDIF
	 Krot=4 ; nrot=8 ; nrl=32 ; rot_solid=.FALSE.
	 ALLOCATE(nom_rot(nrot),rot_min(nrot))

c noms et valeurs non significatives des variables de rotation
	 nom_rot=(/ 'Omeg',' U  ','Psi ','Lamb','Flux','Upsi', 'Phi ',' Pi '/)
	 rot_min=(/ 1.d-7, 1.d-5, 1.d-8, 1.d-7, 1.d-13, 1.d-13, 1.d8, 1.d8 /)

	CASE('cons_loc_mnt_cin')
	 IF(w_rot <= 0.d0)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1110)w_rot ,TRIM(nom_thw)
	  CASE DEFAULT
	   WRITE(*,110)w_rot,TRIM(nom_thw)
	  END SELECT
	  STOP
	 ENDIF
	 Krot=5 ; nrot=1 ; rot_solid=.FALSE. ; nom_des_rot='no_des'
	 ALLOCATE(nom_rot(nrot),rot_min(nrot))
	 nom_rot=(/ 'Omeg' /)

	CASE DEFAULT
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1026)TRIM(nom_thw) ; WRITE(2,1026)TRIM(nom_thw)
1026	  FORMAT('ERROR, in the data file NOM_THW : ',a,/,
	1 'is unknown, known values : no_rot, cte, cons_glob_mnt_cin',
	2 /,'diff_tz97, diff_mz04, cons_loc_mnt_cin')

	 CASE DEFAULT
	  WRITE(*,26)TRIM(nom_thw) ; WRITE(2,26)TRIM(nom_thw)
26	  FORMAT('ERREUR, dans le fichier de données NOM_THW : ',a,/,
	1 'est inconnu, valeurs connues : rot_0, rot_cte, cons_glob_mnt_cin',
	2 /,'diff_tz97, diff_mz04, cons_loc_mnt_cin')
	  STOP
	 END SELECT
	END SELECT

c type de conditions pour U aux limites	ZR/ZC pour Krot=3,4
	SELECT CASE(Krot)
	CASE(3,4)
	 IF(lim_jpz)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1027) ; WRITE(2,1027)
1027	   FORMAT('Boundary condition for U at limits ZR/ZC of JPZh')
	  CASE DEFAULT
	   WRITE(*,27) ; WRITE(2,27)
27	   FORMAT('Conditions limites pour U aux limites ZR/ZC de JPZh')
	  END SELECT
	 ELSE
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1028) ; WRITE(2,1028)
1028	   FORMAT('Boundary condition for U at limits ZR/ZC of PM')
	  CASE DEFAULT
	   WRITE(*,28) ; WRITE(2,28)
28	   FORMAT('Conditions limites pour U aux limites ZR/ZC de PM')
	  END SELECT
	 ENDIF
	 
c avec rotation il vaut mieux utiliser la diffusion microscopique	 	
	 IF(.NOT.diffusion)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1029) ; WRITE(2,1029)
1029	   FORMAT('With rotation, use also microscopic diffusion')
	  CASE DEFAULT
	   WRITE(*,29) ; WRITE(2,29)
29	   FORMAT('Avec rotation il est préférable d''utiliser la diffusion')
	  END SELECT
	 ENDIF 	 
	END SELECT

c w_rot est la valeur initiale de la vitesse angulaire, avec Krot=1, 2
c wrot deviendra la valeur évoluée de la vitesse angulaire dans les
c unités, avec unit=kms/s, w_rot est la vitesse de surface en
c kilomètres par seconde, la vitesse angulaire initiale (rd/s) en sera
c déduite dans la routine cesam dès que le rayon initial sera connu
c pour les modèles initiaux de ZAMS et PMS
	IF(Krot == 0)THEN
	 wrot=w_rot
	ELSE
	 SELECT CASE(unit)
	 CASE('rad/s','kms/s')
	  wrot=w_rot
	 CASE('jours')
	   wrot=2.d0*pi/24.d0/3600.d0/w_rot ; w_rot=wrot
	 CASE DEFAULT
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1012)unit
1012	   FORMAT('ERROR on the input of angular velocity unit : ',a,/,
	1 'known units : rad/s, jours, kms/s')
	  CASE DEFAULT
	   WRITE(*,12)unit
12	   FORMAT('ERREUR d''unité de vitesse angulaire',a,/,
	1 'unités reconnues : rad/s, jours, kms/s')
	  END SELECT
	  STOP
	 END SELECT
	ENDIF

c affectations de Kdes_rot
c no_des : pas de fichier pour dessin
c end_evol : formation du fichier pour dessin à l'issue de l'évolution
c all_mod : formation des fichiers avec le numéro du modèle, conservation de
c tous les fichiers
c all_iter : pour mise au point le fichier de dessin est créé à chaque itération
c end_mod : formation du fichier pour dessin pour chaque modèle
	IF(TRIM(nom_des_rot) == 'no_des')THEN
	 Kdes_rot=0
	ELSEIF(TRIM(nom_des_rot) == 'end_evol')THEN
	 Kdes_rot=1
	ELSEIF(TRIM(nom_des_rot) == 'all_mod')THEN
	 Kdes_rot=2
	ELSEIF(TRIM(nom_des_rot) == 'all_iter')THEN
	 Kdes_rot=3
	ELSEIF(TRIM(nom_des_rot) == 'end_mod')THEN
	 Kdes_rot=4
	ELSE
	 Kdes_rot=0
	ENDIF

c avec n_max < 0 les derniers modèles d'une évolution seront calculés avec
c |n_max| couches
	nc_max = n_max < 0 ; n_max=ABS(n_max)
	
	RETURN

c fichier d'entrée de CESAM2k (mon_modele.don) est incorrect
100	SELECT CASE(langue)
	CASE('english')
	 WRITE(*,1020)TRIM(chain)
1020	 FORMAT('STOP, mistake in the input file : ',a,/,
	1 'see the example in the sub-directory EXPLOIT')
	CASE DEFAULT
	 WRITE(*,20)TRIM(chain)
20	 FORMAT('ARRET, le fichier de données : ',a,/,
	1 'est incorrect, un exemple se trouve dans le sous-directory EXPLOIT')
	END SELECT
	CLOSE(unit=3)

	STOP

	END SUBROUTINE lit_nl
