
c******************************************************************

	SUBROUTINE lit_nl(wrot)

c	routine public des modules mod_donnees et mod_exploit

c	lecture des NAMELISTs du fichier de donn�es *.don

c	le fichier *.don doit se trouver dans le directory o� sont
c	effectu�s les calculs

c	variables num�riques

c	cmtot: masse totale initiale
c	mdot: taux de perte de masse Msol/an
c	tau_max: �paisseur optique au fond de l'atmosph�re
c	agemax: �ge maximun a atteindre
c	dtlist: intervalle de temps min. entre deux listings
c	complets du mod�le
c	x_stop: arr�t si X centre < = xstop
c	log_teff: arr�t si Log(Teff) <> +/-log_teff
c	t_stop: arr�t si T centre > t_stop
c	x0: abondance de X initiale
c	y0: abondance de Y initiale
c	zsx0: rapport Z/X initial zsx0=0 X0 est utilise
c	d_turb: coefficient de diffusion turbulente
c	re_nu: param�tre de diffusivite radiative
c	w_rot: vitesse angulaire initiale
c	pmw : param�tre libre de perte de moment cin�tique
c	alpha: param�tre de longueur de melange
c	cpturb: coefficient de pression turbulente
c	ovshts: coefficient d'overshoot sup�rieur
c	ovshti: coefficient d'overshoot inf�rieur

c variables logiques

c	lim_ro=.true.: condition limite externe de l'atmosph�re sur ro
c	grille_fixe=.true.: utilisation d'une grille fixe en masse
c	pour interpolation de la comp.chim.
c	rot_solid=.true.: la rotation est solide
c	jpz=.true.: overshoot selon JpZ
c	ledoux=.true.: utilisation du crit�re de Ledoux
c	diffusion=.true.: calcul avec diffusion
c	mitler=.true.: effet d'�cran selon Mitler

c variables sous forme de chaines de caract�res

c	precision: niveau de pr�cision requis
c	arret: arr�t sur zams, post, cohe, coca ou autre
c	nom_des: nom de la routine de dessin on line � utiliser
c	nom_ctes_phys: nom la routine de constantes physiques � utiliser
c	nom_pertm: nom la routine de perte de masse � utiliser
c	nom_pertw: nom la routine de perte de moment cin�tique � utiliser
c	nom_tdetau: nom la routine de loit T(tau) � utiliser
c	nom_atm: nom la routine de restitution d'atmosph�re � utiliser
c	nom_conv: nom la routine de convection � utiliser
c	nom_nuc: nom la routine de r�actions nucl�aires � utiliser
c	nom_nuc_cpl: nom la routine de compilation de reac. nuc. � utiliser
c	nom_abon: nom la routine d'abondances initiales � utiliser
c	nom_diffm: nom la routine de diffusion microscopique � utiliser
c	nom_difft: nom la routine de  diffusion turbulente � utiliser
c	nom_etat: nom la routine d'�quation d'�tat � utiliser
c	nom_opa: nom la routine d'opacit� � utiliser

c variables sous forme de chaines de caract�res de noms de fichiers

c	f_eos: noms des fichiers d'�quation d'�tat
c	f_opa: noms des fichiers d'opacit�

c	Auteur: P.Morel, D�partement J.D. Cassini, O.C.A.
c	CESAM2k

c	Modifications:
c	29 09 00 : introduction de zsx0
c	22 12 00 : passage en F95

c---------------------------------------------------------------------

	USE mod_kind
	USE mod_numerique, ONLY : pause
      
	IMPLICIT NONE

	REAL (kind=dp), INTENT(out) :: wrot	
	LOGICAL :: ok

	CHARACTER (len=50) :: chain

c NAMELISTs de CESAM2k avec diffusion du moment cin�tique		
	NAMELIST/nl_cesam/nom_chemin,nom_ctes,nom_des,nom_output,
	1 n_max,precision
	NAMELIST/nl_mass/mtot,nom_pertm,mdot
	NAMELIST/nl_evol/agemax,arret,dtlist,log_teff,nb_max_modeles,
	1 he_core,t_stop,x_stop	
	NAMELIST/nl_chim/grille_fixe,nom_abon,modif_chim,garde_xish,
	1 x0,y0,zsx0
c	NAMELIST/nl_conv/nom_conv,alpha,ovshts,ovshti,beta_cgm,zeta_cgm,f_cgm,
c	1 jpz,cpturb,ledoux	
	NAMELIST/nl_conv/nom_conv,alpha,ovshts,ovshti,
	1 jpz,cpturb,ledoux,lov_ad	
	NAMELIST/nl_diff/diffusion,nom_diffm,nom_difft,d_turb,re_nu,
	1 nom_frad,beta_v,zeta_v,frac_vp
	NAMELIST/nl_rot/rot_solid,w_rot,unit,w_form,nom_diffw,nom_pertw,
	1 p_pertw, pw_extend
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
10	  FORMAT('ARRET, fichier de donn�es inconnu dans lit_nl : ',a)
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
		
	READ(3,nl_cesam, ERR=100, END=100)
	WRITE(*,nl_cesam) ; WRITE(2,nl_cesam)	
	READ(3,nl_mass) ; WRITE(*,nl_mass) ; WRITE(2,nl_mass)	
	READ(3,nl_evol) ; WRITE(*,nl_evol) ; WRITE(2,nl_evol)	
	READ(3,nl_chim) ; WRITE(*,nl_chim) ; WRITE(2,nl_chim)	
	READ(3,nl_conv) ; WRITE(*,nl_conv) ; WRITE(2,nl_conv)	
	READ(3,nl_diff) ; WRITE(*,nl_diff) ; WRITE(2,nl_diff)	
	READ(3,nl_rot)  ; WRITE(*,nl_rot)  ; WRITE(2,nl_rot)	
	READ(3,nl_etat) ; WRITE(*,nl_etat) ; WRITE(2,nl_etat)	
	READ(3,nl_opa)  ; WRITE(*,nl_opa)  ; WRITE(2,nl_opa)	
	READ(3,nl_nuc)  ; WRITE(*,nl_nuc)  ; WRITE(2,nl_nuc)	
	READ(3,nl_atm)  ; WRITE(*,nl_atm)  ; WRITE(2,nl_atm)	
	GOTO 110
	
c fichier d'entr�e de CESAM2k sans diffusion du moment cin�tique, �mulateur 	
100	CALL lit_nl_2korg(ok)
	IF(ok)THEN
	 GOTO 110
	ELSE
	 CALL lit_nl_45(ok)
	 IF(ok)THEN
	  GOTO 110	
	 ELSE
	  PRINT*,'ARRET, le fichier de donn�es n''est pas de type connu'
	  CLOSE(unit=3) ; STOP
	 ENDIF		
	ENDIF
110	CLOSE(unit=3) ; WRITE(*,*); WRITE(2,*)

c tams est synonyme de post

	 IF(arret == 'tams')arret='post'

c arr�t d�s que T au centre d�passe t_stop	 
c t_stop en ln, on �limine les cas invraisemblables
	 t_ajuste=t_stop > 1.d5
	 IF(t_ajuste)THEN
	  lnt_stop=LOG(t_stop)
	 ELSE
	  lnt_stop=LOG(1.d20) 
	 ENDIF

c la composition chimique initiale

c	�limination des inconsistences pour la d�termination de
c	la composition chimique initiale

	IF(y0 < 0.d0)THEN
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1015)y0 ; WRITE(2,1015)y0	 
1015	  FORMAT('STOP,Y0=',es10.3)
	 CASE DEFAULT
	  WRITE(*,15)y0 ; WRITE(2,15)y0	 
15	  FORMAT('ARRET, incoh�rence : Y0=',es10.3)
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
14	  FORMAT('ARRET, incoh�rence : zsx0=',es10.3,', x0=',es10.3)
	 END SELECT
	 STOP
	ENDIF	

c d�termination de la composition chimique initiale
	IF(zsx0 > 0.d0)THEN
	 x0=(1.d0-y0)/(1.d0+zsx0) ; z0=1.d0-x0-y0
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1003)chain,x0,y0,z0 ; WRITE(2,1003)chain,x0,y0,z0	 
1003	  FORMAT('X0 et Z0 are calculated from Y0 et Z/X0, read on the',/,
	1 'file : ',a,/,'X0=',es10.3,', Y0=',es10.3,', Z0=',es10.3)
	 CASE DEFAULT
	  WRITE(*,3)chain,x0,y0,z0 ; WRITE(2,3)chain,x0,y0,z0	 
3	  FORMAT('X0 et Z0 sont d�duits de Y0 et Z/X0, du fichier: ',a,/,
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
13	  FORMAT('X0 et Z0 seront d�duits de Y0=',es10.3,/,
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
9	  FORMAT('M�tallicit� initiale d�duite de X0 et Y0, Z0=',es10.3)
	 END SELECT	
	ENDIF
	
c arr�t d�s que X au centre devient inf�rieur � x_stop
	x_ajuste= x_stop > 0.d0 .AND. x_stop < x0
	
c arr�t d�s que X devient < X_tams=0.01 au niveau m=He_core < Mtot/2
	he_ajuste=he_core > 0.d0 .AND. he_core < mtot/2.d0
	IF(he_ajuste)THEN
	 hhe_core=he_core**(2.d0/3.d0)
	ELSE
	 hhe_core=1.d-30
	ENDIF 	 	
	
c il y a perte ou gain de masse
	lvent=mdot /= 0.d0	
	
c la convection	
	IF(ledoux .AND. .NOT.diffusion)THEN
	 SELECT CASE(langue)
	 CASE('english')	 	
	  WRITE(*,1005)chain ; WRITE(2,1005)chain
1005	  FORMAT('STOP, with ledoux=.TRUE. we must have diffusion=.TRUE.',
	1 /,'Modify the data file : ',a)
	 CASE DEFAULT	
	  WRITE(*,5)chain ; WRITE(2,5)chain
5	  FORMAT('ARRET, avec ledoux=.TRUE. on doit avoir diffusion=.TRUE.',
	1 /,'Modifier le fichier de donn�es : ',a)
	 END SELECT
	 STOP	
	ENDIF
	
c sans acc�l�rations radiatives
	IF(.FALSE.)THEN
c	IF(.TRUE.)THEN
	IF(diffusion .AND. nom_frad /= 'no_frad')THEN	
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1018) ; WRITE(2,1018)
1018	  FORMAT('STOP, the radiative accelerations are not available') 
	 CASE DEFAULT
	  WRITE(*,18) ; WRITE(2,18)
18	  FORMAT('ARRET, les acc�l�rations radiatives sont indisponible')
	 END SELECT	 
	 STOP
	ENDIF
	ENDIF	

c la rotation

	IF(.NOT.rot_solid)THEN      !rotation non solide	 
	 IF(w_rot <= 0.d0)THEN
	  SELECT CASE(langue)
	  CASE('english')	 	
	   WRITE(*,1006)chain ; WRITE(2,1006)chain
1006	   FORMAT('STOP, with non rigid rotation we must have w_rot > 0',
	1  /,'Modify the data file : ',a)
	  CASE DEFAULT		 
	   WRITE(*,6)chain ; WRITE(2,6)chain
6	   FORMAT('ARRET, avec rotation non solide il faut w_rot > 0',/,
	1  'Modifier le fichier de donn�es : ',a)
	  END SELECT
	  STOP
	 ELSEIF(.NOT.diffusion)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1007)chain ; WRITE(2,1007)chain
1007	   FORMAT('STOP, with non rigid rotation we must have',/,
	1  'diffusion=.TRUE.',/,'Modify the data file : ',a)
	  CASE DEFAULT
	   WRITE(*,7)chain ; WRITE(2,7)chain 
7	   FORMAT('ARRET, Avec rotation non solide il faut diffusion=.TRUE.',/,
	1  'Modifier le fichier de donn�es: ',a)
	  END SELECT
	  STOP		
	 ENDIF
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
	1 /,'Modifier le fichier de donn�es: ',a)
	 END SELECT
	 STOP		
	ENDIF
	
c initialisation des constantes

	CALL ini_ctes
	
c initialisation de Krot : flag pour le traitement de la rotation

c  ==>	 si wrot = 0  : pas de rotation on pose Krot=0   
      
c avec rot_solid = .TRUE. :
c	 la vitesse angulaire locale sera la m�me en tout point du mod�le

c  ==>	 si wrot > 0, la vitesse angulaire restera �gale � w_rot pendant
c	 toute l'�volution. Dans ce cas on posera Krot=1

c  ==>	 si wrot < 0, le moment angulaire est globalement conserv�,
c	 la vitesse angulaire locale wrot > 0, la m�me en tout point
c	 du mod�le, variera au cours du temps; elle sera r�ajust�e
c	 au cours des it�rations � l'issue du calcul du moment angulaire
c	 total dont on assurera la conservation en faisant varier wrot.
c	 Dans ce cas on posera Krot=2

c avec rot_solid = .FALSE. :

c   avec diffusion= .FALSE.
c  ==>	 Cas sans objet, non pr�vu

c  avec diffusion= .TRUE.
c	 diffusion du moment angulaire dans les ZR, m�lange dans les ZC.
c	 Dans ce cas on posera Krot=3

c	WRITE(*,2000)w_rot ; PAUSE'lit nl1'	
	
	IF(w_rot == 0.d0)THEN
	 Krot=0 ; nrot=0 ; rot_solid=.TRUE. ; wrot=0.d0 ; RETURN
	ELSEIF(rot_solid)THEN
	 IF(w_rot > 0.d0)THEN
	  Krot=1 ; nrot=0
	 ELSE
	  Krot=2 ; nrot=0
	 ENDIF
	ELSE
	 IF(diffusion)THEN
	  IF(.FALSE.)THEN
c	  IF(.TRUE.)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1017) ; WRITE(2,1017)
1017	   FORMAT('STOP, the diffusion of the angular momentum is not available')
	   STOP 
	  CASE DEFAULT
	   WRITE(*,17) ; WRITE(2,17)
17	   FORMAT('ARRET, la diffusion du moment cin�tique est indisponible')
	   STOP
	  END SELECT
	  ENDIF
	   
	  Krot=3 ; nrot=5
	  ALLOCATE(nom_rot(nrot), rot_min(nrot))
	  nom_rot=(/ 'Omeg',' U  ','teta','lamb','psi ' /)
	  rot_min=(/ 1.d-7, 1.d-5, 1.d-3, 1.d-2, 1.d-2 /)
c	  SELECT CASE(langue)
c	  CASE('english')
c	   WRITE(*,1016) ; WRITE(2,1016)
c1016	   FORMAT('Diff. of angular momentum : grille_fixe=TRUE is needed') 
c	  CASE DEFAULT
c	   WRITE(*,16) ; WRITE(2,16)
c16	   FORMAT('Diff. du moment cin�tique : grille_fixe=TRUE est impos�')
c	  END SELECT	  
 	  
	 ELSE
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1011)rot_solid,diffusion,TRIM(chain)
	   WRITE(2,1011)rot_solid,diffusion,TRIM(chain)
1011	   FORMAT('STOP, solid rotation = ',l1,' with diffusion = ',l1,
	1  'is an unknown case in CESAM',/,'change the input data file ',a)
	  CASE DEFAULT
	   WRITE(*,11)rot_solid,diffusion,TRIM(chain)
	   WRITE(2,11)rot_solid,diffusion,TRIM(chain)
11	   FORMAT('ARRET, rotation solide = ',l1,' avec diffusion = ',l1,
	1  'est un cas non pr�vu dans CESAM',/, 
	2  'il convient de corriger le fichier de donn�es ',a) 
	  END SELECT
	  STOP
	 ENDIF
	ENDIF

c	w_rot est la valeur initiale de la vitesse angulaire
c	avec Krot=1,2
c	wrot deviendra la valeur �volu�e de la vitesse angulaire dans les
c	cas Krot=1 et Krot=2	

c	pmw(1) : q pour w_initial
c	pmw(2) : Omega pour l'input pertw_ptm de moment cin�tique
c	pmw(3) : pw limite pour pertw_ptm de moment cin�tique
c	pmw(4) : a pour pertw_sch (Schumanisch)
c	pmw(5) : gamma pour pert_loc de moment cin�tique

c 	unit�s, avec unit=kms/s, w_rot est la vitesse de surface en
c	kilom�tres par seconde, la vitesse angulaire initiale (rd/s) en sera
c	d�duite dans la routine cesam d�s que le rayon initial sera connu
c	pour les mod�les initiaux de ZAMS et PMS
	
	SELECT CASE(unit)
	CASE('rad/s')
	 wrot=ABS(w_rot)
	CASE('jours')
	 IF(w_rot /= 0.d0)THEN
	  wrot=2.d0*pi/24.d0/3600.d0/ABS(w_rot)
	 ELSE
	  wrot=ABS(w_rot)
	 ENDIF
	 w_rot=SIGN(wrot,w_rot)
	CASE('kms/s')
	 wrot=ABS(w_rot)
	CASE DEFAULT
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1012)unit
1012	  FORMAT('ERROR on the input of angular velocity unit : ',a,/,
	1  'known units : rad/s, jours, kms/s')
	 CASE DEFAULT
	  WRITE(*,12)unit
12	  FORMAT('ERREUR d''unit� de vitesse angulaire',a,/,
	1  'unit�s reconnues : rad/s, jours, kms/s')
	 END SELECT	 
	 STOP	 
	END SELECT
	
	RETURN
	
	CONTAINS

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	 SUBROUTINE lit_nl_2korg(ok)
	
c	    lecture du fichier de donn�es de CESAM2k initial
	
	  LOGICAL, INTENT(out) :: ok
	  
	  REAL (kind=dp) :: dws		
	  CHARACTER (len=20) :: nom_perte
	  LOGICAL :: der_num, f_rad
	
c NAMELISTs de CESAM2k sans diffusion du moment cin�tique, �mulateur	
	 NAMELIST/nl_noms/nom_ctes,nom_des,nom_perte,nom_atm,
	1  nom_tdetau,nom_abon,nom_diffm,nom_difft,nom_conv,nom_etat,
	2  nom_opa,nom_nuc,nom_nuc_cpl,nom_output,nom_chemin
	 NAMELIST/nl_espace/mtot,mdot,der_num,precision,n_max
	 NAMELIST/nl_atmosphere/tau_max,lim_ro
	 NAMELIST/nl_temps/agemax,dtlist,x_stop,log_teff,t_stop,arret,
	1  nb_max_modeles
	 NAMELIST/nl_chim/x0,y0,zsx0,garde_xish,d_turb,re_nu,diffusion,
	1  f_rad,grille_fixe
	 NAMELIST/nl_rot/w_rot,dws,rot_solid  
	 NAMELIST/nl_conv/alpha,cpturb,ovshts,ovshti,jpz,ledoux
	 NAMELIST/nl_etat/f_eos
	 NAMELIST/nl_opa/f_opa
	 NAMELIST/nl_nuc/mitler

c-----------------------------------------------------------------------------------
	
	 REWIND(unit=3)	
	 READ(3,nl_noms, ERR=101, END=101)
	 WRITE(*,nl_noms)      ; WRITE(2,nl_noms)
	 READ(3,nl_espace)     ; WRITE(*,nl_espace)     ; WRITE(2,nl_espace)
	 READ(3,nl_atmosphere) ; WRITE(*,nl_atmosphere) ; WRITE(2,nl_atmosphere)
	 READ(3,nl_temps)      ; WRITE(*,nl_temps)      ; WRITE(2,nl_temps)
	 READ(3,nl_chim)       ; WRITE(*,nl_chim)       ; WRITE(2,nl_chim)
	 READ(3,nl_rot)        ; WRITE(*,nl_rot)        ; WRITE(2,nl_rot)
	 READ(3,nl_conv)       ; WRITE(*,nl_conv)       ; WRITE(2,nl_conv)
	 READ(3,nl_etat)       ; WRITE(*,nl_etat)       ; WRITE(2,nl_etat)
	 READ(3,nl_opa)        ; WRITE(*,nl_opa)        ; WRITE(2,nl_opa)
	 READ(3,nl_nuc)        ; WRITE(*,nl_nuc)        ; WRITE(2,nl_nuc)

	 IF(.NOT.rot_solid)THEN      !rotation non solide
	  SELECT CASE(langue)
	  CASE('english')	 	
	   WRITE(*,1002)chain ; WRITE(2,1002)chain
1002	   FORMAT('STOP, non solid rotation is not operating.',
	1   /,'with this input data file : ',a)
	  CASE DEFAULT		 
	   WRITE(*,2)chain ; WRITE(2,2)chain
2	   FORMAT('ARRET, rotation non solide non impl�ment�e.',/,
	1  'avec ce fichier de donn�es : ',a)
	  END SELECT
	  ok=.FALSE.
	  
c �mulation des donn�es	 
	 ELSE
	  IF(f_rad)THEN
	   nom_frad='alecian2'
	  ELSE
	   nom_frad='no_frad'
	  ENDIF
	  IF(ABS(w_rot) > 1.d0)THEN
	   unit='jours'
	  ELSE
	   unit='rad/s'
	  ENDIF
	  nom_pertm='pertm_ext'	 
	  nom_diffw='diffw_0'
	  nom_pertw='pertw_0'
	  p_pertw=0.d0
	  pw_extend=0.d0	 
	  w_form=0.d0	 
	  ok=.TRUE.
	 ENDIF
	 He_core=-100.d0
	 	 
	 RETURN
	
101	 ok=.FALSE.
	 RETURN	
	 
	END SUBROUTINE lit_nl_2korg
	
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	 SUBROUTINE lit_nl_45(ok)
	
c	    lecture du fichier de donn�es de CESAM4 ou 5

	  REAL (kind=dp) :: d_mw
	
	  LOGICAL, INTENT(out) :: ok
	    
	  CHARACTER (len=80) :: f_nuc
	  
	  LOGICAL :: der_num, f_rad, oui, z_cte
	  
	  NAMELIST/nl_espace/mtot,mdot,der_num,precision,tau_max,lim_ro
	  NAMELIST/nl_temps/agemax,dtlist,x_stop,log_teff,t_stop,arret
	  NAMELIST/nl_chim/x0,y0,zsx0,grille_fixe,z_cte,d_turb,diffusion,
	1   f_rad,re_nu,mitler
	  NAMELIST/nl_rot/w_rot,rot_solid,d_mw	
	  NAMELIST/nl_conv/alpha,cpturb,ovshts,ovshti,jpz,ledoux
	  NAMELIST/nl_etat/f_eos
	  NAMELIST/nl_opa/f_opa
	  NAMELIST/nl_nuc/f_nuc
	  
	  NAMELIST/nl_physique/nom_chemin,nom_ctes,nom_des,nom_output,
	1   nom_pertm,nom_conv,nom_diffm,nom_difft,
	2   nom_etat,nom_opa,nom_nuc,nom_nuc_cpl,nom_atm,nom_tdetau

c--------------------------------------------------------------------------------
	  
	 REWIND(unit=3)	
	 READ(3,nl_espace, ERR=101, END=101)	 	  
	 WRITE(*,nl_espace) ; WRITE(2,nl_espace)
	 READ(3,nl_temps)   ; WRITE(*,nl_temps) ; WRITE(2,nl_temps)
	 READ(3,nl_chim)    ; WRITE(*,nl_chim)  ; WRITE(2,nl_chim)
	 READ(3,nl_rot)     ; WRITE(*,nl_rot)   ; WRITE(2,nl_rot)
	 READ(3,nl_conv)    ; WRITE(*,nl_conv)  ; WRITE(2,nl_conv)
	 READ(3,nl_etat)    ; WRITE(*,nl_etat)  ; WRITE(2,nl_etat)
	 READ(3,nl_opa)     ; WRITE(*,nl_opa)   ; WRITE(2,nl_opa)
	 READ(3,nl_nuc)     ; WRITE(*,nl_nuc)   ; WRITE(2,nl_nuc)
	 
	 IF(.NOT.rot_solid)THEN      !rotation non solide
	  SELECT CASE(langue)
	  CASE('english')	 	
	   WRITE(*,1002)chain ; WRITE(2,1002)chain
1002	   FORMAT('STOP, non solid rotation is not operating.',
	1   /,'with this input data file : ',a)
	  CASE DEFAULT		 
	   WRITE(*,2)chain ; WRITE(2,2)chain
2	   FORMAT('ARRET, rotation non solide non impl�ment�e.',/,
	1  'avec ce fichier de donn�es : ',a)
	  END SELECT
	  ok=.FALSE.
	  
c �mulation des donn�es	 
	 ELSE
	  IF(f_rad)THEN
	   nom_frad='alecian2'
	  ELSE
	   nom_frad='no_frad'
	  ENDIF
	  IF(ABS(w_rot) > 1.d0)THEN
	   unit='jours'
	  ELSE
	   unit='rad/s'
	  ENDIF
	  
	  INQUIRE(file='physique45',exist=oui)	 
	  IF(oui)THEN
	  OPEN(unit=4,form='formatted',status='old',delim='apostrophe',
	1   file='physique45')
	   READ(4,nl_physique) ; WRITE(*,nl_physique) ; WRITE(2,nl_physique)
	   CLOSE(unit=4)	  
	  ELSE	 
	   nom_chemin='~/SUN_STAR_DATA/'
	   nom_ctes='ctes_94'
	   nom_des='des_m'
	   nom_output='no_output'
	   nom_pertm='pertm_ext'
	   nom_conv='conv_jmj'
	   nom_diffm='diffm_mp'
	   nom_difft='difft_nu'
	   nom_etat='etat_eff'
	   nom_opa='opa_yveline'
	   nom_nuc='ppcno9'
	   nom_nuc_cpl='NACRE'
	   nom_atm='lim_atm'
	   nom_tdetau='hopf'
	  ENDIF
	  nom_diffw='diffw_mpz'	  
	  nom_pertw='pertw_0'
 	  nom_abon='solaire_gn'	   
 	  n_max=2000	  
	  nb_max_modeles=2000
 	  modif_chim=.FALSE.
 	  garde_xish=.FALSE.
	  d_turb=10.d0  
	  w_form=0.d0 
	  p_pertw=0.d0
	  pw_extend=0.d0
	  He_core=-100.d0	  	  
	  ok=.TRUE.
	  RETURN
	 ENDIF		 
	
101	 ok=.FALSE.
	 RETURN	
	 
	END SUBROUTINE lit_nl_45
	
      
	END SUBROUTINE lit_nl
