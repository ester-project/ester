
c*******************************************************************

	SUBROUTINE rogerYL(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)

c	routine private du module mod_atm
    
c	interface pour la loi T(tau) de R. cayrel avec [Fe/H] donne
c	appel du fichier 'mod_atm'

c	Auteur: P. Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c       modifie par Y. Lebreton (CESAM4 puis 2k)    !YLD

c entrées :
c	tau : profondeur optique Rosseland
c	teff : température effective
c	grav : gravité

c sorties :
c	t : température
c	dtsd* : dérivées t/ tau, teff, grav
c	dtsd* : dérivées t/ tau, teff, grav
c	dro_** : dérivées ro_ext/ teff, grav
c	f_tau, df_tau, df_tau2 : f, d f / d tau, d2 f / d2 tau

c------------------------------------------------------------------------

	USE mod_donnees, ONLY : langue
	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) :: tau, teff, grav
	REAL (kind=dp), INTENT(out) :: t, dtsdtau, dtsdteff, dtsdg,
     1 ro_ext, dro_grav, dro_teff, f_tau, df_tau, d2f_tau

	REAL (kind=dp) :: logg, tef

	CHARACTER (len=80) :: nom
	CHARACTER (len=80) :: nom_chemin = "/data1/sdeheuve/local/src/cesam2k_v1.1.8_ESTA/SUN_STAR_DATA/"

	LOGICAL, SAVE :: init=.TRUE., kuredd     !YLD

	LOGICAL, SAVE :: initedd                  !YLD
c-----------------------------------------------------------------

	IF(init)THEN
	 init=.FALSE.
	 SELECT CASE(langue)	 
	 CASE('english')
	  WRITE(*,1001) ; WRITE(2,1001)
1001	  FORMAT(/,'T(tau,teff,grav) law, rogerYL, no purely radiative',/,  !YLD
	1 '[Fe/H]=+0.2 interpolations in Teff and Log g',/,
	2 'artificial extend of the data to Teff = 3500K',/,
	3 '3500K < Teff < 7000K, 1 < Log10(g) < 5, tau_min=1.d-4',/)	  
	 CASE DEFAULT	  
	  WRITE(*,1) ; WRITE(2,1)
1	  FORMAT(/,'loi T(tau,teff,grav), roger02, non purement radiative',
	1 /,'[Fe/H]=+0.2 interpolations en Teff et gravité',/,
	2 'table étendue artificiellement à Teff = 3500K',/,
	3 '3500K < Teff < 7000K, 1 < Log10(g) < 5, tau_min=1.d-4',/)
	 END SELECT
	 nom=TRIM(nom_chemin)//'ATM/mod_atm' ; rad=.FALSE.      !YLD
	ENDIF

	logg=LOG10(grav) ; tef=teff     !peut être changé dans trho
	
c       YLD : ajout de la variable kuredd dans trho
	CALL trhoYL(nom,tau,tef,logg,t,dtsdtau,dtsdteff,            !YLD
	1 dtsdg,ro_ext,dro_teff,dro_grav,f_tau,df_tau,d2f_tau,kuredd)  !YLD

	IF(kuredd.eqv.(.TRUE.)) THEN                 !YLD
c        init=.TRUE.                                      !YLD
        initedd=.TRUE.                               !YLD
	CALL edding(tau,teff,t,dtsdtau,dtsdteff,dtsdg,      !YLD
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)    !YLD
	ENDIF                                               !YLD

	RETURN

	END SUBROUTINE rogerYL
