
c*******************************************************************

	SUBROUTINE roger10a(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)

c	routine private du module mod_atm
    
c	interface pour la loi T(tau) de R.Cayrel avec [Fe/H]=-1.0 enhanced

c	Auteur: P. Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

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

	LOGICAL, SAVE :: init=.TRUE.

c-----------------------------------------------------------------

	IF(init)THEN
	 init=.FALSE.
	 SELECT CASE(langue)	 
	 CASE('english')
	  WRITE(*,1001) ; WRITE(2,1001)
1001	  FORMAT(/,'T(tau,teff,grav) law, roger10a, no purely radiative',/,
	1 '[Fe/H]=-1. alpha enhanced, interpolations in Teff and Log g',/,
	2 'artificial extend of the data to Teff = 3500K',/,
	3 '3500K < Teff < 7000K, 1 < Log10(g) < 5, tau_min=1.d-4',/)	  
	 CASE DEFAULT	  
 	  WRITE(*,1) ; WRITE(2,1)
1	  FORMAT(/,'loi T(tau,teff,grav) roger10a, non purement radiative',
	1 /,'[Fe/H]=-1. alpha enhanced interpolations en Teff et Log g',/,
	2 'table étendue artificiellement a Teff = 3500K',/,
	3 '3500K < Teff < 7000K, 1 < Log10(g) < 5, tau_min=1.d-4',/)
	 END SELECT 
	 nom=TRIM(nom_chemin)//'fesh-10a.data' ; rad=.FALSE.
	ENDIF

	logg=LOG10(grav) ; tef=teff     !peut être changé dans trho
	CALL trho(nom,tau,tef,logg,t,dtsdtau,dtsdteff,
	1 dtsdg,ro_ext,dro_teff,dro_grav,f_tau,df_tau,d2f_tau)

	RETURN

	END SUBROUTINE roger10a
