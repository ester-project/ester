	
c**********************************************************************

	SUBROUTINE tdetau(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)

c	routine générique pour le calcul des lois t(tau)    
c	il y a des appels differents suivant suivant nom_tdetau

c	routine public du module mod_atm 

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
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
      
c--------------------------------------------------------------------

 	USE mod_donnees, only :  nom_tdetau
	USE mod_kind
        
	IMPLICIT NONE
        
	REAL (kind=dp), INTENT(in) :: tau, teff, grav
	REAL (kind=dp), INTENT(out) :: t, dtsdtau, dtsdteff, dtsdg,
	1 ro_ext, dro_grav, dro_teff, f_tau, df_tau, d2f_tau
      
c-----------------------------------------------------------------------
      
	SELECT CASE(nom_tdetau)
	CASE ('edding')
	 CALL edding(tau,teff,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)
	CASE ('hopf')
	CALL hopf(tau,teff,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)
	CASE ('kurucz')                                       !YLD
	CALL kurucz(tau,teff,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)
	CASE ('k5750')
	CALL k5750(tau,teff,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)
	CASE ('k5777')
	 CALL k5777(tau,teff,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)
	CASE ('rogerYL')                                      !YLD
	CALL rogerYL(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)
	CASE ('roger00')
	CALL roger00(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)
	CASE ('roger02')
	CALL roger02(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)  
	CASE ('roger05')
	CALL roger05(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)
	CASE ('roger10a')
	CALL roger10a(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)  
	CASE DEFAULT
	PRINT*,'routine de loi T(tau) inconnue: ',nom_tdetau
	PRINT*,'routines connues: edding, hopf, k5750, k5777'
	PRINT*,'roger00, roger02, roger05, roger10a'
	PRINT*,'arrêt' ; STOP
	END SELECT
	
	RETURN

	END SUBROUTINE tdetau
