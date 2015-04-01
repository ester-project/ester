
c**************************************************************************

	SUBROUTINE tdetau(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)
C
c routine générique pour le calcul des lois t(tau)    
c il y a des appels differents suivant suivant nom_tdetau
c Appel à k5777 en cas de sorie de table pour roger et marcs
C
c routine public du module mod_atm 
C
c Auteur: P.Morel, Département J.D. Cassini, O.C.A.,CESAM2k
C
c entrées :
c	tau : profondeur optique Rosseland
c	teff : température effective
c	grav : gravité
c
c sorties :
c	t : température
c	dtsd* : dérivées t/ tau, teff, grav
c	dtsd* : dérivées t/ tau, teff, grav
c	dro_** : dérivées ro_ext/ teff, grav
c	f_tau, df_tau, df_tau2 : f, d f / d tau, d2 f / d2 tau
C
c--------------------------------------------------------------------
C
 	USE mod_donnees, only :  langue, nom_tdetau
	USE mod_kind
C        
	IMPLICIT NONE
C        
	REAL (kind=dp), INTENT(in) :: tau, teff, grav
	REAL (kind=dp), INTENT(out) :: t, dtsdtau, dtsdteff, dtsdg,
	1 ro_ext, dro_grav, dro_teff, f_tau, df_tau, d2f_tau
	
	LOGICAL, SAVE :: init=.TRUE.
	
	Character(Len=5) :: nom_reduit
C
c-----------------------------------------------------------------------
C
      nom_reduit = nom_tdetau(1:5)
            
	SELECT CASE(nom_reduit)

	CASE ('eddin')
	 CALL edding(tau,teff,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)

	CASE ('hopf ')
	 CALL hopf(tau,teff,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)

	CASE ('hsra ')
	 CALL hsra(tau,teff,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)	

	CASE ('k5750')
	 CALL k5750(tau,teff,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)

 	CASE ('k5777')
	 CALL k5777(tau,teff,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)	

	CASE ('MARCS')
	 IF(init)THEN
	  init=.FALSE.
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1001)
1001	   FORMAT('The data for the MARCS''s T(tau) laws are available',/,
	1  'from B.Pichon, OCA')	 
	  CASE DEFAULT
	   WRITE(*,1)
1	   FORMAT('Les données des lois T(tau) de MARCS sont disponibles',/,
	1  'sur demande auprès de B.Pichon, OCA')
	  END SELECT
	 ENDIF
	 CALL marcs(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)
	
	CASE ('roger')
	 CALL roger(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)

	CASE DEFAULT
	 PRINT*,'routine de loi T(tau) inconnue: ',nom_tdetau
	 PRINT*,'routines connues: edding, hopf, hsra, k5750, k5777'
	 PRINT*,'roger****, MARCS****'
	 PRINT*,'arrêt' ; STOP
	END SELECT
	
	RETURN

	END SUBROUTINE tdetau
