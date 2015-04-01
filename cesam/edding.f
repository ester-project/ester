
c****************************************************************

	SUBROUTINE edding(tau,teff,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)

c	routine private du module mod_atm

c	Loi T(tau) Eddington

c	Auteur: P.Morel, Département J.D. cassini, O.c.A.,
c	CESAM2k

c entrées :
c	tau : profondeur optique Rosseland
c	teff : température effective
c	grav : gravité

c sorties :
c	t : température
c	ro_ext: densité externe
c	dtsd* : dérivées t/ tau, teff, grav
c	dro_** : dérivées ro_ext/ teff, grav
c	tau_ext: profondeur optique externe
c	f_tau,df_tau,d2_gtau : g(tau) et dérivées

c-----------------------------------------------------------------------

	USE mod_donnees, ONLY : tau_max
	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) :: tau, teff
	REAL (kind=dp), INTENT(out) :: t, dtsdtau, dtsdteff, dtsdg,
	1 ro_ext, dro_grav, dro_teff, f_tau, df_tau, d2f_tau

	REAL (kind=dp), PARAMETER :: ro_ext0=3.55d-9
	REAL (kind=dp), SAVE :: cte1

	LOGICAL, SAVE :: init=.TRUE.

c-------------------------------------------------------------------

	IF (init) then
	 init=.FALSE. ; rad=.TRUE. ; cte1=(3.d0/4.d0)**0.25d0    
	 WRITE(*,1) ; WRITE(2,1)
1	 FORMAT(/,'loi t(tau,teff,grav) Eddington, loi purement radiative',/,
	1 'tau_min=1.d-4, ro_ext=3.55d-9',/)
	 IF (tau_max < 0.8d0) then
	  WRITE(*,2)tau_max ; WRITE(2,2)tau_max
2	  FORMAT('avec Eddington tau_max >= 0.8, tau_max initial=',es10.3)
	  tau_max=0.8d0
	 ENDIF
	 tau_min=1.d-4
	ENDIF

	f_tau=tau+2.d0/3.d0 ; df_tau=1.d0 ; d2f_tau=0.d0

	t=teff*cte1*f_tau**0.25d0 ;dtsdtau=t/4.d0/f_tau*df_tau
	dtsdteff=t/teff
	dtsdg=0.d0 ; dro_grav=0.d0 ; dro_teff=0.d0 ; ro_ext=ro_ext0

c	PRINT*,t,dtsdtau,dtsdteff,dtsdg,cte1

	RETURN

	END SUBROUTINE edding
