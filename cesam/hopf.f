
c****************************************************************

	SUBROUTINE hopf(tau,teff,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)

c	routine private du module mod_tdetau

c	loi t(tau) de hopf Mihalas stellar atmospheres (3.16) p. 55 et 72

c	Auteur: P.Morel, Département J.D. Cassini, O.c.A.
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

c----------------------------------------------------------------------

	USE mod_kind
	USE mod_numerique, ONLY : bsp1ddn, no_croiss   
      
	IMPLICIT NONE
      
	REAL (kind=dp), INTENT(in) :: tau, teff
	REAL (kind=dp), INTENT(out) :: df_tau, dro_grav, dro_teff, dtsdg,
	1 dtsdtau, dtsdteff,d2f_tau, f_tau, ro_ext, t 

	INTEGER, PARAMETER :: pm=20, mx=4

	REAL (kind=dp), SAVE, DIMENSION(pm+2*mx) :: taut
	REAL (kind=dp), PARAMETER, DIMENSION(pm) :: taus=(/
	1 0.0d0,0.01d0,0.03d0,0.05d0,0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,
	2 0.8d0,1.00d0,1.50d0,2.00d0,2.5d0,3.0d0,3.5d0,4.0d0,5.0d0,1.0d2 /)
	REAL (kind=dp), SAVE, DIMENSION(pm) :: q1=(/
	1 0.577351d0,0.588236d0,0.601242d0,0.610758d0,0.627919d0,
	2 0.649550d0,0.663365d0,0.673090d0,0.680240d0,0.685801d0,
	3 0.693534d0,0.698540d0,0.705130d0,0.707916d0,0.709191d0,
	4 0.709806d0,0.710120d0,0.710270d0,0.710398d0,0.710446d0 /)
	REAL (kind=dp), SAVE, DIMENSION(1,pm) :: qs
	REAL (kind=dp), DIMENSION(mx,0:mx-1) :: fx
	REAL (kind=dp), PARAMETER :: ro_ext0=3.55d-9
c$$$	REAL (kind=dp), PARAMETER :: ro_ext0=9.99d-12
	REAL (kind=dp), SAVE :: cte1

	INTEGER, SAVE :: l=1, knot

	LOGICAL, SAVE :: init=.TRUE.

c------------------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE. ; rad=.TRUE. ; cte1=(3.d0/4.d0)**0.25d0
	 qs=RESHAPE(q1,(/1,pm/)) ; qs(1,:)=qs(1,:)+taus !tau+q(tau)
	 CALL bsp1ddn(1,qs,taus,taut,pm,mx,knot,.FALSE.,taus(1),l,fx)
	 IF(no_croiss)THEN
          PRINT*,'Arrêt 1 dans hopf' ; STOP
         ENDIF
	 tau_min=1.d-4 
	 WRITE(*,1)tau_min,ro_ext0 ; WRITE(2,1)tau_min,ro_ext0
1	 FORMAT(/'loi t(tau,teff,grav), hopf, purement radiative',/,
	1 'tau_min=',es10.3,', ro_ext=',es10.3,/)     
	ENDIF

c	PRINT*,'hopf ',tau,teff,grav ; PAUSE'hopf1'

	CALL bsp1ddn(1,qs,taus,taut,pm,mx,knot,.TRUE.,tau,l,fx)
	f_tau=fx(1,0) ; df_tau=fx(1,1) ; d2f_tau=fx(1,2)

	t=teff*cte1*f_tau**0.25d0 ; dtsdtau=t/4./f_tau*df_tau
	dtsdteff=t/teff ; dtsdg=0.d0 ; ro_ext=ro_ext0
	dro_grav=0.d0 ; dro_teff=0.d0

c	PRINT*,t,dtsdtau,dtsdteff,dtsdg,cte1,f_tau,df_tau,d2f_tau

	RETURN
	
	END SUBROUTINE hopf
