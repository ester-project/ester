
c****************************************************************

	SUBROUTINE hsra(tau,teff,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)

c routine private du module mod_atm

c loi t(tau) HSRA 5780, 4.44, He/h by number .1 .Sol. Phys. 18,347
c mise sous la forme t=Teff(3/4(tau+q(tau)))**1/4
c la gravité est fixée
c le point (100,9756) est ajouté pour cohérence
c en tant que loi non purement radiative on ne dépassera pas tau=20

c Auteur: P. Morel, Département J.D. Cassini, O.C.A., CESAM4 --> CESAM2k

c entrées :
c	tau : profondeur optique Rosseland
c	teff : temperature effective

c sorties :
c	t : température
c	dtsd* : dérivées t/ tau, teff, grav
c	dtsd* : dérivées t/ tau, teff, grav
c	dro_** : dérivées ro_ext/ teff, grav
c	f_tau,df_tau,d2_gtau : g(tau) et dérivées

c------------------------------------------------------------------

	USE mod_kind
	USE mod_numerique, ONLY : bsp1ddn, no_croiss
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) :: tau, teff
	REAL (kind=dp), INTENT(out) :: df_tau, dro_grav, dro_teff, dtsdg,
	1 dtsdtau, dtsdteff, f_tau, d2f_tau, ro_ext, t

	INTEGER, PARAMETER :: mx=4, nx=56
	REAL (kind=dp), PARAMETER, DIMENSION(nx) :: taus=(/
	1 1.000d-4, 1.259d-4, 1.585d-4, 1.995d-4, 2.512d-4, 3.162d-4, 3.981d-4,
	2 5.012d-4, 6.310d-4, 7.943d-4, 1.000d-3, 1.259d-3, 1.585d-3, 1.995d-3, 
	3 2.512d-3, 3.162d-3, 3.981d-3, 5.012d-3, 6.310d-3, 7.943d-3, 1.000d-2, 
	4 1.259d-2, 1.585d-2, 1.995d-2, 2.512d-2, 3.162d-2, 3.981d-2, 5.012d-2,
	5 6.310d-2, 7.943d-2, 1.000d-1, 1.259d-1, 1.585d-1, 1.995d-1, 2.512d-1, 
	6 3.162d-1, 3.981d-1, 5.012d-1, 6.310d-1, 7.943d-1, 1.000d+0, 1.259d+0, 
	7 1.585d+0, 1.995d+0, 2.512d+0, 3.162d+0, 3.981d+0, 5.012d+0, 6.310d+0,
	8 7.943d+0, 1.000d+1, 1.259d+1, 1.585d+1, 1.995d+1, 2.512d+1, 1.000D+2/)
		
	REAL (kind=dp), PARAMETER, DIMENSION(nx) :: t1=(/
	1 4.170d3, 4.175d3, 4.190d3, 4.205d3, 4.225d3, 4.250d3, 4.280d3,
	2 4.305d3, 4.330d3, 4.355d3, 4.380d3, 4.405d3, 4.430d3, 4.460d3,
	3 4.490d3, 4.525d3, 4.550d3, 4.575d3, 4.600d3, 4.630d3, 4.660d3,
	4 4.690d3, 4.720d3, 4.750d3, 4.790d3, 4.840d3, 4.895d3, 4.950d3,
	5 5.010d3, 5.080d3, 5.160d3, 5.240d3, 5.330d3, 5.430d3, 5.540d3,
	6 5.650d3, 5.765d3, 5.890d3, 6.035d3, 6.200d3, 6.390d3, 6.610d3,
	7 6.860d3, 7.140d3, 7.440d3, 7.750d3, 8.030d3, 8.290d3, 8.520d3,
	8 8.710d3, 8.880d3, 9.050d3, 9.220d3, 9.390d3, 9.560d3, 9.756d3 /)

	REAL (kind=dp), SAVE, DIMENSION(nx+2*mx) :: taut
	REAL (kind=dp), SAVE, DIMENSION(1,nx) :: ts
	REAL (kind=dp), DIMENSION(mx,0:mx-1) :: fx
	REAL (kind=dp), PARAMETER :: teffc=5750.d0, ro_ext0=3.55d-9
	REAL (kind=dp), SAVE :: cte1
	
	INTEGER, SAVE :: knot, l=1
	
	LOGICAL, SAVE :: init=.TRUE.

c-----------------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE. ; cte1=(3.d0/4.d0)**0.25d0 ; ts=RESHAPE(t1,(/1,nx/))
	 ts=(ts/teffc)**4*4.d0/3.d0
	 CALL bsp1ddn(1,ts,taus,taut,nx,mx,knot,.FALSE.,taus(1),l,fx)
         IF(no_croiss)STOP 'Arrêt 1 dans hsra'
	 rad=.FALSE. ; tau_min=1.d-4
	 WRITE(*,1)tau_min,ro_ext0 ; WRITE(2,1)tau_min,ro_ext0
1	 FORMAT(/,'loi t(tau,teff), HSRA, non purement radiative',/,
	1 'tau_min=',es10.3,' ro_ext=',es10.3,/)
	ENDIF

	CALL bsp1ddn(1,ts,taus,taut,nx,mx,knot,.TRUE.,tau,l,fx)
	f_tau=fx(1,0) ; df_tau=fx(1,1) ;d2f_tau=fx(1,2)
	t=teff*cte1*f_tau**0.25d0 ; dtsdtau=t/4.d0/f_tau*df_tau
	dtsdteff=t/teff
	dtsdg=0.d0 ; ro_ext=ro_ext0 ; dro_grav=0.d0 ; dro_teff=0.d0

	RETURN

	END SUBROUTINE hsra
