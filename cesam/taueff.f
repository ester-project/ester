
c**********************************************************************

	SUBROUTINE taueff(teff,grav,tau)

c	subroutine private du module mod_tdetau
c	on détermine la valeur de tau tel que Teff=T(tau,teff,grav)
c	utilisation de la fonction générique tdetau

c entrées
c	teff : température effective
c	grav : gravité
c	tdetau : loi T(tau)

c sorties
c	tau : tel que t=T(tau,teff,grav)

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c----------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) :: grav, teff
	REAL (kind=dp), INTENT(out) :: tau
	REAL (kind=dp) :: epsi, dtsdtau, t, dtsdteff, dtsdg,
	1 ro_ext, dro_grav ,dro_teff, g_tau, dg_tau, d2g_tau

	INTEGER :: tour=0

c------------------------------------------------------------------------

2000	FORMAT(8es10.3)

c	WRITE(*,2000)teff,grav
	tau=2.d0/3.d0 ; epsi=1.d3 ; tour=0
	
	B1: DO
	 tour=tour+1
	 IF(tour > 30)THEN
	  WRITE(*,*)'taueff, non CV ou tau erronné: on fixe taueff=2/3'
	  tau=2.d0/3.d0 ; EXIT B1
	 ELSEIF(ABS(epsi) < 1.d-4)THEN
	  EXIT B1
	 ELSE
	  CALL tdetau(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,g_tau,dg_tau,d2g_tau)
	  epsi=(t-teff)/dtsdtau ; tau=tau-epsi	  
	 ENDIF
	ENDDO B1
	
	RETURN

	END SUBROUTINE taueff
