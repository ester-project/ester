
c*****************************************************************

	SUBROUTINE opa_gong(t,ro,kap,dkapt,dkapro,dkapx)

c	routine private du module mod_opa

c	opacité simplifiée du projet GONG

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c entrées :
c	xchim : composition chimique
c	t : température
c	ro : densité

c sorties :
c	kap : opacité
c	dkapt : d kap / dt
c	dkapro : d kap / dro
c	dkapx : d kap /dx

c---------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) :: t, ro
	REAL (kind=dp), INTENT(out) ::  kap, dkapt, dkapro, dkapx

	REAL (kind=dp), PARAMETER :: ce=1.6236784d-33, kapme=0.407895d0,
	1 kapne=9.28289d0, ci=7.1548412d13, kapmi=0.138316d0,
	2 kapni=-1.97541d0
	REAL (kind=dp) :: dkero, dket, dkiro, dkit, ke, ki

	LOGICAL, SAVE :: init=.TRUE.

c------------------------------------------------------------------

	IF(init)THEN
	 init=.FALSE.
	 WRITE(*,1) ; WRITE(2,1)
1	 FORMAT(' opacité de GONG : formule analytique, Z=0.02',/)
	ENDIF

	ke=ro**kapme*ce*t**kapne  !opacité 3.2
	dkero=ke*kapme/ro ; dket= ke*kapne/t
	ki=ci*ro**kapmi*t**kapni  !3.3
	dkiro=ki*kapmi/ro ; dkit =ki*kapni/t
	kap=ke*ki/(ke+ki) ; dkapx=0.d0      !3.1
	dkapro=(ke**2*dkiro+ki**2*dkero)/(ki+ke)**2
	dkapt= (ke**2*dkit +ki**2*dket )/(ki+ke)**2

	RETURN
	
	END SUBROUTINE opa_gong
