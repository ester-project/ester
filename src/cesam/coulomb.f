
c*********************************************************************
	 
	SUBROUTINE coulomb(zi,zj,thetae,ro,x,t,lnlambij,
	1 lnlambijr,lnlambijx,cij,cijr,cijx)

c	routine private du module mod_evol

c	calcul du logarithme de Coulomb
c	d'après Michaud-Proffit, IAU 137, p.250

c	Auteur: J. Mathias, DASGAL, pour H et He
c	adaptation: P. Morel, Departement J.D. Cassini, O.C.A.
c	cesam2k

c entrées :
c	zi,zj :charges éléments i et j
c	thetae = 1, pour mixture H, He
c	ro : densité
c	x : abondance H1 ou He4
c	t : température

c sorties :
c	lnlambij,cij : coeff de Coulomb
c	lnlambijr,cijr : derivees coeff de Coulomb / ro
c	lnlambijx,cijx : derivees coeff de Coulomb / x

c-----------------------------------------------------------------------

	USE mod_kind
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) :: zi, zj, thetae, ro, x, t	
	REAL (kind=dp), INTENT(out) :: lnlambij, cij, lnlambijx, cijx,
	1 lnlambijr, cijr
	
c------------------------------------------------------------------------
	
2000	FORMAT(8es10.3)

	lnlambij=-19.26d0-log(zi*zj)-0.5d0*log(ro) 
	1 -0.5d0*log(1.d0+((x+1.d0)/2.d0)*thetae)+1.5d0*log(t)
	lnlambijr=-0.5d0/ro		!derivée / ro	
	lnlambijx=-0.5d0/(1.d0+((x+1.d0)/2.d0)*thetae)*thetae/2.d0 !der./x

	cij =log(exp(1.2d0*lnlambij)+1.d0)/1.2d0
	cijx=1.2d0*exp(1.2d0*lnlambij)*lnlambijx
	cijx=cijx/1.2d0/(exp(1.2d0*lnlambij)+1.d0)	!dérivee / x
	cijr=1.2d0*exp(1.2d0*lnlambij)*lnlambijr
	cijr=cijr/1.2d0/(exp(1.2d0*lnlambij)+1.d0)	!dérivee / ro
	
c	PRINT*,cij,lnlambij ; WRITE(*,2000)cij,lnlambij,zi,zj,thetae,ro,t,x

	RETURN

	END SUBROUTINE coulomb
