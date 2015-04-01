		
c*************************************************************************

	SUBROUTINE mu_mol(dxchim,hp,nu,r,ro,t,xchim,dlnmu,dlnmu_x,grad_mu,
	1 dgrad_mux,mu,dmu_x,nel,Zbar)
	
c subroutine PRIVATE du module mod_etat

c calcul de diverses quantités reliées au poids moléculaire moyen

c entrées :
c	dxchim : d xchim / d nu
c	hp : échelle hauteur de pression
c	nu : m^2/3
c	r : r/Rsol
c	ro : densité
c	t : température
c	xchim : composition chimique par mole

c sorties :
c	dlnmu, dlnmu_x : d ln mu / d nu, dérivée / X mole
c	grad_mu, dgrad_mux : dln mu /d ln P, dérivée / X mole
c	mu, dmu_x : poids moléculaire moyen, dérivée / X mole
c	nel : nombre d'électrons libres
c	Zbar : charge moyenne	

c Auteur: P.Morel, Département Cassiopée, O.C.A., Observatoire de Nice
c CESAM2k

c-----------------------------------------------------------------
	
	USE mod_donnees, ONLY : amu, ihe4, msol, mu_saha, nchim, nucleo, pi,
	1 rsol, zi, z0	
	USE mod_kind
	Use mod_variables, ONLY : chim_gram
			
	IMPLICIT NONE
		
	REAL (kind=dp), INTENT(in), DIMENSION(nchim) :: dxchim, xchim	
	REAL (kind=dp), INTENT(in) :: hp, nu, r, ro, t

	REAL (kind=dp), INTENT(out), DIMENSION(nchim) :: dlnmu_x, dgrad_mux,
	1 dmu_x 		
	REAL (kind=dp), INTENT(out) :: dlnmu, grad_mu, mu, nel, Zbar
		
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: ioni
	REAL (kind=dp), DIMENSION(nchim) :: z_bar

	REAL (kind=dp), SAVE ::	cte1
	REAL (kind=dp) :: eta, Y	
		
	LOGICAL, SAVE :: init=.TRUE.	
	
c--------------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN	!initialisations
	 init=.FALSE.
	 cte1=8.d0*pi*rsol**2/3.d0/msol
	 IF(mu_saha)THEN
	  ALLOCATE(ioni(0:NINT(MAXVAL(zi)),nchim))
	 ELSE
	  z_bar=zi
	 ENDIF	  	 
	ENDIF

c cas de PP1 
	IF(nchim == 1)THEN
	 Y=1.d0-xchim(1)*nucleo(1)-z0
	 mu=16.d0/(23.d0*xchim(1)*nucleo(1)+3.d0*Y+9.d0)
	 dmu_x=-mu**2*23.d0/16.d0*nucleo(1)
	 dlnmu=dmu_x(1)/mu*dxchim(1) ; dlnmu_x(1)=dlnmu
	 nel=DOT_PRODUCT(zi,xchim)*ro/amu
	 
c autres réseaux	 
	ELSE	 
c ionisation partielle	 
	 IF(mu_saha)CALL saha(xchim,t,ro,ioni,z_bar,nel,eta)
	 mu=1.d0/DOT_PRODUCT((1.d0+z_bar),xchim)
	 dmu_x=-mu**2*(1.d0+z_bar)
	 dlnmu=-mu*DOT_PRODUCT((1.d0+z_bar),dxchim)
	 dlnmu_x=-dmu_x*DOT_PRODUCT((1.d0+z_bar),dxchim)	 
	ENDIF
	Zbar=DOT_PRODUCT(z_bar,xchim)/SUM(xchim)
	
c grad_mu
	IF(nu <= 0.d0)THEN
	 grad_mu=0.d0	
	ELSE
	
c signe - devant dlnmu			
	 grad_mu=-cte1*r**2*ro*dlnmu*hp/SQRT(ABS(nu))
	 dgrad_mux=-cte1*r**2*ro*dlnmu_x*hp/SQRT(ABS(nu))	 
	ENDIF
		
c grad_mu > 0
	IF(grad_mu <= 0.d0)THEN
	 grad_mu=0.d0 ; dgrad_mux=0.d0
	ENDIF
			
	RETURN
	
	END SUBROUTINE mu_mol	
