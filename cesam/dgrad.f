
c****************************************************************************

	REAL (kind=dp) FUNCTION dgrad(pt,p,t,dlpp,xchim,m,l,r,dxchim,w)

c	routine public du module mod_static

c	calcul de la différence des gradients pour le critère de convection

c	Auteur: P.Morel + N. Audard, Département J.D. Cassini, O.C.A.
c	CESAM2k

c	Modifs
c	09 10 96 : introduction de w
c	09 11 99 : correction critère de Ledoux
c	19 11 99 : suppression de nh1, nhe1, nhe2, lamb
c	30 07 00 : introduction F95

c entrées :
c	pt : pression totale
c	p : pression gazeuse cgs
c	t : température K
c	xchim : composition chimique ATTENTION en 1/mole : *ah pour avoir X
c	dxchim : d xchim/d m pour critère de Ledoux uniquement
c	m : masse/msol
c	l : luminosite/lsol
c	r : rayon / rsol
c	mstar: masse au temps du calcul, avec perte de masse
c	w : rotation
c	dlpp : d ln Pgaz / d ln Ptot

c sortie:
c	dgrad=grad_rad-grad_ad > 0 dans ZR

c----------------------------------------------------------------

	USE mod_donnees, ONLY : aradia, clight, cpturb, g, ihe4,
	1 ledoux, lsol, msol, nchim, pi, rsol, t_inf
	USE mod_etat, ONLY : etat
	USE mod_kind
	USE mod_opa, ONLY : opa
	USE mod_nuc, ONLY : nuc
	USE mod_variables, ONLY : chim_gram
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: dxchim, xchim
	REAL (kind=dp), INTENT(in) :: dlpp, l, m, p, pt, t, r, w

	REAL (kind=dp), SAVE, ALLOCATABLE, DIMENSION(:,:) :: jac	
	REAL (kind=dp), SAVE, ALLOCATABLE, DIMENSION(:) :: depsx, dcomp,
	1 dxchimm, xchi, xchimm
	REAL (kind=dp), SAVE, DIMENSION(5) :: epsilon
	REAL (kind=dp), SAVE :: cte1, cte2, cte7, cte8,
	1 cte9, cte13
	REAL (kind=dp) :: alfa, beta, be7, b8, cp, dcpp, dcpt, dcpx, delta,
	1 deltap, deltat, deltax, depst, depsro, drop, drot, drox, dup,
	2 dut, dux, f17, gamma1, gradad, dgradadp, dgradadt, dgradadx,
	3 gradrad, gravite, hh, hp, kap, dkapro, dkapt, dkapx,
	4 krad, ldx, n13, o15, ro, u

	LOGICAL, SAVE :: init=.TRUE.

c--------------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN	!initialisations
	 init=.FALSE.
	 cte1=4.d0/3.d0*aradia*clight ; cte13=g*msol/rsol/rsol	 	 
	 cte2=2.d0/3.d0*rsol ; cte7=4.d0*pi*rsol**2/msol
	 cte8=lsol/4.d0/pi/rsol/rsol	!de 5.9
	 cte9=3.d0/16.d0/pi/aradia/clight/g
	 
	 ALLOCATE(xchimm(nchim),xchi(nchim),dxchimm(nchim),
	1 depsx(nchim),jac(nchim,nchim),dcomp(nchim))
	 
	ENDIF		!initialisation
	
c	composition chimique /gr
	
	xchimm=abs(xchim) ; dxchimm=dxchim ; xchi=xchimm		 
	
	CALL chim_gram(xchimm,dxchimm)

	CALL etat(p,t,xchimm,.FALSE.,
	1 ro,drop,drot,drox,u,dup,dut,dux,
	2 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	IF(cpturb < 0.d0)gradad=gradad*dlpp

	CALL opa(xchimm,t,ro,kap,dkapt,dkapro,dkapx)
	krad=cte1/kap/ro*t**3		!5.1 conductivite radiative
	
	IF(m*l*r /= 0.d0)THEN		!gradient radiatif
	 gravite=cte13*m/r**2-cte2*w**2*r !gravité effective avec rotation
	 gradrad=cte8*l*p/gravite/ro/r**2/krad/t	!5.9
	 hp=pt/gravite/ro		!echelle de hauteur de pression
	ELSE		!au centre
	 IF(t > t_inf)THEN
	  CALL nuc(t,ro,xchi,dcomp,jac,.FALSE.,3,
	1 epsilon,depst,depsro,depsx,hh,be7,b8,n13,o15,f17)
	 ELSE
	  epsilon(1)=0.d0	!total
	 ENDIF
	 gradrad=cte9*kap*epsilon(1)*p/t**4	!au centre l/m ~ epsilon
	ENDIF
	dgrad=gradrad-gradad*dlpp	!critere de Schwarzschild

	IF(ledoux .AND. ihe4 > 1 .AND. m > 0.d0)THEN
		
c	 formulation si ro(P,T,X) eg. OPAL, MHD
	
	 ldx=-cte7*drox*dxchimm(1)*r**2*hp/delta

c	 formulation 3-30 Cox & Guili p. 276 p = ro R T / mu + a/3 T^4
c	 mu = 4 / (2 + ^ X + Y) ie. tot. ionisé
 
c	 ldx=cte7*r**2*hp*ro/(2.d0+6.*xchimm(1)+xchimm(ihe4))*
c	1	beta/(4.d0-3.d0*beta)*(6.d0*dxchimm(1)+dxchimm(ihe4))	
	
c	 IF(dgrad > 0.d0 .AND. ldx >= dgrad)THEN 
c	  WRITE(*,2000)dxchimm(1),gradrad,gradad,ldx,dgrad,r
c	  PAUSE'dans dgrad'
c	 ENDIF
	 dgrad=dgrad-ldx
	ENDIF
	
c	IF(m*l*r /= 0.d0)WRITE(*,2000)l,p,gravite,ro,r,krad,t,
c	1 gradrad,gradad,dgrad
	
	RETURN

	END FUNCTION dgrad
