
c******************************************************************

	SUBROUTINE lim_tau1 (l,r,xchim,pt,dptl,dptr,t,dtl,dtr,m,dml,
	1 dmr,p,dpl,dpr,teff)

c	routine interne du module lim_atm

c	calcul de la condition limite en tau=1
c	relations Pext(l,r), Text(l,r), Mext(l,r) et dérivées

c	modifs :
c	19 11 99 : suppression de nh1, nhe1, nhe2, lamb
c	16 08 01 : F95

c	Auteur: P.Morel, Departement J.D. Cassini, O.C.A.
c	CESAM2k

c entrées
c	list=.TRUE. : on calcule p_atm,t_atm,r_atm,tau,m_atm	
c	l :luminosité
c	r :rayon
c	xchim : composition chimique par gramme
c	mstar: masse avec perte de masse
c	wrot : vitesse angulaire solide

c sorties
c	pt : Ptot
c	p : Pgaz
c	t : température
c	m : masse
c	dptl,dptr,dtl,dtr,dml,dmr,dpl,dpr : dérivées pt,t,m p / l,r
c	teff : température effective
c	rstar : rayon externe
c	pt_atm,t_atm,r_atm,tau,m_atm,p_atm, : sans objet
c	n_atm : mis a 1

c------------------------------------------------------------------------

	USE mod_donnees, ONLY : aradia, clight, g, Krot, langue, lsol,
	2 msol, mtot, n_atm, pi, pturb, rsol
	USE mod_etat, ONLY : etat
	USE mod_kind
	USE mod_opa, ONLY : opa
	USE mod_variables, ONLY : dim_ch, mstar, rota, rstar, wrot	
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim
	REAL (kind=dp), INTENT(in) :: l, r
	REAL (kind=dp), INTENT(out) :: pt, dptl, dptr, t, dtl, dtr,
	1 m, dml, dmr, p, dpl, dpr, teff

	REAL (kind=dp), PARAMETER :: tau_b=1.d0, preci=1.d-8
	
	REAL (kind=dp), SAVE :: pprec=1.d5, cte1, cte2, cte3
	REAL (kind=dp) :: w, ro, drop, drot, drox, u, dup, dut, dux,
	1 delta, deltap, deltat, deltax, cp, dcpp, dcpt, dcpx,
	2 gradad, dgradadp, dgradadt, dgradadx, alfa, beta, gamma1,
	3 kap, dkapt, dkapro, dkapx, dkapp, corr

	INTEGER :: ntour

	LOGICAL, SAVE :: init=.TRUE.
	
c----------------------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE.
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(*,1002)tau_b ; WRITE(2,2)tau_b
1002	  FORMAT(/,'Simplified external limit, ',
	1 'dP / dtau = Pext = gravity / kappa',/,
	2 'at tau=',es10.3,' : m=Mstar, r=Rstar, t=Teff')	 
	 CASE DEFAULT	 
	  WRITE(*,2)tau_b ; WRITE(2,2)tau_b
2	  FORMAT(/,'Limite externe sans calcul d''atmosphère ',
	1 'd P / d tau = Pext = gravite / kappa',/,
	2 'en tau=',es10.3,' : m=Mstar, r=Rstar, t=Teff')
	 END SELECT		
	 IF(pturb)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1003) ; WRITE(2,1003)
1003	   FORMAT(/,'WARNING, calculations with the turbulent pressure',/,
	1  'needs a full atmosphere, difficulties may occur')	 	 
	  CASE DEFAULT	 
	   WRITE(*,3) ; WRITE(2,3)
3	   FORMAT(/,'ATTENTION, la pression turbulente nécessite une',/,
	1 'atmosphère, risques de difficultes')
	  END SELECT		
	 ENDIF
	 cte1=lsol/rsol/rsol/aradia/clight*4.d0/pi/4.d0	 
	 cte2=g*msol/rsol/rsol*tau_b ; cte3=2.d0/3.d0*rsol*tau_b
	 
c	 allocations des tableaux pt_atm etc..
	 
	 n_atm=0
	 ALLOCATE(bp_atm(0,0),bp_atm_t(0,0),pt_atm(0),t_atm(0),
	1  m_atm(0),tau(0),p_atm(0),r_atm(0),dlpp_atm(0),x_atm(0),x_atm_t(0),
	2  xt_atm(0),xt_atm_t(0))	
	ENDIF

c	résolution de kap p r**2 = G m par itération Newton-Raphson

	t=(cte1*l/r**2)**(1.d0/4.d0)

c	WRITE(*,*)'dans lim_tau1 p,corr,t,r,l,mtot,xchim(1),kap,dkapp'
c	WRITE(*,2000)p,corr,t,r,l,mtot,xchim(1),kap,dkapp

	p=pprec		!valeur provisoire de pext
	
c	Omega

	SELECT CASE(Krot)
	CASE(0,1,2)
	 w=wrot
	CASE(3)
	 w=rota(1,dim_ch)	 
	END SELECT
	w=w**2 

c	initialisations
	
	ntour=0
	B1: DO
	 ntour=ntour+1
	 IF(ntour > 30)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1004)p,corr,t,r,l,mtot,xchim(1),kap,dkapp
	   WRITE(2,1004)p,corr,t,r,l,mtot,xchim(1),kap,dkapp	   
1004	   FORMAT('STOP, no convergence for pext in lim_tau1',/,
	1  'p=',es10.3,', corr=',es10.3,', t=',es10.3,', r=',es10.3,/,
	2  'l=',es10.3,', mtot=',es10.3,', xchim(1)=',es10.3,
	3  ', kap=',es10.3,', dkapp=',es10.3) 
	  CASE DEFAULT	  		 
	   WRITE(*,4)p,corr,t,r,l,mtot,xchim(1),kap,dkapp
	   WRITE(2,4)p,corr,t,r,l,mtot,xchim(1),kap,dkapp
4	   FORMAT('ARRET, pas de convergence pour pext dans lim_tau1',/,
	1  'p=',es10.3,', corr=',es10.3,', t=',es10.3,', r=',es10.3,/,
	2  'l=',es10.3,', mtot=',es10.3,', xchim(1)=',es10.3,
	3  ', kap=',es10.3,', dkapp=',es10.3)	 
	  END SELECT
	  STOP
	 ENDIF

	 CALL etat(p,t,xchim,.TRUE.,
	1 ro,drop,drot,drox,u,dup,dut,dux,
	2 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	
	 CALL opa(xchim,t,ro,kap,dkapt,dkapro,dkapx)
	 dkapp=dkapro*drop

	 corr=(kap*p-cte2*mstar/r**2+cte3*w*r)/(dkapp*p+kap)
	 p=p-corr
c	 WRITE(*,2000)p,corr,t,r,l,mtot,xchim(1),kap,dkapp
	 IF(ABS(corr/p) < preci)EXIT B1
	ENDDO B1

	dtr=-t/r/2.d0 ; dtl=t/l/4.d0
	dpr=-(cte2*mstar*2.d0/r**3+cte3*w
	1 +p*(dkapt+dkapro*drot)*dtr)/kap/(1.d0+p*dkapp/kap)
	dpl=-p*(dkapt+dkapro*drot)*dtl /kap/(1.d0+p*dkapp/kap)

	rstar=r ; teff=t

	m=mstar		!raccord en masse a l'exterieur
c	print*,'mstar',mstar

	dml=0.d0 ; dmr=0.d0

	pprec=p		!pour initialisation appel suivant
	
	pt=p ; dptl=dpl ; dptr=dpr
	
c	IF(list)THEN
c	 pt_atm(1)=p ; t_atm(1)=t ; m_atm(1)=m ; tau(1)=tau_b
c	 r_atm(1)=r ; p_atm(1)=p ; dlpp_atm(1)=1.d0
c	ENDIF

	RETURN
	
	END SUBROUTINE lim_tau1
