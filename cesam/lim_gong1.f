
c****************************************************************

	SUBROUTINE lim_gong1(l,r,xchim,pt,dptl,dptr,t,dtl,dtr,m,dml,
	1 dmr,p,dpl,dpr,teff)

c routine public du module mod_atm

c calcul de la condition limite pour gong cas 1
c relations Pext(l,r), Text(l,r), Mext(l,r) et derivees

c modifs :
c 19 11 99 : suppression de nh1, nhe1, nhe2, lamb

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c entree
c	list=.true. : on calcule p_atm,t_atm,r_atm,tau,m_atm	
c	l :luminosité
c	r :rayon
c	xchim : composition chimique par gramme
c	mstar: masse avec perte de masse

c sortie
c	pt : Ptot
c	p : Pgaz
c	t : température
c	m : masse
c	dptl, dptr, dtl, dtr, dml, dmr, dpl,
c	dpr : dérivées pt,t,m,p / l, r
c	teff : température effective
c	rstar : rayon externe
c	pt_atm,t_atm,r_atm,tau,m_atm,p_atm : sans objet

c----------------------------------------------------------------

	USE mod_donnees, ONLY : aradia, clight, g, Krot, langue, lsol, msol,
	1 n_atm, pi, rsol
	USE mod_etat, ONLY : etat
	USE mod_kind
	USE mod_opa, ONLY : opa
	USE mod_variables, ONLY : dim_ch, mstar, rota, rstar, wrot
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim
	REAL (kind=dp), INTENT(in) :: l, r
	REAL (kind=dp), INTENT(out) :: pt, dptl, dptr, t, dtl, dtr,
	1 m, dml, dmr, p, dpl, dpr, teff		
	REAL (kind=dp), SAVE :: cte1, cte2, pprec, cte20, cte3, mstarp=-1.d0
	REAL (kind=dp) :: ro, corr, w,
	1 dux, drox, dkapx, u, kap, dkapp, dkapt, drop, drot, dup, dut,
	2 dkapro, delta, deltap, deltat, deltax,
	3 cp, dcpp, dcpt, dcpx, gradad, dgradadp, dgradadt, dgradadx,	
	4 alfa, beta, gamma1, lambda
	
	INTEGER :: ntour
	
	LOGICAL, SAVE :: init=.FALSE.
	
c------------------------------------------------------------------------

2000	FORMAT(8es10.3)

c	WRITE(*,*)'entrée lim_gong1'

	IF(.NOT. init)THEN
	 init=.TRUE. ; lambda=6.d0 ; beta=7.22d0
c	 WRITE(*,*)'entrer beta' ; read(5,*)beta
	 WRITE(*,1)lambda,beta
1	 FORMAT('lambda=',es10.3,' beta=',es10.3)
	 WRITE(2,2)lambda,beta
2	 FORMAT(/,t10,'pour GONG 1  lambda=',es10.3,' beta=',es10.3,//)
	 cte1=lsol/rsol/rsol/lambda/aradia/clight*4.d0/pi/4.d0
	 cte20=beta*g*msol/rsol/rsol ; cte3=2.d0/3.d0*rsol*beta	 
	 pprec=1.d7	!valeur provisoire de pext
	 
c allocations des tableaux pt_atm etc..	 
	 n_atm=0
	 ALLOCATE(bp_atm(0,0),bp_atm_t(0,0),pt_atm(0),t_atm(0),
	1 m_atm(0),tau(0),p_atm(0),r_atm(0),dlpp_atm(0),x_atm(0),x_atm_t(0),
	2 xt_atm(0),xt_atm_t(0))
	ENDIF

c	PRINT*,'lim_gong1',mstar
	IF(mstar /= mstarp)THEN
	 mstarp=mstar ; cte2=cte20*mstar
	ENDIF

c résolution de kap p r**2 = G m beta par iteration newton-raphson
	t=(cte1*l/r**2)**0.25d0

c	WRITE(*,*)'dans lim_gong1 pprec,r,l,xchim(1)'
c	WRITE(*,2000)pprec,r,l,xchim(1)
c	WRITE(*,*)'dans lim_gong1 p,corr,t,r,l,mtot,xchim(1),kap,dkapp'
 	 
c Omega
	SELECT CASE(Krot)
	CASE(0,1,2)
	 w=wrot
	CASE(3)
	 w=rota(1,dim_ch)	 
	END SELECT
	w=w**2 

c initialisations	
	p=pprec	; ntour=0
	B1: DO
	 ntour=ntour+1
	 IF(ntour > 10)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1003) ; WRITE(2,1003) ; STOP
1003	   FORMAT('STOP, no convergence for pext in lim_gong1')	  
	  CASE DEFAULT
	   WRITE(*,3) ; WRITE(2,3) ; STOP
3	   FORMAT('ARRET, pas de convergence pour pext dans lim_gong1')
	  END SELECT
	 ENDIF

	 CALL etat(p,t,xchim,.true.,
	1 ro,drop,drot,drox,u,dup,dut,dux,
	2 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	 CALL opa(xchim,t,ro,kap,dkapt,dkapro,dkapx)
	 dkapp=dkapro*drop

	 corr=(kap*p-cte2/r**2+cte3*w*r)/(dkapp*p+kap) ; p=p-corr
c	 WRITE(*,2000)p,corr,t,r,l,xchim(1),kap,dkapp
	 IF(ABS(corr/p) <= 1.d-8)EXIT B1
	ENDDO B1

	dtr=-t/r/2.d0 ; dtl=t/l/4.d0
	dpr=-(cte2*2.d0/r**3+cte3*w
	1 +p*(dkapt+dkapro*drot)*dtr)/kap/(1.d0+p*dkapp/kap)
	dpl=-p*(dkapt+dkapro*drot)*dtl /kap/(1.d0+p*dkapp/kap)
	
	rstar=r ; teff=t

	m=mstar		!raccord en masse a l'exterieur
	
c	PRINT*,'mstar',mstar

	dmr=0.d0 ; dml=0.d0

	pprec=p		!pour initialisation appel suivant
	
	pt=p ; dptl=dpl ; dptr=dpr

c	IF(list)THEN	
c	 pt_atm(1)=p ; t_atm(1)=t ; m_atm(1)=m ; tau(1)=1.d0
c	 r_atm(1)=r ; p_atm(1)=p ; dlpp_atm(1)=1.d0
c	ENDIF

	RETURN

	END SUBROUTINE lim_gong1
