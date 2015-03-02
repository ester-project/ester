
c****************************************************************

	SUBROUTINE conv_cml(r,krad,grav,cp,ro,hp,gradrad,gradad,der,
	1 grad,dgradk,dgradgr,dgradel,dgradcp,dgradro,
	2 dgradhp,dgradtaur,dgradrad,dgradad,
	3 gam,dgamk,dgamgr,dgamdel,dgamcp,dgamro,
	4 dgamhp,dgamtaur,dgamrad,dgamad)

c	routine private du module mod_conv

c	calcul du gradient convectif selon Canuto Mazitelli ApJ 370, 295, 1991
c	la longueur de mélange est égale à la distance au sommet de la ZC

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c entrées :
c	krad : conductivité radiative 4 ac T**3 / 3 kap ro
c	grav : G m / r**2
c	delta : - ( d ln ro / d ln T )p
c	cp : chaleur spécifique
c	ro : densité
c	hp : échelle de hauteur de pression
c	gradrad : gradient radiatif
c	gradad : gradient adiabatique
c	taur : ep. opt. de la bulle
c	der=.true. : calcul des dérivées

c sorties :
c	grad : gradient convectif
c	dgrad* : dérivées
c	gam : le gamma de la convection
c	dgam* : dérivées

c----------------------------------------------------------------

	USE mod_donnees, ONLY: rsol
	USE mod_kind
	USE mod_variables, ONLY: id_conv, if_conv, r_zc_conv

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) ::cp, gradad, gradrad, grav, hp, krad,
	1 r, ro
	LOGICAL, INTENT(in) :: der
	REAL (kind=dp), INTENT(out) :: grad, dgradk, dgradgr, dgradel,
	1 dgradcp, dgradro, dgradhp, dgradtaur, dgradrad, dgradad,
	2 gam, dgamk, dgamgr, dgamdel, dgamcp, dgamro,
	3 dgamhp, dgamtaur, dgamrad, dgamad

	REAL (kind=dp), PARAMETER :: a1=24.868d0, a2=9.766d-2, epsi=1.d-6,
	1 m=0.14972d0, n=0.18931d0, p=1.8503d0

	REAL (kind=dp) :: a, a2s, a2sn, a2s1, b, corr, dacp, dagr, dahp,
	1 dak, daro, da2s, da2scp, da2sad, da2sgr, da2sn, da2shp, da2sk,
	2 da2snad, da2sncp, da2sngr, da2snhp, da2snk, da2snro, da2sro,
	3 dbcp, dbgr, dbhp, dbk, dbro, df, dgams, dgrad, dkicp, dkik,
	4 dkiro, dphi, dphiad, dphicp, dphigr, dphihp, dphik, dphiro,
	5 dsig, dsigad, dsigcp, dsighp, dsiggr, dsigk, dsigro, f, ki,
	6 l, l0, phi=10.d0, phip, sg12, sig

	INTEGER, PARAMETER :: iter_max=30
	INTEGER :: iter, i

	LOGICAL, SAVE :: init=.TRUE.

c--------------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN    !vsal: V/Al
	init=.FALSE. ; WRITE(2,1) ; WRITE(*,1)
1	FORMAT(/,'gradient dans zone convective calculé selon',/,
	1 'Canuto-Mazitelli ApJ 370, 295, 1991;',/,'longueur de mélange ',
	2 'égale à la distance au sommet de la ZC')
	ENDIF
	
c	le point de calcul est au delà de R*
c	ce cas se présente pour le point courant ZC d'atmosphère lors de
c	l'initialisation car R* = Rrac = cte, on utilise alors conv_cm
	
	IF(r > r_zc_conv(if_conv))THEN
	 grad=gradrad ; dgradrad=1.d0 ; gam=0.d0
	 dgradk=0.d0 ; dgradgr=0.d0 ; dgradel=0.d0 ; dgradcp=0.d0
	 dgradro=0.d0 ; dgradhp=0.d0 ; dgradtaur=0.d0 ; dgradad=0.d0
	 dgamk=0.d0 ; dgamgr=0.d0 ; dgamdel=0.d0 ; dgamcp=0.d0
	 dgamro=0.d0 ; dgamhp=0.d0 ; dgamtaur=0.d0 ; dgamrad=0.d0
	 dgamad=0.d0	  
	 RETURN
	ENDIF
		
c	détermination de la longueur de mélange : distance au sommet de
c	la ZC.
c	r_zc_conv(ifin_conv) est parfois la limite externe de l'enveloppe
c	et non pas R*, rayon bolomètrique de l'étoile qui se trouve dans
c	l'atmosphère

	B1: DO i=id_conv,if_conv,2
	 IF(r >= r_zc_conv(i) .AND. r <= r_zc_conv(i+1))THEN
	  l0=r_zc_conv(i+1)-r
c	  PRINT*,i
c	  WRITE(*,2000)r_zc_conv(i),r,r_zc_conv(i+1)
	  EXIT B1
	 ENDIF
	ENDDO B1
	
	l=l0*rsol
	ki=krad/cp/ro       !conductivite thermometrique
	b=2.d0*l**2/9.d0/ki*SQRT(grav/2.d0/hp)      !2 x (6)
	a=b**2          !c'est le 4a**2 de (5CM)
c	WRITE(*,2000)l0,l,hp,ki,b,a ; PAUSE'l0, l, Hp, ki, b ,a'

c	initialisations

	grad=gradad*1.1d0 ; iter=0

c	Newton-Raphson "d" signifie dérivée/ grad 
c	si phi=0, grad=grad_rad cf. eq. 63 de CM. (cas où l est petit) 

	B2: DO
	 iter=iter+1
	 IF(iter > iter_max)THEN
	  WRITE(*,*)'pas de convergence dans conv_cml, iter,iter_max ',
	1 iter,iter_max	
	  PRINT*,id_conv,if_conv
	  WRITE(*,*)'l0,l,r,r_zc_conv(id_conv:if_conv)'
	  WRITE(*,2000)l0,l,r,r_zc_conv(id_conv:if_conv)
	  WRITE(*,*)'r,gradad,gradrad,gradad-gradrad'
	  WRITE(*,*)r,gradad,gradrad,gradad-gradrad
	  PRINT*,'krad,gravité,cp,ro,hp,gradrad,gradad,grad'
	  WRITE(*,2000)krad,grav,cp,ro,hp,gradrad,gradad,grad
	  PRINT*,'ki,a,b,phi,phip,,f,df'
	  WRITE(*,2000)ki,a,b,phi,phip,f,df ; STOP
	 ENDIF
	 phip=phi ; sig=a*(grad-gradad) !(5CM)
	 dsig=a ; a2s=a2*sig            !(32CM)
	 da2s=a2*dsig ; a2s1=1.d0+a2s ; a2sn=a2s1**n
	 da2sn=n*a2sn/a2s1*da2s ; phi=a1*sig**m*(a2sn-1.d0)**p

c	si phi=0, grad=grad_rad cf. eq. 63 de CM. (cas où l est petit) 
	 IF(phi <= 0.d0)THEN
c	  PRINT*,'cas phi=0' ; WRITE(*,2000)phi
	  grad=gradrad ; dgradrad=1.d0 ; gam=0.d0
	  dgradk=0.d0 ; dgradgr=0.d0 ; dgradel=0.d0 ; dgradcp=0.d0
	  dgradro=0.d0 ; dgradhp=0.d0 ; dgradtaur=0.d0 ; dgradad=0.d0
	  dgamk=0.d0 ; dgamgr=0.d0 ; dgamdel=0.d0 ; dgamcp=0.d0
	  dgamro=0.d0 ; dgamhp=0.d0 ; dgamtaur=0.d0 ; dgamrad=0.d0
	  dgamad=0.d0	  
	  RETURN
	 ELSE	 
	  dphi=phi*(m*dsig/sig+p*da2sn/(a2sn-1.d0))
	  f=grad+phi*(grad-gradad)-gradrad       !(63CM)
	  df=1.d0+dphi*(grad-gradad)+phi ; corr=f/df

	  IF(.FALSE.)THEN	 
c	  IF(r >= r_zc_conv(if_conv))THEN
	   WRITE(*,*)'a,a2,n,p,m,grad-gradad,grad-gradrad,iter',iter
	   WRITE(*,2000)a,a2,n,p,m,grad-gradad,grad-gradrad
	   WRITE(*,*)'sig,a2s,a2s1,a2sn,sig**m'
	   WRITE(*,*)sig,a2s,a2s1,a2sn,sig**m
	   WRITE(*,*)'grad-gradad,sig,a2sn-1.d0'
	   WRITE(*,*)grad-gradad,sig,a2sn-1.d0
	   PRINT*,'corr,phi,phip,ABS(phi-phip),f,df'
	   WRITE(*,2000)corr,phi,phip,ABS(phi-phip)/phi,f,df
	  ENDIF
	  grad=grad-corr
c	  IF(ABS((phi-phip)/phi) <= epsi)EXIT B2
	  IF(ABS(corr) <= epsi)EXIT B2
	 ENDIF  
	ENDDO B2

	sg12=SQRT(sig+1.d0) ; gam=(sg12-1.d0)/2.d0

	IF(der)THEN
	 dkik=ki/krad ; dkicp=-ki/cp ; dkiro=-ki/ro

	 dbhp=-b*0.5d0/hp ; dbk= -b/ki*dkik
	 dbcp=-b/ki*dkicp ; dbro=-b/ki*dkiro ; dbgr=b*0.5d0/grav

	 dahp=2.d0*b*dbhp ; dak= 2.d0*b*dbk ; dacp=2.d0*b*dbcp
	 daro=2.d0*b*dbro ; dagr=2.d0*b*dbgr

	 sig=a*(grad-gradad)    !(5CM)
	 dsig=a ; dsighp=sig/a*dahp ; dsigk= sig/a*dak ; dsigcp=sig/a*dacp
	 dsigro=sig/a*daro ; dsiggr=sig/a*dagr ; dsigad=-a

	 a2s=sig*a2     !(32CM)
	 a2s1=1.d0+a2s ; da2s=dsig*a2 ; da2shp=dsighp*a2 ; da2sk= dsigk*a2
	 da2scp=dsigcp*a2 ; da2sro=dsigro*a2 ; da2sgr=dsiggr*a2
	 da2sad=dsigad*a2

	 a2sn=a2s1**n ; da2sn=  n*a2sn/a2s1*da2s
	 da2snhp=n*a2sn/a2s1*da2shp ; da2snk= n*a2sn/a2s1*da2sk
	 da2sncp=n*a2sn/a2s1*da2scp ; da2snro=n*a2sn/a2s1*da2sro
	 da2sngr=n*a2sn/a2s1*da2sgr ; da2snad=n*a2sn/a2s1*da2sad

	 phi=a1*sig**m*(a2sn-1.d0)**p
	 dphi=  phi*(m*dsig  /sig+p*da2sn  /(a2sn-1.d0))
	 dphihp=phi*(m*dsighp/sig+p*da2snhp/(a2sn-1.d0))
	 dphik= phi*(m*dsigk /sig+p*da2snk /(a2sn-1.d0))
	 dphicp=phi*(m*dsigcp/sig+p*da2sncp/(a2sn-1.d0))
	 dphiro=phi*(m*dsigro/sig+p*da2snro/(a2sn-1.d0))
	 dphigr=phi*(m*dsiggr/sig+p*da2sngr/(a2sn-1.d0))
	 dphiad=phi*(m*dsigad/sig+p*da2snad/(a2sn-1.d0))

	 dgrad=grad-gradad

	 dgradhp=-dphihp*dgrad/(1.d0+phi+dphi*dgrad)
	 dgradk=-dphik*dgrad/(1.d0+phi+dphi*dgrad)
	 dgradcp=-dphicp*dgrad/(1.d0+phi+dphi*dgrad)
	 dgradro=-dphiro*dgrad/(1.d0+phi+dphi*dgrad)
	 dgradgr=-dphigr*dgrad/(1.d0+phi+dphi*dgrad)
	 dgradad=-(dphiad*dgrad-phi)/(1.d0+phi+dphi*dgrad)
	 dgradrad=1.d0/(1.d0+phi+dphi*dgrad)     
	 dgradel=0.d0 ; dgradtaur=0.d0

	 dgams=0.25d0/sg12
	 dgamhp=dgams*(dsighp+dsig*dgradhp)
	 dgamk= dgams*(dsigk+dsig*dgradk)
	 dgamcp=dgams*(dsigcp+dsig*dgradcp)
	 dgamro=dgams*(dsigro+dsig*dgradro)
	 dgamgr=dgams*(dsiggr+dsig*dgradgr)
	 dgamad=dgams*(dsigad+dsig*dgradad)
	 dgamrad=dgams*dsig*dgradrad
	 dgamdel=0.d0 ; dgamtaur=0.d0
	ENDIF
	
	RETURN

	END SUBROUTINE conv_cml
