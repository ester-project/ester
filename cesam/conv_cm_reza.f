
c****************************************************************

	SUBROUTINE conv_cm_reza(krad,grav,delta,cp,ro,hp,gradrad,gradad,der,
	1 grad,dgradk,dgradgr,dgradel,dgradcp,dgradro,
	2 dgradhp,dgradtaur,dgradrad,dgradad,
	3 gam,dgamk,dgamgr,dgamdel,dgamcp,dgamro,
	4 dgamhp,dgamtaur,dgamrad,dgamad)

c	routine private du module mod_conv
	
c	calcul du gradient convectif selon Canuto Mazitelli ApJ 370, 295, 1991

c       Updated version (30.04.02) by 
c	Reza Samadi, LESIA , Observatoire de Paris 
c       New Version in which Eq.(6) of CM is modified so as to take into
c       account the quantity delta=- ( d ln ro / d ln T )p
c       
c entree :
c	krad : conductivité radiative 4 ac T**3 / 3 kap ro
c	gravite : G m / r**2
c	delta : - ( d ln ro / d ln T )p
c	cp : chaleur spécifique
c	ro : densité
c	hp : échelle de hauteur de pression
c	gradrad : gradient radiatif
c	gradad : gradient adiabatique
c	taur : ep. opt. de la bulle
c	der=.true. : calcul des dérivées

c	alpha : longueur de melange (dans le common /modele/)

c sortie :
c	gradconv : gradient convectif
c	dgrad* : dérivées
c	gam : le gamma de la convection
c	dgam* : dérivées

c	Auteur: Reza Samadi, LESIA , Observatoire de Paris
c	Adaptation: P. Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c-----------------------------------------------------------------------

	USE mod_donnees, ONLY: alpha
	USE mod_kind	
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: krad, grav, delta, cp, ro, hp,
	1 gradrad, gradad		
	LOGICAL, INTENT(in) :: der
	REAL (kind=dp), INTENT(out) :: grad, dgradk, dgradgr, dgradel,
	1 dgradcp, dgradro, dgradhp, dgradtaur, dgradrad, dgradad,
	2 gam, dgamk, dgamgr, dgamdel, dgamcp, dgamro,
	3 dgamhp, dgamtaur, dgamrad, dgamad	
	
	REAL (kind=dp), PARAMETER :: a1=24.868d0, a2=9.766d-2, m=0.14972d0,
	1 n=0.18931d0, p=1.8503d0	
	
	REAL (kind=dp) :: l, ki, a, corr, b, a2s, da2s, sig, dsig, a2sn,
	1 da2sn, a2s1, phi, phip, dphi, f, df, dkik, dkicp, dkiro,
	2 dbhp, dbk, dbcp, dbro, dbgr, dahp, dak, dacp, daro,
	3 dagr, dsighp, dsigk, dsigcp, dsigro, dsiggr, dsigad,
	4 da2shp, da2sk, da2scp, da2sro, da2sgr, da2sad, da2snhp,
	5 da2snk, da2sncp, da2snro, da2sngr, da2snad, dphihp, dphik,
	6 dphicp, dphiro, dphigr, dphiad, dgrad, sg12, dgams
	
	INTEGER :: iter
		
	LOGICAL, SAVE :: init=.TRUE.
	LOGICAL :: conv 

c--------------------------------------------------------------------
	
2000	FORMAT(8es10.3)

	IF(init)THEN	!vsal: V/Al
	 init=.FALSE. ; WRITE(2,1) ; WRITE(*,1)
1	 FORMAT(/,'gradient dans zone convective calcule selon',/,
	1 ' Canuto-Mazitelli ApJ 370, 295, 1991')
	ENDIF
	
	l=alpha*hp		!longueur de mélange
	ki=krad/cp/ro		!conductivité thermometrique
	b=2.d0*l**2/9.d0/ki*sqrt(grav/2.d0/hp*delta)		!2 x (6)
	! eq.(6) corrected so as to include delta:
	! delta= - (dln rho / ln T)_P                     [RS 30.04.02]
	a=b**2			!c'est le 4a**2 de CM 
	
c	initialisations
		
	conv=.FALSE. ; grad=gradad*1.1d0 ; iter=0 ; phi=1.d30
	
c	Newton-Raphson "d" signifie dérivée/ grad
c	Newton-Raphson "d" means dérivative/ grad	
	
	do while(.NOT.conv)
	 iter=iter+1
	 IF(iter > 30)THEN
	  PRINT*,'pas de convergence dans conv_cm'
	  PRINT*,'données : krad,grav,cp,ro,hp,gradrad,gradad,gradrad-gradad'	  
	  WRITE(*,2000)krad,grav,cp,ro,hp,gradrad,gradad,gradrad-gradad
	  PRINT*,'non convergence : phi,phip,abs(phi-phip)/phi'
	  WRITE(*,2000)phi,phip,abs(phi-phip)/phi 
	  IF(ABS(gradrad-gradad)/gradad < 1.d-3)THEN
	     WRITE(*,2)(gradrad-gradad)/gradad
 2	     FORMAT('convergence forcée (gradrad-gradad)/gradad=',es10.3)
	     EXIT
	  ELSE	   
	     STOP
	  ENDIF
	 ENDIF
	 phip=phi ; sig=a*(grad-gradad)	!(5)
	 dsig=a ; a2s=a2*sig		!(32)
	 da2s=a2*dsig ; a2s1=1.d0+a2s ; a2sn=a2s1**n
	 da2sn=n*a2sn/a2s1*da2s ; phi=a1*sig**m*(a2sn-1.d0)**p
	 dphi=phi*(m*dsig/sig+p*da2sn/(a2sn-1.d0))
	 f=grad+phi*(grad-gradad)-gradrad		!(63)
	 df=1.d0+dphi*(grad-gradad)+phi ; corr=f/df
c	 PRINT*,iter
c	 WRITE(*,2000)corr,abs(phi-phip)/phi
	 grad=grad-corr ; conv=abs(phi-phip)/phi .le. 1.d-10
	ENDDO
c	PAUSE
	
	sg12=sqrt(sig+1.d0) ; gam=(sg12-1.d0)/2.d0
	
	IF(der)THEN	
	 dkik=ki/krad ; dkicp=-ki/cp ; dkiro=-ki/ro
	
	 dbhp=b*(2.d0*alpha/l-0.5d0/hp) ; dbk= -b/ki*dkik
	 dbcp=-b/ki*dkicp ; dbro=-b/ki*dkiro ; dbgr=b*0.5d0/grav
	
	 dahp=2.d0*b*dbhp ; dak= 2.d0*b*dbk ; dacp=2.d0*b*dbcp
	 daro=2.d0*b*dbro ; dagr=2.d0*b*dbgr
	 
	 sig=a*(grad-gradad)	!(5)
	 dsig=a ; dsighp=sig/a*dahp ; dsigk= sig/a*dak ; dsigcp=sig/a*dacp
	 dsigro=sig/a*daro ; dsiggr=sig/a*dagr ; dsigad=-a
	
	 a2s=sig*a2		!(32)
	 a2s1=1.d0+a2s ; da2s=dsig*a2 ; da2shp=dsighp*a2
	 da2sk= dsigk *a2 ; da2scp=dsigcp*a2 ; da2sro=dsigro*a2
	 da2sgr=dsiggr*a2 ; da2sad=dsigad*a2
	 
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
	 dgradrad=1.d0/(1.d0+phi+dphi*dgrad) ; dgradel=0.d0 ; dgradtaur=0.d0
	 
	 dgams=.25d0/sg12 ; dgamhp=dgams*(dsighp+dsig*dgradhp)
	 dgamk= dgams*(dsigk+dsig*dgradk)	 
	 dgamcp=dgams*(dsigcp+dsig*dgradcp)
	 dgamro=dgams*(dsigro+dsig*dgradro)	 
	 dgamgr=dgams*(dsiggr+dsig*dgradgr)	 
	 dgamad=dgams*(dsigad+dsig*dgradad)
	 dgamrad=dgams*dsig*dgradrad ; dgamdel=0.d0 ; dgamtaur=0.d0 
	ENDIF
	
	RETURN
		
	END SUBROUTINE conv_cm_reza
