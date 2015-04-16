
c****************************************************************

	SUBROUTINE conv_cgm_reza(krad,grav,delta,cp,ro,hp,gradrad,
	1    gradad,der,grad,dgradk,dgradgr,dgradel,dgradcp,dgradro,
	2    dgradhp,dgradtaur,dgradrad,dgradad,
	3    gam,dgamk,dgamgr,dgamdel,dgamcp,dgamro,
	4    dgamhp,dgamtaur,dgamrad,dgamad)

c	routine private du module  mod_conv
	
c	calcul du gradient convectif selon CGM : Canuto Goldman Mazitelli
c	ApJ 473, 550-559, 1996 en adoptant la prescription de Bernkopf
c       on se réfère ici au 2 articles suivant :
c       H02 : Heiter et al, A&A, 2002
c       CM  : Canuto Mazitelli ApJ 370, 295, 1991

c       Date:  1/03/2002  version pour CESAM4
c       Reza Samadi, LESIA , Observatoire de Paris 
c	adapte programme conv_cm de P.Morel  Cassini, O.C.A.
c       
c       Updated version (30.04.02) by Reza Samadi
c       New Version in which Eq.(6) of CM is modified so as to take into
c       account the quantity delta=- ( d ln ro / d ln T )p

c	Adaptation:  P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c       entrées :
c	krad : conductivité radiative 4 ac T**3 / 3 kap ro
c	gravité : G m / r**2
c	delta : - ( d ln ro / d ln T )p
c	cp : chaleur spécifique
c	ro : densité
c	hp : échelle de hauteur de pression
c	gradrad : gradient radiatif
c	gradad : gradient adiabatique
c	taur : ep. opt. de la bulle
c	der=.true. : calcul des dérivées

c	alpha : longueur de mélange

c       sorties :
c	gradconv : gradient convectif
c	dgrad* : dérivées
c	gam : le gamma de la convection
c	dgam* : dérivées

c---------------------------------------------------------------------------

	USE mod_donnees, ONLY: alpha
	USE mod_kind	
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: krad, grav, delta, cp, ro, hp,
	1    gradrad, gradad		
	LOGICAL, INTENT(in) :: der
	REAL (kind=dp), INTENT(out) :: grad, dgradk, dgradgr, dgradel,
	1    dgradcp, dgradro, dgradhp, dgradtaur, dgradrad, dgradad,
	3    gam, dgamk, dgamgr, dgamdel, dgamcp, dgamro,
	4    dgamhp, dgamtaur, dgamrad, dgamad	
	
	REAL (kind=dp), PARAMETER :: K0=1.7d0, S_sig=81.d0/2.d0,
	2    ca=10.8654d0, cb=0.00489073d0, ck=0.149888d0, cm=0.189238d0,
	3    cn=1.85011d0, cc=0.0108071d0, cd=0.00301208d0,
	4    ce=0.000334441d0, cf=0.000125d0, cpp=0.72d0, cq=0.92d0, 
	5    cr=1.2d0, ct=1.5d0
	
	REAL (kind=dp) :: l, ki, a, corr, b, sig, dsig,
	1    phi, phip, dphi, f, df, dkik, dkicp, dkiro,
	2    dbhp, dbk, dbcp, dbro, dbgr, dahp, dak, dacp, daro,
	3    dagr, dsighp, dsigk, dsigcp, dsigro, dsiggr, dsigad,
	5    dphihp, dphik,
	6    dphicp, dphiro, dphigr, dphiad, dgrad, sg12, dgams, s,
	7    dS, bS1, bS1m, F1, F2, dF1, dF2, cSp, dSq1, eSr, fSt1

	INTEGER :: iter
	
	LOGICAL, SAVE :: init=.TRUE.
	LOGICAL :: conv 

c--------------------------------------------------------------------
	
 2000	FORMAT(8es10.3)

	IF(init)THEN		!vsal: V/Al
	   init=.FALSE.
	   
	   WRITE(2,1) ; WRITE(*,1)
 1	   FORMAT(/,'gradient dans zone convective calculé selon',/,
	1	'Canuto Goldman Mazitelli (CGM) ApJ 473, 550, 1996') 
	ENDIF
	
	l=alpha*hp		!longueur de mélange
	ki=krad/cp/ro		!conductivité thermometrique
	b=2.d0*l**2/9.d0/ki*sqrt(grav/2.d0/hp*delta) !2 x Eq(6) CM
! eq.(6) corrected so as to include delta:
! delta= - (dln rho / ln T)_P           [RS 30.04.02]
	a=b**2			!c'est le 4a**2 de H02 (Eq. 6)
	
c	initialisations
	
	conv=.FALSE.
	grad=gradad*1.1d0
	iter=0
	phi=1.d30
	
c	Newton-Raphson "d" signifie derivee/ grad	
	
	DO WHILE(.NOT.conv)
	   iter=iter+1
	   IF(iter > 30)THEN
	      PRINT*,'pas de convergence dans conv_cgm'
	      PRINT*,'donnees : krad,grav,cp,ro,hp,gradrad,gradad,gradrad-gradad'	  
	      WRITE(*,2000)krad,grav,cp,ro,hp,gradrad,gradad,gradrad-gradad
	      PRINT*,'non convergence : phi,phip,abs(phi-phip)/phi'
	      WRITE(*,2000)phi,phip,abs(phi-phip)/phi
	      IF(ABS(gradrad-gradad)/gradad < 1.d-3)THEN
		 WRITE(*,2)(gradrad-gradad)/gradad
 2		 FORMAT('convergence forcée (gradrad-gradad)/gradad=',es10.3)
		 EXIT
	      ELSE	   
		 STOP
	      ENDIF
	   ENDIF		
	   phip=phi ; S=S_sig*a*(grad-gradad) ; dS=S_sig*a ! Eq.(6) H02

c       F1  (Eq. 18 H02) :
	   bS1=1d0+cb*S ; bS1m=bS1**cm
	   F1=(K0/1.5d0)**3 * ca * (S**ck) * ((bS1m - 1d0)**cn)
	   dF1=F1* (ck/S + cn*cm*cb / (bS1m-1d0) * (bS1m/bS1) ) ; dF1=dF1* dS
c       F2  (Eq. 20 H02) :
	   cSp=cc*(S**cpp) ; dSq1=1d0+cd*(S**cq) ; eSr= ce*(S**cr)
	   fSt1=1d0+cf*(S**ct) ; F2=1d0 + cSp/dSq1+eSr/fSt1
	   dF2=cpp* cSp/S/dSq1-cSp*cd*cq*(S**(cq-1d0))/(dSq1**2) 
	   dF2=cr * eSr/S/fSt1-eSr*cf*ct*(S**(ct-1d0))/(fSt1**2)+dF2 ; dF2=dF2*dS
c       Phi  (Eq. 17 H02) :
	   phi=F1  * F2 ; dphi= F1*dF2+F2*dF1
	   f=grad+phi*(grad-gradad)-gradrad ! Eq. (63) CM
	   df=1.d0+dphi*(grad-gradad)+phi ; corr=f/df
c       PRINT*,iter ; WRITE(*,2000)corr,abs(phi-phip)/phi
	   grad=grad-corr ; conv=abs(phi-phip)/phi .le. 1.d-10
	ENDDO
c	PAUSE
	sig=a*(grad-gradad)	! Eq.(6) H02
	sg12=sqrt(sig+1.d0) ; gam=(sg12-1.d0)/2.d0
	
	IF(der)THEN	
	   dkik=ki/krad ; dkicp=-ki/cp ; dkiro=-ki/ro
	   
	   dbhp=b*(2.d0*alpha/l-.5d0/hp) ; dbk= -b/ki*dkik
	   dbcp=-b/ki*dkicp ; dbro=-b/ki*dkiro ; dbgr=b*.5/grav
	   
	   dahp=2.d0*b*dbhp ; dak= 2.d0*b*dbk ; dacp=2.d0*b*dbcp
	   daro=2.d0*b*dbro ; dagr=2.d0*b*dbgr
	   
	   sig=a*(grad-gradad)	! Eq.(6) H02
	   dsig=a ; dsighp=sig/a*dahp ; dsigk= sig/a*dak ; dsigcp=sig/a*dacp
	   dsigro=sig/a*daro ; dsiggr=sig/a*dagr ; dsigad=-a

c       F1  (Eq. 18 H02) :
	   bS1=1d0+cb*S ; bS1m=bS1**cm
	   F1=(K0/1.5d0)**3 * ca * (S**ck) * ((bS1m - 1d0)**cn)
	   dF1=F1* (ck/S + cn*cm*cb / (bS1m-1d0) * (bS1m/bS1) )
! ici dF1 = dF1/dS	 
c       F2  (Eq. 20 H02) :
	   cSp=cc*(S**cpp) ; dSq1=1d0+cd*(S**cq) ; eSr= ce*(S**cr)
	   fSt1=1d0+cf*(S**ct) ; F2=1d0 + cSp/dSq1+eSr/fSt1
	   dF2=cpp* cSp/S/dSq1-cSp*cd*cq*(S**(cq-1d0))/(dSq1**2) 
	   dF2=cr * eSr/S/fSt1-eSr*cf*ct*(S**(ct-1d0))/(fSt1**2) + dF2
! ici dF2 = dF2/dS
c       Phi  (Eq. 17 H02) :
	   phi=F1  * F2 ; dphi= F1*dF2+F2*dF1 ! ici dphi = dphi/dS

	   dphihp=dphi* dsighp * S_sig ; dphik =dphi* dsigk * S_sig  
	   dphicp=dphi* dsigcp * S_sig ; dphiro=dphi* dsigro * S_sig
	   dphigr=dphi* dsiggr * S_sig ; dphiad=dphi* dsigad * S_sig 
	   dphi=dphi*dsig*  S_sig ! dphi = dphi/dgrad

	   dgrad=grad-gradad
	   
	   dgradhp=-dphihp*dgrad/(1.d0+phi+dphi*dgrad)
	   dgradk=-dphik*dgrad/(1.d0+phi+dphi*dgrad)
	   dgradcp=-dphicp*dgrad/(1.d0+phi+dphi*dgrad)
	   dgradro=-dphiro*dgrad/(1.d0+phi+dphi*dgrad)
	   dgradgr=-dphigr*dgrad/(1.d0+phi+dphi*dgrad)
	   dgradad=-(dphiad*dgrad-phi)/(1.d0+phi+dphi*dgrad)
	   dgradrad=1.d0/(1.d0+phi+dphi*dgrad)	 
	   dgradel=0.d0 ; dgradtaur=0.d0
	   
	   dgams=.25d0/sg12 ; dgamhp=dgams*(dsighp+dsig*dgradhp)
	   dgamk= dgams*(dsigk+dsig*dgradk) ; dgamcp=dgams*(dsigcp+dsig*dgradcp)
	   dgamro=dgams*(dsigro+dsig*dgradro) ; dgamgr=dgams*(dsiggr+dsig*dgradgr)
	   dgamad=dgams*(dsigad+dsig*dgradad) ; dgamrad=dgams*dsig*dgradrad
	   dgamdel=0.d0 ; dgamtaur=0.d0	 	 	 
	ENDIF
	
	RETURN
	
	END SUBROUTINE conv_cgm_reza
