!****************************************************************

	SUBROUTINE v2_cgm(pr,te,lu,ma,ra,l,xchim,v2)
	
!	calcul du gradient convectif selon CGM : Canuto Goldman Mazitelli
!	ApJ 473, 550-559, 1996 en adoptant la prescription de Bernkopf
!       on se réfère ici au 2 articles suivant :
!       H02 : Heiter et al, A&A, 2002
!       CM  : Canuto Mazitelli ApJ 370, 295, 1991

!       Date:  1/03/2002  version pour CESAM4
!       Reza Samadi, LESIA , Observatoire de Paris 
!	adapte programme conv_cm de P.Morel  Cassini, O.C.A.
!       
!       Updated version (30.04.02) by Reza Samadi
!       New Version in which Eq.(6) of CM is modified so as to take into
!       account the quantity delta=- ( d ln ro / d ln T )p

!	Adaptation:  P.Morel, Département J.D. Cassini, O.C.A.
!	CESAM2k

!       entrées :
!	krad : conductivité radiative 4 ac T**3 / 3 kap ro
!	gravité : G m / r**2
!	delta : - ( d ln ro / d ln T )p
!	cp : chaleur spécifique
!	ro : densité
!	hp : échelle de hauteur de pression
!	gradrad : gradient radiatif
!	gradad : gradient adiabatique

!	alpha : longueur de mélange

!       sorties :
!	gradconv : gradient convectif
!	gam : le gamma de la convection

!---------------------------------------------------------------------------

	USE mod_kind
	USE mod_etat, ONLY: etat
	USE mod_opa, ONLY: opa
	USE mod_donnees, ONLY: lsol,rsol,msol,pi,aradia,clight,g
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: pr, te, lu, ma, ra, l
	REAL (kind=dp), DIMENSION(:), INTENT(in) :: xchim
	REAL (kind=dp), INTENT(out) :: v2
	
	REAL (kind=dp), PARAMETER :: a1=24.868d0, a2=9.766d-2, m=0.14972d0,
	1    n=0.18931d0, p=1.8503d0, K0=1.7d0, S_sig=81.d0/2.d0,
	2    ca=10.8654d0, cb=0.00489073d0, ck=0.149888d0, cm=0.189238d0,
	3    cn=1.85011d0, cc=0.0108071d0, cd=0.00301208d0,
	4    ce=0.000334441d0, cf=0.000125d0, cpp=0.72d0, cq=0.92d0, 
	5    cr=1.2d0, ct=1.5d0
	
	REAL (kind=dp) :: ki, a, corr, b, sig, dsig,
	1    phi, phip, dphi, f, df, dkik, dkicp, dkiro,
	2    dbhp, dbk, dbcp, dbro, dbgr, dahp, dak, dacp, daro,
	3    dagr, dsighp, dsigk, dsigcp, dsigro, dsiggr, dsigad,
	5    dphihp, dphik,grad,gam,pturb,
	6    dphicp, dphiro, dphigr, dphiad, dgrad, sg12, dgams, S,
	7    dS, bS1, bS1m, F1, F2, dF1, dF2, cSp, dSq1, eSr, fSt1,
	8    F3, F4, F5, ro,drop,drot,drox,u,dup,dut,dux,kap,dkapt,dkapro,dkapx,
	1    delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,grav,hp,
	2    gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1,gradrad,krad
	REAL (kind=dp), SAVE :: cte1, cte8, cte13

	INTEGER :: iter
	
	LOGICAL, SAVE :: init=.TRUE.
	LOGICAL :: conv 

!--------------------------------------------------------------------
	
 2000	FORMAT(8es10.3)

	IF(init) THEN
	   init = .FALSE.
	   cte1 = 4.d0/3.d0*aradia*clight
	   cte8 = lsol/4.d0/pi/rsol/rsol
	   cte13 = g*msol/rsol/rsol
	ENDIF

	CALL etat(pr,te,xchim,.FALSE.,ro,drop,drot,drox,u,dup,dut,dux,
	1    delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2    gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)

	CALL opa(xchim,te,ro,kap,dkapt,dkapro,dkapx)


	krad = cte1/kap/ro*te**3	!diffusivite radiative

	IF(ma*lu*ra /= 0.d0)THEN !hors du centre
	   grav = cte13*ma/ra**2	! sem wrot
	   
!       rot=-cte2*w**2*r	!gravite effective avec rotation
!       grav=grav+rot	!remarque de N. Audard
	   
!       echelle de hauteur de pression
	   
	   hp = pr/grav/ro	!echelle de hauteur de pression
	   gradrad = cte8*lu*hp/ra**2/krad/te !gradient radiatif

	ELSE
	   hp = 0.d0 ; grav = 0.d0 ; gradrad = 0.d0
	ENDIF
	
	
	ki = krad/cp/ro		!conductivité thermometrique
	b = 2.d0*l**2/9.d0/ki*SQRT(grav/2.d0/hp*delta) !2 x Eq(6) CM
! eq.(6) corrected so as to include delta:
! delta= - (dln rho / ln T)_P           [RS 30.04.02]
	a = b**2			!!'est le 4a**2 de H02 (Eq. 6)
	
!	initialisations
	
	conv = .FALSE.
	grad = gradad*1.1d0
	iter = 0
	phi = 1.d30
	
!	Newton-Raphson "d" signifie derivee/ grad	
	
	DO WHILE(.NOT.conv)
	   iter = iter + 1
	   IF(iter > 30)THEN
	      PRINT*,'pas de convergence dans v2_cgm'
	      PRINT*,'l/hp =',l/hp
	      PRINT*,'donnees : krad,grav,cp,ro,hp,gradrad,gradad'	  
	      WRITE(*,2000)krad,grav,cp,ro,hp,gradrad,gradad
	      PRINT*,'non convergence : phi,phip,abs(phi-phip)/phi'
	      WRITE(*,2000)phi,phip,abs(phi-phip)/phi ; STOP
	   ENDIF
	   phip = phi
	   S = S_sig*a*(grad - gradad)
	   dS = S_sig*a		! Eq.(6) H02

!       F1  (Eq. 18 H02) :
	   bS1 = 1d0 + cb*S
	   bS1m = bS1**cm
	   F1 = (K0/1.5d0)**3*ca*(S**ck)*((bS1m - 1d0)**cn)
	   dF1 = F1*(ck/S + cn*cm*cb/(bS1m-1d0)*(bS1m/bS1))
	   dF1 = dF1*dS

!       F2  (Eq. 20 H02) :
	   cSp = cc*(S**cpp)
	   dSq1 = 1d0 + cd*(S**cq)
	   eSr = ce*(S**cr)
	   fSt1 = 1d0 + cf*(S**ct)
	   F2 = 1d0 + cSp/dSq1 + eSr/fSt1
	   dF2 = cpp*cSp/S/dSq1 - cSp*cd*cq*(S**(cq - 1.d0))/(dSq1**2) 
	   dF2 = cr*eSr/S/fSt1 - eSr*cf*ct*(S**(ct - 1.d0))/(fSt1**2) + dF2
	   dF2 = dF2*dS

!       Phi  (Eq. 17 H02) :
	   phi = F1*F2
	   dphi = F1*dF2 + F2*dF1
	   f = grad + phi*(grad - gradad) - gradrad ! Eq. (63) CM
	   df = 1.d0 + dphi*(grad - gradad) + phi
	   corr = f/df
!       PRINT*,iter ; WRITE(*,2000)corr,abs(phi-phip)/phi
	   grad = grad - corr
	   conv = ABS(phi - phip)/phi <= 1.d-10
	ENDDO

	sig = a*(grad - gradad)	! Eq.(6) H02
	sg12 = SQRT(sig + 1.d0)
	gam = (sg12 - 1.d0)/2.d0

!	F3  (Eq. 89 CGM)
	F3 = (K0/1.5d0)**3.d0*0.00101392d0*S**2.d0/
	1    (1.d0 + SQRT(1.d0 + 0.000017848d0*S**2.d0))
	
!	F4  (Eq. 90 CGM)
	F4 = 6.39899d0 + 2.256815d0*(-1.d0 + 0.000777055d0*S**0.868589d0)/ 
	1    (1.d0 + 0.000777055d0*S**0.868589d0)

!	F5  (Eq. 92 CGM)
	F5=1.49168d0 + 0.45185d0*(-1.d0 + 0.00111378d0*S**0.868589d0)/ 
	1    (1.d0 + 0.00111378d0*S**0.868589d0)

!	v^2  (Eq. 88 CGM)
	v2 = (ki/l)**2.d0*F3*F4

!	Pturb  (Eq. 91 CGM)
	pturb = ro*(ki/l)**2.d0*F3*F5	
	
	END SUBROUTINE v2_cgm

