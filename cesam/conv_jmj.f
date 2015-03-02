
c****************************************************************

	SUBROUTINE conv_jmj(krad,gravite,delta,cp,ro,hp,taur,gradrad,
	1 gradad,der,gradconv,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
	2 dgradhp,dgradtaur,dgradgrad,dgradgad,
	3 gam,dgamkra,dgamgra,dgamdel,dgamcp,dgamro,
	4 dgamhp,dgamtaur,dgamgrad,dgamgad)

c	subroutine private du module mod_conv
	
c	calcul du gradient convectif

c	formulation de:
c	J.Provost, Département J.D. Cassini, O.C.A., Observatoire de Nice 
c	MJo Goupil DASGAL Observatoire de Meudon
c	Correction erreur sur b S. Brun 30 08 96

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

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

c	alpha : longueur de mélange

c sortie :
c	gradconv : gradient convectif
c	dgrad* : dérivées
c	gam : le gamma de la convection
c	dgam* : dérivées

c----------------------------------------------------------------

	USE mod_donnees, ONLY: alpha, langue
	USE mod_kind
	
	IMPLICIT NONE
		
	REAL (kind=dp), intent(in) :: krad, gravite, delta, cp, ro, hp, taur,
	1 gradrad, gradad		
	LOGICAL, INTENT(in) :: der
	REAL (kind=dp), INTENT(out) :: gradconv, dgradkra, dgradgra, dgradel,
	1 dgradcp, dgradro, dgradhp, dgradtaur, dgradgrad, dgradgad,
	3 gam, dgamkra, dgamgra, dgamdel, dgamcp, dgamro,
	4 dgamhp, dgamtaur, dgamgrad, dgamgad	

	REAL (kind=dp), PARAMETER :: vsal=2.d0/9.d0, phi0=3.d0/2.d0/(3.d0*vsal),
	1 ksi0=1.d0/72.d0*(3.d0*vsal)**2, epsi=1.d-6

	REAL (kind=dp) :: ksi, alpha4, phi, b, dbgra, dbkra, dbhp, dbdel,
	1 dbro, dbcp, gradmd, a0, da0gra, da0kra, da0hp, da0del,
	2 da0ro, da0cp, da0rad, da0ad, gam1, gam2, corr, da0taur,
	3 dbtaur, dal4taur, dksitaur, dvsltaur, vsal32,
	4 dphitaur, dvst, dvstaur, gam3

	INTEGER :: iter
	
	LOGICAL, SAVE :: init=.TRUE.
	
c----------------------------------------------------------------------------
	
2000	FORMAT(8es10.3)

	IF(init)THEN	!vsal: V/Al
	 init=.FALSE.
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1001)vsal,ksi0,phi0 ; WRITE(2,1001)vsal,ksi0,phi0
1001	  FORMAT('Convective gradient computed according to MLT : ',/,
	1 'one takes account of the optical depth of the bubble',/,
	2 'V/Al= ',es10.3,', ksi0=',es10.3,', phi0=',es10.3)	 
	 CASE DEFAULT	 
	  WRITE(2,1)vsal,ksi0,phi0 ; WRITE(*,1)vsal,ksi0,phi0
1	  FORMAT('gradient convectif calcule par MLT : ',/,
	1 'formulation tenant compte de l''ep. op. de la bulle',/,
	2 'V/Al= ',es10.3,', ksi0=',es10.3,', phi0=',es10.3)
	 END SELECT	
	ENDIF

	dvst=1.d0+2.d0/3.d0/vsal/taur**2
	phi=phi0/dvst
	vsal32=dvst**2
	ksi=ksi0*vsal32
	alpha4=ksi*alpha**4
	b=alpha4*gravite/krad**2*hp**3*delta*(ro*cp)**2

c	racine réelle de la cubique par newton raphson
c	search of the real root of the cubic by newton raphson algorithm

c	WRITE(*,*)'ksi,gravité,delta,ro,cp,krad,hp,b,gradrad,gradad'
c	WRITE(*,2000)ksi,gravite,delta,ro,cp,krad,hp,b,
c	1 gradrad,gradad
	 
	gradmd=gradrad-gradad
	a0=b*gradmd			!B*gradrad-gradad

	gam=abs(a0/phi)**(1.d0/3.d0)	!initialisation de la racine
	iter=0
	
	B1: DO
	 iter=iter+1
	 IF(iter > 10)THEN
	  SELECT CASE(langue)
	  CASE('english')	  
	   WRITE(*,1002)iter,corr,gam,abs(corr/gam),ksi,gravite,delta,ro,
	1  cp,krad,hp,b,gradrad,gradad
	   WRITE(2,1002)iter,corr,gam,abs(corr/gam),ksi,gravite,delta,ro,
	1   cp,krad,hp,b,gradrad,gradad	  
1002	   FORMAT('STOP, divergence in conv_jmj after ',i3,' iterations'
	1  ,/,'Delta gamma = ',es10.3,', gamma = ',es10.3,', precision = ',
	2  es10.3,/,'Inputs : ksi = ',es10.3,' gravité = ',es10.3,/,
	3  'delta = ',es10.3,', ro = ',es10.3,', cp = ',es10.3,
	4  ', krad = ',es10.3,', hp = ',es10.3,/,'b = ',es10.3,
	5  ', gradrad = ',es10.3,', gradad = ',es10.3)	  
	  CASE DEFAULT
	   WRITE(*,2)iter,corr,gam,abs(corr/gam),ksi,gravite,delta,ro,cp,
	1  krad,hp,b,gradrad,gradad
	   WRITE(2,2)iter,corr,gam,abs(corr/gam),ksi,gravite,delta,ro,cp,
	1  krad,hp,b,gradrad,gradad	  
2	   FORMAT('ARRET, divergence dans conv_jmj après ',i3,' itérations'
	1  ,/,'Delta gamma = ',es10.3,', gamma = ',es10.3,', précision = ',
	2  es10.3,/,'Données : ksi = ',es10.3,' gravité = ',es10.3,/,
	3  'delta = ',es10.3,', ro = ',es10.3,', cp = ',es10.3,
	4  ', krad = ',es10.3,', hp = ',es10.3,/,'b = ',es10.3,
	5  ', gradrad = ',es10.3,', gradad = ',es10.3)
	  END SELECT	 
	  STOP
	 ENDIF

	 corr=(((phi*gam+1.d0)*gam+1.d0)*gam-a0)/
	1 ((3.d0*phi*gam+2.d0)*gam+1.d0)   !correction de Newton Raphson

c	 WRITE(*,*)'gam,corr,b' ; WRITE(*,2000)gam,corr,b

	 gam=gam-corr	 	!algorithme de Newton Raphson
	 IF(gam == 0.d0)EXIT B1
	 IF(abs(corr/gam) < epsi)EXIT B1		!test arret	  
	ENDDO B1
	  
	IF(gam /= 0.d0)THEN
	 gam1=gam*(gam+1.d0)/b ; gradconv=gradad+gam1
	
	 IF(der)THEN
	  dvstaur=-2.d0*(dvst-1.d0)/taur ; dphitaur=-phi/dvst*dvstaur	
	  dvsltaur=2.d0*dvst*dvstaur ; dksitaur=ksi0*dvsltaur
	  dal4taur=alpha4/ksi*dksitaur ; dbgra=b/gravite
	  dbkra=-b*2.d0/krad ; dbhp=b*3.d0/hp ; dbdel=b/delta
	  dbro=b*2.d0/ro ; dbcp=b*2.d0/cp ; dbtaur=b*dal4taur/alpha4
	  da0gra=dbgra*gradmd ; da0kra=dbkra*gradmd ; da0hp=dbhp*gradmd
	  da0del=dbdel*gradmd ; da0ro=dbro*gradmd ; da0cp=dbcp*gradmd
	  da0rad=b ; da0ad=-b ; da0taur=dbtaur*gradmd
	  
	  gam3=(3.d0*phi*gam+2.d0)*gam+1.d0 ; gam2=(2.d0*gam+1.d0)/gam3
	   	 
	  dgradkra=(gam2*da0kra-gam1*dbkra)/b
	  dgradgra=(gam2*da0gra-gam1*dbgra)/b
	  dgradel= (gam2*da0del-gam1*dbdel)/b
	  dgradhp= (gam2*da0hp -gam1*dbhp )/b
	  dgradro= (gam2*da0ro -gam1*dbro )/b
	  dgradcp= (gam2*da0cp -gam1*dbcp )/b
	  dgradgrad=gam2*da0rad/b ; dgradgad=1.d0+gam2*da0ad/b
	  dgradtaur=(gam2*(da0taur-gam**3*dphitaur) -gam1*dbtaur)/b
	  
	  dgamkra=da0kra/gam3 ; dgamgra=da0gra/gam3 ; dgamdel=da0del/gam3
	  dgamhp=da0hp/gam3 ; dgamro=da0ro/gam3 ; dgamcp=da0cp/gam3
	  dgamtaur=(da0taur-gam**3*dphitaur)/gam3 ; dgamgrad=da0rad/gam3
	  dgamgad=da0ad/gam3
	  
	 ENDIF
	ELSE
	 gradconv=gradad			!gamma=0
	 IF(der)THEN
	  dgradkra=0.d0 ; dgradgra=0.d0 ; dgradel=0.d0 ; dgradhp=0.d0
	  dgradro=0.d0 ; dgradcp=0.d0 ; dgradgrad=0.d0 ; dgradgad=1.d0
	  dgradtaur=0.d0 ; dgamkra=0.d0 ; dgamgra=0.d0 ; dgamdel=0.d0
	  dgamcp=0.d0 ; dgamro=0.d0 ; dgamhp=0.d0 ; dgamtaur=0.d0
	  dgamgrad=0.d0 ; dgamgad=0.d0
	 ENDIF
	ENDIF

c	WRITE(*,2000)gradconv,gradrad,gradad

	RETURN
	
	END SUBROUTINE conv_jmj
