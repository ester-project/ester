
c****************************************************************

	SUBROUTINE conv_a0(krad,gravite,delta,cp,ro,hp,taur,gradrad,
	1 gradad,der,gradconv,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
	2 dgradhp,dgradtaur,dgradgrad,dgradgad,
	3 gam,dgamkra,dgamgra,dgamdel,dgamcp,dgamro,
	4 dgamhp,dgamtaur,dgamgrad,dgamgad)

c	routine private du module mod_conv
	
c	calcul du gradient de température dans le milieu convectif
c	formalisme MLT avec pour longueur de mélange
c	alpha(1-gradad/gradrad) ie, nulle aux limites ZR/ZC

c	formulation de:
c	J. Provost, Département J.D. Cassini, O.C.A., Observatoire de Nice 
c	M. J. Goupil DASGAL Observatoire de Meudon
c	Correction erreur sur b S. Brun 30 08 96

c	Adaptation P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c entrées :
c	krad : conductivité radiative 4 ac T**3 / 3 kap ro
c	gravité : G m / r**2
c	delta : - ( d ln ro / d ln T )p
c	cp : chaleur spécifique
c	ro : densité
c	hp : échelle de hauteur de pression
c	gradrad : gradient radiatif
c	gradad : gradient adiabatique
c	taur : ép. opt. de la bulle
c	der=.true. : calcul des derivees

c	alpha : longueur de mélange

c sortie :
c	gradconv : gradient convectif
c	dgrad* : dérivées
c	gam : le gamma de la convection
c	dgam* : dérivées

c----------------------------------------------------------------

	USE mod_donnees, ONLY : alpha
	USE mod_kind
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: krad, gravite, delta, cp, ro, hp,
	1 taur, gradrad, gradad		
	LOGICAL, INTENT(in) :: der
	REAL (kind=dp), INTENT(out) :: gradconv, dgradkra, dgradgra,
	1 dgradel, dgradcp, dgradro, dgradhp, dgradtaur, dgradgrad,
	2 dgradgad, gam, dgamkra, dgamgra, dgamdel, dgamcp, dgamro,
	3dgamhp, dgamtaur, dgamgrad, dgamgad	
	
	REAL (kind=dp), PARAMETER :: vsal=2.d0/9.d0, phi0=3.d0/2.d0/(3.d0*vsal),
	1ksi0=1.d0/72.d0*(3.d0*vsal)**2
	
	REAL (kind=dp) :: ksi, alpha4, phi, b, dbgra, dbkra, dbhp, dbdel,
	1 dbro, dbcp, gradmd, a0, da0gra, da0kra, da0hp, da0del,
	2 da0ro, da0cp, da0rad, da0ad, gam1, gam2, corr, da0taur,
	3 dbtaur, dal4taur,  dksitaur, dvsltaur, vsal32,
	4 dphitaur, dvst, dvstaur, alp, dalpad, dalprad,
	5 gam3, dalp4ad, dalp4rad, dbad, dbrad
	
	INTEGER :: iter
		
	LOGICAL, SAVE :: init=.TRUE.

c--------------------------------------------------------------------
	
2000	FORMAT(8es10.3)

	IF(init)THEN	!vsal: V/Al
	 init=.FALSE. ; WRITE(2,1)vsal,ksi0,phi0 ; WRITE(*,1)vsal,ksi0,phi0
1	 FORMAT(/,'gradient convectif calcule par MLT : ',/,
	1 'formulation tenant compte de l''ep. op. de la bulle',/,
	2 'alpha=1-gradad/gradrad ie, =0 aux limites ZR/ZC',/,
	3 'V/Al= ',es10.3,' ksi0=',es10.3,' phi0=',es10.3)
	ENDIF

	dvst=1.d0+2.d0/3.d0/vsal/taur**2 ; phi=phi0/dvst ; vsal32=dvst**2
	ksi=ksi0*vsal32 ; alp=alpha*(1.d0-gradad/gradrad)
	alpha4=ksi*alp**4 ; b=alpha4*gravite/krad**2*hp**3*delta*(ro*cp)**2

c	racine réelle de la cubique par newton raphson

c	WRITE(*,*)'ksi,gravite,delta,ro,cp,krad,hp,b,gradrad,gradad'
c	WRITE(*,2000)ksi,gravite,delta,ro,cp,krad,hp,b,gradrad,gradad
	 
	gradmd=gradrad-gradad ; a0=b*gradmd	!B*gradrad-gradad

	gam=ABS(a0/phi)**(1.d0/3.d0)	!initialisation de la racine
	
c	B1 : boucle infinie de convergence, Infinite loop for convergency 	

	iter=0
	B1: DO
	 iter=iter+1
	 IF(iter > 10)THEN
	  WRITE(*,*)'pas de convergence dans conv_a0' ; STOP
	 ELSE 
	  	 
c	  correction de Newton Raphson, Newton Raphson correction

	  corr=(((phi*gam+1.d0)*gam+1.d0)*gam-a0)/
	1 ((3.d0*phi*gam+2.d0)*gam+1.d0)	

c	  WRITE(*,*)'gam,corr,b' ; WRITE(*,2000)gam,corr,b

c	  algorithme de Newton Raphson,  Newton Raphson algorithm
	  
	  gam=gam-corr
	 
	  IF(gam == 0.d0)THEN	 
	   gradconv=gradad			!gamma=0
	   IF(der)THEN	!on a besoin des dérivées, derivatives needed
	    dgradkra=0.d0 ; dgradgra=0.d0 ; dgradel=0.d0 ;dgradhp=0.d0
	    dgradro=0.d0 ; dgradcp=0.d0 ; dgradgrad=0.d0 ; dgradgad=1.d0
	    dgradtaur=0.d0 ; dgamkra=0.d0 ; dgamgra=0.d0 ; dgamdel=0.d0
	    dgamcp=0.d0 ; dgamro=0.d0 ; dgamhp=0.d0 ; dgamtaur=0.d0
	    dgamgrad=0.d0 ; dgamgad=0.d0
	    EXIT B1
	   ENDIF
	  ELSEIF(ABS(corr/gam) < 1.d-10)THEN	!test arret, is it conv.?
	 
	   gam1=gam*(gam+1.d0)/b ; gradconv=gradad+gam1
	
	   IF(der)THEN	!on a besoin des dérivées, derivatives needed
	    dvstaur=-2.d0*(dvst-1.d0)/taur ; dphitaur=-phi/dvst*dvstaur	
	    dvsltaur=2.d0*dvst*dvstaur ; dksitaur=ksi0*dvsltaur
	    dalpad=-alpha/gradrad ; dalprad=alpha*gradad/gradrad**2
	    dalp4ad=4.*alpha4/alp*dalpad ; dalp4rad=4.*alpha4/alp*dalprad
	    dal4taur=alpha4/ksi*dksitaur ; dbgra=b/gravite
	    dbkra=-b*2.d0/krad ; dbhp=b*3.d0/hp ; dbdel=b/delta
	    dbro=b*2.d0/ro ; dbcp=b*2.d0/cp ; dbad=b/alpha4*dalp4ad
	    dbrad=b/alpha4*dalp4rad ; dbtaur=b*dal4taur/alpha4
	    da0gra=dbgra*gradmd ; da0kra=dbkra*gradmd
	    da0hp=dbhp*gradmd ; da0del=dbdel*gradmd ; da0ro=dbro*gradmd
	    da0cp=dbcp*gradmd ; da0rad=b+dbrad*gradmd
	    da0ad=-b+dbad*gradmd ; da0taur=dbtaur*gradmd
	    gam3=(3.d0*phi*gam+2.d0)*gam+1.d0 ; gam2=(2.d0*gam+1.d0)/gam3
	   	 
	    dgradkra=(gam2*da0kra-gam1*dbkra)/b
	    dgradgra=(gam2*da0gra-gam1*dbgra)/b
	    dgradel= (gam2*da0del-gam1*dbdel)/b
	    dgradhp= (gam2*da0hp -gam1*dbhp )/b
	    dgradro= (gam2*da0ro -gam1*dbro )/b
	    dgradcp= (gam2*da0cp -gam1*dbcp )/b
	    dgradgrad=(gam2*da0rad-gam1*dbrad)/b
	    dgradgad= (gam2*da0ad -gam1*dbad )/b+1.d0
	    dgradtaur=(gam2*(da0taur-gam**3*dphitaur) -gam1*dbtaur)/b
	  
	    dgamkra=da0kra/gam3 ; dgamgra=da0gra/gam3 ; dgamdel=da0del/gam3
	    dgamhp=da0hp/gam3 ; dgamro=da0ro/gam3 ; dgamcp=da0cp/gam3
	    dgamtaur=(da0taur-gam**3*dphitaur)/gam3
	    dgamgrad=da0rad/gam3 ; dgamgad=da0ad/gam3
	   ENDIF	!der
	   EXIT B1
	  ENDIF		!gam=0
	 ENDIF  	!iter > 10
	ENDDO B1

c	WRITE(*,2000)gradconv,gradrad,gradad
	
	END SUBROUTINE conv_a0
