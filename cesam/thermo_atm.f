
c*********************************************************************

	SUBROUTINE thermo_atm(pt,p,t,xchim,m,l,r,dlpp,
	1 tau,df_tau,d2f_tau,rstar,ro,drop,drot,kap,dkapp,dkapt,gradad,
	2 dgradadp,dgradadt,grad,dgradpt,dgradp,dgradt,dgradr,dgradrs,
	3 dgradm,dgradtau,dgradlpp,gam,dgampt,dgamp,dgamt,dgamr,dgamrs,
	4 dgamm,dgamtau,dgamlpp,hp,dhppt,dhpp,dhpt,dhpr,dhpm,delta,deltap,
	5 deltat,cp,gradrad,alfa,beta,gamma1,radiatif,deriv)

c	routine public du module mod_atm

c	calcul de la thermodynamique pour l'atmosphère avec raccord
c	du gradient
c	on tient éventuellement compte de la pression turbulente
c	Se differencie de thermo par les points suivants :
c	1 - w, L et Xchim sont fixes
c	2 - pas de réactions thermonucléaires
c	3 - loi T(tau) purement radiative ou non
c	4 - dérivées /rstar (rs) et tau

c	19 11 99 : suppression de nh1, nhe1, nhe2, lamb
c	14 07 01 : mise en place de rotation uniforme avec conservation
c		   du moment angulaire
c	16 08 01 : F95

c	Auteur: P. Morel, Département J.D. Cassini,
c	O.C.A., Observatoire de Nice
c	CESAM2k

c entree :
c	pt : Ptot
c	p : Pgaz (Pther+Prad)
c	t : température K
c	xchim : composition chimique ATTENTION /gramme contrairement
c	à thermo
c	m : masse/msol
c	l : luminosité/lsol
c	r : rayon / rsol
c	dlpp : dln Pgaz / dln Ptot
c	tau : profondeur optique
c	df_tau,d2f_tau : fonction  g(tau) et dérivées/tau
c	rad=.true. : loi T(tau) purement radiative
c	rstar : rayon total
c	deriv=.true. : calcul de toutes les dérivées
c	wrot : vitesse angulaire

c sortie : (drop : d ro /dp etc..)
c	ro : densité cgs
c	grad : gradient
c	gam : gamma de la cubique de la MLT
c	dgrad*, dgam* : dérivées
c	kap : opacite Rosseland cm2/gr
c	gradad, gradrad : notations évidentes
c	radiatif=.true. : on est dans une ZR
C
c-----------------------------------------------------------------------

	USE mod_conv, ONLY : conv
	USE mod_donnees, ONLY : alpha, aradia, clight, cpturb, g, Krot,
	1 lsol, msol, nchim, pi, pturb, rsol, tau_max
	USE mod_etat, ONLY : etat
	USE mod_kind
	USE mod_opa, ONLY : opa
	USE mod_variables, ONLY : dim_ch, dim_rot, rota, wrot
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim
	REAL (kind=dp), INTENT(in) :: df_tau, dlpp, d2f_tau, l,
	1 pt, m, r, rstar, tau
	LOGICAL, INTENT(in) :: deriv

	REAL (kind=dp), INTENT(inout) :: p, t
	
	REAL (kind=dp), INTENT(out) :: ro, drop, drot, kap, dkapp, dkapt,
	1 gradad, dgradadp, dgradadt, grad, dgradpt, dgradp, dgradt,
	2 dgradr, dgradrs, dgradm, dgradtau, dgradlpp,
	3 gam, dgampt, dgamp, dgamt, dgamr, dgamrs, dgamm, dgamtau,
	4 dgamlpp, hp, dhppt, dhpp, dhpt, dhpr, dhpm, delta, deltap,
	5 deltat, cp, gradrad, alfa, beta, gamma1	
	LOGICAL, INTENT(out) ::	radiatif

	REAL (kind=dp), SAVE :: cte1, cte2, cte8, cte13, dxdt
			
	REAL (kind=dp) :: u, dup, dut, drox, dux, dkapx, dgradlpp0,
	1 dtaurpt, dtaurp, dtaurt, dtaurr, dtaurm, x, deltax, dcpx,
	3 dgradadx, dgamkra, dgamgra, dgamdel, dgamcp, dgamro,
	4 dgamhp, dgamtaur, dgamgrad, dgamgad, dcpp, dcpt,
	5 krad, dkradp, dkradt, dgradtau0,
	6 gravite, dgravr, dgravm, dgradkra, dgradgra, dgradel,
	7 dgradcp, dgradro, dgradhp, dgradgrad, dgradgad, dkapro,
	8 d_grad, taur, dgradtaur, rot, w, gradrad0, dgradrs0, 
	9 dgradpt0, dgradp0, dgradt0, dgradr0, dgradm0, graddp,
	1 dgraddpp, dgraddpt

cYLD     ajout de tau_rac,grad_int,grad_atm,grad_trans
	REAL (kind=dp) ::grad_int,grad_atm,grad_trans	!YLD
	REAL (kind=dp), SAVE :: tau_rac                !YLD
		
	LOGICAL, SAVE :: init=.true., der
	
c-----------------------------------------------------------------------

2000	FORMAT(8es10.3)
2001	FORMAT(5es15.8)
2002	FORMAT(10es9.2)

c	WRITE(*,2000)pt,p,t,xchim(1),m,l,r
c	PAUSE'entrée thermo_atm'

	IF(init)THEN	!initialisations
	 init=.false.

c         OPEN(unit=30,file='transition.dat',status='unknown')     !YLD  
	 	 
	 cte1=4.d0/3.d0*aradia*clight ; cte2=2.d0/3.d0*rsol 
	 cte8=lsol/4.d0/pi/rsol/rsol ; cte13=g*msol/rsol/rsol
		 
	 WRITE(*,4) ; WRITE(2,4)
4	 FORMAT('convection dans atmosphère: critère de SCHWARZSCHILD')
	 
c--------------------------------------------------------------------	 
cYLD	 dxdt=1.d0/(tau_max-1.d0)  !pente pour le raccord du gradient
         tau_rac=1.D0                                              !YLD
	 dxdt=1.d0/(tau_max-tau_rac)	!pente pour le raccord du gradient     !YLD
	 WRITE(2,*) 'tau_rac au sommet zone de transition :', tau_rac       !YLD
c----------------------------------------------------------------------
	 der=cpturb < 0.d0

	 IF(pturb)THEN
	  WRITE(*,1) ; WRITE(2,1)	  
1	  FORMAT('avec press. turb. dans atmosphère')
	  IF(.NOT.rad)THEN
	   WRITE(*,3) ; WRITE(2,3)
3	   FORMAT('ATTENTION press. turb. avec loi T(tau) non tot. rad')
	  ENDIF
	
	 ELSE
	  WRITE(*,2) ; WRITE(2,2)	  
2	  FORMAT('sans press. turb. dans atmosphère') 
	 ENDIF		 
	ENDIF		!initialisation

c	pression gazeuse (Pgas + Prad)
	
c	PRINT*,deriv ; WRITE(*,2000)p,t ; WRITE(*,2000)xchim(1:nchim)	
c	PAUSE'avant etat'
	CALL etat(p,t,xchim,deriv,
	1 ro,drop,drot,drox,u,dup,dut,dux,
	2 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
c	CALL pause('apres etat')

c	Omega

	SELECT CASE(Krot)
	CASE(0,1,2)
	 w=wrot
	CASE(3)
	 w=rota(1,dim_rot)	 
	END SELECT

	IF(deriv)THEN	

	 CALL opa(xchim,t,ro,kap,dkapt,dkapro,dkapx)
c	 CALL pause('apres opa')
	 dkapp=dkapro*drop ; dkapt=dkapro*drot+dkapt
	 
c	 diffusivité radiative
	
	 krad=cte1/kap/ro*t**3 ; dkradp=krad*(-drop/ro-dkapp/kap)
	 dkradt=krad*(-drot/ro-dkapt/kap+3.d0/t)

c	 gravité effective

	 gravite=cte13*m/r**2 ; dgravr=-2.*gravite/r ; dgravm= gravite/m
	 
	 rot=-cte2*w**2*r	!gravite effective avec rotation
	 gravite=gravite+rot	!rotation, remarque de N. Audard
	 dgravr=dgravr-3.*rot/r
	
c	 échelle de hauteur de pression
	 
	 hp=pt/gravite/ro		!echelle de hauteur de pression
	 dhppt=hp/pt ; dhpp=-hp*drop/ro ; dhpt=-hp*drot/ro
	 dhpr=-hp*dgravr/gravite ; dhpm=-hp*dgravm/gravite
	
c	 gradient radiatif	 
	 
	 gradrad=cte8*l*hp/r**2/krad/t	!gradient radiatif
	 dgradpt=gradrad*dhppt/hp
	 dgradp=gradrad*(dhpp/hp-dkradp/krad)	 
	 dgradt=gradrad*(dhpt/hp-dkradt/krad-1.d0/t)
	 dgradm=gradrad*dhpm/hp ; dgradr=gradrad*(dhpr/hp-2.d0/r)
	 
c	 épaisseur optique de la bulle	 
	 	
	 taur=kap*ro*alpha*hp		!epaisseur optique de la bulle
	 dtaurpt=taur*dhppt/hp
	 dtaurp=taur*(dkapp/kap+drop/ro+dhpp/hp)	 	 
	 dtaurt=taur*(dkapt/kap+drot/ro+dhpt/hp)
	 dtaurr=taur*dhpr/hp ; dtaurm=taur*dhpm/hp
	 
	ELSE
	 CALL opa(xchim,t,ro,kap,dkapt,dkapro,dkapx)
	 krad=cte1/kap/ro*t**3	
	 
	 gravite=cte13*m/r**2
	 rot=-cte2*w**2*r	!gravité effective avec rotation
	 gravite=gravite+rot	!rotation, remarque de N. Audard
	 hp=pt/gravite/ro		!echelle de hauteur de pression
	 gradrad=cte8*l*hp/r**2/krad/t	!gradient radiatif
	 taur=kap*ro*alpha*hp		!epaisseur optique de la bulle
	ENDIF		!deriv
	
c	test de convection
c	gradrad0 est la valeur du gradient dérivée de la loi T(tau)
c	-- si la loi T(tau) n'est pas purement radiative, on fait un
c		raccord avec la valeur du gradient issue de la convection, 
c	-- si la loi T(tau) est purement radiative, on utilise gradrad0
c		pour valeur du gradient radiatif dans le calcul de la
c		valeur du gradient avec convection
c	le gradient radiatif, dérivé de la loi T(tau) (gradrad0),
c	ne suppose pas
c	l'approximation de diffusion et tient compte, en principe,
c	d'un transfert radiatif complet, il diffère donc, un peu,
c	du gradient
c	radiatif (gradrad) l*hp/r**2/krad/t de la structure interne
c	c'est gradrad qui est, arbitrairement, utilisé pour le critère
c	de convection

	gradrad0=gradrad*df_tau*(r/rstar)**2
	dgradpt0=dgradpt*df_tau*(r/rstar)**2
	dgradp0 =dgradp *df_tau*(r/rstar)**2	
	dgradt0= dgradt *df_tau*(r/rstar)**2
	dgradr0= dgradr *df_tau*(r/rstar)**2+2.*gradrad0/r
	dgradrs0=-2.d0*gradrad0/rstar		!dérivée/Rstar
	dgradm0= dgradm *df_tau*(r/rstar)**2
	dgradtau0=gradrad*d2f_tau*(r/rstar)**2
	
c	différence des gradients
	
	IF(rad)THEN		!loi T(tau) purement radiative	
	 d_grad=gradrad0-gradad*dlpp	!Schwarzschild pour la convection
	ELSE
	 d_grad=gradrad-gradad*dlpp
	ENDIF	
	
c	PRINT*,'d_grad,gradrad,gradad'
c	WRITE(*,2000)d_grad,gradrad,gradad
c	PRINT*,'krad,gravite,delta,cp,ro,hp,taur'
c	WRITE(*,2000)krad,gravite,delta,cp,ro,hp,taur ; PAUSE'krad'
	
	radiatif=d_grad <= 0.d0	!zone radiative grad_rad < grad_ad
	
	IF(radiatif)THEN	
	 grad=gradrad0			!gradient=gradrad * df/dtau
	 dgradpt=dgradpt0 ; dgradp=dgradp0 ; dgradt=dgradt0
	 dgradr=dgradr0 ; dgradrs=dgradrs0 ; dgradm=dgradm0
	 dgradtau=dgradtau0 ; dgradlpp=0.d0
	 
c	 mise a 0 de gam et dérivées
	
	 gam=0.d0 ; dgampt=0.d0 ; dgamp=0.d0 ; dgamt=0.d0 ; dgamm=0.d0
	 dgamr=0.d0 ; dgamrs=0.d0 ; dgamtau=0.d0 ; dgamlpp=0.d0
 	 	 	 	  
	ELSE

	 IF(pturb .AND. der)THEN		!avec pression turbulente
	  graddp=gradad*dlpp	!gradad est * dlnPgaz/dlnPtot
	  dgraddpp=dgradadp*dlpp ; dgraddpt=dgradadt*dlpp
	  dgradlpp0=gradad
	 ELSE
	  graddp=gradad ; dgraddpp=dgradadp ; dgraddpt=dgradadt
	  dgradlpp0=0.d0	  
	 ENDIF
	 
c	 PAUSE'avant conv'
	 
	 IF(rad)THEN		!loi T(tau) purement radiative
	 	 
c	  dans une ZC avec une loi T(tau) purement radiative gradrad
c	  est modifié par la correction d'atmosphère
c	  on prend gradrad0 comme gradient radiatif pour le
c	  calcul du gradient

c	  WRITE(*,*)'krad,gravité,delta,cp,ro,hp,taur,gradrad0,graddp'
c	  WRITE(*,2000)krad,gravite,delta,cp,ro,hp,taur,gradrad0,graddp
c	  PAUSE'avant conv'
	  	  
	  CALL conv(r,krad,gravite,delta,cp,ro,hp,taur,gradrad0,graddp,
	1 deriv,grad,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
	2 dgradhp,dgradtaur,dgradgrad,dgradgad,
	3 gam,dgamkra,dgamgra,dgamdel,dgamcp,dgamro,
	4 dgamhp,dgamtaur,dgamgrad,dgamgad)
c	  WRITE(*,*)'grad,gam' ; WRITE(*,2000)grad,gam
c	  CALL pause('après conv')
	
	  IF(deriv)THEN
	  
c	   dérivées du gradient	  
	 
	   dgradpt=dgradhp*dhppt  +dgradtaur*dtaurpt +dgradgrad*dgradpt0
	
	   dgradp=dgradkra*dkradp  +dgradel*deltap+dgradcp *dcpp 
	1  +dgradro* drop  +dgradhp*dhpp  +dgradtaur*dtaurp
	2  +dgradgrad*dgradp0+dgradgad*dgraddpp
	  
	   dgradt=dgradkra*dkradt   +dgradel*deltat+dgradcp*dcpt
	1  +dgradro*drot   +dgradhp*dhpt  +dgradtaur*dtaurt
	2  +dgradgrad*dgradt0+dgradgad*dgraddpt
	
	   dgradm=dgradgra*dgravm+dgradhp*dhpm+dgradtaur*dtaurm  
	1  +dgradgrad*dgradm0
	
	   dgradr=dgradgra*dgravr+dgradhp*dhpr+dgradtaur*dtaurr
	1  +dgradgrad*dgradr0
	
	   dgradrs=dgradgrad*dgradrs0
	  
	   dgradtau=dgradgrad*dgradtau0		  	

	   dgradlpp= dgradgad*dgradlpp0 	!dérivée/ dln Pgaz/dln Ptot
	   
c	   dérivées du gamma de la cubique

	   dgampt=dgamhp*dhppt  +dgamtaur*dtaurpt +dgamgrad*dgradpt0
	
	   dgamp=dgamkra*dkradp  +dgamdel*deltap+dgamcp *dcpp 
	1   +dgamro* drop  +dgamhp*dhpp  +dgamtaur*dtaurp
	2   +dgamgrad*dgradp0+dgamgad*dgraddpp
	   
	   dgamt=dgamkra*dkradt   +dgamdel*deltat+dgamcp*dcpt
	1  +dgamro*drot   +dgamhp*dhpt  +dgamtaur*dtaurt
	2  +dgamgrad*dgradt0+dgamgad*dgraddpt
	
	   dgamm=dgamgra*dgravm+dgamhp*dhpm+dgamtaur*dtaurm  
	1  +dgamgrad*dgradm0
	
	   dgamr=dgamgra*dgravr+dgamhp*dhpr+dgamtaur*dtaurr
	1  +dgamgrad*dgradr0
	
	   dgamrs=dgamgrad*dgradrs0
	  
	   dgamtau=dgamgrad*dgradtau0		  	

	   dgamlpp= dgamgad*dgradlpp0		!dérivée/ dln Pgaz/dln Ptot
	   
c	   WRITE(*,*)'gam,dgamkra,dgamgra,dgamdel,dgamcp,dgamro,
c	1  dgamhp,dgamtaur,dgamgrad,dgamgad'	   
c	   WRITE(*,2000)gam,dgamkra,dgamgra,dgamdel,dgamcp,dgamro,
c	1  dgamhp,dgamtaur,dgamgrad,dgamgad
c	   WRITE(*,*)'gam,dgamp,dgamv,dgamt,dgamr,dgamrs,dgamm,
c	1  dgamtau,dgamlpp'	   
c	   WRITE(*,2000)gam,dgamp,dgamv,dgamt,dgamr,dgamrs,dgamm,
c	1  dgamtau,dgamlpp	   
c	   PAUSE'les gam'

	  ENDIF
	    
	 ELSE	!loi T(tau) non purement radiative
	 
c	  avec une loi non purement radiative on pondère par
c	  exp(tau-tau_max)
c	  le gradient de la loi T(tau) ie. gradrad*df_tau=gradrad0
c	  et le gradient de la MLT calcule avec gradrad

	  CALL conv(r,krad,gravite,delta,cp,ro,hp,taur,gradrad,
	1 graddp,deriv,grad,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
	2 dgradhp,dgradtaur,dgradgrad,dgradgad,
	3 gam,dgamkra,dgamgra,dgamdel,dgamcp,dgamro,
	4 dgamhp,dgamtaur,dgamgrad,dgamgad)
	
cYLD	  x=max((tau-1.d0)*dxdt,0.d0)
	  x=max((tau-tau_rac)*dxdt,0.d0)   !YLD
	  	
	  IF(deriv)THEN
	   dgradpt=dgradhp*dhppt  +dgradtaur*dtaurpt +dgradgrad*dgradpt0	  
	   dgradpt=dgradpt0*(1.d0-x)+dgradpt*x
	   	  
	   dgradp=dgradkra*dkradp + dgradel*deltap+ dgradcp *dcpp 
	1  +dgradro* drop + dgradhp*dhpp  +dgradtaur*dtaurp
	2  +dgradgrad*dgradp+dgradgad*dgraddpp
	   dgradp=dgradp0*(1.d0-x)+dgradp*x	
	   
	   dgradt=dgradkra*dkradt  + dgradel*deltat+  dgradcp*dcpt
	1  +dgradro*drot  + dgradhp*dhpt  +dgradtaur*dtaurt
	2  +dgradgrad*dgradt+dgradgad*dgraddpt
	   dgradt=dgradt0*(1.d0-x)+dgradt*x	
	
	   dgradm=dgradgra*dgravm+dgradhp*dhpm +dgradtaur*dtaurm  
	1  +dgradgrad*dgradm
	   dgradm=dgradm0*(1.d0-x)+dgradm*x	
	
	   dgradr=dgradgra*dgravr +dgradhp*dhpr+dgradtaur*dtaurr
	1  +dgradgrad*dgradr
	   dgradr=dgradr0*(1.d0-x)+dgradr*x	
	
	   dgradrs=dgradgrad*dgradrs0
	   dgradrs=dgradrs*(1.d0-x)		   
	  
	   dgradlpp= dgradgad*dgradlpp0		!dérivée/ dln Pgaz/dln Ptot
	   dgradlpp=dgradlpp*x
	   
	   dgradtau=dgradtau0*(1.d0-x)+dxdt*(grad-gradrad0)	   	   
	  
c	   WRITE(*,2001)krad,gravite,delta,cp,ro,hp,taur,gradrad,graddp
c	1  ,grad

c	   gam ne servant qu'avec la pression turbulente
c	   il ne peut y avoir raccord correct avec la loi T(tau)
c	   ses dérivées n'ont donc pas d'intérêt

	  dgampt=0.d0 ;  
	  dgamp=0.d0 ; dgamt=0.d0 ; dgamr=0.d0 ; dgamrs=0.d0 ; dgamm=0.d0
	  dgamtau=0.d0 ; dgamlpp=0.d0
	  	  
	  ENDIF		!deriv
	  
c---------------------------------------------------------------------------
cYLD       quantités pour atmosphère : fichier transition.dat
	  grad_int=grad     !YLD
	  grad_atm=gradrad0  !YLD
	  grad=gradrad0*(1.d0-x)+grad*x		!ne pas deplacer
	  grad_trans=grad    !YLD
	  write(30,2001) tau, grad_trans,grad_int,grad_atm,t,p,krad,hp  !YLD
	    
c	  WRITE(*,2000)x,dxdt,dgradtau0,gradrad0,grad,tau,d2f_tau,dgradtau
c	  WRITE(*,2001)tau,gradrad0,grad ; WRITE(*,2001)tau,gradrad0,
c	1 dgradtau0
	 ENDIF
	  
	ENDIF
	
c	WRITE(*,2000)tau,d_grad ; WRITE(*,2000)p,pg,t,m,l,r,rstar,grad
	
	RETURN

	END SUBROUTINE thermo_atm
