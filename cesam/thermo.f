
c*********************************************************************

	SUBROUTINE thermo(pt,p,t,m,l,r,dlpp,xchim,dxchim,
	1    ro,drop,drot,drox,u,dup,dut,dux,grad,dgradpt,dgradp,dgradt,
	2    dgradx,dgradm,dgradl,dgradr,dgradlpp,
	3    gam,dgampt,dgamp,dgamt,dgamx,dgamm,dgaml,dgamr,dgamlpp,
	4    epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
	5    delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	6    gradad,dgradadp,dgradadt,dgradadx,
	7    hp,dhppt,dhpp,dhpt,dhpx,dhpr,dhpm,
	8    gradrad,alfa,beta,gamma1,radiatif)

c	subroutine private du module mod_static
c	calcul de la thermodynamique

c	Auteur: P.Morel, Département J.D. Cassini,
c	O.C.A., Observatoire de Nice
c	CESAM2k

c	MODIF:
c	09 10 96 introduction de la gravité effective
c       dont on tient compte dans Hp
c	10 09 97 correction de dgradr
c	23 10 97 : pression turbulente
c	09 11 99 : correction critère de Ledoux
c	19 11 99 : suppression de nh1, nhe1, nhe2, lamb
c	30 07 00 : introduction F95

c       avec Purb, pt et p sont deux variables indépendantes
c       la pression turbulente affecte l'échelle de hauteur de pression
c       de façon différente que la pression gazeuse,
c       ces différences n'affectent que les dérivées
c       sans pression turbulente la pression totale est la pression gazeuse
c       et leurs dérivées sont identiques	

c       entrées :
c	pt : pression totale Ptot
c	p : pression thermique+radiative Pgas
c	t : température K
c	xchim : composition chimique ATTENTION en 1/mole : *ah pour avoir X
c	dxchim : d xchim/d m pour critère de Ledoux uniquement
c	m : masse/msol
c	l : luminosité/lsol
c	r : rayon / rsol
c	lim: nombre de limites ZR/ZC
c	mstar: masse au temps du calcul, avec perte de masse
c	dlpp : d lnPgas /d lnPtot

c       sortie : (drop : d ro /dp etc..)
c	ro : densité cgs
c	u : énergie interne
c	grad : gradient
c	epsilon : énergie nucléaire cgs
c	kap : opacité Rosseland cm2/gr
c	delta, cp, hp, gradad, gradrad : notations évidentes
c	radiatif=.true. : zone radiative
c	gam : gamma de la convection MLT, sert avec la pression turbulente

c       routines externes :
c	etat, opa, conv, nuc

c-----------------------------------------------------------------------

	USE mod_conv, ONLY: conv
	USE mod_donnees, ONLY : alpha, aradia, clight, cpturb, g,
	1    ihe4, Krot, ledoux, lov_ad, lsol, msol, m_ch, m_rot, nchim, 
	2    nrot, ne, nucleo, ovshts, ovshti, pi, pturb, rsol, tau_max
	USE mod_kind	
	USE mod_etat, ONLY : etat
	USE mod_opa, ONLY : opa
	USE mod_nuc, ONLY: nuc
	USE mod_numerique, ONLY : bsp1dn, no_croiss
	USE mod_variables, ONLY : chim_gram, knotc, knotr, lim, mc, mct,
	1    mrot, mrott, n_ch, n_rot, rota, r_ov, r_zc, wrot

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim, dxchim
	REAL (kind=dp), INTENT(in) :: pt ,p, t, m, l, r, dlpp
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: depsx, epsilon
	REAL (kind=dp), INTENT(out) :: ro, drop, drot, drox, u, dup, dut,
	1    dux, grad, dgradpt, dgradp, dgradt, dgradx, dgradm, dgradl,
	2    dgradr, dgradlpp, gam, dgampt, dgamp, dgamt, dgamx, dgamm,
	3    dgaml, dgamr, dgamlpp, depsp, depst, kap, dkapp, dkapt, dkapx,
	4    delta, deltap, deltat, deltax, cp, dcpp, dcpt, dcpx,
	5    gradad, dgradadp, dgradadt, dgradadx, hp, dhppt, dhpp, dhpt,
	6    dhpx, dhpr, dhpm, gradrad, alfa, beta, gamma1
	LOGICAL, INTENT(out) :: radiatif
	
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION (:,:) :: jac
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION (:) :: xchimm, 
	1    dxchimm, xchiml, dxchiml
	REAL (kind=dp), DIMENSION (nrot) :: dfrot, frot	
	REAL (kind=dp), SAVE :: cte1, cte13, cte2, cte7, cte8, cte9
	REAL (kind=dp) :: hh, be7, b8, n13, o15, f17, taur, dtaurpt,
	1    dtaurp, dtaurt, dtaurx, dtaurr, dtaurm, dgraddpp, dgradlpp0,
	2    dgraddpt, dgraddpx, dgradkra, dgradgra, dgradel, dgradcp,
	3    dgradro, dgradhp, dgradtaur, dgradgrad, dgradgad,
	4    d_grad, graddp, dgradpt0, dgradp0, dgradt0, dgradx0,
	5    dgradm0, dgradr0, dgradl0, krad, dkradp, dkradt, dkradx,
	6    dkapro, depsro, gravite, dgravr, dgravm, rot, w, dw, ldx,
	7    dgamkra, dgamgra, dgamdel, dgamcp, dgamro,
	8    dgamhp, dgamtaur, dgamgrad, dgamgad, m13  
	
	INTEGER :: i=1
	LOGICAL, SAVE :: init=.true.
	LOGICAL	:: ovsht, der
	
c------------------------------------------------------------------------
	
 2000	FORMAT(8es10.3)
 2001	FORMAT(8es11.4)
	
	IF(init)THEN		!initialisations
	   init=.false.	 
	   cte1=4.d0/3.d0*aradia*clight
	   cte2=2.d0/3.d0*rsol	 
	   cte7=4.d0*pi*rsol**2/msol	 
	   cte8=lsol/4.d0/pi/rsol/rsol
c       cte8=cte9*lsol/msol		!calcul de gradrad sans grav. eff.
	   cte9=3.d0/16.d0/pi/aradia/clight/g
	   cte13=g*msol/rsol/rsol

	   WRITE(*,2) ; WRITE(2,2)
 2	   FORMAT(/,'---------THERMO----------',/)	 
	   
	   IF(ledoux)THEN
	      WRITE(*,4) ; WRITE(2,4)	 
 4	      FORMAT('convection: critère de LEDOUX')	 
	   ELSE
	      WRITE(*,5) ; WRITE(2,5)	 
 5	      FORMAT('convection: critère de SCHWARZSCHILD')
	   ENDIF
	   
	   ALLOCATE(xchimm(nchim),dxchimm(nchim),xchiml(nchim),
	1	dxchiml(nchim),jac(nchim,nchim))
	   
	   der=cpturb < 0.d0
	   
	   IF(pturb)THEN
	      WRITE(*,1)abs(cpturb) ; WRITE(2,1)abs(cpturb)	  
 1	      FORMAT('avec pression turbulente C=',es10.3)
	      IF(der)THEN
		 WRITE(*,*)'on tient compte de dln Pgaz/dln Ptot'
		 WRITE(2,*)'on tient compte de dln Pgaz/dln Ptot'	   
	      ELSE
		 WRITE(*,*)'on ne tient pas compte de dln Pgaz/dln Ptot'
		 WRITE(2,*)'on ne tient pas compte de dln Pgaz/dln Ptot'	   
	      ENDIF

	   ELSE
	      WRITE(*,3) ; WRITE(2,3)
 3	      FORMAT('sans pression turbulente')
	   ENDIF
	ENDIF			!initialisations
	
c	composition chimique par gramme

	xchimm=ABS(xchim) ; dxchimm=dxchim ; xchiml=xchim
	CALL chim_gram(xchimm,dxchimm)
	
c	WRITE(*,*)'pt,p,t,xchim,nucleo,xchimm'
c	WRITE(*,2000)pt,p,t,xchim(1),nucleo(1),xchimm(1)
c	WRITE(*,2000)xchimm(1)
c	WRITE(*,2000)(xchimm(i),i=1,nchim)
c	PAUSE'1'

c	vitesse angulaire

	SELECT CASE(Krot)
	CASE(0,1,2)
	   w=wrot ; dw=0.d0
	CASE(3)
	   IF(m == 0.d0)THEN
	      w=rota(1,1) ; dw=0.d0
	   ELSE
	      m13=m**(2.d0/3.d0)
	      CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,
	1	   MIN(m13,mrot(n_rot)),i,frot,dfrot)
	      IF(no_croiss)PRINT*,'Pb. en 1 dans thermo'	
	      w=frot(1) ; dw=dfrot(1)*2.d0/3.d0/m13
	   ENDIF
	END SELECT

c	pression gazeuse (Pgas + Prad)
	
	CALL etat(p,t,xchimm,.TRUE.,
	1    ro,drop,drot,drox,u,dup,dut,dux,
	2    delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3    gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
c	WRITE(*,2001)p,t,xchimm(1),ro,drox ; PAUSE'2'

	drox=drox*nucleo(1)	!par mole
	dux=dux*nucleo(1)
	deltax=deltax*nucleo(1)
	dcpx=dcpx*nucleo(1)
	dgradadx=dgradadx*nucleo(1)
	
c	WRITE(*,*)'ro,drop,drot,drox,u,dup,dut,dux'
c	WRITE(*,2000)ro,drop,drot,drox,u,dup,dut,dux
c	PAUSE'4'

	CALL nuc(t,ro,xchiml,dxchiml,jac,.TRUE.,3, !dxchiml : Vt
	1    epsilon,depst,depsro,depsx,hh,be7,b8,n13,o15,f17)
	depsp=depsro*drop
	depst=depsro*drot+depst
	depsx(1)=depsro*drox+depsx(1)
	
	CALL opa(xchimm,t,ro,kap,dkapt,dkapro,dkapx)
	dkapp=dkapro*drop
	dkapt=dkapro*drot+dkapt	 
	dkapx=dkapro*drox+dkapx*nucleo(1) !par mole

c	WRITE(*,*)'delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
c	1 gradad,dgradadp,dgradadt,dgradadx'
c	WRITE(*,2000)delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
c	1 gradad,dgradadp,dgradadt,dgradadx
c	PAUSE'5'

c	WRITE(*,*)'epsilon,depsp,depst,depsx,p,t,r,l,m'
c	WRITE(*,2000)depsp,depst,p,t,r,l,m,ro
c	WRITE(*,2000)epsilon ; WRITE(*,2000)depsx
c	WRITE(*,*)'kap,dkapp,dkapt,dkapx'
c	WRITE(*,2000)kap,dkapp,dkapt,dkapx
c	PAUSE'6'

c       diffusivite radiative	
	krad=cte1/kap/ro*t**3
	dkradp=krad*(-drop/ro-dkapp/kap)
	dkradt=krad*(-drot/ro-dkapt/kap+3.d0/t)
	dkradx=krad*(-drox/ro-dkapx/kap)

	IF(m*l*r /= 0.d0)THEN	!hors du centre
	   
c       gravité effective
	   gravite=cte13*m/r**2 ; dgravr=-2.*gravite/r ; dgravm= gravite/m	 
	   rot=-cte2*w**2*r	!gravite effective avec rotation
	   gravite=gravite+rot	!remarque de N. Audard
	   dgravr=dgravr-3.*rot/r ; dgravm=dgravm-2.*cte2*w*r*dw
	   
c       échelle de hauteur de pression
	   IF(pturb)THEN
	      hp=pt/gravite/ro	 
	      dhppt=hp/pt ; dhpp=-hp*drop/ro
	   ELSE
	      hp=p/gravite/ro	 
	      dhpp=hp*(1.d0/pt-drop/ro) ; dhppt=dhpp
	   ENDIF	 	 
	   dhpt=-hp*drot/ro ; dhpx=-hp*drox/ro
	   dhpr=-hp*dgravr/gravite ; dhpm=-hp*dgravm/gravite
	   
c       gradient radiatif	 
	   gradrad=cte8*l*hp/r**2/krad/t !gradient radiatif
	   IF(pturb)THEN	 
	      dgradpt=gradrad*dhppt/hp ; dgradp=gradrad*(dhpp/hp-dkradp/krad)
	   ELSE
	      dgradp=gradrad*(dhpp/hp-dkradp/krad) ; dgradpt=dgradp
	   ENDIF	 
	   dgradt=gradrad*(dhpt/hp-dkradt/krad-1.d0/t)
	   dgradx=gradrad*(dhpx/hp-dkradx/krad) ; dgradm=gradrad*dhpm/hp
	   dgradr=gradrad*(dhpr/hp-2.d0/r) ; dgradl=gradrad/l
	   
c       épaisseur optique de la bulle	  	
	   taur=kap*ro*alpha*hp
	   IF(pturb)THEN	 
	      dtaurpt=taur*dhppt/hp ; dtaurp=taur*(dkapp/kap+drop/ro+dhpp/hp)
	   ELSE
	      dtaurp=taur*(dkapp/kap+drop/ro+dhpp/hp) ; dtaurpt=dtaurp
	   ENDIF
	   dtaurt=taur*(dkapt/kap+drot/ro+dhpt/hp)
	   dtaurx=taur*(dkapx/kap+drox/ro+dhpx/hp)
	   dtaurr=taur*dhpr/hp ; dtaurm=taur*dhpm/hp

c       WRITE(*,*)'krad,dkradp,dkradt,dkradx,gravite,dgravr,dgravm,
c	1 hp,dhpp,dhpt,dhpx,dhpr,dhpm,gradrad,dgradp,dgradt,dgradx,
c	2 dgradm,dgradr,dgradl (en dehors du centre)'
c       WRITE(*,2000)krad,dkradp,dkradt,dkradx,gravite,dgravr,dgravm,
c	1 hp,dhpp,dhpt,dhpx,dhpr,dhpm,gradrad,dgradp,dgradt,dgradx,
c	2 dgradm,dgradr,dgradl

	ELSE			!au centre, mais approximativement  
	   gradrad=cte9*kap*epsilon(1)*pt/t**4 !au centre l/m ~ epsilon
	   hp=1.d38 ; gravite=0.d0
	   IF(gradrad /= 0.d0)THEN
	      IF(pturb)THEN	 
		 dgradpt=gradrad/pt ; dgradp=gradrad*(dkapp/kap+depsp/epsilon(1))
	      ELSE	  
		 dgradp=gradrad*(1.d0/p+dkapp/kap+depsp/epsilon(1))
		 dgradpt=dgradp
	      ENDIF
	      dgradt=gradrad*(dkapt/kap+depst/epsilon(1)-4.d0/t)
	      dgradx=gradrad*(dkapx/kap+depsx(1)/epsilon(1))
	   ELSE
	      dgradpt=0.d0 ; dgradp=0.d0 ; dgradt=0.d0 ; dgradx=0.d0
	   ENDIF
	   dgradm=0.d0 ; dgradr=0.d0 ; dgradl=0.d0 ; dgradlpp=0.d0
c       WRITE(*,*)'gradrad,dgradp,dgradt,dgradx (au centre)'
c       WRITE(*,2000)gradrad,dgradp,dgradt,dgradx
	ENDIF
	
c       différence des gradients
	d_grad=gradrad-gradad*dlpp !; PAUSE'7'

c       critère de Ledoux pour la convection, cas ou ro(P,T,X)
c       i.e. mu ne depEND que de X, e.g. EOS MHD ou OPAL

	IF(ledoux .AND. ihe4 > 1 .AND. m > 0.d0)THEN
	   
c       formulation si ro(P,T,X) eg. OPAL, MHD
	   
	   ldx=-cte7*drox*dxchimm(1)*r**2*hp/delta

c       formulation 3-30 Cox & Guili p. 276 p = ro R T / mu + a/3 T^4
c       mu = 4 / (2 + ^ X + Y) ie. tot. ionise
c       ldx=cte7*r**2*hp*ro/(2.d0+6.*xchimm(1)+xchimm(ihe4))*
c	1	beta/(4.d0-3.d0*beta)*(6.d0*dxchimm(1)+dxchimm(ihe4))	
	   
c       IF(d_grad > 0.d0 .AND. ldx >= d_grad)THEN 
c       WRITE(*,2000)dxchimm(1),gradrad,gradad,ldx,d_grad,r
c       PAUSE
c       ENDIF
	   d_grad=d_grad-ldx	 
	ELSE
	   ldx=0.d0
	ENDIF

c	WRITE(*,*)'jpz,lim/xchim(1),d_grad,r,r_zc,r_ov',jpz,lim
c	WRITE(*,2000)xchim(1),d_grad,r,(r_zc(i),r_ov(i),i=1,lim)
c	WRITE(*,2000)d_grad,beta,gradrad,gradad,gradmu,dxchim(1),z,mu

c       test de convection
	radiatif=d_grad <= 0.d0
	
c       overshoot		
	IF(radiatif)THEN	!zone radiative grad_rad < grad_ad
	   
c       mise a 0 de gam et dérivées
	   gam=0.d0 ; dgampt=0.d0 ; dgamp=0.d0 ; dgamt=0.d0 ; dgamx=0.d0
	   dgamm=0.d0 ; dgaml=0.d0 ; dgamr=0.d0 ; dgamlpp=0.d0
	   
c       overshoot 	
	   ovsht=.false.	!.true. : il y aura overshoot
	   DO i=1,lim		!est-on dans une zone d'overshoot?
	      
c       si r_zc a ete fixe et que la limite ZR/ZC s'est legerement deplacee
c       r est considere comme en dehors de la zone d'overshoot
c       on ajoute/retire arbitrairement 0.01
	      
	      IF(r_ov(i) >= 0.d0)THEN
		 IF(ovshts > 0.d0)ovsht=
	1	      ovsht .or. (r <= r_ov(i) .AND. r >= r_zc(i)-0.01d0)
		 IF(ovshti > 0.d0)ovsht=
	1	      ovsht .or. (r >= r_ov(i) .AND. r <= r_zc(i)+0.01d0)
	      ENDIF
	   ENDDO		!lim
c------------------------------------------------------------------------
c	   ovsht=.false.	!pour avoir grad=grad_rad dans zone ovsht
c------------------------------------------------------------------------
c       WRITE(*,*)ovsht,lim
c       WRITE(*,2000)d_grad,r,(r_zc(i),r_ov(i),i=1,lim)
	   IF(ovsht .AND. lov_ad)THEN	!.true. : il y a penetration convective
	      grad=gradad
	      IF(pturb)THEN
		 dgradpt=0.d0 ; dgradp=dgradadp
	      ELSE
		 dgradp=dgradadp ; dgradpt=dgradp
	      ENDIF	  
	      dgradt=dgradadt ; dgradx=dgradadx ; dgradm=0.d0 ; dgradr=0.d0
	      dgradl=0.d0 ; dgradlpp=0.d0
	   ELSE			!pc
	      grad=gradrad	!pas de penetration convective
	      dgradlpp=0.d0	!les dérivées de grad sont celles de gradrad
	   ENDIF		!pc
	   
	ELSE			!zone convective
c       IF(m*l*r == 0.d0 .or. t > 5.d5)THEN !ZC a m=0, ou profonde
	   IF(m*l*r == 0.d0)THEN !ZC a m=0	 
	      grad=gradad
	      IF(pturb)THEN	  
		 dgradpt=0.d0 ; dgradp=dgradadp
	      ELSE	  
		 dgradp=dgradadp ; dgradpt=dgradp
	      ENDIF	  
	      dgradt=dgradadt
	      dgradx=dgradadx ; dgradm=0.d0 ; dgradr=0.d0
	      dgradl=0.d0 ; dgradlpp=0.d0	   
c       WRITE(*,*)'grad,dgradp,dgradt,dgradx (ZC au centre)'
c       WRITE(*,2000)grad,dgradp,dgradt,dgradx
	      
c       zone convective	  
	      
	   ELSE			!zone convective pour m /= 0
	      IF(pturb .AND. der)THEN !avec pression turbulente
		 graddp=gradad*dlpp+ldx	!gradad est * dlnPgaz/dlnPtot
		 dgraddpp=dgradadp*dlpp ; dgraddpt=dgradadt*dlpp
		 dgraddpx=dgradadx*dlpp ; dgradlpp0=gradad	   	   
	      ELSE
		 graddp=gradad+ldx ; dgraddpp=dgradadp ; dgraddpt=dgradadt
		 dgraddpx=dgradadx ; dgradlpp0=0.d0	   	   
	      ENDIF
	      
c       PRINT*,'krad,gravite,delta,cp,ro,hp,taur,gradrad,graddp'
c       WRITE(*,2000)krad,gravite,delta,cp,ro,hp,taur,gradrad,graddp
c       PAUSE'avant conv'	  
	      
	      CALL conv(r,krad,gravite,delta,cp,ro,hp,taur,gradrad,graddp,
	1	   .TRUE.,grad,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
	2	   dgradhp,dgradtaur,dgradgrad,dgradgad,
	3	   gam,dgamkra,dgamgra,dgamdel,dgamcp,dgamro,
	4	   dgamhp,dgamtaur,dgamgrad,dgamgad)
c       PRINT*,'grad,gam'
c       WRITE(*,2000)grad,gam
c       PAUSE'apres conv'
	      
c       pour le gradient de la convection
	      dgradpt0=dgradpt ; dgradt0=dgradt ; dgradr0=dgradr
	      dgradl0=dgradl ; dgradm0=dgradm ; dgradp0=dgradp ; dgradx0=dgradx

	      IF(pturb)THEN	   	   	  
		 dgradpt=dgradhp*dhppt  +dgradtaur*dtaurpt +dgradgrad*dgradpt0
		 dgradp=dgradkra*dkradp + dgradel*deltap + dgradcp *dcpp 
	1	      +dgradro* drop  + dgradhp*dhpp   +dgradtaur*dtaurp
	2	      +dgradgrad*dgradp0+dgradgad*dgraddpp	   
	      ELSE
		 dgradp=dgradkra*dkradp + dgradel*deltap + dgradcp *dcpp 
	1	      +dgradro* drop  + dgradhp*dhpp   +dgradtaur*dtaurp
	2	      +dgradgrad*dgradp0+dgradgad*dgraddpp
		 dgradpt=dgradp
	      ENDIF	  
c       dgradp=dgradkra*dkradp + dgradel*deltap + dgradcp *dcpp 
c	1 +dgradro* drop  + dgradhp*dhpp   +dgradtaur*dtaurp
c	2 +dgradgrad*dgradp0+dgradgad*dgraddpp
	      
	      dgradt=dgradkra*dkradt + dgradel*deltat +  dgradcp*dcpt
	1	   +dgradro*drot  + dgradhp*dhpt   +dgradtaur*dtaurt
	2	   +dgradgrad*dgradt0+dgradgad*dgraddpt
	      
	      dgradx=dgradkra*dkradx +dgradel*deltax  +  dgradcp*dcpx
	1	   +dgradro*drox   +dgradhp*dhpx    +dgradtaur*dtaurx
	2	   +dgradgrad*dgradx0+dgradgad*dgraddpx
	      
	      dgradm=dgradgra*dgravm +dgradhp*dhpm    +dgradtaur*dtaurm  
	1	   +dgradgrad*dgradm0
	      
	      dgradr=dgradgra*dgravr +dgradhp*dhpr    +dgradtaur*dtaurr
	1	   +dgradgrad*dgradr0

	      dgradl=dgradgrad*dgradl0
	      
	      dgradlpp= dgradgad*dgradlpp0 !dérivée/ dln Pgaz/dln Ptot
c       PRINT*,'dgradlpp,dgradgad,dgradlpp0'
c       WRITE(*,2000)dgradlpp,dgradgad,dgradlpp0
	      
c       pour le gamma de la convection
	      IF(pturb)THEN
		 dgampt=dgamhp*dhppt  +dgamtaur*dtaurpt +dgamgrad*dgradpt0
		 dgamp=dgamkra*dkradp + dgamdel*deltap + dgamcp *dcpp 
	1	      +dgamro* drop  + dgamhp*dhpp   +dgamtaur*dtaurp
	2	      +dgamgrad*dgradp0+dgamgad*dgraddpp
	      ELSE	  
		 dgamp=dgamkra*dkradp + dgamdel*deltap + dgamcp *dcpp 
	1	      +dgamro* drop  + dgamhp*dhpp   +dgamtaur*dtaurp
	2	      +dgamgrad*dgradp0+dgamgad*dgraddpp
		 dgampt=dgamp
	      ENDIF	  
c       dgamp=dgamkra*dkradp + dgamdel*deltap + dgamcp *dcpp 
c	1 +dgamro* drop  + dgamhp*dhpp   +dgamtaur*dtaurp
c	2 +dgamgrad*dgradp0+dgamgad*dgraddpp
	      
	      dgamt=dgamkra*dkradt + dgamdel*deltat +  dgamcp*dcpt
	1	   +dgamro*drot + dgamhp*dhpt   +dgamtaur*dtaurt
	2	   +dgamgrad*dgradt0+dgamgad*dgraddpt
	      
	      dgamx=dgamkra*dkradx +dgamdel*deltax  +  dgamcp*dcpx
	1	   +dgamro*drox  +dgamhp*dhpx    +dgamtaur*dtaurx
	2	   +dgamgrad*dgradx0+dgamgad*dgraddpx
	      
	      dgamm=dgamgra*dgravm +dgamhp*dhpm    +dgamtaur*dtaurm  
	1	   +dgamgrad*dgradm0
	      
	      dgamr=dgamgra*dgravr +dgamhp*dhpr    +dgamtaur*dtaurr
	1	   +dgamgrad*dgradr0

	      dgaml=dgamgrad*dgradl0
	      
	      dgamlpp= dgamgad*dgradlpp0 !dérivée/ dln Pgaz/dln Ptot
	      
c       WRITE(*,*)'gam,dgamkra,dgamgra,dgamdel,dgamcp,dgamro,dgamhp,dgamtaur,dgamgrad,dgamgad'	   
c       WRITE(*,2000)gam,dgamkra,dgamgra,dgamdel,dgamcp,dgamro,
c	1 dgamhp,dgamtaur,dgamgrad,dgamgad
c       WRITE(*,*)'gam,dgamp,dgamt,dgamx,dgamm,dgaml,dgamr,dgamlpp'	
c       WRITE(*,2000)gam,dgamp,dgamt,dgamx,dgamm,dgaml,
c	1 dgamr,dgamlpp
c       PAUSE'les gam'
	      
c       WRITE(*,*)'zone convective'
c       WRITE(*,*)'grad,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
c	1 dgradhp,dgradgrad,dgradgad,dgradp,dgradt,dgradx,dgradm,
c	2 dgradr,dgradl'
c       WRITE(*,2000)grad,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
c	1 dgradhp,dgradgrad,dgradgad,dgradp,dgradt,dgradx,dgradm,
c	2 dgradr,dgradl

	   ENDIF		!en m=0	
	ENDIF
c	WRITE(*,*)'jpz,lim/xchim(1),d_grad,r,r_zc,r_ov',jpz,lim
c	WRITE(*,2000)grad,gradad,gradrad,d_grad,r,(r_zc(i),
c	1 r_ov(i),i=1,lim)
c	PAUSE'8'

	RETURN

	END SUBROUTINE thermo
