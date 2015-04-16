
c*************************************************************

	SUBROUTINE coeff_rota3(dt,nui,y,frl,coeff)

c routine private du module mod_evol

c calcul des coefficients pour la rotation
c formalisme de Talon & Zahn 1997, Krot=3

c entrées :
c	dt : pas temporel
c	nu : abscisse en m^2/3
c	y : variables de la rotation

c sorties :
c	frl : coefficients de rotation et dérivées
c	deff, dh, dv : coefficients de diffusion du moment cinétique
c	coeff : OPTIONAL, sortie des coefficients pour dessins

c	Auteur: P. Morel, Département Cassiopée, O.C.A., CESAM2k

c NOTATIONS (hélas incohérentes) pour les développements sur B-splines
c	n_ch : nombre VARIABLE de points élément de mod_variables
c	nch : nombre FIXE de fonctions élément de mod_donnees
c	m_ch : ordre FIXE des splines élément de mod_donnees
c	mch(n_ch) : abscisses VARIABLES élément de mod_variables

c Dans une ZC sont seuls utiles
c frl(1, 2, 3, 4, 8, 9), Deff, Dv, Dh

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : amu, aradia, clight, diffusion, d_conv,
	1 en_masse, g, ihe4, lsol, msol, mu_saha, m_ch, m_ptm, ne, nchim,
	2 nucleo, nrl, nrot, pi, rsol, secon6, zi, ord_qs, ord_rot
	USE mod_etat, ONLY: etat, saha
	USE mod_kind
	USE mod_nuc, ONLY : l_planet, mzc_ext, nuc, nuzc_ext, planetoides
	USE mod_numerique, ONLY : bsp1dn, no_croiss
	USE mod_opa, ONLY : opa
	USE mod_variables, ONLY : bp, bp_t, chim, chim_t, chim_gram,
	1 inter, knot, knotc, knotc_t, knotr, knotr_t, knot_ptm, knot_t,
	2 mc, mct, mc_t, mct_t, mrot, mrot_t, mrott, mrott_t, m23, m23_t,
	3 n_ch, n_ch_t, n_ptm, n_qs, n_qs_t, n_rot, n_rot_t,
	4 old_ptm, q, qt, qt_t, q_t, rota, rota_t, r2, r2_t, x_ptm, xt_ptm

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(nrot,0:1) :: y
	REAL (kind=dp), INTENT(in) :: dt, nui
	REAL (kind=dp), INTENT(out), DIMENSION(nrl) :: frl
	REAL (kind=dp), INTENT(out), OPTIONAL, DIMENSION(ncoeff+nchim) :: coeff
	REAL (kind=dp), DIMENSION(2,0:1) :: ys
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: ioni
	REAL (kind=dp), DIMENSION(0,0) :: jac
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: cte_0
	REAL (kind=dp), DIMENSION(nchim) :: depsx, dxchim, dxchim_t, dxchimg,
	1 dxchimg_t, xchim, xchim_t, xchimg, xchimg_t, z_bar
	REAL (kind=dp), DIMENSION(ne) :: dfqs, fqs

	REAL (kind=dp), DIMENSION(nrot) :: dfrot, y_t
	REAL (kind=dp), SAVE, DIMENSION(15) :: cts_0
	REAL (kind=dp), DIMENSION(15) :: cts
	REAL (kind=dp), DIMENSION(5) :: epsilon
	REAL (kind=dp), DIMENSION(2) :: mw_planet	
	REAL (kind=dp), DIMENSION(1) :: dfm, fm
	REAL (kind=dp), DIMENSION(0) :: dcomp

	REAL (kind=dp), SAVE ::	cte1_0, cte2_0, cte3_0, cte6_0, ln13
	REAL (kind=dp) :: alfa, Abar, beta, be7e, b8e, chi, chi_mu, chi_t,
	1 cp, dcpp, dcpt, dcpx, deff, delta,
	2 deltap, deltat, deltax, depsro, depst, dgrad, dgradadp, dgradadt,
	3 dgradadx, dh, dkapr, dkapt, dkapx, dlldlm, dlnmu, dlnro_nu,
	4 dlnP_nu, dlnro_t, drop, drot, drox, dup, dut, dux, dv, eps_mu,
	5 eps_nuc, eps_t, eta, f17e, gamma1, grad, gradad, grad_mu_shp,
	6 grand_k, grav, hhe, hp, ht, kap, lnP, lnP_t, ln_lambda, lum,
	7 mr, mu, mu_t, mw_dot, nel, nu, nu12, n13e, nu32, nu52,
	8 n2mu, n2t, o15e, p, phi_e, p_t, r, ro, ro_t, rr_t, t, t_t, u, vis_cin,
	9 Zbar

	INTEGER, SAVE :: l=1

	LOGICAL, SAVE :: init=.TRUE.
	LOGICAL :: zc

c------------------------------------------------------------------------

2000	FORMAT(8es10.3)

c définition des constantes
	IF(init)THEN
	 init=.FALSE.

c pour le calcul de grad_mu, chi, grav
	 cte1_0=8.d0*pi*rsol**2/3.d0/msol
	 cte2_0=4.d0*aradia*clight/3.d0
	 cte3_0=g*msol/rsol**2
	 cte6_0=2.d-15
	 ln13=LOG(1.3d4)

c coefficients
	 ALLOCATE(cte_0(nrl))
	 cte_0(1)=1.d0/secon6
	 cte_0(2)=1.d0
	 cte_0(3)=8.d0*pi*rsol**2/15.d0/msol
	 cte_0(4)=64.d0*pi**2*rsol**4/9.d0/msol**2
	 cte_0(5)=4.d0*rsol**3/3.d0/pi/g**2/msol
	 cte_0(6)=msol/lsol
	 cte_0(7)=1.d0
	 cte_0(8)=1.d0/secon6
	 cte_0(9)=40.d0*pi*rsol**2/3.d0/msol
	 cte_0(10)=msol/lsol/secon6
	 cte_0(11:12)=1.d0
	 cte_0(13)=8.d0*pi*rsol**2/3.d0/msol
	 cte_0(14)=16.d0*pi*rsol**2/9.d0/msol
	 cte_0(15)=1.d0		!0.d0
	 cte_0(16:18)=-2.d0/3.d0
	 cte_0(19)=16.d0*pi*rsol**6/9.d0/g/msol**2
	 cte_0(20:21)=1.d0
	 cte_0(22)=1.d0/secon6
	 cte_0(23:24)=1.d0

	 cts_0(1)=1.d0/secon6
	 cts_0(2)=msol/2.d0/pi/rsol**3
	 cts_0(3)=msol/lsol
	 cts_0(4)=2.d0
	 cts_0(5:6)=3.d0*msol/2.d0/pi/rsol**4
	 cts_0(7)=msol/lsol
	 cts_0(8:9)=1.d0
	 cts_0(10)=msol/lsol
	 cts_0(11)=1.d0/secon6
	 cts_0(12)=msol/rsol**2		!M^dot_Omega en Msol
	 cts_0(13)=6.d0/rsol**2
	 cts_0(14)=8.d0*rsol**3/3.d0/g/msol
	 cts_0(15)=msol/lsol/secon6

	 IF(mu_saha)ALLOCATE(ioni(0:NINT(MAXVAL(zi)),nchim))

	ENDIF

c mises à 0
	cts=0.d0 ; frl=0.d0 ; mw_dot=0.d0
	IF(PRESENT(coeff))coeff=0.d0

c nu12 etc...
	nu=MAX(nui,nu_min) ; nu12=SQRT(nu) ; nu32=nu12**3 ; nu52=nu32*nu

c est-on dans une ZC?
	zc=lmix(nu)

c~~~~~~~~~~~~~~~ au temps t ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	IF(en_masse)THEN
	 mr=nu
	ELSE
	 mr=nu32
	ENDIF
	mr=MAX(x_ptm(1),MIN(mr,x_ptm(n_ptm)))	!perte de masse
	CALL bsp1dn(1,old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm,.TRUE.,
	1 mr,l,fm,dfm)	!masse au temps t fm(1)
	IF(no_croiss)PRINT*,'Pb en 2 dans coeff_rota3'
	IF(.NOT.en_masse)fm(1)=fm(1)**(2.d0/3.d0)
	mr=MAX(mrot_t(1),MIN(fm(1),mrot_t(n_rot_t)))
	CALL bsp1dn(nrot,rota_t,mrot_t,mrott_t,n_rot_t,ord_rot,
	1 knotr_t,.TRUE.,mr,l,y_t,dfrot)
	IF(no_croiss)PRINT*,'Pb en 3 dans coeff_rota3'
	mr=MAX(m23_t(1),MIN(fm(1),m23_t(n_qs_t)))
 	CALL inter('m23',bp_t,q_t,qt_t,n_qs_t,knot_t,mr,fqs,dfqs,r2_t,m23_t)
	IF(no_croiss)PRINT*,'Pb. en 4 dans coeff_rota3'
	lnP_t=fqs(1) ; p_t=EXP(lnP_t) ; t_t=EXP(fqs(2))
	IF(en_masse)THEN
	 rr_t=ABS(fqs(3))
	ELSE
	 rr_t=fqs(3)**2
	ENDIF

c composition chimique
	mr=MAX(mc_t(1),MIN(nu,mc_t(n_ch_t)))
	CALL bsp1dn(nchim,chim_t,mc_t,mct_t,n_ch_t,m_ch,knotc_t,.TRUE.,mr,l,
	1 xchim_t,dxchim_t)
	IF(no_croiss)PRINT*,'Pb en 4 dans coeff_rota3'

c composition chimique par gramme, équation d'état
	xchimg_t=ABS(xchim_t) ; dxchimg_t=dxchim_t
	CALL chim_gram(xchimg_t,dxchimg_t)
	CALL etat(p_t,t_t,xchimg_t,.FALSE.,ro_t,drop,drot,drox,u,dup,dut,dux,
	1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)

c poids moléculaire moyen, appel à saha pour les charges moyennes
	IF(mu_saha)THEN
	 CALL saha(xchim_t,t_t,ro_t,ioni,z_bar,nel,eta)
	 mu_t=1.d0/DOT_PRODUCT((1.d0+z_bar),xchim_t)
	ELSE
	 mu_t=16.d0/(23.d0*xchimg_t(1)+3.d0*xchimg_t(ihe4)+9.d0)
	ENDIF

c ~~~~~~~~~~~~~~~~au temps t+dt	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c ATTENTION, dérivées / m^2/3 ou m effectuées dans inter
	CALL inter('m23',bp,q,qt,n_qs,knot,nu,fqs,dfqs,r2,m23)
	IF(no_croiss)PRINT*,'Pb en 0 dans coeff_rota3'
	lnP=fqs(1) ; p=EXP(fqs(1)) ; t=EXP(fqs(2))
	IF(en_masse)THEN
	 r=SQRT(ABS(fqs(3))) ; lum=SQRT(ABS(fqs(4)))**3 ; dlnP_nu=dfqs(1)
	ELSE
	 r=ABS(fqs(3)) ; lum=fqs(4)
	ENDIF
	dlldlm=dfqs(4)*fqs(5)/fqs(4)

c grad = d ln T / d ln P
	grad=dfqs(2)/dfqs(1)

c composition chimique
	mr=MAX(mc(1),MIN(nu,mc(n_ch)))
	CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,knotc,.TRUE.,mr,l,xchim,dxchim)
	IF(no_croiss)PRINT*,'Pb en 1 dans coeff_rota3'

c composition chimique par gramme, équation d'état
	xchimg=ABS(xchim) ; dxchimg=dxchim
	CALL chim_gram(xchimg,dxchimg)
	CALL etat(p,t,xchimg,.FALSE.,ro,drop,drot,drox,u,dup,dut,dux,
	1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)

c hp, ht~~~~~~~~~~~~~
	hp=r**2*p/ro/nu32/cte3_0
	ht=hp/grad

c gravité
	grav=cte3_0*nu32/r**2

c appel à saha pour les charges moyennes, poids moléculaire moyen
	IF(mu_saha)THEN
	 CALL saha(xchim,t,ro,ioni,z_bar,nel,eta)
	 mu=1.d0/DOT_PRODUCT((1.d0+z_bar),xchim)
	ELSE
	 mu=16.d0/(23.d0*xchimg(1)+3.d0*xchimg(ihe4)+9.d0)
	ENDIF

c phi_e : estimation de phi= dln ro/dln mu
	phi_e=delta*beta/(4.d0-3.d0*beta)
	
c d Ln mu
	 IF(mu_saha)THEN
	  dlnmu=-mu*DOT_PRODUCT((1.d0+z_bar),dxchim)
	 ELSE
	  dlnmu=-mu/16.d0*(23.d0*dxchimg(1)+3.d0*dxchimg(ihe4))
	 ENDIF

c d ln ro / dt
	 dlnro_t=(alfa-delta*grad)*(lnP-lnP_t)+phi_e*LOG(mu/mu_t)

c d ln ro / d nu
	 dlnro_nu=(alfa-delta*grad)*dlnP_nu+phi_e*dlnmu

c Grad µ / Hp
	 grad_mu_shp=-cte1_0*r**2*ro/nu12*dlnmu

c appel à nuc pour les epsilon
	 CALL nuc(t,ro,xchim,dcomp,jac,.TRUE.,3,
	1 epsilon,depst,depsro,depsx,hhe,be7e,b8e,n13e,o15e,f17e)
	 eps_nuc=epsilon(1)
	 IF(eps_nuc > 0.d0)THEN
	  eps_t=t/eps_nuc*depst
	  eps_mu=phi_e*ro/eps_nuc*depsro
	 ELSE
	  eps_mu=0.d0 ; eps_t=0.d0
	 ENDIF

c estimation de  vis_cin, la viscosité cinématique, Lang 3.15, 3.25 et 3.28
	 Abar=DOT_PRODUCT(xchim,nucleo)/SUM(xchim)
	 IF(mu_saha)THEN
	  Zbar=DOT_PRODUCT(z_bar,xchim)/SUM(xchim)
	  nel=DOT_PRODUCT(z_bar,xchim)*ro/amu
	 ELSE
	  Zbar=DOT_PRODUCT(zi,xchim)/SUM(xchim)
	  nel=DOT_PRODUCT(zi,xchim)*ro/amu
	 ENDIF
	 ln_lambda=ln13+1.5d0*LOG(t)-LOG(nel)/2.d0
	 vis_cin=cte6_0*SQRT(t)**5*SQRT(Abar)/ro/Zbar**4/ln_lambda

c opacités
	 CALL opa(xchimg,t,ro,kap,dkapt,dkapr,dkapx)

c conductivité radiative chi, dérivées et diffusivité thermique K
	 chi=cte2_0*t**3/ro/kap
	 chi_t=3.d0-t/kap*dkapt+delta*(1.d0+ro/kap*dkapr)
	 chi_mu=-phi_e*(1.d0+ro/kap*dkapr)
	 grand_k=chi/ro/cp

c vaïssälä
	 dgrad=gradad-grad
	 n2t=grav*delta/hp*ABS(dgrad)
	 dgrad=ABS(dgrad+phi_e*grad_mu_shp*hp/delta)
	 n2mu=grav*phi_e*grad_mu_shp

c coefficients de diffusion
	IF(zc)THEN
	 deff=d_conv ; dh=d_conv ; dv=d_conv

c perte/gain de moment cinétique pour nu > nuzc_ext dus aux 
c différentes sources externes de moment cinétique:
c 1 - pertw: perte/gain de masse dû au vent
c 2 - mw_planet(1) apport des planétoïdes de vitesse angulaire non nulle
c 3 - correction affectée du signe -  dmw_planet(2) addition fictive
	 IF(nu >= nuzc_ext)THEN
	  IF(lw_perte)CALL pertw(r,mw_dot)
	  IF(l_planet)THEN
	   CALL planetoides(mw_planet=mw_planet)
	   mw_dot=mw_dot+mw_planet(1)-mw_planet(2)*r**2*y(1,0)*dt
	  ENDIF
	 ENDIF

	ELSE

c coefficients de diffusion du moment cinétique y --> ys
c on utilise Omega et U et leurs dérivées au pas temporel précédent
	 IF(PRESENT(coeff))THEN
	  ys(1:2,:)=y(1:2,:)
	 ELSE
c	  ys(1:2,0)=y_t(1:2) ; ys(1:2,1)=dfrot(1:2)
	  ys(1:2,:)=y(1:2,:)		!totalement implicite
	 ENDIF
	 CALL diffw(dlnro_nu,grand_k,nu,vis_cin,n2mu,n2t,r,ro,ys,deff,dh,dv)
	ENDIF

c ~~~~~~~~~~~~~~~~~les coefficients pour les ZR~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c les C*
	cts(1)=1.d0/dt
	cts(2)=nu12/r**3/ro
	cts(3)=eps_nuc*nu12/lum
	cts(4)=dlldlm/nu
	cts(5)=ht*dh*nu12/r**4/ro/grand_k/delta
	cts(6)=ht*nu12/r**4/ro
	cts(7)=eps_nuc*eps_t*nu12/lum
	cts(8)=dlldlm*(chi_t+1.d0)/nu
	cts(9)=dlldlm*chi_mu/nu
	cts(10)=eps_nuc*eps_mu*nu12/lum
	cts(11)=rr_t/dt*y_t(1)
	cts(12)=-mw_dot
	cts(13)=dh/r**2
	cts(14)=r**3/nu52*(1.d0-dlldlm)
	cts(15)=cp*t*nu12/lum/delta/dt
	cts=cts*cts_0

c les C
	frl(1)=r**2*nu12/dt
	frl(2)=(cts(11)+cts(12))*nu12
	frl(3)=r**4*ro
	frl(4)=r**6*ro**2*dv/nu12
	frl(5)=r**3/ro/nu52
	frl(6)=cp*t*nu12/lum*dgrad/hp
	frl(7)=cts(14)
	frl(8)=nu12/dt*(r**2*(2.d0+dlnro_t)-rr_t)
	frl(9)=r**2*ro/nu12*dv
	frl(10)=cp*t*nu12*y_t(3)/lum/delta/dt
	frl(11)=cts(2)-cts(3)+cts(4)-cts(5)-cts(15)
	frl(12)=cts(6)-cts(7)+cts(8)
	frl(13)=dlldlm*r**2*ro*ht/nu32
	frl(14)=r**2*ro*ht/nu12
	frl(15)=cts(9)-cts(10)
	frl(16)=1.d0+chi_t		!chi_t non nul ici fait diverger
	frl(17)=1.d0
	frl(18)=chi_mu
	frl(19)=r**6*ro/nu**2
	frl(20)=cts(1)+cts(13)
	frl(21)=grad_mu_shp
	frl(22)=y_t(4)/dt
	frl(23)=1.d0/delta
	frl(24)=phi_e/delta
	frl=frl*cte_0

c arguments optionnels
	IF(PRESENT(coeff))THEN
	 coeff(1)=r					!R
	 coeff(2)=nu32					!M
	 coeff(3:8)=y(1:6,0)		!Omega, U, Théta, Lambda, Psi, Tau
	 IF(zc)THEN
	  coeff(4)=0.d0 ; coeff(6:7)=0.d0	!dans ZC U, Lambda, Psi=0
	 ENDIF	 
	 coeff(9)=deff					!Deff
	 coeff(10)=dh					!Dh
	 coeff(11)=dv					!Dv
	 coeff(12)=t					!T
	 coeff(13)=ro					!ro
	 coeff(14)=grad_mu_shp*hp			!grad_mu
	 coeff(15)=cte_0(6)*eps_nuc*nu32/lum/dlldlm	!f_epsilon
	 coeff(16)=eps_T				!epsilon_T
	 coeff(17)=eps_nuc				!epsilon
	 coeff(18)=dlldlm				!dlnL/dlnM
	 coeff(19)=eps_mu				!epsilon_mu
	 coeff(20)=chi					!chi
	 coeff(21)=chi_t				!chi_t
	 coeff(22)=chi_mu				!chi_mu
	 coeff(23)=grand_k				!K
	 coeff(24)=dlnro_nu				!d ln ro / dnu
	 coeff(25)=dlnro_t/dt*cte_0(8)			!d ln ro / dt
	 coeff(26)=y(1,1)				!d Omega / dnu
	 coeff(27)=y(2,1)				!d U / d nu
	 coeff(28)=n2t					!n2t
	 coeff(29)=n2mu					!n2mu
	 coeff(30)=mu					!mu
	 coeff(ncoeff+1:ncoeff+nchim)=xchimg		!Xchim
	 WHERE(ABS(coeff) < 1.d-30)coeff=0.d0
	ENDIF

	RETURN

	END SUBROUTINE coeff_rota3
