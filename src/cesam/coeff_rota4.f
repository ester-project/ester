
c*************************************************************

	SUBROUTINE coeff_rota4(dt,nu,y,frl,coeff)

c routine private du module mod_evol
c les coefficients frl pour la rotation avec Krot=4
c les coefficients coeff, pour les dessins
c le potentiel gravitationnel et sa dérivée sont inclus dans les frl

c entrées :
c	dt : pas temporel
c	nu : abscisse en m^2/3
c	y : variables de la rotation

c sorties :
c	frl  : coefficients de rotation
c	deff, dh, dv : coefficients de diffusion du moment cinétique
c	coeff : OPTIONAL, sortie des coefficients pour dessins

c Auteurs: P. Morel, Département Cassiopée, O.C.A., CESAM2k
c	   A. Moya Institut d'Astrophysique d'Andalousie,

c NOTATIONS (hélas incohérentes) pour les développements sur B-splines
c	n_ch : nombre VARIABLE de points élément de mod_variables
c	nch : nombre FIXE de fonctions élément de mod_donnees
c	m_ch : ordre FIXE des splines élément de mod_donnees
c	mch(n_ch) : abscisses VARIABLES élément de mod_variables

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : amu, aradia, clight, diffusion, d_conv,
	1 g, ihe4, lsol, m_ch, msol, nchim,
	2 nucleo, nrl, nrot, nvth, ord_rot, pi, rsol, secon6, x_tams, zi, z0
	USE mod_etat, ONLY: etat
	USE mod_kind
	USE mod_nuc, ONLY : l_planet, mzc_ext, nuc, nuzc_ext, planetoides
	USE mod_numerique, ONLY : bsp1dn, no_croiss
	USE mod_opa, ONLY : opa
	USE mod_variables, ONLY : chim, chim_gram, knotc, knotc_t, knotr,
	1 knotr_t, mc, mct, mc_t, mct_t, mrot, mrot_t, mrott, mrott_t, n_ch,
	2 n_ch_t, n_ptm, n_qs, n_qs_t, n_rot, n_rot_t, rota, rota_t, vth, vth_t

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(nrot,0:1) :: y
	REAL (kind=dp), INTENT(in) :: dt, nu

	REAL (kind=dp), INTENT(out), DIMENSION(nrl) :: frl
	REAL (kind=dp), INTENT(out), OPTIONAL, DIMENSION(ncoeff+nchim) :: coeff

	REAL (kind=dp), DIMENSION(2,0:1) :: ys
	REAL (kind=dp), DIMENSION(0,0) :: jac
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: cte_0
	REAL (kind=dp), DIMENSION(nchim) :: depsx, dxchim, dxchimg, xchim,
	1 xchimg
	REAL (kind=dp), DIMENSION(nvth) :: dfvth, fvth
	REAL (kind=dp), DIMENSION(nrot) :: dy_t, y_t
	REAL (kind=dp), SAVE, DIMENSION(13) :: cts_0
	REAL (kind=dp), DIMENSION(13) :: cts
	REAL (kind=dp), DIMENSION(5) :: epsilon
	REAL (kind=dp), DIMENSION(2) :: mw_planet
	REAL (kind=dp), DIMENSION(0) :: dcomp
	REAL (kind=dp), SAVE ::	cte2_0, cte3_0, cte4_0, cte5_0, cte6_0, ln13
	REAL (kind=dp) :: Abar, beta, be7e, b8e, chi, chi_mu, chi_t,
	1 cp, deff, delta, depsro, depst, dgrad,
	3  dh, dkapr, dkapt, dkapx, dlngrav, dlldlm, dlnro,
	5 dv, eps_mu, eps_nuc, eps_T, f17e, grad, gradad,
	6 grad_mu, grand_k, grav, hhe, hp, ht, kap, ln_lambda,
	7 lnmu, lnmu_t, lnro, lnro_t, lum, mu_s, mw_dot, nel,
	8 nu12, n13e, nu32, n2mu, n2t, o15e, p, phi, p_phi, r, ro,
	9 r2_t, t, vis_cin, Zbar

	INTEGER, SAVE :: l=1

	LOGICAL, SAVE :: init=.TRUE.
	LOGICAL :: zc

c------------------------------------------------------------------------

2000	FORMAT(8es10.3)

c définition des constantes
	IF(init)THEN
	 init=.FALSE.

c quelques coefficients constants
	 cte2_0=4.d0*aradia*clight/3.d0
	 cte3_0=g*msol/rsol**2
	 cte4_0=rsol**2
	 cte5_0=4.d0*pi*rsol**3/msol
	 cte6_0=2.d-15
	 ln13=LOG(1.3d4)

c les coefficients constants des F(nrl)
	 ALLOCATE(cte_0(nrl))
	 cte_0(1)=-rsol/3.d0/pi/g
	 cte_0(2)=-16.d0*rsol**4/27.d0/msol/g
	 cte_0(3)=2.d0*rsol/3.d0
	 cte_0(4)=8.d0*rsol**3/9.d0/g/msol
	 cte_0(5)=8.d0*rsol**2/3.d0/g/msol
	 cte_0(6)=-1.d0/pi/g/rsol
	 cte_0(7)=128.d0*pi*rsol**5/27.d0/msol**2/g
	 cte_0(8)=-16.d0*rsol**2/9.d0/msol/g
	 cte_0(9)=2.d0/rsol
	 cte_0(10)=-16.d0*pi*rsol**2/3.d0/msol
	 cte_0(11)=-2.d0/3.d0
	 cte_0(12)=-8.d0*pi*rsol**2/3.d0/msol
	 cte_0(13:14)=-1.d0
	 cte_0(15)=-msol/lsol/secon6
	 cte_0(16)=msol/lsol
	 cte_0(17)=16.d0*pi*rsol**2/9.d0/msol
	 cte_0(18)=-2.d0/3.d0
	 cte_0(19)=16.d0*pi*rsol**4/9.d0/msol
	 cte_0(20)=-1.d0
	 cte_0(21)=1.d0
	 cte_0(22)=-1.d0
	 cte_0(23)=1.d0
	 cte_0(24)=-1.d0/secon6
	 cte_0(25)=8.d0*pi*rsol**3/3.d0/msol
	 cte_0(26)=1.d0/secon6
	 cte_0(27)=1.d0
	 cte_0(28)=8.d0*pi*g*rsol**3/3.d0
	 cte_0(29)=4.d0*pi*g*rsol**3/3.d0
	 cte_0(30)=1.d0
	 cte_0(31)=8.d0*pi*rsol**2/15.d0/msol
	 cte_0(32)=64.d0*pi**2*rsol**4/9.d0/msol**2

c les coefficients constants des F*
	 cts_0(1)=msol/lsol
	 cts_0(2)=1.d0
	 cts_0(3)=-msol/lsol/secon6
	 cts_0(4)=-3.d0*msol/2.d0/pi/rsol**4
	 cts_0(5)=msol/lsol
	 cts_0(6)=-1.d0
	 cts_0(7)=-msol/lsol/secon6
	 cts_0(8)=6.d0/rsol**2
	 cts_0(9)=-1.d0/secon6
	 cts_0(10)=9.d0*msol/4.d0/pi/rsol**3
	 cts_0(11)=4.d0*pi*g*rsol
	 cts_0(12)=1.d0/secon6
	 cts_0(13)=-1.d0/rsol**2

	ENDIF

c mises à 0
	cts=0.d0 ; frl=0.d0 ;  mw_dot=0.d0

c nu12 etc...
	nu12=SQRT(nu) ; nu32=nu12**3

c est-on dans une ZC?
	zc=lmix(nu)

c~~~~~~~~~~~~~~~ au temps t ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c utilisation de la tabulation vth
c fvth(1)=lnP, fvth(2)=lnT, fvth(3)=r**2, fvth(4)=l**2/3, fvth(5)=ln ro,
c fvth(6)=cp, fvth(7)=delta, fvth(8)=gamma1, fvth(9)=ln µ, fvth(10)=ln kap
	CALL bsp1dn(nvth,vth_t,mc_t,mct_t,n_ch_t,m_ch,knotc_t,.TRUE.,
	1 MAX(mc_t(1),MIN(nu,mc_t(n_ch_t))),l,fvth,dfvth)
	IF(no_croiss)PRINT*,'coeff_rota4, Pb. en 1'

c quantités directement obtenues de la tabulation
	lnmu_t=fvth(9)
	lnro_t=fvth(5)

c quantités dérivées de la tabulation
	r2_t=ABS(fvth(3))

c la rotation
	CALL bsp1dn(nrot,rota_t,mrot_t,mrott_t,n_rot_t,ord_rot,knotr_t,.TRUE.,
	1 MAX(mrot_t(1),MIN(nu,mrot_t(n_rot_t))),l,y_t,dy_t)
	IF(no_croiss)PRINT*,'coeff_rota4, Pb. en 2'

c~~~~~~~~~~~~~~~~au temps t+dt	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	CALL bsp1dn(nvth,vth,mc,mct,n_ch,m_ch,knotc,.TRUE.,
	1 MAX(mc(1),MIN(nu,mc(n_ch))),l,fvth,dfvth)
	IF(no_croiss)PRINT*,'coeff_rota4, Pb. en 3'

c quantités directement obtenues de la tabulation
	cp=fvth(6)
	delta=fvth(7)
	dlnro=dfvth(5)
	lnmu=fvth(9)
	lnro=fvth(5)

c quantités dérivées de la tabulation
	dlldlm=dfvth(4)*nu/fvth(4)	!d ln L / d ln m
	grad=dfvth(2)/dfvth(1)
	grad_mu=dfvth(9)/dfvth(1)
	lum=SQRT(ABS(fvth(4)))**3
	p=EXP(fvth(1))
	ro=EXP(fvth(5))
	r=SQRT(ABS(fvth(3)))
	t=EXP(fvth(2))

c hp, ht~~~~~~~~~~~~~
	hp=r**2*p/ro/nu32/cte3_0
	ht=hp/grad

c gravité, d ln grav / d ln R
	grav=cte3_0*nu32/r**2 ; dlngrav=cte5_0*r**3*ro/nu32-2.d0

c gradient adiabatique
	gradad=p*delta/ro/cp/t

c approximation de phi
	beta=1.d0-aradia*t**4/p/3.d0
	phi=delta*beta/(4.d0-3.d0*beta)

c composition chimique
	CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,knotc,.TRUE.,
	1 MAX(mc(1),MIN(nu,mc(n_ch))),l,xchim,dxchim)
	IF(no_croiss)PRINT*,'coeff_rota4, Pb. en 5'

c A bar, Zbar, nel complètement ionisé
	Abar=DOT_PRODUCT(nucleo,xchim)/SUM(xchim)
	Zbar=DOT_PRODUCT(zi,xchim)/SUM(xchim)
	nel=DOT_PRODUCT(zi,xchim)*ro/amu

c estimation de  vis_cin, la viscosité cinématique, Lang 3.15, 3.25 et 3.28
	ln_lambda=ln13+1.5d0*LOG(t)-LOG(nel)/2.d0
	vis_cin=cte6_0*SQRT(t)**5*SQRT(Abar)/ro/Zbar**4/ln_lambda

c appel à nuc pour les epsilon
	CALL nuc(t,ro,xchim,dcomp,jac,.TRUE.,3,
	1 epsilon,depst,depsro,depsx,hhe,be7e,b8e,n13e,o15e,f17e)
	eps_nuc=epsilon(1)
	IF(eps_nuc > 0.d0)THEN
	 eps_T=t/eps_nuc*depst
	 eps_mu=phi*ro/eps_nuc*depsro
	ELSE
	 eps_mu=0.d0 ; eps_T=0.d0
	ENDIF

c composition chimique par gramme, opacité et dérivées
	xchimg=ABS(xchim) ; dxchimg=dxchim
	CALL chim_gram(xchimg,dxchimg)
	CALL opa(xchimg,t,ro,kap,dkapt,dkapr,dkapx)

ca.moya conductivité radiative chi, dérivées et diffusivité thermique K
	chi=cte2_0*t**3/ro/kap
	chi_t=3.d0-t/kap*dkapt+delta*(1.d0+ro/kap*dkapr)
	chi_mu=-phi*(1.d0+ro/kap*dkapr)
	grand_k=chi/ro/cp

c a.moya, p_phi = fonction thermodynamique (6.8) de Maeder & Zahn A&A 98
c n'est valable que si Z=cte, i.e. sans diffusion et sur la MS
c pour cohérence il faut utiliser mu simplifié et Z initial
	IF(xchim(1) < x_tams .OR. diffusion)THEN
	 p_phi=0.d0
	ELSE
	 mu_s=16.d0/(20.d0*xchimg(1)+12.d0-3.d0*z0)
	 p_phi=-mu_s*(12.d0*(z0-4.d0)*LOG(4.d0*(16.d0+3.d0*mu_s*(z0-4.d0)))
	1 +(32.d0-23.d0*z0)*LOG(-16.d0+mu_s*(32.d0-23.d0*z0))
	2 +5.d0*z0*LOG(5.d0*z0)
	3 +2.d0*(8.d0+3.d0*z0)*LOG(2.d0*(16.d0+mu_s*(8.d0+3.d0*z0))))/(80.d0*cp)
	ENDIF

c dgrad, n2t, n2mu
	 dgrad=ABS(gradad-grad)
	 n2t=grav*delta/hp*dgrad
	 n2mu=grav*phi*grad_mu/hp

c coefficients de diffusion du moment cinétique y --> ys
	IF(zc)THEN
	 deff=d_conv ; dh=d_conv ; dv=d_conv	 
	ELSE
c	 ys(1:2,0)=y_t(1:2) ; ys(1:2,1)=dy_t(1:2)
	 ys(1:2,:)=y(1:2,:)	!totalement implicite
	 CALL diffw(dlnro,grand_k,nu,vis_cin,n2mu,n2t,r,ro,ys,deff,dh,dv)
	ENDIF

c perte/gain de moment cinétique pour nu > nuzc_ext dus aux 
c différentes sources externes de moment cinétique:
c 1 - pertw: perte/gain de masse dû au vent
c 2 - mw_planet(1) apport des planétoïdes de vitesse angulaire non nulle
c 3 - correction affectée du signe -  dmw_planet(2) addition fictive
	IF(.FALSE.)THEN
c	IF(nu >= nuzc_ext)THEN
	 IF(lw_perte)CALL pertw(r,mw_dot)
	 IF(l_planet)THEN
	  CALL planetoides(mw_planet=mw_planet)
	  mw_dot=mw_dot+mw_planet(1)-cte4_0*mw_planet(2)*r**2*y(1,0)
c	  PRINT*,'mw_dot,cte4_0*mw_planet(2)*r**2*y(1,0)'
c	  WRITE(*,2000)mw_dot,cte4_0*mw_planet(2)*r**2*y(1,0),
c	1 mw_planet(2)*r**2*y(1,0)
c	  PAUSE'mw_dot'
	 ENDIF
	ENDIF  

c ~~~~~~~~~~~~~~~~~les coefficients~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ca.moya les F*
 	cts(1)=eps_nuc*nu12*(eps_T-delta)/lum
	cts(2)=dlldlm/nu*(2.d0*delta-1.d0-chi_t)
	cts(3)=nu12*cp*T/lum/dt
	cts(4)=ht*nu12/r**4/ro*(1.d0+Dh/grand_k)
	cts(5)=eps_nuc*nu12*(eps_mu+phi)/lum
	cts(6)=dlldlm/nu*(2.d0*phi+chi_mu)
	cts(7)=nu12*cp*T*p_phi*(lnmu-lnmu_t)/lum/dt	!bug signalé par ANDY
	cts(8)=dh/r**2
	cts(9)=(lnmu-lnmu_t-1.d0)/dt	
	cts(10)=nu12/r**3/ro
	cts(11)=r*ro/grav*dlnro
	cts(12)=r2_t*y_t(1)/dt		!r2_t est r_t**2
	cts(13)=mw_dot

c produit par les constantes
	cts=cts*cts_0

c les F
	frl(1)=r/ro/nu/grav*(dlngrav-2.d0)
	frl(2)=r**4/nu32/grav*(dlngrav-2.d0)
	frl(3)=r/nu/grav*(1.d0-dlldlm)*(dlngrav-2.d0)
	frl(4)=r**3/nu32
	frl(5)=r**2/grav/nu32
	frl(6)=dlngrav/grav/ro/nu/r
	frl(7)=r**5*ro/grav/nu**2
	frl(8)=r**2/grav/nu32*dlngrav
	frl(9)=(1.d0-dlldlm)/grav/nu/r*dlngrav
	frl(10)=r**2*ro/grav/nu32*(1.d0-dlldlm)	
	frl(11)=phi+chi_mu
	frl(12)=r**2*ro*ht/nu32*dlldlm
	frl(13)=cts(1)+cts(2)+cts(3)+cts(4)
	frl(14)=cts(5)+cts(6)+cts(7)
	frl(15)=nu12*cp*T*y_t(3)/lum/dt
	frl(16)=cp*T/lum/hp*nu12*dgrad
	frl(17)=r**2*ro*ht/nu12
	frl(18)=1.d0-delta+chi_t	!chi_t non nul ici fait diverger
	frl(19)=r**4*ro/grav/nu12
	frl(20)=phi
	frl(21)=delta
	frl(22)=grad_mu/hp
	frl(23)=cts(8)+cts(9)
	frl(24)=y_t(4)/dt
	frl(25)=r**3*ro/nu12
	frl(26)=r**2*nu12/dt
	frl(27)=cts(10)+cts(11)
	frl(28)=r**3*ro/grav
	frl(29)=r**3*ro/grav*dlnro
	frl(30)=(cts(12)+cts(13))*nu12
	frl(31)=r**4*ro
	frl(32)=r**6*ro**2/nu12*dv 

c produit par les constantes
	frl=frl*cte_0

c arguments optionnels
	IF(PRESENT(coeff))THEN
	 coeff(1)=r					!R
	 coeff(2)=nu32					!M
	 coeff(3:6)=y(1:4,0)				!Omega, U, Psi, Lambda
	 coeff(7)=y(7,0)				!Phi 	 
	 coeff(8)=y(5,0)				!Tau
	 coeff(9)=deff					!Deff
	 coeff(10)=dh					!Dh
	 coeff(11)=dv					!Dv
	 coeff(12)=t					!T
	 coeff(13)=ro					!ro
	 coeff(14)=grad_mu				!grad_mu
	 coeff(15)=cte_0(16)*eps_nuc*nu32/lum/dlldlm	!f_epsilon
	 coeff(16)=eps_T				!epsilon_T
	 coeff(17)=eps_nuc				!epsilon
	 coeff(18)=dlldlm				!dlnL/dlnM
	 coeff(19)=eps_mu				!epsilon_mu
	 coeff(20)=chi					!chi
	 coeff(21)=chi_t				!chi_t
	 coeff(22)=chi_mu				!chi_mu
	 coeff(23)=grand_k				!K
	 coeff(24)=dlnro	 			!d ln ro / dnu
	 coeff(25)=(lnro-lnro_t)/dt/secon6		!dln ro / dt
	 coeff(26)=y(1,1)				!d Omega / dnu
	 coeff(27)=y(2,1)				!d U / d nu
	 coeff(27)=y(2,0)				!dans ZC U fictif
	 coeff(28)=n2t					!n2t
	 coeff(29)=n2mu					!n2mu
	 coeff(30)=EXP(lnmu)				!mu
	 coeff(ncoeff+1:ncoeff+nchim)=xchimg		!Xchim
	 WHERE(ABS(coeff) < 1.d-30)coeff=0.d0
	ENDIF

	RETURN

	END SUBROUTINE coeff_rota4
