
c*************************************************************

	SUBROUTINE coeff_rota(dt,idis,ndis)

c	routine private du module mod_evol

c	on construit l'interpolation coef_rl des coefficients lissés
c	pour la rotation

c entrées :
c	dt : pas temporel

c	Auteur: P. Morel, Département Cassiopée, O.C.A., CESAM2k

c NOTATIONS (hélas incohérentes) pour les développements sur B-splines
c	n_ch : nombre VARIABLE de points élément de mod_variables
c	nch : nombre FIXE de fonctions élément de mod_donnees
c	m_ch : ordre FIXE des splines élément de mod_donnees 
c	mch(n_ch) : abscisses VARIABLES élément de mod_variables

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : aradia, clight, d_conv, en_masse,
	1 g, langue, lsol, m_ch, m_rl, m_ptm, m_rot, msol, ne, nchim,
	2 nom_rot, nucleo, nrot, pi, rsol, secon6, zi
	USE mod_etat, ONLY: etat, saha	
	USE mod_kind
	USE mod_nuc, ONLY : nuc
	USE mod_numerique, ONLY : bsp1dn, bsp_dis, no_croiss, sum_n
	USE mod_opa, ONLY : opa	
	USE mod_variables, ONLY : bp, bp_t, chim, chim_t, chim_gram,
	1 inter, knot, knotc, knotr, knotc_t, knotr_t, knot_ptm, knot_t,
	2 mc, mct, mc_t, mct_t, mrot, mrott, mrot_t, mrott_t, m23, m23_t,
	3 n_ch, n_ch_t, n_ptm, n_qs, n_qs_t, n_rot, n_rot_t,
	3 old_ptm, q, qt, qt_t, q_t, rota, rota_t, r2, r2_t, tot_conv,
	4 tot_rad, x_ptm, xt_ptm
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) :: dt
	INTEGER, INTENT(in), DIMENSION(0:) :: idis
	INTEGER, INTENT(in) :: ndis
		
	REAL (kind=dp), DIMENSION(0:NINT(MAXVAL(zi)),nchim) :: ioni
	REAL (kind=dp), DIMENSION(nrl,ndis) :: coef_rl_dis	
	REAL (kind=dp), DIMENSION(0,0) :: jac	
	REAL (kind=dp), DIMENSION(nchim) :: depsx, dxchim, dxchimm,
	1 xchim, xchimm, z_bar		
	REAL (kind=dp), DIMENSION(ne) :: dfqs, fqs
	REAL (kind=dp), SAVE, DIMENSION(nrl) :: cte_0
	REAL (kind=dp), DIMENSION(nrl) :: dfrl, frl	
	REAL (kind=dp), DIMENSION(nrot) :: dfrot, frot, y_t
	REAL (kind=dp), SAVE, DIMENSION(16) :: cts_0
	REAL (kind=dp), DIMENSION(16) :: cts
	REAL (kind=dp), DIMENSION(5) :: epsilon
	REAL (kind=dp), DIMENSION(1) :: dfm, fm	
	REAL (kind=dp), DIMENSION(0) :: dcomp
	REAL (kind=dp), SAVE ::	cte1_0, cte2_0, cte3_0, cte4_0	
	REAL (kind=dp) :: alfa, beta, be7e, b8e, chi, chi_mu, chi_t,
	1 cp, dcpp, dcpt, dcpx, deff, delta,
	2 deltap, deltat, deltax, depsro, depst, dgradadp, dgradadt,
	3 dgradadx, dh, dkapr, dkapt, dkapx, dlldlm, dlnmu, dlnro,
	4 drop, drot, drox, dup, dut, dux, dv, eps_mu, eps_nuc, eps_t,
	5 eta, f17e, gamma1, grad, gradad, grad_mu, grand_k, grav,
	6 hhe, hp, ht, kap, lum, mu, mw_dot, nel, nu, nu12,
	7 n13e, nu32, nu52, n2mu, n2t, o15e, p, phi, r, r_t,
	8 ro, t, u, vais, w0sro0

	INTEGER, ALLOCATABLE ,DIMENSION(:) :: id
	INTEGER, SAVE :: l=1
	INTEGER :: i, jdis
	
	LOGICAL, SAVE :: init=.TRUE.
		
c------------------------------------------------------------------------
	
2000	FORMAT(8es10.3)

c	définition des constantes

	IF(init)THEN
	 init=.FALSE.
	 
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1001)m_rl ; WRITE(2,1001)m_rl
1001	  FORMAT(/,'The coefficients for the angular mom. diffusion are',
	1  /,'smoothed by Bézier functions with B-splines of order=',i3)
	 CASE DEFAULT	 
	  WRITE(*,1)m_rl ; WRITE(2,1)m_rl
1	  FORMAT(/,'Les coefficients pour la diffusion du moment cinétique',
	1  /,'sont lissés par des Béziers avec des B-splines d''ordre',i3)
	 END SELECT
	 
c	pour le calcul de grad_mu, chi, grav	 
	 cte1_0=8.d0*pi*rsol**2/3.d0/msol
	 cte2_0=4.d0*aradia*clight/3.d0
	 cte3_0=g*msol/rsol**2
	 cte4_0=1.d0/2.d0/pi/g
	 
	 cte_0(1)=1.d0/secon6
	 cte_0(2)=1.d0
	 cte_0(3)=8.d0*pi*rsol**2/15.d0/msol
	 cte_0(4)=64.d0*pi**2*rsol**4/9.d0/msol**2
	 cte_0(5)=4.d0*rsol**3/3.d0/pi/g**2/msol
	 cte_0(6)=msol/lsol 
	 cte_0(7)=1.d0
	 cte_0(8)=2.d0*rsol**3/3.d0/g/lsol
	 cte_0(9)=2.d0*rsol**3/3.d0/g/lsol/secon6
	 cte_0(10)=msol/lsol/secon6
	 cte_0(11:12)=1.d0
	 cte_0(13)=8.d0*pi*rsol**2/3.d0/msol
	 cte_0(14)=16.d0*pi*rsol**2/9.d0/msol
	 cte_0(15)=1.d0	 
	 cte_0(16:18)=2.d0/3.d0
	 cte_0(19)=16.d0*pi*rsol**6/9.d0/g/msol**2
	 cte_0(20:21)=1.d0
	 cte_0(22)=1.d0/secon6
	 cte_0(23:25)=1.d0
 
	 cts_0(1)=1.d0/secon6
	 cts_0(2)=msol/2.d0/pi/rsol**3	 
	 cts_0(3)=msol/lsol	 
	 cts_0(4)=2.d0	 
	 cts_0(5:6)=3.d0*msol/2.d0/pi/rsol**4	 	 
	 cts_0(7)=msol/lsol	 
	 cts_0(8:9)=1.d0
	 cts_0(10)=msol/lsol
	 cts_0(11)=1.d0/secon6
	 cts_0(12)=1.d0/rsol**2
	 cts_0(13)=6.d0/rsol**2
	 cts_0(14)=8.d0*rsol**3/3.d0/g/msol
	 cts_0(15)=2.d0*rsol**3/3.d0/g/lsol/secon6
	 cts_0(16)=msol/lsol/secon6
	ENDIF	
	
c	nettoyage
	IF(ALLOCATED(coef_rl))DEALLOCATE(coef_rl,mrl,mrlt)

c facteur correctif Omega_0/2 pi G ro_0
	w0sro0=cte4_0*rota(1,1)/EXP(bp(7,1))

c cas totalement radiatif ou totalement convectif	
	IF(tot_conv .OR. tot_rad)THEN	
	 n_rl=n_rot ; knotrl=n_rl+m_rl
	 ALLOCATE(coef_rl(nrl,n_rl),mrl(n_rl),mrlt(knotrl))	
	 mrl=mrot
	 DO i=1,n_rl
	  nu=MAX(mrl(i),nu_min) ; frl=0.d0
	  IF(tot_conv)THEN
	   CALL z_conv
	  ELSE
	   CALL z_rad
	  ENDIF
	  coef_rl(:,i)=frl(:)
	 ENDDO
	 CALL bsp1dn(nrl,coef_rl,mrl,mrlt,n_rl,m_rl,knotrl,.FALSE.,nu,l,
	1  frl,dfrl,.TRUE.)
	 RETURN
	ENDIF
	
c avec ZR et ZC	
c	définitions et allocations pour l'interpolation spline lissée	
c	idis --> id pour pouvoir réinitialiser id(0:ndis+1) dans noeud_dis
	
	n_rl=n_rot ; knotrl=n_rl+ndis+m_rl	
	ALLOCATE(coef_rl(nrl,n_rl+ndis),id(SIZE(idis)),mrl(n_rl),
	1 mrlt(knotrl))	
	mrl=mrot ; id=idis		
	jdis=1
	DO i=1,n_rl
	 nu=MAX(mrl(i),nu_min) ; frl=0.d0
	 	
c est-on dans une zone de mélange?	 
	 CALL bsp1dn(nmix,fmix,x_mix,x_mixt,n_mix,m_mix,knot_mix,.TRUE.,
	1  nu,l,fm,dfm)	
	 IF(fm(1) > 0.d0)THEN
	  CALL z_rad
	  coef_rl(:,i)=frl(:)
	  
c est-on sur une limite? si oui, on ne peut être que convectif à droite
	  IF(i == idis(jdis))THEN
	   nu=un_eps*mrl(i)+eps*mrl(i+1)
	   CALL z_conv
	   coef_rl_dis(:,jdis)=frl(:)
	   jdis=MIN(jdis+1,ndis)
	  ENDIF	  
	 ELSE
	  CALL z_conv
	  coef_rl(:,i)=frl(:)
	  
c est-on sur une limite? si oui, on ne peut être que radiatif à droite
	  IF(i == idis(jdis))THEN
	   nu=un_eps*mrl(i)+eps*mrl(i+1)	  
	   CALL z_rad
	   coef_rl_dis(:,jdis)=frl(:)
	   jdis=MIN(jdis+1,ndis)
	  ENDIF
	 ENDIF  	  
	ENDDO

c	création du développement en B-splines avec discontinuités aux
c	limites ZR/ZC
c	lissage par contour avec .TRUE, interpolation bsp1dn avec .FALSE.
	
       CALL bsp_dis(nrl,mrl,coef_rl,ndis,id,coef_rl_dis,eps,n_rl,m_rl,mrlt,
	1 knotrl,.TRUE.)
	
c        CALL bsp_dis(nrl,mrl,coef_rl,ndis,id,coef_rl_dis,eps,n_rl,m_rl,
c	1 mrlt,knotrl)
	
	DEALLOCATE(id)
	
c	IF(.TRUE.)THEN
	IF(.FALSE.)THEN

	 PRINT*,'test interpolation'	
	 jdis=1
	 DO i=1,n_rl
	  nu=mrl(i)
	  CALL bsp1dn(nrl,coef_rl,mrl,mrlt,n_rl,m_rl,knotrl,.TRUE.,nu,l,
	1   frl,dfrl)
	  CALL bsp1dn(nmix,fmix,x_mix,x_mixt,n_mix,m_mix,knot_mix,.TRUE.,
	1   nu,l,fm,dfm)	
	  WRITE(*,2000)nu,fm(1),frl(4:6),frl(11),frl(25)
	  IF(i == idis(jdis))THEN
	   PRINT*,'limite ZR/ZC j=',jdis,idis(jdis),i
	   	  
	   nu=un_eps*mrl(i)+eps*mrl(i+1)
	   
	   PRINT*,mrl(i-1),mrl(i),nu,mrl(i+1)
	   
	   jdis=MIN(jdis+1,ndis)
	   CALL bsp1dn(nrl,coef_rl,mrl,mrlt,n_rl,m_rl,knotrl,.TRUE.,nu,l,
	1    frl,dfrl)
	   CALL bsp1dn(nmix,fmix,x_mix,x_mixt,n_mix,m_mix,knot_mix,.TRUE.,
	1    nu,l,fm,dfm)	
	   WRITE(*,2000)nu,fm(1),frl(4:6),frl(11),frl(25)
	  ENDIF
	 ENDDO
	 PAUSE'test interpolation fin'
	ENDIF	
	
	RETURN

	CONTAINS
	
c******************************************************************
	
	 SUBROUTINE z_conv
	 
c	 routine  subordonnée de coeff_rota
c	 les variables au temps t+dt dans une zone convective

c-----------------------------------------------------------------------------

2000	 FORMAT(8es10.3)

	 nu12=SQRT(nu)
	 CALL inter('m23',bp,q,qt,n_qs,knot,nu,fqs,dfqs,r2,m23)
	 ro=EXP(fqs(7))
	 IF(en_masse)THEN
	  r=SQRT(ABS(fqs(3)))
	 ELSE
	  r=ABS(fqs(3))
	 ENDIF
	  
c	 les variables au temps t
	  
	 CALL bsp1dn(1,old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm,.TRUE.,
	1 MAX(x_ptm(1),MIN(nu,x_ptm(n_ptm))),l,fm,dfm)	!masse au temps t f(1)
	 IF(no_croiss)PRINT*,'Pb. en 2 dans coeff_rota'	
	 IF(.NOT.en_masse)fm(1)=fm(1)**(2.d0/3.d0)
	 fm(1)=MAX(nu_min,MIN(fm(1),mrot_t(n_rot_t)))
	 CALL bsp1dn(nrot,rota_t,mrot_t,mrott_t,n_rot_t,m_rot,
	1 knotr_t,.TRUE.,MAX(mrot_t(1),MIN(fm(1),mrot_t(n_rot_t))),l,y_t,
	2 dfrot)
	 IF(no_croiss)PRINT*,'Pb en 3 dans coeff_rota'
 	 CALL inter('m23',bp_t,q_t,qt_t,n_qs_t,knot_t,
	1  MAX(m23_t(1),MIN(fm(1),m23_t(n_qs_t))),fqs,dfqs,r2_t,m23_t)
	 IF(no_croiss)PRINT*,'Pb. en 4 dans coeff_rota'
	 IF(en_masse)THEN
	  r_t=SQRT(ABS(fqs(3)))
	 ELSE
	  r_t=ABS(fqs(3))	 	 
	 ENDIF
	  
c Omega, U, theta, Lambda, psi au temps t+dt
	 CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,nu,l,
	1 frot,dfrot)	  
c perte de moment cinétique au temps t+dt 	 
	 CALL pertw(nu,frot(1),r,mw_dot)
	  
c	 les coefficients /= 0 dans une zone de mélange

	 deff=d_conv ; dv=d_conv
	 cts(11)=r**2/dt*y_t(1)*cts_0(11)
	 cts(12)=mw_dot/ro*cts_0(12)	 
	 frl(1)=r*nu12/dt*(3.d0*r-2.d0*r_t)	 
	 frl(2)=(cts(11)+cts(12))*nu12
	 frl(3)=r**4*ro
	 frl(4)=r**6*ro**2*dv/nu12	 	 	 
	 frl(1:4)=frl(1:4)*cte_0(1:4)	 	 
	 frl(25)=deff
	 
c	 WRITE(*,2000)nu,r,r_t,ro,ro_t,cts(11:12),y_t(1)
	  
	 RETURN 
	  
	 END SUBROUTINE z_conv

c*****************************************************************
	 	
	 SUBROUTINE z_rad
	 
c	 routine subordonnée de coeff_rota	 
c	 les variables au temps t+dt dans une zone radiative

c-----------------------------------------------------------------------

2000	 FORMAT(8es10.3)

	 nu12=SQRT(nu) ; nu32=nu12**3 ; nu52=nu32*nu
	 
	 CALL inter('m23',bp,q,qt,n_qs,knot,nu,fqs,dfqs,r2,m23)
	 p=EXP(fqs(1)) ; t=EXP(fqs(2))
	 IF(en_masse)THEN
	  r=SQRT(ABS(fqs(3))) ; lum=SQRT(ABS(fqs(4)))**3
	  dlnro=dfqs(7)/dfqs(5)
	 ELSE
	  r=ABS(fqs(3)) ; lum=fqs(4)
	  dlnro=3.d0/2.d0*dfqs(7)/dfqs(5)*nu12
	 ENDIF
	 dlldlm=dfqs(4)/dfqs(5)*fqs(5)/fqs(4)

c	 grad = d ln T / d ln P
 
	 grad=dfqs(2)/dfqs(1)

c	 hp, ht, grav
	
	 hp=-rsol*dfqs(3)/dfqs(1) ; IF(en_masse)hp=0.5d0*hp/r
	 ht=hp/grad ; grav=cte3_0*nu32/r**2	
	
c	 composition chimique
	
	 CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,knotc,.TRUE.,nu,
	1 l,xchim,dxchim)
	 IF(no_croiss)PRINT*,'Pb en 1 dans coeff_rota'	 
	
c	 composition chimique par gramme, équation d'état

	 xchimm=ABS(xchim) ; dxchimm=dxchim
	 CALL chim_gram(xchimm,dxchimm)
	 CALL etat(p,t,xchimm,.FALSE.,ro,drop,drot,drox,u,dup,dut,dux,
	1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)

c	 les variables au temps t

	 CALL bsp1dn(1,old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm,.TRUE.,
	1 nu,l,fm,dfm)	!masse au temps t f(1)
	 IF(no_croiss)PRINT*,'Pb en 2 dans coeff_rota'	
	 IF(.NOT.en_masse)fm(1)=fm(1)**(2.d0/3.d0)
	 fm(1)=MAX(nu_min,MIN(fm(1),mrot_t(n_rot_t)))
	 CALL bsp1dn(nrot,rota_t,mrot_t,mrott_t,n_rot_t,m_rot,
	1  knotr_t,.TRUE.,MIN(mrot_t(1),MAX(fm(1),mrot_t(n_rot_t))),l,
	1  y_t,dfrot)
	 IF(no_croiss)PRINT*,'Pb en 3 dans coeff_rota'	
 	 CALL inter('m23',bp_t,q_t,qt_t,n_qs_t,knot_t,
	1  MAX(m23_t(1),MIN(fm(1),m23_t(n_qs_t))),fqs,dfqs,r2_t,m23_t)
	 IF(no_croiss)PRINT*,'Pb. en 4 dans coeff_rota'	
	 IF(en_masse)THEN
	  r_t=SQRT(ABS(fqs(3)))
	 ELSE
	  r_t=ABS(fqs(3))	 	 
	 ENDIF

c opacité
	 CALL opa(xchimm,t,ro,kap,dkapt,dkapr,dkapx)	
	
c appel à saha pour les charges moyennes
	 CALL saha(xchim,t,ro,ioni,z_bar,nel,eta)

c poids moléculaire moyen
	 mu=1.d0/DOT_PRODUCT((1.d0+z_bar),xchim)
	 dlnmu=-mu*DOT_PRODUCT((1.d0+z_bar),dxchim)
	 grad_mu=-cte1_0*r**2*ro*hp/nu12*dlnmu

c	 phi = d ln ro / d ln mu phi étant mal connu, on pose

	 phi=1.d0
	 
c vaïssälä
c	 n2t=grav*delta/hp*(gradad-grad) ; n2mu=grav*phi/hp*grad_mu
c	 vais=ABS(gradad-grad+phi/delta*grad_mu)

	 n2t=grav*delta/hp*(gradad-grad)
	 n2mu=0.d0 ; vais=ABS(gradad-grad)	 

c conductivité radiative chi, dérivées et diffusivité thermique K
	 chi=cte2_0*t**3/ro/kap
	 chi_t=3.d0-t/kap*dkapt+delta
	 chi_mu=-nucleo(1)/mu/kap/(1.d0+z_bar(1))*dkapx-phi
	 grand_k=chi/ro/cp
		
c appel à nuc pour les epsilon
	 CALL nuc(t,ro,xchim,dcomp,jac,.TRUE.,3,
	1 epsilon,depst,depsro,depsx,hhe,be7e,b8e,n13e,o15e,f17e)
	 eps_nuc=epsilon(1)
	 IF(epsilon(1) > 0.d0)THEN
	  eps_mu=-SUM(depsx(1:nchim)/(1.d0+z_bar(1:nchim)))/mu/epsilon(1)
	  eps_t=t/epsilon(1)*depst
	 ELSE
	  eps_mu=0.d0 ; eps_t=0.d0
	 ENDIF

c Omega, U, theta, Lambda, psi au temps t+dt
	 CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,nu,l,
	1 frot,dfrot)
	 IF(no_croiss)PRINT*,'Pb en 4 dans coeff_rota'
	 	 
c Dh, Dv, Deff	 
	 CALL diffw(dlnro,grand_k,nu,n2mu,n2t,r,ro,dfrot,frot,deff,dh,dv)
	 
c perte de moment cinétique 	 
	 CALL pertw(nu,frot(1),r,mw_dot)
	
c les C*
	 cts(1)=1.d0/dt
	 cts(2)=nu12/r**3/ro	 
	 cts(3)=eps_nuc*nu12/lum
	 cts(4)=dlldlm/nu	 
	 cts(5)=ht*nu12*dh/r**4/ro/grand_k/delta
	 cts(6)=ht*nu12/r**4/ro	 
	 cts(7)=eps_nuc*eps_t*nu12/lum	 
	 cts(8)=dlldlm*(chi_t+1.d0)/nu
	 cts(9)=dlldlm*chi_mu/nu	 
	 cts(10)=eps_nuc*eps_mu*nu12/lum 
	 cts(11)=r**2/dt*y_t(1)
	 cts(12)=mw_dot/ro
	 cts(13)=dh/r**2
	 cts(14)=r**3/nu52*(1.d0+w0sro0-dlldlm)
	 cts(15)=r**3*cp*t*y_t(3)/lum/nu/delta/dt
	 cts(16)=cp*t*nu12/lum/delta/dt
	 	 
	 cts(:)=cts(:)*cts_0(:)
	 
c les C	  
	 frl(1)=r*nu12/dt*(3.d0*r-2.d0*r_t)	 
	 frl(2)=(cts(11)+cts(12))*nu12
	 frl(3)=r**4*ro
	 frl(4)=r**6*ro**2*dv/nu12
	 frl(5)=r**3/ro/nu52	 
	 frl(6)=cp*t*nu12/lum*vais/hp  
	 frl(7)=cts(14)+cts(15)  
	 frl(8)=r**3*cp*t*vais/lum/hp/nu	 
	 frl(9)=r**3*cp*t/lum/nu/delta/dt	 
	 frl(10)=cp*t*nu12*y_t(3)/lum/delta/dt	 
	 frl(11)=cts(2)-cts(3)+cts(4)-cts(5)-cts(16)	 
	 frl(12)=cts(6)-cts(7)+cts(8)	 
	 frl(13)=dlldlm*r**2*ro*ht/nu32
	 frl(14)=r**2*ro*ht/nu12	 
	 frl(15)=cts(9)-cts(10)	 
	 frl(16)=1.d0+chi_t	 
	 frl(17)=-1.d0
c	 frl(17)=chi_t*(1.d0-1.d0/delta)-1.d0 	 	 
	 frl(18)=chi_mu
	 frl(19)=r**6*ro/nu**2	 
	 frl(20)=cts(1)+cts(13)	 
	 frl(21)=grad_mu/hp
	 frl(22)=y_t(4)/dt	 
	 frl(23)=1.d0/delta
	 frl(24)=phi/delta
	 
	 frl(25)=deff
	 
	 frl(:)=frl*cte_0

	 RETURN
	 
	 END SUBROUTINE z_rad	
	END SUBROUTINE coeff_rota
