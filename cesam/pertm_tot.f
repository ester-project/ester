
c**********************************************************************

	SUBROUTINE pertm_tot(dt)

c	routine private du module mod_static
c	interpolation m(t+dt)--->m(t) en tenant compte 	
c	des variations de masse dues a E=mc**2
c	de la perte de masse (mdot>0 : gain de masse,
c	mdot<0 : perte de masse)
c	la perte de masse externe est supposée concentrée dans la
c	couche n_ptm-1 n_ptm

c	on entre avec mstar=masse au temps t'
c	on sort avec mstar=masse au temps t+dt

c	on intègre les variations dues a E=mc**2
c	d'où la variation de masse totale et la nouvelle masse
c	totale mstar
c	calcule la masse qui correspondait au temps t a la masse
c	en chaque point
c	tabule l'ancienne masse en fonction de la nouvelle (en m**2/3)
c	on néglige la variation de masse de l'atmosphère entre les
c	instants t+dt et t

c	utilisation par bsp1dn
c	en lagrangien m^2/3 ---> m^2/3 ancien
c	en eulérien m ---> m ancien

c	avec pression turbulente 7 inconnues 
c	sans pression turbulente 6 inconnues, Ptot=Pgaz

c	25 08 97 : mise en place des variables eulériennes
c	30 07 00 : introduction F95

c	Auteur: P. Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c entrées
c	bp,q,n_qs,qt,knot,chim,mc,nc,mct,knotc : modele au temps t
c	dt : pas temporel

c entrées/sorties
c	mstar : masse totale avec perte de masse
c	en entrée au temps t, en sortie au temps t+dt 

c sorties
c	old_ptm,x_ptm,xt_ptm,n_ptm,knot_ptm : interpolation de
c	l'ancienne masse
c	en fonction de la nouvelle masse en m^2/3 en lagrangien
c	en m en eulérien
C
c----------------------------------------------------------------

	USE mod_donnees, ONLY : clight, en_masse, eve,
	1 nchim, mdot, m_ch, m_ptm, m_qs, ne, nom_pertm, pturb,
	2 ord_qs, secon6, t_inf
	USE mod_etat, ONLY : etat
	USE mod_kind	
	USE mod_nuc, ONLY : nuc
	USE mod_numerique, ONLY : bsp1dn, no_croiss, sum_n
	USE mod_variables, ONLY : bp, chim, chim_gram, knot, knotc,
	1 knot_ptm,
	1 mc, mct, mstar, mstar_t, n_ch, n_ptm, n_qs, old_ptm, q, qt,
	2 x_ptm, xt_ptm

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) :: dt

	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: jac
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: dm	
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: df, dxchim,
	1 depsx, f, xchim, xchimm
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: vt1, vt2, m, mt	
	REAL (kind=dp), SAVE, DIMENSION(5) :: epsilon
	
	REAL (kind=dp), SAVE :: cte1, ehh, ebe7, eb8, ef17, en13,
	1 eo15, mev
	REAL (kind=dp) :: p, t, cte2, ro, drop, drot, drox,
	1 u, dup, dut, dux, depst, depsro, hh, be7, b8, n13, o15, f17,
	2 bid, delta, deltap, deltat, deltax, cp, dcpp, dcpt, dcpx,
	3 gradad, dgradadp, dgradadt, dgradadx, alfa, beta, gamma1

	INTEGER, SAVE :: np_ptm=-100

	INTEGER :: i, l
	
	LOGICAL, SAVE :: init=.TRUE.
	
c------------------------------------------------------------------------
	
2000	FORMAT(8es10.3)
	
	IF(init)THEN
	 init=.FALSE. 	 
	 cte1=secon6/clight/clight
	 mev=1.d6*eve	!energie des neutrinos selon Clayton p.380 et 392
	 ehh=mev*0.263d0 ; ebe7=mev*0.80d0 ; eb8=mev*7.2d0
	 en13=mev*0.71d0 ; eo15=mev*1.d0 ; ef17=mev*0.94d0	 
	 WRITE(*,1) ; WRITE(2,1)
1	 FORMAT(/,'on tient compte du delta Masse due à E=mc**2',/)

	 ALLOCATE(df(ne),f(ne),xchim(nchim),dxchim(nchim),xchimm(nchim),
	1 depsx(nchim))

	ENDIF	 
	
	ALLOCATE(dm(1,n_qs),m(n_qs),mt(n_qs+6),vt2(n_qs))
	
	cte2=cte1*dt

c	extraction des masses
c	calcul des réactions thermonucléaires puis des pertes par neutrinos
c	m(i)=delta m du a e=mc2 pendant dt en fraction de Mstar x Msol
	
	dm=0.d0 ; m=0.d0 ; x_ptm=0.d0
	DO i=2,n_qs
	 CALL bsp1dn(ne,bp,q,qt,n_qs,m_qs,knot,.TRUE.,q(i),l,f,df)
c	 WRITE(*,2000)f
	 IF(pturb)THEN		!avec pression turbulente 7 inconnues
	  p=exp(f(7))
	 ELSE			!sans pression turbulente 6 inconnues
	  p=exp(f(1))
	 ENDIF
	 t=exp(f(2)) ; vt2(i)=f(5)	!soit en m23 ou en m
	 IF(en_masse)THEN
	  bid=abs(f(5)) ; m(i)=sqrt(bid)**3	!m est en Msol
	 ELSE
	  m(i)=abs(f(5)) ; bid=m(i)**(2.d0/3.d0)
	 ENDIF
c	 WRITE(*,2000)m(i),f(5)
	 CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1 knotc,.TRUE.,min(bid,mc(n_ch)),l,xchim,dxchim)
	 IF(t > t_inf)THEN
	  xchimm=xchim
	  CALL chim_gram(xchimm,dxchim)	  
	  CALL etat(p,t,xchimm,.FALSE.,
	1 ro,drop,drot,drox,u,dup,dut,dux,
	2 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	  CALL nuc(t,ro,xchim,dxchim,jac,.FALSE.,3,	!3 : epsilon
	1 epsilon,depst,depsro,depsx,hh,be7,b8,n13,o15,f17)
	
c	  hh,be7,b8,n13,o15,f17 : nombre de neutrinos par unité de masse
c	  et de temps pour les diverses réactions
	
	  CALL nuc(t,ro,xchim,dxchim,jac,.FALSE.,4,	!4 : neutrinos
	1 epsilon,depst,depsro,depsx,hh,be7,b8,n13,o15,f17)
	  dm(1,i)=(epsilon(1)+hh*ehh+be7*ebe7+b8*eb8+n13*en13+
	1 o15*eo15+f17*ef17)*cte2
	 ENDIF	 
c	 WRITE(*,2000)m(i),f(5),dm(1,i)
	ENDDO
	
c	WRITE(*,2000)dt,cte2,mstar_t ; WRITE(*,2000)dm(1,:)
c	WRITE(*,2000)m

c	suppression des inversions, elles sont anormales avec perte
c	de masse mais peuvent se produire avec gain de masse

	n_ptm=n_qs ; i=n_ptm-1
	DO WHILE(i > 1)
	 IF(m(i) >= m(i+1))THEN
	  n_ptm=n_ptm-1	 
c	  PRINT*,'perte_tot, inversion en',i,n_ptm,m(i),m(i+1) ; PAUSE
	  DO l=i,n_ptm
	   m(l)=m(l+1) ; dm(1,l)=dm(1,l+1) ; vt2(l)=vt2(l+1)
	  ENDDO
	  i=min(i,n_ptm-1)
	 ELSE
	  i=i-1
	 ENDIF
	ENDDO	
		
c	spline de dm en fonction de m pour intégration
	
	CALL bsp1dn(1,dm,m,mt,n_ptm,m_ptm,knot_ptm,.FALSE.,m(1),l,xchim,
	1 dxchim)
        IF(no_croiss)THEN
         PRINT*,'Arrêt 1 dans pertm_tot' ; STOP
        ENDIF
	
c	integration mdot<0, pour perte de masse	
c	mstar(t+dt)=mstar(t)+m_dot*dt-defaut de masse
c	mstar contient l'atmosphere alors que m(n_ptm) ne la contient pas
	
	CALL sum_n(1,dm,mt,m_ptm,knot_ptm,.FALSE.,m(1),m(n_ptm),f)
	
	mstar=mstar_t+mdot*1.d6*dt-f(1)
	
c	on ajoute la masse perdue par perte de masse due a E=mc2
c	old_ptm est initialisé avec la masse au temps t+dt

	ALLOCATE(vt1(n_ptm))
	DO i=1,n_ptm
	 CALL sum_n(1,dm,mt,m_ptm,knot_ptm,.TRUE.,m(1),m(i),f)
	 vt1(i)=m(i)+f(1)
c	 WRITE(*,2000)f(1)
	ENDDO
	
c	en lagrangien on repasse en m**2/3

	IF(en_masse)THEN
	 DO i=1,n_ptm
	  vt1(i)=vt1(i)**(2.d0/3.d0)
	 ENDDO
	ENDIF
	
c	suppression des masses négatives anormales

	vt1(1)=0.d0 ; vt2(1)=0.d0
	i=2
	DO WHILE(i <= n_ptm)
	 IF(vt1(i) <= 0.d0 .OR. vt2(i) <= 0.d0)THEN
	  n_ptm=n_ptm-1
	  DO l=i,n_ptm
	   vt1(l)=vt1(l+1) ; vt2(l)=vt2(l+1)
	  ENDDO
	 ELSE
	  i=i+1
	 ENDIF
	ENDDO
	
c	suppression des inversions, elles sont anormales avec
c	perte de masse mais peuvent se produire avec gain de masse

	i=n_ptm-1
	DO WHILE(i > 1)
	 IF(vt1(i) >= vt1(i+1) .OR. vt2(i) >= vt2(i+1))THEN
c	  PRINT*,'perte_ext, inversion en',i,n_ptm,
c	1 vt1(i),vt1(i+1),vt2(i),vt2(i+1)
c	  PAUSE
	  n_ptm=n_ptm-1
	  DO l=i,n_ptm
	   vt1(l)=vt1(l+1) ; vt2(l)=vt2(l+1)
	  ENDDO
	  i=min(i,n_ptm-1)
	 ELSE
	  i=i-1
	 ENDIF
	ENDDO

c	allocations
	
	IF(ALLOCATED(old_ptm))THEN
	 IF(n_ptm /= np_ptm)THEN
	  DEALLOCATE(old_ptm,x_ptm,xt_ptm)
	  ALLOCATE(old_ptm(1,n_ptm),x_ptm(n_ptm),xt_ptm(n_ptm+m_ptm))
	 ENDIF
	ELSE
	 ALLOCATE(old_ptm(1,n_ptm),x_ptm(n_ptm),xt_ptm(n_ptm+m_ptm))
	ENDIF
	old_ptm(1,1:n_ptm)=vt1(1:n_ptm) ; x_ptm(1:n_ptm)=vt2(1:n_ptm)
		
c	tabulation

	CALL bsp1dn(1,old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm,.FALSE.,
	1 x_ptm(1),l,f,df)
        IF(no_croiss)THEN
         PRINT*,'Arrêt 2 dans pertm_tot' ; STOP
        ENDIF	
c	PRINT*,mstar,dt ; PRINT*,'mtot-mstar/mtot,mstar'	
c	WRITE(*,2000)mtot-mstar ; PRINT*,mtot,mstar
c	PAUSE'sortie pertm_tot'

	np_ptm=n_ptm
	
	DEALLOCATE(dm,m,mt,vt1,vt2)

	RETURN

	END SUBROUTINE pertm_tot

