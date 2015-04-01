
c***********************************************************************

	SUBROUTINE eq_diff_chim(nu,y,dt,ad,as,bd,bs)	
	
c routine private du module mod_evol

c formation des équations à résoudre dans diffus par la méthode de Galerkin
c les dérivées/xchim sont calculées dans la routine de calcul
c du coefficient de diffusion
c on n'a pas a tenir compte de la pression turbulente puisqu'elle
c est nulle dans les zones radiatives où il y a diffusion
c microscopique et qu'il y a mélange dans les zones convectives

c on suppose que la perte de masse pour l'élément X est X Mdot
c en faisant porter l'intégrale de Petrov-Galerkin sur [0, M(t+dt)]
c cette perte de masse est automatiquement prise en compte
c pour une autre valeur de la perte de masse il faudrait intégrer
c sur [0, M(t)] en tenant compte du fait que sur [M(t+dt),M(t)]
c les abondances sont nulles... donc modifier les produits scalaires.

c Deff coefficient de diffusionn turbulente induit par la rotation est
c calculé directement, pour fixer Deff=0, voir lignes 232 et 285

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c entrées
c	nu: point de calcul en m**2/3
c	dt: pas temporel
c	y: solution intermédiaire

c sorties	
c	as, ad: coefficients de la spline(S)/spline dérivée(dS)
c	bs, bd: seconds membres

c_______________________________________________________________________

c pour chaque spline Nj le système est constitué de nchim produits scalaires
c les équations:

c	< 3/2 nu^1/2 [ Xi(t+1)-Xi(t) - Psi dt ] , Fi_j > +
c	+ < 32/(3 nu^1/2) pi^2 r^4 ro^2 dt
c	[Di1 dX1/dnu+...+Din dXn/dnu],dFi_j/dnu >
c	- < 4 pi r^2 ro dt Vi Xi , dFi_j/dnu >=0

c	nu: m**2/3

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : ah, amu, aradia, clight, en_masse, g, ihe4, Ipg,
	1 Krot, lsol, mdot, mu_saha, m_ch, m_ptm, m_qs, msol, nchim,
	2 ne, nom_difft, nom_elem, nrot, nucleo, ord_rot, pi, pturb, rsol,
	3 secon6, zi
	USE mod_etat, ONLY : etat, saha
	USE mod_kind
	USE mod_nuc, ONLY : l_planet, l_vent, nuc, nuzc_ext, planetoides, vent
	USE mod_numerique, ONLY : bsp1dn, no_croiss
	USE mod_opa, ONLY : opa
	USE mod_variables, ONLY : bp, chim_gram, chim_t, inter, knot,
	1 knotc, knotc_t, knotr, knot_ptm, mc, mct, mc_t, mrot, mrott,
	2 mct_t, m23, n_ch, n_ch_t, n_ptm, n_qs, n_rot, old_ptm, q, qt,
	3 rota, r2, xt_ptm, x_ptm, wrot
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(nchim,0:1) :: y
	REAL (kind=dp), INTENT(in) :: nu, dt
	REAL (kind=dp), INTENT(out), DIMENSION(nchim,nchim,0:1) :: ad, as
	REAL (kind=dp), INTENT(out), DIMENSION(nchim) :: bd, bs
	
	REAL (kind=dp), DIMENSION(nchim,nchim,nchim) :: dd
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: ioni		
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: yd, jac0
	REAL (kind=dp), DIMENSION(nrot,0:1) :: ys
	REAL (kind=dp), DIMENSION(nchim,nchim) :: d, dv, jac		
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: bds, bss, dcomp0
	REAL (kind=dp), DIMENSION(nchim) :: comp_t, dcomp,
	1 dcomp_t, depsx, dxchim, dychim, v, xchim, ychim, z_bar
	REAL (kind=dp), DIMENSION(ne) :: dfqs, fqs
	REAL (kind=dp), DIMENSION(nrot) :: dfrot, frot		
	REAL (kind=dp), DIMENSION(5) :: epsilon
	REAL (kind=dp), DIMENSION(1) :: dfm, fm	
	REAL (kind=dp), PARAMETER :: dx=1.d-5, unpdx=1.d0+dx
	REAL (kind=dp), SAVE :: cte1_0, cte2_0, cte3_0, cte4_0, cte5_0, cte6_0,
	1 cte8_0, cte9_0, ln13
	REAL (kind=dp) :: alfa, Abar, beta, be7e, bid, b8e, chi, cp,
	1 cte1, cte2, cte3, dcpp, dcpt, dcpx, dcte2x, dcte3x, deff, delta,
	2 deltap, deltat, deltax, depsro, depst, grad, gradad, gradrad,
	3 dgradadx, dgradadp, dgradadt, dgradradx, dh, dkapro, dkapt,
	3 dkapx, dlnmu, dlnp, dlnro, drop, drot, drox, dstor, dup, dut,
	5 dux, dvv, eta, f17e, gamma1, grad_mu_shp, grand_k, grav, hhe, hp, kap,
	6 l, ln_lambda, mu, nel, nu12, nu32, n13e, n2mu, n2t,
	7 vis_cin, o15e, p, phi, r, ro, u, stor0, stor, t, w, Zbar

	INTEGER, SAVE :: lx=1	
	INTEGER	:: i, j, k, kd	

	LOGICAL, SAVE :: init=.TRUE.
	LOGICAL :: deriv=.FALSE., melange
		
c----------------------------------------------------------------------

c le nombre d'inconnues est nchim
c y(1:nchim,0:der): dérivée d'ordre der
c (der=0 pour la fonction, der=1 pour dérivée première)
c de la variable var y : 1, .. nchim

2000	FORMAT(8es10.3)
2001	FORMAT(10es8.1)
2002	FORMAT(9es9.2)

c	WRITE(*,*)fait,nu ; WRITE(*,2000)dt ; WRITE(*,2001)y(1:nchim,0)
c	PAUSE'entree eq_diff_chim'

	IF(init)THEN
	 init=.FALSE.
	 cte1_0=1.5d0 
	 cte2_0=2.d0/3.d0*(4.d0*pi*rsol**2/msol)**2*secon6
	 cte3_0=4.d0*pi*rsol**2/msol*secon6
	 cte4_0=3.d0/16.d0/pi/aradia/clight/g*lsol/msol
	 cte5_0=4.d0*aradia*clight/3.d0
	 cte6_0=g*msol/rsol**2
	 cte8_0=8.d0*pi*rsol**2/3.d0/msol
	 cte9_0=2.d-15
	 ln13=LOG(1.3d4)

c si on utilise saha pour les taux d'ionisation
	 IF(mu_saha)ALLOCATE(ioni(0:NINT(MAXVAL(zi)),nchim))
	ENDIF

c les conditions limites sont incluses dans les parties intégrées
c et sont nulles dXi/dnu=0 au centre et à la surface
c M, P, T, R, L au point nu=m**2/3	
	nu12=SQRT(nu) ; nu32=nu12**3
	IF(en_masse)THEN
	 bid=MIN(nu,m23(n_qs))
	ELSE
	 bid=MIN(nu32,m23(n_qs))
	ENDIF

c variables quasi-statiques au point courant, les dérivées sont / bid		 
	CALL inter('m23',bp,q,qt,n_qs,knot,bid,fqs,dfqs,r2,m23)
	IF(pturb)THEN		!avec pression turbulente
	 p=EXP(fqs(Ipg)) 
	ELSE			!sans pression turbulente
	 p=EXP(fqs(1))	 
	ENDIF
	t=EXP(fqs(2))
	IF(en_masse)THEN
	 r=SQRT(ABS(fqs(3))) ; l=SQRT(ABS(fqs(4)))**3 ; dlnp=dfqs(1)
	ELSE
	 r=abs(fqs(3)) ; l=fqs(4) ; dlnp=3.d0/2.d0*dfqs(1)*nu12
	ENDIF
	
c hp, grad	
	hp=-rsol*dfqs(3)/dfqs(1) ; IF(en_masse)hp=0.5d0*hp/r
	grad=dfqs(2)/dfqs(1)
	 		
c	PRINT*,'r',r ; WRITE(*,2000)p,t,r,l,nu32

c est-on dans une zone mélangée ?
	melange=lmix(nu)
	
c avec difft_gab il y a mélange si T < t_gab
	IF(nom_difft == 'difft_gab')melange=melange .OR. t < t_gab

c la composition chimique au temps t
	IF(en_masse)THEN
	 u=MIN(nu,x_ptm(n_ptm))	!u:vT
	ELSE
	 u=MIN(nu32,x_ptm(n_ptm))	!u:vT
	ENDIF	
	CALL bsp1dn(1,old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm,.TRUE.,
	1 u,lx,fm,dfm)	!masse au temps t f(1)
	IF(no_croiss)PRINT*,'Pb en 2 dans eq_diff_chim'	
	IF(.NOT.en_masse)fm(1)=fm(1)**(2.d0/3.d0)
	fm(1)=MAX(0.d0,MIN(fm(1),mc_t(n_ch_t)))
	CALL bsp1dn(nchim,chim_t,mc_t,mct_t,n_ch_t,m_ch,
	1 knotc_t,.TRUE.,fm(1),lx,comp_t,dcomp_t)
	IF(no_croiss)PRINT*,'Pb en 3 dans eq_diff_chim'	
c	WRITE(*,2000)comp_t,y(:,0)

c ychim composition chimique par mole ---> xchim composition chimique par gramme
	ychim(:)=y(:,0) ; dychim=y(:,1) ; xchim=ychim ; dxchim=dychim
	CALL chim_gram(xchim,dxchim)		
	
c la densité, le gradient adiabatique et leurs dérivées/X
	CALL etat(p,t,xchim,.TRUE.,ro,drop,drot,drox,u,dup,dut,dux,
	1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	
c l'opacité et le gradient radiatif
	CALL opa(xchim,t,ro,kap,dkapt,dkapro,dkapx)	
	gradrad=cte4_0*l*p*kap/nu32/t**4
	dgradradx=gradrad/kap*(dkapx+dkapro*drox)

c initialisation des coefficients de diffusion microscopique
	v=0.d0 ; dv=0.d0 ; d=0.d0 ; dd=0.d0

c les coefficients de diffusion microscopique et Deff
	IF(melange)THEN
	
c le coefficient de diffusion turbulente sera d_conv dans difft	
	 deff=0.d0
	ELSE
	
c la vitesse angulaire w et le coefficient de diffusion Deff
	 SELECT CASE(Krot)
	 CASE(0,1,2)	 	 
	  w=wrot ; deff=0.d0
	  
c diffusion du moment cinétique	Deff  
	 CASE(3,4)	!commenter  pour supprimer le calcul de Deff
c	 CASE DEFAULT	!décommenter  pour supprimer le calcul de Deff	 
	 
	  CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,nu,
	1 lx,frot,dfrot)
	  IF(no_croiss)PRINT*,'Pb en 5 dans eq_diff_chim' 	
	  w=frot(1) ; ys(:,0)=frot ; ys(:,1)=dfrot
c	  WRITE(*,2000)nu,w

c conductivité radiative chi, dérivées et diffusivité thermique K
	  chi=cte5_0*t**3/ro/kap
	  grand_k=chi/ro/cp
	 
c appel à saha pour les charges moyennes, ychim par mole, xchim par gramme
c poids moléculaire moyen
	  IF(mu_saha)THEN
	   CALL saha(ychim,t,ro,ioni,z_bar,nel,eta)
	   mu=1.d0/DOT_PRODUCT((1.d0+z_bar),ychim)
	   dlnmu=-mu*DOT_PRODUCT((1.d0+z_bar),dychim)
	  ELSE	  
	   mu=16.d0/(23.d0*xchim(1)+3.d0*xchim(ihe4)+9.d0)
	   dlnmu=-mu/16.d0*(23.d0*dxchim(1)+3.d0*dxchim(ihe4))	   	  
	  ENDIF	  
	  grad_mu_shp=-cte8_0*r**2*ro/nu12*dlnmu
	  
c estimation de  vis_cin, la viscosité cinématique, Lang 3.15, 3.25 et 3.28
	  Abar=DOT_PRODUCT(ychim,nucleo)/SUM(ychim)
	  IF(mu_saha)THEN
	   Zbar=DOT_PRODUCT(z_bar,ychim)/SUM(ychim)
	   nel=DOT_PRODUCT(z_bar,ychim)*ro/amu	 
	  ELSE
	   Zbar=DOT_PRODUCT(zi,ychim)/SUM(ychim)
	   nel=DOT_PRODUCT(zi,ychim)*ro/amu	 
	  ENDIF	 	
	  ln_lambda=ln13+1.5d0*LOG(t)-LOG(nel)/2.d0
	  vis_cin=cte9_0*SQRT(t)**5*SQRT(Abar)/ro/Zbar**4/ln_lambda
	 
c gravité	
	  grav=cte6_0*nu32/r**2
	  
c phi = d ln ro / d ln mu [approx. de Cox & Guili phi/delta=beta/(4-3beta)]
	  phi=delta*beta/(4.d0-3.d0*beta)
	 	 
c n2t, n2mu
	  n2t=grav*delta/hp*ABS(gradad-grad) ; n2mu=grav*phi*grad_mu_shp

c d ln ro
	  dlnro=(alfa-delta*grad)*dlnp+phi*dlnmu

c calcul de Dh, Dv, Deff
	  CALL diffw(dlnro,grand_k,nu,vis_cin,n2mu,n2t,r,ro,ys,deff,dh,dvv)

c vitesse angulaire
	 CASE(5)	!commenter  pour supprimer le calcul de Deff
c	 CASE(3,4,5)	!décommenter  pour supprimer le calcul de Deff
	 	  
	  CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,nu,
	1 lx,frot,dfrot)
	  IF(no_croiss)PRINT*,'Pb en 6 dans eq_diff_chim' 	
	  w=frot(1) ; deff=0.d0
	 END SELECT

c le coefficient de diffusion microscopique
c inactif s'il y a mélange ou initialisation de la rotation
	 CALL diffm(p,t,r,l,nu32,ro,drox,kap,dkapx,w,gradrad,dgradradx,ychim,
	1 d,dd,v,dv)
	ENDIF		!melange

c le coefficient de diffusion turbulente + moment angulaire
	CALL difft(l,melange,p,t,r,ro,drox,kap,dkapro,dkapx,deff,gradad,gradrad,
	1 nu32,y,d,dd)

c	PRINT*,'p,t,r,nu32' ; WRITE(*,2000)p,t,r,nu32
c	PRINT*,'v' ; WRITE(*,2000)v
c	DO i=1,nchim
c	 PRINT*,'pour ',nom_elem(i)
c	 WRITE(*,2000)d(i,:)
c	ENDDO
c	PAUSE'eq_diff_chim'
	
c les réactions thermonucléaires Xdot
	CALL nuc(t,ro,ychim,dcomp,jac,.TRUE.,2,
	1 epsilon,depst,depsro,depsx,hhe,be7e,b8e,n13e,o15e,f17e)

c s'il y a perte/gain de masse
	IF(nu >= nuzc_ext)THEN
	 IF(melange .AND. l_vent)CALL vent(ychim,dcomp,jac)
	 IF(melange .AND. l_planet)CALL planetoides(ychim,dcomp,jac)
	ENDIF	
	 
c formation des équations
	ad=0.d0 ; as=0.d0 ; bd=0.d0 ; bs=0.d0
	 
c < 3/2 nu^1/2 [ Xi(t+1)-Xi(t) - Psi dt ] , Fi_j > +
c + < 32/(3 nu^1/2) pi^2 r^4 ro^2 dt
c [Di1 dX1/dnu+...+Din dXn/dnu],dFi_j/dnu >
c  - < 4 pi r^2 ro dt Vi Xi , dFi_j/dnu >=0
	cte1=cte1_0*SQRT(nu)
	cte2=cte2_0*r**4*ro**2/SQRT(nu)*dt
	cte3=cte3_0*r**2*ro*dt
	drox=drox*ah/ro		!d ro / d x_1 à ro pres
	dcte2x=cte2*2.d0*drox ; dcte3x=cte3*drox	

c pour bs(spline):	3/2 nu^1/2 [ Xi(t+1)-Xi(t) - Psi dt ]
c pour bd(dérivée de la spline): 32/(3 nu^1/2) pi^2 r^4 ro^2 dt *
c			*  (Dt+Di+Dm) dXi/dnu  -4 pi r^2 ro dt Vi Xi

c bs : coefficient de la spline du produit scalaire
c bd : coefficient de la dérivée de la spline du produit scalaire dérivées
c y(var,der) := dérivée d'ordre der=0,1
	 
c < bs . S > = 0
 
c bs=cte1*(y(:,0)-comp_t-dt*dcomp)
	 bs(1:nchim)=cte1*(y(1:nchim,0)-comp_t(1:nchim)-dt*dcomp(1:nchim))

c as(i,j,k)=coefficient de la k-ieme dérivée
c de la j-ieme variable de la i-ieme équation
	DO i=1,nchim
	 as(i,i,0)=1.d0		!der / x_i à cte1 près
	 DO j=1,nchim	
	  as(i,j,0)=as(i,j,0)-dt*jac(i,j)	!der / x_j jacobien
	 ENDDO
	 as(i,:,0)=as(i,:,0)*cte1	  
	ENDDO
	
c < bd . dS > = 0

c bd=-cte3*v*y(:,0) !coeff. de Xi: vitesse
	bd(1:nchim)=-cte3*v(1:nchim)*y(1:nchim,0) !coeff. de Xi: vitesse 
	DO j=1,nchim		!produit scalaire D . dx/dm	
	 bd(1:nchim)=bd(1:nchim)+cte2*d(1:nchim,j)*y(j,1)	  
	ENDDO	!j

c ad(i,j,k)=coefficient de la k-ieme dérivée
c de la j-ieme variable de la i-ieme équation <.> d/dx Ni
c dérivée dans bd de cte3 v_i x_i + SUM d_ij x_j' / x_k
c ad(:,1,0)=-dcte3x*v*y(:,0)  !der / x_1
	ad(1:nchim,1,0)=-dcte3x*v(1:nchim)*y(1:nchim,0)  !der / x_1	 
	DO i=1,nchim	  
	 ad(i,i,0)=ad(i,i,0)-cte3*v(i)	!der / x_i
	 DO j=1,nchim
	  ad(i,j,0)=ad(i,j,0)-cte3*dv(i,j)*y(i,0)	!der / x_i	  
	 ENDDO
	 DO j=1,nchim
	  ad(i,j,1)=cte2*d(i,j)	!der / x_i'
	  ad(i,1,0)=ad(i,1,0)+dcte2x*d(i,j)*y(j,1)	!der / x_1
	  DO k=1,nchim	   
	   ad(i,k,0)=ad(i,k,0)+cte2*dd(i,j,k)*y(j,1)	!der / x_k
	  ENDDO
	 ENDDO
	ENDDO

c vérification des dérivées
c	 deriv=.TRUE.
c	 deriv=.NOT.melange
c	 deriv=fait == 2
	IF(deriv)THEN
	 PRINT*,'nchim=',nchim ; PRINT*,'p,t,r,l,nu32,nu,ro,cte2'	
	 WRITE(*,2000)p,t,r,l,nu32,nu,ro,cte2
	 PRINT*,'xchim' ; WRITE(*,2000)y(:,0)
	 PRINT*,'dxchim' ; WRITE(*,2000)y(:,1)
	 PRINT*,'bs' ; WRITE(*,2000)bs
	 PRINT*,'bd' ; WRITE(*,2000)bd
c	 les équations et les dérivées numériques en Xi
	  	 
	 ALLOCATE(bss(nchim),bds(nchim),yd(nchim,0:1),
	1 jac0(nchim,nchim),dcomp0(nchim))
	 yd=y ; jac0=jac ; dcomp0=dcomp
	 DO kd=0,1		!dérivées 0 et '
	  DO k=1,nchim
	   stor0=yd(k,kd) ; stor=stor0*unpdx ; IF(stor == 0.d0)stor=dx
	   yd(k,kd)=stor ; dstor=stor-stor0

c xchim composition chimique par gramme, ychim composition chimique par mole
	   ychim(:)=yd(:,0) ; xchim=ychim
c	   CALL chim_gram(xchim,dcomp)	!dcomp VT

c la densité et le gradient adiabatique
c	    PRINT*,k ; WRITE(*,2000)xchim ; PAUSE'ici'

c	    CALL etat(p,t,xchim,.FALSE.,
c	1   ro,drop,drot,drox,u,dup,dut,dux,
c	2   delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
c	3   gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)

c	    l'opacité et le gradient radiatif

c	    CALL opa(xchim,t,ro,kap,dkapt,dkapro,dkapx)
c	    gradrad=cte4_0*l*p*kap/nu32/t**4
	 
c le coefficient de diffusion

	   v=0.d0 ; dv=0.d0 ; d=0.d0 ; dd=0.d0

c le coefficient de diffusion microscopique
	   IF(.NOT.melange)CALL diffm(p,t,r,l,nu32,ro,drox,kap,dkapx,w,
	1   gradrad,dgradradx,ychim,d,dd,v,dv)

c le coefficient de diffusion turbulente + moment angulaire
	   CALL difft(l,melange,p,t,r,ro,drox,kap,dkapro,dkapx,deff,gradad,
	1  gradrad,nu32,y,d,dd)

c les réactions thermonucléaires Xdot
	   CALL nuc(t,ro,ychim,dcomp,jac,.FALSE.,2,
	1   epsilon,depst,depsro,depsx,hhe,be7e,b8e,n13e,o15e,f17e)

c	    WRITE(*,2000)ro,d(1:nchim,1),dcomp(1:nchim)

c s'il y a perte/gain de masse
	   IF(melange .AND. l_vent)CALL vent(ychim,dcomp,jac)
	   IF(melange .AND. l_planet)CALL planetoides(ychim,dcomp,jac)
	   bss=0.d0 ; bds=0.d0

	   DO i=1,nchim
	    bss(i)=cte1*(yd(i,0)-comp_t(i)-dt*dcomp(i))
	    bds(i)=-cte3*v(i)*yd(i,0)	!coeff. de Xi: vitesse
	    DO j=1,nchim		!produit scalaire D . dx/dm
	     bds(i)=bds(i)+cte2*d(i,j)*yd(j,1)
	    ENDDO	!j
	   ENDDO	!i
	     
	   SELECT CASE(kd)
	   CASE(0)	     
	    WRITE(*,1)nom_elem(k)
1	    FORMAT('dérivée /',a5)	      
	    PRINT*,'dérivée jacobien, analytique - numérique'
	    WRITE(*,2000)(jac0(i,k),(dcomp(i)-dcomp0(i))/dstor,i=1,nchim)     
	   CASE(1)
	    WRITE(*,2)nom_elem(k)	     
2	    FORMAT('dérivée / d',a5)
	   END SELECT	     
	   PRINT*,'dérivée as, analytique - numérique'
	   WRITE(*,2000)(as(i,k,kd),(bss(i)-bs(i))/dstor,i=1,nchim)
	   PRINT*,'dérivée ad, analytique - numérique'
	   WRITE(*,2000)(ad(i,k,kd),(bds(i)-bd(i))/dstor,i=1,nchim)  
	   PRINT*
	   yd(k,kd)=stor0
	  ENDDO	!k
	 ENDDO	!kd
	 DEALLOCATE(bss,bds,yd,jac0,dcomp0)
	 PAUSE'point suivant'
	ENDIF		!test
	 
c------------------fin des équations--------------------------

	RETURN

	END SUBROUTINE eq_diff_chim
