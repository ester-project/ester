	
c***********************************************************************	

	SUBROUTINE coeff_vth(n,xx,y)	 

c subroutine PRIVATE du module mod_evol

c détermination des coefficients pour l'approximation vth
c il y a nvth=10 fonctions, interpolation sur la base mc, mct, knotc
c y(1)=lnP, y(2)=lnT, y(3)=r**2, y(4)=l**2/3, y(5)=ln ro, 
c y(6)=cp, y(7)=delta, y(8)=gamma1, y(9)=ln µ, y(10)=ln kappa

c entrées:
c	n : nombre de variables nvth=10
c	xx : abscisse en m^2/3

c sorties
c	y(n)

c Auteur: P.Morel, Département Cassiopée, O.C.A., CESAM2k

c----------------------------------------------------------------------------

	USE mod_donnees, ONLY : diffusion, en_masse, lisse, mu_saha, m_ch,
	1 nchim, ne, zi
	USE mod_etat, ONLY : etat, saha
	USE mod_kind	
	USE mod_numerique, ONLY : bsp1dn
	USE mod_opa, ONLY : opa
	USE mod_variables, ONLY : bp, chim, chim_gram, inter, knot, knotc,
	1 mc, mct, m23, n_ch, n_qs, q, qt, r2, vth
	
	IMPLICIT NONE
	
	REAL(kind=dp), INTENT(in) :: xx	   
	INTEGER, INTENT(in) :: n
	REAL(kind=dp), INTENT(out), DIMENSION(n) :: y

	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: ioni
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: z_bar		
	REAL (kind=dp), DIMENSION(nchim) :: dxchim, dxchimg, xchim, xchimg
	REAL (kind=dp), DIMENSION(ne) :: dfdx, f	
	REAL (kind=dp) :: alfa, beta, cp, dcpp, dcpt, dcpx, delta,
	1 deltap, deltat, deltax, drop, drot, drox, dup, dut, dux, eta,	
	2 gamma1, gradad, dgradadp, dgradadt, dgradadx, dkapro, dkapt, dkapx,
	3 kap, nel, p, ro, t, u, xxx
	
	INTEGER, SAVE :: ll

	LOGICAL, SAVE :: init=.TRUE.

c----------------------------------------------------------------------------

2000	FORMAT(8es10.3)

c initialisations
	IF(init)THEN
	 init=.FALSE.
	 ALLOCATE(z_bar(nchim))
	 IF(mu_saha)THEN
	  ALLOCATE(ioni(0:NINT(MAXVAL(zi)),nchim))
	 ELSE
	  z_bar=zi
	 ENDIF	
	ENDIF
	
c m ou m**2/3 pour inter	
	IF(en_masse)THEN
	 xxx=xx			!m**2/3
	ELSE
	 xxx=SQRT(xx**3)	!m
	ENDIF	
		
c variables de structure aux xx		
	CALL inter('m23',bp,q,qt,n_qs,knot,xxx,f,dfdx,r2,m23)
	y(1:2)=f(1:2)	!lnP, lnT
	p=EXP(y(1)) ; t=EXP(y(2))
	IF(en_masse)THEN		
	 y(3:4)=f(3:4)	!r**2, l**2/3
	ELSE
	 y(3)=f(3)**2 ; y(4)=(f(4)**2)**(1.d0/3.d0)	!r**2, l**2/3
	ENDIF	

c composition chimique	 
	CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1 knotc,.TRUE.,xx,ll,xchim,dxchim)

c densité	 
	xchimg=xchim ; dxchimg=dxchim ; CALL chim_gram(xchimg,dxchimg)	 
	CALL etat(p,t,xchimg,.FALSE.,
	1 ro,drop,drot,drox,u,dup,dut,dux,
	2 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)	
	y(5)=LOG(ro)

c y(6)=ln cp
	y(6)=cp

c y(7)=ln delta			
	y(7)=delta	

c y(8)=ln gamma1			
	y(8)=gamma1	
	 	 
c calcul de ln mu, ionisation partielle	 
	IF(mu_saha)CALL saha(xchim,t,ro,ioni,z_bar,nel,eta)

c y(9)=ln mu	 
	y(9)=-LOG(DOT_PRODUCT((1.d0+z_bar),xchim))
	
c opacité	
	CALL opa(xchimg,t,ro,kap,dkapt,dkapro,dkapx)
	
c y(10)=ln kappa	 
	y(10)=LOG(kap)
		
	RETURN
	
	END SUBROUTINE coeff_vth

	
