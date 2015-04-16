
c---------------------------------------------------------------
	
	SUBROUTINE difft_smc(melange,deff,lum,m,p,r,t,y,d,dd)
	
c routine PRIVATE du module mod_evol
c formation du coefficient de diffusion turbulente en tenant compte de la
c semi-convection Langer et al. A&A 126,207, 1983
c d_melange + d_turb + nu_rad + Deff

c nu_rad, suivant Morel & Thévenin évite la sédimentation de l'hélium
c sauf dans les ZC mélangées

c entrées
c	melange
c	deff : diffusion turbulente de rotation
c	lum : luminosité/Lsol
c	m : masse m/Msol
c	p : pression
c	r : rayon/Rsol
c	t : température K
c	y(nchim,0:1) : composition chimique, dérivée / m^2/3

c sorties
c	d, dd : coefficients D_ij de diffsion de l'élément x_i dû à x_j,
c	 dérivées / x_k
c	Dimensions et initialisations dans le programme appelant
c	d(nchim,nchim), dd(nchim,nchim,nchim)

c les coefficient grad_ad, grad_rad, ro etc... sont recalculés

c Auteur: P.Morel, OCA, CESAM2k

c-------------------------------------------------------------------------

	USE mod_conv, ONLY : conv
	USE mod_donnees, ONLY : alpha, aradia, clight, d_conv, d_turb, g, Krot,
	1 langue, ledoux, lsol, msol, nchim, nrot, ord_rot, pi, nucleo,
	2 re_nu, rsol
	USE mod_etat, ONLY : etat, mu_mol
	USE mod_kind	
	USE mod_numerique, ONLY : bsp1dn, no_croiss
	USE mod_opa, ONLY : opa
	USE mod_variables, ONLY : chim_gram, knotc, knotr, mc, mct,
	1 mrot, mrott, n_ch, n_rot, rota, wrot
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(nchim,0:1) :: y		
	REAL (kind=dp), INTENT(in) :: deff, lum, m, p, t, r
	LOGICAL, INTENT(in) :: melange	
	REAL (kind=dp), INTENT(inout), DIMENSION(nchim,nchim,nchim) :: dd
	REAL (kind=dp), INTENT(inout), DIMENSION(nchim,nchim) :: d
		
	REAL (kind=dp), DIMENSION(nchim) :: dgrad_mux,
	1 dlnmu_x, dmu_x, dxchim, dxchimm, xchim, xchimm
	REAL (kind=dp), DIMENSION(nrot) :: dfrot, frot	
	REAL (kind=dp), SAVE :: alfa_l, cte1, cte2, cte4, cte8, cte13
	REAL (kind=dp) :: alfa, beta, cp, delta, deltap, deltat, deltax, dconv,
	1 dconvx, dcpp, dcpt, dcpx, dgradadp, dgradadt, dgradadx, dgrad_ldx, 
	2 dgrad_radx, dhpx, dkapt, dkapro,  dkapx, dkradx, dlnmu,
	3 dnu_radx, drop, drot, drox,  dup, dut, dux,
	4 gamma1, grad_ad, grad_ldx, grad_mu, grad_rad,
	5 grav, hp, kap, krad, ldx, mu, nel, nu, nu_rad, ro,  u, w, Zbar

	INTEGER, SAVE :: l=1
	INTEGER :: i

	LOGICAL, SAVE :: init=.TRUE.

c--------------------------------------------------------------------------

2000	FORMAT(8es10.3)

c------------------initialisations----------------
	IF(init)THEN
	 init=.FALSE.
	 alfa_l=1.d0/6.d0	 
	 cte1=4.d0/3.d0*aradia*clight
	 cte2=2.d0/3.d0*rsol
	 cte4=re_nu*aradia/clight*4.d0/15.d0	  
	 cte8=lsol/4.d0/pi/rsol/rsol
	 cte13=g*msol/rsol/rsol
	 
	 IF(ledoux)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1001) ; WRITE(2,1001)
1001	   FORMAT('Ledoux=.TRUE. is unconsistent with NOM_DIFFT=difft_smc',/,
	1  'modify the input file, STOP')	 
	  CASE DEFAULT	 
	   WRITE(*,1) ; WRITE(2,1)
1	   FORMAT('Ledoux=.TRUE. est incompatible avec NOM_DIFFT=difft_smc',/,
	1  'modifier le fichier de données, ARRET')
	   STOP
	  END SELECT
	 ENDIF
	 	 
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1010)d_turb,re_nu
	  WRITE(2,1010)d_turb,re_nu
1010	  FORMAT('Convective mixing, semi-convective, turbulent : dturb=',
	1 es10.3,', radiative : Re_nu=',es10.3,' + Deff with rotation')	 
	 CASE DEFAULT	 
	  WRITE(*,10)d_turb,re_nu ; WRITE(2,10)d_turb,re_nu
10	  FORMAT('Mélange convectif, semi-convectif, turbulent : dturb=',
	1 es10.3,', radiatif : Re_nu=',es10.3,' + Deff avec rotation')
	 END SELECT
	ENDIF

c---------------diffusivité radiative------------------------

c nu=m^2/3
	nu=m**(2.d0/3.d0)

c composition chimique / mole et gramme	
	xchim(:)=y(:,0) ; dxchim(:)=y(:,1)
	xchimm=xchim ; dxchimm =dxchim 
	CALL chim_gram(xchimm,dxchimm)		

c pression gazeuse (Pgas + Prad)
	CALL etat(p,t,xchimm,.TRUE.,
	1 ro,drop,drot,drox,u,dup,dut,dux,
	2 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3 grad_ad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)

	drox=drox*nucleo(1)	!par mole
	dcpx=dcpx*nucleo(1)
	dgradadx=dgradadx*nucleo(1)
	
	CALL opa(xchimm,t,ro,kap,dkapt,dkapro,dkapx)
	dkapx=dkapro*drox+dkapx*nucleo(1)	!par mole

c diffusivité radiative	
	krad=cte1/kap/ro*t**3 ;	dkradx=-krad*(drox/ro+dkapx/kap)
	
c coefficient de viscosité radiative, dérivée par mole
	nu_rad=cte4*t**4/kap/ro**2 ; dnu_radx=-nu_rad*(dkapx/kap+2.d0*drox/ro)

c---------------------les gradients-------------------------
	
c vitesse angulaire
	SELECT CASE(Krot)
	CASE(0,1,2)
	 w=wrot
	CASE(3,4,5)
	 IF(m == 0.d0)THEN
	  w=rota(1,1)
	 ELSE
	  CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,
	1 MIN(nu,mrot(n_rot)),l,frot,dfrot)
	  IF(no_croiss)PRINT*,'Pb. en 1 dans thermo'	
	  w=frot(1)
	 ENDIF	 
	END SELECT	
	
c gravité effective
	grav=cte13*m/r**2-cte2*w**2*r
	 
c échelle de hauteur de pression
	hp=p/grav/ro ; dhpx=-hp*drox/ro
	  
c gradient radiatif	 
	grad_rad=cte8*lum*hp/r**2/krad/t
	dgrad_radx=grad_rad*(dhpx/hp-dkradx/krad)

c dans une zone	mélangée	
	IF(melange)THEN
	
c on peut être radiatif dans une extension d'une zone convective (overshoot)
	 IF(grad_rad <= grad_ad)THEN	
	  dconv=d_conv ; dconvx=0.d0
c	  PRINT*,'partie radiative avec mélange'

c zone de mélange
	 ELSE		 
	
c calcul de grad_mu		
	  CALL mu_mol(dxchim,hp,nu,r,ro,t,xchim,dlnmu,dlnmu_x,grad_mu,dgrad_mux,
	1 mu,dmu_x,nel,Zbar)	
		
c gradient de Ledoux
	  ldx=beta/(4.d0-3.d0*beta)*grad_mu
	  grad_ldx=grad_ad+ldx
	  dgrad_ldx=dgradadx+beta/(4.d0-3.d0*beta)*dgrad_mux(1)
	  
c----semi-convection : grad_ad < grad_rad < grad_ldx ---------

	  IF(MAX(grad_ad,MIN(grad_rad,grad_ldx)) == grad_rad)THEN 
	   dconv=alfa_l*krad/cp/ro*(grad_rad-grad_ad)/(grad_ldx-grad_rad)
	   IF(dconv >= d_conv)THEN
	    dconv=d_conv ; dconvx=0.d0
	   ELSE
	    dconvx=dconv*(dkradx/krad-dcpx/cp-drox/ro+(dgrad_radx-dgradadx)/
	1  (grad_rad-grad_ad)-(dgrad_ldx-dgrad_radx)/(grad_ldx-grad_rad))
	   ENDIF	   
	
c convection
	  ELSE
	   dconv=d_conv ; dconvx=0.d0
	  ENDIF	   	 
	 ENDIF
	
c zone radiative sans mélange	
	ELSE
	 dconv=0.d0 ; dconvx=0.d0
c	 PRINT*,'radiatif sans mélange'
	ENDIF	
		 
c contributions des diverses diffusivités turbulentes
	DO i=1,nchim
	 d(i,i)=d(i,i)+dconv+d_turb+nu_rad+deff
	 dd(i,i,1)=dd(i,i,1)+dconvx+dnu_radx
	ENDDO
	 
	RETURN 
	 
	END SUBROUTINE difft_smc
