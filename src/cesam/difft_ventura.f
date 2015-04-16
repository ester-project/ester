
	SUBROUTINE difft_ventura(melange,nu,t,ro,drox,kap,dkapx,deff,d,dd)

!	routine private du module mod_evol

!	calcul du coefficient de diffusion turbulente, d_turb + nu_rad
!	et du Moment Angulaire, sauf dans les ZC mélangées

!	Dimensions et initialisations dans le programme appelant
!	d(nbelem,nbelem), dd(nbelem,nbelem,nbelem), v(nbelem),
!	dv(nbelem,nbelem)

!	convention de notation :
!	équation de diffusion dXi/dt=dFi/dm + nuclear, i=1,nbelem
!	Fi=4 pi r**2 ro (4 pi r**2 ro D.dX/dm - Vi Xi)

!	d=D=(di1, di2,... din) avec Dij coefficient de d Xj / dm
!	dans le produit scalaire D.dX/dm=sum d_ij d Xj / dm

!	pour ligne d'indice i
!	v(i) coefficient de x_i,
!	dv(i,k)=dv(nbelem*(k-1)+i)=dérivée v_i / x_k
!	seule la première colonne de dv
!	est non nulle (pas de dérivées / Xi, i .ne. 1)
!	d(i,j)=coefficient d_ij de d x_j / dm
!	dd(i,j,k)= dérivée de d_ij / x_k

!	Auteur: P.Morel, OCA
!	CESAM2k

!       entrées
!	melange=.TRUE.: on est dans une ZC
!	p, t, r, l, m, ro: données au point de calcul
!	xi: composition chimique, par mole
!	kap: opacité 
!	gradad, gradrad: gradients
!	terminaisons x : dérivées/ X1 (ie H) par gramme
!	mstar: masse avec perte de masse
!	m_zc, r_zc, lim : masses, rayons, nombre de limites ZR/ZC
!	age, gamma1, cp, delta: notations evidentes

!       sorties
!	d0, dd : coefficients d_ij de d x_j / d m et dérivées / x_k
!	v0, dv : coefficients v_i de x_i et dérivées / x_k

!-------------------------------------------------------------------------

	USE mod_donnees, ONLY : aradia, clight, d_turb, en_masse, Ipg, Krot, ne, 
	1    nchim, pturb, re_nu, m_ch, rsol, msol, g, beta_v, zeta_v, frac_vp, ord_qs
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn, no_croiss
	USE mod_variables, ONLY: r_zc_conv, hp_zc_conv, pr_zc_conv, id_conv, 
	1    u_zc_conv, chim_t, bp, chim_gram,jlim,tot_conv, n_qs, m23, q, qt, r2, 
	2    knot, mc, n_ch, chim, mct, knotc, inter
	USE mod_etat, ONLY: etat

	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: nu, deff, kap, dkapx, t
	REAL (kind=dp), INTENT(inout) :: drox, ro
	LOGICAL, INTENT(in) :: melange
	REAL (kind=dp), INTENT(inout), DIMENSION(:,:,:) :: dd
	REAL (kind=dp), INTENT(inout), DIMENSION(:,:) :: d

	REAL (kind=dp), SAVE :: cte2, d_conv, w_diff, cte1, cte13
	REAL (kind=dp) :: dnu_radx, nu_rad, l, u, ub, pb, zup, zlow, hpup,
	1    krad, grav, hp, fthick, pc
	REAL (kind=dp), DIMENSION(2) :: ui, pri, tei, mai, rai, lui, li, hpi, roi
	REAL (kind=dp), DIMENSION(ne) :: fqs, dfqs
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: ychim, dychim, uint
	REAL (kind=dp) :: ro1,drop1,drot1,drox1,uu1,duup1,duut1,duux1,ra,ma,
	1    delta1,deltap1,deltat1,deltax1,cp1,dcpp1,dcpt1,dcpx1,bid,lu,pr,te,
	2    gradad1,dgradadp1,dgradadt1,dgradadx1,alfa1,beta_v1,gamma11,drop,drot,
	3    uu,duup,duut,duux,u_conv,
	4    delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	5    gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:), SAVE :: xchim, dxchim, 
	1    xchimi, dxchimi

	INTEGER :: i, j, j1, j2
	INTEGER, SAVE :: ll, lc

	LOGICAL, SAVE :: init=.TRUE.

!--------------------------------------------------------------------------

 2000	FORMAT(8es10.3)

	IF(init)THEN
	   init = .FALSE.
	   PRINT* ; WRITE(2,*)
	   cte1 = 4.d0/3.d0*aradia*clight
	   cte13 = g*msol/rsol/rsol
	   cte2 = re_nu*aradia/clight*4.d0/15.d0 ; d_conv=1.d13

	   ALLOCATE(xchim(nchim),dxchim(nchim),xchimi(nchim),dxchimi(nchim))
	   WRITE(*,10)d_conv,d_turb,re_nu ; WRITE(2,10)d_conv,d_turb,re_nu
 10	   FORMAT('Diffusion turbulente: Dturb + Dradiative',/,
	1	'Dans ZC, Dconv=',es10.3,'. Dans ZR, Dturb=',es10.3,
	2	', Re_nu=',es10.3)
	   IF(Krot > 0)THEN
	      w_diff=0.d0
!       PRINT*,'entrer w_diff' ; read*,w_diff
	      WRITE(*,11)w_diff ; WRITE(2,11)w_diff  
 11	      FORMAT('avec diffusion constante du moment angulaire w_diff =',
	1	   es10.3,//)
	   ENDIF
	ENDIF

	IF(en_masse)THEN
	   bid = MIN(nu,m23(n_qs))
	   ma = nu**(1.5d0)
	ELSE
	   bid = MIN(nu**1.5d0,m23(n_qs))
	   ma = nu**(1.5d0)
	ENDIF	

	CALL inter('m23',bp,q,qt,n_qs,knot,bid,fqs,dfqs,r2,m23)

c	variables quasi-statiques au point courant

	IF(pturb)THEN		!avec pression turbulente
	   pr = EXP(fqs(Ipg))	 
	ELSE			!sans pression turbulente
	   pr = EXP(fqs(1))	 
	ENDIF
	te = EXP(fqs(2))
	IF(en_masse)THEN
	   ra = SQRT(ABS(fqs(3)))
	   lu = SQRT(ABS(fqs(4)))**3
	ELSE
	   ra = ABS(fqs(3))
	   lu = fqs(4)
	ENDIF

	IF(en_masse)THEN
	   bid = MAX(mc(1),MIN(nu,mc(n_ch)))
	ELSE
	   bid = MAX(mc(1),MIN(nu**1.5d0,mc(n_ch)))
	ENDIF	
	CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,knotc,.TRUE.,bid,ll,xchim,dxchim)
	CALL chim_gram(xchim,dxchim)
	CALL etat(pr,te,xchim,.FALSE.,ro,drop,drot,drox,uu,duup,duut,duux,
	1    delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2    gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)

	grav=cte13*ma/ra**2	! sem wrot
	      
!       rot=-cte2*w**2*r	!gravite effective avec rotation
!       grav=grav+rot	!remarque de N. Audard
	      
!       echelle de hauteur de pression
	      
	hp = pr/grav/ro		!echelle de hauteur de pression
	
	IF(melange)THEN ! we are in a convective zone
	   IF(id_conv==0 .AND. ra<=r_zc_conv(1)) THEN ! central CZ
	      pc = pr_zc_conv(0)
	      pb = pr_zc_conv(1)
	      hpup = hp_zc_conv(1)
	      zup = (r_zc_conv(1)-ra)*rsol
	      zlow = ra*rsol
	      l = (zup + beta_v*hpup)*zlow/(zup + zlow) ! mixing length

	      IF((pr>(1.+frac_vp)*pb .AND. pr<(1.-frac_vp)*pc) .OR. jlim(1)<5) THEN
!       if we are not near the border of the CZ, and the CZ is not too small...
		 IF(pr==pc .OR. pr==pb) THEN
		    u = 0.d0
		 ELSE
		    CALL v2_cgm(pr,te,lu,ma,ra,l,xchim,u)
		    u = SQRT(u)	! turbulent velocity
		 ENDIF

		 DO i=1,nchim
		    d(i,i) = u*l/3.d0 ! diffusion coefficient
		 ENDDO

	      ELSEIF(pr<(1.d0+frac_vp)*pb) THEN
!       if we are near the top of the CZ, we extrapolate the conv. velocities
!       otherwise, u goes to zero as we approach the border of the CZ
		 DO j=jlim(1),1,-1 ! we go inward from the top of the CZ
c$$$		    CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(j),lc,
c$$$	1		 fqs,dfqs)
		    CALL inter('m23',bp,q,qt,n_qs,knot,mc(j),fqs,dfqs,r2,m23)

		    pri(1) = EXP(fqs(1))
		    IF(pri(1)>(1.+frac_vp)*pb) THEN ! we passed the critical pr.
		       j1 = j
		       j2 = j + 1	! the 2 points for extrapolation.
		       IF(pri(2)>=pc) THEN
			  j1 = 1
			  j2 = 2
		       ENDIF

c$$$		       CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(j1),lc,
c$$$	1		    fqs,dfqs)
		       CALL inter('m23',bp,q,qt,n_qs,knot,mc(j1),fqs,dfqs,r2,m23)
		       pri(1) = EXP(fqs(1))
		       tei(1) = EXP(fqs(2))
		       lui(1) = fqs(4)**1.5
		       rai(1) = SQRT(fqs(3))
		       mai(1) = fqs(5)**1.5
		       IF(en_masse)THEN
			  bid = MAX(mc(1),MIN(fqs(5),mc(n_ch)))
		       ELSE
			  bid = MAX(mc(1),MIN(fqs(5)**1.5d0,mc(n_ch)))
		       ENDIF	
		       CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,knotc,.TRUE.,bid,ll,
	1		    xchimi,dxchimi)
		       CALL chim_gram(xchimi,dxchimi)

		       CALL etat(pri(1),tei(1),xchimi,.FALSE.,roi(1),drop1,drot1,
	1		    drox1,uu1,duup1,duut1,duux1,delta1,deltap1,deltat1,
	2		    deltax1,cp1,dcpp1,dcpt1,dcpx1,gradad1,dgradadp1,
	3		    dgradadt1,dgradadx1,alfa1,beta_v1,gamma11)

		       hpup = hp_zc_conv(1) !; PRINT*,'hpup =',hpup
		       zup = (r_zc_conv(1) - rai(1))*rsol
		       zlow = rai(1)*rsol
		       li(1) = (zup + beta_v*hpup)*zlow/(zup + zlow) ! mixing length
!       PRINT*,'j1, j2', j1, j2
		       CALL v2_cgm(pri(1),tei(1),lui(1),mai(1),rai(1),li(1),xchimi,
	1		    ui(1))
		       ui(1) = SQRT(ui(1))

c$$$		       CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(j2),lc,
c$$$	1		    fqs,dfqs)
		       CALL inter('m23',bp,q,qt,n_qs,knot,mc(j2),fqs,dfqs,r2,m23)

		       pri(2) = EXP(fqs(1))
		       tei(2) = EXP(fqs(2))
		       lui(2) = fqs(4)**1.5
		       rai(2) = SQRT(fqs(3))
		       mai(2) = fqs(5)**1.5

		       IF(en_masse)THEN
			  bid = MAX(mc(1),MIN(fqs(5),mc(n_ch)))
		       ELSE
			  bid = MAX(mc(1),MIN(fqs(5)**1.5d0,mc(n_ch)))
		       ENDIF	
		       CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,knotc,.TRUE.,bid,ll,
	1		    xchimi,dxchimi)
		       CALL chim_gram(xchimi,dxchimi)

		       CALL etat(pri(2),tei(2),xchimi,.FALSE.,roi(2),drop1,drot1,
	1		    drox1,uu1,duup1,duut1,duux1,delta1,deltap1,deltat1,
	2		    deltax1,cp1,dcpp1,dcpt1,dcpx1,gradad1,dgradadp1,
	3		    dgradadt1,dgradadx1,alfa1,beta_v1,gamma11)

		       zup = (r_zc_conv(1) - rai(2))*rsol
		       zlow = rai(2)*rsol
		       li(2) = (zup + beta_v*hpup)*zlow/(zup + zlow) ! mixing length
!       PRINT*,'j1, j2', j1, j2
		       CALL v2_cgm(pri(2),tei(2),lui(2),mai(2),rai(2),li(2),xchimi,
	1		    ui(2))
		       ui(2) = SQRT(ui(2))

!       extrapolated velocity
		       ui = LOG(ui) ; pri = LOG(pri)
		       u = ui(1) + (ui(2) - ui(1))/(pri(2) - pri(1))*
	1		    (LOG(pr) - pri(1))
		       u = EXP(u)

c$$$		       print"(20es22.14)",EXP(ui(1)),EXP(ui(2)),u,
c$$$	1		    EXP(pri(1)),EXP(pri(2)),pr

		       DO i=1,nchim
			  d(i,i) = u*l/3.d0	! diffusion coefficient
		       ENDDO

		       EXIT ! since we reached the desired point
		    ENDIF
		 ENDDO ! going inwards from the top of the CZ

	      ELSE
!       if we are near the bottom of the CZ, we extrapolate the conv. velocities,
!       otherwise, u goes to zero as we approach the border of the CZ
		 DO j=1,jlim(1) ! we go outwards from the center
		    CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(j),lc,
	1		 fqs,dfqs)

		    pri(1) = EXP(fqs(1))
		    IF(pri(1)<(1.d0-frac_vp)*pc) THEN ! we passed the critical pr.
		       j1 = j - 1
		       j2 = j	! the 2 points for extrapolation.
		       IF(pri(1)<=pb) THEN
			  j1 = j - 2
			  j2 = j - 1
		       ENDIF

c$$$		       CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(j1),lc,
c$$$	1		    fqs,dfqs)
		       CALL inter('m23',bp,q,qt,n_qs,knot,mc(j1),fqs,dfqs,r2,m23)
		       pri(1) = EXP(fqs(1))
		       tei(1) = EXP(fqs(2))
		       lui(1) = fqs(4)**1.5
		       rai(1) = SQRT(fqs(3))
		       mai(1) = fqs(5)**1.5
		       IF(en_masse)THEN
			  bid = MAX(mc(1),MIN(fqs(5),mc(n_ch)))
		       ELSE
			  bid = MAX(mc(1),MIN(fqs(5)**1.5d0,mc(n_ch)))
		       ENDIF	
		       CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,knotc,.TRUE.,bid,ll,
	1		    xchimi,dxchimi)
		       CALL chim_gram(xchimi,dxchimi)

		       CALL etat(pri(1),tei(1),xchimi,.FALSE.,roi(1),drop1,drot1,
	1		    drox1,uu1,duup1,duut1,duux1,delta1,deltap1,deltat1,
	2		    deltax1,cp1,dcpp1,dcpt1,dcpx1,gradad1,dgradadp1,
	3		    dgradadt1,dgradadx1,alfa1,beta_v1,gamma11)

		       hpup = hp_zc_conv(1) !; PRINT*,'hpup =',hpup
		       zup = (r_zc_conv(1) - rai(1))*rsol
		       zlow = rai(1)*rsol
		       li(1) = (zup + beta_v*hpup)*zlow/(zup + zlow) ! mixing length
!       PRINT*,'j1, j2', j1, j2
		       CALL v2_cgm(pri(1),tei(1),lui(1),mai(1),rai(1),li(1),xchimi,
	1		    ui(1))
		       ui(1) = SQRT(ui(1))

c$$$		       CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(j2),lc,
c$$$	1		    fqs,dfqs)
		       CALL inter('m23',bp,q,qt,n_qs,knot,mc(j2),fqs,dfqs,r2,m23)

		       pri(2) = EXP(fqs(1))
		       tei(2) = EXP(fqs(2))
		       lui(2) = fqs(4)**1.5
		       rai(2) = SQRT(fqs(3))
		       mai(2) = fqs(5)**1.5

		       IF(en_masse)THEN
			  bid = MAX(mc(1),MIN(fqs(5),mc(n_ch)))
		       ELSE
			  bid = MAX(mc(1),MIN(fqs(5)**1.5d0,mc(n_ch)))
		       ENDIF	
		       CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,knotc,.TRUE.,bid,ll,
	1		    xchimi,dxchimi)
		       CALL chim_gram(xchimi,dxchimi)

		       CALL etat(pri(2),tei(2),xchimi,.FALSE.,roi(2),drop1,drot1,
	1		    drox1,uu1,duup1,duut1,duux1,delta1,deltap1,deltat1,
	2		    deltax1,cp1,dcpp1,dcpt1,dcpx1,gradad1,dgradadp1,
	3		    dgradadt1,dgradadx1,alfa1,beta_v1,gamma11)

		       zup = (r_zc_conv(1) - rai(2))*rsol
		       zlow = rai(2)*rsol
		       li(2) = (zup + beta_v*hpup)*zlow/(zup + zlow) ! mixing length
!       PRINT*,'j1, j2', j1, j2
		       CALL v2_cgm(pri(2),tei(2),lui(2),mai(2),rai(2),li(2),xchimi,
	1		    ui(2))
		       ui(2) = SQRT(ui(2))

!       extrapolated velocity
		       ui = LOG(ui) ; pri = LOG(pri)
		       u = ui(1) + (ui(2) - ui(1))/(pri(2) - pri(1))*
	1		    (LOG(pr) - pri(1))
		       u = EXP(u)

c$$$		       print"(50es22.14)",EXP(ui(1)),EXP(ui(2)),u,
c$$$	1		    EXP(pri(1)),EXP(pri(2)),pr,mai(1),mai(2),ma,
c$$$	2		    rai(1),rai(2),ra,lui(1),lui(2),lu

		       DO i=1,nchim
			  d(i,i) = u*l/3.d0	! diffusion coefficient
		       ENDDO

		       EXIT ! since we reached the desired point
		    ENDIF
		 ENDDO
		 
	      ENDIF
	   ELSE ! not a central CZ; intermediate or surface CZ
	      DO i=1,nchim
		 d(i,i)=d_conv
	      ENDDO
	   ENDIF ! central CZ

	ELSE ! we are in a radiative zone; diffusive overshoot 

	   IF(id_conv==0) THEN 	! there is a central CZ; OV
	      grav=cte13*ma/ra**2 ! sem wrot
	      
!       rot=-cte2*w**2*r	!gravite effective avec rotation
!       grav=grav+rot	!remarque de N. Audard
	      
!       echelle de hauteur de pression
	      
	      hp = pr/grav/ro	!echelle de hauteur de pression
	      pb = pr_zc_conv(1) ; fthick = MIN(r_zc_conv(1)*rsol/hp, 1.d0)
	      pc = pr_zc_conv(id_conv)

	      IF(jlim(1)>=5) THEN ! central CZ not too small
		 l = beta_v*hp
!	extrapolação!!!
		 
		 DO j=jlim(1),1,-1
c$$$		    CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(j),lc,
c$$$	1		 fqs,dfqs)
		    CALL inter('m23',bp,q,qt,n_qs,knot,mc(j),fqs,dfqs,r2,m23)
		    
		    pri(1) = EXP(fqs(1))
		    IF(pri(1)>(1.+frac_vp)*pb) THEN ! we passed the critical pr.
		       j1 = j
		       j2 = j + 1 ! the 2 points for extrapolation.
		       IF(pri(2)>=pc) THEN
			  j1 = 1
			  j2 = 2
		       ENDIF
		       
c$$$		       CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(j1),lc,
c$$$	1		    fqs,dfqs)
		       CALL inter('m23',bp,q,qt,n_qs,knot,mc(j1),fqs,dfqs,r2,m23)
		       pri(1) = EXP(fqs(1))
		       tei(1) = EXP(fqs(2))
		       lui(1) = fqs(4)**1.5
		       rai(1) = SQRT(fqs(3))
		       mai(1) = fqs(5)**1.5
		       IF(en_masse)THEN
			  bid = MAX(mc(1),MIN(fqs(5),mc(n_ch)))
		       ELSE
			  bid = MAX(mc(1),MIN(fqs(5)**1.5d0,mc(n_ch)))
		       ENDIF	
		       CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,knotc,.TRUE.,bid,ll,
	1		    xchimi,dxchimi)
		       CALL chim_gram(xchimi,dxchimi)

		       CALL etat(pri(1),tei(1),xchimi,.FALSE.,roi(1),drop1,drot1,
	1		    drox1,uu1,duup1,duut1,duux1,delta1,deltap1,deltat1,
	2		    deltax1,cp1,dcpp1,dcpt1,dcpx1,gradad1,dgradadp1,
	3		    dgradadt1,dgradadx1,alfa1,beta_v1,gamma11)

		       hpup = hp_zc_conv(1) !; PRINT*,'hpup =',hpup
		       zup = (r_zc_conv(1) - rai(1))*rsol
		       zlow = rai(1)*rsol
		       li(1) = (zup + beta_v*hpup)*zlow/(zup + zlow) ! mixing length
!       PRINT*,'j1, j2', j1, j2
		       CALL v2_cgm(pri(1),tei(1),lui(1),mai(1),rai(1),li(1),xchimi,
	1		    ui(1))
		       ui(1) = SQRT(ui(1))

c$$$		       CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(j2),lc,
c$$$	1		    fqs,dfqs)
		       CALL inter('m23',bp,q,qt,n_qs,knot,mc(j2),fqs,dfqs,r2,m23)

		       pri(2) = EXP(fqs(1))
		       tei(2) = EXP(fqs(2))
		       lui(2) = fqs(4)**1.5
		       rai(2) = SQRT(fqs(3))
		       mai(2) = fqs(5)**1.5

		       IF(en_masse)THEN
			  bid = MAX(mc(1),MIN(fqs(5),mc(n_ch)))
		       ELSE
			  bid = MAX(mc(1),MIN(fqs(5)**1.5d0,mc(n_ch)))
		       ENDIF	
		       CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,knotc,.TRUE.,bid,ll,
	1		    xchimi,dxchimi)
		       CALL chim_gram(xchimi,dxchimi)

		       CALL etat(pri(2),tei(2),xchimi,.FALSE.,roi(2),drop1,drot1,
	1		    drox1,uu1,duup1,duut1,duux1,delta1,deltap1,deltat1,
	2		    deltax1,cp1,dcpp1,dcpt1,dcpx1,gradad1,dgradadp1,
	3		    dgradadt1,dgradadx1,alfa1,beta_v1,gamma11)

		       zup = (r_zc_conv(1) - rai(2))*rsol
		       zlow = rai(2)*rsol
		       li(2) = (zup + beta_v*hpup)*zlow/(zup + zlow) ! mixing length
!       PRINT*,'j1, j2', j1, j2
		       CALL v2_cgm(pri(2),tei(2),lui(2),mai(2),rai(2),li(2),xchimi,
	1		    ui(2))
		       ui(2) = SQRT(ui(2))

!       extrapolated velocity at the border of the CZ
		       ui = LOG(ui) ; pri = LOG(pri)
		       ub = ui(1) + (ui(2) - ui(1))/(pri(2) - pri(1))*
	1		    (LOG(pb) - pri(1))
		       ub = EXP(ub)
c$$$		       ub = ui(1) + (ui(2) - ui(1))/(pri(2) - pri(1))*(pb - pri(1))
!       turbulent velocity (exponentially decaying outside the CZ)
		       u = ub*EXP(LOG(pr/pb)/(zeta_v*fthick))
		       
		       DO i=1,nchim
			  d(i,i)=u*l/3.	! diffusion coefficient
		       ENDDO
!       PRINT'(4es12.5)',ra, u, l, u*l/3.
		       
!       PRINT*,'ub, u: ', ub, u, '3'
!       PRINT*,'pb, pc: ', pb, pc
!       PRINT*,'hp, fthick: ', hp, fthick
!       PRINT*,'p1, p2: ', pri(1), pri(2)
!       PRINT*,'u1, u2: ', ui(1), ui(2)
		       EXIT
		    ENDIF
		 ENDDO 
	      ENDIF ! central CZ not too small
	   ENDIF ! there is a central CZ
	   
!       coefficient de diffusivité radiative
	   
	   nu_rad=cte2*te**4/kap/ro**2
	   dnu_radx=-nu_rad*(dkapx/kap+2.d0*drox/ro)
	   DO i=1,nchim
	      d(i,i)=d(i,i)+d_turb+nu_rad ; dd(i,i,1)=dd(i,i,1)+dnu_radx
	   ENDDO
	   u_conv = 0.d0
	   
	ENDIF ! convective or radiative zone (melange)

c$$$	WRITE(999,"(50es22.14)"),nu**1.5,ra,d(1,1),u, u_conv,l,l/hp,pr
	
	RETURN

	END SUBROUTINE difft_ventura
