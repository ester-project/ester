
c********************************************************************

	SUBROUTINE evol(compt,dt,dts,reprend)

c	subroutine public du module mod_evol
	
c	gestion de l'évolution temporelle de la composition chimique
c	les points d'intégration en comp. chim. sont les points de raccord

c	sans diffusion, les ZC sont mélangées, dans ZR on conserve W
c	avec diffusion, W est diffusé

c	on ne tient pas compte de la différence Ptot/Pgaz négligeable
c	dans les régions ou les réactions nucléaires sont actives

c	on reconstruit le vecteur de composition chimique avec
c	discontinuité des xchim aux limites ZR/ZC

c	avec pression turbulente 8 inconnues 
c	sans pression turbulente 7 inconnues, Ptot=Pgaz

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c	23 09 96 : suppression du mouvement des ZR/ZC avec diffusion
c	09 10 96 : modif de la gestion de somme dm/somme r2 dm dans ZC
c	26 06 97 : remplacement du moment angulaire par la v. angulaire
c	25 08 97 : mise en place des variables eulériennes
c	20 10 99 : ajout des variables de structure au temps t dans diffus
c	19 11 99 : suppression de nh1, nhe1, nhe2, lamb
c	24 02 00 : mc_max=(1.d0-1.d-9)*mstar
c	18 04 00 ; coeff_diff ---> diffm, difft, age dans diffus
c	30 07 00 : introduction F95

c entrées
c	compt=0: première intégration t-->dt, .FALSE. : améliorations
c	dt: pas temporel (fixe dans UPDATE)
c	bp,q,qt,knot: solution temps t+dt

c sortie
c	reprend: il faut réduire le pas temporel à cause de TdS ou de non CV
c	estim : estimation de l'erreur max
c	dts : estimation du pas temporel à utiliser pour le pas suivant
c	kmax : indice de couche de l'erreur max

c	modèle totalement convectif: lim=1,jlim(1)=n,lconv(1)=.FALSE.
c	modèle totalement radiatif: lim=0,jlim(i)=-100,lconv(i)=.FALSE.
c	les m_zc sont en m/Mstar

c----------------------------------------------------------------

	USE mod_donnees, ONLY : ab_min, ab_ini, agemax, diffusion, dtmax, dtmin,
	1 en_masse, grille_fixe, ihe4, Ipg, Krot, langue, lvent, mdot,
	2 mvt_dis, m_ch, m_ptm, nchim, ne, nrot, nucleo, pi, pnzc,
	3 pturb, p_vent, ordre, ord_qs, precit
	USE mod_etat, ONLY : etat
	USE mod_kind
	USE mod_nuc, ONLY : vent
	USE mod_numerique, ONLY : bsp1dn, bsp_dis, linf, no_croiss, shell
	USE mod_variables, ONLY : age, bp, bp_t, chim, chim_gram, chim_t,
	1 dim_ch, inter, jlim, jlim_t, knot, knotc, knotc_t, knot_ptm,
	2 knot_t, lconv, lconv_t, lim, lim_t, mc, mct, mct_t, mc_fixe,
	3 mc_t, mstar, m_zc, m_zc_t, m23, m23_t, nc_fixe,
	4 n_ch, n_ch_t, n_ptm, n_qs, n_qs_t, old_ptm, q, qt, qt_t, q_t,
	5 r_zc, r2, r2_t, rota, sortie, tot_conv, tot_rad, xt_ptm, x_ptm

	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: dt
	INTEGER, INTENT(in) :: compt
	REAL (kind=dp), INTENT(out) :: dts
	LOGICAL, INTENT(out) :: reprend

	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: chimd, compx,
	1 chim_zc, chim_zr
	REAL (kind=dp), SAVE, ALLOCATABLE, DIMENSION(:) :: ab_max, bidd,
	1 comp1, comp2, compy, dfdx, esti, estim, f, mc_maxi
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: dm, mc_tmp, ro, ro_t,
	1 t, t_t
     	REAL (kind=dp), DIMENSION(lim) :: m_zc23

c plus de lim. ZR/ZC au delà de Mstar X dmc_max
c plus mvt de lim. ZR/ZC au delà de Mstar X dmc_mvt
c épaisseur des couches pour la comp.chim > Mstar X pas_min		
	REAL (kind=dp), PARAMETER :: dmc_max=1.d-7, dmc_mvt=dmc_max
	REAL (kind=dp), SAVE :: pas_min=dmc_max*1.d-6, mc_max, mc_mvt
	REAL (kind=dp) :: a, alfa, b, beta, bid, cp,  dcpp, dcpt, dcpx,
	1 delta, deltap, deltat, deltax, dgradadp, dgradadt, dgradadx,
	2 drop, drot, drox, dtn, dtnew, dup, dut, dux, est, gamma1, gradad,
	3 mass_zc, mk, mk32, mk32_t, mk_t, mstar23, m_zc23t,
	4 norm, p, pas, p_t, r_ext, u

c nb. de points minimum dans ZR, ZC
	INTEGER, PARAMETER :: n_min_zr=10, n_min_zc=10
	INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: convd, convf, kmaxi
	INTEGER, ALLOCATABLE, DIMENSION(:) :: idis
	INTEGER, SAVE :: ll			
	INTEGER :: i, it, izc, j, jfin, k, kmax, nadd, nc_tmp,
	1 ndis, n_pt_zc, nzc

	LOGICAL, SAVE :: init=.TRUE., demi=.TRUE.
	LOGICAL :: ok, core_reg, tot_conv_t, zc_aug, z_vent

c------------------------------------------------------------------------

2000	FORMAT(8es10.3)
2001	FORMAT(es10.3,3es22.15)

c	PRINT*,'n,nchim,nchim,n_ch,knotc',n_qs,nchim,nchim,n_ch,knotc
c	PRINT*,'dt,t_inf,lim,jlim,l_conv',dt,t_inf,lim,
c	1 (jlim(i),i=1,lim),(lconv(i),i=1,lim)
c	PAUSE'entrée evol'
	
	IF(init)THEN
	 init=.FALSE.
	 
c pour m >  mc_max on ignore les limites ZR/ZC pour le mélange convectif
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1001)dmc_max,dmc_mvt,pas_min
	  WRITE(2,1001)dmc_max,dmc_mvt,pas_min
1001	  FORMAT(/,'--Parameters used in evolution of chem. species--',/,
	2 'no limit RZ/CZ futher than Mstar X ',es10.3,/,
	1 'no change of limit RZ/CZ futher than Mstar X ',es10.3,/,
	2 'thickness of the shells for the chem. comp. > Mstar X ',es10.3) 
	 CASE DEFAULT	 
	  WRITE(*,1)dmc_max,dmc_mvt,pas_min
	  WRITE(2,1)dmc_max,dmc_mvt,pas_min
1	  FORMAT(/,'--Paramètres utilisés pour évoluer la comp. chim.--',/,
	1 'plus de lim. ZR/ZC au delà de Mstar X ',es10.3,/,
	1 'plus mvt de lim. ZR/ZC au delà de Mstar X ',es10.3,/,
	2 'épaisseur des couches pour la comp.chim > Mstar X ',es10.3)
	 END SELECT
	 bid=mstar**(2.d0/3.d0)	; mc_max=(1.d0-dmc_max)*bid
	 mc_mvt=(1.d0-dmc_mvt)*bid
	 pas_min=pas_min*bid	 
	 IF(pturb)THEN
	  WRITE(*,11) ; WRITE(2,11)	  
11	  FORMAT(/,'Evolution temporelle de la composition chimique',/,
	1 'on ne tient pas compte de la différence Ptot-Pgaz',/,
	2 'qui est négligeable dans les régions où les',/,
	3 'réactions thermonucléaires sont actives',/)
	 ENDIF
	 
	 kmax=MAX(nchim,ne)
	 ALLOCATE(ab_max(nchim),bidd(nchim),comp1(nchim),comp2(nchim),
	1 compy(nchim),dfdx(kmax),esti(nchim),estim(nchim),
	2 f(kmax),mc_maxi(nchim),kmaxi(nchim),convd(pnzc+1),
	3 convf(0:pnzc))

c initialisation des paramètres du vent	
	 IF(lvent)THEN
	  ALLOCATE(compx(0,0))
	  CALL vent(comp1,comp2,compx)
	  DEALLOCATE(compx)
	 ENDIF
	ENDIF
	
c initialisations des tableaux de précision
	est=-100.d0 ; estim=-100.d0 ; ab_max=-100.d0
	mc_maxi=-100.d0 ; kmaxi=-100

c sortie si dt < dtmin
	IF(dt < dtmin)THEN
	 SELECT CASE(langue)	
	 CASE('english')
	  WRITE(*,1010)dt,dtmin ; WRITE(2,1010)dt,dtmin	 
1010	  FORMAT('STOP : beginning of evol, dt=',es10.3,' < dtmin=',es10.3)	 	
	 CASE DEFAULT	
	  WRITE(*,10)dt,dtmin ; WRITE(2,10)dt,dtmin	 
10	  FORMAT('ARRET : Entrée de evol, dt=',es10.3,' < dtmin=',es10.3)
	 END SELECT
	 CALL sortie
	ENDIF

c détermination de r_ext
	CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(n_qs),ll,f,dfdx)
	IF(en_masse)THEN
	 r_ext=SQRT(f(3))
	ELSE
	 r_ext=f(3)
	ENDIF	
c 	WRITE(*,2000)r_ext ; PAUSE'r_ext'
	
c masse externe		
	IF(en_masse)THEN
	 mstar23=MIN(mstar**(2.d0/3.d0),m23(n_qs))
	ELSE
	 mstar23=MIN(mstar,m23(n_qs))**(2.d0/3.d0)
	ENDIF

c	modèle totalement convectif: lim=1, jlim(1)=n_qs, lconv(1)=.FALSE.
c	nzc=1, convd(1)=1, convf(1)=nc_tmp
c	modèle totalement radiatif: lim=0, jlim(i)=-100, lconv(i)=.FALSE.
c	nzc=0, convd(1)=nc_tmp, convf(0)=1

	tot_conv=lim == 1 .AND. (jlim(1) == n_qs) .AND. .NOT.lconv(1)
	tot_rad=lim == 0
	
c	PRINT*,tot_conv,lim,jlim(1),.NOT.lconv(1) ; PAUSE'tot_conv'

c	en eulérien et en lagrangien les m_zc et r_zc sont en Rsol et Msol
c	élimination des limites ZR/ZC trop externes

	IF(.NOT.tot_conv .AND. .NOT.tot_rad)THEN
	 i=1
	 DO WHILE(i <= lim)
c	  PRINT*,i,lim,lconv(i),lconv(i-1),m_zc(i),mc_max ; PAUSE'evol 179'
	  IF(m_zc(i) > mc_max)THEN
	   IF(lim == 1)THEN
	    IF(.NOT.lconv(1))THEN 	!fin de l'unique ZC
	     jlim(1)=n_qs	!que l'on déplace à la fin du modèle
	     SELECT CASE(langue)
	     CASE('english')
	      WRITE(*,1008)
1008	      FORMAT(/,'fully mixed model')
	     CASE DEFAULT	     	     	     
	      WRITE(*,8)
8	      FORMAT(/,'Modèle complètement mélangé')
	     END SELECT
	     tot_conv=.TRUE.	!qui devient totalement convectif
	     lim=1 ; jlim(1)=n_qs ; lconv(1)=.FALSE.	     
	    ENDIF
c	   ELSEIF(lconv(i))THEN  !on est sur un début de ZC externe
c	     WRITE(*,7)i
c7	     FORMAT('suppression de la ZC externe',i3)
c	     IF(lim > i)WRITE(*,6)
c	     lim=max(i-1,0)
	   ELSEIF(lconv(i-1))THEN 	!on est sur une fin de ZC
	    jlim(i)=n_qs	!que l'on déplace à la fin du modèle
	    SELECT CASE(langue)
	    CASE('english')	    
	     WRITE(*,1005)i-1	    
1005	     FORMAT(/,'extend to the surface of the CZ #',i3)
	    CASE DEFAULT	    	    
	     WRITE(*,5)i-1	    
5	     FORMAT('prolongement à la surface de la ZC',i3)
	    END SELECT
c	    IF(lim > i)WRITE(*,6)
c6	    FORMAT('et des limites suivantes')	  
	    lim=max(i-1,0) 
	   ENDIF	!lim=1	   
	  ENDIF		!m_zc(i) > mc_max
	  i=i+1
	 ENDDO		!while
	ENDIF
	SELECT CASE(langue)
	CASE('english')	    
	 WRITE(*,1006)
1006	 FORMAT(/,'limits RZ/CZ used for the convective mixing :')	
	CASE DEFAULT
	 WRITE(*,6)
6	 FORMAT(/,'limites ZR/ZC utilisées pour le mélange convectif :')
	END SELECT	
	IF(tot_rad)THEN
	 SELECT CASE(langue)
	 CASE('english')	    
	  WRITE(*,1007)
1007	  FORMAT('fully radiative model')
	 CASE DEFAULT
	  WRITE(*,7)
7	  FORMAT(/,'modèle complètement radiatif')
	 END SELECT	 
	 lim=0 ; jlim=-100 ; lconv=.FALSE.
	ELSEIF(tot_conv)THEN
	 SELECT CASE(langue)
	 CASE('english')	    
	  WRITE(*,1009)
1009	  FORMAT(/'fully convective model')	
	 CASE DEFAULT
	  WRITE(*,9)
9	  FORMAT(/'modèle complètement convectif')
	 END SELECT	 
	ELSE
	 DO i=1,lim
	  IF(lconv(i))THEN
	   SELECT CASE(langue)	  
	   CASE('english')	    
	    WRITE(*,1012)jlim(i),m_zc(i),1.d0-m_zc(i)/mstar,r_zc(i),
	1   r_zc(i)/r_ext
1012	    FORMAT('beginning of CZ shell : ',i4,/,'mass=',es10.3,
	1   ', 1-m/Mstar=',es10.3,', radius=',es10.3,', radius/Rstar=',
	2   es10.3)	  
	   CASE DEFAULT	  
	    WRITE(*,12)jlim(i),m_zc(i),1.d0-m_zc(i)/mstar,r_zc(i),
	1   r_zc(i)/r_ext	  	  
12	    FORMAT('début de ZC couche : ',i4,/,'masse=',es10.3,
	1  ', 1-m/Mstar=',es10.3,', rayon=',es10.3,', rayon/Rstar=',es10.3)
	   END SELECT		   
	  ELSE
	   SELECT CASE(langue)	  
	   CASE('english')	    
	    WRITE(*,1017)jlim(i),m_zc(i),1.d0-m_zc(i)/mstar,r_zc(i),
	1     r_zc(i)/r_ext
1017	    FORMAT('end of CZ shell : ',i4,/,'mass=',es10.3,
	1    ', 1-m/Mstar=',es10.3,', radius=',es10.3,', radius/Rstar=',
	2    es10.3)
	   CASE DEFAULT	   
	    WRITE(*,17)jlim(i),m_zc(i),1.d0-m_zc(i)/mstar,r_zc(i),
	1     r_zc(i)/r_ext
17	    FORMAT('fin de ZC couche : ',i4,/,'masse=',es10.3,
	1    ', 1-m/Mstar=',es10.3,', rayon=',es10.3,', rayon/Rstar=',es10.3)
	   END SELECT
	  ENDIF
	 ENDDO
	ENDIF

c limites ZR/ZC en m^(2/3)	
	m_zc23(1:lim)=m_zc(1:lim)**(2.d0/3.d0)	

c mc_tmp : masses (m^2/3) temporaires pour tabul. de la comp.chim.
c	nc_temp : nombre temporaire de masses

	IF(ALLOCATED(mc_tmp))DEALLOCATE(mc_tmp)
	ALLOCATE(mc_tmp(n_qs+n_qs_t+2*pnzc))	!allocation généreuse

c on prend la répartition déduite de bp ou la grille fixe
c	puis on ajoute les limites ZR/ZC enfin, on s'assure que
c	les ZC et ZR retenues ont, au moins, n_min_zc et n_min_zr points	

	IF(grille_fixe)THEN
	 mc_tmp(1:nc_fixe)=mc_fixe ; nc_tmp=nc_fixe
	ELSE
	 DO i=1,n_qs
	  CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),ll,f,dfdx)
	  IF(en_masse)THEN
	   mc_tmp(i)=f(5)
	  ELSE
	   mc_tmp(i)=f(5)**(2.d0/3.d0)	 
	  ENDIF
	 ENDDO
	 nc_tmp=n_qs
	ENDIF	
c	PRINT*,nc_tmp ; WRITE(*,2000)mc_tmp(1:nc_tmp) ; PAUSE'1'
	
c un pas temporel sur 2 on prend les mc au milieu
	IF(compt == 0)demi=.NOT.demi
	IF(grille_fixe)demi=.FALSE.	!sauf avec la grille fixe
c	PRINT*,demi, nc_tmp ; WRITE(*,*)mc_tmp(nc_tmp-5:nc_tmp)	
	IF(demi)THEN
	 DO i=nc_tmp-1,2,-1
	  mc_tmp(i)=(mc_tmp(i+1)+mc_tmp(i))/2.d0
	 ENDDO
	ENDIF	
c	PRINT*,nc_tmp ; WRITE(*,2000)mc_tmp(1:nc_tmp) ; PAUSE'2'
c	PRINT*,nc_tmp ; WRITE(*,*)mc_tmp(nc_tmp-4:nc_tmp) ; PAUSE'2'

c élimination des couches plus serrées que pas_min
	DO i=2,nc_tmp-1
	 IF(mc_tmp(i)-mc_tmp(i-1) < pas_min)mc_tmp(i)=mc_tmp(i-1)
	ENDDO	
	k=1			!suppression des doubles
	DO i=2,nc_tmp
	 IF(mc_tmp(k) /= mc_tmp(i))THEN
	  k=k+1 ; mc_tmp(k)=mc_tmp(i)
	 ENDIF
	ENDDO
	nc_tmp=k

c couches par ordre croissant	
	CALL shell(nc_tmp,mc_tmp)
c	WRITE(*,*)mc_tmp(nc_tmp-4:nc_tmp) ; PAUSE'21'	
	
c	mc_tmp(nc_tmp) = mstar23	
c on s'assure que mstar23 est à l'extérieur	
c	IF(mc_tmp(nc_tmp) < mstar23)THEN
c	 nc_tmp=nc_tmp+1 ; mc_tmp(nc_tmp) = mstar23
c	ENDIF

c on introduit les limites ZR/ZC
	IF(.NOT.(tot_rad .OR. tot_conv))THEN
	 j=1
	 B4: DO i=2,nc_tmp-1
	  b=mc_tmp(i)-m_zc23(j) ; a=m_zc23(j)-mc_tmp(i-1)
	  IF(a*b < 0.d0)THEN
	   CYCLE B4
	  ELSEIF(b <= a)THEN
	   mc_tmp(i)=m_zc23(j) ; j=j+1	   
	   IF(j > lim)EXIT B4
	  ELSE   	    
	   mc_tmp(i-1)=m_zc23(j) ; j=j+1	   
	   IF(j > lim)EXIT B4	   	  
	  ENDIF
	 ENDDO B4
	ENDIF

c	écritures
	IF(demi)THEN
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(*,1003)
1003	  FORMAT('nodes at middle of the shells')
	 CASE DEFAULT	
	  WRITE(*,3)
3	  FORMAT('couches aux points semi entiers')
	 END SELECT	
	ELSE
	 SELECT CASE(langue)
	 CASE('english')	
	  WRITE(*,1004)
1004	  FORMAT('nodes on limits of shells')
	 CASE DEFAULT	
	  WRITE(*,4)
4	  FORMAT('couches aux points entiers')
	 END SELECT	
	ENDIF

c	délimitation des zones convectives
c	il y a nzc zones convectives chacune entre convd(.) et convf(.)
	convd=-100 ; convf=-100
	IF(tot_conv)THEN
	 nzc=1 ; convd(1)=1 ; convf(1)=nc_tmp
	ELSEIF(tot_rad)THEN
	 nzc=0 ; convd(1)=nc_tmp ; convf(0)=1
	ELSE
	 nzc=0		!nombre de ZC
	 ll=1		!indice des limites
	 DO i=1,nc_tmp
	  IF(m_zc23(ll) == mc_tmp(i))THEN
	   IF(lconv(ll))THEN		!début de ZC
	    nzc=nzc+1			!une ZC en plus
	    convd(nzc)=i		!indice du début de la ZC à droite
	   ELSE				!fin de ZC
	    IF(nzc == 0)THEN		!la zone externe est ZC
	     nzc=1 ; convd(nzc)=1
	    ENDIF
	    convf(nzc)=i		!indice de fin de la ZC à gauche
	   ENDIF
	   ll=MIN(lim,ll+1)		!limite suivante
	  ENDIF
	 ENDDO	!i, pas de limite en n_qs
	 IF(nzc >= 1 .AND. convf(nzc) < 0)convf(nzc)=nc_tmp
	ENDIF
	
c------------------------Pour vérifications-------------------------------

	IF(.FALSE.)THEN
c	IF(.TRUE.)THEN
	 PRINT*,lim,nzc,(convd(i),convf(i),i=1,nzc)
	 DO i=1,nzc
	  PRINT*,i,convd(i),convf(i),mc_tmp(convd(i)),mc_tmp(convf(i))
	 ENDDO	!i
	 PRINT*
	 j=1
	 DO i=1,nc_tmp
	  IF(i == convd(j))THEN	   
	   IF(i /= 1)PRINT*,mc_tmp(i),i
	   PRINT*,'début de ZC',j,i ; PRINT*,mc_tmp(i),i
	  ELSEIF(i == convf(j))THEN
	   PRINT*,mc_tmp(i),i ; PRINT*,'fin de ZC',j,i
	   IF(i /= nc_tmp)PRINT*,mc_tmp(i),i	   
	   j=MIN(j+1,nzc)
	  ELSE
	   PRINT*,mc_tmp(i),i
	  ENDIF
	 ENDDO
	 PAUSE'avant ajustement des ZR et ZC'
	ENDIF
	
c------------------------------------------------------------------------

c les ZR doivent avoir au moins m_ch+1 couches
c	n_min_zr : nombre de couches min. dans ZR
c	n_min_zc : nombre de couches min. dans ZC

c	PRINT*,(mc_tmp(i),i=1,nc_tmp) ; PAUSE'répartition initiale'	   

	IF(nzc > 0)THEN
	
c pour une ZR centrale	
	 IF(convd(1) > 1)THEN
	  nadd=n_min_zr+1-convd(1)
	  IF(nadd > 0)THEN
           WRITE(*,14)nadd
14         FORMAT('adjonction de',i3,' couches pour la ZR centrale')
	   IF(nadd > n_min_zr+1)THEN
	    PRINT*,'erreur nadd > n_min_zr+1' ; PAUSE'ABANDON'
	    CALL sortie
	   ENDIF
	   DO i=nc_tmp,convd(1),-1		!decalages
	    mc_tmp(i+nadd)=mc_tmp(i)
	   ENDDO
	   DO i=1,nzc
	    convd(i)=convd(i)+nadd ; convf(i)=convf(i)+nadd
	   ENDDO
	   nc_tmp=nc_tmp+nadd	!adjonction de masses
	   pas=(mc_tmp(convd(1))-mc_tmp(1))/(n_min_zr+1)
	   DO i=2,n_min_zr+1
	    mc_tmp(i)=mc_tmp(1)+pas*i
	   ENDDO
	  ENDIF
	 ENDIF

c pour les ZR internes
	 DO k=1,nzc-1
	  nadd=n_min_zr+1+convf(k)-convd(k+1)
	  IF(nadd > 0)THEN
         WRITE(*,16)nadd,k
16       FORMAT('adjonction de',i3,' couches pour la ZR qui suit la ZC',i3)
	   IF(nadd > n_min_zr+1)THEN
	    PRINT*,'erreur nadd > n_min_zr+1' ;  PAUSE'ABANDON'
	    CALL sortie
	   ENDIF
	   DO i=nc_tmp,convd(k+1),-1	!decalages
	    mc_tmp(i+nadd)=mc_tmp(i)
	   ENDDO
	   DO i=k+1,nzc
	    convd(i)=convd(i)+nadd ; convf(i)=convf(i)+nadd
	   ENDDO
	   nc_tmp=nc_tmp+nadd	!adjonction de masses
	   pas=(mc_tmp(convd(k+1))-mc_tmp(convf(k)))/(n_min_zr+1)
	   DO i=1,n_min_zr
	    mc_tmp(convf(k)+i)=mc_tmp(convf(k))+pas*i
	   ENDDO
	  ENDIF
	 ENDDO
	 
c pour une ZR externe
	 nadd=0
	 IF(convf(nzc) /= nc_tmp)nadd=n_min_zr+1+convf(nzc)-nc_tmp
	 IF(nadd > 0)THEN
	  WRITE(*,15)nadd
15	  FORMAT('adjonction de',i3,' couches pour une ZR externe')
	  IF(nadd > n_min_zr+1)THEN
	   PRINT*,'erreur nadd > n_min_zr+1' ; PAUSE'ABANDON'
	   CALL sortie
	  ENDIF
c	  PRINT*,'avant',nc_tmp,k,mc_tmp(convf(nzc)),convf(nzc)
	  pas=(mc_tmp(nc_tmp)-mc_tmp(convf(nzc)))/(n_min_zr+1)
	  nc_tmp=nc_tmp+nadd ; mc_tmp(nc_tmp)=mstar23	  
c	  PRINT*,'après',nc_tmp
	  DO i=1,n_min_zr
	   mc_tmp(convf(nzc)+i)=mc_tmp(convf(nzc))+pas*i
c	   PRINT*,convf(nzc)+i,i,mc_tmp(convf(nzc)+i),mc_tmp(convf(nzc)),
c	1  pas*i
	  ENDDO
c	  PRINT*,(mc_tmp(i),i=convf(nzc),nc_tmp)
c	  PAUSE'nouvelle repartition'
	 ENDIF
	
c les ZC doivent avoir aussi plus de n_min_zc+1 couches
	 DO izc=1,nzc
	  nadd=n_min_zc+1+convd(izc)-convf(izc)
c	  PRINT*,'izc,convd(izc),convf(izc),n_min_zc,nadd'	  
c	  PRINT*,izc,convd(izc),convf(izc),n_min_zc,nadd ; PAUSE'izc'	  
	  IF(nadd > 0)THEN
	   WRITE(*,13)nadd,izc         
13         FORMAT('adjonction de',i3,' couches dans la ZC',i3)
	   IF(nadd > n_min_zc+1)THEN
	    PRINT*,'erreur nadd > n_min_zc+1' ; PAUSE'ABANDON'
	    CALL sortie
	   ENDIF
c	   PRINT*,'avant',nc_tmp,izc,mc_tmp(convd(izc)),convd(izc),
c	1	mc_tmp(convf(izc)),convf(izc)
	   DO i=nc_tmp,convf(izc),-1	!decalages
	    mc_tmp(i+nadd)=mc_tmp(i)
	   ENDDO
	   convf(izc)=convf(izc)+nadd
	   DO i=izc+1,nzc
	    convd(i)=convd(i)+nadd ; convf(i)=convf(i)+nadd
	   ENDDO
	   nc_tmp=nc_tmp+nadd	!adjonction de masses	   
c	   PRINT*,'après',nc_tmp,izc,mc_tmp(convd(izc)),convd(izc),
c	1	mc_tmp(convf(izc)),convf(izc)	   
	   pas=(mc_tmp(convf(izc))-mc_tmp(convd(izc)))/(n_min_zc+1)
	   DO i=1,n_min_zc
	    mc_tmp(convd(izc)+i)=mc_tmp(convd(izc))+pas*i
c	    PRINT*,convd(izc)+i,i,izc,mc_tmp(convd(izc)+i),
c	1	mc_tmp(convd(izc)),pas*i
	   ENDDO
c	   PRINT*,convd(izc),convf(izc) ; PRINT*,(mc_tmp(i),i=1,nc_tmp)
c	   PAUSE'nouvelle répartition'	   
	  ENDIF
	 ENDDO
	ENDIF		!sur nzc
	
c	avec diffusion il y aura mélange au delà de mc_max
c	on ajoute / retire des ZM selon les besoins
c	IF(.NOT.(tot_conv .OR. tot_rad) .AND. diffusion)THEN
c	 IF(mc_tmp(convd(nzc)) > mc_max)THEN
c	  B1: DO i=convd(nzc),1,-1
c	   IF(mc_tmp(i) > mc_max)CYCLE B1
c	   IF(mc_tmp(i) <= mc_tmp(convf(nzc-1)))THEN
c	    convf(nzc-1)=convf(nzc) ; nzc=nzc-1
c	   ELSE
c	    convd(nzc)=i
c	   ENDIF
c	   EXIT B1	 
c	  ENDDO B1
c	 ELSEIF(convf(nzc) /= nc_tmp)THEN
c	  B2: DO i=nc_tmp,1,-1
c	   IF(mc_tmp(i) > mc_max)CYCLE B2
c	   nzc=nzc-1 ; convf(nzc)=mc_tmp(nc_tmp) ; convd(nzc)=i
c	   IF(convd(nzc) < convf(nzc-1))THEN
c	    convf(nzc-1)=convf(nzc) ; nzc=nzc-1
c	   ENDIF
c	  ENDDO B2 
c	 ENDIF	
c	ENDIF	

c limites fictives pour faciliter les algorithmes
	IF(convd(1) /= 1)THEN
	 convf(0)=1
	ELSE
	 convf(0)=10000
	ENDIF
	IF(convf(nzc) /= nc_tmp)THEN
	 convd(nzc+1)=nc_tmp
	ELSE
	 convd(nzc+1)=-10000
	ENDIF
	SELECT CASE(langue)
	CASE('english')
	 WRITE(*,1002)nc_tmp
1002	 FORMAT('number of shells for the interpolation of chemicals:',i5)		
	CASE DEFAULT	
	 WRITE(*,2)nc_tmp
2	 FORMAT('nombre de couches utilisées pour la comp. chim.:',i5)
	END SELECT

c masse à la limite inférieure de la zone affectée par le vent
	IF(lvent)mlim_vent=(mstar*(1.d0-p_vent))**(2.d0/3.d0)
	
	mc_tmp(1)=0.d0
	
c---------------------Pour vérifications-----------------------------

c	IF(.TRUE.)THEN
	IF(.FALSE.)THEN
	 PRINT*,'nzc',nzc,nc_tmp
	 PRINT*,'convd',convd(1:nzc+1)
	 PRINT*,'convf',convf(0:nzc)
	 PRINT*,'mlim_vent',mlim_vent
	 PRINT*,'mc_max',mc_max	 
	 DO i=1,nzc
	  PRINT*,i,convd(i),convf(i),mc_tmp(convd(i)),mc_tmp(convf(i))	  
	 ENDDO
	 PRINT*
	 j=1
	 DO i=1,nc_tmp
	  IF(i == convd(j))THEN
	   IF(i /= 1)PRINT*,mc_tmp(i),i,mc_tmp(i)**(3.d0/2.d0)	  
	   PRINT*,'début de ZC',j,i
	   PRINT*,mc_tmp(i),i,mc_tmp(i)**(3.d0/2.d0)
	  ELSEIF(i == convf(j))THEN
	   PRINT*,mc_tmp(i),i,mc_tmp(i)**(3.d0/2.d0)
	   PRINT*,'fin de ZC',j,i
	   IF(i /= nc_tmp)PRINT*,mc_tmp(i),i,mc_tmp(i)**(3.d0/2.d0)	   
	   j=MIN(j+1,nzc)
	  ELSE
	   PRINT*,mc_tmp(i),i,mc_tmp(i)**(3.d0/2.d0)
	  ENDIF
	 ENDDO
	 PAUSE'après ajustement des ZR et des ZC'	 
	ENDIF

c-----------------------  diffusion --------------------------------------

	IF(diffusion)THEN	
c	IF(diffusion .AND. .NOT.tot_conv)THEN	
	
c	 DO i=1,n_ch
c	  CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
c	1 knotc,.TRUE.,mc_tmp(i),ll,comp1,bidd)
c	  WRITE(*,2002)mc_tmp(i),(comp1(j),j=1,nchim)
2002	  FORMAT(13es8.1)
c	 ENDDO

	 CALL diffus(ok,dt,convd,convf,nzc,mc_tmp,nc_tmp)
	 reprend=.NOT.ok	 
	 IF(reprend)THEN
	 SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1018) ; WRITE(2,1018)
1018	   FORMAT('Reinitialize and dt-->dt/2')	 	
	  CASE DEFAULT	
	   WRITE(*,18) ; WRITE(2,18)
18	   FORMAT('Réinitialisation et dt-->dt/2')  	  
	  END SELECT	  
	  IF(dt > dtmin)THEN
	   RETURN
	  ELSE
	   SELECT CASE(langue)
	   CASE('english')	  
	    WRITE(*,1019)dt,dtmin ; WRITE(2,1019)dt,dtmin
1019	    FORMAT('STOP, dt =,'es10.3,' <',es10.3,' = dtmin')  	  	  	
	   CASE DEFAULT	
	    WRITE(*,19)dt,dtmin ; WRITE(2,19)dt,dtmin
19	    FORMAT('ARRET, dt =,'es10.3,' <',es10.3,' = dtmin')  	  
	   END SELECT
	   CALL sortie
	  ENDIF	
	 ENDIF	 

c estimation de la précision, à cause de sauts ZR --> ZC avec variation
c	éventuellement forte de comp.chim. on élimine les ZC du test

c	 WRITE(*,2002)(mct_t(i),i=1,knotc_t)
c	 PAUSE'estimation de la précision'

	 IF(tot_conv)THEN
	  it=5 ; jfin=it
	 ELSE
	  it=1 ; jfin=n_ch
	 ENDIF
	 
c	 DO i=1,n_ch
	 DO i=it,jfin
	  CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,knotc,.TRUE.,mc(i),
	1   ll,comp1,bidd)
	  IF(en_masse)THEN
	   bid=MIN(mc(i),x_ptm(n_ptm))
	  ELSE
	   bid=MIN(SQRT(ABS(mc(i)))**3,x_ptm(n_ptm))
	  ENDIF	  
	  CALL bsp1dn(1,old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm,.TRUE.,
	1   bid,ll,f,dfdx)	!masse au temps t mk_t
	  IF(.NOT.en_masse)f(1)=ABS(f(1))**(2.d0/3.d0)
	  mk_t=MIN(MAX(mc_t(1),f(1)),mc_t(n_ch_t))
	  CALL bsp1dn(nchim,chim_t,mc_t,mct_t,n_ch_t,m_ch,	   
	1   knotc_t,.TRUE.,mk_t,ll,comp2,bidd)
	
c estimation de la précision
	  DO j=1,nchim
	   bid=MAX(ABS(comp1(j)),ABS(comp2(j)))
	   IF(bid > 0.d0)THEN
	    esti(j)=ABS(comp1(j)-comp2(j))/bid
	   ELSE
	    esti(j)=0.d0
	   ENDIF	   
	   IF(ABS(comp1(j)) > ab_min(j))THEN
	    IF(j == 1 .OR. j == ihe4)THEN
	     est=MAX(est,esti(j))
	     IF(est == esti(j))kmax=i
	    ENDIF
	    IF(estim(j) <= esti(j))THEN		!précision par élément	    
c	     PRINT*,i,j ; WRITE(*,2000)estim(j),esti(j),est	     
	     estim(j)=esti(j) ; kmaxi(j)=i
	     ab_max(j)=comp1(j) ; mc_maxi(j)=mc(i)
	    ENDIF
	   ENDIF
	  ENDDO
	 ENDDO
	  
c	 PAUSE'après diffusion'

c----------------------------- sans diffusion ------------------------------

	ELSE	!intégration et mélange classiques

c	 PRINT*,'knotc',knotc,nchim ; WRITE(*,2000)dt

c	 PRINT*,'avant intégr.: n_ch_t,knotc_t,compt',n_ch_t,knotc_t,compt
c	 PRINT*,'mc_t',mc_t(n_ch_t) ; WRITE(*,2000)(mc_t(i),i=1,nc_tmp)
c	 PRINT*,'mct_t',mct_t(knotc_t)
c	 WRITE(*,2000)(mct_t(i),i=1,knotc_t) ; PRINT*,'xchim'
c	 DO i=1,nc_t
c	  CALL bsp1dn(nchim,chim_t,mc_t,mct_t,n_ch_t,m_ch,
c	1	knotc_t,.TRUE.,mc_t(i),ll,compx,compy)
c	  WRITE(*,2000)mc_t(i),(compx(j)*nucleo(j),j=1,MIN(7,nchim))
c	  WRITE(*,2002)mc_t(i),(compx(j)*nucleo(j),j=1,nchim)
c	 ENDDO
c	 PAUSE'avant intégration'

c on intègre par rapport à t en mélangeant les ZM
c	 d'abord les ZC puis les ZR
c	 on garde dans chim_zc(nchim,nzc) les valeurs des ZM
c	 la variable en masse est mc=m**2/3
c	 convf(0)=1 si convd(1) /= 1, sinon convf(0)=10000
c	 convd(nzc+1)=nc_tmp si convf(nzc) /= nc_tmp,
c	 sinon convd(nzc+1)=-10000
c	 dans ZM on mélange MW, dans ZR on conserve MW
c	 on reconstruit le vecteur de composition chimique avec
c	 discontinuité des xchim aux LMR.
	 
c les ZC d'abord

	 ALLOCATE(chim_zc(nchim,nzc))
	 DO izc=1,nzc
	  mass_zc=0.d0
	  
c n_pt_zc : nombre de points dans la ZC	  
	  n_pt_zc=convf(izc)-convd(izc)
	  
	  ALLOCATE(compx(nchim,n_pt_zc),ro(n_pt_zc),ro_t(n_pt_zc),
	1 t(n_pt_zc),t_t(n_pt_zc),dm(n_pt_zc))
	  DO i=convd(izc),convf(izc)-1
	   k=i-convd(izc)+1
	   mk=(mc_tmp(i)+mc_tmp(i+1))/2.d0		!point milieu
	   mk32=SQRT(ABS(mk))**3 ; dm(k)=(mc_tmp(i+1)-mc_tmp(i))*SQRT(mk)
	   mass_zc=mass_zc+dm(k)		!masse de la ZC
c	   WRITE(*,2000)mass_zc,dm(k),mc_tmp(i),mc_tmp(i+1),SQRT(mk),mk

c composition chimique en mk_t au temps t, et en mk au temps t+dt
c	   on tient compte de la perte/gain de masse

	   IF(en_masse)THEN
	    bid=MIN(mk,x_ptm(n_ptm))
	   ELSE
	    bid=MIN(mk32,x_ptm(n_ptm))
	   ENDIF
	   CALL bsp1dn(1,old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm,.TRUE.,
	1  bid,ll,f,dfdx)	!masse au temps t mk_t
	   IF(en_masse)THEN
	    bid=f(1)
	   ELSE
	    bid=ABS(f(1))**(2.d0/3.d0)
	   ENDIF	    
	   mk_t=MAX(mc_t(1),bid)
	   mk32_t=SQRT(ABS(mk_t))**3
	   IF(mk_t>mc_t(n_ch_t)) THEN
	    comp1=ab_ini/nucleo
	   ELSE
	    CALL bsp1dn(nchim,chim_t,mc_t,mct_t,n_ch_t,m_ch,
	1   knotc_t,.TRUE.,mk_t,ll,comp1,bidd)
	   ENDIF
	   
c	   PRINT*,'après bsp1dn 1, izc,i,k',izc,i,k
c	   WRITE(*,2000)mk,(comp1(j),j=1,MIN(7,nchim))
c	   PAUSE'après bsp1dn 1'

c	   comp1, compx : abondances au temps t
c	   comp2, compy : abondances au temps t+dt
c	   compx, compy : abondances par mole, comp1, comp2 par gramme

	   comp1=MAX(comp1,1.d-100) ; compx(:,k)=comp1(:)

	   IF(compt == 0)THEN		!comp. chim. en t+dt	    
	    comp2=comp1
	   ELSE
	    CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1   knotc,.TRUE.,MIN(mk,mct(knotc)),ll,comp2,bidd)
c	    WRITE(*,2000)mk,(compy(j),j=1,MIN(7,nchim))
c	    PAUSE'après bsp1dn 2'
	    comp2=MAX(comp2,1.d-100) 
	   ENDIF		!new
	   
c	   comp1 et comp2 par gramme pour appel à l'EOS
	   
	   CALL chim_gram(comp1,bidd) ; CALL chim_gram(comp2,bidd)	   
	   
c	   T et ro en t+dt par interpolation inverse en m23 par inter

c	   PRINT*,'avant, izc,i,k',izc,i,k

	   IF(en_masse)THEN
	    CALL inter('m23',bp,q,qt,n_qs,knot,MIN(mk,m23(n_qs)),
	1   f,dfdx,r2,m23)	    
	   ELSE
	    CALL inter('m23',bp,q,qt,n_qs,knot,MIN(mk32,m23(n_qs)),
	1   f,dfdx,r2,m23)	    
	   ENDIF
	   IF(pturb)THEN	!avec pression turbulente
	    p=EXP(f(Ipg))	 
	   ELSE			!sans pression turbulente
	    p=EXP(f(1))
	   ENDIF
	   t(k)=EXP(f(2))
	   CALL etat(p,t(k),comp2,.FALSE.,			!les ro
	1  ro(k),drop,drot,drox,u,dup,dut,dux,
	2  delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3  gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	
c	   WRITE(*,2000)p,t(k),ro(k),mk
	   
c T et ro en t par interpolation inverse en m23	    
	   IF(en_masse)THEN
	   CALL inter('m23',bp_t,q_t,qt_t,n_qs_t,knot_t,
	1  MIN(mk_t,m23_t(n_qs_t)),f,dfdx,r2_t,m23_t)	    
	   ELSE
	    bid=mk32_t
	   CALL inter('m23',bp_t,q_t,qt_t,n_qs_t,knot_t,
	1  MIN(mk32_t,m23_t(n_qs_t)),f,dfdx,r2_t,m23_t)  
	   ENDIF		    
	   IF(pturb)THEN	!avec pression turbulente
	    p_t=EXP(f(Ipg))	 
	   ELSE			!sans pression turbulente
	    p_t=EXP(f(1))
	   ENDIF
	   t_t(k)=EXP(f(2))
	   CALL etat(p_t,t_t(k),comp1,.FALSE.,
	1  ro_t(k),drop,drot,drox,u,dup,dut,dux,
	2  delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3  gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	
c	   WRITE(*,2000)p_t,t_t(k),ro_t(k),mk
	  ENDDO		!i dans chaque ZC

c	  PRINT*,k,n_pt_zc ; PAUSE'evol: n_pt_zc'
c	  WRITE(*,2000)compx(1,1:k)  	   
c	  WRITE(*,2000)dm ; WRITE(*,2000)t
c	  WRITE(*,2000)ro ; WRITE(*,2000)mass_zc ; PAUSE'après les dm(k)'

c normalisation des dm
	  dm=dm/mass_zc
	  
c est-on dans la zone affectée par le vent ?	  
	  z_vent=mc_tmp(i) >= mlim_vent
	  
c intégration temporelle  
	  CALL rk_imps(t_t,ro_t,compx,t,ro,compy,dt,esti,ok,n_pt_zc,
	1 z_vent,dm)
c	  IF(z_vent)THEN
c	   WRITE(*,2000)compx(1:nchim,1) ; WRITE(*,2000)compy(1:nchim)
c	   WRITE(*,2000)SUM(compy*nucleo)	
c	   PAUSE'solution z_vent'
c	  ENDIF
	
c	  PRINT*,ok ; WRITE(*,2000)compy ; PAUSE'après rk_imps'
	
c 	  on garde, provisoirement, les abondances dans chim_zc
	   
	  IF(ok)THEN		!il y a convergence
	   chim_zc(:,izc)=compy(:) !chim_zc: comp. chim. de la ZC
	   
c 	   gestion de la précision locale
	   
	   DO j=1,nchim
	    IF(compy(j) > ab_min(j) .AND. estim(j) <= esti(j))THEN
	     estim(j)=esti(j)	!précision par élément
	     kmaxi(j)=convd(izc) ; ab_max(j)=compy(j)
	     mc_maxi(j)=mc_tmp(convd(izc))
	    ENDIF
c	    PRINT*,'ZC',convd(izc),j,estim(j),esti(j),est
	   ENDDO		!j=1,nchim
c	   WRITE(*,2000)estim ; WRITE(*,2000)esti ; WRITE(*,2000)mc_maxi
c	   WRITE(*,2000)ab_max ; WRITE(*,2000)compy
c	   PRINT*,kmaxi ; PAUSE'ZC'     

c	   précision globale estimee sur He4 ou H1	   
	   
	   IF(ihe4 > 1)THEN
	    est=MAX(est,esti(ihe4)) ; IF(est == esti(ihe4))kmax=convd(izc)
	   ELSE 
	    est=MAX(est,esti(1)) ; IF(est == esti(1))kmax=convd(izc)
	   ENDIF
	  ELSE
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,2020)izc ; WRITE(2,2020)izc
2020	    FORMAT('EVOL/rk_imps Pb. at ZC #',i3,'Reboot and dt-->dt/2') 
	   CASE DEFAULT
	    WRITE(*,20)izc ; WRITE(2,20)izc
20	    FORMAT('EVOL/rk_imps Pb. à la ZC #',i3,
	1     'Réinitialisation et dt-->dt/2') 
	   END SELECT
	   IF(dt > dtmin)THEN
	    reprend=.TRUE. ; RETURN
	   ELSE
	    SELECT CASE(langue)
	    CASE('english')	  
	     WRITE(*,1019)dt,dtmin ; WRITE(2,1019)dt,dtmin  	  	
	    CASE DEFAULT	
	     WRITE(*,19)dt,dtmin ; WRITE(2,19)dt,dtmin	  
	    END SELECT
	    CALL sortie
	   ENDIF	
	  ENDIF

c	  PRINT*,'zone mélangee',izc
c	  WRITE(*,2000)mc_tmp(convd(izc)),mc_tmp(convf(izc))
c	  WRITE(*,2002)(chim_zc(1:nchim,izc) ; PAUSE'ZM'

	  DEALLOCATE(compx,ro,ro_t,t,t_t,dm)	  
	 ENDDO	!izc
	 
c ensuite les ZR
	 ALLOCATE(compx(nchim,1),ro(1),ro_t(1),t(1),t_t(1),dm(1),
	1 chim_zr(nchim,nc_tmp))
	 dm=1.d0
	 DO izc=1,nzc+1
c	  PRINT*,'ZR',izc,convf(izc-1),convd(izc)	  	  
	  DO i=convf(izc-1),convd(izc)

c	   la comp.chim. en mk_t au temps t et en mc_tmp(i) au temps t+dt
c	   les mc peuvent différer des mct

	   IF(en_masse)THEN
	    bid=mc_tmp(i)
	   ELSE
	    bid=SQRT(ABS(mc_tmp(i)))**3
	   ENDIF
	   CALL inter('m23',bp,q,qt,n_qs,knot,MIN(bid,m23(n_qs)),f,
	1  dfdx,r2,m23)
	   bid=MIN(bid,x_ptm(n_ptm)) 
	   CALL bsp1dn(1,old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm,.TRUE.,
	1  bid,ll,f,dfdx)	!masse au temps t
	   IF(.NOT.en_masse)f(1)=ABS(f(1))**(2.d0/3.d0)
	   mk_t=MAX(mc_t(1),f(1)) ; mk32_t=SQRT(mk_t)**3
	   IF(mk_t>mc_t(n_ch_t)) THEN
	    comp1=ab_ini/nucleo
	   ELSE
	    CALL bsp1dn(nchim,chim_t,mc_t,mct_t,n_ch_t,m_ch,
	1   knotc_t,.TRUE.,mk_t,ll,comp1,bidd)
	   ENDIF
	
c	   PRINT*,'après bsp1dn 1, izc,i',izc,i
c	   WRITE(*,2000)mc_tmp(i),(comp1(j),j=1,MIN(7,nchim))
c	   PAUSE'après bsp1dn 1 radiatif'

	   comp1=MAX(comp1,1.d-100) ; compx(:,1)=comp1	    

	   IF(compt == 0)THEN		!comp. chim. en t+dt
	    comp2=comp1
	   ELSE
	    CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1   knotc,.TRUE.,MIN(mc_tmp(i),mct(knotc)),ll,comp2,bidd)
c	    WRITE(*,2000)compy ; PAUSE'après bsp1dn 2'
	    comp2=MAX(comp2,1.d-100)
	   ENDIF	!compt=0
	   CALL chim_gram(comp1,bidd) ; CALL chim_gram(comp2,bidd)
	   
c	   les T et ro en t+dt, par interpolation inverse en m23

c	   PRINT*,'avant, izc,i',izc,i

	   IF(en_masse)THEN
	    bid=mc_tmp(i)	    
	   ELSE
	    bid=SQRT(ABS(mc_tmp(i)))**3	    
	   ENDIF
	   CALL inter('m23',bp,q,qt,n_qs,knot,MIN(bid,m23(n_qs)),f,
	1  dfdx,r2,m23)
	   IF(pturb)THEN	!avec pression turbulente
	    p=EXP(f(Ipg))	 
	   ELSE			!sans pression turbulente
	    p=EXP(f(1))	 
	   ENDIF
	   t(1)=EXP(f(2))
	   CALL etat(p,t(1),comp2,.FALSE.,			!les ro
	1  ro(1),drop,drot,drox,u,dup,dut,dux,
	2  delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3  gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
c	   WRITE(*,2000)p,t(1),ro(1),mc_tmp(i)
	   
c	   les T et ro en t, par interpolation inverse en m23

	   CALL bsp1dn(1,old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm,.TRUE.,
	1  MIN(bid,x_ptm(n_ptm)),ll,f,dfdx)
	   mk_t=f(1)	   
	   CALL inter('m23',bp_t,q_t,qt_t,n_qs_t,knot_t,
	1  MIN(mk_t,m23_t(n_qs_t)),f,dfdx,r2_t,m23_t)
	   IF(pturb)THEN	!avec pression turbulente
	    p_t=EXP(f(Ipg))	 
	   ELSE			!sans pression turbulente
	    p_t=EXP(f(1))	 
	   ENDIF
	   t_t(1)=EXP(f(2))
	   CALL etat(p_t,t_t(1),comp1,.FALSE.,
	1  ro_t(1),drop,drot,drox,u,dup,dut,dux,
	2  delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3  gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
c	   WRITE(*,2000)p_t,t_t(1),ro_t(1),mc_tmp(i)

c est-on dans la zone affectée par le vent ?	  
	   z_vent=mc_tmp(i) >= mlim_vent
	   
c intégration temporelle	   
	   CALL rk_imps(t_t,ro_t,compx,t,ro,compy,dt,esti,ok,1,z_vent,dm)
	   compy=ABS(compy)
c	   WRITE(*,2000)mc_tmp(i),t,ro,compy(1),compx(1,1)
c	   PRINT*,'dans ZR',i
	   
c	   on place les abondances dans chim_zr
	   
	   IF(ok)THEN		!il y a eu convergence
	    chim_zr(:,i)=compy
	   
c	    gestion de la précision locale	   

	    DO j=1,nchim
	     IF(compy(j) > ab_min(j) .AND. estim(j) <= esti(j))THEN
	      estim(j)=esti(j)	!précision par élément
	      kmaxi(j)=i ; ab_max(j)=compy(j)
	      mc_maxi(j)=mc_tmp(i)
	     ENDIF
c	     PRINT*,'ZR, couche',i,j,estim(j),esti(j),est
	    ENDDO		!j=1,nchim	     
c	    WRITE(*,2000)estim ; WRITE(*,2000)esti ; WRITE(*,2000)mc_maxi
c	    WRITE(*,2000)ab_max ; PRINT*,kmaxi ; WRITE(*,2000)compy
c	    PAUSE'ZR'     

c	    précision globale estimée sur He4 ou H1	   
	   
	    IF(ihe4 > 1)THEN
	     est=MAX(est,esti(ihe4)) ; IF(est == esti(ihe4))kmax=convd(izc)
	    ELSE 
	     est=MAX(est,esti(1)) ; IF(est == esti(1))kmax=convd(izc)
	    ENDIF
	   ELSE			!ok
	   
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,2021)i ; WRITE(2,2021)i
2021	     FORMAT('EVOL/rk_imps Pb. at the shell #',i3,
	1     'Reboot and dt-->dt/2') 
	    CASE DEFAULT
	     WRITE(*,21)i ; WRITE(2,21)i
21	     FORMAT('EVOL/rk_imps Pb. à la couche #',i3,
	1     'Réinitialisation et dt-->dt/2') 
	    END SELECT
	    IF(dt > dtmin)THEN
	     reprend=.TRUE. ; RETURN
	    ELSE
	     SELECT CASE(langue)
	     CASE('english')	  
	      WRITE(*,1019)dt,dtmin ; WRITE(2,1019)dt,dtmin  	  	
	     CASE DEFAULT	
	      WRITE(*,19)dt,dtmin ; WRITE(2,19)dt,dtmin	  
	     END SELECT
	     CALL sortie
	    ENDIF	
	   ENDIF
	  ENDDO	!i
	 ENDDO		!izc

c	 PAUSE'après intégration'
	 
c	 cas sans diffusion, mouvements des limites ZR/ZC
c	 mvt_dis=.TRUE. on modifie la comp. chim. dans la ZR	 

	 IF(.NOT.(tot_conv .OR. tot_rad) .AND. mvt_dis)THEN
	  tot_conv_t=lim_t == 1 .AND. jlim_t(1) == n_qs_t 
	1 .AND. .NOT.lconv_t(1)
	 
c	  cas de la régression d'un coeur convectif :
c	  on interpole linéairement la comp.chim. sur la zone de régression

	  core_reg=.NOT.lconv(1) .AND. .NOT.lconv_t(1)
	1 .AND. m_zc(1) < m_zc_t(1) .AND. m_zc(1) < mstar
	
c	  PRINT*,core_reg,tot_conv,tot_rad,lconv(1),lconv_t(1),
c	1 m_zc(1),m_zc_t(1),convd,convf ; PAUSE
	
c	  A droite de lim. on met la comp. chim. de la ZC

	  IF(core_reg)THEN  
	   chim_zr(:,convf(1))=chim_zc(:,1) 
	   m_zc23t=m_zc_t(1)**(2.d0/3.d0)
	   CALL linf(m_zc23t,mc_tmp,nc_tmp,ll)
c	   PRINT*,convf(1),ll,mc_tmp(convf(1)),m_zc23t
	   
c	   la ZC a regresse de m_zc23t à mc_tmp(convf(1)) et on a :
c	   mc_tmp(ll+1) > m_zc23t >= mc_tmp(ll) >= mc_tmp(convf(1))
c	   interp. lin. de la comp. chim. pour mc entre
c	   mc_tmp(ll+1) et mc_tmp(convf(1))
	  
	   DO i=convf(1)+1,ll
c	    PRINT*,i,mc_tmp(i)
	    bid=(mc_tmp(i)-mc_tmp(convf(1)))/(mc_tmp(ll+1)
	1   -mc_tmp(convf(1)))
	    chim_zr(:,i)=chim_zr(:,convf(1))+(chim_zr(:,ll)
	1   -chim_zr(:,convf(1)))*bid
	   ENDDO	 
	  ENDIF
	  
c	  cas de l'augmentation de l'abscisse lagrangienne (masse) de la
c	  base d'une ZC suffisamment interne
c	  on interp. linéairement la comp.chim. sur la zone d'augmentation

	  DO izc=1,nzc
	   zc_aug=.FALSE.
	   IF(mc_tmp(convd(izc)) <= mc_mvt)THEN
	    IF(tot_conv_t .AND. convd(izc) > 1)THEN
	     zc_aug=.TRUE. ; m_zc23t=0.d0
	    ELSEIF(lim_t == 1)THEN	!au temps t : une seule ZC
	     m_zc23t=m_zc_t(1)**(2.d0/3.d0)
	     zc_aug=m_zc23t < mc_tmp(convd(izc)) .AND. lconv_t(1)
	    ELSE		!plusieurs ZC au temps t
	     jfin=lim_t-1 ; it=1
	     DO WHILE(it <= jfin)
	      zc_aug=m_zc_t(it) < mc_tmp(convd(izc)) .AND. lconv_t(it)
	1     .AND. m_zc_t(it+1)**(2.d0/3.d0) > mc_tmp(convd(izc))
	      IF(zc_aug)THEN
	       jfin=-4 ; m_zc23t= m_zc_t(it)	!il y a encadrement
	      ELSE
	       it=it+1
	      ENDIF
	     ENDDO
	    ENDIF
   
c	    PRINT*,zc_aug,izc,convd(izc),it,tot_conv_t,lim_t,
c	1   jlim_t(1),n_qs_t,lconv_t(1)
c	    WRITE(*,2000)m_zc23t,mc_tmp(convd(izc)),
c	    m_zc_t(it+1)**(2.d0/3.d0)

	    IF(zc_aug)THEN
	     CALL linf(m_zc23t,mc_tmp,nc_tmp,ll)
c	     PRINT*,ll
c	     la lim. de la ZC est passée de m_zc23t=mc_tmp(ll) a
c	     mc_tmp(convd(iz))
c	     interp. lin. de la comp. chim. pour mc entre mc_tmp(ll)
c	     et mc_tmp(convd(izc))

	     DO i=ll+1,convd(izc)
c	      PRINT*,i,mc_tmp(i),mc_tmp(ll)
	      bid=(mc_tmp(i)-mc_tmp(convd(izc)))/(mc_tmp(ll)
	1     -mc_tmp(convd(izc)))
	      chim_zr(:,i)=chim_zc(:,izc)+(chim_zr(:,ll)
	1     -chim_zc(:,izc))*bid
	     ENDDO
c	     PRINT*,convd(izc)
	    ENDIF
	   ENDIF
	  ENDDO
	 ENDIF

c	 PAUSE'après mouvements des ZR/ZC'

c on modifie le vecteur de composition chimique

c	 convf(0)=1 si convd(1) /= 1, sinon convf(0)=10000
c	 convd(nzc+1)=nc_tmp si convf(nzc) /= nc_tmp,
c	 sinon convd(nzc+1)=-10000

c	 PRINT*,nzc;  PRINT*,(convf(i),i=0,nzc)
c	 PRINT*,(convd(i),i=1,nzc+1) ; PAUSE'nzc'

c pour allocations, estimation du nombre de discontinuités ndis

	 ndis=0		!nombre de discontinuités
	 DO izc=1,nzc+1
	  DO i=convf(izc-1),convd(izc)	!zone radiative  ZR
	   IF(i == convf(izc-1) .AND. i > 1)ndis=ndis+1
	  ENDDO
	  IF(izc <= nzc)THEN
	   DO i=convd(izc),convf(izc)
	    IF(i == convd(izc) .AND. i > 1)ndis=ndis+1
	   ENDDO
	  ENDIF
	 ENDDO

c redéfinition du vecteur de composition chimique nc_tmp-->n_ch
c	 chimd : abondances aux discontinuités, idis(ndis): indices discon.
c	 la seconde dim. du vecteur de comp.chim. chim(nchim,n_ch+ndis)
c	 tient compte des discontinuités, la dimension du vecteur nodal
c	 correspondant sera knotc=n_ch+m_ch+ndis, dim_ch=n_ch+ndis

	 n_ch=nc_tmp ; DEALLOCATE(chim) 
	 ALLOCATE(chim(nchim,n_ch+ndis),chimd(nchim,ndis),idis(0:ndis+1))
	 idis=-100
	 	
c	 j est l'indice de la discontinuité, à la fin de la séquence j=ndis

	 j =0 
	 DO izc=1,nzc+1
	  DO i=convf(izc-1),convd(izc)	!zone radiative  ZR
	   IF(i == convf(izc-1) .AND. i > 1)THEN
	    j=j+1 ; chimd(:,j)=chim_zr(:,i) ; idis(j)=i
c	    PRINT*,'j,i',j,i ; PRINT*,mc(i)
c	    PAUSE'début de ZR à droite de ZM'
	   ELSE
	    chim(:,i)=chim_zr(:,i)
	   ENDIF
	  ENDDO
	  
	  IF(izc <= nzc)THEN		!pour ZM	   
	   DO i=convd(izc),convf(izc)
	    IF(i == convd(izc) .AND. i > 1)THEN
	     j=j+1 ; chimd(:,j)=chim_zc(:,izc) ; idis(j)=i
c	     PRINT*,'j,i,izc',j,i,izc
c	     WRITE(*,2000)chim_zc(1,izc),chimd(1,j)
c	     PAUSE'début de ZM à droite de ZR'
	    ELSE    
	     chim(:,i)=chim_zc(:,izc)
	    ENDIF	!i == convd
	   ENDDO	!i=convd , convf
	  ENDIF		!izc <= nzc
	 ENDDO		!nzc
	 
c	 vérification: on doit avoir j=ndis

	 IF(j /= ndis)THEN
	  PRINT*,'evol: j=',j,' /= ndis=',ndis ; PRINT*,'ARRET'
	  CALL sortie
	 ENDIF
	 
c	 PRINT*,nc,nzc,ndis,idis
c	 ll=1
c	 DO i=1,n_ch
c	  WRITE(*,2000)mc(i),chim(1:MIN(4,nchim),i)	 
c	  IF(i == idis(ll))THEN
c	   PRINT*,idis(ll) ; WRITE(*,2000)mc,chimd(1:MIN(4,nchim),ll)
c	   ll=MIN(ll+1,ndis)
c	  ENDIF
c	 ENDDO
c	 DO i=1,ndis
c	  WRITE(*,2000)chimd(1:MIN(4,nchim),i)
c	 ENDDO
c	 PAUSE'après ecritures'

c	 transfert des masses temporaires
c	 la dimension du vecteur dual mct est knotc=n_ch+m_ch+ndis

	 DEALLOCATE(mc,mct) ; ALLOCATE(mc(n_ch),mct(n_ch+m_ch+ndis))
	 mc(1:n_ch)=mc_tmp(1:n_ch)

c	 tabulation de la composition chimique, bsp_dis envisage les cas
c	 avec et sans discontinuité

c	 PRINT*,nchim,ndis,idis,n_ch,m_ch,mvt_dis
c	 WRITE(*,2000)mc ; WRITE(*,2000)chim(1,1:n_ch)
c	 WRITE(*,2000)chimd(1,1:ndis) ; PAUSE'avant bsp_dis'

	 CALL bsp_dis(nchim,mc,chim,ndis,idis,chimd,eps,n_ch,m_ch,mct,knotc)
         IF(no_croiss)THEN
          PRINT*,'Arrêt 1 dans evol' ; CALL sortie
         ENDIF
	 dim_ch=knotc-m_ch
	 
c	 PRINT*,'après intégration : nc,knotc',nc,knotc
c	 PRINT*,'mc',mc(nc) ; WRITE(*,2000)mc
c	 PRINT*,'mct',mct(knotc) ; WRITE(*,2000)mct ; PRINT*,'xchim'
c	 DO i=1,nc
c	  CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
c	1	knotc,.TRUE.,mc(i),ll,compx,compy)
c	  WRITE(*,2000)mc(i),(compx(j)*nucleo(j),j=1,MIN(7,nchim))
c	  WRITE(*,2002)mc(i),(compx(j)*nucleo(j),j=1,nchim)
c	  WRITE(*,2002)mc(i),(compx(j),j=1,nchim)	 
c	 ENDDO
c	 PAUSE'après intégration'
c
c	 PAUSE'test rota'
c	 DO i=1,nc
c	  CALL inter('m23',bp,q,qt,n_qs,knot,MIN(m23(i),m23(n_qs)),
c	1	f,dfdx,r2,m23)
c	  CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
c	1	knotc,.TRUE.,m23(i),ll,compx,compy)
c	  IF(r2(i) > 0.d0)PRINT*,m23(i),r2(i),,compx(1)	
c	  IF(f(3) > 0.d0)THEN
c	   PRINT*,m23(i),f(3)
c	  ELSE
c	   PRINT*,m23(i),f(3)
c	  ENDIF
c	 ENDDO
c	 PAUSE'après test'

	 DEALLOCATE(chimd,compx,chim_zc,chim_zr,dm,ro,ro_t,t,t_t,idis)
	ENDIF	!diffusion ou intégration et mélange classiques
	DEALLOCATE(mc_tmp)
		
c	normalisation: afin d'obtenir Somme Xi=1	
c	on compense le fait que les nucleo i.e. les masses des éléments
c	chimiques ne sont pas des nombres entiers	
c	pour chaque élément, dans la base des splines, on forme
c	norn = somme des coefficients X nucleo, norm est proche de 1
c	on divise ensuite les coefficients par norm
c	ce qui est valable parceque la somme des B-splines normalisées = 1

	IF(nchim > 1)THEN
	 DO i=1,dim_ch	!dim_ch=knotc-m_ch: dimension de la base
	  norm=SUM(chim(1:nchim,i)*nucleo(1:nchim))
	  chim(1:nchim,i)=chim(1:nchim,i)/norm
	 ENDDO
c	 PAUSE
	ENDIF	

c	estimation du pas temporel suivant

	IF(est /= 0.d0)THEN
	 IF(diffusion)THEN
	  dtnew=0.95d0*dt*precit/est
	 ELSE
	  dtnew=0.9d0*dt*(precit/est)**(1.d0/(ordre+1))	 
	 ENDIF
	ELSE
	 dtnew=1.2d0*dt
	ENDIF
	dtn=MAX(0.8d0*dt,MIN(dtnew,1.2d0*dt))
	mc_maxi=SQRT(ABS(mc_maxi))**3	!les masses avec erreur maximale
	mc_maxi=mc_maxi/mstar
	
c	écritures	

	WRITE(*,111)kmax,est,dt,dtn
	WRITE(*,112)mc_maxi(1:MIN(10,nchim)) 
	WRITE(*,112)(ab_max(i)*nucleo(i),i=1,MIN(10,nchim))
	WRITE(*,112)estim(1:MIN(10,nchim))
	WRITE(*,113)kmaxi(1:MIN(10,nchim)) 
111	FORMAT(/,'EVOL couche ',i4,', var. max. de H ou He4:',
	1 es10.3,/,'dt=',es10.3,', dt optimal=',es10.3,
	2 ' masse/abondance/variation/couche:')
112	FORMAT(12es8.1)
113	FORMAT(12(2x,i4,2x))

	dts=MIN(dtmax,dtn,agemax-age-dt) ; reprend=.FALSE.
c	PRINT*,'fin evol dt,dts,dtn,dtmax,agemax,age,agemax-age-dt'
c	WRITE(*,2000)dt,dts,dtn,dtmax,agemax,age,agemax-age-dt

c       IF(.TRUE.)PAUSE'fin de evol'
c	IF(.FALSE.)PAUSE'fin de evol'

	RETURN

	END SUBROUTINE evol
