
c********************************************************************

	SUBROUTINE evol(compt,dt,dts,reprend)

c subroutine public du module mod_evol

c gestion de l évolution temporelle de la composition chimique
c les points d intégration en comp. chim. sont les points de raccord

c sans diffusion, les ZC sont mélangées, dans ZR on conserve W
c avec diffusion, W est diffusé

c on ne tient pas compte de la différence Ptot/Pgaz négligeable
c dans les régions ou les réactions nucléaires sont actives

c on reconstruit le vecteur de composition chimique avec
c discontinuité des xchim aux limites ZR/ZC

c avec pression turbulente 8 inconnues
c sans pression turbulente 7 inconnues, Ptot=Pgaz

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c	23 09 96 : suppression du mouvement des ZR/ZC avec diffusion
c	09 10 96 : modif de la gestion de somme dm/somme r2 dm dans ZC
c	26 06 97 : remplacement du moment angulaire par la v. angulaire
c	25 08 97 : mise en place des variables eulériennes
c	20 10 99 : ajout des variables de structure au temps t dans diffus
c	19 11 99 : suppression de nh1, nhe1, nhe2, lamb
c	18 04 00 ; coeff_diff ---> diffm, difft, age dans diffus
c	30 07 00 : introduction F95

c entrées
c	compt=0: compteur des itérations globales, compt=0 pour
c	la première intération globale de chaque nouveau pas temporel
c	dt: pas temporel (fixé dans UPDATE)

c sorties
c	reprend: il faut réduire le pas temporel à cause de TdS ou de non CV
c	dts : estimation du pas temporel à utiliser pour le pas suivant

c	modèle totalement convectif: lim=1,jlim(1)=n,lconv(1)=.FALSE.
c	modèle totalement radiatif: lim=0,jlim(i)=-100,lconv(i)=.FALSE.
c	les m_zc sont en m/Mstar

c----------------------------------------------------------------

	USE mod_donnees, ONLY : ab_min, agemax, baratine, diffusion, dlntc,
	1 dtmax, dtmin, en_masse, grille_fixe, ihe4, Ipg, Kdes_rot,
	2 Krot, langue, lisse, l_demi, mdot, msol, m_rot, mvt_dis,
	3 m_ch, m_ptm, nchim, ne, nom_elem, nom_fich2, nrot, nucleo,
	4 pi, pnzc, pturb, ordre, ord_qs, ord_rot, precit, rsol, r_qs,
	5 secon6, thw, mtot
	USE mod_etat, ONLY : etat
	USE mod_kind
	USE mod_nuc, ONLY : mzc_ext, nuzc_ext
	USE mod_numerique, ONLY : bsp1dn, bsp_dis, coll, linf, newspl_gal,
	1 noeud, no_croiss, shell
	USE mod_variables, ONLY : age, bp, bp_t, chim, chim_gram, chim_t,
	1 dim_ch, dim_rot, inter, jlim, jlim_t, knot, knotc, knotc_t, knot_ptm,
	2 knot_t, knotr, knotr_t, lconv, lconv_t, lim, lim_t, mc, mct, mct_t,
	3 mc_fixe, mc_t, mrot, mrot_t, mrott, mrott_t, mstar, m_zc, m_zc_t, m23,
	4 m23_t, nc_fixe, n_ch, n_ch_t, n_ptm, n_qs, n_qs_t, n_rot,
	5 n_rot_t, old_ptm, q, qt, qt_t, q_t, r_zc, r_zc_t, r2, r2_t, rota,
	6 rota_t, sortie, tot_conv, tot_rad, xt_ptm, x_ptm

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) :: dt
	INTEGER, INTENT(in) :: compt
	REAL (kind=dp), INTENT(out) :: dts
	LOGICAL, INTENT(out) :: reprend

	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: chimd, compx,
	1 chim_zc, chim_zr, rotad, rota_tmp
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: dm, mc_tmp, mrott_tmp,
	1 Omega_zc, Omega_zr, ro, ro_t, t, t_t
	REAL (kind=dp), DIMENSION(nchim) :: ab_max, bidd,
	1 comp1, comp2, compy, dfdx, esti, estim, f, mc_maxi	
     	REAL (kind=dp), DIMENSION(lim) :: m_zc23
	REAL (kind=dp), DIMENSION(nrot) :: omega_t

	REAL (kind=dp) :: a, alfa, b, beta, bid, cp,  dcpp, dcpt, dcpx,
	1 delta, deltap, deltat, deltax, den, dgradadp, dgradadt, dgradadx,
	2 drop, drot, drox, dtn, dtnew, dtt, dup, dut, dux, est, gamma1, gradad,
	3 mass_zc, mk, mk32, mk32_t, mk_t, mnt_t, m_zc23t,
	4 norm, nu_inf, nu_sup, p, precit_eff, p_t, r, r_ext, r_t, u

	INTEGER, PARAMETER :: n_min_zc=10
	INTEGER, DIMENSION(nchim) :: kmaxi
	INTEGER, DIMENSION(1) :: ich_max
	INTEGER, SAVE :: ll, kmax
	INTEGER :: i, ich_abon_max, it, izc, j, jfin, k, knotr_tmp,
	1 nadd, nc_tmp, n_pt_zc

	LOGICAL, SAVE :: demi=.FALSE., init=.TRUE.
	LOGICAL :: discon, core_reg, l_nzc, ok, tot_conv_t

c------------------------------------------------------------------------

2000	FORMAT(8es10.3)
2001	FORMAT(es10.3,3es22.15)
2002	FORMAT(13es8.1)

c	PRINT*,'n,nchim,nchim,n_ch,knotc',n_qs,nchim,nchim,n_ch,knotc
c	PRINT*,'dt,t_inf,lim,jlim,l_conv',dt,t_inf,lim,
c	1 (jlim(i),i=1,lim),(lconv(i),i=1,lim)
c	PAUSE'entrée evol'

	IF(init)THEN
	 init=.FALSE.

	 IF(pturb)THEN
	  WRITE(*,11) ; WRITE(2,11)
11	  FORMAT(/,'Evolution temporelle de la composition chimique',/,
	1 'on ne tient pas compte de la différence Ptot-Pgaz',/,
	2 'qui est négligeable dans les régions où les',/,
	3 'réactions thermonucléaires sont actives',/)
	 ENDIF

	 SELECT CASE(Krot)
	 CASE(3,4)
	  SELECT CASE(langue)
	  CASE('english')
           WRITE(*,1022)TRIM(thw(Krot)) ; WRITE(2,1022)TRIM(thw(Krot))
1022       FORMAT(/,'Diffusion of ang. mom. with collocation',/
	1  'Formalism : ',a)
           WRITE(*,1023)m_rot ; WRITE(2,1023)m_rot
1023       FORMAT('order of B-splines:',i3)
	  CASE DEFAULT
           WRITE(*,22)TRIM(thw(Krot)) ; WRITE(2,22)TRIM(thw(Krot))
22          FORMAT(/,'Diffusion du moment cinétique par collocation.',/
	1   'Formalisme : ',a)
           WRITE(*,23)ord_rot ; WRITE(2,23)ord_rot
23          FORMAT('ordre des B-splines:',i3)
	  END SELECT
	  IF(.NOT.ALLOCATED(xcoll_rot))ALLOCATE(xcoll_rot(0))	!initialisation
	 END SELECT
	 kmax=MAX(nchim,ne)

c baratine=.FALSE. permet de dérouter sur le fichier mon_modele_103 les
c informations concernant le déroulement des calculs de l évolution de la
c composition chimique et de la vitesse angulaire	
	 IF(baratine)THEN
	  usl_evol=6
	 ELSE	  	 
	  usl_evol=103
	  OPEN(unit=103,form='formatted',status='unknown',!access='append',
	1 file=TRIM(nom_fich2)//'_evol')  
	 ENDIF
	 
c l_demi est initialisé	dans cesam.f
	 IF(grille_fixe)l_demi=.FALSE.

	ENDIF	!init

c le modèle au temps t était-il totalement convectif?
	tot_conv_t=lim_t == 1 .AND. jlim_t(1) == n_qs_t .AND. .NOT.lconv_t(1)

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

c modèle totalement convectif: lim=1, jlim(1)=n_qs, lconv(1)=.FALSE.
c nzc=1, convd(1)=1, convf(1)=nc_tmp
c modèle totalement radiatif: lim=0, jlim(i)=-100, lconv(i)=.FALSE.
c nzc=0, convd(1)=nc_tmp, convf(0)=1
	tot_conv=lim == 1 .AND. (jlim(1) == n_qs) .AND. .NOT.lconv(1)
	tot_rad=lim == 0

c adaptation de la limite de la zone externe TOUJOURS convective
	I0: IF(.NOT.tot_conv .AND. .NOT.tot_rad)THEN
	 I2: IF(.NOT.lconv(lim))THEN	!dernière limite fin de ZC
	  I1: IF(lim == 1)THEN		!limite unique
	   SELECT CASE(langue)	!que l'on déplace à la fin du modèle
	   CASE('english')	!qui devient totalement convectif
	    WRITE(usl_evol,1008)
1008	    FORMAT(/,'fully mixed model')
	   CASE DEFAULT
	    WRITE(usl_evol,8)
8	    FORMAT(/,'Modèle complètement mélangé')
	   END SELECT
	   tot_conv=.TRUE.	!rappel de la définition de tot_conv
	   lim=1 ; jlim(1)=n_qs ; lconv(1)=.FALSE.
	  ELSE I1

c la dernière limite, fin de ZC est supprimée et lim > 1
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(usl_evol,1005)lim
1005	    FORMAT(/,'extend to the surface of the outer ZC, limit #',i3)
	   CASE DEFAULT
	    WRITE(usl_evol,5)lim
5	    FORMAT(/,'prolongement à la surface de la ZC externe, limite',i3)
	   END SELECT
	   lim=lim-1
	  ENDIF I1
	 ENDIF I2 
	ENDIF I0

c limites ZR/ZC effectives
	SELECT CASE(langue)
	CASE('english')
	 WRITE(usl_evol,1006)
1006	 FORMAT(/,'limits RZ/CZ used for the convective mixing :')
	CASE DEFAULT
	 WRITE(usl_evol,6)
6	 FORMAT(/,'limites ZR/ZC utilisées pour le mélange convectif :')
	END SELECT

c cas totalement radiatif, masse fictive de la ZC externe	
	IF(tot_rad)THEN
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(usl_evol,1007)
1007	  FORMAT('Unusual condition : fully radiative model')
	 CASE DEFAULT
	  WRITE(usl_evol,7)
7	  FORMAT(/,'Situation anormale : modèle complètement radiatif')
	 END SELECT
	 lim=0 ; jlim=-100 ; lconv=.FALSE.
	 mzc_ext=1.d-20 ;  nuzc_ext=mstar**(2.d0/3.d0)

c cas totalement convectif, masse de la ZC = mstar
	ELSEIF(tot_conv)THEN
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(usl_evol,1009)
1009	  FORMAT(/'fully convective model')
	 CASE DEFAULT
	  WRITE(usl_evol,9)
9	  FORMAT(/'modèle complètement convectif')
	  mzc_ext=mstar ; nuzc_ext=0.d0
	 END SELECT

	ELSE
	 DO i=1,lim
	  IF(lconv(i))THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(usl_evol,1012)jlim(i),m_zc(i),1.d0-m_zc(i)/mstar,r_zc(i),
	1   r_zc(i)/r_ext
1012	    FORMAT('beginning of CZ shell : ',i4,/,'mass=',es10.3,
	1   ', 1-M/Mstar=',es10.3,', radius=',es10.3,', radius/Rstar=',es10.3)
	   CASE DEFAULT
	    WRITE(usl_evol,12)jlim(i),m_zc(i),1.d0-m_zc(i)/mstar,r_zc(i),
	1   r_zc(i)/r_ext
12	    FORMAT('début de ZC couche : ',i4,/,'masse=',es10.3,
	1  ', 1-m/Mstar=',es10.3,', rayon=',es10.3,', rayon/Rstar=',es10.3)
	   END SELECT
	  ELSE
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(usl_evol,1017)jlim(i),m_zc(i),1.d0-m_zc(i)/mstar,r_zc(i),
	1   r_zc(i)/r_ext
1017	    FORMAT('end of CZ shell : ',i4,/,'mass=',es10.3,
	1   ', 1-m/Mstar=',es10.3,', radius=',es10.3,', radius/Rstar=',es10.3)
	   CASE DEFAULT
	    WRITE(usl_evol,17)jlim(i),m_zc(i),1.d0-m_zc(i)/mstar,r_zc(i),
	1   r_zc(i)/r_ext
17	    FORMAT('fin de ZC couche : ',i4,/,'masse=',es10.3,
	1   ', 1-m/Mstar=',es10.3,', rayon=',es10.3,', rayon/Rstar=',es10.3)
	   END SELECT
	  ENDIF
	 ENDDO		 
	ENDIF	
c	PAUSE'lim,lconv(lim)'
	
c limites ZR/ZC en m^(2/3)
	m_zc23(1:lim)=m_zc(1:lim)**(2.d0/3.d0)

c mc_tmp : masses (m^2/3) temporaires pour tabul. de la comp.chim.
c nc_temp : nombre temporaire de masses
	IF(ALLOCATED(mc_tmp))DEALLOCATE(mc_tmp)
	ALLOCATE(mc_tmp(n_qs+n_qs_t+2*pnzc))	!allocation généreuse

c on prend la répartition déduite de bp ou la grille fixe
	IF(grille_fixe)THEN
	 mc_tmp(1:nc_fixe)=mc_fixe ; nc_tmp=nc_fixe
	ELSE
	 DO i=1,n_qs
	  CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),ll,f,dfdx)
	  IF(en_masse)THEN
	   mc_tmp(i)=f(5)
	  ELSE
	   mc_tmp(i)=ABS(f(5))**(2.d0/3.d0)
	  ENDIF
	 ENDDO
	 nc_tmp=n_qs
	ENDIF
c	PRINT*,nc_tmp ; WRITE(*,2000)mc_tmp(1:nc_tmp) ; PAUSE'1'

c un pas temporel sur 2 on prend les mc_tmp au milieu, sauf avec grille_fixe
	IF(l_demi .AND. compt == 0)demi=.NOT.demi
	IF(.NOT.grille_fixe .AND. demi)THEN
	 DO i=nc_tmp-1,2,-1
	  mc_tmp(i)=(mc_tmp(i+1)+mc_tmp(i))/2.d0
	 ENDDO
	ENDIF
c	PRINT*,nc_tmp ; PRINT*,mc_tmp(nc_tmp-4:nc_tmp) ; PAUSE'evol23'

c couches par ordre croissant
	CALL shell(nc_tmp,mc_tmp)
c	WRITE(*,*)mc_tmp(nc_tmp-4:nc_tmp) ; PAUSE'21'

c suppression des multiplicités, actualisation de nc_tmp
	k=1
	DO i=2,nc_tmp
	 IF(mc_tmp(k) /= mc_tmp(i))THEN
	  k=k+1 ; mc_tmp(k)=mc_tmp(i)
	 ENDIF
	ENDDO
	nc_tmp=k
c	PRINT*,nc_tmp ; PRINT*,mc_tmp(nc_tmp-4:nc_tmp) ; PAUSE'evol22'	

c on introduit les limites ZR/ZC en remplacement du point le plus proche
c délimitation des zones convectives
c il y a nzc zones convectives chacune entre convd(.) et convf(.)
	convd=-100 ; convf=-100
	IF(tot_conv)THEN
	 nzc=1 ; convd(1)=1 ; convf(1)=nc_tmp
	ELSEIF(tot_rad)THEN
	 nzc=0 ; convd(1)=nc_tmp ; convf(0)=1
	ELSE
	 nzc=0		!nombre de ZC
	 j=1		!indice des limites
	 B4: DO i=2,nc_tmp-1
	  b=mc_tmp(i)-m_zc23(j) ; a=m_zc23(j)-mc_tmp(i-1)
	  IF(a*b < 0.d0)CYCLE B4	   
	  IF(b < a)THEN		!point le plus proche	
	   mc_tmp(i)=m_zc23(j)
	  ELSE
	   mc_tmp(i-1)=m_zc23(j)	!si égalité coté gauche	   
	  ENDIF	   
	  IF(lconv(j))THEN		!début de ZC
	   nzc=nzc+1			!une ZC en plus
	   convd(nzc)=i
	  ELSE				!fin de ZC
	   IF(nzc == 0)THEN	!début de la ZC en 1 car pas encore identifié
	    nzc=1 ; convd(nzc)=1	!la zone externe est ZC
	   ENDIF
	   convf(nzc)=i
	  ENDIF	   
	  j=j+1
	  IF(j > lim)EXIT B4
	 ENDDO B4
	 IF(nzc >= 1 .AND. convf(nzc) < 0)convf(nzc)=nc_tmp	 
	ENDIF
	
c la ZC externe doit avoir plus de n_min_zc couches
	 nadd=n_min_zc+1+convd(nzc)-convf(nzc)
	 IF(nadd > 0)THEN
	  convd(nzc)=nc_tmp-n_min_zc
	  IF(nzc > 1 .AND. convf(nzc-1) >= convd(nzc))THEN 
	   convf(nzc-1)=convf(nzc) ; nzc=nzc-1
	  ENDIF
	 ENDIF
	 
c masse de la ZC externe en Msol
	IF(.NOT.(tot_conv .OR. tot_rad))THEN
	 nuzc_ext=mc_tmp(convd(nzc)) ; mzc_ext=mstar-m_zc(lim)
c	 WRITE(*,2000)mzc_ext,nuzc_ext,m_zc(lim)**(2.d0/3.d0) ; PAUSE'mzc_ext'
	ENDIF	 
	
c------------------------Pour vérifications-------------------------------

	IF(.FALSE.)THEN
c	IF(.TRUE.)THEN
	 PRINT*,lim,nzc,nc_tmp,mc_tmp(nc_tmp)
	 PRINT*,convd(1:nzc),convf(1:nzc),mc_tmp(convd(1:nzc)),
	1 mc_tmp(convf(1:nzc))
	 PAUSE'convd,convf'
	 ELSEIF(.FALSE.)THEN
c	 ELSEIF(.TRUE.)THEN	 
	 j=1
	 DO i=1,nc_tmp
c	  ok=lmix(mc_tmp(i))
	  IF(i == convd(j))THEN
	   IF(i /= 1)PRINT*,mc_tmp(i),i!,ok
	   PRINT*,'début de ZC',j,i ; PRINT*,mc_tmp(i),i!,ok
	  ELSEIF(i == convf(j))THEN
	   PRINT*,mc_tmp(i),i!,ok
	   PRINT*,'fin de ZC',j,i
	   IF(i /= nc_tmp)PRINT*,mc_tmp(i),i!,ok
	   j=MIN(j+1,nzc)
	  ELSE
	   PRINT*,mc_tmp(i),i!,ok
	  ENDIF
	 ENDDO
	 PAUSE'avant ajustement des ZR et ZC'
	ENDIF

c------------------------------------------------------------------------		

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
	 WRITE(usl_evol,1002)nc_tmp
1002	 FORMAT('number of shells for the interpolation of chemicals:',i5)
	CASE DEFAULT
	 WRITE(usl_evol,2)nc_tmp
2	 FORMAT('nombre de couches utilisées pour la comp. chim.:',i5)
	END SELECT

c on fixe le centre
	mc_tmp(1)=0.d0

c le vecteur de mélange pour la composition chimique
c le mélange convectif doit affecter toute la couche contenant une limite ZR/ZC
c avec rotation on fixe la grille à celle de la 1-ère itération
c sans rotation on adapte la grille
	IF(convd(1) == 1)THEN
	 n_mix=2*nzc ; j=0
	ELSE
	 n_mix=2*nzc+1 ; j=1
	 x_mix(1)=mc_tmp(1) ; mix(1)=.FALSE.
	ENDIF
	DO i=1,nzc
	 j=j+1
	 x_mix(j)=mc_tmp(convd(i)) ; mix(j)=.TRUE.
	 j=j+1
	 x_mix(j)=mc_tmp(convf(i)) ; mix(j)=.FALSE.
	ENDDO

c estimation de ndis et idis : nombre et tableau des discontinuités
	idis=-HUGE(1)
	IF(tot_conv)THEN
	 ndis=0
	ELSE
	 ndis=2*nzc
	 IF(convd(1) == 1)ndis=ndis-1
	 IF(convf(nzc) == nc_tmp)ndis=ndis-1
	 idis(0)=1 ; idis(ndis+1)=nc_tmp
	 j=0
	 DO i=1,nzc
	  IF(convd(i) /= 1)THEN
	   j=j+1 ; idis(j)=convd(i)
	  ENDIF
	  IF(convf(i) /= n_mix)THEN
	   j=j+1 ; idis(j)=convf(i)
	  ENDIF
	 ENDDO
	ENDIF
c	PRINT*,nzc,ndis,idis ; PAUSE'ndis'

c---------------------Pour vérifications-----------------------------

c	IF(.TRUE.)THEN
	IF(.FALSE.)THEN
	 PRINT*,'nzc',nzc,nc_tmp
	 PRINT*,'convd',convd(1:nzc+1)
	 PRINT*,'convf',convf(0:nzc)
	 DO i=1,nzc
	  PRINT*,i,convd(i),convf(i),mc_tmp(convd(i)),mc_tmp(convf(i))
	 ENDDO
	 PRINT*
	 j=1
	 DO i=1,nc_tmp
	  ok=lmix(mc_tmp(i))	  
	  IF(i == convd(j))THEN
	   IF(i /= 1)PRINT*,mc_tmp(i),i,mc_tmp(i)**(3.d0/2.d0),ok
	   PRINT*,'début de ZC',j,i
	   PRINT*,mc_tmp(i),i,mc_tmp(i)**(3.d0/2.d0),ok
	  ELSEIF(i == convf(j))THEN
	   PRINT*,mc_tmp(i),i,mc_tmp(i)**(3.d0/2.d0),ok
	   PRINT*,'fin de ZC',j,i
	   IF(i /= nc_tmp)PRINT*,mc_tmp(i),i,mc_tmp(i)**(3.d0/2.d0),ok
	   j=MIN(j+1,nzc)
	  ELSE
	   PRINT*,mc_tmp(i),i,mc_tmp(i)**(3.d0/2.d0),ok
	  ENDIF
	 ENDDO
	 PAUSE'après ajustement des ZR et des ZC'
	ENDIF

c-----------------------  diffusion des éléments chimiques-----------

	IF(diffusion)THEN
	 CALL diffus(ok,dt,mc_tmp,nc_tmp)
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
c éventuellement forte de comp.chim. on élimine les ZC du test
c	 WRITE(*,2002)(mct_t(i),i=1,knotc_t)
c	 PAUSE'estimation de la précision'

c cas totalement convectif, on se place en un point test (en mc(5))
	 IF(tot_conv)THEN
	  it=5 ; jfin=it

c sinon tout le modèle est balayé
	 ELSE
	  it=1 ; jfin=n_ch
	 ENDIF

	 DO i=it,jfin
	  CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,knotc,.TRUE.,mc(i),
	1 ll,comp1,bidd)
	  IF(en_masse)THEN
	   bid=MIN(mc(i),x_ptm(n_ptm))
	  ELSE
	   bid=MIN(SQRT(ABS(mc(i)))**3,x_ptm(n_ptm))
	  ENDIF
	  CALL bsp1dn(1,old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm,.TRUE.,
	1 bid,ll,f,dfdx)	!masse au temps t mk_t
	  IF(.NOT.en_masse)f(1)=ABS(f(1))**(2.d0/3.d0)
	  mk_t=MIN(MAX(mc_t(1),f(1)),mc_t(n_ch_t))
	  CALL bsp1dn(nchim,chim_t,mc_t,mct_t,n_ch_t,m_ch,
	1 knotc_t,.TRUE.,mk_t,ll,comp2,bidd)

c indice de l élément le plus abondant en mole
	  ich_max=MAXLOC(comp2)

c estimation de la variation maximale de l élément localement le plus abondant
	  DO j=1,nchim
	   bid=MAX(ABS(comp1(j)),ABS(comp2(j)))
	   IF(bid > 0.d0)THEN
	    esti(j)=ABS(comp1(j)-comp2(j))/bid
	   ELSE
	    esti(j)=0.d0
	   ENDIF
	   IF(estim(j) <= esti(j))THEN		!précision par élément
	    estim(j)=esti(j) ; kmaxi(j)=i
	    ab_max(j)=comp1(j) ; mc_maxi(j)=mc(i)
	   ENDIF	   
	  ENDDO

c précision globale estimée sur l élément le plus abondant en mole
	  est=MAX(est,esti(ich_max(1)))
	  IF(est == esti(ich_max(1)))THEN
	   kmax=i ; ich_abon_max=ich_max(1)
	  ENDIF
	 ENDDO

c	 PAUSE'après diffusion'

c----------------sans diffusion mélange classique début------------------------------

	ELSE

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
c	 d abord les ZC puis les ZR
c	 on garde dans chim_zc(nchim,nzc) les valeurs des ZM
c	 la variable en masse est mc=m**2/3
c	 convf(0)=1 si convd(1) /= 1, sinon convf(0)=10000
c	 convd(nzc+1)=nc_tmp si convf(nzc) /= nc_tmp,
c	 sinon convd(nzc+1)=-10000
c	 dans ZM on mélange MW, dans ZR on conserve MW
c	 on reconstruit le vecteur de composition chimique avec
c	 discontinuité des xchim aux LMR.

c les ZC d abord
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
c on tient compte de la perte/gain de masse
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
	   mk_t=MIN(MAX(mc_t(1),bid),mc_t(n_ch_t))
	   mk32_t=SQRT(ABS(mk_t))**3
	   CALL bsp1dn(nchim,chim_t,mc_t,mct_t,n_ch_t,m_ch,
	1  knotc_t,.TRUE.,mk_t,ll,comp1,bidd)

c	   PRINT*,'après bsp1dn 1, izc,i,k',izc,i,k
c	   WRITE(*,2000)mk,(comp1(j),j=1,MIN(7,nchim))
c	   PAUSE'après bsp1dn 1'

c comp1, compx : abondances au temps t
c comp2, compy : abondances au temps t+dt
c compx, compy : abondances par mole, comp1, comp2 par gramme
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

c comp1 et comp2 par gramme pour appel à l EOS
	   CALL chim_gram(comp1,bidd) ; CALL chim_gram(comp2,bidd)

c T et ro en t+dt par interpolation inverse en m23 par inter
c	   PRINT*,'avant, izc,i,k',izc,i,k
	   IF(en_masse)THEN
	    CALL inter('m23',bp,q,qt,n_qs,knot,MIN(mk,m23(n_qs)),f,dfdx,r2,m23)
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
	1   MIN(mk_t,m23_t(n_qs_t)),f,dfdx,r2_t,m23_t)
	   ELSE
	    CALL inter('m23',bp_t,q_t,qt_t,n_qs_t,knot_t,
	1   MIN(mk32_t,m23_t(n_qs_t)),f,dfdx,r2_t,m23_t)
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

c l_nzc : on est dans la ZC externe
	   l_nzc= izc == nzc
	  ENDDO		!i dans chaque ZC

c	  PRINT*,k,n_pt_zc ; PAUSE'evol: n_pt_zc'
c	  WRITE(*,2000)compx(1,1:k)
c	  WRITE(*,2000)dm ; WRITE(*,2000)t
c	  WRITE(*,2000)ro ; WRITE(*,2000)mass_zc ; PAUSE'après les dm(k)'

c normalisation des dm
	  dm=dm/mass_zc

c intégration temporelle
	  CALL rk_imps(t_t,ro_t,compx,t,ro,compy,dt,esti,ok,n_pt_zc,dm,l_nzc)
c	  PRINT*,ok ; WRITE(*,2000)compy ; PAUSE'après rk_imps'

c  on garde, provisoirement, les abondances dans chim_zc
	  IF(ok)THEN		!il y a convergence
	   chim_zc(:,izc)=compy(:) !chim_zc: comp. chim. de la ZC

c gestion de la précision locale
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

c indice de l élément le plus abondant en mole
	   ich_max=MAXLOC(compy)

c précision globale estimée sur l élément le plus abondant en mole
	   est=MAX(est,esti(ich_max(1)))
	   IF(est == esti(ich_max(1)))THEN
	    kmax=convd(izc) ; ich_abon_max=ich_max(1)
	   ENDIF
	    	  	   
c difficulté	   
	  ELSE
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,2020)izc ; WRITE(2,2020)izc
2020	    FORMAT('EVOL/rk_imps Pb. at ZC #',i4,' Reboot and dt-->dt/2',/)
	   CASE DEFAULT
	    WRITE(*,20)izc ; WRITE(2,20)izc
20	    FORMAT('EVOL/rk_imps Pb. à la ZC #',i4,
	1   ' Réinitialisation et dt-->dt/2',/)
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
	 dm=1.d0 ; l_nzc=.FALSE.
	 DO izc=1,nzc+1
c	  PRINT*,'ZR',izc,convf(izc-1),convd(izc)
	  DO i=convf(izc-1),convd(izc)

c la comp.chim. en mk_t au temps t et en mc_tmp(i) au temps t+dt
c les mc peuvent différer des mct
	   IF(en_masse)THEN
	    bid=mc_tmp(i)
	   ELSE
	    bid=SQRT(ABS(mc_tmp(i)))**3
	   ENDIF
	   CALL inter('m23',bp,q,qt,n_qs,knot,MIN(bid,m23(n_qs)),f,dfdx,r2,m23)
	   bid=MIN(bid,x_ptm(n_ptm))
	   CALL bsp1dn(1,old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm,.TRUE.,
	1  bid,ll,f,dfdx)	!masse au temps t
	   IF(.NOT.en_masse)f(1)=ABS(f(1))**(2.d0/3.d0)
	   mk_t=MIN(MAX(mc_t(1),f(1)),mc_t(n_ch_t)) ; mk32_t=SQRT(mk_t)**3
	   CALL bsp1dn(nchim,chim_t,mc_t,mct_t,n_ch_t,m_ch,
	1  knotc_t,.TRUE.,mk_t,ll,comp1,bidd)	!perte de Mw

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

c les T et ro en t+dt, par interpolation inverse en m23
c	   PRINT*,'avant, izc,i',izc,i
	   IF(en_masse)THEN
	    bid=mc_tmp(i)
	   ELSE
	    bid=SQRT(ABS(mc_tmp(i)))**3
	   ENDIF
	   CALL inter('m23',bp,q,qt,n_qs,knot,MIN(bid,m23(n_qs)),f,dfdx,r2,m23)
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

c les T et ro en t, par interpolation inverse en m23
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

c intégration temporelle
	   CALL rk_imps(t_t,ro_t,compx,t,ro,compy,dt,esti,ok,1,dm,l_nzc)
	   compy=ABS(compy)
c	   WRITE(*,2000)mc_tmp(i),t,ro,compy(1),compx(1,1)
c	   PRINT*,'dans ZR',i

c on place les abondances dans chim_zr
	   IF(ok)THEN		!il y a eu convergence
	    chim_zr(:,i)=compy

c gestion de la précision locale
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

c indice de l'élément le plus abondant en mole
	    ich_max=MAXLOC(compy)

c précision globale estimée sur l'élément le plus abondant en mole
	    est=MAX(est,esti(ich_max(1)))
	    IF(est == esti(ich_max(1)))THEN
	     kmax=i ; ich_abon_max=ich_max(1)
	    ENDIF

c difficulté
	   ELSE			!ok
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,2021)i ; WRITE(2,2021)i
2021	     FORMAT('EVOL/rk_imps Pb. at the shell #',i4,
	1     ', Reboot and dt-->dt/2')
	    CASE DEFAULT
	     WRITE(*,21)i ; WRITE(2,21)i
21	     FORMAT('EVOL/rk_imps Pb. à la couche #',i4,
	1     ', Réinitialisation et dt-->dt/2')
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

c------- mouvements des limites ZR/ZC si mvt_dis=.TRUE. ----------------

c le modèle et le précédent ne peuvent être totalement convectifs
	 discon=.NOT.(tot_conv_t .AND. tot_conv) .AND. mvt_dis 

c cas de la régression d'un coeur convectif :
c sur la 1-ière limite, ie. en jlim(1), il y passage de
c ZC à ZR, donc lconv(1) est FALSE
	 core_reg=discon .AND. (.NOT.lconv(1) .AND. .NOT.lconv_t(1)
	1 .AND. m_zc(1) < m_zc_t(1) .AND. m_zc(1) < mstar)

c sans rotation	
	 core_reg=core_reg .AND. Krot == 0

c le noyau convectif a régressé de nu_sup=m_zc23t à nu_inf=mc_tmp(convf(1)) 
	 IF(core_reg)THEN 
	  m_zc23t=m_zc_t(1)**(2.d0/3.d0)
	  nu_inf=mc_tmp(convf(1)) ; nu_sup=m_zc23t

c écritures
	 SELECT CASE(langue)
	  CASE('english')
	   WRITE(usl_evol,1025)SQRT(nu_sup/mc_t(n_ch_t))**3,
	1  SQRT(nu_inf/mc_t(n_ch_t))**3
1025	   FORMAT('The convective core recedes from',es10.3' to',es10.3,
	1  ' Mstar')
	  CASE DEFAULT
	   WRITE(usl_evol,25)SQRT(nu_sup/mc_t(n_ch_t))**3,
	1  SQRT(nu_inf/mc_t(n_ch_t))**3
25	   FORMAT('Le noyau convectif a régressé de',es10.3,' à' ,
	1  es10.3,' Mstar')
	  END SELECT

c avec régression, à gauche de nu_inf on met la comp. chim. de la ZC
	  chim_zr(:,convf(1))=chim_zc(:,1)	  
	  CALL linf(nu_sup,mc_tmp,nc_tmp,ll)

c la ZC a regressé de m_zc23t à mc_tmp(convf(1)) et on a :
c mc_tmp(ll+1) > m_zc23t >= mc_tmp(ll) >= mc_tmp(convf(1))
c on se place en ll+1 pour être à droite de la limite située dans [ll,ll+1[
	  DO i=convf(1),ll+1
	   bid=(mc_tmp(i)-nu_inf)/(mc_tmp(ll+1)-nu_inf)
	   chim_zr(:,i)=(chim_zr(:,ll+1)-chim_zr(:,convf(1)))*bid
	1  +chim_zr(:,convf(1))
c	    PRINT*,convf(1),i,ll+1
c	    WRITE(*,2000)nu_inf,mc_tmp(i),mc_tmp(ll+1)
c	    WRITE(*,2000)chim_zr(1,convf(1)),chim_zr(1,i),chim_zr(1,ll+1)	
	  ENDDO
	 ENDIF

c-------------------zc_aug éliminé---------------------------------
c cas de l'augmentation de l'abscisse lagrangienne (masse) de la limite
c inférieure de la ZC externe, interp. linéaire la comp.chim. sur la zone
c d'accroissement. Il reste une discontinuité, au moins formelle, de comp.chim.
c	  m_zc23t=m_zc_t(lim_t)**(2.d0/3.d0)
c	  zc_aug=discon .AND. m_zc23t < mc_tmp(convd(nzc))  	  	  	  
c	  IF(zc_aug)THEN
c	   CALL linf(m_zc23t,mc_tmp,nc_tmp,ll)
c mc_tmp(ll) <= m_zc23t < mc_tmp(convd(nzc))
c la lim. de la ZC est passée de m_zc23t=mc_tmp(ll) à mc_tmp(convd(nzc))
c interp. lin. de la comp. chim. pour mc entre mc_tmp(ll) et mc_tmp(convd(nzc))
c	   DO i=ll,convd(nzc)
c	    bid=(mc_tmp(i)-mc_tmp(convd(nzc)))/(mc_tmp(ll)-mc_tmp(convd(nzc)))
c	    chim_zr(:,i)=(chim_zr(:,ll)-chim_zr(:,convd(nzc)))*bid
c	1   +chim_zr(:,convd(nzc))
c	   ENDDO
c	  ENDIF	!zc_aug
c-------------------zc_aug éliminé---------------------------------

c si on tient compte du mouvement des limites ZR/ZC
	 Imvt: IF(discon)THEN

c on modifie le vecteur de composition chimique
c	  convf(0)=1 si convd(1) /= 1, sinon convf(0)=10000
c	  convd(nzc+1)=nc_tmp si convf(nzc) /= nc_tmp,
c	  sinon convd(nzc+1)=-10000

c	  PRINT*,nzc;  PRINT*,(convf(i),i=0,nzc)
c	  PRINT*,(convd(i),i=1,nzc+1) ; PAUSE'nzc'

c redéfinition du vecteur de composition chimique nc_tmp-->n_ch
c chimd : abondances aux discontinuités, idis(ndis): indices discon.
c la seconde dim. du vecteur de comp.chim. chim(nchim,n_ch+ndis)
c tient compte des discontinuités, la dimension du vecteur nodal
c correspondant sera knotc=n_ch+m_ch+ndis, dim_ch=n_ch+ndis
c transfert des masses temporaires
c la dimension du vecteur dual mct est knotc=n_ch+m_ch+ndis
	  n_ch=nc_tmp ; DEALLOCATE(chim,mc,mct) ; ALLOCATE(mc(n_ch))
	  mc(1:n_ch)=mc_tmp(1:n_ch)
	  ALLOCATE(chim(nchim,n_ch+ndis),chimd(nchim,ndis),mct(n_ch+m_ch+ndis))

c j est l'indice de la discontinuité, à la fin de la séquence j=ndis
c coté gauche on place dans chim, suivant le cas, chim_zr ou chim_zc
c coté droit on place dans chimd, suivant le cas chim_zc ou chim_zr 
	  j=0
	  DO izc=1,nzc+1	  
c	   PRINT*,'izc,convf(izc-1),convd(izc)',izc,convf(izc-1),convd(izc)
	   DO i=convf(izc-1),convd(izc)	!zone radiative  ZR
	    IF(i == convf(izc-1) .AND. i > 1)THEN 
	     IF(i == idis(j+1))THEN
	      j=j+1 ; chimd(:,j)=chim_zr(:,i)
c	      PRINT*,'j,i',j,i ; PRINT*,mc(i)
c	      PAUSE'début de ZR à droite de ZM'
	     ENDIF
	    ELSE
	     chim(:,i)=chim_zr(:,i)
	    ENDIF
	   ENDDO

	   IF(izc <= nzc)THEN		!pour ZM	   
c	    PRINT*,'izc,convd(izc),convf(izc)',izc,convd(izc),convf(izc)
	    DO i=convd(izc),convf(izc)
	     IF(i == convd(izc) .AND. i > 1)THEN
	      IF(i == idis(j+1))THEN	     	     
	       j=j+1 ; chimd(:,j)=chim_zc(:,izc)	      
c	       PRINT*,'j,i,izc',j,i,izc
c	       WRITE(*,2000)chim_zc(1,izc),chimd(1,j)
c	       PAUSE'début de ZM à droite de ZR' 
	      ENDIF
	     ELSE
	      chim(:,i)=chim_zc(:,izc)
	     ENDIF	!i == convd
	    ENDDO	!i=convd , convf
	   ENDIF	!izc <= nzc
	  ENDDO		!nzc  	   	  

c vérification: on doit avoir j=ndis
	  IF(j /= ndis)THEN
	   PRINT*,'evol: j=',j,' /= ndis=',ndis ; PRINT*,'ARRET'
	   CALL sortie
	  ENDIF

c	  PRINT*,n_ch,nzc,ndis,idis
c	  ll=1
c	  DO i=1,n_ch
c	   WRITE(*,2000)mc(i),chim(1:MIN(4,nchim),i)
c	   IF(i == idis(ll))THEN
c	    PRINT*,idis(ll) ; WRITE(*,2000)mc,chimd(1:MIN(4,nchim),ll)
c	    ll=MIN(ll+1,ndis)
c	   ENDIF
c	  ENDDO
c	  DO i=1,ndis
c	   WRITE(*,2000)chimd(1:MIN(4,nchim),i)
c	  ENDDO
c	  PAUSE'après écritures'

c tabulation de la composition chimique, bsp_dis envisage les cas
c avec et sans discontinuité, avec lisse=.TRUE. lissage par contour
c	  PRINT*,nchim,ndis,idis,n_ch,m_ch,mvt_dis
c	  PRINT*,'mc',mc(n_ch)
c	  WRITE(*,2000)mc ; WRITE(*,2000)chim(1,1:n_ch)
c	  WRITE(*,2000)chimd(1,1:ndis) ;	 WRITE(*,2000)mc(300:350)

c	  PRINT*,convf(1),ll
c	  WRITE(*,2000)chim(1,convf(1)-1:ll+1) ; PAUSE'bsp_dis'

	  CALL bsp_dis(nchim,mc,chim,ndis,idis,chimd,eps,n_ch,m_ch,mct,knotc,
	1 lisse)

          IF(no_croiss)THEN
           PRINT*,'Arrêt 1 dans evol' ; CALL sortie
          ENDIF
	  dim_ch=knotc-m_ch ; DEALLOCATE(chimd)

c	  PRINT*,'après intégration : n_ch,knotc',n_ch,knotc
c	  PRINT*,'mc',mc(n_ch) ; WRITE(*,2000)mc
c	  PRINT*,'mct',mct(knotc),mc(n_ch)
c	  WRITE(*,2000)mct ; PRINT*,'xchim'
c	  DO i=1,nc
c	   CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
c	1  knotc,.TRUE.,mc(i),ll,compx,compy)
c	   WRITE(*,2000)mc(i),(compx(j)*nucleo(j),j=1,MIN(7,nchim))
c	   WRITE(*,2002)mc(i),(compx(j)*nucleo(j),j=1,nchim)
c	   WRITE(*,2002)mc(i),(compx(j),j=1,nchim)
c	  ENDDO
c	  PAUSE'après intégration'

c	  PAUSE'test rota'
c	  DO i=1,nc
c	   CALL inter('m23',bp,q,qt,n_qs,knot,MIN(m23(i),m23(n_qs)),
c	1  f,dfdx,r2,m23)
c	   CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
c	1  knotc,.TRUE.,m23(i),ll,compx,compy)
c	   IF(r2(i) > 0.d0)PRINT*,m23(i),r2(i),,compx(1)
c	   IF(f(3) > 0.d0)THEN
c	    PRINT*,m23(i),f(3)
c	   ELSE
c	    PRINT*,m23(i),f(3)
c	   ENDIF
c	  ENDDO
c	  PAUSE'après test'

	 ELSE Imvt	!sans discontinuités

c on place les abondances des ZR et des ZC dans chim sans tenir compte des
c discontinuités
c [convd(i-1) ZC [convf(i-1) ZR [convd(i) ZC [convf(i) ZR [convd(i+1)....
c .... ZR [convd(nzc) ZC convf(nzc)=nc_tmp], i=1,nzc
c Si coeur convectif conf(0)=10000, sinon conf(0)=1
	  n_ch=nc_tmp ; DEALLOCATE(chim,mc,mct) ; ALLOCATE(mc(n_ch))
	  mc(1:n_ch)=mc_tmp(1:n_ch)
	  ALLOCATE(chim(nchim,n_ch),mct(n_ch+m_ch))
	  DO izc=1,nzc+1
	   DO i=convf(izc-1),convd(izc)-1	!zone radiative
	    chim(:,i)=chim_zr(:,i)		![convf(i-1) ZR [convd(i)
	   ENDDO

	   IF(izc <= nzc)THEN
	    DO i=convd(izc),convf(izc)-1	!zone convective
	     chim(:,i)=chim_zc(:,izc)		![convd(i) ZC [convf(i)
	    ENDDO	    
	   ENDIF	!izc <= nzc
	  ENDDO		!nzc
	  chim(:,n_ch)=chim_zc(:,nzc)

c interpolation avec, éventuellement, lissage
	  CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,knotc,.FALSE.,mc(1),ll,
	1 comp1,bidd,lisse)
	  dim_ch=knotc-m_ch
	 ENDIF	Imvt	!discontinuités

	 DEALLOCATE(compx,chim_zc,chim_zr,dm,ro,ro_t,t,t_t)
	ENDIF
	
c---------diffusion ou intégration et mélange classiques fin------------

c normalisation: afin d'obtenir Somme Xi=1
c on compense le fait que les nucleo i.e. les masses des éléments
c chimiques ne sont pas des nombres entiers
c pour chaque élément, dans la base des splines, on forme
c norn = somme des coefficients X nucleo, norm est proche de 1
c on divise ensuite les coefficients par norm
c ce qui est valable parceque la somme des B-splines normalisées = 1
	IF(nchim > 1)THEN
	 DO i=1,dim_ch	!dim_ch=knotc-m_ch: dimension de la base
	  norm=SUM(chim(1:nchim,i)*nucleo(1:nchim))
	  chim(1:nchim,i)=chim(1:nchim,i)/norm
	 ENDDO
c	 PAUSE
	ENDIF
		
c~~~~le vecteur nodal pour les équations différentielles de la rotation~~~~~~

c formation des points de collocation, qui serviront pour rota et poisson
c pour un ordre d'interpolation m_rot+r_qs=ord_rot, il y a m_rot points de
c collocation dans chaque intervalle, la multiplicité des noeuds étant mrot
c il y a ncoll_rot=(n_rot-1)*m_rot points de collocation et r_qs points limite
c r_qs ordre des équations différentielles
c ord_rot=m_rot+r_qs ordre d'interpolation
c avec nd discontinuités il convient d'ajouter nzc*r_qs points limite
c soit (nd+1)*r_qs points limite
c sans discontinuité interne dimension du vecteur nodal knotr=dim_rot+m_rot+r_qs
c avec discontinuités internes knotr=dim_rot+m_rot+r_qs(nzc+1)
c dim_rot, knotr,  éléments de mod_variables, définis dans base_rota
c ord_rot est initialisé dans ini_rota,
c m_rot, r_qs sont des éléments de mod_donnees
c ce n'est qu'au second passage que les modifications de masse totale
c (planétoïdes, perte de masse...) sont pris en compte
	
c-------------------- diffusion du moment cinétique---------------------

	SELECT CASE(Krot)
	CASE(3,4)
	 
c interpolation de variables thermodynamiques (calcul de dro/dR)
	 CALL tab_vth	 

c mrot --> mrot_tmp
	 ALLOCATE(mrott_tmp(SIZE(mrott)),rota_tmp(nrot,SIZE(rota,2)))
	 rota_tmp=rota ; mrott_tmp=mrott ; knotr_tmp=knotr ; n_rot=nc_tmp

c mc --> mrot
	 DEALLOCATE(mrot) ; ALLOCATE(mrot(n_rot))
	 mrot=mc_tmp(1:n_rot) ; mrot(1)=0.d0

c génération de la base de de Boor pour la rotation
	 CALL base_rota

c spline sur le nouveau vecteur nodal mrott
	 DEALLOCATE(rota) ;  ALLOCATE(rota(nrot,dim_rot))

c	 CALL newspl(nrot,mrot,mrott_tmp,knotr_tmp,ord_rot,mrott,knotr,ord_rot,
c	1 rota_tmp,rota,.TRUE.)	
	 CALL newspl_gal(nrot,mrot,mrott_tmp,knotr_tmp,ord_rot,rota_tmp,
	1 mrott,knotr,ord_rot,rota)	
	 IF(no_croiss)PRINT*,'Pb. en 2 dans evol'

c les points de collocation
	  ncoll_rot=(n_rot-1)*m_rot
	  DEALLOCATE(xcoll_rot) ; ALLOCATE(xcoll_rot(ncoll_rot))
	  CALL coll(mrot,n_rot,m_rot,xcoll_rot)
	 
c nettoyage
	  DEALLOCATE(mrott_tmp,rota_tmp)

c intégration
	 PRINT*
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(usl_evol,*)'---- diffusion of angular momentun (begin) ----'
	 CASE DEFAULT
	  WRITE(usl_evol,*)'---- diffusion du moment cinétique (début) ----'
	 END SELECT

c résolution des équations
	CALL resout_rota(dt,ok)	
	reprend=.NOT.ok ; IF(reprend)RETURN	 

c écriture des variables de la rotation off line
	 IF(Kdes_rot == 3)CALL ecrit_rota(dt)
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(usl_evol,*)'---- diffusion of angular momentun (end) ----'
	 CASE DEFAULT
	  WRITE(usl_evol,*)'---- diffusion du moment cinétique (fin) ----'
	 END SELECT

c ~~~~~~~~~~~conservation locale du moment cinétique~~~~~~~~~~~~~~~~

	CASE(5)

c mélange dans les ZC
	 ALLOCATE(Omega_zc(nzc))
	 DO izc=1,nzc
	  den=0.d0 ; mnt_t=0.d0

c n_pt_zc : nombre de points dans la ZC
	  n_pt_zc=convf(izc)-convd(izc) ; ALLOCATE(dm(n_pt_zc))

	  DO i=convd(izc),convf(izc)-1
	   k=i-convd(izc)+1
	   mk=(mc_tmp(i)+mc_tmp(i+1))/2.d0		!point milieu
	   mk32=SQRT(ABS(mk))**3 ; dm(k)=(mc_tmp(i+1)-mc_tmp(i))*SQRT(mk)

c omega_t en mk_t au temps t
c on tient compte de la perte/gain de masse
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
	   mk_t=MIN(MAX(mrot_t(1),bid),mrot_t(n_rot_t))
	   mk32_t=SQRT(ABS(mk_t))**3
	   CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,
	1  mk_t,ll,omega_t,dfdx)

c r_t en mk_t au temps t par interpolation inverse en m23
	   IF(en_masse)THEN
	    CALL inter('m23',bp_t,q_t,qt_t,n_qs_t,knot_t,
	1   MIN(mk_t,m23_t(n_qs_t)),f,dfdx,r2_t,m23_t)
	    r_t=SQRT(f(3))
	   ELSE
	    CALL inter('m23',bp_t,q_t,qt_t,n_qs_t,knot_t,
	1   MIN(mk32_t,m23_t(n_qs_t)),f,dfdx,r2_t,m23_t)
	    r_t=f(3)
	   ENDIF
	   mnt_t=mnt_t+r_t**2*omega_t(1)*dm(k)

c r^2 en mk au temps t+dt par interpolation inverse en m23
	   IF(en_masse)THEN
	    CALL inter('m23',bp,q,qt,n_qs,knot,MIN(mk,m23(n_qs)),f,dfdx,r2,m23)
	    r=SQRT(f(3))
	   ELSE
	    CALL inter('m23',bp,q,qt,n_qs,knot,MIN(mk32,m23(n_qs)),f,dfdx,
	1   r2,m23)
	    r=f(3)
	   ENDIF
	   den=den+r*(3.d0*r-2.d0*r_t)*dm(k)
	  ENDDO			!i dans chaque ZC

c Oméga moyen pour la ZC
	  Omega_zc(izc)=mnt_t/den ; DEALLOCATE(dm)
	 ENDDO

c ensuite les ZR
	 ALLOCATE(Omega_zr(nc_tmp)) ; Omega_zr=0.d0

	 DO izc=1,nzc+1
	  DO i=convf(izc-1),convd(izc)

c centre radiatif
	   IF(i == 1)THEN
	    Omega_zr(i)=-100.d0

c ailleurs
	   ELSE
	    IF(en_masse)THEN
	     bid=mc_tmp(i)
	    ELSE
	     bid=SQRT(ABS(mc_tmp(i)))**3
	    ENDIF
	    CALL inter('m23',bp,q,qt,n_qs,knot,MIN(bid,m23(n_qs)),f,dfdx,r2,m23)
	    IF(en_masse)THEN
	     r=SQRT(f(3))
	    ELSE
	     r=f(3)
	    ENDIF
	    bid=MIN(bid,x_ptm(n_ptm))
	    CALL bsp1dn(1,old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm,.TRUE.,
	1   bid,ll,f,dfdx)	!masse au temps t

c omega_t au temps t
	    IF(.NOT.en_masse)f(1)=ABS(f(1))**(2.d0/3.d0)
	    mk_t=MIN(MAX(mc_t(1),f(1)),mc_t(n_ch_t))
	    CALL bsp1dn(nrot,rota_t,mrot_t,mrott_t,n_rot_t,ord_rot,knotr_t,
	1   .TRUE.,mk_t,ll,omega_t,dfdx)

c r_t au temps t
	    CALL inter('m23',bp_t,q_t,qt_t,n_qs_t,knot_t,MIN(bid,m23_t(n_qs_t)),
	1   f,dfdx,r2_t,m23_t)
	    IF(en_masse)THEN
	     r_t=SQRT(f(3))
	    ELSE
	     r_t=f(3)
	    ENDIF

c Oméga en mc_tmp(i)
	    Omega_zr(i)=(r*Omega_t(1))/(3.d0*r-2.d0*r_t)
	   ENDIF
	  ENDDO
	 ENDDO

c cas du centre radiatif
	 IF(Omega_zr(1) < 0.d0)Omega_zr(1)=Omega_zr(2)

c	 WRITE(*,2000)Omega_zc ; WRITE(*,2000)Omega_zr(1:4) ;  PAUSE'Krot=5'

c interpolation de rota
	 n_rot=nc_tmp ; DEALLOCATE(rota)
	 ALLOCATE(rota(nrot,n_rot+ndis),rotad(nrot,ndis))

c j est l'indice de la discontinuité, à la fin de la séquence j=ndis
	 j =0
	 DO izc=1,nzc+1
	  DO i=convf(izc-1),convd(izc)	!zone radiative  ZR
	   IF(i == convf(izc-1) .AND. i > 1)THEN
	    j=j+1
	    rotad(1,j)=Omega_zr(i)
	   ELSE
	    rota(1,i)=Omega_zr(i)
	   ENDIF
	  ENDDO

	  IF(izc <= nzc)THEN		!pour ZM
	   DO i=convd(izc),convf(izc)
	    IF(i == convd(izc) .AND. i > 1)THEN
	     j=j+1
	     rotad(1,j)=Omega_zc(izc)
	    ELSE
	     rota(1,i)=Omega_zc(izc)
	    ENDIF	!i == convd
	   ENDDO	!i=convd , convf
	  ENDIF		!izc <= nzc
	 ENDDO		!nzc

	 DEALLOCATE(mrot,mrott) ; ALLOCATE(mrot(n_rot),mrott(n_rot+ord_rot+ndis))
	 mrot(1:n_rot)=mc_tmp(1:n_ch)

	 CALL bsp_dis(nrot,mrot,rota,ndis,idis,rotad,eps,n_rot,ord_rot,mrott,
	1 knotr,mvt_dis)
	 dim_rot=knotr-ord_rot

	 DEALLOCATE(omega_zc,omega_zr,rotad)

	END SELECT

	DEALLOCATE(mc_tmp)

c---------------------estimation du pas temporel suivant-------------------

c approche de la fin de la MS, dt < 0.01
	I4: IF(chim_t(1,1) > ab_min(1) .AND. chim(1,1) < ab_min(1))THEN
	 reprend=chim_t(1,1)-chim(1,1) > 0.05d0	 
	 I5: IF(reprend)THEN
	  SELECT CASE(langue)
	   CASE('english')
	   WRITE(*,1030) ; WRITE(2,1030)
1030	   FORMAT('Too fast approach of the end of the MS, dt-->dt/2')
	  CASE DEFAULT
	   WRITE(*,30) ; WRITE(2,30)
30	   FORMAT('Approche trop rapide de la fin de la MS, dt-->dt/2')
	  END SELECT
	  I6: IF(dt < dtmin)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1031)dt,dtmin ; WRITE(2,1031)dt,dtmin
1031	    FORMAT('STOP, dt =,'es10.3,' <',es10.3,' = dtmin')
	   CASE DEFAULT
	    WRITE(*,31)dt,dtmin ; WRITE(2,31)dt,dtmin
31	    FORMAT('ARRET, dt =,'es10.3,' <',es10.3,' = dtmin')
	   END SELECT
	   CALL sortie
	   	  
c on reprend avec dt/2
	  ELSE I6
	   RETURN
	  ENDIF	I6  	  
	 ENDIF I5
	ENDIF I4

c lors de l hélium burning 0.1 < T9 < 0.2 au centre on limite precit < 10%
c avant la combustion du carbone 0.2 < T9 < 0.6, precit > 30%
	precit_eff=precit ; bid=EXP(bp(2,1))
	IF(bid > 1.d8)THEN
	 IF(bid < 2.d8)THEN
	  precit_eff=MIN(0.1d0,precit)
	 ELSEIF(bid < 5.d8)THEN
	  precit_eff=MAX(0.3d0,precit)
	 ELSEIF(bid < 1.d9)THEN
	  precit_eff=MAX(1.d0,precit)
	 ELSE
	  precit_eff=MAX(1.5d0,precit)
	 ENDIF
	ELSEIF(bid > 2.5d7)THEN
	 precit_eff=MAX(0.3d0,precit)
	ENDIF

	IF(est /= 0.d0)THEN
	 IF(diffusion)THEN
	  dtnew=0.95d0*dt*precit_eff/est
	 ELSE
	  dtnew=dt*(precit_eff/est)**(1.d0/(ordre+1))
	 ENDIF
	ELSE
	 dtnew=1.2d0*dt
	ENDIF
	
c on impose une variation maximale de dlntc~10% à Tc
	dtt=ABS(bp(2,1)-bp_t(2,1))/bp(2,1)
	IF(dtt > dlntc)dtnew=MIN(dtnew,dt*dlntc/dtt)

c pas temporel estimé
	SELECT CASE(Krot)
	CASE(3,4)
	 dtn=MAX(0.8d0*dt,MIN(dtnew,1.1d0*dt))
	CASE DEFAULT
	 dtn=MAX(0.8d0*dt,MIN(dtnew,1.2d0*dt))
	END SELECT

c les masses pour l erreur maximale
	mc_maxi=SQRT(ABS(mc_maxi))**3
	mc_maxi=mc_maxi/mstar

c écritures
	WRITE(usl_evol,111)kmax,est,nom_elem(ich_abon_max),dt,dtn
111	FORMAT(/,'EVOL couche ',i4,', var. max. :',
	1 es10.3,' pour ',a4,/,'dt=',es10.3,', dt optimal=',es10.3,
	2 ' masse/abondance/variation/couche:')	
	WRITE(usl_evol,114)nom_elem(1:MIN(10,nchim))
114	FORMAT(10(2x,a4,2x))
	WRITE(usl_evol,112)mc_maxi(1:MIN(10,nchim))
	WRITE(usl_evol,112)(ab_max(i)*nucleo(i),i=1,MIN(10,nchim))
	WRITE(usl_evol,112)estim(1:MIN(10,nchim))
112	FORMAT(10es8.1)	
	WRITE(usl_evol,113)kmaxi(1:MIN(10,nchim))
113	FORMAT(10(2x,i4,2x))

c pas temporel estimé
	dts=MIN(dtmax,dtn,agemax-age-dt) ; reprend=.FALSE.
c	PRINT*,'fin evol dt,dts,dtn,dtmax,agemax,age,agemax-age-dt',Krot
c	WRITE(*,2000)dt,dts,dtn,dtmax,agemax,age,agemax-age-dt
c	PAUSE'fin de evol'

	RETURN

	CONTAINS
	 INCLUDE 'base_rota.f'
	 
	END SUBROUTINE evol
