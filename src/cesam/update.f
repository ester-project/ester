      
c*******************************************************************

	SUBROUTINE update(next,dt,dts)

c	routine private du module mod_static	
c	passage au pas temporel suivant ou réinitialisation

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c entrées
c   next=.TRUE.	 ***--->***_t, dts-->dt : on decale d'un dt

c   next=.FALSE. ***_t--->*** reinitialisation

c solution, age=t+dt : 
c	bp,q,n,qt,knot,chim,mc,nc,mct,knotc,r2,m23,dts,
c	m_zc,r_zc,r_ov,lim,jlim,lconv,mstar
c solution, age=t :
c	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,r2_t,m23_t,dt,
c	m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,mstar_t
c	entre t+dt et t

c	dts: pas temporel estimé pour intégration entre age et age+1
c	dt : pas temporel utilisé entre age-1 et age
c	lconv=.true. : début de ZC

c sorties
c	bp: solution
c	dt: nouveau pas temporel
c	lnp,lnt,r2,l23,m23,q_t,qt_t,knot_t,nt,mu,dt: solution en t
c	m_t,mt_t,n_t,knot_t,xchims_t: idem en comp. chim.
c	nt_,m_,mt_,n_,knot_,xchims_: idem en t-1
c	m_zct : masses des limites ZR/ZC au temps t
c	lim_t : nombre de limites ZR/ZC au temps t
c	jlim_t : limites ZR/ZC au temps t
c	lconv_t=.true. : debut de ZC temps t

c Au debut de resout, il faut transferrer tds ---> tds_t qui sert à
c initialiser les itérations NR du modèle quasi-statique
c en cas d'échec des itérations il est inutile de
c remettre tds dans tds_t qui n'a pas été modifié
c sur le disque on écrit tds et non pas tds_t

c NOTATIONS (hélas incohérentes)
c	n_ch : nombre VARIABLE de points élément de mod_variables
c	nch : nombre FIXE de fonctions élément de mod_donnees
c	m_ch : ordre FIXE des splines élément de mod_donnees 
c	mch(n_ch) : abscisses VARIABLES élément de mod_variables

c----------------------------------------------------------------

	USE mod_atm, ONLY : bp_atm, bp_atm_t, x_atm, x_atm_t, xt_atm,
	1 xt_atm_t
	USE mod_donnees, ONLY : agemax, dtmax, dtmin, Krot, langue, 
	1 m_ch, m_rot, nchim, ne, nrot, nvth, ord_qs
	USE mod_kind
	USE mod_nuc, ONLY : age_deb, age_fin, dt_planet, l_planet
	USE mod_variables, ONLY: age, bp, bp_t, chim, chim_t, dim_ch,
	1 dim_qs, dim_rot, jlim, jlim_t, knotc, knotc_t, knot, knotr, knotr_t,
	2 knot_t, knot_tds, knot_tds_t, lconv, lconv_t, lim, lim_t,
	3 mc, mc_t, mct, mct_t, mrot, mrott, mrott_t, mrot_t, mstar,
	4 mstar_t, mw_tot, mw_tot_t, m_zc, m_zc_t, n_ch, n_ch_t, n_qs,
	5 n_qs_t, n_rot, n_rot_t, n_tds, n_tds_t, m23, m23_t, q, q_t,
	6 qt, qt_t, rota, rota_t, r2, r2_t, r_ov, r_ov_t, r_zc,
	7 r_zc_t, tds, vth, vth_t, wrot, wrot_t, x_tds, xt_tds, tds_t, x_tds_t,
	8 xt_tds_t
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(inout) :: dt, dts
	
	LOGICAL, INTENT(in) :: next
	
	LOGICAL :: reduc
	
c----------------------------------------------------------------------

2000	FORMAT(8es10.3)

c	PRINT*,'update', dt

c	en_masse = .true. variables lagrangiennes m23=m**23, r2=r**2
c	en_masse = .false. variables eulériennes m23=m, r2=r

	IF(next)THEN	!*** ---> ***_t
	
c on garde les vth (interpolation de variables thermodynamiques)
	 IF(ALLOCATED(vth_t))DEALLOCATE(vth_t)
	 ALLOCATE(vth_t(nvth,SIZE(vth,2))) ; vth_t=vth	

c les tableaux au temps age-dt sont alloués mais	 
c le nombre n_qs_t de couches du modele quasi-statique
c au temps age-dt, en est-il le même n_qs qu'au temps t ? 	  
	 IF(ALLOCATED(bp_t))DEALLOCATE(bp_t,m23_t,q_t,qt_t,r2_t)
	 ALLOCATE(bp_t(ne,dim_qs),m23_t(n_qs),q_t(n_qs),qt_t(knot),
	1 r2_t(n_qs))	  
	 IF(ALLOCATED(tds_t))DEALLOCATE(tds_t,x_tds_t,xt_tds_t)
	 ALLOCATE(tds_t(1,n_tds),x_tds_t(n_tds),xt_tds_t(knot_tds))

c le nombre n_ch de couches du vecteur de comp. chim. au temps t+dt
c a peut être changé	 
	 IF(ALLOCATED(chim_t))DEALLOCATE(chim_t,mc_t,mct_t,mrot_t,mrott_t,
	1 rota_t)	
	 ALLOCATE(chim_t(nchim,dim_ch),mc_t(n_ch),mct_t(knotc),mrot_t(n_rot),
	1 mrott_t(knotr),rota_t(nrot,SIZE(rota,2)))	
	
c transferts *** ---> ***_t
	 bp_t=bp ; chim_t=chim ; jlim_t=jlim ; knotc_t=knotc ; knot_t=knot
	 knot_tds_t=knot_tds
	 lconv_t=lconv ; lim_t=lim ; mc_t=mc ; mct_t=mct ; mstar_t=mstar
	 m_zc_t=m_zc ; n_ch_t=n_ch ; n_qs_t=n_qs ; n_tds_t=n_tds
	 q_t=q ; qt_t=qt ; mw_tot_t=mw_tot ; m23_t=m23
	 r_zc_t=r_zc ; r_ov_t=r_ov ; r2_t=r2
	 tds_t(1,1:n_tds)=tds(1,1:n_tds) ; x_tds_t=x_tds ; xt_tds_t=xt_tds
	 wrot_t=wrot
	 bp_atm_t=bp_atm ; x_atm_t=x_atm ; xt_atm_t=xt_atm	 
	 n_rot_t=n_rot ; mrot_t=mrot ; mrott_t=mrott
	 knotr_t=knotr ; rota_t=rota
	 
c nouveau dt	 	 
	 dt=MAX(dtmin,MIN(dtmax,dts,dt*1.2d0,agemax-age))

c contraintes sur le pas temporel pendant les chutes de planétoïdes
	 IF(l_planet)THEN
	  reduc=.FALSE.
	  IF(age < age_deb)THEN	  
	   IF(age+dt >= age_deb)THEN
	    dt=MIN(age_deb-age,dtmax) ; reduc=.TRUE.
	   ENDIF
	  ELSEIF(age >= age_deb-dtmin .AND. age <= age_fin)THEN
	   dt=MIN(dt,dt_planet,dtmax) ; reduc=.TRUE.
	  ELSEIF(age < age_fin)THEN
	   dt=MIN(age_fin-age,dtmax) ; reduc=.TRUE.	   
	  ENDIF
	  
c écriture	  
	  IF(reduc)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1001) ; WRITE(2,1001)
1001	    FORMAT('time step adjusted during the infall of planetoids',/)	
	   CASE DEFAULT	
	    WRITE(*,1) ; WRITE(2,1)
1	    FORMAT('chute de planétoïdes en cours, pas temporel ajusté',/)
	   END SELECT
	  ENDIF
	 ENDIF

c	 WRITE(*,2000)dtmin,dtmax,dts,age,agemax,dt ; WRITE(*,2000)q
c	 PAUSE'sortie update *** ---> ***_t'
 	 
c retour ***_t ---> ***
	ELSE
c	 WRITE(*,*)'UPDATE(2) dtmin,dtmax,dts,age,agemax,dt'
c	 WRITE(*,2000)dtmin,dtmax,dts,age,agemax,dt

	 dt=MAX(dtmin,MIN(dtmax,dt,agemax-age))	!réinitialisation
	 
c il est inutile de transferrer tds_t dans tds puisque seul
c tds_t est utilisé au temps t+dt et il est déjà connu	 
	 IF(n_qs_t /= n_qs)THEN
	  DEALLOCATE(bp,q,qt,r2,m23)
	  n_qs=n_qs_t ; knot=knot_t
	  n_tds=n_tds_t ; dim_qs=knot-ord_qs
	  ALLOCATE(bp(ne,dim_qs),m23(n_qs),q(n_qs),qt(knot),r2(n_qs))
	 ENDIF

c les dimensions des fichiers de comp. chim. ont-elles changé? 
	 IF(SIZE(chim_t,2) /= SIZE(chim,2))THEN
	  DEALLOCATE(chim,mc,mct)
	  ALLOCATE(chim(nchim,SIZE(chim_t,2)),mc(SIZE(mc_t)),mct(SIZE(mct_t)))
	 ENDIF
	 
c les dimensions des fichiers de rotation ont-elles changé?	
	 IF(SIZE(rota_t,2) /= SIZE(rota,2))THEN	
	  DEALLOCATE(mrot,mrott,rota)
	  ALLOCATE(mrot(SIZE(mrot_t)),mrott(SIZE(mrott_t)),
	1 rota(nrot,SIZE(rota_t,2)))
	 ENDIF
	 	 
	 bp=bp_t ; chim=chim_t ; jlim=jlim_t ; knot=knot_t
	 lconv=lconv_t ; lim=lim_t ; mstar=mstar_t
	 knotc=knotc_t ; dim_ch=knotc-m_ch ; n_ch=n_ch_t 
	 mc=mc_t ; mct=mct_t ; m_zc=m_zc_t ; n_qs=n_qs_t
	 q=q_t ; qt=qt_t ; mw_tot=mw_tot_t ; m23=m23_t
	 r_zc=r_zc_t ; r_ov=r_ov_t ; r2=r2_t ; wrot=wrot_t
	 bp_atm=bp_atm_t ; x_atm=x_atm_t ; xt_atm=xt_atm_t
	 n_rot=n_rot_t ; mrot=mrot_t ; mrott=mrott_t ; knotr=knotr_t
	 rota=rota_t ; dim_rot=knotr-m_rot
	 	 	 
	ENDIF
	
c	PRINT*,'sortie UPDATE next,dt',next,dt

	RETURN

	END SUBROUTINE update
