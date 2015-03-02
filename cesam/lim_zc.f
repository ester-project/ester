
c********************************************************************

	SUBROUTINE lim_zc(no_rep,new_nqs)

c	routine private du module mod_static

c	répartition des couches et détermination des limites ZR/ZC.
c	Le nombre est modifié si new_n /= n

c	avec pression turbulente 8 inconnues 
c	sans pression turbulente 7 inconnues, Ptot=Pgaz

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c	17 09 96 : correction bug dans déclaration de vt; signalé par
c	D. Cordier
c	25 08 97 : mise en place des variables eulériennes
c	19 11 99 : suppression de nh1, nhe1, nhe2, lamb
c	02 02 00 : introduction de no_supp
c	20 02 00 : écriture si non affinement de la loc. des lim. des ZR/ZC

c       entrées
c	no_rep : il n'y a pas eu de reprise
c	new_nqs : nouveau nombre de couches
c	mc,mct,nc,knotc,chim : pour interpolation
c	mstar: masse au temps t+dt, avec perte de masse

c       entrées/sorties
c	bp, q, qt, knot : solution en B-spline, var. pub. de mod_static
c	n_qs : nombre de couches, var. pub. de mod_quasi_static
c	fac : ctes de répar. par couche, var. priv. de mod_quasi_static
c	nl: nombre de lignes du Jacobien
c	dim : dimension de la base
c	r2, m23 : r**2, m**2/3 en lagrangien, r et m en eulérien

c       sorties
c	jlim(i) : plus proche indice de la i-ieme limite ZR/ZC
c	lim : nombre de limites ZR/ZC
c	lconv(i)=.TRUE. : début d'une ZC à la i-ieme limite ZR/ZC
c	derxx, lderxx : coefficients pour les splines
c	xx, xl : points de collocation et limites
c	ni(i) : nombre de limite en i
c	nl : nombre de lignes du Jacobien
c	m_zc : masses en Msol des limites ZR/ZC éventuellement overshootées
c	r_zc : rayons en Rsol des limites ZR/ZC sans overshoot
c	r_ov : rayons en Rsol des limites ZR/ZC éventuellement overshootées
c	fac(i) : poids de la couche [i, i+1]

c----------------------------------------------------------------

	USE mod_donnees, ONLY: alpha, aradia, clight, cpturb,
	1    ctel, ctem, ctep, cter, ctet, diffusion, dn_fixe, en_masse,
	2    g, grille_fixe, ihe4, Ipg, jpz, Krot, langue, loc_zc, ledoux,
	3    lsol, msol, m_ch, m_qs, m_rot, nchim, ne, nrot, ord_qs, ovshti,
	4    ovshts, pi, pnzc, pturb, rsol, r_qs, t_inf
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn, newspl, noedif, no_croiss, pause,
	1    zoning
	USE mod_variables, ONLY: bp, chim, dim_qs, hp_zc, hp_zc_conv,
	1    id_conv, if_conv, inter, jlim, lconv, knot, knotc, knotr, lim,
	2    mc, mct, mc_fixe, mrot, mrott, mstar, m_zc, m23, nc_fixe, n_ch,
	3    n_qs, n_rot, pr_zc, pr_zc_conv, q, qt, rota, rstar, r2, r_ov, r_zc,
	4    r_zc_conv, wrot, tot_conv, tot_rad

	IMPLICIT NONE

	INTEGER, INTENT(in) :: new_nqs
	LOGICAL, INTENT(inout) :: no_rep

	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: new_b, new_esp
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: comp, dcomp,
	1    depsx, dfdq, dfrot, d_grad, f, frot, rap
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: new_q, new_qt, esp, new
	REAL (kind=dp), DIMENSION(5) :: epsilon

	REAL (kind=dp) :: pt, dlpp, p, t, r, l, m, m13, w, rn, grp, gr,
	1    j_new, j1, j2, j3, kip, drot, ro, dkapt, kap, gradad, drop,
	2    dkapp, hp, drox, u, dup, dut, dux, grad, dgradpt, dgradp, dgradt,
	3    dgradx, dgradm, dgradl, dgradr, dgradlpp, gam, dgampt, dgamp,
	4    dgamt, dgamx, dgamm, dgaml, dgamr, dgamlpp, depsp,
	5    depst, dkapx, delta, deltap, deltat, deltax, cp, dcpp,
	6    dcpt, dcpx, dgradadp, dgradadt, dgradadx, dhppt, dhpp, dhpt,
	7    dhpx, dhpr, dhpm, gradrad, alfa, beta, gamma1

	INTEGER, SAVE :: dlt_n=60, lc=1, ncouches_zc=5  
	INTEGER :: new_knot, i, j, lq, compt

	LOGICAL, SAVE :: der, ovsht, init=.TRUE.
	LOGICAL :: no_supp, radiatif

c-------------------------------------------------------------------

 2000	FORMAT(8es10.3)

	IF(init)THEN
	   init=.FALSE.
	   der=cpturb < 0.d0	!der=T. on tient compte de dln Pgaz/dln Ptot
	   SELECT CASE(langue)
	   CASE('english')
	      WRITE(2,1020)loc_zc ; WRITE(*,1020)loc_zc
 1020	      FORMAT('distance of the limit RZ/CZ to the closest node :',
	1	   es10.3,/)	 
	   CASE DEFAULT
	      WRITE(2,20)loc_zc ; WRITE(*,20)loc_zc
 20	      FORMAT('écart des limites ZR/ZC au plus proche noeud :',es10.3,/)
	   END SELECT
	   ovsht=ovshti > 0.d0 .OR. ovshts > 0.d0
	   IF(ovsht)THEN
	      IF(jpz)THEN
		 IF(ovshti > 0.d0)THEN
		    SELECT CASE(langue)
	            CASE('english')
		       WRITE(*,1009)INT(ovshti*100) ; WRITE(2,1009)INT(ovshti*100)
 1009		       FORMAT('overshooting JPZ beneath the CZ : ',i3,'% Hp/Kip',/)
		    CASE DEFAULT	   
		       WRITE(*,9)INT(ovshti*100) ; WRITE(2,9)INT(ovshti*100)
 9		       FORMAT('overshooting inférieur de JPZ des ZC :',i3,'% Hp/Kip')
		    END SELECT
		 ENDIF	   
		 IF(ovshts > 0.d0)THEN
		    SELECT CASE(langue)
		    CASE('english')
		       WRITE(*,1006)INT(ovshts*100) ; WRITE(2,1006)INT(ovshts*100)
 1006		       FORMAT('overshooting JPZ above the CZ : :',i3,'% R noyau',/)
		    CASE DEFAULT	   
		       WRITE(*,6)INT(ovshts*100) ; WRITE(2,6)INT(ovshts*100)
 6		       FORMAT('overshooting supérieur de JPZ des ZC :',i3,
	1		    '% R noyau',/)
		    END SELECT
		 ENDIF
	      ELSE
		 IF(ovshti > 0.d0)THEN
		    SELECT CASE(langue)
		    CASE('english')
		       WRITE(*,1007)INT(ovshts*100) ; WRITE(2,1007)INT(ovshts*100)
 1007		       FORMAT('overshooting beneath the CZ :',i3,'%  Hp')	    
		    CASE DEFAULT
		       WRITE(*,7)INT(ovshts*100) ; WRITE(2,7)INT(ovshts*100)
 7		       FORMAT('overshooting inférieur des ZC :',i3,'%  Hp',/)
		    END SELECT   
		 ENDIF
		 IF(ovshts > 0.d0)THEN
		    SELECT CASE(langue)
		    CASE('english')
		       WRITE(*,1008)INT(ovshts*100) ; WRITE(2,1008)INT(ovshts*100) 
 1008		       FORMAT('overshooting above the CZ :',i3,'%  Hp',/)	    
		    CASE DEFAULT	   
		       WRITE(*,8)INT(ovshts*100) ; WRITE(2,8)INT(ovshts*100)
 8		       FORMAT('overshooting supérieur des ZC :',i3,'%  Hp',/)
		    END SELECT 
		 ENDIF
	      ENDIF
	   ELSE
	      SELECT CASE(langue)
	      CASE('english')
		 WRITE(*,1005) ; WRITE(2,1005)
 1005		 FORMAT('Model without overshooting')	 
	      CASE DEFAULT	 
		 WRITE(*,5) ; WRITE(2,5)
 5		 FORMAT('Modèle sans overshooting')
	      END SELECT
	      
	      IF(Krot == 3)ALLOCATE(dfrot(nrot),frot(nrot))
	   ENDIF
	   
c       allocations initiales
	   
	   ALLOCATE(fac(n_qs-1),d_grad(n_qs),f(ne),dfdq(ne),comp(nchim),
	1	dcomp(nchim),depsx(nchim),rap(2*pnzc))
	ENDIF			!initialisation
	
c	on redéfinit la répartition avant la recherche des limites ZR/ZC
c	si le modèle n'est pas estimé "bien répartit" ou si le nb.
c	de couches du modèle quasi. stat. doit changer
	
	IF(.NOT.no_rep .OR. new_nqs /= n_qs) THEN
	   
c       les vecteurs fac et d_grad n'ont plus la même dimension
	   
	   DEALLOCATE(fac,d_grad) ; no_rep=.TRUE.

c       on met dans esp la répartition précédente

	   ALLOCATE(esp(n_qs))
	   DO i=1,n_qs
	      CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),lq,f,dfdq)
	      IF(no_croiss)PRINT*,'Problème localisé en 1 dans lim_zc'
	      esp(i)=ctep*f(1)+ctet*f(2)+cter*f(3)+ctel*f(4)+ctem*f(5)
	      IF(i > 1)THEN
		 DO WHILE(esp(i) < esp(i-1))
		    esp(i)=esp(i)+0.01d0
		 ENDDO
	      ENDIF
c       WRITE(*,2000)esp(1:ne,i)
	   ENDDO
c       CALL pause('nouvel esp')

c       on dispose les couches pour assurer une répar. approx.
c       uniforme de la fonction de répartition sur new_nqs couches 

	   ALLOCATE(new_q(new_nqs))
	   CALL zoning(n_qs,q,esp,new_nqs,new_q) !choix des nouveaux q

c       esp est désormais inutile

	   DEALLOCATE(esp)

c       PRINT*,'anciens q',n_qs,new_nqs ; WRITE(*,2000)q
c       PRINT*,'esp' ; WRITE(*,2000)esp ; PRINT*,'nouveaux q'
c       WRITE(*,2000)new_q ; CALL pause('après zoning')

c       dans new_b on se place aux new_nqs points new_q
c       on calcule la nouvelle fonction de répartition

	   ALLOCATE(new_b(ne,new_nqs),new_esp(1,new_nqs))
	   DO i=1,new_nqs
	      CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,new_q(i),lq,f,
	1	   dfdq)	  
	      IF(no_croiss)PRINT*,'Problème localisé en 2 dans lim_zc'
	      new_b(1:ne,i)=f(1:ne)
c       WRITE(*,2000)new_b(1:ne,i)
	      new_esp(1,i)=ctep*f(1)+ctet*f(2)+cter*f(3)+ctel*f(4)+ctem*f(5)
c       WRITE(*,2000)q(i),(f(j),j=1,ne),new_esp(1,i)
	      IF(i > 1)THEN
		 DO WHILE(new_esp(1,i) < new_esp(1,i-1))
		    new_esp(1,i)=new_esp(1,i)+0.01d0
		 ENDDO
	      ENDIF
	   ENDDO
c       CALL pause('new_b')

c       spline d'ordre 2 sur les nouveaux q sur new_nqs points

	   new_q=(/ (i, i=1,new_nqs) /)
	   ALLOCATE(new_qt(new_nqs+ord_qs))
	   CALL bsp1dn(1,new_esp,new_q,new_qt,new_nqs,2,new_knot,.FALSE.,
	1	new_q(1),lq,f,dfdq)  
	   IF(no_croiss)THEN
	      PRINT*,'Arrêt 1 dans lim_zc' ; STOP
	   ENDIF

c       dérivée de la fonction de répartition dans new_b(6,:)
	   
	   DO i=1,new_nqs 
	      CALL bsp1dn(1,new_esp,new_q,new_qt,new_nqs,2,new_knot,.TRUE.,
	1	   new_q(i),lq,f,dfdq)
	      IF(no_croiss)PRINT*,'Problème localisé en 3 dans lim_zc'	
	      new_b(6,i)=dfdq(1)
	   ENDDO
c       WRITE(*,2000)new_b(6,:) ; CALL pause('new_b(6,:)')

c       on spline new_b

	   new_q=(/ (i, i=1,new_nqs) /)
	   CALL bsp1dn(ne,new_b,new_q,new_qt,new_nqs,ord_qs,new_knot,.FALSE.,
	1	new_q(1),lq,f,dfdq)
	   IF(no_croiss)THEN
	      PRINT*,'Arrêt 2 dans lim_zc' ; STOP
	   ENDIF

c       DO i=1,n_qs
c       CALL bsp1dn(ne,new_b,new_q,new_qt,new_nqs,ord_qs,new_knot,.TRUE.,
c	1  new_q(i),lq,f,dfdq)
c       WRITE(*,2000)f(1:ne)
c       ENDDO
c       CALL pause('new_b')

c       dimension de la base splines pour l'intégration équi.quasi.stat. 

	   n_qs=new_nqs ; dim_qs=(n_qs-1)*m_qs+r_qs

c       vecteur nodal pour l'intégration sur n_qs points

	   DEALLOCATE(q,qt)
	   ALLOCATE(q(n_qs),qt(dim_qs+ord_qs))
	   q=(/ (i, i=1,n_qs) /)   
	   CALL noedif(q,n_qs,m_qs,r_qs,qt,knot)
	   IF(no_croiss)THEN
	      PRINT*,'Arrêt 3 dans lim_zc' ; STOP
	   ENDIF
	   
c       PRINT*,ord_qs,knot,dim_qs
c       WRITE(*,2000)qt(1:knot) ; CALL pause('noedif')

c       on place la solution initiale dans la base de q, bp_t ==> bp

	   DEALLOCATE(bp)
	   ALLOCATE(bp(ne,dim_qs))
	   CALL newspl(ne,new_q,new_qt,new_knot,ord_qs,qt,knot,ord_qs,
	1	new_b,bp)
c       CALL pause('newspl')

c       new_b, new_q, new_qt, new_esp sont désormais inutiles
c       r2 et m23 sont à modifier

	   DEALLOCATE(new_b,new_q,new_qt,new_esp,r2,m23)

c       adaptation de fac, d_grad, r2 et m23 à la nouvelle répartition

	   ALLOCATE(fac(n_qs-1),d_grad(n_qs),r2(n_qs),m23(n_qs))

c       tabulation de r2 et m23 qui accélèrent la recherche dans
c       le ssp. inter

	   DO i=1,n_qs 
c       CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),lc,f,dfdq)
	      r2(i)=bp(3,m_qs*(i-1)+1) ; m23(i)=bp(5,m_qs*(i-1)+1)
	   ENDDO
	   
c       test pour changement de grille fixe pour la composition chimique
c       il y a changement de grille fixe quand si
c       new_nqs diffère de nc_fixe de plus de dn_fixe (5% environ)

	   IF(grille_fixe)THEN
	      IF(ABS(REAL(new_nqs-nc_fixe)/new_nqs) > dn_fixe)THEN
		 WRITE(*,30)nc_fixe,n_qs ; WRITE(2,30)nc_fixe,n_qs	  
 30		 FORMAT(/,'changement de grille fixe: nc_fixe=',i4,' ==> ',i4,/)
		 DEALLOCATE(mc_fixe) ; ALLOCATE(mc_fixe(n_qs))	  
		 nc_fixe=n_qs ; mc_fixe=m23
		 IF(.NOT.en_masse)mc_fixe=mc_fixe**(2.d0/3.d0)
	      ENDIF
	   ENDIF
	ENDIF
	bp(3:5,1)=0.d0		!au centre R, L, M =0

c------------------recherche des limites ZR/ZC----------------------------

c	recherche les limites zones convectives / zones radiatives
c	en fonction de q, détermination des facteurs de répartition
c	overshooting inférieur et supérieur : extension de la zone mélangée
c	ou overshooting inférieur PC de JPZ avec saut du gradient

c	modèle totalement convectif: lim=1, jlim(1)=n, lconv(1)=.FALSE.
c	modèle totalement radiatif: lim=0, jlim(i)=-100, lconv(i)=.FALSE.

c	CALL pause('début recherche limites ZR/ZC')

	SELECT CASE(langue)
	CASE('english')
	   WRITE(*,1010)
 1010	   FORMAT(/,'----- Search of the limits RZ/CZ (begin) --------')		
	CASE DEFAULT
	   WRITE(*,10)
 10	   FORMAT(/,'------ Recherche des limites ZR/ZC (début) --------')
	END SELECT

c	mise à zero des r_zc et r_ov par securite (remarque de D. Cordier)

	r_ov=-100.d0 ;  r_zc=-100.d0 ; m_zc=-100.d0

	tot_conv=.FALSE. ; tot_rad =.FALSE. ; ovsht=MAX(ovshts,ovshti) > 0.d0

c	calcul de grad-grad

c	PRINT*,n_ch,m_ch,knotc,n_qs,nchim ; WRITE(*,2000)mc
c	WRITE(*,2000)mct ; PRINT*,mstar ; CALL pause('début lim_zc')

	DO i=1,n_qs		!calcul de grad-grad
	   CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),lq,f,dfdq)
	   IF(no_croiss)PRINT*,'Problème localisé en 4 dans lim_zc'	 
c       WRITE(*,2000)f
	   pt=EXP(f(1))
	   dlpp=1.d0
	   IF(pturb)THEN  	!avec pression turbulente
	      p=EXP(f(Ipg))
	      IF(der)dlpp=dfdq(Ipg)/dfdq(1) !dlpp=dln Pgaz/dln Ptot
	   ELSE			!sans pression turbulente
	      p=pt
	   ENDIF
	   t=EXP(f(2))   
	   IF(en_masse)THEN
	      r=SQRT(ABS(f(3))) ; l=SQRT(ABS(f(4)))**3 ; m=SQRT(ABS(f(5)))**3
	   ELSE
	      r=ABS(f(3)) ; l=f(4) ; m=ABS(f(5)) ; f(5)=m**(2.d0/3.d0)
	   ENDIF
	   m13=MAX(m**(1.d0/3.d0),1.d-5) !en 0 décalage pour calcul dérivée
	   CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1	knotc,.TRUE.,MAX(mc(1),MIN(f(5),mc(n_ch))),lc,comp,dcomp)
	   IF(no_croiss)PRINT*,'Problème localisé en 5 dans lim_zc'	
c       WRITE(*,2000)p,t,r,l,m,comp(1),q(i)    !,mc(n_ch),f(5)
c       PRINT*,mc(1),mc(n_ch),f(5),comp(1) ; CALL pause('boucle lim_zc')
	   DO j=1,nchim		!dX/dm
	      dcomp(j)=dcomp(j)*2.d0/3.d0/m13
	   ENDDO   

c       vitesse angulaire

	   SELECT CASE(Krot)
	   CASE(0,1,2)
	      w=wrot
	   CASE(3)
	      CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,
	1	   MAX(mrot(1),MIN(f(5),mrot(n_rot))),lc,frot,dfrot)
	      IF(no_croiss)PRINT*,'Problème localisé en 6 dans lim_zc'	
	      w=frot(1)
	   END SELECT 
	   
	   d_grad(i)=dgrad(pt,p,t,dlpp,comp,m,l,r,dcomp,w)
c       WRITE(*,2000)p,t,m,l,r,comp(1),comp(ihe4),d_grad(i)	 
	ENDDO			!i
c	WRITE(*,2000)d_grad ; CALL pause('d_grad')

	rn=r			!r à l'extérieur

	grp=d_grad(1)		!ZC avec zero exact en 1
	IF(grp == 0.d0 .AND. d_grad(2) >= 0.d0)grp=-1
	lim=0			!initialisation 
	B1: DO j=2,n_qs		!recherche des limites ZR/ZC
	   gr=d_grad(j)
	   IF(grp /= 0.d0 .AND. gr*grp <= 0.d0)THEN
	      lim=lim+1
	      IF(lim > 2*pnzc-1)THEN
		 PRINT*,'nombre de limites ZR/ZC supérieur à',pnzc*2-1
		 PRINT*,'limites ZR/ZC imprécises, floues ?'
		 PRINT*,'difficulté d''initialisation ? erreur ? '
		 WRITE(*,"('poursuite en limitant à',i2,' le nb. de ZC')")pnzc   
		 lim=lim-1  
		 EXIT B1    
c       tot_rad=.TRUE. ; lim=0
c       DO i=1,pnzc
c       jlim(i)=-100 ; lconv(i)=.FALSE.
c       ENDDO
c       EXIT B1
	      ENDIF
	      lconv(lim)=grp < 0.d0 !.TRUE. passage ZR-->ZC
	      rap(lim)=ABS(grp/gr)
	      rap(lim)=rap(lim)/(1.d0+rap(lim))
	      jlim(lim)=j	!limite ZR/ZC entre j et j-1
	      IF(lconv(lim))THEN
		 SELECT CASE(langue)
	         CASE('english')
		    WRITE(*,1003)INT(rap(lim)*100.d0),j-1,j,lim
 1003		    FORMAT('a CZ begins at',i3,'% between the shells #',
	1		 i5,' and #',i5,', limit RZ/CZ #',i3)
		 CASE DEFAULT
		    WRITE(*,3)INT(rap(lim)*100.d0),j-1,j,lim
 3		    FORMAT('début ZC localisé à',i3,'% entre les couches #',
	1		 i5,' et #',i5,', limite ZR/ZC #',i3)
		 END SELECT
	      ELSE
		 SELECT CASE(langue)
	         CASE('english')
		    WRITE(*,1004)INT(rap(lim)*100.),j-1,j,lim
 1004		    FORMAT('a CZ ends at ',i3,'% between the shells #',
	1		 i5,' and #',i5,', limit RZ/CZ #',i3)
		 CASE DEFAULT	  
		    WRITE(*,4)INT(rap(lim)*100.),j-1,j,lim
 4		    FORMAT('fin de ZC localisée à',i3,'% entre les couches #',
	1		 i5,' et #',i5,', limite ZR/ZC #',i3)
		 END SELECT
	      ENDIF
	   ENDIF
	   grp=gr
	ENDDO B1	!j

c	PRINT*,'sortie du tri, lim',lconv(1:lim),lim
c	WRITE(*,2000)rap(1:lim) ; CALL pause('tri')

c	indices des points de grille
c	IF (ALLOCATED(new)) DEALLOCATE(new)
	ALLOCATE(new(n_qs)) ; new=q

c	il y a une limite ZR/ZC entre jlim(i)-1 et jlim(i), i=1,lim
c	on détermine l'incide new(jlim(i)) à la limite ZR/ZC
c	on fixe la limite soit à droite soit à gauche
c	jlim(i) : indice le plus proche de la limite

	DO i=1,lim		!j_new: indice de la lim. ZR/ZC
	   j_new=new(jlim(i)-1)+rap(i)*(new(jlim(i))-new(jlim(i)-1))
c       PRINT*,i,j_new,new(jlim(i)-1),rap(i),new(jlim(i))
	   jlim(i)=NINT(j_new)    
	   new(jlim(i))=j_new	!nouvel indice pour les limites ZR/ZC
c       PRINT*,i,j_new,jlim(i)
	ENDDO

c	PRINT*,'indice pour les limites ZR/ZC' ; PRINT*,lim
c	PRINT*,(new(jlim(i)),i=1,lim) ; PRINT*,(jlim(i),i=1,lim)    

c	suppression des limites trop proches
c	no_supp=.FALSE. : suppression de limites séparées par moins de
c	2 couches, pas d'affinement de la position des limites ZR/ZC

	no_supp=.TRUE.
	i=1
	DO WHILE (i <= lim)
	   IF(jlim(i) <= 3)THEN	!limites sur premieres couches
	      PRINT*,'suppression d''une limite ZR/ZC trop centrale #',i
	      lim=lim-1		!supprimees
	      DO j=1,lim
		 jlim(j)=jlim(j+1) ; new(jlim(j))=new(jlim(j+1))
		 lconv(j)=lconv(j+1)
	      ENDDO
	   ELSE
	      i=i+1
	   ENDIF
	ENDDO

c	PRINT*,'limites sur premières couches' ; PRINT*,lim
c	PRINT*,(new(jlim(i)),i=1,lim) ; PRINT*,(jlim(i),i=1,lim)    

	i=2
	DO WHILE (i <= lim)
	   IF(jlim(i)-jlim(i-1) <= ncouches_zc)THEN !limites trop proches
	      PRINT*,'limites trop proches, suppression des limites #',i-1,i
	      no_supp=no_supp .AND. jlim(i)-jlim(i-1) >= 2    
	      lim=lim-2		!suppression de i-1 et i
	      DO j=i-1,lim
		 jlim(j)=jlim(j+2) ; new(jlim(j))=new(jlim(j+2))
		 lconv(j)=lconv(j+2)
	      ENDDO		!no_supp=.FALSE. : pas d'affinement des limites
	   ELSE
	      i=i+1
	   ENDIF		!pour la suppression de limites trop rapprochées
	ENDDO

c	PRINT*,'limites trop rapprochées'
c	PRINT*,lim
c	PRINT*,(new(jlim(i)),i=1,lim)
c	PRINT*,(jlim(i),i=1,lim)    

	IF(.NOT.diffusion)THEN  !diffusion no suppression de lim externes 
	   i=1
	   DO WHILE (i <= lim)
	      IF(jlim(i) >= n_qs-ncouches_zc)THEN
		 PRINT*,'suppression d''une limite ZR/ZC trop externe #',lim
		 new(jlim(i))=jlim(i) ; lim=lim-1
	      ELSE
		 i=i+1
	      ENDIF
	   ENDDO
	ENDIF

c	PRINT*,'limite ZR/ZC trop externe' ; PRINT*,lim
c	PRINT*,(new(jlim(i)),i=1,lim) ; PRINT*,(jlim(i),i=1,lim)

c	nature complètement convectif/radiatif du modèle 

	IF(lim == 0 .AND. jlim(1) /= -100)THEN
	   IF(d_grad(4) >= 0.d0)THEN
	      PRINT*,'modèle complètement convectif'
	      tot_conv=.TRUE. ; lim=1 ; jlim(1)=n_qs
	      lconv(1)=.FALSE. ; m_zc(1)=mstar ; r_zc(1)=rn

c******************************novo ***********************************
	      CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,REAL(n_qs,dp),
	1	   lq,f,dfdq)
c       j3=REAL(jlim(i),dp)
c       CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,j3,lq,f,dfdq)	
	      IF(no_croiss)PRINT*,'Problème localisé en 23 dans lim_zc'
c----------------NEW --------- NEW --------- NEW ------------------------	
	      pt=EXP(f(1)) ; dlpp=1.d0
	      IF(pturb)THEN	!avec pression turbulente
		 p=EXP(f(Ipg))
		 IF(der)dlpp=dfdq(Ipg)/dfdq(1) !dlpp=dln Pgaz/dln Ptot
	      ELSE		!sans pression turbulente
		 p=pt
	      ENDIF
	      t=EXP(f(2))
	      IF(en_masse)THEN
		 l=SQRT(ABS(f(4)))**3
	      ELSE
		 l=f(4)
		 f(5)=m_zc(1)**(2.d0/3.d0)
	      ENDIF  
	      m13=MAX(m_zc(1)**(1.d0/3.d0),1.d-5) !en 0 decalage
	      CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1	   knotc,.TRUE.,MAX(mc(1),MIN(f(5),mc(n_ch))),lc,comp,dcomp)
	      IF(no_croiss)PRINT*,'Problème localisé en 25 dans lim_zc'
	      DO j=1,nchim	!dX/dm
		 dcomp(j)=dcomp(j)*2.d0/3.d0/m13
	      ENDDO
	      IF(no_croiss)PRINT*,'Problème localisé en 25 dans lim_zc'
	      DO j=1,nchim	!dX/dm
		 dcomp(j)=dcomp(j)*2.d0/3.d0/m13
	      ENDDO
c       PRINT*,'m_zc,r_zc,new(jlim(i))',lconv(i),i
c       WRITE(*,2000)m_zc(i),r_zc(i),new(jlim(i))
c       WRITE(*,2000)p,t,r_zc(i),l,m_zc(i),comp(1)

c       vitesse angulaire

	      SELECT CASE(Krot)
	      CASE(0,1,2)
		 w=wrot
	      CASE(3)
		 CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,
	1	      MAX(mrot(1),MIN(f(5),mrot(n_rot))),lc,frot,dfrot)
		 IF(no_croiss)PRINT*,'Problème localisé en 26 dans lim_zc'	
		 w=frot(1)
	      END SELECT
	      
	      CALL thermo(pt,p,t,m_zc(1),l,r_zc(1),dlpp,comp,dcomp,
	1	   ro,drop,drot,drox,u,dup,dut,dux,
	2	   grad,dgradpt,dgradp,dgradt,dgradx,dgradm,dgradl,dgradr,dgradlpp,
	3	   gam,dgampt,dgamp,dgamt,dgamx,dgamm,dgaml,dgamr,dgamlpp,
	4	   epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
	5	   delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	6	   gradad,dgradadp,dgradadt,dgradadx,
	7	   hp,dhppt,dhpp,dhpt,dhpx,dhpr,dhpm,
	8	   gradrad,alfa,beta,gamma1,radiatif)
	      
	      pr_zc(1)=p ; hp_zc(1)=hp
c******************end of new ******************************
	   ELSE
	      PRINT*,'modèle complètement radiatif'
	      tot_rad=.TRUE. ; lim=0 ; jlim=-100 ; lconv=.FALSE. ; m_zc=-100.d0
	   ENDIF
	   
	ELSE			!il y a des limites ZR/ZC

c--------on affine, par dichotomie, la position des limites retenues

c       PRINT*,no_supp ; WRITE(*,2000)new(jlim(1:lim))

	   IF(no_supp)THEN
	      DO i=1,lim
		 j1=new(jlim(i))
		 CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,j1,lq,f,dfdq)
		 IF(no_croiss)PRINT*,'Problème localisé en 7 dans lim_zc'	   
c       WRITE(*,2000)f
		 pt=EXP(f(1)) ; dlpp=1.d0   
		 IF(pturb)THEN	!avec pression turbulente 
		    p=EXP(f(Ipg))
		    IF(der)dlpp=dfdq(Ipg)/dfdq(1) !dlpp=dln Pgaz/dln Ptot
		 ELSE		!sans pression turbulente
		    p=pt
		 ENDIF
		 t=EXP(f(2))
		 IF(en_masse)THEN
		    r=SQRT(ABS(f(3))) ; l=SQRT(ABS(f(4)))**3 ; m=SQRT(ABS(f(5)))**3
		 ELSE
		    r=ABS(f(3)) ; l=f(4) ; m=ABS(f(5)) ; f(5)=m**(2.d0/3.d0)  
		 ENDIF
		 m13=MAX(m**(1.d0/3.d0),1.d-5) !en 0 decalage pour calcul dérivée
		 CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1	      knotc,.TRUE.,MAX(mc(1),MIN(f(5),mc(n_ch))),lc,comp,dcomp)
		 IF(no_croiss)PRINT*,'Problème localisé en 8 dans lim_zc'	
		 DO j=1,nchim   !dX/dm
		    dcomp(j)=dcomp(j)*2.d0/3.d0/m13
		 ENDDO

c       vitesse angulaire

		 SELECT CASE(Krot)
	         CASE(0,1,2)
		    w=wrot
		 CASE(3)
		    CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,
	1		 MAX(mrot(1),MIN(f(5),mrot(n_rot))),lc,frot,dfrot)
		    IF(no_croiss)PRINT*,'Problème localisé en 9 dans lim_zc'	
		    w=frot(1)
		 END SELECT
		 
		 d_grad(1)=dgrad(pt,p,t,dlpp,comp,m,l,r,dcomp,w)
		 
		 j2=jlim(i)
		 CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,j2,lq,f,dfdq)
		 IF(no_croiss)PRINT*,'Problème localisé en 10 dans lim_zc'
c       WRITE(*,2000)f
		 pt=EXP(f(1))
		 dlpp=1.d0 
		 IF(pturb)THEN	!avec pression turbulente
		    p=EXP(f(Ipg))
		    IF(der) dlpp=dfdq(Ipg)/dfdq(1) !dlpp=dln Pgaz/dln Ptot
		 ELSE		!sans pression turbulente
		    p=pt
		 ENDIF
		 t=EXP(f(2))
		 IF(en_masse)THEN
		    r=SQRT(ABS(f(3))) ; l=SQRT(ABS(f(4)))**3 ; m=SQRT(ABS(f(5)))**3
		 ELSE
		    r=ABS(f(3)) ; l=f(4) ; m=ABS(f(5)) ; f(5)=m**(2.d0/3.d0)  
		 ENDIF
		 m13=MAX(m**(1.d0/3.d0),1.d-5) !en 0 decalage pour calcul der.
		 CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1	      knotc,.TRUE.,MAX(mc(1),MIN(f(5),mc(n_ch))),lc,comp,dcomp)
		 IF(no_croiss)PRINT*,'Problème localisé en 11 dans lim_zc'	
		 DO j=1,nchim	!dX/dm
		    dcomp(j)=dcomp(j)*2.d0/3.d0/m13
		 ENDDO
		 
c       vitesse angulaire
		 
		 SELECT CASE(Krot)
	         CASE(0,1,2)
		    w=wrot
		 CASE(3)
		    CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,
	1		 MAX(mrot(1),MIN(f(5),mrot(n_rot))),lc,frot,dfrot)
		    IF(no_croiss)PRINT*,'Problème localisé en 12 dans lim_zc'	
		    w=frot(1)
		 END SELECT
		 
		 d_grad(2)=dgrad(pt,p,t,dlpp,comp,m,l,r,dcomp,w)    
		 
		 IF(d_grad(1)*d_grad(2) > 0.d0)THEN
		    j2=jlim(i)-1
		    CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,j2,lq,f,dfdq)
		    IF(no_croiss)PRINT*,'Problème localisé en 13 dans lim_zc'	    
c       WRITE(*,2000)f
		    pt=EXP(f(1)) ; dlpp=1.d0
		    IF(pturb)THEN !avec pression turbulente
		       p=EXP(f(Ipg))
		       IF(der)dlpp=dfdq(Ipg)/dfdq(1) !dlpp=dln Pgaz/dln Ptot
		    ELSE	!sans pression turbulente
		       p=pt
		    ENDIF
		    t=EXP(f(2))
		    IF(en_masse)THEN
		       r=SQRT(ABS(f(3))) ; l=SQRT(ABS(f(4)))**3 ; m=SQRT(ABS(f(5)))**3
		    ELSE
		       r=ABS(f(3)) ; l=f(4) ; m=ABS(f(5)) ; f(5)=m**(2.d0/3.d0)
		    ENDIF
		    m13=MAX(m**(1.d0/3.d0),1.d-5) !en 0 decalage pour calcul dérivée
		    CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1		 knotc,.TRUE.,MAX(mc(1),MIN(f(5),mc(n_ch))),lc,comp,dcomp)
		    IF(no_croiss)PRINT*,'Problème localisé en 14 dans lim_zc'	
		    DO j=1,nchim !dX/dm
		       dcomp(j)=dcomp(j)*2.d0/3.d0/m13
		    ENDDO
		    
c       vitesse angulaire
		    
		    SELECT CASE(Krot)
		    CASE(0,1,2)
		       w=wrot
		    CASE(3)
		       CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,
	1		    MAX(mrot(1),MIN(f(5),mrot(n_rot))),lc,frot,dfrot)
		       IF(no_croiss)PRINT*,'Problème localisé en 15 dans lim_zc'	
		       w=frot(1)
		    END SELECT
		    
		    d_grad(2)=dgrad(pt,p,t,dlpp,comp,m,l,r,dcomp,w)
		    
		    IF(d_grad(1)*d_grad(2) > 0.d0)THEN  
		       j2=jlim(i)+1
		       CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,j2,lq,f,dfdq)
		       IF(no_croiss)PRINT*,'Problème localisé en 16 dans lim_zc'	     
c       WRITE(*,2000)f
		       pt=EXP(f(1)) ; dlpp=1.d0 
		       IF(pturb)THEN !avec pression turbulente
			  p=EXP(f(Ipg))
			  IF(der) dlpp=dfdq(Ipg)/dfdq(1) !dlpp=dln Pgaz/dln Ptot
		       ELSE	!sans pression turbulente
			  p=pt
		       ENDIF
		       t=EXP(f(2))
		       IF(en_masse)THEN
			  r=SQRT(ABS(f(3))) ; l=SQRT(ABS(f(4)))**3
			  m=SQRT(ABS(f(5)))**3
		       ELSE
			  r=ABS(f(3)) ; l=f(4) ; m=ABS(f(5)) ; f(5)=m**(2.d0/3.d0)
		       ENDIF
		       m13=MAX(m**(1.d0/3.d0),1.d-5) !en 0 décalage pour calcul dér.
		       CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1		    knotc,.TRUE.,MAX(mc(1),MIN(f(5),mc(n_ch))),lc,comp,dcomp)
		       IF(no_croiss)PRINT*,'Problème localisé en 17 dans lim_zc'	
		       DO j=1,nchim !dX/dm
			  dcomp(j)=dcomp(j)*2.d0/3.d0/m13
		       ENDDO
		       
c       vitesse angulaire
		       
		       SELECT CASE(Krot)
		       CASE(0,1,2)
			  w=wrot
		       CASE(3)
			  CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,
	1		       MAX(mrot(1),MIN(f(5),mrot(n_rot))),lc,frot,dfrot)
			  IF(no_croiss)PRINT*,'Problème localisé en 18 dans lim_zc'	
			  w=frot(1)
		       END SELECT
		       
		       d_grad(2)=dgrad(pt,p,t,dlpp,comp,m,l,r,dcomp,w)
		       
		       IF(d_grad(1)*d_grad(2) > 0.d0)THEN
			  PRINT*,'lim_zc : erreur dans la dichotomie, limite ZR/ZC',i
			  PRINT*,'j1, jlim(i)',j1,jlim(i) ; PRINT*,'on arrète'
c       CALL pause('on arrète')
			  STOP
		       ENDIF	!jlim(i)+1
		    ENDIF	!jlim(i)-1
		 ENDIF		!on a encadré la limite entre j1 et j2
		 compt=0
		 j3=(j1+j2)/2.d0
		 DO WHILE(ABS(j1-j2) > loc_zc)
		    CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,j3,lq,f,dfdq)
		    IF(no_croiss)PRINT*,'Problème localisé en 19 dans lim_zc'	    
c       WRITE(*,2000)f
		    pt=EXP(f(1)) ; dlpp=1.d0
		    IF(pturb)THEN !avec pression turbulente
		       p=EXP(f(Ipg))
		       IF(der) dlpp=dfdq(Ipg)/dfdq(1) !dlpp=dln Pgaz/dln Ptot
		    ELSE	!sans pression turbulente
		       p=pt
		    ENDIF
		    t=EXP(f(2))
		    IF(en_masse)THEN
		       r=SQRT(ABS(f(3))) ; l=SQRT(ABS(f(4)))**3
		       m=SQRT(ABS(f(5)))**3
		    ELSE
		       r=ABS(f(3)) ; l=f(4) ; m=ABS(f(5)) ; f(5)=m**(2.d0/3.d0)    
		    ENDIF
		    m13=MAX(m**(1.d0/3.d0),1.d-5) !en 0 décalage calcul de dérivée
		    CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1		 knotc,.TRUE.,MAX(mc(1),MIN(f(5),mc(n_ch))),lc,comp,dcomp)
		    IF(no_croiss)PRINT*,'Problème localisé en 20 dans lim_zc'	
		    DO j=1,nchim !dX/dm
		       dcomp(j)=dcomp(j)*2.d0/3.d0/m13
		    ENDDO
		    
c       vitesse angulaire
		    
		    SELECT CASE(Krot)
		    CASE(0,1,2)
		       w=wrot
		    CASE(3)
		       CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,
	1		    MAX(mrot(1),MIN(f(5),mrot(n_rot))),lc,frot,dfrot)	
		       IF(no_croiss)PRINT*,'Problème localisé en 21 dans lim_zc'	     
		       w=frot(1)
		    END SELECT
		    
		    d_grad(3)=dgrad(pt,p,t,dlpp,comp,m,l,r,dcomp,w)    
		    
c       PRINT*,j1,j2,j3,compt ; WRITE(*,2000)d_grad
		    
		    IF(d_grad(2)*d_grad(3) < 0.d0)THEN
		       j1=j3 ; d_grad(1)=d_grad(3)
		    ELSEIF(d_grad(1)*d_grad(3) <= 0.d0)THEN   
		       j2=j3 ; d_grad(2)=d_grad(3)
		    ENDIF
		    j3=(j1+j2)/2.d0 ; compt=compt+1
		    IF(compt > dlt_n)THEN
		       PRINT*,'dichotomie dans lim_zc, compt>',compt
		       PRINT*,j1,d_grad(1),j2,d_grad(2)
		       CALL pause('abandon')
		       STOP
		    ENDIF
		 ENDDO		!WHILE
		 jlim(i)=NINT(j3) !déplacement après de ENDDO le 02 02 00    
		 new(jlim(i))=j3
c       WRITE(*,2000)new(jlim(1:lim)),ABS(j1-j2),loc_zc
	      ENDDO		!i
	      
c       WRITE(*,2000)new(jlim(1:lim)) ; CALL pause('new jlim')
	      
c---------écritures-------------------------------  
	      
	      PRINT*
	      SELECT CASE(langue)
	      CASE('english')
		 PRINT*,'Limits RZ/CZ after the refinement of their locations :'	  
	      CASE DEFAULT
		 PRINT*,'Limites ZR/ZC après affinement des localisations :'
	      END SELECT	  
	      
	      DO i=1,lim	!masses des limites
		 CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,new(jlim(i)),lq,
	1	      f,dfdq)
c       j3=REAL(jlim(i),dp)
c       CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,j3,lq,f,dfdq)	
		 IF(no_croiss)PRINT*,'Problème localisé en 22 dans lim_zc'

c----------------NEW --- NEW ---- NEW ------------------------------------	
		 pt=EXP(f(1)) ; dlpp=1.d0
		 IF(pturb)THEN	!avec pression turbulente
		    p=EXP(f(Ipg))
		    IF(der)dlpp=dfdq(Ipg)/dfdq(1) !dlpp=dln Pgaz/dln Ptot
		 ELSE		!sans pression turbulente
		    p=pt
		 ENDIF
		 t=EXP(f(2))
		 IF(en_masse)THEN
		    l=SQRT(ABS(f(4)))**3 ; r_zc(i)=SQRT(ABS(f(3)))
		    m_zc(i)=SQRT(ABS(f(5)))**3
		 ELSE
		    l=f(4) ; r_zc(i)=ABS(f(3)) ; m_zc(i)=ABS(f(5))
		    f(5)=m_zc(i)**(2.d0/3.d0)
		 ENDIF  
		 m13=MAX(m_zc(i)**(1.d0/3.d0),1.d-5) !en 0 decalage
		 CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1	      knotc,.TRUE.,MAX(mc(1),MIN(f(5),mc(n_ch))),lc,comp,dcomp)
		 IF(no_croiss)PRINT*,'Problème localisé en 25 dans lim_zc'
		 DO j=1,nchim	!dX/dm
		    dcomp(j)=dcomp(j)*2.d0/3.d0/m13
		 ENDDO
		 IF(no_croiss)PRINT*,'Problème localisé en 25 dans lim_zc'
		 DO j=1,nchim	!dX/dm
		    dcomp(j)=dcomp(j)*2.d0/3.d0/m13
		 ENDDO
c       PRINT*,'m_zc,r_zc,new(jlim(i))',lconv(i),i
c       WRITE(*,2000)m_zc(i),r_zc(i),new(jlim(i))
c       WRITE(*,2000)p,t,r_zc(i),l,m_zc(i),comp(1)
		 
c       vitesse angulaire
		 
		 SELECT CASE(Krot)
	         CASE(0,1,2)
		    w=wrot
		 CASE(3)
		    CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,
	1		 MAX(mrot(1),MIN(f(5),mrot(n_rot))),lc,frot,dfrot)
		    IF(no_croiss)PRINT*,'Problème localisé en 26 dans lim_zc'	
		    w=frot(1)
		 END SELECT
		 
		 CALL thermo(pt,p,t,m_zc(i),l,r_zc(i),dlpp,comp,dcomp,
	1	      ro,drop,drot,drox,u,dup,dut,dux,
	2	      grad,dgradpt,dgradp,dgradt,dgradx,dgradm,dgradl,dgradr,dgradlpp,
	3	      gam,dgampt,dgamp,dgamt,dgamx,dgamm,dgaml,dgamr,dgamlpp,
	4	      epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
	5	      delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	6	      gradad,dgradadp,dgradadt,dgradadx,
	7	      hp,dhppt,dhpp,dhpt,dhpx,dhpr,dhpm,
	8	      gradrad,alfa,beta,gamma1,radiatif)
		 
		 pr_zc(i)=p ; hp_zc(i)=hp
c-----------------------------------------------------------------------
c       IF(en_masse)THEN
c       r_zc(i)=SQRT(ABS(f(3))) ; m_zc(i)=SQRT(ABS(f(5)))**3
c       ELSE
c       r_zc(i)=ABS(f(3)) ; m_zc(i)=ABS(f(5))
c       ENDIF
c-----------------------------------------------------------------------
		 
		 IF(lconv(i))THEN !début de ZC
		    SELECT CASE(langue)
		    CASE('english')
		       WRITE(*,1024)i,INT((new(jlim(i))-INT(new(jlim(i))))*100.),
	1		    INT(new(jlim(i))),INT(new(jlim(i)))+1,r_zc(i)/rn,
	2		    m_zc(i)/mstar,1.d0-m_zc(i)/mstar
 1024		       FORMAT('for the limit RZ/CZ #',i4,' at ',i3,'% between shell',
	1		    i5,' and',i5,', beginning of CZ',/,'R_zc/Rtot=',es10.3,
	2		    ', M_zc/Mstar=',es10.3,', 1-M_zc/Mstar=',es10.3)	    
		    CASE DEFAULT	   
		       WRITE(*,24)i,INT((new(jlim(i))-INT(new(jlim(i))))*100.),
	1		    INT(new(jlim(i))),INT(new(jlim(i)))+1,r_zc(i)/rn,
	2		    m_zc(i)/mstar,1.d0-m_zc(i)/mstar
 24		       FORMAT('pour limite ZR/ZC #',i4,' à ',i3,'% entre couches',i5,
	1		    ' et',i5,', début de ZC',/,'R_zc/Rtot=',es10.3,
	2		    ', M_zc/Mstar=',es10.3,', 1-M_zc/Mstar=',es10.3)
		    END SELECT
		 ELSE
		    SELECT CASE(langue)
		    CASE('english')
		       WRITE(*,1025)i,INT((new(jlim(i))-INT(new(jlim(i))))*100.),
	1		    INT(new(jlim(i))),INT(new(jlim(i)))+1,r_zc(i)/rn,
	2		    m_zc(i)/mstar,1.d0-m_zc(i)/mstar
 1025		       FORMAT('for the limit RZ/CZ #',i4,' at ',i3,'% between shell',
	1		    i5,' and #',i5,', end of CZ',/,'R_zc/Rtot=',es10.3,
	2		    ', M_zc/Mstar=',es10.3,', 1-M_zc/Mstar=',es10.3)
		    CASE DEFAULT	   	   
		       WRITE(*,25)i,INT((new(jlim(i))-INT(new(jlim(i))))*100.),
	1		    INT(new(jlim(i))),INT(new(jlim(i)))+1,r_zc(i)/rn,
	2		    m_zc(i)/mstar,1.d0-m_zc(i)/mstar
 25		       FORMAT('pour limite ZR/ZC #',i4,' à ',i3,'% entre couches #',
	1		    i5,' et #',i5,', fin de ZC',/,'R_zc/Rtot=',es10.3,
	2		    ', M_zc/Mstar=',es10.3,', 1-M_zc/Mstar=',es10.3)
		    END SELECT	
		 ENDIF
	      ENDDO		!i
	      
	   ELSE			!no_supp=.FALSE
	      SELECT CASE(langue)
	      CASE('english')
		 PRINT*,'no refinement of the localization of the limits CZ/RZ'
		 PRINT*,'because too closed limits have been suppressed'    
		 PRINT*,'limits CZ/RZ used:'	 
	      CASE DEFAULT
		 PRINT*,'Pas d''affinement de la localisation des limites ZR/ZC'
		 PRINT*,'à cause de la suppression de limites ZR/ZC trop proches'    
		 PRINT*,'Limites ZR/ZC retenues:'
	      END SELECT
	      
	      DO i=1,lim	!masses des limites
		 CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,new(jlim(i)),
	1	      lq,f,dfdq)
c       j3=REAL(jlim(i),dp)
c       CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,j3,lq,f,dfdq)	
		 IF(no_croiss)PRINT*,'Problème localisé en 23 dans lim_zc'
c----------------NEW --------- NEW --------- NEW ------------------------	
		 pt=EXP(f(1)) ; dlpp=1.d0
		 IF(pturb)THEN	!avec pression turbulente
		    p=EXP(f(Ipg))
		    IF(der)dlpp=dfdq(Ipg)/dfdq(1) !dlpp=dln Pgaz/dln Ptot
		 ELSE		!sans pression turbulente
		    p=pt
		 ENDIF
		 t=EXP(f(2))
		 IF(en_masse)THEN
		    l=SQRT(ABS(f(4)))**3 ; r_zc(i)=SQRT(ABS(f(3)))
		    m_zc(i)=SQRT(ABS(f(5)))**3
		 ELSE
		    l=f(4) ; r_zc(i)=ABS(f(3)) ; m_zc(i)=ABS(f(5))
		    f(5)=m_zc(i)**(2.d0/3.d0)
		 ENDIF  
		 m13=MAX(m_zc(i)**(1.d0/3.d0),1.d-5) !en 0 decalage
		 CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1	      knotc,.TRUE.,MAX(mc(1),MIN(f(5),mc(n_ch))),lc,comp,dcomp)
		 IF(no_croiss)PRINT*,'Problème localisé en 25 dans lim_zc'
		 DO j=1,nchim	!dX/dm
		    dcomp(j)=dcomp(j)*2.d0/3.d0/m13
		 ENDDO
		 IF(no_croiss)PRINT*,'Problème localisé en 25 dans lim_zc'
		 DO j=1,nchim	!dX/dm
		    dcomp(j)=dcomp(j)*2.d0/3.d0/m13
		 ENDDO
c       PRINT*,'m_zc,r_zc,new(jlim(i))',lconv(i),i
c       WRITE(*,2000)m_zc(i),r_zc(i),new(jlim(i))
c       WRITE(*,2000)p,t,r_zc(i),l,m_zc(i),comp(1)
		 
c       vitesse angulaire
		 
		 SELECT CASE(Krot)
	         CASE(0,1,2)
		    w=wrot
		 CASE(3)
		    CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,
	1		 MAX(mrot(1),MIN(f(5),mrot(n_rot))),lc,frot,dfrot)
		    IF(no_croiss)PRINT*,'Problème localisé en 26 dans lim_zc'	
		    w=frot(1)
		 END SELECT
		 
		 CALL thermo(pt,p,t,m_zc(i),l,r_zc(i),dlpp,comp,dcomp,
	1	      ro,drop,drot,drox,u,dup,dut,dux,
	2	      grad,dgradpt,dgradp,dgradt,dgradx,dgradm,dgradl,dgradr,dgradlpp,
	3	      gam,dgampt,dgamp,dgamt,dgamx,dgamm,dgaml,dgamr,dgamlpp,
	4	      epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
	5	      delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	6	      gradad,dgradadp,dgradadt,dgradadx,
	7	      hp,dhppt,dhpp,dhpt,dhpx,dhpr,dhpm,
	8	      gradrad,alfa,beta,gamma1,radiatif)
		 
		 pr_zc(i)=p ; hp_zc(i)=hp
c-----------------end of new  -----------------------------
c       IF(en_masse)THEN
c       r_zc(i)=SQRT(ABS(f(3))) ; m_zc(i)=SQRT(ABS(f(5)))**3
c       ELSE
c       r_zc(i)=ABS(f(3)) ; m_zc(i)=ABS(f(5))
c       ENDIF
c       PRINT*,'masse r_zc,m_zc',i,jlim(i),new(jlim(i))
		 
		 IF(lconv(i))THEN !début de ZC
		    WRITE(*,24)i,INT((new(jlim(i))-INT(new(jlim(i))))*100.),
	1		 INT(new(jlim(i))),INT(new(jlim(i)))+1,r_zc(i)/rn,
	2		 m_zc(i)/mstar,1.d0-m_zc(i)/mstar
		 ELSE
		    WRITE(*,25)i,INT((new(jlim(i))-INT(new(jlim(i))))*100.),
	1		 INT(new(jlim(i))),INT(new(jlim(i)))+1,r_zc(i)/rn,
	2		 m_zc(i)/mstar,1.d0-m_zc(i)/mstar  
		 ENDIF
	      ENDDO		!i
	   ENDIF		!no_supp
	   
c       PRINT*,lconv(1:lim) ; CALL pause('lconv')
	   
c       calcul de r_zc et de r_ov si overshooting ou PC
	   
	   IF(ovsht)THEN	!s'il y a overshooting ou PC
	      DO i=1,lim
c       PRINT*,'limite ',i,' ',lconv(i)
		 r_zc(i)=-100d0 ; r_ov(i)=-100d0
		 CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,new(jlim(i)),
	1	      lq,f,dfdq)
c       j3=REAL(jlim(i),dp)
c       CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,j3,lq,f,dfdq)
		 IF(no_croiss)PRINT*,'Problème localisé en 24 dans lim_zc'
		 pt=EXP(f(1)) ; dlpp=1.d0
		 IF(pturb)THEN	!avec pression turbulente
		    p=EXP(f(Ipg))
		    IF(der)dlpp=dfdq(Ipg)/dfdq(1) !dlpp=dln Pgaz/dln Ptot
		 ELSE		!sans pression turbulente
		    p=pt
		 ENDIF
		 t=EXP(f(2))
		 IF(en_masse)THEN
		    l=SQRT(ABS(f(4)))**3 ; r_zc(i)=SQRT(ABS(f(3)))
		    m_zc(i)=SQRT(ABS(f(5)))**3
		 ELSE
		    l=f(4) ; r_zc(i)=ABS(f(3)) ; m_zc(i)=ABS(f(5))
		    f(5)=m**(2.d0/3.d0)
		 ENDIF  
		 m13=MAX(m_zc(i)**(1.d0/3.d0),1.d-5) !en 0 decalage
		 CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1	      knotc,.TRUE.,MAX(mc(1),MIN(f(5),mc(n_ch))),lc,comp,dcomp)
		 IF(no_croiss)PRINT*,'Problème localisé en 25 dans lim_zc'
		 DO j=1,nchim	!dX/dm
		    dcomp(j)=dcomp(j)*2.d0/3.d0/m13
		 ENDDO
c       PRINT*,'m_zc,r_zc,new(jlim(i))',lconv(i),i
c       WRITE(*,2000)m_zc(i),r_zc(i),new(jlim(i))
c       WRITE(*,2000)p,t,r_zc(i),l,m_zc(i),comp(1)
		 
c       vitesse angulaire
		 
		 SELECT CASE(Krot)
	         CASE(0,1,2)
		    w=wrot
		 CASE(3)
		    CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,
	1		 MAX(mrot(1),MIN(f(5),mrot(n_rot))),lc,frot,dfrot)
		    IF(no_croiss)PRINT*,'Problème localisé en 26 dans lim_zc'	
		    w=frot(1)
		 END SELECT
		 
		 CALL thermo(pt,p,t,m_zc(i),l,r_zc(i),dlpp,comp,dcomp,
	1	      ro,drop,drot,drox,u,dup,dut,dux,
	2	      grad,dgradpt,dgradp,dgradt,dgradx,dgradm,dgradl,dgradr,dgradlpp,
	3	      gam,dgampt,dgamp,dgamt,dgamx,dgamm,dgaml,dgamr,dgamlpp,
	4	      epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
	5	      delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	6	      gradad,dgradadp,dgradadt,dgradadx,
	7	      hp,dhppt,dhpp,dhpt,dhpx,dhpr,dhpm,
	8	      gradrad,alfa,beta,gamma1,radiatif)
		 
c       overshoot: ovshti(s)*Hp ou ovshti*kip de JPZ si T>5d5 à la limite
		 
		 IF(t > 5.d5)THEN
		    IF(lconv(i))THEN !début de ZC, ext. vers le bas
		       IF(ovshti > 0.d0)THEN
			  IF(jpz)THEN
			     kip=(3.d0-t*(drot/ro+dkapt/kap))*gradad-p*(drop/ro+dkapp/kap)
			     r_ov(i)=MAX(r_zc(i)-hp*ovshti/rsol/kip,0.d0)
c       PRINT*,'i/r_zc(i),hp*ovshti/rsol/kip,hp,ovshti,rsol,kip',i
c       WRITE(*,2000)r_zc(i),hp*ovshti/rsol/kip,hp,ovshti,rsol,kip
c       WRITE(*,2000)r_ov(i)
			  ELSE	!ext. vers le bas: 0.<r_ov<r_zc
			     r_ov(i)=MAX(r_zc(i)-hp*ovshti/rsol,0.d0)
			  ENDIF	!jpz
c       PRINT*,'overshoot inférieur'
c       WRITE(*,2000)r_ov(i),r_zc(i),hp
		       ENDIF	!ovshti
		    ELSE	!fin de ZC overshoot supérieur: rext>r_ov>r_zc 
		       IF(ovshts > 0.d0)THEN !il y a ovsht sup.
			  IF(jpz)THEN
			     r_ov(i)=r_zc(i)*(1.d0+ovshts) !overshoot sup.= R_noyauXovshts
			  ELSE
			     r_ov(i)=MIN(r_zc(i)+hp*ovshts/rsol,r_zc(i)*(1.d0+ovshts))
			  ENDIF	!jpz
			  r_ov(i)=MIN(r_ov(i),rn)
c       PRINT*,'ovshts' ; WRITE(*,2000)r_ov(i),r_zc(i),hp,rn
		       ENDIF	!ovshts
		    ENDIF	!lconv
c       CALL pause('après if sur lconv')
		    
c       masse d'overshooting ou de pénétration
c       ATTENTION 'r2 ' et non 'r2'
		    
		    IF(r_ov(i) >= 0.d0)THEN !r_ov=-100: no overshoot en r_zc(i)
		       IF(en_masse)THEN
			  CALL inter('r2 ',bp,q,qt,n_qs,knot,r_ov(i)**2,f,dfdq,r2,m23) 
			  m_zc(i)=SQRT(ABS(f(5)))**3
		       ELSE
			  CALL inter('r2 ',bp,q,qt,n_qs,knot,r_ov(i),f,dfdq,r2,m23)
			  m_zc(i)=ABS(f(5)) 
		       ENDIF
		       j_new=f(6) !astuce pour passer l'indice
c       WRITE(*,2000)f(1:ne) 
		       
c       l'ancien new(jlim(i)) est désormais un point ordinaire
c       qui redevient l'indice entier jlim(i)
c       la limite est désormais en jlim(i)= + proche entier de j_new
c       PRINT*,'jlim(i),new(jlim(i)),j_new',jlim(i),new(jlim(i)),j_new
		       
		       jlim(i)=NINT(j_new)    
		       new(jlim(i))=j_new !nouvel indice pour les limites ZR/ZC
		       
c       PRINT*,'jlim(i),new(jlim(i)),j_new',jlim(i),new(jlim(i)),j_new
c       PRINT*,'j_new/M_ov,R_ov',j_new ; WRITE(*,2000)m_zc(i),r_ov(i)
c       CALL pause('après j_new')
		    ENDIF	!r_ov .ne. 0
		 ELSE
		    r_ov(i)=r_zc(i)
c       CALL pause('t > 1.d5')
		 ENDIF		!t > 1.d5
	      ENDDO		!i sur lim
	      
c       il ne doit pas y avoir de chevauchement après les extensions
	      
c       PRINT*,lim,ovsht,(jlim(i),i=1,lim),(lconv(i),i=1,lim)
c       WRITE(*,2000)(r_ov(j),j=1,lim)
c       WRITE(*,2000)(r_zc(j),j=1,lim)
c       WRITE(*,2000)(m_zc(j),j=1,lim)
	      
	      IF(lconv(1))THEN	!début de ZC à la première limite
		 i=2
	      ELSE
		 i=1
	      ENDIF
	      DO WHILE (i < lim)
		 IF(jlim(i) >= jlim(i+1))THEN
		    DO j=i,lim-2 !chevauchement: suppression d'une ZR décalage
		       jlim(j)=jlim(j+2) ; m_zc(j)=m_zc(j+2) ; lconv(j)=lconv(j+2)
		       r_zc(j)=r_zc(j+2) ; r_ov(j)=r_ov(j+2)
		    ENDDO
		    lim=lim-2
		 ELSE
		    i=i+2
		 ENDIF
	      ENDDO		!WHILE
	      
c       à cause d'un overshoot, le modèle devient totalement convectif
	      
	      IF(lim == 1 .AND. jlim(1) == 1 .AND. lconv(1))THEN
		 jlim(1)=n_qs ; lconv(1)=.FALSE. ; tot_conv=.TRUE. ; r_zc(1)=rn
	      ELSE
		 tot_conv=.FALSE.
		 
c       pas de limite ZR/ZC aux extrémités
		 
		 IF(jlim(1) == 1)THEN !en 1
		    DO i=1,lim-1
		       jlim(i)=jlim(i+1) ; lconv(i)=lconv(i+1) ; m_zc(i)=m_zc(i+1)
		       r_zc(i)=r_zc(i+1) ; r_ov(i)=r_ov(i+1)
		    ENDDO
		    lim=lim-1
		 ENDIF
		 IF(jlim(lim) == n_qs)lim=lim-1 !en n
	      ENDIF
	      SELECT CASE(langue)
	      CASE('english')
		 PRINT*,'limits RZ/CZ after examination of overshooting :'
	      CASE DEFAULT   
		 PRINT*,'limites ZR/ZC après examen de l''overshooting :'
	      END SELECT
	      
c       PRINT*,lim,jlim(1:lim)
	      
	      IF(tot_conv)THEN
		 SELECT CASE(langue)
	         CASE('english')
		    WRITE(*,1033)
 1033		    FORMAT('fully convective model',/)	  
		 CASE DEFAULT	  
		    WRITE(*,33)
 33		    FORMAT('modèle totalement convectif',/)
		 END SELECT
	      ELSE
		 DO j=1,lim
		    IF(lconv(j))THEN !début de ZC
		       IF(ovshti > 0.d0)THEN
			  IF(jpz)THEN
			     SELECT CASE(langue)
			     CASE('english')
				WRITE(*,1031)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
 1031				FORMAT('Convective penetration of JPZ limit CZ/RZ #',i4,
	1			     ', shell #',i5,/,3x,'radius CZ/Rtot=',es10.3,
	2			     ', reduced radius/Rtot=',es10.3,/,
	3			     'm/Mstar at the limit of the convective penetration=',
	4			     es10.3,/)	      
			     CASE DEFAULT		      
				WRITE(*,31)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
 31				FORMAT('Pénétration Convective JPZ limite ZC/ZR #',i4,
	1			     ', couche #',i5,/,3x,'rayon ZC/Rtot=',es10.3,
	2			     ', rayon reduit/Rtot=',es10.3,/,
	3			     'm/Mstar à la limite de la pénétration convective=',
	4			     es10.3,/)
			     END SELECT
			  ELSE	!jpz
			     SELECT CASE(langue)
			     CASE('english')
				WRITE(*,1002)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
 1002				FORMAT('negative overshooting  limit CZ/RZ #',i4,
	1			     ', shell #',i5,/,'rayon ZC/Rtot=',es10.3,
	2			     ', reduced raius/Rtot=',es10.3,/,
	3			     'm/Mstar at the limit of overshooting=',es10.3,/)
			     CASE DEFAULT		      
				WRITE(*,2)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
 2				FORMAT('Overshooting inférieur limite ZC/ZR #',i4,
	1			     ', couche #',i5,/,'rayon ZC/Rtot=',es10.3,
	2			     ', rayon reduit/Rtot=',es10.3,/,
	3			     'm/Mstar à la limite de l''overshooting=',es10.3,/)
			     END SELECT		
			  ENDIF	!jpz
		       ELSE	!sans overshoot
			  SELECT CASE(langue)	     
		          CASE('english')
			     WRITE(*,1011)j,jlim(j),r_zc(j)/rn,m_zc(j)/mstar
 1011			     FORMAT('Limit without overshooting nor PC, limit RZ/CZ #:',
	1			  i4,', shell #',i5,/,'radius CZ/Rtot=',es10.3,
	2			  ', m/Mstar at the end of the CZ=',es10.3,/)	      
			  CASE DEFAULT		      
			     WRITE(*,11)j,jlim(j),r_zc(j)/rn,m_zc(j)/mstar
 11			     FORMAT('Limite sans overshooting ni PC, limite ZR/ZC #:',i4,
	1			  ', couche #',i5,/,'rayon ZC/Rtot=',es10.3,
	2			  ', m/Mstar à la fin de la ZC=',es10.3,/)
			  END SELECT	
		       ENDIF	!ovshti
		    ELSE	!fin de ZC
		       IF(ovshts > 0.d0)THEN
			  IF(jpz)THEN
			     SELECT CASE(langue)	      
			     CASE('english')
				WRITE(*,1012)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
 1012				FORMAT('positive overshoot * R core, limit CZ/RZ #',i4,
	1			     ', shell #',i5,/,3x,'radius CZ/Rtot=',es10.3,
	2			     ', extented radius/Rtot=',es10.3,/,
	3			     'm/Mstar at the limit of overshoot=',es10.3,/)	      
			     CASE DEFAULT		      
				WRITE(*,12)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
 12				FORMAT('Overshoot supérieur * R noyau, limite ZC/ZR #',i4,
	1			     ', couche #',i5,/,3x,'rayon ZC/Rtot=',es10.3,
	2			     ', rayon etendu/Rtot=',es10.3,/,
	3			     'm/Mstar à la limite de l''overshoot=',es10.3,/)
			     END SELECT
			  ELSE	!jpz	    
			     IF(r_ov(j) == r_zc(j))THEN
				SELECT CASE(langue)	       
			        CASE('english')
				   WRITE(*,1013)j,jlim(j),r_zc(j)/rn,m_zc(j)/mstar
 1013				   FORMAT('Limit CZ/RZ without positive overshoot #',i4,
	1				', shell #',i5,/,'because T < 5.d5, radius CZ/Rtot=',
	2				es10.3,', m/Mstar on the limit =',es10.3,/)	      
				CASE DEFAULT		      
				   WRITE(*,13)j,jlim(j),r_zc(j)/rn,m_zc(j)/mstar
 13				   FORMAT('Limite ZC/ZR sans overshooting supérieur #',i4,
	1				', couche #',i5,/,'car T < 5.d5, rayon ZC/Rtot=',es10.3,
	2				', m/Mstar à la limite =',es10.3,/)
				END SELECT
			     ELSE
				SELECT CASE(langue)	       
			        CASE('english')
				   WRITE(*,1001)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
 1001				   FORMAT('Overshooting positive * Hp, limit CZ/RZ #',i4,
	1				', shell #',i5,/,'radius CZ/Rtot=',es10.3,
	2				', extended radius/Rtot=',es10.3,/,
	3				'm/Mstar at the limit of overshooting=',es10.3,/)	      
				CASE DEFAULT		      
				   WRITE(*,1)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
 1				   FORMAT('Overshooting supérieur * Hp, limite ZC/ZR #',i4,
	1				', couche #',i5,/,'rayon ZC/Rtot=',es10.3,
	2				', rayon etendu/Rtot=',es10.3,/,
	3				'm/Mstar à la limite de l''overshooting=',es10.3,/)
				END SELECT
			     ENDIF
			  ENDIF	!jpz
		       ELSE	!pas d'overshoot
			  SELECT CASE(langue)	     
		          CASE('english')
			     WRITE(*,1021)j,jlim(j),r_zc(j)/rn,m_zc(j)/mstar
 1021			     FORMAT('Limit without overshooting nor PC, limit RZ/CZ #',
	1			  i4,', shell #',i5,/,'radius CZ/Rtot=',es10.3,
	2			  ', m/Mstar at the beginning of the CZ=',es10.3,/)	      
			  CASE DEFAULT		      
			     WRITE(*,21)j,jlim(j),r_zc(j)/rn,m_zc(j)/mstar
 21			     FORMAT('Limite sans overshooting ni PC, limite ZR/ZC #',i4,
	1			  ', couche #',i5,/,'rayon ZC/Rtot=',es10.3,
	2			  ', m/Mstar au début de la ZC=',es10.3,/)
			  END SELECT		
		       ENDIF	!ovshts
		    ENDIF	!lconv
		 ENDDO		!j limites ZR/ZC
	      ENDIF		!tot_conv
	   ENDIF		!ovsht
	   
c       CALL pause('après les overshoot')
	   
c       détermination de la fonction de poids fac 
	   
c       interpolation aux nouveaux q de P, T, R, L, M pour le calcul de la
c       fonction de répartition aux new
c       PRINT*,'ctel,ctep' ; WRITE(*,2000)ctel,ctep,ctem,cter,ctet
c       PRINT* ; PRINT*,'lim/new/jlim',lim ; PRINT*,(new(jlim(i)),i=1,lim)
c       PRINT*,(jlim(i),i=1,lim) ; WRITE(*,2000)new    
	   
	   IF(.NOT.(tot_rad .OR. tot_conv))THEN
	      
	      ALLOCATE(esp(n_qs))
	      DO i=1,n_qs
		 CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,new(i),lq,f,dfdq)
		 IF(no_croiss)PRINT*,'Problème localisé en 27 dans lim_zc'
		 esp(i)=ctep*f(1)+ctet*f(2)+cter*f(3)+ctel*f(4)+ctem*f(5)
c       WRITE(*,2000)(f(j),j=1,5),esp(i),new(i)
		 IF(i > 1)THEN
		    DO WHILE(esp(i) <= esp(i-1))
		       esp(i)=esp(i)+0.01d0
		    ENDDO
		 ENDIF
	      ENDDO
	      
c       esp(i) est la valeur de la fonction de répartition aux nouveaux q
c       déterminée avec les anciens facteurs fac(i)
c       calcul des nouveaux fac(i) de facon à amener un noeud sur chaque
c       nouveau q
	      
	      DO i=1,n_qs-1
		 fac(i)=bp(6,1)/(esp(i+1)-esp(i))
c       WRITE(*,2000)REAL(i),fac(i) 
	      ENDDO
c       WRITE(*,2000)bp(6,1) ; PRINT*,n_qs
c       WRITE(*,2000)(fac(i),new(i),i=1,n_qs-1)
c       PRINT*,(lim,jlim(i),i=1,lim)
c       PRINT*,bp(6,1),(fac(jlim(i)),new(jlim(i)),i=1,lim)
c       PRINT*,(esp(jlim(i)-1),esp(jlim(i)),esp(jlim(i)+1),i=1,lim)
c       CALL pause('esp')
	      DEALLOCATE(new,esp)
	      
c       formation des limites ZR/ZC pour la convection
c       avec un coeur convectif on prend l=Rcoeur-r, ce qui sera
c       réalisé en prenant r_zc_conv(id_conv)=-r_zc(1)
	      
	      IF(lconv(1))THEN
		 id_conv=1 ; r_zc_conv(id_conv)=r_zc(1) !coeur pas convectif
		 pr_zc_conv(id_conv)=pr_zc(1) ; hp_zc_conv(id_conv)=hp_zc(1)
	      ELSE
		 id_conv=0 ; r_zc_conv(id_conv)=0.d0 !r_zc_conv(id_conv)=-r_zc(1)   !coeur convectif
		 pr_zc_conv(id_conv)=EXP(bp(1,1)) ; hp_zc_conv(id_conv)=-hp_zc(1)
	      ENDIF
	      DO i=1,lim
		 r_zc_conv(i+id_conv)=r_zc(i)
		 pr_zc_conv(i+id_conv)=pr_zc(i) ; hp_zc_conv(i+id_conv)=hp_zc(i)
	      ENDDO
	      IF(lconv(lim))THEN
		 if_conv=lim+1	!enveloppe totalement convective
		 r_zc_conv(if_conv)=MAX(rstar,r_zc(lim))
c       hp_zc_conv(if_conv)=hp_zc(lim)
c       pr_zc_conv(if_conv)=0.d0

		 
c******************************novo ***********************************
		 CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,REAL(n_qs-1,dp),
	1	      lq,f,dfdq)
c       j3=REAL(jlim(i),dp)
c       CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,j3,lq,f,dfdq)	
		 IF(no_croiss)PRINT*,'Problème localisé en 23 dans lim_zc'
c----------------NEW --------- NEW --------- NEW ------------------------	
		 pt=EXP(f(1)) ; dlpp=1.d0
		 IF(pturb)THEN	!avec pression turbulente
		    p=EXP(f(Ipg))
		    IF(der)dlpp=dfdq(Ipg)/dfdq(1) !dlpp=dln Pgaz/dln Ptot
		 ELSE		!sans pression turbulente
		    p=pt
		 ENDIF
		 t=EXP(f(2))
		 IF(en_masse)THEN
		    l=SQRT(ABS(f(4)))**3
		 ELSE
		    l=f(4)
		    f(5)=f(5)**(2.d0/3.d0)
		 ENDIF  
		 m13=MAX(SQRT(f(5)),1.d-5) !en 0 decalage
		 CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1	      knotc,.TRUE.,MAX(mc(1),MIN(f(5),mc(n_ch))),lc,comp,dcomp)
		 IF(no_croiss)PRINT*,'Problème localisé en 25 dans lim_zc'
		 DO j=1,nchim	!dX/dm
		    dcomp(j)=dcomp(j)*2.d0/3.d0/m13
		 ENDDO
		 IF(no_croiss)PRINT*,'Problème localisé en 25 dans lim_zc'
		 DO j=1,nchim	!dX/dm
		    dcomp(j)=dcomp(j)*2.d0/3.d0/m13
		 ENDDO
c       PRINT*,'m_zc,r_zc,new(jlim(i))',lconv(i),i
c       WRITE(*,2000)m_zc(i),r_zc(i),new(jlim(i))
c       WRITE(*,2000)p,t,r_zc(i),l,m_zc(i),comp(1)
		 
c       vitesse angulaire
		 
		 SELECT CASE(Krot)
	         CASE(0,1,2)
		    w=wrot
		 CASE(3)
		    CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,
	1		 MAX(mrot(1),MIN(f(5),mrot(n_rot))),lc,frot,dfrot)
		    IF(no_croiss)PRINT*,'Problème localisé en 26 dans lim_zc'	
		    w=frot(1)
		 END SELECT
		 
		 CALL thermo(pt,p,t,mstar,l,r_zc_conv(if_conv),dlpp,comp,dcomp,
	1	      ro,drop,drot,drox,u,dup,dut,dux,
	2	      grad,dgradpt,dgradp,dgradt,dgradx,dgradm,dgradl,dgradr,dgradlpp,
	3	      gam,dgampt,dgamp,dgamt,dgamx,dgamm,dgaml,dgamr,dgamlpp,
	4	      epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
	5	      delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	6	      gradad,dgradadp,dgradadt,dgradadx,
	7	      hp,dhppt,dhpp,dhpt,dhpx,dhpr,dhpm,
	8	      gradrad,alfa,beta,gamma1,radiatif)
		 
		 pr_zc_conv(if_conv)=p ; hp_zc_conv(if_conv)=hp
c******************end of new ******************************
		 
		 
	      ELSE
		 if_conv=lim	!enveloppe non totalement convective
	      ENDIF
	   ENDIF
	ENDIF			!lim=0

	IF(tot_rad .OR. tot_conv)THEN
	   fac=1.d0 ; id_conv=0 ; r_zc_conv(id_conv)=0.d0
	   if_conv=1 ; r_zc_conv(if_conv)=rstar ; hp_zc_conv(if_conv)=hp_zc(1)
	   pr_zc_conv(id_conv)=EXP(bp(1,1)) ; pr_zc_conv(if_conv)=pr_zc(1)
	ENDIF

c	PRINT*,n_qs ; WRITE(*,2000)fac ; CALL pause('sortie lim_zc')

	SELECT CASE(langue)
	CASE('english')
	   WRITE(*,1022)
 1022	   FORMAT('---------Search of the limits ZR/ZC(end)-----------')
	CASE DEFAULT
	   WRITE(*,22)
 22	   FORMAT('---------Recherche de limites ZR/ZC(fin)-----------')
	END SELECT
	
	RETURN
	
	END SUBROUTINE lim_zc
