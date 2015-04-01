
c********************************************************************

	SUBROUTINE lim_zc(no_rep,new_nqs)

c routine private du module mod_static

c répartition des couches et détermination des limites ZR/ZC.
c Le nombre est modifié si new_n /= n

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c 17 09 96 : correction bug dans déclaration de vt; signalé par D. Cordier
c 25 08 97 : mise en place des variables eulériennes
c 19 11 99 : suppression de nh1, nhe1, nhe2, lamb
c 02 02 00 : introduction de no_supp
c 20 02 00 : écriture si non affinement de la loc. des lim. des ZR/ZC

c entrées
c	no_rep : il n'y a pas eu de reprise
c	new_nqs : nouveau nombre de couches

c----------------------------------------------------------------

	USE mod_donnees, ONLY: alpha, aradia, clight, cpturb,
	1 ctel, ctem, ctep, cter, ctet, diffusion, dn_fixe, en_masse,
	2 g, grille_fixe, ihe4, Ipg, jpz, Krot, langue, loc_zc, ledoux,
	3 lsol, msol, m_ch, m_qs, nchim, ne, nrot, ord_qs, ord_rot, ovshti,
	4 ovshts, pi, pnzc, precision, pturb, rsol, r_qs, t_inf
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn, newspl, noedif, no_croiss, pause,
	1 zoning
	USE mod_variables, ONLY: bp, chim, dim_qs,
	1 id_conv, if_conv, inter, jlim, lconv, knot, knotc, knotr, lim,
	2 mc, mct, mc_fixe, mrot, mrott, mstar, m_zc, m23, nc_fixe, n_ch,
	3 n_qs, n_rot, q, qt, rota, rstar, r2, r_ov, r_zc, r_zc_conv, wrot,
	4 tot_conv, tot_rad

	IMPLICIT NONE

	INTEGER, INTENT(in) :: new_nqs
	LOGICAL, INTENT(inout) :: no_rep

	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: new_b, new_esp
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: dfrot, d_grad, frot
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: new_q, new_qt, esp, new
	REAL (kind=dp), DIMENSION(2*pnzc) :: rap
	REAL (kind=dp), DIMENSION(nchim) :: comp, dcomp, depsx
	REAL (kind=dp), DIMENSION(ne) :: dfdq, f
	REAL (kind=dp), DIMENSION(5) :: epsilon

	REAL (kind=dp) :: alfa, beta, pt, dlpp, p, t, r, l, m, w, rn, grp, gr,
	1 j_new, j1, j2, j3, kip, drot, ro, dkapt, kap, gradad, drop,
	2 dkapp, hp, drox, u, dup, dut, dux, grad, dgradpt, dgradp, dgradt,
	3 dgradx, dgradm, dgradl, dgradr, dgradlpp, gam, dgampt, dgamp,
	4 dgamt, dgamx, dgamm, dgaml, dgamr, dgamlpp, depsp,
	5 depst, dkapx, delta, deltap, deltat, deltax, cp, dcpp,
	6 dcpt, dcpx, dgradadp, dgradadt, dgradadx, dhppt, dhpp, dhpt,
	7 dhpx, dhpr, dhpm, gradrad, grad_mu, gamma1

	INTEGER, SAVE :: dlt_n=60, lc=1, ncouches_zrc=5  
	INTEGER :: new_knot, i, im, j, lq, compt

	LOGICAL, SAVE :: der, ovsht, init=.TRUE.
	LOGICAL :: no_supp, radiatif

c-------------------------------------------------------------------

2000	FORMAT(8es10.3)

c-----------------------initialisations début-----------------------
	IF(init)THEN
	 init=.FALSE.
	 der=cpturb < 0.d0  !der=T. on tient compte de dln Pgaz/dln Ptot
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(2,1020)loc_zc ; WRITE(*,1020)loc_zc
1020	  FORMAT('distance of the limit RZ/CZ to the closest node :',1 es10.3,/)
	 CASE DEFAULT
	  WRITE(2,20)loc_zc ; WRITE(*,20)loc_zc
20	  FORMAT('écart des limites ZR/ZC au plus proche noeud :',es10.3,/)
	 END SELECT
	 IF(Krot >= 3)ALLOCATE(dfrot(nrot),frot(nrot))
	 ovsht=ovshti > 0.d0 .OR. ovshts > 0.d0
	 IF(ovsht)THEN	 
	  IF(jpz)THEN
	   IF(ovshti > 0.d0)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1009)INT(ovshti*100) ; WRITE(2,1009)INT(ovshti*100)
1009	     FORMAT('overshooting JPZ beneath the CZ : ',i3,'% Hp/Kip',/)
	    CASE DEFAULT	   
	     WRITE(*,9)INT(ovshti*100) ; WRITE(2,9)INT(ovshti*100)
9	     FORMAT('overshooting inférieur de JPZ des ZC :',i3,'% Hp/Kip')
	    END SELECT
	   ENDIF	   
	   IF(ovshts > 0.d0)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1006)INT(ovshts*100) ; WRITE(2,1006)INT(ovshts*100)
1006	     FORMAT('overshooting JPZ above the CZ : :',i3,'% R noyau',/)
	    CASE DEFAULT	   
	     WRITE(*,6)INT(ovshts*100) ; WRITE(2,6)INT(ovshts*100)
6	     FORMAT('overshooting supérieur de JPZ des ZC :',i3,
	1    '% R noyau',/)
	    END SELECT
	   ENDIF
	  ELSE
	   IF(ovshti > 0.d0)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1007)INT(ovshti*100) ; WRITE(2,1007)INT(ovshti*100)
1007	     FORMAT('overshooting beneath the CZ :',i3,'%  Hp')	    
	    CASE DEFAULT
	     WRITE(*,7)INT(ovshti*100) ; WRITE(2,7)INT(ovshti*100)
7	     FORMAT('overshooting inférieur des ZC :',i3,'%  Hp',/)
	    END SELECT   
	   ENDIF
	   IF(ovshts > 0.d0)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1008)INT(ovshts*100) ; WRITE(2,1008)INT(ovshts*100) 
1008	     FORMAT('overshooting above the CZ :',i3,'%  Hp',/)	    
	    CASE DEFAULT	   
	     WRITE(*,8)INT(ovshts*100) ; WRITE(2,8)INT(ovshts*100)
8	     FORMAT('overshooting supérieur des ZC :',i3,'%  Hp',/)
	    END SELECT 
	   ENDIF
	  ENDIF
	 ELSE
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1005) ; WRITE(2,1005)
1005	   FORMAT('Model without overshooting')	 
	  CASE DEFAULT	 
	   WRITE(*,5) ; WRITE(2,5)
5	   FORMAT('Modèle sans overshooting')
	  END SELECT
	 ENDIF

c pour éviter des ZC internes trop étroites pour les stades avancés	  
	 IF(precision == 'av' .OR. precision == 'ar')ncouches_zrc=15

c allocations initiales
	 ALLOCATE(fac(n_qs-1),d_grad(n_qs))
	ENDIF   
c----------------fin des initialisations----------------------------

c on redéfinit la répartition avant la recherche des limites ZR/ZC
c si le modèle n'est pas estimé "bien réparti" ou si le nb.
c de couches du modèle quasi. stat. doit changer
	I1: IF(.NOT.no_rep .OR. new_nqs /= n_qs) THEN

c les vecteurs fac et d_grad n'ont plus la même dimension
	 DEALLOCATE(fac,d_grad) ; no_rep=.TRUE.

c on met dans esp la répartition précédente
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
c	  WRITE(*,2000)esp(1:ne,i)
	 ENDDO
c	 CALL pause('nouvel esp')

c on dispose les couches pour assurer une répar. approx.
c uniforme de la fonction de répartition sur new_nqs couches 
	 ALLOCATE(new_q(new_nqs))
	 CALL zoning(n_qs,q,esp,new_nqs,new_q)  !choix des nouveaux q

c esp est désormais inutile
	 DEALLOCATE(esp)

c	 PRINT*,'anciens q',n_qs,new_nqs ; WRITE(*,2000)q
c	 PRINT*,'esp' ; WRITE(*,2000)esp ; PRINT*,'nouveaux q'
c	 WRITE(*,2000)new_q ; CALL pause('après zoning')

c dans new_b on se place aux new_nqs points new_q
c on calcule la nouvelle fonction de répartition
	 ALLOCATE(new_b(ne,new_nqs),new_esp(1,new_nqs))
	 DO i=1,new_nqs
	  CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,new_q(i),lq,f,dfdq)	  
	  IF(no_croiss)PRINT*,'Problème localisé en 2 dans lim_zc'
	  new_b(1:ne,i)=f(1:ne)
c	  WRITE(*,2000)new_b(1:ne,i)
	  new_esp(1,i)=ctep*f(1)+ctet*f(2)+cter*f(3)+ctel*f(4)+ctem*f(5)
c	  WRITE(*,2000)q(i),(f(j),j=1,ne),new_esp(1,i)
	  IF(i > 1)THEN
	   DO WHILE(new_esp(1,i) < new_esp(1,i-1))
	    new_esp(1,i)=new_esp(1,i)+0.01d0
	   ENDDO
	  ENDIF
	 ENDDO
c	 CALL pause('new_b')

c spline d'ordre 2 sur les nouveaux q sur new_nqs points
	 new_q=(/ (i, i=1,new_nqs) /)
	 ALLOCATE(new_qt(new_nqs+ord_qs))
	 CALL bsp1dn(1,new_esp,new_q,new_qt,new_nqs,2,new_knot,.FALSE.,
	1 new_q(1),lq,f,dfdq)  
	 IF(no_croiss)THEN
	  PRINT*,'Arrêt 1 dans lim_zc' ; STOP
	 ENDIF

c dérivée de la fonction de répartition dans new_b(6,:)
	 DO i=1,new_nqs 
	 CALL bsp1dn(1,new_esp,new_q,new_qt,new_nqs,2,new_knot,.TRUE.,
	1 new_q(i),lq,f,dfdq)
	 IF(no_croiss)PRINT*,'Problème localisé en 3 dans lim_zc'	
	 new_b(6,i)=dfdq(1)
	 ENDDO
c	 WRITE(*,2000)new_b(6,:) ; CALL pause('new_b(6,:)')

c on spline new_b
	 new_q=(/ (i, i=1,new_nqs) /)
	 CALL bsp1dn(ne,new_b,new_q,new_qt,new_nqs,ord_qs,new_knot,.FALSE.,
	1 new_q(1),lq,f,dfdq)
	 IF(no_croiss)THEN
	  PRINT*,'Arrêt 2 dans lim_zc' ; STOP
	 ENDIF

c	 DO i=1,n_qs
c	  CALL bsp1dn(ne,new_b,new_q,new_qt,new_nqs,ord_qs,new_knot,.TRUE.,
c	1  new_q(i),lq,f,dfdq)
c	  WRITE(*,2000)f(1:ne)
c	 ENDDO
c	 CALL pause('new_b')

c commentaires
	 IF(n_qs /=new_nqs)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(2,1015)n_qs,new_nqs ; WRITE(*,1015)n_qs,new_nqs 
1015	   FORMAT('Change of the number of shells :',i4,'-->',i4,/)
	  CASE DEFAULT
	   WRITE(2,15)n_qs,new_nqs  ; WRITE(*,15)n_qs,new_nqs 
15	   FORMAT('Modification du nombre de couches :',i4,'-->',i4,/)
	  END SELECT
	 ENDIF

c dimension de la base splines pour l'intégration équi.quasi.stat. 
	 n_qs=new_nqs ; dim_qs=(n_qs-1)*m_qs+r_qs

c vecteur nodal pour l'intégration sur n_qs points
	 DEALLOCATE(q,qt)
	 ALLOCATE(q(n_qs),qt(dim_qs+ord_qs))
	 q=(/ (i, i=1,n_qs) /)   
	 CALL noedif(q,n_qs,m_qs,r_qs,qt,knot)
	 IF(no_croiss)THEN
	  PRINT*,'Arrêt 3 dans lim_zc' ; STOP
	 ENDIF
 
c	 PRINT*,ord_qs,knot,dim_qs
c	 WRITE(*,2000)qt(1:knot) ; CALL pause('noedif')

c on place la solution initiale dans la base de q, bp_t ==> bp
	 DEALLOCATE(bp)
	 ALLOCATE(bp(ne,dim_qs))
	 CALL newspl(ne,new_q,new_qt,new_knot,ord_qs,qt,knot,ord_qs,new_b,bp) 
c	 CALL pause('newspl')

c new_b, new_q, new_qt, new_esp inutiles, r2 et m23 sont à modifier
	 DEALLOCATE(new_b,new_q,new_qt,new_esp,r2,m23)

c adaptation de fac, d_grad, r2 et m23 à la nouvelle répartition
	 ALLOCATE(fac(n_qs-1),d_grad(n_qs),r2(n_qs),m23(n_qs))

c tabulation de r2 et m23 qui accélèrent la recherche dans le ssp. inter
	 DO i=1,n_qs 
c	  CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),lc,f,dfdq)
	  r2(i)=bp(3,m_qs*(i-1)+1) ; m23(i)=bp(5,m_qs*(i-1)+1)
	 ENDDO
 
c test pour changement de grille fixe pour la composition chimique
c il y a changement de grille fixe quand
c new_nqs diffère de nc_fixe de plus dn_fixe (5% environ fixé dans cesam)
	 IF(grille_fixe)THEN
	  IF(ABS(REAL(new_nqs-nc_fixe,dp)/REAL(new_nqs,dp)) > dn_fixe)THEN
	   WRITE(*,30)nc_fixe,n_qs ; WRITE(2,30)nc_fixe,n_qs	  
30	   FORMAT(/,'changement de grille fixe: nc_fixe=',i4,' ==> ',i4,/)
	   DEALLOCATE(mc_fixe) ; ALLOCATE(mc_fixe(n_qs))	  
	   nc_fixe=n_qs ; mc_fixe=m23
	   IF(.NOT.en_masse)mc_fixe=mc_fixe**(2.d0/3.d0)	   
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(usl_static,1050)
1050	    FORMAT(/,'Change of the fixed grid')		
	   CASE DEFAULT
	    WRITE(usl_static,50)
50	    FORMAT(/,'Changement de grille fixe')
	   END SELECT	   
	  ENDIF
	 ENDIF
	ENDIF I1
	
c au centre R, L, M =0	
	bp(3:5,1)=0.d0

c------------------recherche des limites ZR/ZC----------------------------

c recherche les limites zones convectives / zones radiatives
c en fonction de q, détermination des facteurs de répartition
c overshooting inférieur et supérieur : extension de la zone mélangée
c ou overshooting inférieur PC de JPZ avec saut du gradient
c modèle totalement convectif: lim=1, jlim(1)=n, lconv(1)=.FALSE.
c modèle totalement radiatif: lim=0, jlim(i)=-100, lconv(i)=.FALSE.
c lconv(i)=..TRUE. en r_zc(i) on passe de ZC à ZR

c	CALL pause('début recherche limites ZR/ZC')

	SELECT CASE(langue)
	CASE('english')
	 WRITE(usl_static,1010)
1010	 FORMAT(/,'----- Search of the limits RZ/CZ (begin) --------')		
	CASE DEFAULT
	 WRITE(usl_static,10)
10	 FORMAT(/,'------ Recherche des limites ZR/ZC (début) --------')
	END SELECT

c initialisations
	r_ov=-HUGE(1.d0) ;  r_zc=-HUGE(1.d0) ; m_zc=-HUGE(1.d0)	
	lim=0 ; jlim=-100 ; lconv=.TRUE.
	tot_conv=.FALSE. ; tot_rad =.FALSE. ; ovsht=MAX(ovshts,ovshti) > 0.d0

c calcul de grad-grad
c	PRINT*,n_ch,m_ch,knotc,n_qs,nchim ; WRITE(*,2000)mc
c	WRITE(*,2000)mct ; PRINT*,mstar ; CALL pause('début lim_zc')

	DO i=1,n_qs	!calcul de grad-grad
	 CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),lq,f,dfdq)
	 IF(no_croiss)PRINT*,'Problème localisé en 4 dans lim_zc'	 
c	 WRITE(*,2000)f
	 pt=EXP(f(1))
	 dlpp=1.d0
	 IF(pturb)THEN  	!avec pression turbulente
	  p=EXP(f(Ipg))
	  IF(der)dlpp=dfdq(Ipg)/dfdq(1)   !dlpp=dln Pgaz/dln Ptot
	 ELSE 			!sans pression turbulente
	  p=pt
	 ENDIF
	 t=EXP(f(2))   
	 IF(en_masse)THEN
	  r=SQRT(ABS(f(3))) ; l=SQRT(ABS(f(4)))**3 ; m=SQRT(ABS(f(5)))**3
	 ELSE
	  r=ABS(f(3)) ; l=f(4) ; m=ABS(f(5)) ; f(5)=m**(2.d0/3.d0)
	 ENDIF
	 CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1 knotc,.TRUE.,MAX(mc(1),MIN(f(5),mc(n_ch))),lc,comp,dcomp)
	 IF(no_croiss)THEN
	  PRINT*,'Problème localisé en 5 dans lim_zc',i,knotc,n_ch
	  WRITE(*,2000)q(i),mc(1),f(5),mc(n_ch),mct(knotc) ; PAUSE'en 5'
	 ENDIF

c vitesse angulaire
	 SELECT CASE(Krot)
	 CASE(0,1,2)
	  w=wrot
	 CASE(3,4,5)
	  CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,
	1 MAX(mrot(1),MIN(f(5),mrot(n_rot))),lc,frot,dfrot)
	  IF(no_croiss)PRINT*,'Problème localisé en 6 dans lim_zc'	
	  w=frot(1)
	 END SELECT 
	 d_grad(i)=dgrad(pt,p,t,dlpp,comp,m,l,r,dcomp,w)
	ENDDO   !i
c	WRITE(*,2000)d_grad ; PAUSE'd_grad'

c r à l'extérieur
	rn=r

c ZC avec zero exact en 1
	grp=d_grad(1)
	IF(grp == 0.d0 .AND. d_grad(2) >= 0.d0)grp=-1
	lim=0  !initialisation 
	B1: DO j=2,n_qs    !recherche des limites ZR/ZC
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
	  ENDIF
	  lconv(lim)=grp < 0.d0 !.TRUE. passage ZR-->ZC
	  rap(lim)=ABS(grp/gr)
	  rap(lim)=rap(lim)/(1.d0+rap(lim))
	  jlim(lim)=j  !limite ZR/ZC entre j et j-1
	  IF(lconv(lim))THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(usl_static,1003)INT(rap(lim)*100.d0),j-1,j,lim
1003	    FORMAT('a CZ begins at',i3,'% between the shells #',
	1   i5,' and #',i5,', limit RZ/CZ #',i3)	  
	   CASE DEFAULT		  
	    WRITE(usl_static,3)INT(rap(lim)*100.d0),j-1,j,lim
3	    FORMAT('début ZC localisé à',i3,'% entre les couches #',
	1   i5,' et #',i5,', limite ZR/ZC #',i3)	
	   END SELECT	
	  ELSE
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(usl_static,1004)INT(rap(lim)*100.),j-1,j,lim
1004	    FORMAT('a CZ ends at ',i3,'% between the shells #',
	1   i5,' and #',i5,', limit RZ/CZ #',i3)	   
	   CASE DEFAULT	  
	    WRITE(usl_static,4)INT(rap(lim)*100.),j-1,j,lim
4	    FORMAT('fin de ZC localisée à',i3,'% entre les couches #',
	1   i5,' et #',i5,', limite ZR/ZC #',i3)
	   END SELECT
	  ENDIF
	 ENDIF
	 grp=gr
	ENDDO B1  !j

c	PRINT*,lim,jlim,lconv
c	PAUSE'sortie du tri, lim, jlim, lconv'
c	WRITE(*,2000)rap(1:lim) ; CALL pause('tri')

c indices des points de grille
	ALLOCATE(new(n_qs)) ; new=q

c il y a une limite ZR/ZC entre jlim(i)-1 et jlim(i), i=1,lim
c on détermine l'incide new(jlim(i)) à la limite ZR/ZC
c on fixe la limite soit à droite soit à gauche
c jlim(i) : indice le plus proche de la limite
	DO i=1,lim !j_new: indice de la lim. ZR/ZC
	 j_new=new(jlim(i)-1)+rap(i)*(new(jlim(i))-new(jlim(i)-1))
c	 PRINT*,i,j_new,new(jlim(i)-1),rap(i),new(jlim(i))
	 jlim(i)=NINT(j_new)    
	 new(jlim(i))=j_new!nouvel indice pour les limites ZR/ZC
c	 PRINT*,i,j_new,jlim(i)
	ENDDO

c	PRINT*,'indice pour les limites ZR/ZC' ; PRINT*,lim
c	PRINT*,(new(jlim(i)),i=1,lim) ; PRINT*,(jlim(i),i=1,lim)    

c suppression des limites trop proches
c no_supp=.FALSE. : suppression de limites séparées par moins de
c ncouches_zrc couches, pas d'affinement de la position des limites ZR/ZC
	no_supp=.TRUE.
	i=1
	DO WHILE (i <= lim)
	 IF(jlim(i) <= ncouches_zrc)THEN  !limites sur premières couches
	  WRITE(usl_static,*)'suppression d''une limite ZR/ZC trop centrale #',i
	  lim=lim-1 
	  DO j=1,lim
	   jlim(j)=jlim(j+1) ; new(jlim(j))=new(jlim(j+1))
	   lconv(j)=lconv(j+1)
	  ENDDO	  
c	  PRINT*,'lim,jlim',lim,jlim
	  
	 ELSE
	  i=i+1
	 ENDIF
	ENDDO
c	PRINT*,'limites sur premières couches',lim
c	PRINT*,new(jlim(1:lim)) ; PRINT*,jlim(1:lim) ; PAUSE'lim'   
	i=2
	DO WHILE (i <= lim)
	 IF(jlim(i)-jlim(i-1) <= ncouches_zrc)THEN !limites trop proches
	  no_supp=no_supp .AND. jlim(i)-jlim(i-1) >= 2 
	  IF(i == lim .AND. lconv(i-1))THEN
	   im=1 ; lim=lim-im
	   WRITE(usl_static,*)'limites trop proches, suppression limite ',i
	  ELSE
	   im=2 ; lim=lim-im    !suppression de i-1 et i
	  WRITE(usl_static,*)'limites trop proches, suppression limites',i-1,i
	  ENDIF
	  DO j=i-1,lim
	   jlim(j)=jlim(j+im) ; new(jlim(j))=new(jlim(j+im))
	   lconv(j)=lconv(j+im)
	  ENDDO 	!no_supp=.FALSE. : pas d'affinement des limites
	 ELSE
	  i=i+1
	 ENDIF !pour la suppression de limites trop rapprochées
	ENDDO
c	PRINT*,'limites trop rapprochées'
c	PRINT*,lim
c	PRINT*,(new(jlim(i)),i=1,lim)
c	PRINT*,lim,jlim
c	PAUSE'lim,jlim' 

c suppression d'une limite trop externe
c	IF(lim > 0 .AND. jlim(lim) >= n_qs-ncouches_zrc)THEN
c	 WRITE(usl_static,*)'suppression d''une limite ZR/ZC trop externe #',lim
c	 lim=lim-1
c	ENDIF

c nature complètement convectif/radiatif du modèle
c si lim=0 et que la différence des gradients est positive sur la
c couche externe le modèle est complètement convectif
	IF(lim == 0)THEN	
	 IF(d_grad(n_qs) >= 0.d0)THEN
	  tot_conv=.TRUE. ; lim=1 ; jlim(1)=n_qs
	  lconv(1)=.FALSE. ; m_zc(1)=mstar ; r_zc(1)=rn	 
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(usl_static,1040)
1040	   FORMAT('Fully convective model',/)	   
	  CASE DEFAULT	  
	   WRITE(usl_static,40)
40	   FORMAT('Modèle totalement convectif',/)
	  END SELECT	 	  
	 ELSE
	  tot_rad=.TRUE. ; lim=0 ; jlim=-100 ; lconv=.FALSE. ; m_zc=-100.d0 
	  SELECT CASE(langue)
	  CASE('english')	   	 
	   WRITE(usl_static,1041)
1041	   FORMAT('WARNING, Unphysical fully radiative model',/)	   
	  CASE DEFAULT	  
	   WRITE(usl_static,41)
41	   FORMAT('ATTENTION, Cas non physique, modèle totalement radiatif',/)
	  END SELECT
	 ENDIF
	 
c----il y a des limites ZR/ZC	 
	ELSE

c--------on affine, par dichotomie, la position des limites retenues
c pour chaque limite, encadrement entre q=j1 et q=j2
	 IF(no_supp)THEN
	  DO i=1,lim
	  
c j1 --> dgrad(1)	  
	   j1=new(jlim(i))
	   CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,j1,lq,f,dfdq)
	   IF(no_croiss)PRINT*,'Problème localisé en 7 dans lim_zc'	   
c	   WRITE(*,2000)f
	   pt=EXP(f(1)) ; dlpp=1.d0   
	   IF(pturb)THEN   		!avec pression turbulente 
	    p=EXP(f(Ipg))
	    IF(der)dlpp=dfdq(Ipg)/dfdq(1) !dlpp=dln Pgaz/dln Ptot
	   ELSE    			!sans pression turbulente
	    p=pt
	   ENDIF
	   t=EXP(f(2))
	   IF(en_masse)THEN
	    r=SQRT(ABS(f(3))) ; l=SQRT(ABS(f(4)))**3 ; m=SQRT(ABS(f(5)))**3
	   ELSE
	    r=ABS(f(3)) ; l=f(4) ; m=ABS(f(5)) ; f(5)=m**(2.d0/3.d0)  
	   ENDIF
	   CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1   knotc,.TRUE.,MAX(mc(1),MIN(f(5),mc(n_ch))),lc,comp,dcomp)
	   IF(no_croiss)PRINT*,'Problème localisé en 8 dans lim_zc'	

c vitesse angulaire
	   SELECT CASE(Krot)
	   CASE(0,1,2)
	    w=wrot
	   CASE(3,4,5)
	    CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,
	1   MAX(mrot(1),MIN(f(5),mrot(n_rot))),lc,frot,dfrot)
	    IF(no_croiss)PRINT*,'Problème localisé en 9 dans lim_zc'	
	    w=frot(1)
	   END SELECT
	   d_grad(1)=dgrad(pt,p,t,dlpp,comp,m,l,r,dcomp,w)

c j2 en jlim(i) ---> dgrad(2)
	   j2=jlim(i)
	   CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,j2,lq,f,dfdq)
	   IF(no_croiss)PRINT*,'Problème localisé en 10 dans lim_zc'
c	   WRITE(*,2000)f
	   pt=EXP(f(1))
	   dlpp=1.d0 
	   IF(pturb)THEN   		!avec pression turbulente
	    p=EXP(f(Ipg))
	    IF(der) dlpp=dfdq(Ipg)/dfdq(1)    !dlpp=dln Pgaz/dln Ptot
	   ELSE    			!sans pression turbulente
	    p=pt
	   ENDIF
	   t=EXP(f(2))
	   IF(en_masse)THEN
	    r=SQRT(ABS(f(3))) ; l=SQRT(ABS(f(4)))**3 ; m=SQRT(ABS(f(5)))**3
	   ELSE
	    r=ABS(f(3)) ; l=f(4) ; m=ABS(f(5)) ; f(5)=m**(2.d0/3.d0)  
	   ENDIF
	   CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1  knotc,.TRUE.,MAX(mc(1),MIN(f(5),mc(n_ch))),lc,comp,dcomp)
	   IF(no_croiss)PRINT*,'Problème localisé en 11 dans lim_zc'	
   
c vitesse angulaire
	   SELECT CASE(Krot)
	   CASE(0,1,2)
	    w=wrot
	   CASE(3,4,5)
	    CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,
	1   MAX(mrot(1),MIN(f(5),mrot(n_rot))),lc,frot,dfrot)
	    IF(no_croiss)PRINT*,'Problème localisé en 12 dans lim_zc'	
	    w=frot(1)
	   END SELECT
	   d_grad(2)=dgrad(pt,p,t,dlpp,comp,m,l,r,dcomp,w)
	   
c pas d'inversion entre j1 et j2=jlim(i) --> j2 en jlim(i)-1	       
	   IF(d_grad(1)*d_grad(2) > 0.d0)THEN
	    j2=jlim(i)-1
	    CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,j2,lq,f,dfdq)
	    IF(no_croiss)PRINT*,'Problème localisé en 13 dans lim_zc'	    
c	    WRITE(*,2000)f
	    pt=EXP(f(1)) ; dlpp=1.d0
	    IF(pturb)THEN  		!avec pression turbulente
 	     p=EXP(f(Ipg))
 	     IF(der)dlpp=dfdq(Ipg)/dfdq(1)    !dlpp=dln Pgaz/dln Ptot
	    ELSE  			!sans pression turbulente
	     p=pt
	    ENDIF
	    t=EXP(f(2))
	    IF(en_masse)THEN
	     r=SQRT(ABS(f(3))) ; l=SQRT(ABS(f(4)))**3 ; m=SQRT(ABS(f(5)))**3
	    ELSE
	     r=ABS(f(3)) ; l=f(4) ; m=ABS(f(5)) ; f(5)=m**(2.d0/3.d0)
	    ENDIF
	    CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1   knotc,.TRUE.,MAX(mc(1),MIN(f(5),mc(n_ch))),lc,comp,dcomp)
	    IF(no_croiss)PRINT*,'Problème localisé en 14 dans lim_zc'	

c vitesse angulaire
	    SELECT CASE(Krot)
	    CASE(0,1,2)
	     w=wrot
	    CASE(3,4,5)
	     CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,
	1    MAX(mrot(1),MIN(f(5),mrot(n_rot))),lc,frot,dfrot)
	     IF(no_croiss)PRINT*,'Problème localisé en 15 dans lim_zc'	
	     w=frot(1)
	    END SELECT
	    d_grad(2)=dgrad(pt,p,t,dlpp,comp,m,l,r,dcomp,w)
	    
c pas d'inversion entre j1 et j2=jlim(i)-1 ---> j2 en jlim(i)+1	    
	    IF(d_grad(1)*d_grad(2) > 0.d0)THEN  
	     j2=jlim(i)+1
	     CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,j2,lq,f,dfdq)
	     IF(no_croiss)PRINT*,'Problème localisé en 16 dans lim_zc'	     
c	     WRITE(*,2000)f
	     pt=EXP(f(1)) ; dlpp=1.d0 
	     IF(pturb)THEN 		!avec pression turbulente
	      p=EXP(f(Ipg))
	      IF(der) dlpp=dfdq(Ipg)/dfdq(1)  !dlpp=dln Pgaz/dln Ptot
	     ELSE 			!sans pression turbulente
	      p=pt
	     ENDIF
	     t=EXP(f(2))
	     IF(en_masse)THEN
	      r=SQRT(ABS(f(3))) ; l=SQRT(ABS(f(4)))**3
	      m=SQRT(ABS(f(5)))**3
	     ELSE
	      r=ABS(f(3)) ; l=f(4) ; m=ABS(f(5)) ; f(5)=m**(2.d0/3.d0)
	     ENDIF
	     CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1    knotc,.TRUE.,MAX(mc(1),MIN(f(5),mc(n_ch))),lc,comp,dcomp)
	     IF(no_croiss)PRINT*,'Problème localisé en 17 dans lim_zc'	

c vitesse angulaire
	     SELECT CASE(Krot)
	     CASE(0,1,2)
	      w=wrot
	     CASE(3,4,5)
	      CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,
	1     MAX(mrot(1),MIN(f(5),mrot(n_rot))),lc,frot,dfrot)
	      IF(no_croiss)PRINT*,'Problème localisé en 18 dans lim_zc'	
	      w=frot(1)
	     END SELECT
	     d_grad(2)=dgrad(pt,p,t,dlpp,comp,m,l,r,dcomp,w)
	     IF(d_grad(1)*d_grad(2) > 0.d0)THEN
	      PRINT*,'lim_zc : erreur dans la dichotomie, limite ZR/ZC',i
	      PRINT*,'j1, jlim(i)',j1,jlim(i) ; PRINT*,'on arrète'
c	      CALL pause('on arrète')
	      STOP
	     ENDIF  !jlim(i)+1
	    ENDIF   !jlim(i)-1
	   ENDIF    !on a encadré la limite entre j1 et j2
	   
c limite encadrée entre j1 et j2, dichotomie avec j3	   
	   compt=0
	   j3=(j1+j2)/2.d0
	   DO WHILE(ABS(j1-j2) > loc_zc)
	    CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,j3,lq,f,dfdq)
	    IF(no_croiss)PRINT*,'Problème localisé en 19 dans lim_zc'	    
c	    WRITE(*,2000)f
	    pt=EXP(f(1)) ; dlpp=1.d0
	    IF(pturb)THEN  		!avec pression turbulente
	     p=EXP(f(Ipg))
	     IF(der) dlpp=dfdq(Ipg)/dfdq(1)   !dlpp=dln Pgaz/dln Ptot
	    ELSE  			!sans pression turbulente
	     p=pt
	    ENDIF
	    t=EXP(f(2))
	    IF(en_masse)THEN
	     r=SQRT(ABS(f(3))) ; l=SQRT(ABS(f(4)))**3
	     m=SQRT(ABS(f(5)))**3
	    ELSE
	     r=ABS(f(3)) ; l=f(4) ; m=ABS(f(5)) ; f(5)=m**(2.d0/3.d0)    
	    ENDIF
	    CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1   knotc,.TRUE.,MAX(mc(1),MIN(f(5),mc(n_ch))),lc,comp,dcomp)
	    IF(no_croiss)PRINT*,'Problème localisé en 20 dans lim_zc'	

c vitesse angulaire
	    SELECT CASE(Krot)
	    CASE(0,1,2)
	     w=wrot
	    CASE(3,4,5)
	     CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,
	1    MAX(mrot(1),MIN(f(5),mrot(n_rot))),lc,frot,dfrot)	
	     IF(no_croiss)PRINT*,'Problème localisé en 21 dans lim_zc'	     
	     w=frot(1)
	    END SELECT
	    d_grad(3)=dgrad(pt,p,t,dlpp,comp,m,l,r,dcomp,w)    
c	    PRINT*,j1,j2,j3,compt ; WRITE(*,2000)d_grad

c limite entre j2 et j3
	    IF(d_grad(2)*d_grad(3) < 0.d0)THEN
	     j1=j3 ; d_grad(1)=d_grad(3)
	     
c limite entre j1 et j3	     
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
	   ENDDO   !WHILE
	   jlim(i)=NINT(j3) ; new(jlim(i))=j3 
	  ENDDO			!i

c---------écritures-------------------------------

	  PRINT*
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(usl_static,*)'Limits RZ/CZ after the refinement of their locations :'	  
	  CASE DEFAULT
	   WRITE(usl_static,*)'Limites ZR/ZC après affinement des localisations :'
	  END SELECT	  
	  
	  DO i=1,lim  !masses des limites
	   CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,REAL(jlim(i),dp),
	1  lq,f,dfdq)
	   IF(no_croiss)PRINT*,'Problème localisé en 22 dans lim_zc'
	   IF(en_masse)THEN
	    r_zc(i)=SQRT(ABS(f(3))) ; m_zc(i)=SQRT(ABS(f(5)))**3
	   ELSE
	    r_zc(i)=ABS(f(3)) ; m_zc(i)=ABS(f(5))
	   ENDIF	   

	   IF(lconv(i))THEN !début de ZC
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(usl_static,1024)i,INT((new(jlim(i))-INT(new(jlim(i))))*100.),
	1    INT(new(jlim(i))),INT(new(jlim(i)))+1,r_zc(i)/rn,
	2    m_zc(i)/mstar,1.d0-m_zc(i)/mstar
1024	     FORMAT('for the limit RZ/CZ #',i4,' at ',i3,'% between shell',
	1    i5,' and',i5,', beginning of CZ',/,'R_zc/Rtot=',es10.3,
	2    ', M_zc/Mstar=',es10.3,', 1-M_zc/Mstar=',es10.3)	    
	    CASE DEFAULT	   
	     WRITE(usl_static,24)i,INT((new(jlim(i))-INT(new(jlim(i))))*100.),
	1    INT(new(jlim(i))),INT(new(jlim(i)))+1,r_zc(i)/rn,
	2    m_zc(i)/mstar,1.d0-m_zc(i)/mstar
24	     FORMAT('pour limite ZR/ZC #',i4,' à ',i3,'% entre couches',i5,
	1    ' et',i5,', début de ZC',/,'R_zc/Rtot=',es10.3,
	2    ', M_zc/Mstar=',es10.3,', 1-M_zc/Mstar=',es10.3)
	    END SELECT
	   ELSE
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(usl_static,1025)i,INT((new(jlim(i))-INT(new(jlim(i))))*100.),
	1    INT(new(jlim(i))),INT(new(jlim(i)))+1,r_zc(i)/rn,
	2    m_zc(i)/mstar,1.d0-m_zc(i)/mstar
1025	     FORMAT('for the limit RZ/CZ #',i4,' at ',i3,'% between shell',
	1    i5,' and #',i5,', end of CZ',/,'R_zc/Rtot=',es10.3,
	2    ', M_zc/Mstar=',es10.3,', 1-M_zc/Mstar=',es10.3)
	    CASE DEFAULT	   	   
	     WRITE(usl_static,25)i,INT((new(jlim(i))-INT(new(jlim(i))))*100.),
	1    INT(new(jlim(i))),INT(new(jlim(i)))+1,r_zc(i)/rn,
	2    m_zc(i)/mstar,1.d0-m_zc(i)/mstar
25	     FORMAT('pour limite ZR/ZC #',i4,' à ',i3,'% entre couches #',
	1    i5,' et #',i5,', fin de ZC',/,'R_zc/Rtot=',es10.3,
	2    ', M_zc/Mstar=',es10.3,', 1-M_zc/Mstar=',es10.3)
	    END SELECT	
	   ENDIF
	  ENDDO 	!i

	 ELSE		!no_supp=.FALSE
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(usl_static,*)'no refinement of the localization of the limits CZ/RZ'
	   WRITE(usl_static,*)'because too closed limits have been suppressed'    
	   WRITE(usl_static,*)'limits CZ/RZ used:'	 
	  CASE DEFAULT
	   WRITE(usl_static,*)'Pas d''affinement de la localisation des limites ZR/ZC'
	   WRITE(usl_static,*)'à cause de la suppression de limites ZR/ZC trop proches'    
	   WRITE(usl_static,*)'Limites ZR/ZC retenues:'
	  END SELECT
	    
	  DO i=1,lim   !masses des limites
	   CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,new(jlim(i)),
	1    lq,f,dfdq)
c	   j3=REAL(jlim(i),dp)
c	   CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,j3,lq,f,dfdq)	
	   IF(no_croiss)PRINT*,'Problème localisé en 23 dans lim_zc'	
	   IF(en_masse)THEN
	    r_zc(i)=SQRT(ABS(f(3))) ; m_zc(i)=SQRT(ABS(f(5)))**3
	   ELSE
	    r_zc(i)=ABS(f(3)) ; m_zc(i)=ABS(f(5))
	   ENDIF
c	   PRINT*,'masse r_zc,m_zc',i,jlim(i),new(jlim(i))

	   IF(lconv(i))THEN !début de ZC
	    WRITE(usl_static,24)i,INT((new(jlim(i))-INT(new(jlim(i))))*100.),
	1   INT(new(jlim(i))),INT(new(jlim(i)))+1,r_zc(i)/rn,
	2   m_zc(i)/mstar,1.d0-m_zc(i)/mstar
	   ELSE
	    WRITE(usl_static,25)i,INT((new(jlim(i))-INT(new(jlim(i))))*100.),
	1   INT(new(jlim(i))),INT(new(jlim(i)))+1,r_zc(i)/rn,
	2   m_zc(i)/mstar,1.d0-m_zc(i)/mstar  
	   ENDIF
	  ENDDO 	!i
	 ENDIF  !no_supp

c	 PRINT*,lconv(1:lim) ; CALL pause('lconv')

c----------overshoot-----------------

c calcul de r_zc et de r_ov si overshooting ou PC
	 IF(ovsht)THEN !s'il y a overshooting ou PC
	  DO i=1,lim
c	   PRINT*,'limite ',i,' ',lconv(i)
	   r_zc(i)=-100d0 ; r_ov(i)=-100d0
	   CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,new(jlim(i)),
	1  lq,f,dfdq)
c	   j3=REAL(jlim(i),dp)
c	   CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,j3,lq,f,dfdq)
	   IF(no_croiss)PRINT*,'Problème localisé en 24 dans lim_zc'
	   pt=EXP(f(1)) ; dlpp=1.d0
	   IF(pturb)THEN   			!avec pression turbulente
	    p=EXP(f(Ipg))
	    IF(der)dlpp=dfdq(Ipg)/dfdq(1) !dlpp=dln Pgaz/dln Ptot
	   ELSE    				!sans pression turbulente
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
	   CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1  knotc,.TRUE.,MAX(mc(1),MIN(f(5),mc(n_ch))),lc,comp,dcomp)
	   IF(no_croiss)PRINT*,'Problème localisé en 25 dans lim_zc'
c	   PRINT*,'m_zc,r_zc,new(jlim(i))',lconv(i),i
c	   WRITE(*,2000)m_zc(i),r_zc(i),new(jlim(i))
c	   WRITE(*,2000)p,t,r_zc(i),l,m_zc(i),comp(1)

c vitesse angulaire
	   SELECT CASE(Krot)
	   CASE(0,1,2)
	    w=wrot
	   CASE(3,4,5)
	    CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,
	1   MAX(mrot(1),MIN(f(5),mrot(n_rot))),lc,frot,dfrot)
	    IF(no_croiss)PRINT*,'Problème localisé en 26 dans lim_zc'	
	    w=frot(1)
	   END SELECT

	   CALL thermo(pt,p,t,m_zc(i),l,r_zc(i),dlpp,comp,dcomp,
	1  ro,drop,drot,drox,u,dup,dut,dux,
	2  grad,dgradpt,dgradp,dgradt,dgradx,dgradm,dgradl,dgradr,dgradlpp,
	3  gam,dgampt,dgamp,dgamt,dgamx,dgamm,dgaml,dgamr,dgamlpp,
	4  epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
	5  delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	6  gradad,dgradadp,dgradadt,dgradadx,
	7  hp,dhppt,dhpp,dhpt,dhpx,dhpr,dhpm,grad_mu,
	8  gradrad,alfa,beta,gamma1,radiatif)

c overshoot: ovshti(s)*Hp ou ovshti*kip de JPZ si T > 5d5 à la limite
	   IF(t > 5.d5)THEN
	    IF(lconv(i))THEN !début de ZC, ext. vers le bas
	     IF(ovshti > 0.d0)THEN
	      IF(jpz)THEN
	       kip=(3.d0-t*(drot/ro+dkapt/kap))*gradad-p*(drop/ro+dkapp/kap)
	       r_ov(i)=MAX(r_zc(i)-hp*ovshti/rsol/kip,0.d0)
c	       PRINT*,'i/r_zc(i),hp*ovshti/rsol/kip,hp,ovshti,rsol,kip',i
c	       WRITE(*,2000)r_zc(i),hp*ovshti/rsol/kip,hp,ovshti,rsol,kip
c	       WRITE(*,2000)r_ov(i)
	      ELSE !ext. vers le bas: 0.<r_ov<r_zc
	       r_ov(i)=MAX(r_zc(i)-hp*ovshti/rsol,0.d0)
	      ENDIF !jpz
c	      PRINT*,'overshoot inférieur'
c	      WRITE(*,2000)r_ov(i),r_zc(i),hp
	     ENDIF !ovshti
	    ELSE    !fin de ZC overshoot supérieur: rext>r_ov>r_zc 
	     IF(ovshts > 0.d0)THEN !il y a ovsht sup.
	      IF(jpz)THEN
	       r_ov(i)=r_zc(i)*(1.d0+ovshts) !overshoot sup.= R_noyauXovshts
	      ELSE
	       r_ov(i)=MIN(r_zc(i)+hp*ovshts/rsol,r_zc(i)*(1.d0+ovshts))
	      ENDIF !jpz
	      r_ov(i)=MIN(r_ov(i),rn)
c	      PRINT*,'ovshts' ; WRITE(*,2000)r_ov(i),r_zc(i),hp,rn
	     ENDIF !ovshts
	    ENDIF  !lconv
c	    CALL pause('après if sur lconv')

c masse d'overshooting ou de pénétration, ATTENTION 'r2 ' et non 'r2'
	    IF(r_ov(i) >= 0.d0)THEN !r_ov=-100: no overshoot en r_zc(i)
	     IF(en_masse)THEN
	      CALL inter('r2 ',bp,q,qt,n_qs,knot,r_ov(i)**2,f,dfdq,r2,m23) 
	      m_zc(i)=SQRT(ABS(f(5)))**3
	     ELSE
	      CALL inter('r2 ',bp,q,qt,n_qs,knot,r_ov(i),f,dfdq,r2,m23)
	      m_zc(i)=ABS(f(5)) 
	     ENDIF
	     j_new=f(6)!astuce pour passer l'indice
c	     WRITE(*,2000)f(1:ne) 

c l'ancien new(jlim(i)) est désormais un point ordinaire
c qui redevient l'indice entier jlim(i)
c la limite est désormais en jlim(i)= + proche entier de j_new
c	     PRINT*,'jlim(i),new(jlim(i)),j_new',jlim(i),new(jlim(i)),j_new
	     jlim(i)=NINT(j_new)    
	     new(jlim(i))=j_new	!nouvel indice pour les limites ZR/ZC
	    
c	     PRINT*,'jlim(i),new(jlim(i)),j_new',jlim(i),new(jlim(i)),j_new
c	     PRINT*,'j_new/M_ov,R_ov',j_new ; WRITE(*,2000)m_zc(i),r_ov(i)
c	     CALL pause('après j_new')
	    ENDIF  !r_ov /= 0
	   ELSE
	    r_ov(i)=r_zc(i)
c	    CALL pause('t > 1.d5')
	   ENDIF   !t > 1.d5
	  ENDDO		!i sur lim

c il ne doit pas y avoir de chevauchement après les extensions
c	  PRINT*,lim,ovsht,(jlim(i),i=1,lim),(lconv(i),i=1,lim)
c	  WRITE(*,2000)(r_ov(j),j=1,lim)
c	  WRITE(*,2000)(r_zc(j),j=1,lim)
c	  WRITE(*,2000)(m_zc(j),j=1,lim)

	  IF(lconv(1))THEN  !début de ZC à la première limite
	   i=2
	  ELSE
	   i=1
	  ENDIF
	  DO WHILE (i < lim)
	   IF(jlim(i) >= jlim(i+1))THEN
	    DO j=i,lim-2    !chevauchement: suppression d'une ZR décalage
	     jlim(j)=jlim(j+2) ; m_zc(j)=m_zc(j+2) ; lconv(j)=lconv(j+2)
	     r_zc(j)=r_zc(j+2) ; r_ov(j)=r_ov(j+2)
	    ENDDO
	    lim=lim-2
	   ELSE
	    i=i+2
	   ENDIF
	  ENDDO		!WHILE
   
c à cause d'un overshoot, le modèle devient totalement convectif
	  IF(lim == 1 .AND. jlim(1) == 1 .AND. lconv(1))THEN
	   jlim(1)=n_qs ; lconv(1)=.FALSE. ; tot_conv=.TRUE. ; r_zc(1)=rn
	  ELSE
	   tot_conv=.FALSE.
    
c pas de limite ZR/ZC aux extrémités
	   IF(jlim(1) == 1)THEN			!en 1
	    DO i=1,lim-1
	     jlim(i)=jlim(i+1) ; lconv(i)=lconv(i+1) ; m_zc(i)=m_zc(i+1)
	     r_zc(i)=r_zc(i+1) ; r_ov(i)=r_ov(i+1)
	    ENDDO
	    lim=lim-1
	   ENDIF
	   IF(jlim(lim) == n_qs)lim=lim-1	 !en n
	  ENDIF
	  PRINT*
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(usl_static,*)'limits RZ/CZ after examination of overshooting :'
	  CASE DEFAULT   
	   WRITE(usl_static,*)'limites ZR/ZC après examen de l''overshooting :'
	  END SELECT

	  IF(tot_conv)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(usl_static,1033)
1033	    FORMAT('fully convective model',/)	  
	   CASE DEFAULT	  
	    WRITE(usl_static,33)
33	    FORMAT('modèle totalement convectif',/)
	   END SELECT
	  ELSE
	   DO j=1,lim
	    IF(lconv(j))THEN   !début de ZC
	     IF(ovshti > 0.d0)THEN
	      IF(jpz)THEN
	       SELECT CASE(langue)
	       CASE('english')
	        WRITE(usl_static,1031)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
1031	        FORMAT('Convective penetration of JPZ limit CZ/RZ #',i4,
	1       ', shell #',i5,/,3x,'radius CZ/Rtot=',es10.3,
	2       ', reduced radius/Rtot=',es10.3,/,
	3       'm/Mstar at the limit of the convective penetration=',
	4       es10.3,/)	      
	       CASE DEFAULT		      
	        WRITE(usl_static,31)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
31	        FORMAT('Pénétration Convective JPZ limite ZC/ZR #',i4,
	1       ', couche #',i5,/,3x,'rayon ZC/Rtot=',es10.3,
	2       ', rayon reduit/Rtot=',es10.3,/,
	3       'm/Mstar à la limite de la pénétration convective=',
	4       es10.3,/)
	       END SELECT
	      ELSE !jpz
	       SELECT CASE(langue)
	       CASE('english')
	        WRITE(usl_static,1002)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
1002	        FORMAT('negative overshooting  limit CZ/RZ #',i4,
	1       ', shell #',i5,/,'rayon ZC/Rtot=',es10.3,
	2       ', reduced raius/Rtot=',es10.3,/,
	3       'm/Mstar at the limit of overshooting=',es10.3,/)
	       CASE DEFAULT		      
	        WRITE(usl_static,2)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
2	        FORMAT('Overshooting inférieur limite ZC/ZR #',i4,
	1       ', couche #',i5,/,'rayon ZC/Rtot=',es10.3,
	2       ', rayon reduit/Rtot=',es10.3,/,
	3       'm/Mstar à la limite de l''overshooting=',es10.3,/)
	       END SELECT		
	      ENDIF!jpz
	     ELSE  !sans overshoot
	      SELECT CASE(langue)	     
	      CASE('english')
	       WRITE(usl_static,1011)j,jlim(j),r_zc(j)/rn,m_zc(j)/mstar
1011	       FORMAT('Limit without overshooting nor PC, limit RZ/CZ #:',
	1      i4,', shell #',i5,/,'radius CZ/Rtot=',es10.3,
	2      ', m/Mstar at the end of the CZ=',es10.3,/)	      
	      CASE DEFAULT		      
	       WRITE(usl_static,11)j,jlim(j),r_zc(j)/rn,m_zc(j)/mstar
11	       FORMAT('Limite sans overshooting ni PC, limite ZR/ZC #:',i4,
	1      ', couche #',i5,/,'rayon ZC/Rtot=',es10.3,
	2      ', m/Mstar à la fin de la ZC=',es10.3,/)
	      END SELECT	
	     ENDIF !ovshti
	    ELSE   !fin de ZC
	     IF(ovshts > 0.d0)THEN
	      IF(jpz)THEN
	       SELECT CASE(langue)	      
	       CASE('english')
	        WRITE(usl_static,1012)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
1012	        FORMAT('positive overshoot * R core, limit CZ/RZ #',i4,
	1       ', shell #',i5,/,3x,'radius CZ/Rtot=',es10.3,
	2       ', extented radius/Rtot=',es10.3,/,
	3       'm/Mstar at the limit of overshoot=',es10.3,/)	      
	       CASE DEFAULT		      
	        WRITE(usl_static,12)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
12	        FORMAT('Overshoot supérieur * R noyau, limite ZC/ZR #',i4,
	1       ', couche #',i5,/,3x,'rayon ZC/Rtot=',es10.3,
	2       ', rayon etendu/Rtot=',es10.3,/,
	3       'm/Mstar à la limite de l''overshoot=',es10.3,/)
	       END SELECT
	      ELSE !jpz	    
	       IF(r_ov(j) == r_zc(j))THEN
	        SELECT CASE(langue)	       
	        CASE('english')
	         WRITE(usl_static,1013)j,jlim(j),r_zc(j)/rn,m_zc(j)/mstar
1013	         FORMAT('Limit CZ/RZ without positive overshoot #',i4,
	1        ', shell #',i5,/,'because T < 5.d5, radius CZ/Rtot=',
	2	 es10.3,', m/Mstar on the limit =',es10.3,/)	      
	        CASE DEFAULT		      
	         WRITE(usl_static,13)j,jlim(j),r_zc(j)/rn,m_zc(j)/mstar
13	         FORMAT('Limite ZC/ZR sans overshooting supérieur #',i4,
	1        ', couche #',i5,/,'car T < 5.d5, rayon ZC/Rtot=',es10.3,
	2        ', m/Mstar à la limite =',es10.3,/)
	        END SELECT
	       ELSE
	        SELECT CASE(langue)	       
	        CASE('english')
	         WRITE(usl_static,1001)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
1001	         FORMAT('Overshooting positive * Hp, limit CZ/RZ #',i4,
	1        ', shell #',i5,/,'radius CZ/Rtot=',es10.3,
	2        ', extended radius/Rtot=',es10.3,/,
	3        'm/Mstar at the limit of overshooting=',es10.3,/)	      
	        CASE DEFAULT		      
	         WRITE(usl_static,1)j,jlim(j),r_zc(j)/rn,r_ov(j)/rn,m_zc(j)/mstar
1	         FORMAT('Overshooting supérieur * Hp, limite ZC/ZR #',i4,
	1        ', couche #',i5,/,'rayon ZC/Rtot=',es10.3,
	2        ', rayon etendu/Rtot=',es10.3,/,
	3        'm/Mstar à la limite de l''overshooting=',es10.3,/)
	        END SELECT
	       ENDIF
	      ENDIF!jpz
	     ELSE  !pas d'overshoot
	      SELECT CASE(langue)	     
	      CASE('english')
	       WRITE(usl_static,1021)j,jlim(j),r_zc(j)/rn,m_zc(j)/mstar
1021	       FORMAT('Limit without overshooting nor PC, limit RZ/CZ #',
	1      i4,', shell #',i5,/,'radius CZ/Rtot=',es10.3,
	2      ', m/Mstar at the beginning of the CZ=',es10.3,/)	      
	      CASE DEFAULT		      
	       WRITE(usl_static,21)j,jlim(j),r_zc(j)/rn,m_zc(j)/mstar
21	       FORMAT('Limite sans overshooting ni PC, limite ZR/ZC #',i4,
	1      ', couche #',i5,/,'rayon ZC/Rtot=',es10.3,
	2      ', m/Mstar au début de la ZC=',es10.3,/)
	      END SELECT		
	     ENDIF 	!ovshts
	    ENDIF  	!lconv
	   ENDDO    	!j limites ZR/ZC
	  ENDIF		!tot_conv
	 ENDIF		!ovsht

c	 CALL pause('après les overshoot')

c détermination de la fonction de poids fac 
c interpolation aux nouveaux q de P, T, R, L, M pour le calcul de la
c fonction de répartition aux new
c	 PRINT*,'ctel,ctep' ; WRITE(*,2000)ctel,ctep,ctem,cter,ctet
c	 PRINT* ; PRINT*,'lim/new/jlim',lim ; PRINT*,(new(jlim(i)),i=1,lim)
c	 PRINT*,(jlim(i),i=1,lim) ; WRITE(*,2000)new    

	 IF(.NOT.(tot_rad .OR. tot_conv))THEN
	  ALLOCATE(esp(n_qs))
	  DO i=1,n_qs
	   CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,new(i),lq,f,dfdq)
	   IF(no_croiss)PRINT*,'Problème localisé en 27 dans lim_zc'
	   esp(i)=ctep*f(1)+ctet*f(2)+cter*f(3)+ctel*f(4)+ctem*f(5)
c	   WRITE(*,2000)(f(j),j=1,5),esp(i),new(i)
	   IF(i > 1)THEN
	    DO WHILE(esp(i) <= esp(i-1))
	     esp(i)=esp(i)+0.01d0
	    ENDDO
	   ENDIF
	  ENDDO

c esp(i) est la valeur de la fonction de répartition aux nouveaux q
c déterminée avec les anciens facteurs fac(i)
c calcul des nouveaux fac(i) de facon à amener un noeud sur chaque nouveau q
	  DO i=1,n_qs-1
	   fac(i)=bp(6,1)/(esp(i+1)-esp(i))
c	   WRITE(*,2000)REAL(i),fac(i) 
	  ENDDO
c	  WRITE(*,2000)bp(6,1) ; PRINT*,n_qs
c	  WRITE(*,2000)(fac(i),new(i),i=1,n_qs-1)
c	  PRINT*,(lim,jlim(i),i=1,lim)
c	  PRINT*,bp(6,1),(fac(jlim(i)),new(jlim(i)),i=1,lim)
c	  PRINT*,(esp(jlim(i)-1),esp(jlim(i)),esp(jlim(i)+1),i=1,lim)
c	  CALL pause('esp')
	  DEALLOCATE(esp)
	  
c test fac=1.d0

c formation des limites ZR/ZC pour la convection
c avec un coeur convectif on prend l=Rcoeur-r, ce qui sera
c réalisé en prenant r_zc_conv(id_conv)=-r_zc(1)
	  IF(lconv(1))THEN
	   id_conv=1 ; r_zc_conv(id_conv)=r_zc(1)    !coeur pas convectif
	  ELSE
	   id_conv=0 ; r_zc_conv(id_conv)=-r_zc(1)   !coeur convectif
	  ENDIF
	  DO i=1,lim
	   r_zc_conv(i+id_conv)=r_zc(i)
	  ENDDO
	  IF(lconv(lim))THEN
	   if_conv=lim+1    !enveloppe totalement convective
	   r_zc_conv(if_conv)=MAX(rstar,r_zc(lim))
	  ELSE
	   if_conv=lim !enveloppe non totalement convective
	  ENDIF
	 ENDIF	 
	ENDIF  !lim=0
	
c new est toujours alloué	
	DEALLOCATE(new)

	IF(tot_rad .OR. tot_conv)THEN
	 fac=1.d0 ; id_conv=0 ; r_zc_conv(id_conv)=0.d0
	 if_conv=1 ; r_zc_conv(if_conv)=rstar
	ENDIF

	SELECT CASE(langue)
	CASE('english')
	 WRITE(usl_static,1022)
1022	 FORMAT('---------Search of the limits ZR/ZC(end)-----------')
	CASE DEFAULT
	 WRITE(usl_static,22)
22	 FORMAT('---------Recherche de limites ZR/ZC(fin)-----------')
	END SELECT
	 
	RETURN
 
	END SUBROUTINE lim_zc
