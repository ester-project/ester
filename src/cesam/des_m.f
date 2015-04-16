
c*************************************************************************

	SUBROUTINE des_m(fin,dt,teff)

c	routine private du module mod_cesam

c	dessin de variables en cours d'évolution
c	Auteur: P.Morel, Departement J.D. Cassini, O.C.A.
c	CESAM2k

c	dans le module mod_donnees ou dans le fichier device
c	adapter h, dh, ld, dl suivant l'écran dont on dispose

c	25 08 97 : mise en place des variables eulériennes

c entrées, significations évidentes
c	p,t,m,l,r,ro,grad_ad,grad_mj,alfa,delta,kap,cp,teff,age
c	n: nombre de couches
c	chim,mc,mct,nc,knotc: pour interpolation de la comp. chim.
c	lim: nombre de limites ZR/ZC
c	lconv: .TRUE. si début de ZC
c	m_zc: masse aux limites ZC/ZR
c	mstar: masse totale au temps du dessin

c NOTATIONS (hélas incohérentes) pour les développements sur B-splines
c	n_ch : nombre VARIABLE de points élément de mod_variables
c	nch : nombre FIXE de fonctions élément de mod_donnees
c	m_ch : ordre FIXE des splines élément de mod_donnees 
c	mch(n_ch) : abscisses VARIABLES élément de mod_variables

c------------------------------------------------------------
 
	USE mod_donnees, ONLY : device, dfesh_des, dl, dh, diffusion, dl_des,
	1 dteff_des, en_masse, fesh_des, h, ihe4, Krot, ld, ln10, logl_max,
	2 logl_min, logteff_max, logteff_min, l_des, mtot, m_ch,
	3 nb_max_modeles, nchim, ne, nom_elem, nom_fich2, nrot, ord_qs, ord_rot,
	4 precision, rsol, teff_des, xleft, ybot, zoom_l, zoom_t, zsx_sol
	USE mod_kind
	USE mod_numerique, ONLY : box, bsp1dn, no_croiss
	USE mod_variables, ONLY : age, bp, chim, chim_gram, knot, knotc, 
	1 knotr, lconv, lim, mc, mct, model_num, mstar, mrot, mrott, m_zc,
	2 n_ch, n_qs, n_rot, q, qt, rota
    
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in) :: dt, teff	
	
	LOGICAL, INTENT(in) :: fin

	REAL (kind=dp), SAVE, ALLOCATABLE, DIMENSION(:) :: dfdq, f, dxchim,
	1 xchim
	REAL (kind=dp), DIMENSION(nrot) :: dfrot, frot
	INTEGER, PARAMETER :: n_des=1000           
	REAL (kind=dp), SAVE, DIMENSION(n_des) :: qx
	REAL (kind=dp) :: dc, dq, nu, xdbl
      
	REAL (kind=sp), SAVE, ALLOCATABLE, DIMENSION(:,:) :: y
	REAL (kind=sp), SAVE, ALLOCATABLE, DIMENSION(:) :: yc, yloc, ymaxj, yw
	REAL (kind=sp), SAVE, DIMENSION(n_des) :: fesh, x, yd, yl, ym,
	1 yp, yr, yt
	REAL (kind=sp), SAVE, DIMENSION(4) :: feshf, xf
	REAL (kind=sp), SAVE, DIMENSION(5) :: ylocp	
	REAL (kind=sp), SAVE, DIMENSION(2) :: x4        
	REAL (kind=sp), PARAMETER :: inch=2.54
	REAL (kind=sp), SAVE ::  fesh0, llext, lmax,
	1 lteff, llextp, lteffp, pas=-1., pmax, rmax, teffs, tmax, vext, xloc,
	2 xlocc, xmax, xmax1, xmax3, xmax4, xmin4, xmin, xmin1, xmin3, x_age,
	3 ylocw, ymax, ymax1, ymax3, ymax4, ymin, ymin1, ymin3, ymin4,
	4 y_age, wmax
	
	REAL (kind=sp) :: coord, disp, fjust, xtick, ytick

	INTEGER :: i, ic, j, lq, nxsub, nysub
			
	LOGICAL, SAVE :: des_boite, des_fesh, init=.TRUE., lim_hr, zoom
   
	CHARACTER (len=10), SAVE :: side, xopt, yopt
	CHARACTER (len=15), SAVE, ALLOCATABLE, DIMENSION(:) :: htext
	CHARACTER (len=15), SAVE :: char      
	CHARACTER (len=23), SAVE :: ftext, ltex, ptex, rtex, tefftex, ttex, wtex
	CHARACTER (len=23), SAVE, ALLOCATABLE, DIMENSION(:) :: tex  
	CHARACTER (len=80), SAVE :: text, text_age
	CHARACTER (len=100), SAVE :: texptrl
	
c--------------------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE.
       
c	 lib=2.5e6        !température de disparition du Li

c les boites
	 des_boite=teff_des > 0. .AND. l_des > 0.
	 IF(des_boite)THEN
	  dteff_des=dteff_des/teff_des/ln10 ; teff_des=LOG10(teff_des)      
	  dl_des=dl_des/l_des/ln10 ; l_des=LOG10(l_des)	 
c	  WRITE(*,2000)teff_des,dteff_des,l_des,dl_des ; PAUSE'teff_des'
	 ENDIF
	 des_fesh=fesh_des < 100. .AND. des_boite
	 fesh0=0.5-fesh_des ; fesh0=0.1*ANINT(10.*fesh0)
	 fesh_des=fesh_des+fesh0
	
	 zoom=MAXVAL(ABS(zoom_t)) /= 0. .AND. MAXVAL(ABS(zoom_l)) /= 0.
	1 .AND. des_boite
	 lim_hr=logteff_max > 0.
	 
c	 WRITE(*,2000)logteff_max,logteff_min,logl_max,logl_min
c	 PRINT*,lim_hr ; PAUSE'des_m'

	 texptrl='P/P\dmax\u, T/T\dmax\u, L/L\dmax\u, R/R\dmax\u'	      

	 ALLOCATE(y(nchim,n_des),xchim(nchim),dxchim(nchim),yloc(nchim),
	1 tex(nchim),yc(nchim),ymaxj(nchim),f(ne),dfdq(ne),htext(nchim))
	
	 IF(Krot >= 3)THEN
	  ALLOCATE(yw(n_des)) ; ylocw=1.-0.08*6
	 ENDIF	

	 CALL pgbegin(0,device,1,1)
	 CALL pgscf(2)          !roman font 
c	 CALL pgscr(4,.3,.5,1.) !aménagement du bleu

c	 dans le module mod_donnees ou dans le fichier device
c	 adapter h, dh, ld, dl, y_age suivant l'écran dont on dispose

	 h=h/inch 		!hauteur du cadre      
	 dh=dh/inch		!saut entre cadres     
	 ld=ld/inch  		!largeur du cadre     
	 dl=dl/inch 		!espace en largeur     
	 xleft=xleft/inch ; ybot=ybot/inch
	 
c	 WRITE(*,2000)xleft,ld,ybot,h ; PAUSE'init'
	 
	 IF(nchim > 10)THEN
	  yloc=(/ (1.2-0.08*i, i=1,nchim) /)
	 ELSE
	  yloc=(/ (1.-0.08*i, i=1,nchim) /)
	 ENDIF
	 ylocp=(/ (1.-0.08*i, i=1,5) /)
   
c premier cadre: composition chimique
	 xmin=-0.01 ; xmax=mtot*1.01 ; xlocc=xmax*0.65 ; ymax=1.05 ; ymin=0.
	 CALL pgvsize(xleft,xleft+ld,ybot,ybot+h)
	 CALL pgwindow(xmin,xmax,ymin,ymax)
	 xopt='bcnst' ; xtick=0. ; nxsub=0 ; yopt='bcnst' ; ytick=0.
	 nysub=0
	 CALL pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
	 side='b'   !ABScisses
	 disp=2 ; coord=0.5 ; fjust=0.5
	 text='M en M\d\(2281)\u'
	 CALL pgmtext(side,disp,coord,fjust,text)
	 side='l'   !ordonnees
	 disp=2 ; coord=0.5 ; fjust=0.5 ; text='X\di\u(t)'
	 CALL pgmtext(side,disp,coord,fjust,text)
	 
c tracé des limites de [Z/H] 
	 IF(des_fesh)THEN
	  feshf(1)=fesh_des+dfesh_des(1) ; feshf(2)=feshf(1)
	  feshf(3)=fesh_des+dfesh_des(2) ; feshf(4)=feshf(3)
	  xf(1)=mtot*0.95 ; xf(2)=mtot ; xf(3)=xf(2) ; xf(4)=xf(1)
c	  CALL pgsfs(4) ; CALL pgpoly(4,xf,feshf) !; PAUSE'poly'	 
c	  WRITE(*,2000)xf,feshf ; PAUSE'poly1'
	 ENDIF		

c second cadre: P, T, R, L, M, (Omega si Krot=3) en f(m)
	 xmin1=-0.01 ; xmax1=xmax ; xloc=xmax1*0.5 ; ymin1=-0.01 ; ymax1=1.1
	 IF(Krot >= 3)THEN
	  x_age=(4.*xmin1+xmax1)/5.
	 ELSE	 
	  x_age=(2.*xmin1+xmax1)/3.
	 ENDIF
	 y_age=ymax1+0.2
	 
	 CALL pgvsize(xleft,xleft+ld,ybot+h+dh,ybot+h+dh+h)
	 CALL pgwindow(xmin1,xmax1,ymin1,ymax1)
	 xopt='bcnst' ; xtick=0. ; nxsub=0 ; yopt='bcnst' ; ytick=0.
	 nysub=0
	 CALL pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
	 side='b'   !abscisses
	 disp=2 ; coord=0.5 ; fjust=0.5 ; text='M en M\d\(2281)\u'
	 CALL pgmtext(side,disp,coord,fjust,text)
	 side='l'   !ordonnees
	 disp=2 ; coord=0.5 ; fjust=0.5
	 CALL pgsch(0.9)	 
	 CALL pgmtext(side,disp,coord,fjust,texptrl)
	 side='t'   !sommet
	 disp=1 ; coord=0.5 ; fjust=0.5
	 CALL pgmtext(side,disp,coord,fjust,nom_fich2)
      
c troisième cadre: diagramme HR
	 lteffp=LOG10(teff) ; lq=n_qs
	 CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(n_qs),lq,f,dfdq)
	 IF(no_croiss)PRINT*,'Pb. en 1 dans des_m'
	 IF(en_masse)THEN 
	  llextp=LOG10(f(4))*1.5
	 ELSE
	  llextp=LOG10(f(4))     
	 ENDIF
	 
	 IF(lim_hr)THEN
	  xmin3=logteff_max ; xmax3=logteff_min
	  ymax3=logl_max ; ymin3=logl_min
	 ELSEIF(zoom)THEN	
	  xmin3=teff_des+zoom_t(1) ; xmax3=teff_des+zoom_t(2)
	  ymax3=l_des+zoom_l(1) ; ymin3=l_des+zoom_l(2)
	 ELSEIF(mtot <= 0.80d0)THEN
	  xmin3=3.85 ; xmax3=3.5 ; ymax3=2. ; ymin3=-1.5	  	  
	 ELSEIF(mtot <= 0.95d0)THEN
	  xmin3=3.75 ; xmax3=3.5 ; ymax3=3. ; ymin3=-1.
	 ELSEIF(mtot <= 1.5d0)THEN
	  xmin3=3.90 ; xmax3=3.5 ; ymax3=3. ; ymin3=-0.5             
	 ELSEIF(mtot <= 2.5d0)THEN
	  xmin3=4.05 ; xmax3=3.5 ; ymax3=3.5 ; ymin3=0.                    
	 ELSEIF(mtot <= 6.0d0)THEN
	  xmin3=4.3 ; xmax3=3.5 ; ymax3=5. ; ymin3=1.5
	 ELSE
	  xmin3=4.6 ; xmax3=3.5 ; ymax3=6. ; ymin3=3.
	 ENDIF        
	 CALL pgvsize(xleft+ld+dl,xleft+2*ld+dl,ybot,ybot+h)
	 CALL pgwindow(xmin3,xmax3,ymin3,ymax3)
	 xopt='bcnst' ; xtick=0. ; nxsub=0 ; yopt='bcnst' ; ytick=0.
	 nysub=0
	 CALL pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
	 side='b'   !abscisses
	 disp=2 ; coord=0.5 ; fjust=0.5
	 text='HR --- Log\d10\u T\deff\u'
	 CALL pgmtext(side,disp,coord,fjust,text)
	 side='l'   !ordonnees
	 disp=2 ; coord=0.5 ; fjust=0.5
	 text='HR --- Log\d10\u L/L\d\(2281)\u'
	 CALL pgmtext(side,disp,coord,fjust,text)
	
	 ! IF(des_boite)CALL box(teff_des,dteff_des,l_des,dl_des)
	 
c quatrième cadre: ZC
	 xmin4=-3. ; xmax4=MAX(nb_max_modeles+10,50)
	 ymin4=-.05 ; ymax4=mtot*1.05  
	 CALL pgvsize(xleft+ld+dl,xleft+2*ld+dl,ybot+h+dh,ybot+h+dh+h)
	 CALL pgwindow(xmin4,xmax4,ymin4,ymax4)
	 xopt='bcnst'
	 xtick=0. ; nxsub=0 ; yopt='bcnst' ; ytick=0. ; nysub=0
	 CALL pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
	 side='b'   !ABScisses
	 disp=2 ; coord=0.5 ; fjust=0.5
	 text='pas temporel'
	 CALL pgmtext(side,disp,coord,fjust,text)
	 side='l'   !ordonnees
	 disp=2 ; coord=0.5 ; fjust=0.5
	 text='ZC ---- M/M\d\(2281)\u'
	 CALL pgmtext(side,disp,coord,fjust,text) 
c	 PAUSE
	 side='t'   !sommet
	 disp=1 ; coord=0.5 ; fjust=0.5
	 CALL pgmtext(side,disp,coord,fjust,nom_fich2)
c	 CALL pgslw(1)          !epaisseur du trait

	ELSE

c à partir du second passage on efface
c	 WRITE(*,2000)xleft,ld,ybot,h ; PAUSE'efface'
	 CALL pgsci(0) ; CALL pgsls(1)
	 CALL pgsch(0.8)  !hauteur des caractères  
	 CALL pgvsize(xleft,xleft+ld,ybot,ybot+h) 
	 CALL pgwindow(xmin,xmax,ymin,ymax) 
	 DO j=1,nchim
	  yd=y(j,:)
	  CALL pgline(n_des,x,yd) ; CALL pgtext(xlocc,yloc(j),tex(j))
	 ENDDO
	 IF(des_fesh)THEN
	  CALL pgsls(4) ; CALL pgline(n_des,x,fesh)
	 ENDIF 
	 CALL pgsls(1)       
	 CALL pgvsize(xleft,xleft+ld,ybot+h+dh,ybot+h+dh+h)
	 CALL pgwindow(xmin1,xmax1,ymin1,ymax1)
	 CALL pgline(n_des,ym,yp) ; CALL pgline(n_des,ym,yt)
	 CALL pgline(n_des,ym,yl) ; CALL pgline(n_des,ym,yr)
	 CALL pgtext(xloc,ylocp(1),ptex) ; CALL pgtext(xloc,ylocp(2),ttex)
	 CALL pgtext(xloc,ylocp(3),ltex) ; CALL pgtext(xloc,ylocp(4),rtex)
	 CALL pgtext(xloc,ylocp(5),tefftex)
	 IF(Krot >= 3)THEN
	  CALL pgline(n_des,ym,yw) ;  CALL pgtext(xloc,ylocw,wtex)
	 ENDIF

	 CALL pgsch(1.) ; CALL pgtext(x_age,y_age,text_age)
c	 PAUSE'efface '//text_age

	 CALL pgsci(1)  
	ENDIF       !fin d'initialisation et d'effacage
      
	dc=mc(n_ch)/(n_des-1)
	DO i=1,n_des
	 xdbl=dc*(i-1) ; xdbl=MIN(MAX(xdbl,mc(1)),mc(n_ch))
	 x(i)=SQRT(xdbl)**3
	 CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,knotc,.TRUE.,xdbl,lq,
	1 xchim,dxchim)
	 IF(no_croiss)PRINT*,'Pb. en 2 dans des_m'	
	 CALL chim_gram(xchim,dxchim)
	 y(:,i)=MAX(xchim(:),1.d-30)
c	 write(*,2000)y(:,i)
	 IF(y(1,i) > 1.d-4 .AND. ihe4 > 1)THEN 		!calcul de [Z/H]
	  fesh(i)=LOG10(SUM(y(ihe4+1:nchim,i))/y(1,i)/zsx_sol)
	 ELSE
	  fesh(i)=10.
	 ENDIF 
	ENDDO
	fesh=fesh+fesh0
c	PAUSE'y(:,i)'
   
c normalisation des abondances, sauf H et He4
	DO j=1,nchim
	 ymaxj(j)=MAXVAL(y(j,:)) ; yc(j)=y(j,1)
	 IF(j /= 1 .AND. j /= ihe4)y(j,:)=y(j,:)/ymaxj(j)	            
	ENDDO
	CALL pgsci(1) ; CALL pgsls(1)
	CALL pgvsize(xleft,xleft+ld,ybot,ybot+h)    
	CALL pgwindow(xmin,xmax,ymin,ymax)      
	CALL pgsch(0.8)   !hauteur des caractères
	ic=1
	DO j=1,nchim	
	 IF(precision == 'av' .OR. precision == 'ar')THEN
	  WRITE(htext(j),10)yc(j)
10	  FORMAT(es9.2)	  
	 ELSE
	  WRITE(htext(j),10)ymaxj(j)
	 ENDIF	 
	 tex(j)=nom_elem(j)//'='//htext(j)
	 ic=ic+1 ; IF(ic > 13)ic=2 ; IF(ic == 4)ic=5 ; CALL pgsci(ic)
	 yd=y(j,:) ; CALL pgline(n_des,x,yd)    
	 CALL pgtext(xlocc,yloc(j),tex(j))
	ENDDO 
	IF(des_fesh)THEN		!trace [Z/H]
	 CALL pgsfs(4) ; CALL pgpoly(4,xf,feshf) ;
	 CALL pgsls(4) ; CALL pgline(n_des,x,fesh) ; CALL pgsls(1)
	 CALL pgsch(0.6)
	 IF(fesh0 > 0.)THEN
	  WRITE(ftext,1)fesh0 
1	  FORMAT('[Z/X]+',f4.1)
	 ELSEIF(fesh0 < 0.)THEN
	  WRITE(ftext,2)fesh0 
2	  FORMAT('[Z/X]',f4.1)
	 ELSE
	  ftext='[Z/X]'
	 ENDIF	 	 
	 CALL pgtext(xf(1)-0.08,feshf(1)+0.02,TRIM(ftext))	 
c	 WRITE(*,2000)xf,feshf,fesh0 ; PRINT*,ftext ; PAUSE'poly2'
	ENDIF			
	CALL pgsci(1) ; CALL pgsch(1.)    !hauteur des caractères

c P, T, R, L, R, (Omega si Krot >= 3) en fonction de M
	CALL pgvsize(xleft,xleft+ld,ybot+h+dh,ybot+h+dh+h)
	CALL pgwindow(xmin1,xmax1,ymin1,ymax1)
c	CALL pgsch(.8)                !hauteur des caractères
c	CALL pgtext(xloc,ylocp(1),ptex) ; CALL pgtext(xloc,ylocp(2),ttex)
c	CALL pgtext(xloc,ylocp(3),ltex) ; CALL pgtext(xloc,ylocp(4),rtex)
c	CALL pgsch(1.)    !hauteur des caractères  
c	CALL pgsls(1)

	dq=(n_qs-1.d0)/dble(n_des-1)  !pour les variables      
	qx=(/ (dq*(i-1)+1.d0, i=1,n_des) /)    
	DO i=1,n_des
	 CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,min(qx(i),q(n_qs)),
	1 lq,f,dfdq)
	 IF(no_croiss)PRINT*,'Pb. en 3 dans des_m'
	 yp(i)=f(1) ; yt(i)=f(2)
	 IF(en_masse)THEN 
	  yr(i)=SQRT(ABS(f(3))) ; yl(i)=SQRT(ABS(f(4)))**3
	  ym(i)=SQRT(ABS(f(5)))**3
	 ELSE
	  yr(i)=ABS(f(3)) ; yl(i)=f(4) ; ym(i)=ABS(f(5))
	 ENDIF                    
	ENDDO
	
	IF(Krot >= 3)THEN
	 DO i=1,n_des
	  nu=ym(i)**(2./3.)
	  CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,
	1 MAX(mrot(1),MIN(nu,mrot(n_rot))),lq,frot,dfrot)
	  IF(no_croiss)PRINT*,'Pb. en 4 dans des_m'
	  yw(i)=frot(1)
	 ENDDO
	 vext=yw(n_des)*yr(n_des)*rsol*1.d-5
	 wmax=MAXVAL(yw) ; yw=yw/wmax
	ENDIF

c écriture de l'entête	
	IF(Krot >= 3)THEN
	 WRITE(text_age,3)model_num,age,dt,vext
3	 FORMAT('modèle ',i4.4,', âge en My: ',es10.3,', dt=',es10.3,
	1 ', V\dext\u=',es10.3,'Km/s')
	 CALL pgtext(x_age,y_age,text_age)
	ELSE
	 WRITE(text_age,9)model_num,age,dt
9	 FORMAT('modèle ',i4.4,', âge en My: ',es10.3,', dt=',es10.3)
	 CALL pgtext(x_age,y_age,text_age)
	ENDIF
	    
c normalisations
	yp=EXP(yp) ; yt=EXP(yt)	
	pmax=MAXVAL(yp) ; tmax=MAXVAL(yt)   
	rmax=MAXVAL(yr); lmax=MAXVAL(yl)
	yp=yp/pmax ; yt=yt/tmax ; yr=yr/rmax ; yl=yl/lmax

c	CALL pgvsize(xleft,xleft+ld,ybot+h+dh,ybot+h+dh+h)
c	CALL pgwindow(xmin,xmax,ymin,ymax)
      
	CALL pgsch(0.8)   !hauteur des caractères 
11	FORMAT(es10.3)      
	WRITE(char,11)pmax ; ptex='P\dmax\u='//char ; CALL pgsci(2)
	CALL pgline(n_des,ym,yp) ; CALL pgtext(xloc,ylocp(1),ptex)
	WRITE(char,11)tmax ; ttex='T\dmax\u='//char ; CALL pgsci(3)
	CALL pgline(n_des,ym,yt) ; CALL pgtext(xloc,ylocp(2),ttex)   
	WRITE(char,11)lmax ; ltex='L\dmax\u='//char ; CALL pgsci(6)
	CALL pgline(n_des,ym,yl) ; CALL pgtext(xloc,ylocp(3),ltex)
	WRITE(char,11)rmax ; rtex='R\dmax\u='//char ; CALL pgsci(7)
	CALL pgline(n_des,ym,yr) ; CALL pgtext(xloc,ylocp(4),rtex)
	teffs=teff ; WRITE(char,11)teffs ; tefftex='T\deff\u='//char
	CALL pgsci(9) ;	CALL pgtext(xloc,ylocp(5),tefftex)
	IF(Krot >= 3)THEN
	 WRITE(char,10)wmax ; wtex='\gW\dmax\u='//char ; CALL pgsci(8)
	 CALL pgline(n_des,ym,yw) ; CALL pgtext(xloc,ylocw,wtex)
	ENDIF	   
	CALL pgsch(1.)    !hauteur des caractères
c	PAUSE 
	CALL pgsci(1)

c diagramme HR
	lteff=LOG10(teff) ; lq=n_qs
	CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(n_qs),lq,f,dfdq)
	IF(no_croiss)PRINT*,'Pb. en 5 dans des_m'
	IF(en_masse)THEN  
	 llext=LOG10(f(4))*1.5
	ELSE
	 llext=LOG10(f(4))
	ENDIF
	CALL pgvsize(xleft+ld+dl,xleft+2*ld+dl,ybot,ybot+h)
	CALL pgwindow(xmin3,xmax3,ymin3,ymax3)
	CALL pgmove(lteffp,llextp) ; CALL pgdraw(lteff,llext)
	lteffp=lteff ; llextp=llext
      
c zones radiatives, convectives
c	PRINT*,lim,(lconv(i),i=1,lim) ; PRINT*,(m_zc(i),i=1,lim),age
	pas=pas+1. ; x4(1)=pas ; x4(2)=pas

	IF(lim > 0)THEN   !lim=0: modele totalement radiatif pas de dessin      
	 CALL pgvsize(xleft+ld+dl,xleft+2.*ld+dl,ybot+h+dh,ybot+2.*h+dh)

c	 PRINT*,lim,(lconv(i),i=1,lim)
c	 WRITE(*,2000)(m_zc(i),i=1,lim),mstar,mtot

	 CALL pgwindow(xmin4,xmax4,ymin4,ymax4)
	 IF(lim == 1 .and. .not.lconv(1) .and.
	1 m_zc(1)/mstar .ge. .99)THEN!modele totalement convectIF
	  yd(1)=0. ; yd(2)=mstar ; CALL pgline(2,x4,yd)
	 ELSE
	  i=1
	  DO WHILE(i <= lim)
	   IF(lconv(i))THEN     !debut de ZC
  	    yd(1)=m_zc(i)
	    IF(i == lim)THEN          !ZC externe
	     yd(2)=mstar
	    ELSE
	     yd(2)=m_zc(i+1) ; i=i+1
	    ENDIF      
	    CALL pgline(2,x4,yd)
	   ELSE
	    IF(i == 1)THEN                  !fin de ZC centrale
	     yd(1)=0. ; yd(2)=m_zc(1) ; CALL pgline(2,x4,yd)
	    ELSE
	     PRINT*,'problème' ; CALL pgend ; STOP
	    ENDIF
	   ENDIF
	   i=i+1
	  ENDDO     !WHILE
	 ENDIF
	ENDIF

	IF(fin)CALL pgend
  
	RETURN
      
	END SUBROUTINE des_m
