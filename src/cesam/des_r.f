
c*************************************************************************

      SUBROUTINE des_r(fin,dt,teff)

c routine du module mod_cesam      
c     dessin de variables en cours d'evolution
c     Auteur: P. Morel, Departement J.D. Cassini, O.C.A.
c     CESAM2k

c     25 08 97 : mise en place des variables eulériennes

c entrées, significations évidentes
c     p,t,m,l,r,ro,grad_ad,grad_mj,alfa,delta,kap,cp,teff,age
c     n: nombre de couches
c     chim,mc,mct,n_ch,knotc: pour interpolation de la comp. chim.
c     lim: nombre de limites ZR/ZC
c     lconv: .TRUE. si debut de ZC
c     m_zc: masse aux limites ZC/ZR
c     mstar: masse totale au temps du dessin

c NOTATIONS (hélas incohérentes) pour les développements sur B-splines
c	n_ch : nombre VARIABLE de points, élément de mod_variables
c	nch : nombre FIXE de fonctions, élément de mod_donnees
c	m_ch : ordre FIXE des splines, élément de mod_donnees 
c	mch(n_ch) : abscisses VARIABLES, élément de mod_variables
      
c------------------------------------------------------------

	USE mod_donnees, ONLY : device, dl, dh, diffusion, en_masse, h, ld,
	1 ln10, ihe4, Krot, mtot, m_ch, nchim, nb_max_modeles, ne,
	2 nom_elem, nom_fich2, nrot, ord_qs, ord_rot, rsol, xleft, ybot
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn, no_croiss
	USE mod_variables, ONLY : age, bp, chim, chim_gram, knot, knotc, 
	1 knotr, lconv, lim, mc, mct, model_num, mstar, mrot, mrott, m_zc, m23,
	2 n_ch, n_qs, n_rot, q, qt, r2, inter, rota
 
	IMPLICIT NONE
      
	REAL (kind=dp), INTENT(in) :: dt, teff	
	LOGICAL, INTENT(in) :: fin

	REAL (kind=dp), SAVE, ALLOCATABLE, DIMENSION(:) :: dfdq, dxchim, f,
	1 xchim
	REAL (kind=dp), DIMENSION(nrot) :: dfrot, frot
	INTEGER, PARAMETER :: n_des=1000           
	REAL (kind=dp), SAVE, DIMENSION(n_des) :: qx
	REAL (kind=dp) :: dc, dq, nu, xdbl

	REAL (kind=sp), SAVE, ALLOCATABLE, DIMENSION(:,:) :: y
	REAL (kind=sp), SAVE, ALLOCATABLE, DIMENSION(:) :: yloc, ymaxj, yw    
	REAL (kind=sp), SAVE, DIMENSION(n_des) :: x, yd, yl, ym, yp, yr, yt
	REAL (kind=sp), SAVE, DIMENSION(5) :: ylocp
	REAL (kind=sp), SAVE, DIMENSION(2) :: rt, tli, xt, x4
	REAL (kind=sp), PARAMETER :: fac=4./3., fac1=0.9*fac,
	1 fac2=0.75, inch=2.54, lib=2.5e6
	REAL (kind=sp), SAVE :: llext, llextp, lteff, lteffp, pas=-1.,
	1 lmax, pmax, rtot, rtotp,teffs, tmax,  vext, xloc, xmax, xmax3, xmax4,
	2 xmin, xmin3, xmin4, x_age, ylocw, ymax, ymax3, ymax4, ymin,
	3 ymin3, ymin4, y_age
	REAL (kind=sp):: coord, disp, fjust, mmax, wmax, xtick, ytick

	INTEGER :: i, ic, j, lq, nxsub, nysub
                
	LOGICAL, SAVE :: init=.TRUE.
	LOGICAL :: change

	CHARACTER (len=10) :: side, xopt, yopt
	CHARACTER (len=15), SAVE, ALLOCATABLE, DIMENSION(:) :: htext
	CHARACTER (len=15) :: char          
	CHARACTER (len=80), SAVE :: ltex, mtex, ptex, tefftex, ttex, wtex
	CHARACTER (len=80), SAVE, ALLOCATABLE, DIMENSION(:) :: tex
	CHARACTER (len=80) :: text, text_age
	CHARACTER (len=100), SAVE :: texptrl
	
c----------------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE.

c      fac facteur en r, si rtot > fac rtotp, on change d'échelle en R
c      lib=2.5e6        !température de disparition du Li

	 ALLOCATE(y(nchim,n_des),xchim(nchim),dxchim(nchim),yloc(nchim),
	1 tex(nchim),ymaxj(nchim),f(ne),dfdq(ne),htext(nchim))

	 texptrl='P/P\dmax\u, T/T\dmax\u, L/L\dmax\u, M/M\dmax\u'	      

	 IF(Krot >= 3)THEN
	  ALLOCATE(yw(n_des)) ; ylocw=1.-0.08*6
	 ENDIF	

c	sous window : device=/ws, /NULL ou ?
c	sous linux/unix : device=/xw, /NULL ou ?

	 CALL pgbegin(0,device,1,1)
	 CALL pgscf(2)          	!roman font
	 CALL pgsch(1.)  	!hauteur des caractères
	 CALL pgslw(1)            !épaisseur du trait     
c	 CALL pgscr(4,.3,.5,1.) 	!aménagement du bleu     

	 dc=mtot/(n_des-1)

c	 dans le module mod_donnees ou dans le fichier device
c	 adapter h, dh, ld, dl, y_age suivant l'écran dont on dispose

	 h=h/inch 		!hauteur du cadre      
	 dh=dh/inch		!saut entre cadres     
	 ld=ld/inch  	!largeur du cadre     
	 dl=dl/inch 		!espace en largeur     
	 xleft=xleft/inch ; ybot=ybot/inch
	 
c	 WRITE(*,2000)xleft,ld,ybot,h ; PAUSE'init'	 
   
	 IF(nchim > 10)THEN
	  yloc=(/ (1.1-0.08*i, i=1,nchim) /)
	 ELSE
	  yloc=(/ (1.-0.08*i, i=1,nchim) /)
	 ENDIF
	 ylocp=(/ (1.-0.08*i, i=1,5) /)

c Rtot
	 CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(n_qs),lq,f,dfdq)
	 IF(no_croiss)PRINT*,'Pb. en 1 dans des_r'
	 IF(en_masse)THEN
	  rtot=SQRT(f(3))
	 ELSE
	  rtot=f(3)
	 ENDIF
	 rtotp=rtot ; xt(1)=rtot ; xt(2)=rtot ; rt(1)=-1. ; rt(2)=1.4

c premier cadre: abondances
	xmin=-0.01 ; xmax=rtotp*fac ; xloc=xmax*0.6
	ymax=1.05 ; ymin=-0.05       
	CALL pgvsize(xleft,xleft+ld,ybot,ybot+h)
	CALL pgwindow(xmin,xmax,ymin,ymax)
	xopt='bcnst' ; xtick=0. ; nxsub=0 ; yopt='bcnst' ; ytick=0.
	nysub=0
	CALL pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
	side='b'   !abscisses
	disp=2 ; coord=0.5 ; fjust=0.5
	text='R/R\d\(2281)\u'
	CALL pgmtext(side,disp,coord,fjust,text)
	side='l'   !ordonnees
	disp=2 ; coord=0.5 ; fjust=0.5
	text='X\di\u(t)'
	CALL pgmtext(side,disp,coord,fjust,text) 

c second cadre: P, T, R, L, M, (Omega si Krot>=3) en f(r)
       CALL pgvsize(xleft,xleft+ld,ybot+h+dh,ybot+h+dh+h)
       CALL pgwindow(xmin,xmax,ymin,ymax)
       xopt='bcnst' ; xtick=0. ; nxsub=0 ; yopt='bcnst' ; ytick=0. ; nysub=0
       CALL pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
       side='b'   !abscisses
       disp=2 ; coord=0.5 ; fjust=0.5 ; text='R/R\d\(2281)\u'
       CALL pgmtext(side,disp,coord,fjust,text)
       side='l'   !ordonnees
       disp=2 ; coord=0.5 ; fjust=0.5
       CALL pgsch(0.9)             
       CALL pgmtext(side,disp,coord,fjust,texptrl)
       CALL pgsch(1.)       
       side='t'   !sommet
       disp=1 ; coord=0.5 ; fjust=0.5
       CALL pgmtext(side,disp,coord,fjust,nom_fich2)

c  troisième cadre: diagramme HR
       lteffp=LOG10(teff) ; lq=n_qs
       CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(n_qs),lq,f,dfdq)
       IF(no_croiss)PRINT*,'Pb. en 2 dans des_r'
       IF(en_masse)THEN 
        llextp=LOG10(f(4))*1.5
       ELSE
        llextp=LOG10(f(4))
       ENDIF
       IF(mtot <= 0.80d0)THEN
	xmin3=3.85 ; xmax3=3.68 ; ymax3=2. ; ymin3=-1.
       ELSEIF(mtot <= 0.95d0)THEN
        xmin3=3.75 ; xmax3=3.5 ; ymax3=3. ; ymin3=-1.
       ELSEIF(mtot <= 1.5d0)THEN
        xmin3=3.90 ; xmax3=3.5 ; ymax3=3. ; ymin3=-0.5
c	xmin3=3.855 ; xmax3=3.80 ; ymax3=0.82 ; ymin3=0.54
       ELSEIF(mtot <= 2.5d0)THEN
        xmin3=4.05 ; xmax3=3.5 ; ymax3=3.5 ; ymin3=0.
       ELSEIF(mtot <= 6.0d0)THEN
        xmin3=4.3 ; xmax3=3.5 ; ymax3=5. ; ymin3=2.
       ELSE
        xmin3=4.6 ; xmax3=3.5 ; ymax3=6. ; ymin3=3.
       ENDIF
       CALL pgvsize(xleft+ld+dl,xleft+2*ld+dl,ybot,ybot+h)
       CALL pgwindow(xmin3,xmax3,ymin3,ymax3)
       xopt='bcnst' ; xtick=0. ; nxsub=0 ; yopt='bcnst'
       ytick=0. ; nysub=0
       CALL pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
       side='b'   !abscisses
       disp=2 ; coord=0.5 ; fjust=0.5
       text='HR --- Log\d10\u T\deff\u'
       CALL pgmtext(side,disp,coord,fjust,text)
       side='l'   !ordonnees
       disp=2 ; coord=0.5 ; fjust=0.5
       text='HR --- Log\d10\u L/L\d\(2281)\u'
       CALL pgmtext(side,disp,coord,fjust,text)

c quatrième cadre: ZC
       xmin4=-3. ; xmax4=MAX(nb_max_modeles+10,50)
       ymin4=-0.05 ; ymax4=mtot*1.05
       CALL pgvsize(xleft+ld+dl,xleft+2*ld+dl,ybot+h+dh,ybot+h+dh+h)
       CALL pgwindow(xmin4,xmax4,ymin4,ymax4)
       xopt='bcnst' ; xtick=0. ; nxsub=0 ; yopt='bcnst' ; ytick=0. ; nysub=0
       CALL pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
       side='b'   !abscisses
       disp=2 ; coord=0.5 ; fjust=0.5
       text='pas temporel'
       CALL pgmtext(side,disp,coord,fjust,text)
       side='l'   !ordonnees
       disp=2 ; coord=0.5 ; fjust=0.5
       text='ZC ---- M/M\d\(2281)\u'
       CALL pgmtext(side,disp,coord,fjust,text) 
       side='t'   !sommet
       disp=1 ; coord=0.5 ; fjust=0.5
       CALL pgmtext(side,disp,coord,fjust,nom_fich2)
c       CALL pgslw(1)          !epaisseur du trait
       init=.FALSE.

	ELSE

c à partir du second passage on efface les courbes
c	 WRITE(*,2000)xleft,ld,ybot,h ; PAUSE'efface'
	 
	 CALL pgsci(0) ; CALL pgsls(1)
	 CALL pgsch(0.8)  !hauteur des caractères
	 CALL pgvsize(xleft,xleft+ld,ybot,ybot+h) !les abondances
	 CALL pgwindow(xmin,xmax,ymin,ymax)
	 DO j=1,nchim
	  yd=y(j,:)	!WRITE(*,2000)yd(1:n_des) ; PAUSE'4,0'
	  CALL pgline(n_des,x,yd) ; CALL pgtext(xloc,yloc(j),tex(j))
	 ENDDO
	 CALL pgvsize(xleft,xleft+ld,ybot+h+dh,ybot+h+dh+h)   !P, T, L, M
c	  PRINT*,n_des ; WRITE(*,2000)xmin,xmax,ymin,ymax
c	  WRITE(*,2000)yr(1:n_des) ; WRITE(*,2000)yp(1:n_des) ; PAUSE'4,11'
	 CALL pgwindow(xmin,xmax,ymin,ymax) ; CALL pgline(n_des,yr,yp)
	 CALL pgline(n_des,yr,yt) ; CALL pgline(n_des,yr,yl)
	 CALL pgline(n_des,yr,ym)
	 CALL pgtext(xloc,ylocp(1),ptex) ; CALL pgtext(xloc,ylocp(2),ttex)
	 CALL pgtext(xloc,ylocp(3),ltex) ; CALL pgtext(xloc,ylocp(4),mtex)
	 CALL pgtext(xloc,ylocp(5),tefftex)	
	 IF(Krot >= 3)THEN
	  CALL pgline(n_des,yr,yw) ;  CALL pgtext(xloc,ylocw,wtex)
	 ENDIF          
	 CALL pgsch(1.)
	 CALL pgtext(x_age,y_age,text_age)
	 CALL pgsci(1) 
	ENDIF       !fin d'initialisation et d'effacage 
      
c changement d'échelle  
	CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(n_qs),lq,f,dfdq)
	IF(no_croiss)PRINT*,'Pb. en 3 dans des_r'
	IF(en_masse)THEN
	rtot=SQRT(f(3))
	ELSE
	rtot=f(3)
	ENDIF
	change=rtot < rtotp*fac2 .or. rtot > rtotp*fac1
	IF(change)THEN
	CALL pgsci(0)    !effacement des cadres et Rtot
	CALL pgsch(1.)  	!hauteur des caractères

c effacement du cadre des abondances               
       CALL pgvsize(xleft,xleft+ld,ybot,ybot+h)
       CALL pgwindow(xmin,xmax,ymin,ymax)
       xopt='bcnst' ; xtick=0. ; nxsub=0 ; yopt='bcnst'
       ytick=0. ; nysub=0
       CALL pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
       side='b'   !abscisses
       disp=2 ; coord=0.5 ; fjust=0.5
       text='R/R\d\(2281)\u'
       CALL pgmtext(side,disp,coord,fjust,text)
       side='l'   !ordonnees
       disp=2 ; coord=0.5 ; fjust=0.5 ; text='X\di\u(t)'
       CALL pgmtext(side,disp,coord,fjust,text)

c effacement de l'ancien rayon 
       CALL pgsls(2) ; CALL pgline(2,xt,rt)
       
c effacement du cadre  P, T, L, M
       CALL pgsls(1)     
       CALL pgvsize(xleft,xleft+ld,ybot+h+dh,ybot+h+dh+h)
       CALL pgwindow(xmin,xmax,ymin,ymax)
       xopt='bcnst' ; xtick=0. ; nxsub=0 ; yopt='bcnst' ; ytick=0.
       nysub=0
       CALL pgsci(0)     
       CALL pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
       side='b'   !abscisses
       disp=2 ; coord=0.5 ; fjust=0.5
       text='R/R\d\(2281)\u'
       CALL pgmtext(side,disp,coord,fjust,text)
       side='l'   !ordonnees
       disp=2 ; coord=0.5 ; fjust=0.5
       CALL pgsch(0.9)       
       CALL pgmtext(side,disp,coord,fjust,texptrl)
       CALL pgsch(1.)
       side='t'   !sommet
       disp=1 ; coord=0.5 ; fjust=0.5
       CALL pgmtext(side,disp,coord,fjust,nom_fich2)
       CALL pgsls(2)          !l'ancien rayon
       CALL pgline(2,xt,rt) ; CALL pgsls(1)
       
c nouvelle échelle
       xmax=rtot*fac ; xloc=xmax*.65 ; xt(1)=rtot ; xt(2)=rtot
       rt(1)=-1. ; rt(2)=1.4 ; rtotp=rtot
c      WRITE(*,2000)xt(1),xt(2),rt(1),rt(2),rtot,rtotp ; PAUSE'change rtot'
       CALL pgsci(1)

       CALL pgvsize(xleft,xleft+ld,ybot,ybot+h) !abondances
       CALL pgwindow(xmin,xmax,ymin,ymax)
       xopt='bcnst' ; xtick=0. ; nxsub=0 ; yopt='bcnst' ; ytick=0. ; nysub=0
       CALL pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
       side='b'   !abscisses
       disp=2 ; coord=0.5 ; fjust=0.5
       text='R/R\d\(2281)\u'
       CALL pgmtext(side,disp,coord,fjust,text)
       side='l'   !ordonnees
       disp=2 ; coord=0.5 ; fjust=0.5 ; text='X\di\u(t)'
       CALL pgmtext(side,disp,coord,fjust,text)
      ENDIF

      IF(en_masse)THEN  
       CALL inter('m23',bp,q,qt,n_qs,knot,mc(n_ch),f,dfdq,r2,m23)
      ELSE
       CALL inter('m23',bp,q,qt,n_qs,knot,SQRT(mc(n_ch))**3,f,dfdq,r2,m23)
      ENDIF 
      dc=f(3)/(n_des-1) ; ymaxj=0.
      DO i=1,n_des
       xdbl=dc*(i-1)
       IF(en_masse)THEN
        x(i)=SQRT(xdbl)
       ELSE
        x(i)=xdbl
       ENDIF
       CALL inter('r2 ',bp,q,qt,n_qs,knot,xdbl,f,dfdq,r2,m23)
       IF(en_masse)THEN
        xdbl=f(5)
       ELSE
        xdbl=f(5)**(2./3.)
       ENDIF
       xdbl=MIN(MAX(xdbl,mc(1)),mc(n_ch))
       CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,knotc,.TRUE.,xdbl,lq,
     1  xchim,dxchim)
       IF(no_croiss)PRINT*,'Pb. en 4 dans des_r'    
       CALL chim_gram(xchim,dxchim)
       y(:,i)=MAX(xchim(:),1.d-30)
      ENDDO

c normalisation des abondances, sauf H et He4
      DO j=1,nchim
       ymaxj(j)=MAXVAL(y(j,:))
       IF(j /= 1 .AND. j /= ihe4)y(j,:)=y(j,:)/ymaxj(j)
      ENDDO
      CALL pgsls(1) ; CALL pgvsize(xleft,xleft+ld,ybot,ybot+h)
      CALL pgwindow(xmin,xmax,ymin,ymax)

      CALL pgsch(0.8)   !hauteur des caractères
      ic=1
      DO j=1,nchim
       WRITE(htext(j),10)ymaxj(j)
10     FORMAT(es10.3)   
       tex(j)=nom_elem(j)//'='//htext(j)
       ic=ic+1 ; IF(ic > 13)ic=2 ; IF(ic == 4)ic=5 ; CALL pgsci(ic)
       yd(:)=y(j,:) ; CALL pgline(n_des,x,yd)
       CALL pgtext(xloc,yloc(j),tex(j))
      ENDDO
      CALL pgsci(1)
      CALL pgsch(1.)    !hauteur des caractères

c tracé de Rtotp
      CALL pgsls(2) ; CALL pgline(2,xt,rt) ; CALL pgsls(1)
c     WRITE(*,2000)xt(1),xt(2),rt(1),rt(2) ; PAUSE'rtot1'

      IF(change)THEN    !trace du cadre P, T, L, M
       CALL pgvsize(xleft,xleft+ld,ybot+h+dh,ybot+h+dh+h)
       CALL pgwindow(xmin,xmax,ymin,ymax)
       xopt='bcnst' ; xtick=0. ; nxsub=0 ; yopt='bcnst' ; ytick=0. ; nysub=0
       CALL pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
       side='b'   !abscisses
       disp=2 ; coord=0.5 ; fjust=0.5 ; text='R/R\d\(2281)\u'
       CALL pgmtext(side,disp,coord,fjust,text)
       side='l'   !ordonnees
       disp=2 ; coord=0.5 ; fjust=0.5
       CALL pgsch(0.9)
       CALL pgmtext(side,disp,coord,fjust,texptrl)
       CALL pgsch(1.) 
       side='t'   !sommet
       disp=1 ; coord=0.5 ; fjust=0.5 ;
       CALL pgmtext(side,disp,coord,fjust,nom_fich2)
      ENDIF
      
c P, T, R, L, R, (Omega si Krot>=3) en fonction de R
	dq=(n_qs-1.d0)/dble(n_des-1)  !pour les variables
	qx=(/ (dq*(i-1)+1.d0, i=1,n_des) /)

	DO i=1,n_des
	 CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,MIN(qx(i),q(n_qs)),
	1 lq,f,dfdq)
	  IF(no_croiss)PRINT*,'Pb. en 5 dans des_r'
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
	1 MAX(mrott(1),MIN(nu,mrott(knotr))),lq,frot,dfrot)
	  IF(no_croiss)PRINT*,'Pb. en 6 dans des_r'	
	  yw(i)=frot(1)
	 ENDDO
	 vext=yw(n_des)*yr(n_des)*rsol*1.d-5
	 wmax=MAXVAL(yw) ; yw=yw/wmax
	ENDIF

c normalisations
	yp=EXP(yp) ; yt=EXP(yt)
	pmax=MAXVAL(yp) ; tmax=MAXVAL(yt)
	lmax=MAXVAL(yl) ; mmax=MAXVAL(ym)
	yp=yp/pmax ; yt=yt/tmax ; yl=yl/lmax ; ym=ym/mmax

	tli(1)=LOG(lib)/tmax ; tli(2)=tli(1)

	CALL pgvsize(xleft,xleft+ld,ybot+h+dh,ybot+h+dh+h)    !P, T, L, M
	CALL pgwindow(xmin,xmax,ymin,ymax)

	IF(Krot >= 3)THEN
	 x_age=(4.*xmin+xmax)/5.
	ELSE	 
	 x_age=(2.*xmin+xmax)/3.
	ENDIF
	y_age=ymax+0.2

c écriture de l'entête	
	IF(Krot >= 3)THEN
	 WRITE(text_age,3)model_num,age,dt,vext
	 CALL pgtext(x_age,y_age,text_age)
3	 FORMAT('modèle ',i4.4,', âge en My: ',es10.3,', dt=',es10.3,
	1 ', V\dext\u=',es10.3,'Km/s')
	 CALL pgtext(x_age,y_age,text_age)	
	ELSE
	 WRITE(text_age,9)model_num,age,dt
	 CALL pgtext(x_age,y_age,text_age)
9	 FORMAT('modèle ',i4.4,', âge en My: ',es10.3,', dt=',es10.3)
	 CALL pgtext(x_age,y_age,text_age)
	ENDIF

	CALL pgsch(0.8)   !hauteur des caractères
	WRITE(char,10)pmax ; ptex='P\dmax\u='//char ; CALL pgsci(2)
	CALL pgline(n_des,yr,yp) ; CALL pgtext(xloc,ylocp(1),ptex)
	WRITE(char,10)tmax ; ttex='T\dmax\u='//char ; CALL pgsci(3)
	CALL pgline(n_des,yr,yt) ; CALL pgtext(xloc,ylocp(2),ttex)
	WRITE(char,10)lmax ; ltex='L\dmax\u='//char ; CALL pgsci(6)
	CALL pgline(n_des,yr,yl) ; CALL pgtext(xloc,ylocp(3),ltex)
	WRITE(char,10)mmax ; mtex='M\dmax\u='//char ; CALL pgsci(7)
	CALL pgline(n_des,yr,ym) ; CALL pgtext(xloc,ylocp(4),mtex)
	teffs=teff ; WRITE(char,10)teffs ; tefftex='T\deff\u='//char
	CALL pgsci(9) ;	CALL pgtext(xloc,ylocp(5),tefftex)	
	IF(Krot >= 3)THEN
	 WRITE(char,10)wmax ; wtex='\gW\dmax\u='//char ; CALL pgsci(8)
	 CALL pgline(n_des,yr,yw) ; CALL pgtext(xloc,ylocw,wtex)
	ENDIF	   
	
	CALL pgsci(1) ; CALL pgsch(1.)    !hauteur des caractères

c tracé de Rtotp
	CALL pgsls(2) ; CALL pgline(2,xt,rt) ; CALL pgsls(1)
c	WRITE(6,2000)xt(1),xt(2),rt(1),rt(2) ; PAUSE'rtot2'

c diagramme HR
	lteff=LOG10(teff) ; lq=n_qs
	CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(n_qs),lq,f,dfdq)
	IF(no_croiss)PRINT*,'Pb. en 7 dans des_r'      
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
c     PRINT*,lim,(lconv(i),i=1,lim)
c     PRINT*,(m_zc(i),i=1,lim),age

      pas=pas+1. ; x4(1)=pas ; x4(2)=pas

      IF(lim > 0)THEN   !lim=0: modele totalement radiatif pas de dessin      
       CALL pgvsize(xleft+ld+dl,xleft+2.*ld+dl,ybot+h+dh,ybot+2.*h+dh)

c      PRINT*,lim,(lconv(i),i=1,lim)
c      WRITE(6,2000)(m_zc(i),i=1,lim),mstar,mtot

       CALL pgwindow(xmin4,xmax4,ymin4,ymax4)
       IF(lim == 1 .AND. .NOT.lconv(1) .AND.
     1 m_zc(1)/mstar .ge. .99)THEN       !modele totalement convectif
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
           yd(1)=0 ; yd(2)=m_zc(1) ; CALL pgline(2,x4,yd)
          ELSE
           PRINT*,'probleme' ; CALL pgend ; STOP
          ENDIF
         ENDIF
         i=i+1
        ENDDO     !WHILE
       ENDIF
      ENDIF

      IF(fin)CALL pgend

      RETURN

      END SUBROUTINE des_r
