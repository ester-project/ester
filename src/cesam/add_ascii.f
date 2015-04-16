	
c********************************************************************

	SUBROUTINE add_ascii(var,glob,itot,ivar)
			
c routine private des modules mod_cesam et mod_exploit

c complète le tableau des quantités entrant éventuellement dans
c les fichiers ASCII de type output, par exemple pour les
c oscillations adiabatiques

c entrées :
c	glob :  globales
c	itot : nombre total de points
c	ivar : nombre de variables	

c entrées/sorties :
c	var : tableau des variables

c	les quantités sont les suivantes:

c	var(1,i): Rayon
c	var(2,i): Ln M/Mtot
c	var(3,i): Température
c	var(4,i): Pression
c	var(5,i): Densité
c	var(6,i): Gradient
c	var(7,i): Gradient
c	var(8,i): Luminosité
c	var(9,i): Opacité
c	var(10,i): Energie nuc+grav
c	var(11,i): Grand Gamma1
c	var(12,i): Gradient adiabatique
c	var(13,i): Delta
c	var(14,i): Cp
c	var(15,i): Mue^(-1)
c	var(16,i): Mu      
c	var(17,i): Vaissala
c	var(18,i): Omega
c	var(19,i): dln kappa/dln T
c	var(20,i): dln kappa/dln ro
c	var(21,i): d epsilon(nuc) / d ln T
c	var(22,i): d epsilon(nuc) / d ln ro
c	var(23,i): !Ptot/Pgaz ou grad_mu
c	var(24,i): !Gradient radiatif       
c	var(25,i): d Gamma1 / d lnP (TY)
c	var(26,i): d Gamma1 / d lnT (PY)
c	var(27,i): d Gamma1 / dY	(PT)
c	var(28,i): dP / dro (TX)
c	var(29,i): dP / dT (roX)
c	var(30,i): dP / dX (Tro)
c	var(31,i): du / dro (TX)
c	var(32,i): du / dT (roX)
c	var(33,i): du / dX(Tro)
c	var(34,i): énergie interne
c	var(35,i): d^2P / dro^2 (TX)
c	var(36,i): d^2P / dro dT (X)	 
c	var(37,i): d^2P / dT^2(roX)		 	 
c	var(38,i): d^2U / dro^2 (TX)	 
c	var(39,i): d^2U / dro dT (X) 
c	var(40,i): d^2U / dT^2 (X) 	  
c	var(41,i): dK / dX
c	var(42,i): d^2K / dT^2	
c	var(43,i): d epsi / dX
c	var(44,i): dX / dR
c	var(45,i): J-B	  
c	var(46,i): Edding. facteur
c	var(ivar+j,i): xchim1g(j=1,nchim) Abondances / gramme

c	i=1 correspond à la surface, i=itot au centre

c	Dans les routines osc_adia, invers, etc... on constitue un 
c	fichier ASCII en prenant les quantités désirées parmi les
c	précédentes

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : aradia, clight, g, nchim, nucleo
	USE mod_etat, ONLY : etat, df_rotx
	USE mod_kind
	USE mod_nuc, ONLY : nuc
	USE mod_numerique, ONLY : bsp1dn, no_croiss, sum_n
	USE mod_opa, ONLY : opa
	
	IMPLICIT NONE	

	INTEGER, INTENT(in) :: itot, ivar
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: glob
	REAL (kind=dp), INTENT(inout), DIMENSION(:,:) :: var
		
	REAL (kind=dp), PARAMETER :: dd=1.d-5, unpdd=1.d0+dd
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: frad,  jac, x	
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: comp, dcomp,
	1 df, epsilon, ex, f, rx, rxt, sum, xchim, xchimd
	
	REAL (kind=dp) :: p, t, ro, drodp_tx, drodt_px, drodx_pt,
	1 u, dudp_tx, dudt_px, dudx_pt, m, d_p, cte1, cte2, cte3,
	2 delta, deltap, deltat, deltax, cp, dcpp, dcpt, dcpx,
	3 gradad, dgradadp, dgradadt, dgradadx, alfa, beta,
	4 gamma1, gamma, dgamdlp, dgamdlt, f_edd, b_rad, j_rad,
	5 dlrodx, dgamdx, d2pdt2, dpdxt, dudxt, dudtro,
	6 d2udt2, dudtt, dpdtt, dudrot, dx, d_t, t1, dpdrot,
	7 dudtp, d2pdro2, d2udro2, alfa1, dpdtp, dudrop, delta1, u1,
	8 dpdrop, dudx, p1, ro1, dpdx, dudro, dudt, dpdro, dpdt,
	9 kappa, dkapdt, dkapdr, dkapdx, hh, be7 ,b8,
	1 n13, o15, f17, et, ero, dkapdx1, dkapdt1, d2kdt2,
	2 dpdxp, dudxp, d2pdrot, d2udrot
	 
	INTEGER :: i, knot, l, mx, j
    
c------------------------------------------------------------------------

2000	FORMAT(8es10.3)

	cte1=4.d0*aradia*clight/3.d0 ; cte2=cte1*g
	cte3=aradia*clight/4.d0

	ALLOCATE(xchim(nchim),xchimd(nchim),comp(nchim),dcomp(nchim),
	1 epsilon(4),ex(nchim),jac(nchim,nchim))	
	
c Tabulation de X=f(r) pour calcul de la dérivée dX/dr
c de x(2)=Frad=f(r) pour calcul de la dérivée dFrad/dr=B-J
c on identifie Fs=Frad pour le calcul de la
c somme[0,r] 3/4k ro Frad dr=fJ(r) car Frad est modifié
c mx=2 ie. interpolation linéaire
	mx=2		 !; mx=3 !; mx=4
	ALLOCATE(x(2,itot),rx(itot),rxt(itot+2*mx-2),f(2),df(2),
	1 frad(1,itot),sum(1))	
	DO i=1,itot
	 j=itot-i+1
	 IF(var(2,i) > -1.d30)THEN
	  m=exp(var(2,i))*glob(1)
	  frad(1,j)=cte2/var(8,i)/var(4,i)/var(1,i)**2*var(3,i)**4*
	1 m*var(6,i)
	 ELSE
	  m=0.d0 ; frad(1,j)=0.d0  
	 ENDIF
	 x(1,i)=var(ivar+1,j)		!pour X	
	 x(2,i)=frad(1,j)		!pour Frad
	 rx(i)=var(1,j)			!pour R
	ENDDO
	
c	DO i=1,itot
c	 j=itot-i+1
c	 WRITE(*,2000)rx(i),var(1,j),x(1,i),var(ivar+1,j),x(2,i),frad(i)
c	ENDDO
c	PAUSE'x'	

	CALL bsp1dn(2,x,rx,rxt,itot,mx,knot,.FALSE.,rx(1),l,f,df)
	IF(no_croiss)THEN
	 PRINT*,'Arrêt 1 dans add_ascii' ; STOP
	ENDIF

c Tabulation de Fsum=Frad=f(r) pour calcul de la dérivée dFrad/dr=B-J
c et de la somme[0,r] 3/4 k ro Frad dr=fJ(r) on utilise le même vecteur nodal
	DO i=1,itot
	 j=itot-i+1
	 frad(1,j)=frad(1,j)*3.d0/4.d0*var(8,i)*var(5,i)
	ENDDO
	CALL bsp1dn(1,frad,rx,rxt,itot,mx,knot,.FALSE.,rx(1),l,f,df)
	IF(no_croiss)THEN
	 PRINT*,'Arrêt 2 dans add_ascii' ; STOP
	ENDIF		 
	CALL sum_n(1,frad,rxt,mx,knot,.FALSE.,rx(itot),rx(itot-10),sum)

c le fichier d'oscillations augmenté
	B1: DO i=1,itot
c	B1: DO i=itot-5,itot	
c	B1: DO i=10,10
c	B1: DO i=430,450	
	 p=var(4,i) ; t=var(3,i)
	 xchim(1:nchim)=var(ivar+1:ivar+nchim,i)
c	 WRITE(*,2000)xchim ; PAUSE'xchim'
	 
	 CALL etat(p,t,xchim,.TRUE.,
	1 ro,drodp_tx,drodt_px,drodx_pt,u,dudp_tx,dudt_px,dudx_pt,
	2 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma)
	
	 CALL df_rotx(p,ro,t,alfa,delta,1.d0,0.d0,0.d0,
	1 drodx_pt,dpdro,dpdt,dpdx)
	
	 CALL df_rotx(p,ro,t,alfa,delta,dudp_tx,dudt_px,dudx_pt,
	1 drodx_pt,dudro,dudt,dudx)

c dérivées / p
	 p1=p*unpdd ; d_p=p1-p
	 CALL etat(p1,t,xchim,.TRUE.,
	1 ro1,drodp_tx,drodt_px,drodx_pt,u1,dudp_tx,dudt_px,dudx_pt,
	2 delta1,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3 gradad,dgradadp,dgradadt,dgradadx,alfa1,beta,gamma1)
	
	 dgamdlp=(gamma1-gamma)/d_p*p
	 
c dérivées secondes
	 CALL df_rotx(p1,ro1,t,alfa1,delta1,1.d0,0.d0,0.d0,drodx_pt,
	1 dpdrop,dpdtp,dpdxp)
	 dpdrop=(dpdrop-dpdro)/d_p	![dP/dro{P+dP  - dP/dro{P] /dP
	 dpdtp=(dpdtp-dpdt)/d_p	![dP/dT{P+dP  - dP/dT{P] /dP
	 CALL df_rotx(p,ro,t,alfa,delta,dpdrop,dpdtp,0.d0,drodx_pt,
	1 d2pdro2,d2pdrot,dpdxp)
	  
	 CALL df_rotx(p1,ro1,t,alfa1,delta1,dudp_tx,dudt_px,dudx_pt,
	1 drodx_pt,dudrop,dudtp,dudxp)
	 dudrop=(dudrop-dudro)/d_p	
	 dudtro=(dudtp-dudt)/d_p
	 CALL df_rotx(p,ro,t,alfa,delta,dudrop,dudtro,0.d0,drodx_pt,
	1 d2udro2,d2udrot,dudxp)	  

c dérivées / t
	 t1=t*unpdd ; d_t=t1-t
	 CALL etat(p,t1,xchim,.TRUE.,
	1 ro1,drodp_tx,drodt_px,drodx_pt,u1,dudp_tx,dudt_px,dudx_pt,
	2 delta1,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3 gradad,dgradadp,dgradadt,dgradadx,alfa1,beta,gamma1)
	 dgamdlt=(gamma1-gamma)/d_t*t
	 
c dérivées secondes
	 CALL df_rotx(p,ro1,t1,alfa1,delta1,1.d0,0.d0,0.d0,
	1 drodx_pt,dpdrot,dpdtt,dpdxt)
	 dpdrot=(dpdrot-dpdro)/d_t	
	 dpdtt=(dpdtt-dpdt)/d_t
	 CALL df_rotx(p,ro,t,alfa,delta,dpdrot,dpdtt,0.d0,
	1 drodx_pt,d2pdt2,d2pdrot,dpdxp)
	  
	 CALL df_rotx(p,ro1,t,alfa1,delta1,dudp_tx,dudt_px,dudx_pt,
	1 drodx_pt,dudrot,dudtt,dudxt)
	 dudrot=(dudrot-dudro)/d_t	
	 dudtt=(dudtt-dudt)/d_t
	 CALL df_rotx(p,ro,t,alfa,delta,dudrot,dudtt,0.d0,
	1 drodx_pt,d2udrot,d2udt2,dudxt)

c X+Y+Z=1 ==> dX+dY+dZ=0		
c dérivées(Z=cte) / Y = - dérivées(Z=cte) / X	
c dérivées(Y=cte) / Z = - dérivées(Y=cte) / X
c dgam_dy(dZ=0)=-dgam_dx(dZ=0)=dgam_dz(dY=0)
	 xchimd=xchim ; xchimd(1)=xchimd(1)*unpdd ; dx=xchimd(1)-xchim(1)
	 CALL etat(p,t,xchimd,.TRUE.,
	1 ro1,drodp_tx,drodt_px,drodx_pt,u1,dudp_tx,dudt_px,dudx_pt,
	2 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	
	 dgamdx=(gamma1-gamma)/dx ; dlrodx=(ro1-ro)/dx/ro

c les dérivées du Gamma1
	 var(25,i)=dgamdlp+alfa/delta*dgamdlt	 !d Gamma1 / d lnP (TY)
	 var(26,i)=-dgamdlt/delta		 !d Gamma1 / d lnT (PY)
	 var(27,i)=-dgamdx-dlrodx*dgamdlt/delta !d Gamma1 / dY	(PT)
	 
c	 WRITE(*,2000)var(25:27,i)

c opacités 
	 CALL opa(xchim,t,var(5,i),kappa,dkapdt,dkapdr,dkapdx)
	 t1=t*unpdd ; d_t=t1-t
	 CALL opa(xchim,t1,var(5,i),kappa,dkapdt1,dkapdr,dkapdx1)	 
	 d2kdt2=(dkapdt1-dkapdt)/d_t

c réactions thermonucléaires
	 comp(1:nchim)=xchim(1:nchim)/nucleo(1:nchim)
c	 WRITE(*,2000)t,var(5,i),comp,xchim,nucleo
	 CALL nuc(t,var(5,i),comp,dcomp,jac,.TRUE.,3,
	1 epsilon,et,ero,ex,hh,be7,b8,n13,o15,f17)
c	 WRITE(*,2000)et,ero
c	 WRITE(*,2000)ex

c facteurs d'Eddington
	 CALL bsp1dn(2,x,rx,rxt,itot,mx,knot,.TRUE.,var(1,i),l,f,df)
	 IF(no_croiss)PRINT*,'En 3 dans add_ascii'
	 b_rad=cte3*var(3,i)**4
	 IF(var(8,i) > 0.d0)THEN 
	  CALL sum_n(1,frad,rxt,mx,knot,.TRUE.,var(1,1),var(1,i),sum)
	  j_rad=b_rad-df(2)/4.d0/var(5,i)/var(8,i)
	  f_edd=-sum(1)/j_rad		!MIN(1.d0,-sum(1)/j_rad)
	 ELSE
	  j_rad=b_rad ; f_edd=1.d0		!au centre
	 ENDIF

c les autres quantités 
	 var(28,i)=dpdro	!dP / dro (TX)
	 var(29,i)=dpdt 	!dP / dT (roX)
	 var(30,i)=dpdx 	!dP / dX (Tro)
	 var(31,i)=dudro	!du / dro (TX)
	 var(32,i)=dudt		!du / dT (roX)
	 var(33,i)=dudx	 	!du / dX(Tro)
	 var(34,i)=u 		!énergie interne
	 var(35,i)=d2pdro2 	!d^2P / dro^2 (TX)
	 var(36,i)=d2pdrot	!d^2P / dro dT (X)	 
	 var(37,i)=d2pdt2 	!d^2P / dT^2(roX)		 	 
	 var(38,i)=d2udro2 	!d^2U / dro^2 (TX)	 
	 var(39,i)=d2udrot	!d^2U / dro dT (X) 
	 var(40,i)=d2udt2	!d^2U / dT^2 (X) 	  
	 var(41,i)=dkapdx	!dK / dX
	 var(42,i)=d2kdt2	!d^2K / dT^2	
	 var(43,i)=ex(1)	!d epsi / dX
	 var(44,i)=df(1)	!dX / dR
	 var(45,i)=j_rad-b_rad	!J-B	  
	 var(46,i)=MAX(0.d0,MIN(1.d0,f_edd))	!Edding. facteur

c	 PRINT*,i,itot-i+1
c	 WRITE(*,2000)var(1,1),var(1,i),f(1),df(1),f(2),df(2)
c	 WRITE(*,2000)sum,var(43,i),var(44,i),b_rad,j_rad,f_edd
	 
c	 WRITE(*,2000)var(1,i),var(23:25,i)
c	 WRITE(*,2000)var(26:32,i)
c	 WRITE(*,2000)var(1,i),var(33:38,i)
c	 WRITE(*,2000)var(39:44,i)	 	 	 
c	 PAUSE'interpolations'	 	 
	ENDDO B1
	
	DEALLOCATE(x,xchim,xchimd,comp,dcomp,epsilon,ex,jac,rx,rxt,
	1 f,df,frad)
	
	RETURN

	END SUBROUTINE add_ascii
