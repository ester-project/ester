
c***********************************************************************

	SUBROUTINE etat_eff(p,tt,xchim,ro,drop,drot,drox,u,dup,dut,dux,
	1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)

c       subroutine private du module mod_etat

c       adaptation du logiciel eff écrit par jorgen christensen_dalsgaard
c       à partir de p et t, il calcule ro et u et toutes les dérivées
c       jusqu'au deuxième ordre
c       Annie Baglin le 1 12 1990

c       modifs:
c	14 07 90 : suppression de IMPLICIT REAL*8 (a-h,o-z)
c		  déclaration explicite des variables, Michel Auvergne.
c	25 12 90 : ajout de lamb=5, suppression des variables inutiles
c	01 10 91 : appel à etat_gong2 en cas d'échec de convergence
c	 dans eqstp
c	19 11 99 : suppression des nh1, nh2, nhe2, lamb

c	Adaptation F95: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c entree :
c	p : pression
c	tt : température
c	xchim : composition chimique par gramme

c sortie :
c	ro : densité et dérivées
c	u : énergie interne et dérivées

c----------------------------------------------------------------

	USE mod_donnees, ONLY : aradia, ah, amu, ahe4, ihe4, nchim, Lz0 => z0
	USE mod_kind

	IMPLICIT NONE
		
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim	
	REAL (kind=dp), INTENT(in) :: p, tt	
	REAL (kind=dp), INTENT(out) :: ro,drop,drot,drox,u,dup,dut,dux,
	1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1

c~~~~	
c	initialisations pour remplacement du COMMON potetc initialise
c	par le BLOCKDATA bleqst	
	
	INTEGER :: i
	REAL (kind=dp), PARAMETER, DIMENSION(374) :: chi=(/
	1 11.26d0, 24.38d0, 47.86d0, 64.48d0, 391.99d0,489.84d0,     
	2 14.54d0, 29.60d0, 47.43d0, 77.45d0, 97.86d0, 551.92d0, 666.83d0,
	3 13.61d0, 35.15d0, 54.93d0, 77.39d0, 113.87d0,138.08d0, 739.11d0,
	4 871.12d0,21.56d0, 41.07d0, 63.5d0 , 97.16d0, 126.4d0,  157.91d0,
	5 207.3d0, 239.d0,  1196.d0, 1360.d0, 5.14d0,  47.29d0,  71.65d0,
	6 98.88d0, 138.60d0,172.36d0,208.44d0,264.15d0,299.78d0, 1465.d0,
	7 1646.d0, 7.64d0,  15.03d0 ,80.12d0 ,109.29d0,141.23d0, 186.86d0,
	8 225.31d0,265.96d0,327.90d0,367.36d0,1761.2d0,2085.d0,  5.98d0,
	9 18.82d0, 28.44d0, 119.96d0,153.77d0,190.42d0,241.93d0, 285.13d0,
	1 330.1d0, 398.5d0, 441.9d0, 2085.5d0,2299.d0, 8.15d0,   16.34d0,
	2 33.46d0, 45.13d0, 166.73d0,205.11d0,246.41d0,303.87d0, 351.83d0,
	3 401.3d0, 476.d0,  523.2d0, 2436.d0, 2666.d0, 15.75d0,  27.62d0,
	4 40.90d0, 59.79d0, 75.d0,   91.3d0,  124.d0,  143.46d0, 421.d0,
	5 480.d0,  539.5d0, 621.1d0, 688.5d0, 755.5d0, 854.4d0,  1.d20,
	6 1.d20,   1.d20,   7.90d0,  16.18d0, 30.64d0, 56.d0,    79.d0,
	7 105.d0,  133.d0,  151.d0,  235.d0,  262.d0,  290.d0,   321.d0,
	8 355.d0,  390.d0, 457.d0, (1.d20, i=115,125), (88888.d0,i=126,374)
	9 /)
	REAL (kind=dp), PARAMETER, DIMENSION(23) :: am=(/ 12.d0, 14.01d0,
	1 16.d0, 20.17d0, 22.99d0, 24.31d0, 26.98d0, 28.08d0,       
	2 39.94d0, 55.84d0, (99999.d0,i=11,23) /)		
	INTEGER, PARAMETER, DIMENSION(23) :: iz=(/ 6, 7, 8, 10, 11, 12,
	1 13, 14, 18, 26, (9999,i=11,23) /)
	CHARACTER (len=4), PARAMETER, DIMENSION(23) :: name=(/ '   c',
	1 '   n','   o','  ne','  na','  mg','  al','  si','   a','  fe',
	2 ('bido',i=11,23) /)
c~~~~~	

	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: Lxchim
	REAL (kind=dp) :: drott,drotp,drotx,dropl,dropn,drotn,droxn,drofn,
	1 f,dfpn,dftn,dfxn,dfpl,dftl,dfxl,dptn,Lp,Lt,
	2 droffn,droftn,drottn,drotxn,drofxn,drotpn,
	3 dpffn,dptxn,dpfxn,betgam,betbet,dpttn,dpxn, p10,t10,
	4 unsro,ro2,psro2,psro3,cpp,dutn,duxn,dupn,dutt,dutx,dutp,
	5 duttn,dutpn,dutxn,dhtn,dhxn,dhtpn,dhtxn,dhttn,
	6 dhffn,dhftn,dhfxn,dhpn,dhfn,dpftn,aradias3	

	LOGICAL :: ok
	COMMON /marche/ok

	LOGICAL, SAVE :: init=.TRUE.

	REAL (kind=dp) :: av,ah_eff,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
	1 ct,crho,cpe,che,cao3,caa,ckh,car,ergev
	COMMON/consts/av,ah_eff,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
	1 ct,crho,cpe,che,cao3,caa,ckh,car,ergev
	
	REAL (kind=dp) :: xii1(4),ane(10),rho(20),ht(20),pt(20),
	1 cp_eff(4),dad(4), dlt(4),gm1,tprh,trph,rhxp
	COMMON/eqstd/xii1,ane,rho,ht,pt,cp_eff,dad,
	1 dlt,gm1,tprh,trph,rhxp
	
	INTEGER ihvz
	REAL (kind=dp) :: anh0
	COMMON/eqscnt/anh0,ihvz
	
	REAL (kind=dp) :: amm,amm2,amm3
	COMMON/ln10/amm,amm2,amm3
	
c----------------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(init)THEN
	 CALL setcns		!initialisation des ctes
	 anh0=0.5d0		!initialisation ???????
	 init=.FALSE. ; WRITE(2,1) ; WRITE(*,1)
1	 FORMAT('équation d''état EFF')
	 aradias3=aradia/3.d0 ; ALLOCATE(Lxchim(nchim))
	ENDIF
	
	Lxchim(1:nchim)=xchim(1:nchim) ; ok=.TRUE. ; p10=LOG10(p) ; t10=LOG10(tt)
c	PRINT*,ihe4
c	WRITE(*,2000)p10,t10,xchim(1),xchim(ihe4),x0,y0,Lz0
	IF(ihe4 > 1)THEN
	 CALL theffp(p10,t10,Lxchim(1),Lxchim(ihe4)+Lxchim(ihe4-1),f)
	ELSE                
	 CALL theffp(p10,t10,Lxchim(1),1.d0-Lxchim(1)-Lz0,f)
        ENDIF
         
	IF(ok)THEN

c        sorties

c        ro et ses dérivées

         ro=rho(1)

c        dérivées premieres de rho

	 drofn=rho(2)*rho(1)/f
         dropl=rho(2)/pt(2)
         dropn=(rho(1)/pt(1))*dropl
         drotn=(rho(1)/tt)*(rho(3)-dropl*pt(3))
         droxn=rho(1)*(rho(4)-dropl*pt(4))*amm
         drop=dropn
         drot=drotn
         drox=droxn

c        dérivées secondes de rho et p

c        1. dérivées secondes naturelles dans le systeme ftx

	 dpffn=(pt(1)/(f*f))*((pt(5)/amm)+pt(2)*(pt(2)-1.d0))
         droffn=(rho(1)/(f*f))*((rho(5)/amm)+rho(2)*(rho(2)-1.d0))
         drottn=(rho(1)/tt**2)*((rho(8)/amm)+rho(3)*(rho(3)-1.d0))
	 dpttn=(pt(1)/tt**2)*((pt(8)/amm+pt(3)*(pt(3)-1.d0)))
         droftn=(rho(1)/(f*tt))*(rho(2)*rho(3)+rho(6)/amm)
         drotxn=(rho(1)/tt)*(rho(9)+rho(4)*rho(3)*amm)
	 dptxn=(pt(1)/tt)*(pt(9)+pt(4)*pt(3)*amm)
         drofxn=(rho(1)/f)*(rho(7)+rho(4)*rho(2)*amm)
	 dpfxn=(pt(1)/f)*(pt(7)+pt(4)*pt(2)*amm)
         dpftn=(pt(1)/(f*tt))*(pt(2)*pt(3)+pt(6)/amm)

c        2. dérivées premieres naturelles de f dans le systeme ptx

         dfpl=1.d0/pt(2) ; dftl=-pt(3)/pt(2) ; dfxl=-pt(4)/pt(2)
         dfpn=(f/pt(1))*dfpl ; dftn=(f/tt)*dftl ; dpxn=pt(4)*pt(1)*amm
	 dfxn=-dfpn*dpxn ; dfxn=f*dfxl ; dptn=pt(3)*pt(1)/tt

c        3. dérivées secondes naturelles dans le systeme ptx

         drotpn=dfpn*(droftn-dfpn*dptn*droffn+drofn*dfpn*
	1 (dpffn*dptn*dfpn-dpftn))
	 betbet=-(dpttn+2.d0*dftn*dpftn+dftn*dftn*dpffn)
	 betgam=-(dptxn+dfxn*dpftn+dftn*dpfxn+dftn*dfxn*dpffn)
	 drotxn=drotxn+dftn*drofxn+dfxn*droftn+dftn*dfxn*droffn+betgam*dropn
	 drottn=drottn+2.d0*dftn*droftn+dftn*dftn*droffn+betbet*dropn
         drott=drottn ; drotp=drotpn ; drotx=drotxn

c        energie interne  u   ! attention jorgen calcule l'enthalpie

         unsro=1.d0/ro ; u=ht(1)-pt(1)*unsro

c        dérivées premieres de u

         ro2=ro*ro ; psro2=pt(1)/ro2 ; dhfn=ht(2)/(f*amm) ; dhpn=dhfn*dfpn
	 dupn=dhpn-unsro+psro2*dropn ; dhtn=(ht(3)-ht(2)*pt(3)/pt(2))/(tt*amm)
	 dutn=dhtn+psro2*drotn ; dhxn=ht(4)-ht(2)*pt(4)/pt(2)
	 duxn=dhxn+psro2*droxn ; dux=duxn ; dut=dutn ; dup=dupn	

c        dérivées secondes de u

c	 1- dérivées naturelles de h dans le systeme ftx

	 dhttn=((ht(8)/amm)-ht(3))/(tt**2*amm) ; dhtxn=ht(9)/(tt*amm)
	 dhffn=((ht(5)/amm)-ht(2))/(f*f*amm) ; dhftn=ht(6)/(f*tt*amm2)
	 dhfxn=ht(7)/(f*amm)

c	 2- dérivées naturelles de h dans le systeme ptx

	 dhttn=dhttn+2.d0*dftn*dhftn+dftn*dftn*dhffn+betbet*dhpn
	 dhtpn=dfpn*(dhftn-dfpn*dptn*dhffn+dhfn*dfpn*
	1 (dpffn*dptn*dfpn-dpftn))
	 dhtxn=dhtxn+dftn*dhfxn+dfxn*dhftn+dftn*dfxn*dhffn+betgam*dhpn

c	 3- dérivées de u dans le systeme ptx

	 psro3=2.d0*psro2/rho(1)
	 dutpn=dhtpn+drotn/ro2-psro3*dropn*drotn+psro2*drotpn
	 duttn=dhttn-psro3*drotn*drotn+psro2*drottn
	 dutxn=dhtxn-psro3*drotn*droxn+psro2*drotxn+drotn*dptn/ro2
	 dutt=duttn ; dutp=dutpn ; dutx=dutxn

c        degres d'ionisation

c        nh1=xii1(1)
c        nhe1=xii1(2)
c        nhe2=xii1(3)

	 delta=-tt/ro*drot ; deltap=delta*(-drop/ro+drotp/drot)
	 deltat=delta*(-drot/ro+drott/drot+1.d0/tt)
	 deltax=delta*(-drox/ro+drotx/drot)	 

	 cpp=p/ro/tt*delta ; cp=dut+cpp	
	 dcpp=dutp+cpp*(-drop/ro+deltap/delta+1.d0/p)
	 dcpt=dutt+cpp*(-drot/ro+deltat/delta-1.d0/tt)
	 dcpx=dutx+cpp*(-drox/ro+deltax/delta)

	 gradad=p/ro*delta/cp/tt		!gradient adiabatique
	 dgradadp=gradad*(-drop/ro+deltap/delta-dcpp/cp+1.d0/p)
	 dgradadt=gradad*(-drot/ro+deltat/delta-dcpt/cp-1.d0/tt)
	 dgradadx=gradad*(-drox/ro+deltax/delta-dcpx/cp)
	
	 alfa=p/ro*drop ; beta=1.d0-aradias3*tt**4/p
	 gamma1=1.d0/(alfa-delta*gradad)
	 
	ELSE	!appel a ETAT_GONG2 car il y a pb avec EFF pour les cas extremes
	 Lp=p ; Lt=tt
	 PRINT*,'appel a ETAT_GONG2 car defaillance de EFF'	 
         CALL etat_gong2(Lp,Lt,xchim,
	1 ro,drop,drot,drox,u,dup,dut,dux,
	2 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	ENDIF
	
        RETURN
	
	CONTAINS

c*******************************************************************

	SUBROUTINE bilin(x0,y0,z0,x1,y1,z1,x2,y2,z2,x,y,z)  
       
c	performs bilinear interpolation of the fonction f,          
c	given on three arbitray points (x0,y0),(x1,y1),(x2,y2)      
c	the respective fonction values are z0,z1,z2.     
 
c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P.Morel, B.Pichon
c	CESAM2k

c--------------------------------------------------------------------

	USE mod_kind

	REAL (kind=dp), INTENT(in) :: x0, y0, z0, x1, y1, z1, x2, y2, z2,
	1 x, y
	REAL (kind=dp), INTENT(out) :: z
		
	REAL (kind=dp) :: x10, x20, y10, y20, z10, z20, det, dzdx, dzdy
	 
c--------------------------------------------------------------------	
        
	x10=x1-x0 ; x20=x2-x0 ; y10=y1-y0 ; y20=y2-y0 ; z10=z1-z0
	z20=z2-z0      
	det=x10*y20-y10*x20
	IF(det /= 0.d0)THEN
	 dzdx=(z10*y20-z20*y10)/det ; dzdy=(z20*x10-z10*x20)/det  
	 z=z0+(x-x0)*dzdx+(y-y0)*dzdy     
         RETURN
        ELSE
	 WRITE(*,1)x10,x20,y10,y20
1	 FORMAT(/,' collinear points in s/r bilin. error stop.',         
	1 ' x1-x0,x2-x0,y1-y0,y2-y0 = ',/1x,1p4g15.6/)       
	 STOP
	ENDIF
	 
	END SUBROUTINE bilin

c*******************************************************************

	SUBROUTINE eqstf(fl,tl,x,y,z,nosd,notd)             

c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel, B.Pichon
c	CESAM2k

	USE mod_kind

	IMPLICIT REAL*8 (a-h,o-z)

	LOGICAL :: tstl,nosd,notd,cmplio        
	INTEGER :: ihvz,idmu,l,n,ider,iders,k,ia,i,ii,l1,m,m1,jj,j,k1,k2,
	1 kd,imx,i10,lm,nu,nb,kl,nuu,ln,lt,lf,mn
  
      DIMENSION :: phi(30),hst(30),ex(3),ea(30),xi(30),dxt(4),dxt1(4),  
     .  dxt2(4),ane(10),rho(20),dne(10),dph(20),ddmu(10),he(20),pe(10),         
     .   hi(20),pi(10),hh(20),ph(10),hr(20),pr(10),ht(20),pt(20),cp(4),         
     .   dad(4),dlt(4),xii(30),anuh(10),ueh(10),anuhr(23),uehr(23)   
      COMMON/eqscnt/ anh0,ihvz    
      COMMON/ln10/ amm,amm2,amm3  
      COMMON/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,  
     .   ct,crho,cpe,che,ca03,caa,ckh,car,ergev           
      COMMON/dmuder/ dmu,dmuf,dmut,dmux,dmuff,dmuft,dmufx,dmutt,dmutx,          
     .   dmuxx,idmu    
      COMMON/eqsout/ ea,xii,dne,dph,he,pe,hi,pi,hh,ph,hr,pr          
      COMMON/eqstd/ xii1(4),ane,rho,ht,pt,cp,dad,dlt,gm1,tprh,trhp,  
     .   rhxp          
      COMMON/eqsaux/psi           
c         
      EQUIVALENCE(ex(1),exh),(ddmu(1),dmu)   
c         
      fxp(a)= EXP(MIN(8.d1,MAX(a,-1.d2)))             
      nuu(l,n)=((l-1)*(8-l))/2+n-l+5         
c         
c  number of derivatives          
      ider=20          
      IF(notd) ider=10            
      IF(nosd) ider=4  
      iders=MIN(ider,10)         
c         
      f=1.e1**fl ; t=1.e1**tl       
c         
      CALL phder(fl,tl,phi,hst,nosd,notd)    
c         
c         
c  psi and derivatives            
      wf= SQRT(1+f)    
      psi=2*wf+LOG(f/((1+wf)*(1+wf)))       
      psif=amm*wf      
      psff=amm2*f/(2*wf)          
c         
      zf=x+2*y+6*z     
      zf2=zf*zf        
      zf3=zf*zf2       
c         
      ak0=ckh*zf2      
c  k*t, in ev          
      tk=ck1*t         
c  a0**3*ne/(k*t)      
      aa=caa*phi(1)/(tk*zf3)      
c         
c  delta mu and derivatives       
c         
      bmu=tk+20*ak0    
      dmu=aa*bmu       
c         
      ref=phi(2)       
      ret=phi(3)       
      dmuf=dmu*ref*amm            
      dmut=aa*(bmu*ret-20*ak0)*amm           
      dmux=aa*(3*tk+20*ak0)/zf    
c         
      IF(nosd) GOTO 10           
      dmuff=(dmuf*ref+dmu*phi(4)*amm)*amm-psff            
      dmuft=(dmut*ref+dmu*phi(5)*amm)*amm    
      dmufx=dmux*ref*amm          
      dmutt=aa*(20*ak0*(1-2*ret)+bmu*(phi(6)+ret*ret))*amm2          
      dmutx=aa*(ret*(3*tk+20*ak0)-20*ak0)*amm/zf          
      dmuxx=aa*(12*tk+40*ak0)/zf2 
   10 dmu=dmu-psi      
      dmuf=dmuf-psif   
      idmu=1           
c  test for complete ionization of h and he  
      cmplio=(dmu-54.4/tk) > 19  
      IF(cmplio) GOTO 31         
c         
c  e h, e he, e he+ and derivatives          
c         
      k=-10            
      DO 25 ia=1,3     
      k=k+10           
      ext=ex(ia)/tk    
      eea=dmu-ext      
      IF(eea+30 <= 0.d0)THEN
       DO i=1,10	!no ionization
        ea(k+i)=0.d0
       ENDDO     
       GOTO 25
      ENDIF
      eea=fxp(eea)  
      IF(ia-2)22,23,22          
   22 eea=eea/2.d0        
      GOTO 24         
   23 eea=2*eea        
   24 ea(k+1)=eea      
c  first derivatives   
      ea(k+2)=dmuf*eea            
      ea(k+3)=(amm*ext+dmut)*eea  
      ea(k+4)=dmux*eea            
      IF(nosd) GOTO 25           
c  second derivatives  
      ea(k+5)=dmuff*eea+dmuf*ea(k+2)         
      ea(k+6)=dmuft*eea+dmuf*ea(k+3)         
      ea(k+7)=dmufx*eea+dmux*ea(k+2)         
      ea(k+8)=(dmutt-amm2*ext)*eea+(amm*ext+dmut)*ea(k+3)            
      ea(k+9)=dmutx*eea+dmux*ea(k+3)         
      ea(k+10)=dmuxx*eea+dmux*ea(k+4)        
   25 CONTINUE         
c  reset e he+         
      eea=ea(21)       
      IF(eea == 0.d0) GOTO 35       
      eea1=ea(11)      
      ea(21)=eea*eea1  
      ii=24            
      IF(nosd) GOTO 27           
      DO 26 l=1,3      
      l1=l+21          
      DO 26 m=l,3      
      ii=ii+1          
      m1=m+21          
   26 ea(ii)=eea*ea(ii-10)+ea(l1)*ea(m1-10)+ea(l1-10)*ea(m1)+        
     .   eea1*ea(ii)   
   27 DO 28 i=22,24    
   28 ea(i)=ea(i)*eea1+eea*ea(i-10)          
   30 CONTINUE         
      GOTO 35         
c         
c  x h+, x he+, x he++ and derivatives       
c         
c  complete ionization            
   31 xi(1) =1         
      xi(11)=0         
      xi(21)=1         
      DO 33 i=2,22,10  
      jj=i+8           
      DO 33 j=i,jj     
   33 xi(j)=0          
      GOTO 50         
c  partial ionization  
   35 dnm=1+ea(1)      
      xi(1)=ea(1)/dnm  
      ii=4             
      dnm2=dnm*dnm     
      dnm3=dnm*dnm2    
      DO 40 l=1,3      
      l1=l+1           
      eal=ea(l1)       
      xi(l1)=eal/dnm2  
      IF(nosd) GOTO 40           
      DO 38 m=l,3      
      m1=m+1           
      ii=ii+1          
   38 xi(ii)=ea(ii)/dnm2-2*eal*ea(m1)/dnm3   
   40 CONTINUE         
c         
      dnm=1+ea(11)+ea(21)         
      dnm2=dnm*dnm     
      dnm3=dnm*dnm2    
      k1=11            
      k2=21            
   45 kd=k1-k2         
      ii=k1+3          
      ea1=ea(k1)       
      ea2=ea(k2)       
      xi(k1)=ea1/dnm   
      DO 48 l=1,3      
      l1=k1+l          
      eal1=ea(l1)      
      eal2=ea(k2+l)    
      anm=(1+ea2)*eal1-ea1*eal2   
      xi(l1)=anm/dnm2  
      IF(nosd) GOTO 48           
      DO 46 m=l,3      
      ii=ii+1          
      eam1=ea(k1+m)    
      eam2=ea(k2+m)    
   46 xi(ii)=(eal1*eam2+(1+ea2)*ea(ii)-eal2*eam1-ea1*ea(ii-kd))/dnm2 
     .   -2*anm*(eam1+eam2)/dnm3  
   48 CONTINUE         
      IF(k1==21) GOTO 50       
      k1=21            
      k2=11            
      GOTO 45         
c  ionization of heavy elements   
   50 IF(ihvz) 51,51,52           
   51 anuh(1)=anh0     
      anuh(2:10) = 0.    !! BP   CALL zero(anuh(2),9)        
      ueh(:10) = 0.      !! BP   CALL zero(ueh,10)           
      GOTO 53         
   52 CALL hviona(fl,tl,x,y,z,nosd,notd,anuh,ueh,anuhr,uehr)         
c         
c  ne and derivatives  
c         
c  combine he fractions for an in xi(20+.), and for ionization       
c  energy in xi(10+.)  
   53 exhc=exhe+exhep  
c         
      imx=10           
      IF(nosd) imx=4   
      DO 55 i=1,21,10  
   55 CALL store(xi(i),xii(i),imx)           
      xia=xi(1)        
      imx=imx+20       
c         
      DO 56 i=21,imx   
      i10=i-10         
      xio=exhe*xi(i10)+exhc*xi(i) 
      xi(i)=2*xi(i)+xi(i10)       
   56 xi(i10)=xio      
c  terms needed for x-derivatives 
      DO 57 l=1,4      
   57 dxt(l)=av*(xi(l)/ah-xi(20+l)/ahe)      
c         
      ane(1)=av*(xi(1)*x/ah+xi(21)*y/ahe+z*anuh(1))       
c         
      ii=4             
      DO 60 l=1,3      
      l1=l+1           
      anel=(xi(l1)*x/ah+xi(l+21)*y/ahe+z*anuh(l1))*av     
      tstl=l==3      
      IF(tstl) anel=anel+dxt(1)   
      ane(l1)=anel     
      IF(nosd) GOTO 60           
      DO 58 m=l,3      
      m1=m+1           
      ii=ii+1          
      anelm=(xi(ii)*x/ah+xi(20+ii)*y/ahe+z*anuh(ii))*av   
      IF(tstl) anelm=anelm+dxt(m1)           
      IF(m==3) anelm=anelm+dxt(l1)         
   58 ane(ii)=anelm    
   60 CONTINUE         
c         
c  the density and derivatives (note that rho(2) = dlog rho/dlog f,  
c  and so on).         
c         
      anee=ane(1)      
      rho(1)=crho*phi(1)/anee     
      ii=4             
      jj=10            
      DO 63 l=1,3      
      l1=l+1           
      rhol=-ane(l1)/amm/anee      
      tstl=l <= 2      
      IF(tstl) rhol=rhol+phi(l1)  
      rho(l1)=rhol     
      IF(nosd) GOTO 63           
      DO 62 m=l,3      
      ii=ii+1          
      m1=m+1           
      lm=l+m           
      rholm=(ane(l1)*ane(m1)/anee-ane(ii))/anee/amm       
      IF(tstl.and.m <= 2) rholm=rholm+amm*phi(2+lm)       
      IF(notd) GOTO 62           
      DO 61 n=m,3      
      jj=jj+1          
      rhd=-2*ane(l1)*ane(m1)*ane(n+1)/amm/anee**3         
      IF(l < 3.and.m < 3.and.n < 3) rhd=rhd+amm2*phi(4+lm+n)      
   61 rho(jj)=rhd      
   62 rho(ii)=rholm    
   63 CONTINUE         
      pt(:20) = 0.    !! BP   CALL zero(pt,20)            
      ht(:20) = 0.    !! BP   CALL zero(ht,20)            
c         
c         
c  delta p, delta h and derivatives.         
c         
c  test for complete ionization   
      IF(cmplio) GOTO 80         
c         
c  delta ne**2         
c         
c  reset xi(1) and xi(21) (note that xi(2-10) and xi(22-30) are still           
c  consistent).        
      xi(1)=-1/(ea(1)+1)          
      xi(21)=-(2+ea(11))/(1+ea(11)+ea(21))   
c         
      DO 65 l=1,4      
      dxt1(l)=xi(l+20)/ahe-xi(l)/ah          
      dxtl=-xi(l+20)*y/ahe-xi(l)*x/ah-anuh(l)*z           
      IF(l==1) dxtl=dxtl+anh0*z 
      dxtl2=ane(l)     
      IF(l <= 3) GOTO 64         
      dxtl=dxtl+dxt1(1)           
      dxtl2=dxtl2+avda            
   64 dxt(l)=dxtl*av   
   65 dxt2(l)=dxtl2    
c         
      ann=ane(1)+(x/ah+2*y/ahe+z*anh0)*av    
      xtt=dxt(1)       
      dne(1)=ann*xtt   
c         
      ii=4             
      DO 70 l=1,3      
      l1=l+1           
      dne(l1)=dxt2(l1)*xtt+ann*dxt(l1)       
      IF(nosd) GOTO 70           
      tstl=l==3      
      DO 68 m=l,3      
      m1=m+1           
      ii=ii+1          
      dnlm=-xi(20+ii)*y/ahe-xi(ii)*x/ah      
      IF(tstl) dnlm=dnlm+dxt1(m1) 
      IF(m==3) dnlm=dnlm+dxt1(l1)          
   68 dne(ii)=ane(ii)*xtt+dxt2(l1)*dxt(m1)+dxt2(m1)*dxt(l1)          
     .   +ann*dnlm*av  
   70 CONTINUE         
c  quantities COMMON to delta p and delta h (dph(15-20) is used as   
c  intermediate storage).         
      a03=ca03/zf3/2   
      dph(15:20) = 0.  !! BP   CALL zero(dph(15),6)        
      dxt1(1)=0        
      c1=amm*tk        
      dxt1(2)=c1       
      dph(18)=amm2*tk  
      dph(19)=3*c1/zf  
      dnee=dne(1)      
      a03=a03*rho(1)*rho(1)       
      nu=2             
      nb=20            
      k1=1             
   75 DO 77 l=2,4      
   77 dxt(l-1)=dne(l)+nu*amm*dnee*rho(l)     
c         
      dxt1(3)=(3*tk+nb*ak0)/zf    
      bmu=tk+nb*ak0    
      dph(k1)=a03*bmu*dnee        
      ii=k1+3          
      jj=4             
      DO 79 l=1,3      
      l1=l+1           
      dph(k1+l)=a03*(bmu*dxt(l)+dnee*dxt1(l))             
      IF(nosd) GOTO 79           
      DO 78 m=l,3      
      ii=ii+1          
      jj=jj+1          
      m1=m+1           
      dphlm=bmu*(dne(jj)+nu*amm*(dne(l1)*rho(m1)+dne(m1)*rho(l1)     
     .   +dnee*(rho(jj)+nu*amm*rho(m1)*rho(l1))))+dxt(l)*dxt1(m)     
     .   +dxt(m)*dxt1(l)+dnee*dph(10+jj)     
      IF(m==3.and.l==3) dphlm=dphlm+dnee*(12*tk+2*nb*ak0)/zf2    
   78 dph(ii)=a03*dphlm           
   79 CONTINUE         
      IF(nu==1) GOTO 90        
      nu=1             
      nb=40            
      a03=a03/rho(1)   
      k1=11            
      GOTO 75         
c  complete ionization            
   80 dne(:10) = 0.   !! BP  CALL zero(dne,10)           
      dph(:20) = 0.   !! BP  CALL zero(dph,20)           
      GOTO 100        
   90 CONTINUE         
      DO 95 i=1,iders  
      pt(i)=pt(i)+dph(i)          
   95 ht(i)=ht(i)+dph(i+10)       
c  electron pressure and enthalpy 
  100 ii=4             
      jj=10            
      pee=cpe*phi(11)  
      pt(1)=pt(1)+pee  
      hsst=hst(1)      
      hee=che*anee*hsst           
      ht(1)=ht(1)+hee  
      pe(1)=pee        
      he(1)=hee        
      DO 110 l=1,3     
      l1=l+1           
      hll=0            
      hel=hsst*ane(l1)            
      pel=0            
      tstl=l==3      
      IF(tstl) GOTO 102          
      pel=amm*pee*phi(11+l)       
      hll=hst(l1)      
      hel=hel+anee*hll            
      pt(l1)=pt(l1)+pel           
  102 pe(l1)=pel       
      hel=che*hel      
      ht(l1)=ht(l1)+hel           
      he(l1)=hel       
      IF(nosd) GOTO 110          
      DO 108 m=l,3     
      ii=ii+1          
      m1=m+1           
      pelm=0           
      helm=ane(ii)*hsst+hll*ane(m1)          
      tstl=tstl.or.m==3         
      IF(tstl) GOTO 104          
      lm=l+m           
      pelm=amm2*pee*(phi(11+m)*phi(11+l)+phi(12+lm))      
      pt(ii)=pt(ii)+pelm          
      helm=helm+ane(l1)*hst(m1)+anee*hst(2+lm)            
  104 pe(ii)=pelm      
      helm=che*helm    
      ht(ii)=ht(ii)+helm          
      he(ii)=helm      
      IF(notd) GOTO 108          
      DO 106 n=m,3     
      jj=jj+1          
      helm=0           
      IF(tstl) GOTO 106          
      helm=ane(n+1)*hst(2+lm)     
      IF(n==3) GOTO 105        
      helm=helm+ane(l1)*hst(2+m+n)+ane(m1)*hst(2+l+n)     
     .   +anee*hst(4+lm+n)        
      pt(jj)=pt(jj)+amm3*pee*(phi(11+l)*phi(11+m)*phi(11+n)          
     .   +phi(12+lm)*phi(11+n)+phi(12+l+n)*phi(11+m)      
     .   +phi(12+m+n)*phi(11+l)+phi(14+l+m+n))            
  105 helm=che*helm    
      ht(jj)=ht(jj)+helm          
  106 he(jj)=helm      
  108 CONTINUE         
  110 CONTINUE         
c  ionization enthalpy            
c  (pi is introduced for consistency)        
      pi(:10) = 0.   !! BP  CALL zero(pi,10)            
      hi(11:20) = 0. !! BP  CALL zero(hi(11),10)        
c         
      xi(1)=xia        
      averg=av*ergev   
      DO 112 l=1,4     
  112 dxt(l)=exh*xi(l)/ah-xi(10+l)/ahe       
c         
      hi1=averg*(exh*xia*x/ah+xi(11)*y/ahe+z*ueh(1))      
      hi(1)=hi1        
      ht(1)=ht(1)+hi1  
      ii=4             
      DO 115 l=2,4     
      dhi=exh*xi(l)*x/ah+xi(l+10)*y/ahe+z*ueh(l)          
      tstl=l==4      
      IF(tstl) dhi=dhi+dxt(1)     
      dhi=averg*dhi    
      ht(l)=ht(l)+dhi  
      hi(l)=dhi        
      IF(nosd) GOTO 115          
      DO 114 m=l,4     
      ii=ii+1          
      dhi=exh*xi(ii)*x/ah+xi(10+ii)*y/ahe+z*ueh(ii)       
      IF(tstl) dhi=dhi+dxt(m)     
      IF(m==4) dhi=dhi+dxt(l)   
      dhi=averg*dhi    
      ht(ii)=ht(ii)+dhi           
  114 hi(ii)=dhi       
  115 CONTINUE         
c  pressure and enthalpy of heavy particles  
  120 anh=(x/ah+y/ahe+z/az)*av    
c  k*t, in ergs        
      tk=ck2*t         
      rhh=rho(1)       
      phh=anh*rhh*tk   
      ph(1)=phh        
      ph(2)=amm*phh*rho(2)        
      drht=1+rho(3)    
      ph(3)=amm*phh*drht          
      ph(4)=amm*phh*rho(4)+rhh*tk*avd1       
      DO 121 i=1,4     
  121 pt(i)=pt(i)+ph(i)           
      IF(nosd) GOTO 125          
      ph(10)=amm*(phh*rho(10)+rho(4)*(ph(4)+rhh*tk*avd1))            
      DO 122 k=1,3     
      k1=k+1           
      ph(k+4)=amm*(ph(k1)*rho(2)+phh*rho(k+4))            
      IF(k > 1) ph(k+6)=amm*(ph(k1)*drht+phh*rho(k+6))   
  122 CONTINUE         
      DO 123 i=5,10    
  123 pt(i)=pt(i)+ph(i)           
c         
      IF(notd) GOTO 125          
      DO 124 k=1,3     
      k1=k+1           
      DO 124 l=k,3     
      kl=nuu(k,l)      
      pt(6+kl)=pt(6+kl)+amm*(ph(kl)*rho(2)+ph(k1)*rho(l+4)+          
     .   ph(l+1)*rho(k+4)+phh*rho(kl+6))     
      IF(k > 1) pt(9+kl)=pt(9+kl)+amm*(ph(kl)*drht+ph(k1)*          
     .   rho(6+l)+ph(l+1)*rho(6+k)+phh*rho(9+kl))         
  124 CONTINUE         
      pt(20)=pt(20)+amm*(ph(10)*rho(4)+2*ph(4)*rho(10)    
     .   +phh*rho(20)+rhh*tk*avd1*(amm*rho(4)*rho(4)+rho(10)))       
c         
  125 hhh=2.5*anh*tk   
      hh(1)=hhh        
      hh(2)=0          
      hh(3)=amm*hhh    
      hh(4)=2.5*tk*avd1           
      hh(5)=0          
      hh(6)=0          
      hh(7)=0          
      hh(8)=amm2*hhh   
      hh(9)=amm*hh(4)  
      hh(10)=0         
      hh(11:20) = 0.  !! BP  CALL zero(hh(11),10)        
      hh(17)=amm*hh(8)            
      hh(18)=amm*hh(9)            
      ht(1)=ht(1)+hhh  
      ht(3)=ht(3)+amm*hhh         
      hh4=2.5*tk*avd1  
      ht(4)=ht(4)+hh4  
      IF(nosd) GOTO 130          
      ht(8)=ht(8)+amm2*hhh        
      ht(9)=ht(9)+amm*hh4         
      IF(notd) GOTO 130          
      ht(17)=ht(17)+amm3*hhh      
      ht(18)=ht(18)+amm2*hh4      
c  pressure and enthalpy of radiation        
  130 t4=t*t*t*t       
      prr=car*t4/3     
      pr(1)=prr        
      pr(2:10) = 0.   !! BP  CALL zero(pr(2),9)          
      pr(3)=4*amm*pr(1)           
      pr(8)=4*amm*pr(3)           
      pt(1)=pt(1)+prr  
      pt(3)=pt(3)+4*amm*prr       
      IF(nosd) GOTO 135          
      pt(8)=pt(8)+16*amm2*prr     
      IF(notd) GOTO 135          
      pt(17)=pt(17)+64*amm3*prr   
c         
  135 hrr=4*pr(1)/rhh  
      hr(1)=hrr        
      DO 136 i=2,4     
  136 dxt(i)=-rho(i)   
      dxt(3)=4+dxt(3)  
      ii=4             
      jj=10            
      DO 140 l=1,3     
      l1=l+1           
      hr(l1)=amm*hrr*dxt(l1)      
      IF(nosd) GOTO 140          
      DO 138 m=l,3     
      ii=ii+1          
      m1=m+1           
      hr(ii)=amm*(hr(l1)*dxt(m1)-hrr*rho(ii))             
      IF(notd) GOTO 138          
      DO 137 n=m,3     
      jj=jj+1          
      ln=nuu(l,n)      
      mn=nuu(m,n)      
  137 hr(jj)=amm*(hr(ii)*dxt(n+1)-amm*hrr*dxt(m1)*rho(ln)-hr(l1)*rho(mn)        
     .   -hrr*rho(jj))            
  138 CONTINUE         
  140 CONTINUE         
      DO 145 i=1,ider  
  145 ht(i)=ht(i)+hr(i)           
c  change to derivatives of log p 
      ptt=pt(1)        
      IF(notd) GOTO 155          
      jj=10            
      DO 152 l=1,3     
      DO 152 m=l,3     
      DO 152 n=m,3     
      lm=nuu(l,m)      
      ln=nuu(l,n)      
      mn=nuu(m,n)      
      jj=jj+1          
  152 pt(jj)=(pt(jj)+(-(pt(lm)*pt(n+1)+pt(ln)*pt(m+1)+pt(mn)*pt(l+1))           
     .   +2*pt(m+1)*pt(l+1)*pt(n+1)/ptt)/ptt)/ptt/amm     
  155 ii=4             
      DO 158 l=2,4     
      ptl=pt(l)/ptt/amm           
      IF(nosd) GOTO 158          
      DO 156 m=l,4     
      ii=ii+1          
  156 pt(ii)=(pt(ii)-pt(l)*pt(m)/ptt)/ptt/amm             
  158 pt(l)=ptl        
c         
c  cp and dad          
c         
  160 pf=pt(2)         
      hf=ht(2)         
      dxtt=ht(3)*pf-hf*pt(3)      
      lt=6             
      DO 165 l=1,3     
      lf=l+4           
      IF(l > 1) lt=6+l           
  165 dxt(l+1)=ht(lt)*pf+ht(3)*pt(lf)-ht(lf)*pt(3)-hf*pt(lt)         
c         
      cpp=dxtt/amm/t/pf           
      cp(1)=cpp        
      IF(nosd) GOTO 173          
      DO 170 l=2,4     
      dcp=dxt(l)/amm/t/pf-cpp*pt(l+3)/pf     
      IF(l==3) dcp=dcp-cpp*amm  
  170 cp(l)=dcp        
c         
  173 prh=amm*ptt/rhh  
      anum=pf*prh-hf   
      dad(1)=anum/dxtt            
      IF(nosd) GOTO 177          
      DO 175 l=2,4     
      lf=l+3           
  175 dad(l)=(prh*(amm*(pt(l)-rho(l))*pf+pt(lf))-ht(lf)-anum*dxt(l)  
     .   /dxtt)/dxtt   
c  further derivatives            
  177 rhf=rho(2)       
      dxtt=pt(3)*rhf-pf*rho(3)    
      gm1=pf/(rhf-dad(1)*dxtt)    
      tprh=rhf/dxtt    
      trhp=-pf/dxtt    
      rhxp=amm*x*(rho(4)-rhf*pt(4)/pf)       
c  delta and derivatives          
      delta=-1/trhp    
      dlt(1)=delta     
      lt=6             
      IF(nosd) GOTO 190          
      DO 180 l=2,4     
      lf=l+3           
      IF(l > 2) lt=l+5           
  180 dlt(l)=(pt(lt)*rhf+pt(3)*rho(lf)-pt(lf)*rho(3)-pf*rho(lt))/pf  
     .   -delta*pt(lf)/pf         
  190 xii1(1)=xii(1)   
      xii1(2)=xii(11)  
      xii1(3)=xii(21)  
      xii1(4)=anuh(1)/anh0        
C
	END SUBROUTINE eqstf

c**********************************************************************
   
	SUBROUTINE eqstp(pl,tl,x,y,z,fl)       
C
c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel, B.Pichon
c	CESAM2k
         
c  this routine iterates s/r eqstf in order to find the appropriate  
c  fl value. at the END, eqstf is CALLed to compute td-quantities.   
c         
	IMPLICIT REAL*8 (a-h,o-z)
   
      LOGICAL nosdit,notdit       
      COMMON/eqstd/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),         
     *  dlt(4),gm1,tprh,trhp,rhxp 
      COMMON/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,  
     *  ct,crho,cpe,che,ca03,caa,ckh,car,ergev         
      DIMENSION flini(2),wt(2)
	LOGICAL ok
	COMMON /marche/ok 

	INTEGER iter,icase,iextr,i

c         
c         
      iter = 15        
      eps = 1.d-9      
      epsf = 1.d-3     
      carsav=car/3.    
      beta=1.-carsav*10.**(4.*tl-pl)         
      beta=MAX(beta,1.0d-6)
      pgl0=pl + LOG10(beta)      
c...... initial guess of fl       
c         
c================ icase=1 for gas pressure as argument    
      icase=1          
      CALL inteff(tl,pgl0,rlbid,flini,wt,icase,iextr)     
      flin=wt(1)*flini(1)+wt(2)*flini(2) ; fl=flin          
      nosdit = .TRUE. ; notdit = .TRUE. ; car    = 0.      
c........note that for the iteration no second derivatives are needed.          
c........car is switched off temporarily to account for no radiation.           
      CALL eqstf(fl,tl,x,y,z,nosdit,notdit)  
c         
      plt=LOG10(pt(1))           
      s=pgl0-plt       
      dfl=1.d1*sign(epsf,s)       
c         
c  iterate to desired pressure    
      fl=fl+dfl        
      DO 10 i=1,iter
      CALL eqstf(fl,tl,x,y,z,nosdit,notdit)  
      pln=LOG10(pt(1))           
      IF(i==1) dlpdlf=(pln-plt)/dfl        
      IF(abs(dfl) < epsf) GOTO 5          
      dlpdlf=(pln-plt)/dfl        
      plt=pln          
    5 dfl=(pgl0-pln)/dlpdlf       
      fl=fl+dfl        
      rdif=abs(dfl/fl)           
      IF(rdif < eps) GOTO 20    
   10 CONTINUE         
      PRINT 110,pgl0,tl,rdif
	ok=.FALSE.      
 20   nosdit = .FALSE.            
      notdit = .TRUE.  
c........this CALL is here to provide the second derivatives.        
c........for this car is again switched on to account for radiation. 
      car = carsav*3.  
      CALL eqstf(fl,tl,x,y,z,nosdit,notdit)  

110   FORMAT(' >>>>>>>>> warning from eqstp: number of iterations',  
	1	' not sufficient. pgl0,tl,rdif = ',2f8.4,e12.3,
	2	' utilisation de ETAT_GONG2')    

	END SUBROUTINE eqstp
C
c*************************************************************************
C
	SUBROUTINE hviona(fl,tl,x,y,z,nosd,notd,anu,ue,anur,uer)       

c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel, B.Pichon
c	CESAM2k

c-----------------------------------------------------------------

	IMPLICIT REAL*8 (a-h,o-z)

      LOGICAL nosd,notd           
      DIMENSION anu(*),ue(*),anur(*),uer(*)  
      DIMENSION eir(10),dr(10),hr(10),gr(10),chir(27),ir(6),           
     1 izr(6),abr(6),amr(6),xi(29),phi(30),hst(30)
      INTEGER :: iab,ihvz,idmu,icount,iwrite,lj,k,izj,i,im,irab,
     1 ir,jfirst,j,ljr,izi,izrj,izr,idmx,izoj,lj0,imax,ii,l   
c      INTEGER*4 iz,name     
c      am(23),name(23),chi(374)
c      COMMON/potetc/ chi,am,iz,name    
      COMMON/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,  
     1 ct,crho,cpe,che,ca03,caa,ckh,car,ergev                     
      COMMON/hvabnd/ ab(10),iab   
      COMMON/eqscnt/ anz0,ihvz    
      COMMON/dmuder/ dmup(10),idmu    
      COMMON/ln10/ amm,amm2,amm3       

      DATA icount/0/,iwrite/3/
         
c--------------------------------------------------------------------- 
        
      IF(icount > 0) GOTO 10    
      icount=1         
      IF(iwrite /= 2) GOTO 61    
      WRITE(*,200) ; lj=0             
      DO k=1,10      
       izj=iz(k)        
       DO i=1,izj     
        lj=lj+1 ; xi(i)=chi(lj)
       ENDDO    
       im=MIN(izj,13) ; WRITE(*,210)name(k),izj,(xi(i),i=1,im)
      ENDDO                
  200 FORMAT(///,' ionization potentials:',/)  
  210 FORMAT(1x,a4,i4,13f9.2)
       
c     reset to restricted set of variables 
     
   61 irab=3 ; ir(1)=1 ; ir(2)=3 ; ir(3)=10
            
c     element with lowest ionization potential (here fe)
     
      jfirst=3 ; j=1 ; lj=0 ; ljr=0            
      DO i=1,10      
       izi=iz(i)        
       IF(i /= ir(j))GOTO 5      
       izrj=izi         
       IF(j > 2)izrj=1           
       izr(j)=izrj ; amr(j)=am(i)     
       DO k=1,izrj    
        lj=lj+1 ; ljr=ljr+1 ; chir(ljr)=chi(lj) 
       ENDDO          
       izi=izi-izrj ; j=j+1            
    5  lj=lj+izi
      ENDDO       
      abr(1)=1.02d0*ab(1) ; abr(2)=1.35d0*ab(3) ; abr(3)=3.71d0*ab(10)
                
c     reset anz0 
         
      anz0=0.d0           
      DO j=1,irab   
       anz0=anz0+abr(j)*izr(j)/amr(j)
      ENDDO
      IF(iwrite /= 1)GOTO 8
      WRITE(*,220) ; lj=0             
      DO k=1,irab    
       izrj=izr(k)      
       DO i=1,izrj    
        lj=lj+1 ; xi(i)=chir(lj)
       ENDDO   
       im=MIN(izrj,13)            
       WRITE(*,210) name(ir(k)),izrj,(xi(i),i=1,im)
      ENDDO      
  220 FORMAT(///,' ionization potentials after resetting:',//)         
      WRITE(*,230) anz0           
  230 FORMAT(//,' now anz0 =',f10.5//)
          
c     change to summed ionization potentials 
   
    8 lj=0             
      DO k=1,irab    
       sum=0.d0 ; izj=izr(k)       
       DO i=1,izj     
        lj=lj+1 ; sum=sum+chir(lj) ; chir(lj)=sum
       ENDDO    
      ENDDO
   10 f=1.d1**fl ; t=1.d1**tl
             
c     k*t, in ev 
        
      tk=ck1*t         
      IF(idmu == 1)GOTO 15
      CALL phder(fl,tl,phi,hst,nosd,notd)
         
c     psi and derivatives

      wf= SQRT(1.d0+f) ; psi=2.d0*wf+LOG(f/(1.d0+wf)**2)       
      psif=amm*wf ; psff=amm2*f/2.d0/wf
      zf=x+2.d0*y+6.d0*z ; zf2=zf*zf ; zf3=zf*zf2 ; ak0=ckh*zf2
      
c     a0**3*ne/(k*t)
      
      aa=caa*phi(1)/(tk*zf3)      
         
c     delta mu and derivatives       
         
      bmu=tk+20.d0*ak0 ; dmu=aa*bmu ; dmps=dmu-psi ; dmup(1)=dmps     
      ref=phi(2) ; ret=phi(3) ; dmuf=dmu*ref*amm            
      dmut=aa*(bmu*ret-20.d0*ak0)*amm ; dmux=aa*(3.d0*tk+20.d0*ak0)/zf    
      dmup(2)=dmuf-psif ; dmup(3)=dmut ; dmup(4)=dmux     
         
      IF(nosd)GOTO 18           
      dmup(5)=(dmuf*ref+dmu*phi(4)*amm)*amm-psff          
      dmup(6)=(dmut*ref+dmu*phi(5)*amm)*amm  
      dmup(7)=dmux*ref*amm        
      dmup(8)=aa*(20.d0*ak0*(1.d0-2.d0*ret)+bmu*(phi(6)+ret*ret))*amm2        
      dmup(9)=aa*(ret*(3.d0*tk+20.d0*ak0)-20.d0*ak0)*amm/zf        
      dmup(10)=aa*(12.d0*tk+40.d0*ak0)/zf2
      GOTO 18         
   15 dmps=dmup(1)
   18 idmx=10          
      IF(nosd)idmx=4  
      DO i=1,10     
       anu(i)=0.d0 ; ue(i)=0.d0 ; anur(i)=0.d0 ; uer(i)=0.d0
      ENDDO  
      lj=0 ; tki=1.d0/tk      
      b50: DO j=1,irab   
       izoj=iz(ir(j)) ; izj=izr(j) ; lj0=lj           
       DO i=1,izj    
        lj=lj+1; dxp=i*dmps-chir(lj)*tki ; xi(i)=dxp
       ENDDO 
              
c      complete or no ionization?     
         
       IF(izj-1)21,21,25
      
c      only one level
      
   21  phm=xi(1)        
c       IF(phm)22,22,23
       IF(phm <= 0.d0)THEN
        imax=0 ; phmh=phm ; phm=0.d0
       ELSE         
        imax=1 ; phmh=0.d0
       ENDIF
       GOTO 30
               
c      more than one level
            
   25  ii=izj+1 ; xi(ii)=0.d0 ; imax=1 ; phm=xi(1)
              
c      largest xi 
         
       b26: DO i=2,ii     
        IF(xi(i) <= phm)CYCLE b26   
        imax=i ; phm=xi(i)        
       ENDDO b26
                
c      largest but one xi
 
       phmh=-100.d0        
       DO i=1,ii     
        IF(i /= imax)phmh=MAX(phmh,xi(i))      
       ENDDO         
         
       IF(imax == ii)imax=0
      
c      test for complete ionization 
  
   30  dphm=phm-phmh ; dptst=19.d0         
       IF(imax < izj .OR. dphm <= dptst)GOTO 32
                 
c      complete ionization
           
       fct1=abr(j)/amr(j) ; anur(j)=fct1*izj ; uer(j)=fct1*chir(lj)        
       anu(1)=anu(1)+anur(j) ; ue(1)=ue(1)+uer(j) ; CYCLE b50         
               
c      test for no ionization

   32  IF(imax > 0)GOTO 34      
       IF(j == jfirst)dptst=160.d0   
       IF(dphm > dptst)CYCLE b50
        
c      general case
        
   34  DO i=1,idmx   
        dr(i)=0.d0 ; hr(i)=0.d0 ; gr(i)=0.d0 
       ENDDO      
       IF(phm <= dptst)dr(1)=omega(0,izoj)*EXP(-phm)     
       b40: DO i=1,izj    
        cchi=chir(lj0+i) ; dxp=xi(i)-phm    
        IF(dxp < -dptst)CYCLE b40 
        dxp=omega(i,izoj)*EXP(dxp) ; eir(1)=1         
        DO k=2,4      
         eir(k)=i*dmup(k)
        ENDDO            
        eir(3)=eir(3)+cchi*amm*tki ; ii=4             
        IF(nosd)GOTO 38           
        DO k=2,4      
         DO l=k,4      
          ii=ii+1 ; eir(ii)=eir(k)*eir(l)+i*dmup(ii)
         ENDDO
        ENDDO    
        eir(8)=eir(8)-amm2*cchi*tki
   38   DO k=1,idmx   
         eeir=dxp*eir(k) ; dr(k)=dr(k)+eeir ; hr(k)=hr(k)+i*eeir          
         gr(k)=gr(k)+cchi*eeir
        ENDDO
       ENDDO b40
       dr1i=1.d0/dr(1) ; fct1=abr(j)/(amr(j)*dr(1)) ; hrdr=hr(1)*dr1i
       grdr=gr(1)*dr1i ; anur(j)=fct1*hr(1) ; uer(j)=fct1*gr(1)
       anu(1)=anu(1)+anur(j) ; ue(1)=ue(1)+uer(j)
                
c      derivatives
         
       ii=4             
       b48: DO k=2,4      
        anu(k)=anu(k)+fct1*(hr(k)-dr(k)*hrdr)  
        ue(k)=ue(k)+fct1*(gr(k)-dr(k)*grdr)    
        IF(nosd)CYCLE b48           
        DO l=k,4      
         ii=ii+1          
         anu(ii)=anu(ii)+fct1*(hr(ii)-(dr(k)*hr(l)+dr(l)*hr(k)+(dr(ii)  
     1   -2*dr(k)*dr(l)*dr1i)*hr(1))*dr1i)   
         ue(ii)=ue(ii)+fct1*(gr(ii)-(dr(k)*gr(l)+dr(l)*gr(k)+(dr(ii)    
     1   -2*dr(k)*dr(l)*dr1i)*gr(1))*dr1i)
        ENDDO    
       ENDDO b48       
      ENDDO b50         

      END SUBROUTINE hviona

c************************************************************************

	SUBROUTINE inteff(tl,pgl,rl,flini,wt,icase,iextr)   

c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel, B.Pichon
c	CESAM2k

	IMPLICIT REAL*8 (a-h,o-z)
  
c>>>>> finds starting value of flini for eff routines if the argument           
c>>>>> is either LOG10(gas pressure) or LOG10(density). iextr is     
c>>>>> set to 1 if extrapolation outside the table had to be made.   
c>>>>>    
c>>>>> based on eff-pre-computation with x=.73,y=.25,z=.02           
c>>>>>    
c===== icase = 1: input=(tl,pgl),  rl ignored             
c===== icase = 2: input=(tl,rl ), pgl ignored             
c 

	INTEGER nm,nml,nmh,ntl,nth
       	PARAMETER (nm=11, nml=5, nmh=6, ntl=nml*nm, nth=nmh*nm) 
c         
      DIMENSION tLOG(nm),rLOG(nm,nm),pgLOG(nm,nm),fLOG(nm,nm)        
      DIMENSION rll(ntl),rlh(nth),pgll(ntl),pglh(nth),fll(ntl),flh(nth)         
      DIMENSION flini(2),wt(2)    

	INTEGER icase,iextr,it,ip,ir,i,it1,ip1,ir1

c         
      EQUIVALENCE(rlog (1,1),rll (1)) ,(rlog (1,nmh),rlh (1))        
      EQUIVALENCE(pgLOG(1,1),pgll(1)) ,(pgLOG(1,nmh),pglh(1))        
      EQUIVALENCE(flog (1,1),fll (1)) ,(flog (1,nmh),flh (1))        
c         
      DATA tlog/       
     .    3.40000,   3.80000,   4.20000,   4.60000,   5.00000,       
     .    5.40000,   5.80000,   6.20000,   6.60000,   7.00000,       
     .    7.40000/     
c         
      DATA rll/        
     .  -15.93806, -14.59077, -13.11393, -11.71751, -10.38835,       
     .   -8.97767,  -7.59375,  -6.16023,  -4.74661,  -3.34048,       
     .   -1.93629,     
     .  -14.39806, -13.08988, -11.78384, -10.46730,  -9.15835,       
     .   -7.85767,  -6.54375,  -5.24023,  -3.92661,  -2.62048,       
     .   -1.31627,     
     .  -12.79806, -11.57078, -10.36165,  -9.14343,  -7.92835,       
     .   -6.71767,  -5.50375,  -4.28023,  -3.06660,  -1.85045,       
     .   -0.63617,     
     .  -11.19806, -10.06214,  -8.95086,  -7.82325,  -6.69833,       
     .   -5.57766,  -4.45374,  -3.33021,  -2.20654,  -1.08028,       
     .    0.04407,     
     .   -9.59806,  -8.56014,  -7.52659,  -6.49503,  -5.46791,       
     .   -4.43760,  -3.40366,  -2.37003,  -1.33614,  -0.30961,       
     .    0.72372/     
      DATA rlh/        
     .   -7.99806,  -7.05242,  -6.11177,  -5.17967,  -4.23200,       
     .   -3.29677,  -2.35276,  -1.41861,  -0.47419,   0.46015,       
     .    1.40286,     
     .   -6.39806,  -5.55147,  -4.69828,  -3.84877,  -3.00142,       
     .   -2.15675,  -1.30499,  -0.46287,   0.38404,   1.23722,       
     .    2.08971,     
     .   -4.79807,  -4.04144,  -3.28452,  -2.53132,  -1.76802,       
     .   -1.01603,  -0.25685,   0.49360,   1.25430,   2.00696,       
     .    2.76640,     
     .   -3.19831,  -2.53173,  -1.86318,  -1.20584,  -0.53947,       
     .    0.12553,   0.78489,   1.45261,   2.11985,   2.77988,       
     .    3.44720,     
     .   -1.59748,  -1.02039,  -0.44723,   0.12189,   0.69032,       
     .    1.26763,   1.83971,   2.41156,   2.98379,   3.55092,       
     .    4.12641,     
     .    0.00202,   0.48728,   0.96517,   1.44458,   1.92401,       
     .    2.40470,   2.88056,   3.36421,   3.84159,   4.32331,       
     .    4.80062/     
c         
      DATA pgll/       
     .   -4.71634,  -2.70110,  -0.79422,   1.01898,   2.74815,       
     .    4.55883,   6.34275,   8.17627,   9.98989,  11.79601,       
     .   13.60020,     
     .   -3.17634,  -1.30952,   0.53582,   2.26909,   3.97815,       
     .    5.67883,   7.39275,   9.09627,  10.80989,  12.51602,       
     .   14.22022,     
     .   -1.57634,   0.09093,   1.95691,   3.59094,   5.20815,       
     .    6.81883,   8.43275,  10.05627,  11.66990,  13.28603,       
     .   14.90028,     
     .    0.02366,   1.56549,   3.35740,   4.90076,   6.43816,       
     .    7.95883,   9.48275,  11.00628,  12.52992,  14.05612,       
     .   15.58052,     
     .    1.62366,   3.06194,   4.76523,   6.22483,   7.66836,       
     .    9.09886,  10.53279,  11.96637,  13.40014,  14.82664,       
     .   16.26059/     
      DATA pglh/       
     .    3.22366,   4.56932,   6.08520,   7.53751,   8.90120,       
     .   10.23927,  11.58325,  12.91715,  14.26164,  15.59761,       
     .   16.94023,     
     .    4.82375,   6.07052,   7.37716,   8.83489,  10.11675,       
     .   11.37463,  12.62823,  13.87314,  15.12352,  16.37610,       
     .   17.63041,     
     .    6.42746,   7.58888,   8.76598,  10.05827,  11.32914,       
     .   12.52233,  13.69271,  14.83754,  16.00047,  17.15703,       
     .   18.32207,     
     .    8.15380,   9.30678,  10.50206,  11.73099,  12.68725,       
     .   13.70433,  14.76934,  15.84437,  16.91982,  17.98796,       
     .   19.06506,     
     .   10.79939,  11.94400,  13.08449,  13.22093,  14.17459,       
     .   15.14145,  16.10082,  17.05997,  18.01959,  18.97151,       
     .   19.93379,     
     .   13.96914,  14.93979,  14.51362,  15.31529,  16.11750,       
     .   16.92229,  17.71977,  18.53028,  19.33067,  20.13707,       
     .   20.93387/     
c         
      DATA fll/        
     .  -15.20875, -12.62213, -11.68335, -10.85416, -10.12500,       
     .   -9.31434,  -8.53047,  -7.69708,  -6.88379,  -6.07850,       
     .   -5.27639,     
     .  -13.66875, -11.40213, -10.35335,  -9.60416,  -8.89500,       
     .   -8.19434,  -7.48047,  -6.77708,  -6.06379,  -5.35850,       
     .   -4.65639,     
     .  -12.06875, -10.50213,  -8.93335,  -8.28416,  -7.66500,       
     .   -7.05434,  -6.44047,  -5.81708,  -5.20379,  -4.58850,       
     .   -3.97639,     
     .  -10.46875,  -9.61213,  -7.54335,  -6.98416,  -6.43500,       
     .   -5.91434,  -5.39047,  -4.86708,  -4.34379,  -3.81850,       
     .   -3.29639,     
     .   -8.86875,  -8.40213,  -6.15335,  -5.66416,  -5.20500,       
     .   -4.77434,  -4.34047,  -3.90708,  -3.47379,  -3.04850,       
     .   -2.61639/     
      DATA flh/        
     .   -7.26875,  -6.92213,  -4.97335,  -4.35416,  -3.97500,       
     .   -3.63434,  -3.29047,  -2.95708,  -2.61379,  -2.27850,       
     .   -1.93639,     
     .   -5.66875,  -5.42213,  -4.11335,  -3.09416,  -2.77500,       
     .   -2.50434,  -2.25047,  -2.00708,  -1.75379,  -1.49850,       
     .   -1.24639,     
     .   -4.06875,  -3.91213,  -3.30335,  -2.10416,  -1.64500,       
     .   -1.40434,  -1.20047,  -1.03708,  -0.87379,  -0.71850,       
     .   -0.55639,     
     .   -2.46875,  -2.40213,  -2.26335,  -1.25416,  -0.26500,       
     .   -0.17434,  -0.11047,  -0.03708,   0.03621,   0.10150,       
     .    0.17361,     
     .   -0.85875,  -0.88213,  -0.90335,   1.24584,   1.20500,       
     .    1.17566,   1.13953,   1.10292,   1.06621,   1.02150,       
     .    0.98361,     
     .    0.92125,   0.77787,   3.15665,   2.99584,   2.83500,       
     .    2.67566,   2.50953,   2.35292,   2.18621,   2.02150,       
     .    1.84361/     
c         
      it=-999          
      ip=-999          
      ir=-999          
      iextr=0          
c========== select isotherms for linear inter(extra)polation         
      DO 11 i=1,nm     
      IF(tl < tLOG(i)) GOTO 20   
      it=i             
 11   CONTINUE         
         iextr=1       
         it=nm-1       
c         
  20  IF(it==-999) THEN         
         iextr=1       
         it=1          
      ENDIF           
c========== tl-part of arguments  
      it1=it+1         
      x0 =tLOG(it)     
      x1 =tLOG(it1)    
      x2 =x0           
      x3 =x1           
      x  =  tl         
c         
c................................. icase = 1 ............................       
      IF(icase==2) GOTO 200     
c         
c========== select pressure points for linear inter(extra)polation   
      DO 21 i=1,nm     
      IF(pgl < pgLOG(it,i)) GOTO 30         
      ip=i             
  21  CONTINUE         
         iextr=1       
         ip=nm-1       
c         
  30  IF(ip==-999) THEN         
         iextr=1       
         ip=1          
      ENDIF           
c         
c========== define the three fonction values for inter(extra)polation           
      ip1=ip+1         
      y0 =pgLOG(it,ip)            
      y1 =pgLOG(it1,ip)           
      y2 =pgLOG(it,ip1)           
      y3 =pgLOG(it1,ip1)          
      z0 =fLOG(it,ip)  
      z1 =fLOG(it1,ip)            
      z2 =fLOG(it,ip1)            
      z3 =fLOG(it1,ip1)           
c         
c========== define arguments      
      y  = pgl         
      GOTO 1000        
c................................. icase = 2 ............................       
c         
c========== select density points for linear inter(extra)polation    
 200  DO 31 i=1,nm     
      IF(rl < rLOG(it,i)) GOTO 40           
      ir=i             
  31  CONTINUE         
         iextr=1       
         ir=nm-1       
c         
  40  IF(ir==-999) THEN         
         iextr=1       
         ir=1          
      ENDIF           
c         
c========== define the three fonction values for inter(extra)polation           
      ir1=ir+1         
      y0 =rLOG(it,ir)  
      y1 =rLOG(it1,ir)            
      y2 =rLOG(it,ir1)            
      y3 =rLOG(it1,ir1)           
      z0 =fLOG(it,ir)  
      z1 =fLOG(it1,ir)            
      z2 =fLOG(it,ir1)            
      z3 =fLOG(it1,ir1)           
c         
c========== define argument       
      y  =  rl         
c         
c========== CALL bilinear interpolation      
c         
1000  CONTINUE         
c         
c....... lower triangle           
      CALL  bilin (x0,y0,z0,x1,y1,z1,x2,y2,z2,x,y,z)      
      flini(1) = z     
c         
c....... upper triangle           
      CALL  bilin (x3,y3,z3,x2,y2,z2,x1,y1,z1,x,y,z)      
      flini(2) = z     
c         
c....... weights (quite arbitrary)           
      wlow = 1./( (x-x0)**2 + (y-y0)**2 + 1.e-5 )         
      whig = 1./( (x-x3)**2 + (y-y3)**2 + 1.e-5 )         
      wtot = wlow + whig          
      wlow = wlow/wtot            
      whig = whig/wtot            
      wt(1) = wlow     
      wt(2) = whig     
c         
      RETURN

	END SUBROUTINE inteff

c*************************************************************************

	function omega(i,iz)   

c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel, B.Pichon
c	CESAM2k

      IMPLICIT REAL*8 (a-h,o-z)

c  calculates statistical weight of i-th ionization stage of element 
c  with number iz      
      DIMENSION iom(26),iom1(20)

	REAL*8 omega
	INTEGER i,iz,iom,iom1

      DATA iom/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,10,21,28,25,6,25, 
     .   28,21/        
      DATA iom1/2,1,1,2,10,15,21,28,28,25,7,6,6,7,25,30,28,21,21,10/ 
c         
      IF(i <= 1.and.iz.ge.19) GOTO 20       
      IF(i==iz) GOTO 15        
      omega=iom(iz-i)  
      RETURN           
   15 omega=15         
      RETURN           
   20 omega=iom1(2*(iz-19)+i+1)   
      RETURN           

	END FUNCTION omega
C
c*************************************************************************
C
	SUBROUTINE phder(fl,tl,phi,hst,nosd,notd)  

c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel, B.Pichon
c	CESAM2k
        
	IMPLICIT REAL*8 (a-h,o-z)

      LOGICAL nosd,notd           
      DIMENSION phi(30),hst(30)           
      DIMENSION sig(10),tmn(10),ff(4),gg(4),c(48)           
      COMMON/phdsms/ s0,sf,sg,sff,sfg,sgg,sfff,sffg,sfgg,sggg,       
     .   cfg,tf,tg,tff,tfg,tgg,tfff,tffg,tfgg,tggg        
      COMMON/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,  
     .   ct,crho,cpe,che,ca03,caa,ckh,car,ergev    
      COMMON/ln10/ amm,amm2,amm3
	INTEGER ic,imax,i,icase,ic2,l,mio,kk,l0,k,in,im,ig,ik

c  LOGICAL EQUIVALENCE            
      EQUIVALENCE(s0,sig(1)),(cfg,tmn(1))    
c     ***********************************    
c  coefficients        
      DATA ic,c/4,2.315472,7.128660,7.504998,2.665350,7.837752,23.507934        
     .,23.311317,7.987465,9.215560,26.834068,25.082745,8.020509,3.693280        
     .,10.333176,9.168960,2.668248,2.315472,6.748104,6.564912,2.132280,         
     .7.837752,21.439740,19.080088,5.478100,9.215560,23.551504,19.015888        
     .,4.679944,3.693280,8.859868,6.500712,1.334124,1.157736,3.770676,          
     .4.015224,1.402284,8.283420,26.184486,28.211372,10.310306,14.755480        
     .,45.031658,46.909420,16.633242,7.386560,22.159680,22.438048,   
     .7.664928/        
c         
c  number of sums      
      imax=10          
      IF(notd) imax=6  
      IF(nosd) imax=3  
c  f,g    
      f=1.d1**fl       
      t=1.d1**tl       
      ts=ct*t          
      wf=SQRT(1.d0+f)    
      g1=ts*wf          
c         
c  1/(1+f), 1/(1+g1) etc.          
c         
      vf=1.d0/(1.d0+f)       
      vg=1.d0+g1
      fdf=vg*g1         
      fdf=vf*fdf*SQRT(fdf)       
      vg=1.d0/vg          
      vfg=vf*vg        
      ug=g1*vg          
      uf=f*vf          
c         
      ug2=ug*vg        
      ug3=vg*ug2       
c  powers of f and g   
      ff(1)=f          
      gg(1)=1.d0          
      DO 10 i=2,ic     
      ff(i)=f*ff(i-1)  
      gg(i)=g1*gg(i-1)  
   10 fdf=fdf*vfg      
c  test on size of f and g        
      icase=1          
      IF(.not.notd)GOTO12      
      IF(f < 1.d-4) icase=icase+1           
      IF(g1 < 1.d-4) icase=icase+2           
c         
   12 ic2=ic*ic        
c         
c  calculate phi* and derivatives 
c         
      l=1   
      anu=1.5d0        
      mio=ic           
      an32=2.5d0-ic    
      kk=-10           
c         
      l0=1             
      DO 50 k=1,3      
      kk=kk+10         
      IF(k-2) 18,16,17            
c  reset fdf for k=2 and 3        
   16 anu=2.5d0        
      fdf=g1*fdf        
      GOTO18         
   17 mio=mio+1        
      fdf=fdf*vf       
   18 DO 19 i=1,imax   
   19 sig(i)=0.d0      
      annu=anu-1       
c         
c  the summation       
c         
      l=l0             
      GOTO(20,25,30,35), icase  
c  the general case    
   20 DO 23 in=1,ic    
      annu=annu+1      
      DO 23 im=1,ic    
c  phimn*(f**(m+1))*(g1**n)        
      cfg=c(l)*ff(im)*gg(in)      
c         
      tg=annu*cfg      
      tf=im*cfg        
      IF(nosd) GOTO 21           
c  second derivatives  
      tgg=annu*tg      
      tfg=im*tg        
      tff=im*tf        
      IF(notd) GOTO 21           
c  third derivatives   
      tggg=annu*tgg    
      tfgg=im*tgg      
      tffg=im*tfg      
      tfff=im*tff      
c  summing             
   21 DO 22 i=1,imax   
   22 sig(i)=sig(i)+tmn(i)        
   23 l=l+1            
c  the summation is finished      
      IF(nosd) GOTO 45           
c         
c  the sigma tilde (cf (22.2)) are stored in the corresponding sigma.           
c  this is o.k. provided that we go backwards.            
c         
      s02=s0*s0        
      sg2=sg*sg        
      sf2=sf*sf        
      IF(notd) GOTO 24           
c  third derivatives   
      s03=s02*s0       
      sfff=(sfff*s02-sf*(3*sff*s0-2*sf2))/s03             
      sffg=(sffg*s02-sff*sg*s0-2*sf*(sfg*s0-sg*sf))/s03   
      sfgg=(sfgg*s02-sgg*sf*s0-2*sg*(sfg*s0-sg*sf))/s03   
      sggg=(sggg*s02-sg*(3*sgg*s0-2*sg2))/s03             
c  second derivatives  
   24 sff=(sff*s0-sf2)/s02        
      sfg=(sfg*s0-sf*sg)/s02      
      sgg=(sgg*s0-sg2)/s02        
      GOTO 45         
c  f is small          
   25 DO 28 in=1,ic    
      annu=annu+1      
      DO 27 im=1,2     
      cfg=c(l)*ff(im)*gg(in)      
      sig(im)=sig(im)+cfg         
      tg=annu*cfg      
      sg=sg+tg         
      IF(nosd) GOTO 27           
      sgg=sgg+annu*tg  
      sfg=sfg+im*tg    
   27 l=l+1            
   28 l=l+2            
c  the summation is finished. set precursors for sigma tilde         
      sff=s0*sf        
      s0=s0+sf         
      sf=s0+sf         
      IF(nosd) GOTO 45           
      sgg=sgg*s0-sg*sg            
      GOTO 40         
c  g1 is small          
   30 ig=1             
      DO 33 in=1,2     
      annu=annu+1      
      DO 32 im=1,4     
      cfg=c(l)*ff(im)*gg(in)      
      sig(ig)=sig(ig)+cfg         
      tf=im*cfg        
      sf=sf+tf         
      IF(nosd) GOTO 32           
      sff=sff+im*tf    
      sfg=sfg+annu*tf  
   32 l=l+1            
   33 ig=3             
c  the summation is finished. set precursors for sigma tilde.        
      sgg=s0*sg        
      s0=s0+sg         
      sg=anu*s0+sg     
      IF(nosd) GOTO 45           
      sff=sff*s0-sf*sf            
      GOTO 40         
c  both f and g1 are small         
   35 ig=3             
c  in this case we must also zero sfg        
      sfg=0.e0         
      DO 38 in=1,2     
      annu=annu+1      
      DO 37 im=1,2     
      cfg=c(l)*ff(im)*gg(in)      
      sig(im)=sig(im)+cfg         
      sig(ig)=sig(ig)+cfg         
   37 l=l+1            
      ig=5             
   38 l=l+2            
c  the summation is finished. set precursors for sigma tilde.        
      sff=s0*sf        
      s0=s0+sf         
      sf=s0+sf         
      sgg=sg*sfg       
      sg=anu*s0+sfg    
      IF(nosd) GOTO 45           
c  set final values of the sigma tilde.      
   40 s02=s0*s0        
      sff=sff/s02      
      sgg=sgg/s02      
      IF(f*g1 < 1.00001d-8) GOTO 42         
      sfg=(sfg*s0-sf*sg)/s02      
      GOTO 45         
c  approximate expression for sfg (may need fixing up, if f = o(1)   
c  or g1 = o(1))        
   42 sfg=f*g1*(c(l0+5)-c(l0+1)*c(l0+4)/c(l0))/c(l0)       
c         
c  phi* and first derivatives     
c         
   45 phi(kk+1)=fdf*s0            
      pht=an32*ug+sg/s0           
      phi(kk+3)=pht    
      phi(kk+2)=(pht/2-mio)*uf+sf/s0         
      IF(nosd) GOTO 50           
c         
c  second derivatives of phi*.    
c         
      phtt=an32*ug2+sgg           
      phi(kk+6)=phtt   
      phi(kk+5)=phtt*uf/2+sfg     
      phi(kk+4)=sff+uf*(sfg+vf*(pht/2-mio+f*phtt/4))      
c         
      IF(notd) GOTO 50           
c  third derivatives   
      phttt=an32*ug3*(1.d0-g1)+sggg   
      phi(kk+10)=phttt            
      phi(kk+9)=sfgg+uf*phttt/2   
      phfft=sffg+uf*(sfgg+vf*(phtt+f*phttt/2)/2)
      phi(kk+8)=phfft  
      phi(kk+7)=sfff+uf*(sffg+phfft/2+vf*(1.5*sfg+f*sfgg/4           
     .   +vf*((1-f)*(pht/2-mio)+f*phtt/2)))  
   50 l0=l0+ic2        
c         
c  h* and its derivatives (pp 23-25)         
c         
      DO 55 i=2,imax   
      ik=20+i          
   55 phi(ik)=phi(ik)-phi(i)      
c         
      hs=phi(21)/phi(1)           
      wft1=2.d0*g1         
      hst(1)=hs+wft1   
c         
      hf=phi(22)       
      ht1=phi(23)       
      wft2=ts*f/wf
      hst(2)=hs*hf+wft2           
      hst(3)=hs*ht1+wft1   
      IF(nosd) GOTO 58           
c  second derivatives  
      hff=phi(24)      
      hft=phi(25)      
      htt=phi(26)      
      wft3=uf*(1+f/2)*ts/wf       
      hst(4)=hs*(hf*hf+hff)+wft3  
      hst(5)=hs*(hf*ht1+hft)+wft2  
      hst(6)=hs*(ht1*ht1+htt)+wft1  
      IF(notd) GOTO 58           
c  third derivatives   
      hst(7)=hs*(hf*(hf*hf+3*hff)+phi(27))+uf*vf*(1+f*(2+f)/4)*ts/wf 
      hst(8)=hs*(hf*(hf*ht1+2*hft)+ht1*hff+phi(28))+wft3    
      hst(9)=hs*(ht1*(ht1*hf+2*hft)+hf*htt+phi(29))+wft2    
      hst(10)=hs*(ht1*(ht1*ht1+3*htt)+phi(30))+wft1          
c  change to derivatives wrt LOG10 f and LOG10 t          
   58 fct=amm          
      DO 60 i=2,imax   
      IF(i==4 .or. i==7) fct=fct*amm       
   60 hst(i)=hst(i)*fct
             
	END SUBROUTINE phder
   
c***************************************************************************

	SUBROUTINE setcns

c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel, B.Pichon
c	CESAM2k

c	Les valeurs des constantes ont ete modifiees pour rendre
c	leur valeurs coherentes avec le ss-pg ctes du code de Pierre.
c	Michel Juillet 90.
c	les constantes modifies sont signalees par !*
c  sets the physical and mathematical constants used in the program
          
	IMPLICIT REAL*8 (a-h,o-z)
	INTEGER :: idmu

	COMMON/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,  
	1 ct,crho,cpe,che,ca03,caa,ckh,car,ergev           
	COMMON/ln10/ amm,amm2,amm3  
	COMMON/eqstd/ ccc1(90)      
	COMMON/eqsout/ ccc2(210)    
	COMMON/dmuder/ ccc3(10),idmu
                 
c  avogadro's number   
      av=6.0221366d23		!*
c  atomic weights of h and he     
      ah=1.007825d0		!*
      ahe=4.002603d0		!*
      az=17.8d0      
      avda=av*(1/ah-2/ahe)
      avd1=av*(1/ah-1/ahe)
c  boltzmann's constant in ev/deg 
      ck1=8.6170837d-5
c  the same, in ergs/deg
      ck2=1.380658d-16	!*
c  ionization potentials
      exh=13.595d0
      exhe=24.580d0
      exhep=54.403d0
c  constants for transition to starred variables          
      ct=1.686304d-10
      crho=1.759547d30
      cpe =1.440588d24
      che =8.187265d-7
c  constants for pressure ionization
      ca03=2.147d-24
      caa=1.759547d30*ca03
      ckh=13.5d0
c  change from ev to ergs         
      ca03=1.60217733d-12*ca03	!*
c  the radiation constant         
      car=7.5659122d-15		!*
c  number of ergs in 1 ev         
      ergev=1.60217733d-12		!*
c  ln 10  
      amm=LOG(1.d1)   
      amm2=amm*amm     
      amm3=amm2*amm    
c  set COMMONs from s/r eqstf to zero        
      ccc1 = 0.   !! BP  CALL zero(ccc1,90)          
      ccc2 = 0.   !! BP  CALL zero(ccc2,210)         
      ccc3 = 0.   !! BP  CALL zero(ccc3,10)          

	END SUBROUTINE setcns
 
c************************************************************************

	SUBROUTINE store(a,b,n)     
 
c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel, B.Pichon
c	CESAM2k

	IMPLICIT REAL*8 (a-h,o-z)
      
c     stores first n elements of a into b    
c        
      DIMENSION a(1),b(1)   
	INTEGER i,n    
         
   10 DO 11 i=1,n      
   11 b(i)=a(i)        

	END SUBROUTINE store

c****************************************************************************

	SUBROUTINE theffp(p10,t10,chem,ychem,f)   

c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel, B.Pichon
c	CESAM2k
       
	IMPLICIT REAL*8 (a-h,o-z)

	INTEGER ihvz       
      COMMON/eqscnt/ anh0,ihvz    
      COMMON/eqstd/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp 
      COMMON/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,  
     &  ct,crho,cpe,che,ca03,caa,ckh,car,ergev     
      COMMON/eqsaux/psi           
         
c================= output =========================================  
      COMMON/eos/ald8,cp8,psi8,rhp8,rht8,vlro8,vmol8,     
     *vmy8,vmyp8,vmyt8,vna8,      
     *xio8(3),xmol8,e8,uint8,dpdrs8,scgs8,cv8             
c==================================================================  
         
      umod  = LOG(10.d0)
      xc=chem          
      yc=ychem         
      zc=1.-xc-yc      
         
      CALL eqstp(p10,t10,xc,yc,zc,flnew)     
         
      psi8   =psi      
      vlro8  = LOG10(rho(1))     
      rhp8   =rho(2)/pt(2)        
      rht8   = -dlt(1)            
      vna8   =dad(1)   
      xio8(1)=xii1(1)  
      xio8(2)=xii1(2)  
      xio8(3)=xii1(3)  
      xmol8  = 0.      
c>>>>>>>>>>>>>>>>>> no molecules in eff      
      rmu    = xc/ah+yc/ahe+zc/az 
      vmol8  = 1./rmu  
      cp88   =cp(1)    
      rgas0  = 8.31434 e 7        
      cp8    =cp88*vmol8/rgas0    
      emax   =(xc/ah+2.*yc/ahe+zc*anh0)*av   
      e8     =ane(1)/emax         
      vmy8   =vmol8/(1.+e8)       
      fact   =1.0 + e8            
      xelp   =(ane(2)/pt(2))/emax 
      xelt   =(-ane(2)*pt(3)/pt(2) + ane(3))/emax         
      vmyp8  = - xelp/(fact*umod) 
      vmyt8  = - xelt/(fact*umod) 
      dpdrs8 = gm1
      f=10.d0**flnew
      
      RETURN

      END SUBROUTINE theffp

	END SUBROUTINE etat_eff
