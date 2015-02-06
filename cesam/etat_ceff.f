
c***************************************************************************

	SUBROUTINE etat_ceff(p,t,xchim,
	1 ro,drop,drot,drox,u,dup,dut,dux,
	2 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	
c	routine public du module mod_etat

c       adaptation du logiciel ceff écrit par Jorgen Christensen_Dalsgaard
c       à partir de P, T & X calcul  de ro & u et de toutes les dérivées
c       jusqu'au deuxième ordre
c       annie baglin le 1 09 1991
c
c	MODIF:
c	ajout d'un appel a EFF puis ETAT_GONG2 en cas de Pb.

c	P. Morel octobre 91
c	appel systématique a EFF si T <3000 ou si X <.1

c       réinitialisation de icoulm lorsqu'on doit appeler EFF
c       appel a EFF seulement lorsque X < 1.D-50
c       Y. Lebreton 1996 (pour cesam2)
c	21 04 00 : adaptation a CESAM4.2.0.0 (P.Morel)
c	suppression des nh1, nh2, nhe2, lamb des appels a EFF (P.Morel)

c	Adaptation F95 B.Pichon & P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c---------------------------------------------------------

	USE mod_donnees, only : Lz0 => z0, aradia, ihe4, nchim
	USE mod_kind

	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim	
	REAL (kind=dp), INTENT(in) :: p, t
	REAL (kind=dp), INTENT(out) :: ro,drop,drot,drox,u,dup,dut,dux,
	1 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	2 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1

	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: Lxchim
	
	REAL (kind=dp), SAVE :: aradias3

	REAL (kind=dp) :: drott, drotp, drotx, Lt, Lp, cpp,
	1 dutt, dutp, dutx,
	2 pl, tl, fl, f, drofn, dropl, dropn, drotn, droxn,
	3 dpffn,droffn,drottn,droftn,drotxn,dptxn,drofxn,
	4 dpfxn,dpftn,dfpl,dftl,dfxl,dfpn,dpttn,dftn,dpxn,dfxn,
	5 dptn,drotpn,betbet,betgam,unsro,ro2,psro2,
	6 dhfn,dhpn,dupn,dhtn,dutn,dhxn,duxn,dhttn,dhtxn,dhffn,
	7 dhftn,dhfxn,dhtpn,psro3,dutpn,duttn,dutxn

	INTEGER :: nit
			
	LOGICAL :: nosd, notd, pass=.true.,init=.true.
	
	REAL (kind=dp) :: amm,amm2,amm3
	COMMON/ln10c/amm,amm2,amm3

	REAL (kind=dp) :: av,ah_ceff,ahe,az,avda,avd1,ck1,ck2,exh,
	1 exhe,exhep,ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm    
	COMMON/constsc/av,ah_ceff,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
	1 ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm

	REAL (kind=dp), DIMENSION(20) :: ht, pt, rho
	REAL (kind=dp), DIMENSION(10) :: ane	
	REAL (kind=dp), DIMENSION(4) ::	dad, dlt, cp_ceff, gmm1, xii1
	REAL (kind=dp) :: gm1, rhxp, tprh, trhp 	
	COMMON/eqstdc/xii1,ane,rho,ht,pt,cp_ceff,dad,dlt,gm1,tprh,trhp,
	1 rhxp,gmm1

	REAL (kind=dp), DIMENSION(30) :: east, xii
	REAL (kind=dp), DIMENSION(20) :: dph, he, hi, hh, hr
	REAL (kind=dp), DIMENSION(10) :: dne, hcoul, pcoul, pe, pi_ceff,
	1 ph, pr	
	COMMON/eqsoutc/east,xii,dne,dph,he,pe,hi,pi_ceff,hh,ph,hr,pr,
	1 pcoul,hcoul
    
	INTEGER :: idgbcs,idgrhs,idgeos,idgopc,idgeng
	COMMON/diagns/idgbcs,idgrhs,idgeos,idgopc,idgeng

	INTEGER :: istdou,istdpr
	COMMON/cstdio/ istdou, istdpr

	LOGICAL :: ok
	COMMON /marche/ok
	
        REAL (kind=dp) :: anh0, anhe0
	INTEGER :: ihvz, iprrad, ihmin	
	COMMON/eqscntc/ anh0,anhe0,ihvz,iprrad,ihmin
	
	REAL (kind=dp), DIMENSION(10) :: ab
        INTEGER :: iab	      
        COMMON/hvabndc/ ab,iab

        REAL (kind=dp) :: frhi, bdcoh, bdcoz
        INTEGER :: idpco	     
	COMMON/eqdpco/ frhi,bdcoh,bdcoz,idpco
	
	REAL (kind=dp), DIMENSION(48) :: c
        INTEGER :: ic	      
	COMMON/eqphcs/ c,ic
      
	REAL (kind=dp), DIMENSION(125) :: chi
	REAL (kind=dp), DIMENSION(10) :: am
	INTEGER, DIMENSION(10) :: iz
	COMMON/potetcc/ chi, am, iz
      
	CHARACTER (len=4), DIMENSION(10) :: name
	COMMON/hvname/ name
      
	REAL (kind=dp) :: dptst0,dptst1    
        INTEGER :: icnthv,iwrthv       
	COMMON /hvcntl/ icnthv,iwrthv,dptst0,dptst1
      
	INTEGER, DIMENSION(26) :: iom 
	INTEGER, DIMENSION(20) :: iom1	     
	COMMON/hvomeg/ iom,iom1
	
        INTEGER :: iomfll      
	COMMON/hvomcl/ iomfll
	
	REAL (kind=dp) :: epsdmu, epssdr
        INTEGER ::  icoulm, iclmit, iclmsd, nitdmu	        
	COMMON/ccoulm/ epsdmu, icoulm, iclmit, iclmsd, nitdmu, epssdr
	
c---------------------------------------------------------
	
2000	FORMAT(8es10.3)

	IF(init)THEN
	
c	 initialisation des variables du BLOCKDATA blstio

	 istdou=6 ; istdpr=6
	 
c	 initialisation des variables du BLOCKDATA bleqstc	 

c	 COMMON /eqscntc/

	 anh0=0.5d0 ; anhe0=6.d0 ; ihvz=1 ; iprrad=1 ; ihmin=0
	
c	 COMMON	/hvabndc/

	 ab=(/ 0.2254d0, 0.0549d0, 0.4987d0, 0.03350d0, 0.00197d0,
	1      0.0436d0, 0.00403d0,0.0565d0, 0.00180d0, 0.0795d0 /) 	 
	 iab=10

c	 COMMON/eqdpco/

	 frhi=1.d0 ; bdcoh=1.d0 ; bdcoz=1.d0 ; idpco=0
	 
c	 COMMON/eqphcs/

	 c=(/
	1 2.315472d0, 7.128660d0, 7.504998d0, 2.665350d0, 7.837752d0,
	2 23.507934d0,23.311317d0,7.987465d0, 9.215560d0, 26.834068d0,
	3 25.082745d0,8.020509d0, 3.693280d0, 10.333176d0,9.168960d0,
	4 2.668248d0, 2.315472d0, 6.748104d0, 6.564912d0, 2.132280d0,
	5 7.837752d0, 21.439740d0,19.080088d0,5.478100d0, 9.215560d0,
	6 23.551504d0,19.015888d0,4.679944d0, 3.693280d0, 8.859868d0,
	7 6.500712d0, 1.334124d0, 1.157736d0, 3.770676d0, 4.015224d0,
	8 1.402284d0, 8.283420d0, 26.184486d0,28.211372d0,10.310306d0,
	9 14.755480d0,45.031658d0,46.909420d0,16.633242d0,7.386560d0,
	1 22.159680d0,22.438048d0,7.664928d0 /)
	 ic=4		 
	 
c	 COMMON /potetcc/
	 
	 chi=(/
	1	11.26d0, 	24.38d0,	47.86d0,	64.48d0,
	2 	391.99d0,	489.84d0,	14.54d0,	29.60d0,
	3	47.43d0,	77.45d0,	97.86d0,	551.92d0,
	4	666.83d0,	13.61d0,	35.15d0,	54.93d0,
	5	77.39d0,	113.87d0,	138.08d0,	739.11d0,
	6	871.12d0,	21.56d0,	41.07d0,	63.5d0,
	7	97.16d0,	126.4d0,	157.91d0,	207.3d0,
	8	239.d0,		1196.d0,	1360.d0,	5.14d0,
	9	47.29d0,	71.65d0,	98.88d0,	138.60d0,
	1	172.36d0,	208.44d0,	264.15d0,	299.78d0,
	2	1465.d0,	1646.d0,	7.64d0,		15.03d0,
	3	80.12d0,	109.29d0,	141.23d0,	186.86d0,
	4	225.31d0,	265.96d0,	327.90d0,	367.36d0,
	5	1761.2d0,	2085.d0,	5.98d0,		18.82d0,
	6	28.44d0,	119.96d0,	153.77d0,	190.42d0,
	7	241.93d0,	285.13d0,	330.1d0,	398.5d0,
	8	441.9d0,	2085.5d0,	2299.d0,	8.15d0,
	9	16.34d0,	33.46d0,	45.13d0,	166.73d0,
	1	205.11d0,	246.41d0,	303.87d0,	351.83d0,
	2	401.3d0,	476.0d0,	523.2d0,	2436.d0,
	3	2666.d0,	15.75d0,	27.62d0,	40.90d0,
	4	59.79d0,	75.0d0,		91.3d0,		124.d0,
	5	143.46d0,	421.d0,		480.d0,		539.5d0,
	6	21.1d0,		688.5d0,	755.5d0,	854.4d0,
	7	918.d0,		4121.d0,	4426.d0,	7.90d0,
	8	16.18d0,	30.64d0,	56.d0,		79.d0,
	9	105.d0,		133.d0,		151.d0,		235.d0,
	1	262.d0,		290.d0,		321.d0,		355.d0,
	2	390.d0,		457.d0,		489.d0,		1266.d0,
	3	1358.d0,	1456.d0,	1582.d0,	1689.d0,
	4	1799.d0,	1950.d0,	2045.d0,	8828.d0,
	5	9278.d0 /)	 
	 am=(/ 12.d0,   14.01d0, 16.d0,   20.17d0, 22.99d0,
	1       24.31d0, 26.98d0, 28.08d0, 39.94d0, 55.84d0 /)
	 iz=(/ 6, 7, 8, 10, 11, 12, 13, 14, 18, 26 /)	 
	 
c	 COMMON /hvname/

	 name=(/ '   C','   N','   O','  Ne','  Na','  Mg','  Al',
	1 '  Si','  Ar','  Fe' /)
	 
c	 COMMON /hvcntl/ icnthv,iwrthv,dptst0,dptst1	 
	 
	 icnthv=0 ; iwrthv=0 ; dptst0=85.d0 ; dptst1=19.d0
	 
c	 COMMON/hvomeg/

	 iom=(/ 2, 1, 2, 1, 6, 9, 4, 9, 6, 1, 2, 1, 6, 9, 4, 9,
	1 6, 1, 10, 21, 28, 25, 6, 25, 28, 21 /)
	 iom1=(/ 2, 1, 1, 2, 10, 15, 21, 28, 28, 25, 7, 6, 6, 7, 25,
	2 30, 28, 21, 21, 10 /)

c	 COMMON/hvomcl/

	 iomfll=1
	
c	 COMMON/ccoulm/ epsdmu, icoulm, iclmit, iclmsd, nitdmu, epssdr

	 epsdmu=1.d-12 ; icoulm=0 ; iclmit=1 ; iclmsd=1 ; epssdr=1.d-3

c~~~~~~~~~~~~~~~~~~~~~~~~
	 
	 CALL setcnsc		!initialisation des ctes
	 init=.FALSE.
	 WRITE(2,1) ; WRITE(*,1)
1	 FORMAT(//,5x,'equation d''etat CEFF',//)
	 aradias3=aradia/3.d0
	 ALLOCATE(Lxchim(nchim))
	ENDIF

	Lxchim(1:nchim)=xchim(1:nchim) ; Lp=p ; Lt=t

c	IF(t < 3000.d0 .OR. xchim(1) < 0.2d0 .OR. p < 1.d2)THEN
	IF(t < 3000.d0 .OR. xchim(1) < 1.d-50 .OR. p < 1.d2) THEN   !yveline 040100
	 IF(pass)WRITE(*,20)t,xchim(1),p
	 pass=.FALSE.
20	 FORMAT('divers appels a EFF car 3000 > t=',es10.3,
	1 ' ou Xlim > X(1)=',es10.3, ' ou 1.d2 > P=',es10.3)
	 CALL etat_eff(Lp,Lt,Lxchim,
	1	ro,drop,drot,drox,u,dup,dut,dux,
	2	delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3	gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	 RETURN
	ENDIF

c	index divers

	idgeos=3		!lie a l'appel de setf4?
	ihmin=0			!0 sans  h-, 1 avec h-
c	ihmin=1	
	iomfll=1		!utilise dans le calcul de omegac
	amm2=amm*amm
c
	ok=.true.
        pl=log10(p)
        tl=log10(t)
	nosd=.FALSE.
	notd=.true.

	IF(ihe4 > 1)THEN
c         PRINT*,'on passe ici',nosd,notd
c	  WRITE(*,2000)pl,tl,xchim(1),xchim(ihe4)+xchim(ihe4-1),
c	1 1.d0-xchim(1)-xchim(ihe4)-xchim(ihe4-1)
          CALL eqstpc(pl,tl,Lxchim(1),Lxchim(ihe4)+Lxchim(ihe4-1),
	1 1.d0-Lxchim(1)-Lxchim(ihe4)-Lxchim(ihe4-1),nosd,notd,fl,nit)
c         PRINT*,'on repasse ici',nit
c         WRITE(*,2000)pl,tl,xchim(1)
        ELSE		!X uniquement
         CALL eqstpc(pl,tl,Lxchim(1),1.d0-Lxchim(1)-Lz0,
	1 Lz0,nosd,notd,fl,nit)
        ENDIF
  	IF(ok)THEN
c        sorties

c        ro et ses dérivées

         ro=rho(1)
	 f=10.d0**fl
c        dérivées premieres de rho

	 drofn=rho(2)*rho(1)/f
         dropl=rho(2)/pt(2)
         dropn=(rho(1)/pt(1))*dropl
         drotn=(rho(1)/t)*(rho(3)-dropl*pt(3))
         droxn=rho(1)*(rho(4)-dropl*pt(4))*amm
         drop=dropn
         drot=drotn
         drox=droxn

c        dérivées secondes de rho et p

c        1. dérivées secondes naturelles dans le systeme ftx

	 dpffn=(pt(1)/(f*f))*((pt(5)/amm)+pt(2)*(pt(2)-1.d0))
         droffn=(rho(1)/(f*f))*((rho(5)/amm)+rho(2)*(rho(2)-1.d0))
         drottn=(rho(1)/(t*t))*((rho(8)/amm)+rho(3)*(rho(3)-1.d0))
	 dpttn=(pt(1)/(t*t))*((pt(8)/amm+pt(3)*(pt(3)-1.d0)))
         droftn=(rho(1)/(f*t))*(rho(2)*rho(3)+rho(6)/amm)
         drotxn=(rho(1)/t)*(rho(9)+rho(4)*rho(3)*amm)
	 dptxn=(pt(1)/t)*(pt(9)+pt(4)*pt(3)*amm)
         drofxn=(rho(1)/f)*(rho(7)+rho(4)*rho(2)*amm)
	 dpfxn=(pt(1)/f)*(pt(7)+pt(4)*pt(2)*amm)
         dpftn=(pt(1)/(f*t))*(pt(2)*pt(3)+pt(6)/amm)

c        2. dérivées premieres naturelles de f dans le systeme ptx

         dfpl=1.d0/pt(2)
         dftl=-pt(3)/pt(2)
         dfxl=-pt(4)/pt(2)
         dfpn=(f/pt(1))*dfpl
         dftn=(f/t)*dftl
	 dpxn=pt(4)*pt(1)*amm
	 dfxn=-dfpn*dpxn
         dfxn=f*dfxl
	 dptn=pt(3)*pt(1)/t

c        3. dérivées secondes naturelles dans le systeme ptx

         drotpn=dfpn*(droftn-dfpn*dptn*droffn+drofn*dfpn*
	1	(dpffn*dptn*dfpn-dpftn))
	 betbet=-(dpttn+2.d0*dftn*dpftn+dftn*dftn*dpffn)
	 betgam=-(dptxn+dfxn*dpftn+dftn*dpfxn+dftn*dfxn*dpffn)
	 drotxn=drotxn+dftn*drofxn+dfxn*droftn+dftn*dfxn*droffn+betgam*dropn
	 drottn=drottn+2.d0*dftn*droftn+dftn*dftn*droffn+betbet*dropn
         drott=drottn
         drotp=drotpn
         drotx=drotxn

c        énergie interne  u   ! attention jorgen calcule l'enthalpie

         unsro=1.d0/ro
	 u=ht(1)-pt(1)*unsro

c        dérivées premieres de u

         ro2=ro*ro
         psro2=pt(1)/ro2
	 dhfn=ht(2)/(f*amm)
         dhpn=dhfn*dfpn
	 dupn=dhpn-unsro+psro2*dropn
	 dhtn=(ht(3)-ht(2)*pt(3)/pt(2))/(t*amm)
	 dutn=dhtn+psro2*drotn
	 dhxn=ht(4)-ht(2)*pt(4)/pt(2)
	 duxn=dhxn+psro2*droxn
	 dux=duxn
	 dut=dutn
	 dup=dupn	

c        dérivées secondes de u

c	 1- dérivées naturelles de h dans le systeme ftx

	 dhttn=((ht(8)/amm)-ht(3))/(t*t*amm)
	 dhtxn=ht(9)/(t*amm)
	 dhffn=((ht(5)/amm)-ht(2))/(f*f*amm)
	 dhftn=ht(6)/(f*t*amm2)
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
	 dutt=duttn
	 dutp=dutpn
	 dutx=dutxn

	 delta=-t/ro*drot
	 deltap=delta*(-drop/ro+drotp/drot)
	 deltat=delta*(-drot/ro+drott/drot+1.d0/t)
	 deltax=delta*(-drox/ro+drotx/drot)	 

	 cpp=p/ro/t*delta
	 cp=dut+cpp	
	 dcpp=dutp+cpp*(-drop/ro+deltap/delta+1.d0/p)
	 dcpt=dutt+cpp*(-drot/ro+deltat/delta-1.d0/t)
	 dcpx=dutx+cpp*(-drox/ro+deltax/delta)

	 gradad=p/ro*delta/cp/t		!gradient adiabatique
	 dgradadp=gradad*(-drop/ro+deltap/delta-dcpp/cp+1.d0/p)
	 dgradadt=gradad*(-drot/ro+deltat/delta-dcpt/cp-1.d0/t)
	 dgradadx=gradad*(-drox/ro+deltax/delta-dcpx/cp)
	
	 alfa=p/ro*drop
	 beta=1.d0-aradias3*t**4.d0/p
	 gamma1=1.d0/(alfa-delta*gradad)

	ELSE
	 WRITE(*,*)'appel a etat_eff car defaillance de ceff'
	 CALL etat_eff(p,t,xchim,
	1 ro,drop,drot,drox,u,dup,dut,dux,
	2 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3 gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)

	ENDIF

	END SUBROUTINE etat_ceff

C*****************************************************************************

	SUBROUTINE bilinc(x0,y0,z0,x1,y1,z1,x2,y2,z2,x,y,z)

c------- performs bilincear interpolation of the fonction f,
c------- given on three arbitray points (x0,y0),(x1,y1),(x2,y2),
c------- where the respective fonction values are z0,z1,z2.

	IMPLICIT DOUBLE PRECISION(a-h,o-z)
	IMPLICIT INTEGER(i-n)
	
      x10=x1-x0
      x20=x2-x0
      y10=y1-y0
      y20=y2-y0
      z10=z1-z0
      z20=z2-z0
c
      det=x10*y20-y10*x20
c
      IF(det == 0.d0) GOTO 999
c
      dzdx=(z10*y20-z20*y10)/det
      dzdy=(z20*x10-z10*x20)/det
c
      z = z0 + (x-x0)*dzdx + (y-y0)*dzdy
c
      GOTO 1000
 999  PRINT 8000,x10,x20,y10,y20
      STOP
1000  RETURN
8000  FORMAT(/,'collinear points in s/r bilinc. error STOP.',
     . ' x1-x0,x2-x0,y1-y0,y2-y0 = ',/1x,1p4g15.6/)
     
      END SUBROUTINE bilinc

c*****************************************************************************

	SUBROUTINE clmnew(east,eahat,dmuc,ane,x,y,ea,npar,nitdmu)
c
c  finds Newton iteration correction to ionization fractions
c  Input: east: current value of ionization fractions
c         eahat: Coulomb-indepENDent part of ionization fractions
c         dmuc: Coulomb corrections (dmuc(i,1) contains correction,
c               dmuc(i,2) derivative wrt log rho, dmuc(i, 2+j)
c               derivative wrt ionization fraction no. j.
c         ane: electron number density (per unit mass)
c         x and y: hydrogen and helium abundances
c  Returns updated ionization fractions in ea.
c
c  Original version: 10/5/90
c
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
	IMPLICIT INTEGER(i-n)
      PARAMETER(nparw=3, nparw1=nparw+1)
      DIMENSION east(30), eahat(30), ea(30), dmuc(npar,*),
     *  rhld(nparw),w(nparw,nparw1)

C ajout Bernard PICHON de Janvier 2003

      DIMENSION WA(nparw,nparw), WB(nparw,1)
C
      COMMON/constsc/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
      COMMON/ln10c/ amm,amm2,amm3
      COMMON/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
c  set derivatives of log rho wrt ionization fractions
c
      anxfct=x*av/(amm*ane*ah)
      anyfct=y*av/(amm*ane*ahe)
c
c..      WRITE(*,*) 'east:',east(1),east(11),east(21)
c..      WRITE(*,*) 'eahat:',eahat(1),eahat(11),eahat(21)
      denom=1+east(1)
      denom=denom*denom
      rhld(1)=-anxfct/denom
c
      denom=1+east(11)*(1+east(21))
      denom=denom*denom
      rhld(2)=-anyfct*(1+2*east(21))/denom
      rhld(3)=-anyfct*east(11)*(2+east(11))/denom
c..      WRITE(*,*) 'rhld:', rhld
c..      WRITE(*,*) 'dmuc(1,.):',(dmuc(1,j),j=1,5)
c..      WRITE(*,*) 'dmuc(2,.):',(dmuc(2,j),j=1,5)
c..      WRITE(*,*) 'dmuc(3,.):',(dmuc(3,j),j=1,5)
c
c  set up linear equations
c
      DO 20 i=1,npar
      is=1+10*(i-1)
      eea=eahat(is)*exp(dmuc(i,1))
      w(i,nparw1)=eea-east(is)
c
      DO 15 j=1,npar
   15 w(i,j)=-eea*(dmuc(i,j+2)+dmuc(i,2)*rhld(j))
c
   20 w(i,i)=1+w(i,i)
c
      IF(idgeos>=4) THEN
        WRITE(*,110)
        DO 22 i=1,npar
   22   WRITE(*,120) (w(i,j),j=1,npar),w(i,nparw1)
      ENDIF
c
c  solve linear equations
c
C correction d'une erreur courante en f77 
C      CALL leq(w,w(1,nparw1),npar,1,nparw,nparw,err)
C Bernard PICHON en janvier 2003
C
      WA(:,:npar) = W(:,:npar)
      WB = W(:,(/nparw1/))
      CALL LEQ_BP(WA,WB,npar,1,nparw,nparw,err)
      W(:,:npar) = WA(:,:npar)
      W(:,(/nparw1/)) = WB 
c
c  as pure ad hoc fudge, halve correction if there are convergence
c  problems
c
      IF(nitdmu > 20) THEN
        DO 24 i=1,npar
   24   w(i,nparw1)=0.5*w(i,nparw1)
      ENDIF
c
c  set updated ea
c
      ea(:30) = east(:30)    ! BP 2000-11-14    CALL storec(east,ea,30)
      DO 30 i=1,npar
      is=1+10*(i-1)
   30 ea(is)=ea(is)+w(i,nparw1)
c
      IF(idgeos>=4) THEN
        WRITE(*,*) 'Corr. to ea', (w(i,nparw1),i=1,npar)
      ENDIF
c
      RETURN
  110 FORMAT(' equations in clmnew:')
  120 FORMAT(1p5e13.5)
      END SUBROUTINE clmnew
c
c*******************************************************************************
c
	SUBROUTINE dmpeqs
c  dumps COMMONs from s/r eqstfc
	IMPLICIT DOUBLE PRECISION (a-h, o-z)
	IMPLICIT INTEGER(i-n)
	COMMON/eqstdc/ c1(94)
	COMMON/eqsoutc/ c2(230)
	COMMON/dmudec/ c3(10),idmu

c
c  COMMON defining standard input and output
c
      COMMON/cstdio/ istdou, istdpr
c
      WRITE(istdpr,105)
      WRITE(istdpr,110) c1
      WRITE(istdpr,115) c2
      WRITE(istdpr,120) c3
      RETURN
  105 FORMAT(///' output from dmpeqs:'/1x,20(1h*))
  110 FORMAT(/' COMMON/eqstdc/:'/1p4e13.5/(10e13.5))
  115 FORMAT(/' COMMON/eqsoutc/:'/(1p10e13.5))
  120 FORMAT(/' COMMON/dmudec/:'/1p10e13.5)
      END SUBROUTINE dmpeqs
c
c******************************************************************************
c
	SUBROUTINE derive(x,y,dydx,nn,iy,idy,ipy,ipdy)

c     derive sets dydx(1,n) = first  derivative of y(1,n) w.r.t. x(n)
c     deriv2 sets dydx(1,n) = second derivative of y(1,n) w.r.t. x(n)
c                             n=1,nn
c
c     iy  is the first DIMENSION of y    in the CALLing programme
c     idy is the first DIMENSION of dydx in the CALLing programme
c
c     second order accuracy differences are used at interior points
c     at END points third order differences are used for first
c     derivative.  second derivative is obtained by quadratic
c     interpolation from the interior
c
c
c  revised on 6/9 1984 to include scaling by interval length,
c  to avoid underflow and consequent divide errors.
c
c                      ****************************************

c     notes on precision PARAMETERs:

c   The arguments ipy and ipdy are kept for consistency with previous
c   versions of routine. However they have no effect.

c  .............................................................................

c  Double precision version.
c  ++++++++++++++++++++++++

c  Dated 10/3/90.

c  Note: this DOUBLE PRECISION version of the routine has same name
c  as single precision version, but is distinguished by its file name

      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL second
      DIMENSION x(*), y(*),dydx(*)
      DATA epsil/1.d-8/
c
c
      second=.FALSE.
      GOTO 10
c
      entry deriv2(x,y,dydx,nn,iy,idy,ipy,ipdy)
      second=.true.
c
c
   10 ky=iy*ipy
      kdy=idy*ipdy
      IF(ipdy==1) GOTO 30
c
c     if dydx is REAL*8, set insignificant half to zero
      j=1-kdy
   20 DO 21 n=1,nn
      j=j+kdy
   21 dydx(j)=0.
c
   30 n1=nn-1
      n=0
      i=1
      k=1+ky
      hn =x(2)-x(1)
      xn= ABS(x(1))
      hna=ABS(hn)
      e=y(k)-y(1)
      IF(second) GOTO 37
c
c     first derivative at END points
      hn1=x(3)-x(2)
      h3=x(4)-x(3)
      nx1=0
      xn1= ABS(x(4))
c
      hna1=ABS(hn1)
      ha3=ABS(h3)
c
c  rescale differences
c
      IF(hn == 0.d0) GOTO 361
      hnin=1./hn
      hn1=hnin*hn1
      h3=hnin*h3
c
   31 hn2=1.+hn1
      hn3=hn2+h3
c     test size of intervals
      xxa=(xn+xn1)*epsil
      IF(hna < xxa .OR. hna1 < xxa .OR. ha3 < xxa) GOTO 361
      c=hn2*hn3*(hn2*((hn1+h3)*(hn3+1.d0)-2.d0*h3-hn1*hn3)-
     1 (h3*h3+hn1*hn3))
      IF(second) GOTO 34
   32 d=hn2*hn2
      dysc=hnin
      f=1.
      b=f*h3
      f=f*hn1
      hn1=hn3*hn3
      a=d*hn1*h3
      b=-hn1*(b+f)
      d=d*f
      IF(n) 35,35,33
   33 c=-c
      GOTO 36
   34 a=-hn2*hn3*h3*(hn2+hn3)
      b=hn3*(hn1+h3)*(hn3+1.)
      d=-hn1*hn2*(1.+hn2)
      c=0.5*c
      dysc=hnin*hnin
      IF(n) 35,35,36
   35 j=1
      l=k+ky*2
   36 dydx(i)=dysc*(a*e+b*(y(k+ky)-y(j))+d*(y(l)-y(j)))/c
      GOTO 362
c     zero interval. WRITE diagnostics
  361 nj1=nx1+1
      nj2=nx1+4
      WRITE(*,1000) nj1,nj2,(x(nx1+j),j=1,4)
      dydx(i)=0.
c
  362 IF(n /= 0) RETURN
c
c     derivative at interior points
   37 DO 42 n=2,n1
      i=i+kdy
      j=k
      k=k+ky
      xn1=xn
      xn= ABS(x(n))
      xxa=(xn1+xn)*epsil
      d=e
      e=y(k)-y(j)
      hn1=hn
      hn=x(n+1)-x(n)
      hna1=hna
      hna=ABS(hn)
c     test size of intervals
      IF(hna>=xxa.and.hna1>=xxa) GOTO 371
      dydx(i)=dydx(i-kdy)
      nj1=n-1
      nj2=n+1
      WRITE(*,1000) nj1,nj2,x(nj1),x(n),x(nj2)
      GOTO 42
c
c  rescale differences
c
  371 hnin=1./hn
      hn1=hnin*hn1
c
      c=hn1*(1.+hn1)
      IF(second) GOTO 39
   38 a=hn1*hn1
      b=1.
      dysc=hnin
      GOTO 40
   39 a=hn1
      b=-1.
      c=0.5*c
      dysc=hnin*hnin
   40 dydx(i)=dysc*(a*e+b*d)/c
   42 CONTINUE
c
      h3=x(nn-2)-x(nn-3)
c
      ha3=ABS(h3)
      h3=hnin*h3
c
      xn= ABS(x(nn))
      xn1= ABS(x(nn-3))
      nx1=nn-4
      IF(second) GOTO 50
c
c     storage indices for first derivative at last point
      i=i+kdy
      j=k
      k=k-ky*3
      l=k
      e=-e
      GOTO 31
c
c     second derivative at END points
   50 j=i
      k=i-kdy
      l=k-kdy
      i=i+kdy
   51 a=1.+hn1
      b=a+h3
      c=hn1+h3
      xxa=(xn1+xn)*epsil
c     test size of intervals
      IF(hna>=xxa.and.hna1>=xxa.and.ha3>=xxa) GOTO 52
51100 nj1=nx1+1
      nj2=nx1+4
      WRITE(*,1000) nj1,nj2,(x(nx1+j),j=1,4)
      dydx(i)=0.
      GOTO 53
   52 dydx(i)=(a*b/(hn1*c))*dydx(j)-(b/(hn1*h3))*dydx(k)
     1 +(a/(c*h3))*dydx(l)
   53 IF(i == 1) RETURN
      i=1
      j=i+kdy
      k=j+kdy
      l=k+kdy
      hn=x(2)-x(1)
      hn1=x(3)-x(2)
      h3=x(4)-x(3)
      xn= ABS(x(1))
      xn1= ABS(x(4))
      nx1=0
      hna=ABS(hn)
      hna1=ABS(hn1)
      ha3=ABS(h3)
c
c  rescale differences
c
      IF(hn == 0.d0) GOTO 51100
      hnin=1./hn
      hn1=hnin*hn1
      h3=hnin*h3
c
      GOTO 51
c
 1000 FORMAT(' **** from derive: degeneracy among x(',i5,') - x(',
     .  i5,') = ',4es16.8)
      END SUBROUTINE derive
c
c*******************************************************************
c
	SUBROUTINE eostst(rhl, tl, xh, yh, zh, nosd, notd, fl, nit,
     *  ids1, ids2)
c
c  routine to CALL EOS at given log rho = rl and log T = tl, with
c  output of various quantities to file, on units ids1 and ids2
c
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
	IMPLICIT INTEGER(i-n)
       PARAMETER(nspe=6, npar = 3)
	LOGICAL nosd, notd
      DIMENSION sn(nspe), df4(npar), d2f4i(npar,npar), d2f4f(npar,npar),
     .          d2f4x(npar),d2f4r(npar),d2f4t(npar),
     .          dp4i(npar),dh4i(npar)
c
      COMMON/eqscntc/ anh0,anhe0,ihvz,iprrad,ihmin
      COMMON/ln10c/ amm,amm2,amm3
      COMMON/constsc/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
	COMMON/dmudec/ dmu,dmuf,dmut,dmux,dmuff,dmuft,dmufx,dmutt,dmutx,
	1	dmuxx,idmu
      COMMON/eqsoutc/ east(30),xii(30),dne(10),dph(20),he(20),pe(10),
     *  hi(20),pi(10),hh(20),ph(10),hr(20),pr(10),pcoul(10), hcoul(10)
      COMMON/eqstdc/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
      DATA idsp1, idsp2 /-1, -1/
c
c  test for writing header
c
      IF(ids1 /= idsp1) THEN
        WRITE(ids1,110)
        idsp1=ids1
      ENDIF
c
      IF(ids2 /= idsp2) THEN
        WRITE(ids2,115)
        idsp2=ids2
      ENDIF
c
      CALL eqstrh(rhl, tl, xh, yh, zh, nosd, notd, fl, nit)
c
c  test for failure of convergence
c
      IF(fl < -1.e10) RETURN
c
c  prepare for CALLing WD f4 routine
c
      anx=xh*av/ah
      any=yh*av/ahe
c
      sn(1)=(1-xii1(1))*anx
      sn(2)=xii1(1)*anx
      sn(3)=(1-xii1(2)-xii1(3))*any
      sn(4)=xii1(2)*any
      sn(5)=xii1(3)*any
      sn(6)=ane(1)
c
      CALL f4der(rhl,tl,sn,nspe,f4,e4,p4,h4,d2f4r2,d2f4t2,
     1 d2f4rt,df4,d2f4i,d2f4f,d2f4x,d2f4r,d2f4t,
     2 dp4i,dp4dr,dp4dt,dp4dx,dh4i,dh4dr,dh4dt,dh4dx,npar)
c
      tk=ck2*10.d0**tl
c
      dmuc1=df4(1)/tk
      dmuc2=df4(2)/tk
      dmuc3=df4(3)/tk
c
      WRITE(3,120) rhl,tl,(xii1(i),i=1,3),dmu,dmuc1,dmuc2,dmuc3
c
c  set comparison quantities
c
      chit=pt(3)-pt(2)*rho(3)/rho(2)
      chirho=pt(2)/rho(2)
c
      WRITE(4,125) rhl, tl, pt(1), chit, chirho, cp(1)
c..      WRITE(*,130) tl,xii1,dph(1),pe(1),pr(1),pcoul(1), pt(1)
  130 FORMAT('extra output:',f8.3,1p9e10.3)
c
      RETURN
  110 FORMAT('# log rho, log T, xi(1-3), dmu, dmu1 - dmu3:'/'#')
  115 FORMAT('# log rho, log T, log p, chi_T, chi_rho, c_p:'/'#')
  120 FORMAT(2f10.5,1p3e11.3,4e12.4)
  125 FORMAT(2f12.7,1pe14.7,0p2f12.7,1pe14.7)
  
      END SUBROUTINE eostst
c
c******************************************************************************
c
	SUBROUTINE eqstfc(fl,tl,x,y,z,nosd,notd)
c
c  This version of the EFF routine includes
c  effects of Coulomb interaction in the Debye-Hueckel approximation,
c  using the F4 from the MHD package.
c
c  Modifications initiated 19/4/90.
c
c  jc-d equation of state updated to be compatible with  dog equation
c  of state routines. notice in particular that COMMON/eqscntc/ has
c  been changed, by the inclusion of anhe0 after anh0. in addition  a
c  dummy sbroutine seteqs has been added. furthermore the DIMENSIONs
c  in COMMON /potetc/ have been changed, to allow for only up to
c  10 heavy elements.
c
c  controls:
c  ********
c
c  in argument list:
c
c  LOGICAL nosd and notd are switches for no calculation of second and
c  third derivatives, respectively.
c
c  in COMMON/eqscntc/:
c
c  ihvz determines treatment of heavy elements.
c  ihvz = 0: assume heavy elements to be fully ionized everywhere
c       = 1: complete treatment of ionization of c and o; first level
c            of ionization of fe. reset abundances to fit results of
c            full treatment of all 10 elements considered.
c       = 2: first level of ionization of all 10 elements.
c       = 3: complete treatment of c, n and o; first level of
c            ionization of remaining elements.
c       = 4: complete treatment of ionization of all 10 elements
c            included.
c
c  iprrad  /=  0: include pressure and enthalpy of radiation (in
c                 diffusion approximation)
c  iprrad  =   0: DO not include pressure and enthalpy of radiation.
c
c  ihmin   =   1: include h-
c          /=  1: DO not include h-
c
c  these are initialized  in module mod_bleqstc to ihvz = 1,
c  iprrad = 1, ihmin = 0.
c
c  in COMMON/eqdpco/, relating to fudge departure coefficients to
c  mimick nlte effects:
c
c  idpco = 0: DO not include departure coefficients
c        = 1: determine departure coefficient for h (RETURNed in
c             bdcoh) such that n(h+)/n(h neutr.) = frhi.
c        = 2: treat h as for idpco = 1. take departure coefficient
c             for he and heavy elements to be bdcoz.
c        = 3: take departure coefficient for h to be bdcoh, and
c             departure coefficient for he and heavy elements to be
c             bdcoz.
c
c  frhi, bdcoh and bdcoz are initialized to 1., and idpco to 0.
c
c  addition 5/1/84: flag iomfll in COMMON/hvomcl/.
c    if iomfll = 0 statistical weight omegac for fully ionized state
c    of heavy elements is set to 15.
c    this corresponds to the situation before 5/1/84. the
c    origin of this error is unclear.
c    if iomfll  /=  0 omegac for fully ionized state is set to 1,
c    as it should be.
c
c  Controls in COMMON /ccoulm/
c    epsdmu: convergence criterion for iteration for Coulomb
c    effect in Saha equation. (Default: 1.e-12)
c    icoulm: determines how Coulomb terms are included:
c      icoulm = -1: DO not include Coulomb effects (should correspond
c         to basic eqstf).
c      icoulm = 0: include Coulomb effects fully.
c      icoulm = 1: include Coulomb effects in Saha equations, but not
c         explicit contributions to pressure and enthalpy.
c      icoulm = 2: include explicit contributions to pressure and
c         enthalpy from Coulomb effects, but not effect
c         in Saha equations.
c      (Default: 0)
c    iclmit: determines type of Coulomb iteration
c      iclmit = 0: backsubstitution
c      iclmit = 1: Newton iteration
c      (Default: 1)
c    iclmsd; if iclmsd = 1, set include second derivatives of
c      Coulomb terms
c      (Default: 0)
c
c
c  modification 15/9/86: SAVE statement added in all sbroutines.
c     (needed in f77, apparently)
c
c  modification 3/1/88: DIMENSION of dxt1 increased to 20.
c
c  Modification 6/8/88: second derivative of He ionization fractions
c     corrected. Also in commom/eqsout/ is now stored ea before rescaling,
c     for test of derivatives.
c
c  Modification 17/5/90: include possibility of including detailed
c     ionization of heavy elements, for ihvz = 4. Note that this
c     involves changing size of arrays in COMMONs /hvredu/ and /xhv/
c
c  Modified 22/5/90 to include all levels of ionization of Ar and Fe.
c     Previously only the first 15 levels were included, the remainder
c     being forced to be unionized
c
c  Modification 4/6/90: include average degree of ionization for
c     each element in array xhvmn(10) at the beginning of
c     COMMON/hviond/
c
c  Modified 4/6/90, adding array gmm1(4) in COMMON/ eqstd/ containing
c     Gamma1 and derivatives.
c
c  Modified 4/6/90, to include numerical second derivatives of
c     Coulomb terms.
c
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 9/3/90
c
c            *********************************

c  results RETURNed in COMMON/eqstdc/. in addition value of fl
c  is set into COMMON/eqstfl/.

	IMPLICIT DOUBLE PRECISION (a-h,o-z)
	IMPLICIT INTEGER(i-n)
	INTEGER, PARAMETER :: nspe = 6, npar = 3
	LOGICAL :: tstl,nosd,notd,cmplio, secder, thrder
	DIMENSION phi(30),hst(10),ex(3),ea(30),xi(30),dxt(4),dxt1(20),
	1 dxt2(4),anuh(10),ueh(10),anuhr(23),uehr(23),
	2 ehm(10),aneder(10), sn(nspe), dmucpr(npar),
	3 eahat(30), dmucc(npar,4), dmucp(npar,4,6),
	4 pcoulc(4), pcoulp(4,6), hcoulc(4), hcoulp(4,6)
	COMMON/dmudec/ dmu,dmuf,dmut,dmux,dmuff,dmuft,dmufx,dmutt,dmutx,
	1 dmuxx,idmu
	COMMON/eqscntc/ anh0,anhe0,ihvz,iprrad,ihmin
	COMMON/ln10c/ amm,amm2,amm3
	COMMON/constsc/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
	1 ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm

c  note: second DIMENSION of dmuc must be at least as great as npar + 2

	COMMON/dmucdr/ dmuc(npar,10)
	COMMON/df4der/ df4(npar), d2f4i(npar,npar), d2f4f(npar,npar),
	1 d2f4x(npar),d2f4r(npar),d2f4t(npar),
	2 dp4i(npar),dh4i(npar), p4, h4, dp4dr, dp4dt, dp4dx,
	3 dh4dr, dh4dt, dh4dx
      COMMON/eqsoutc/ east(30),xii(30),dne(10),dph(20),he(20),pe(10),
     *  hi(20),pi(10),hh(20),ph(10),hr(20),pr(10), pcoul(10), hcoul(10)
      COMMON/eqstdc/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
      COMMON/eqstfl/ flcm
      COMMON/eqdpco/ frhi,bdcoh,bdcoz,idpco
      COMMON/xhminc/ xhm(10)

c  controls for Coulomb effect.
c  epsdmu is accuracy requirement in backsubstitution
c  icoulm is used for selectively switching off parts of Coulomb effect

      COMMON/ccoulm/ epsdmu, icoulm, iclmit, iclmsd, nitdmu, epssdr
c
c  COMMON defining standard input and output
c
      COMMON/cstdio/ istdou, istdpr
      COMMON/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      COMMON/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng

      EQUIVALENCE(ex(1),exh)
c
      DATA flprev, tlprev, xprev /-1.d10, -1.d10, -1.d10/
c
      SAVE
c
      fxp(a)=exp(min(85.d0,max(a,-85.d0)))
      nuu(l,n)=((l-1)*(8-l))/2+n-l+5
c
c  *****************************************************************
c
c  set equation of state version number, depENDing on icoulm
c
      IF(icoulm==0) THEN
        ivreos=1
      ELSE IF(icoulm==2) THEN
        ivreos=2
      ENDIF
c
c  *****************************************************************
c
c  for changing internal logics, set LOGICAL variables for setting
c  second and third derivatives
c
      secder = .not. nosd
      thrder = .not. notd
c
c  set y
      y=1-x-z
c  store fl in COMMON
      flcm=fl
c  number of derivatives
      IF(thrder) THEN
        ider=20
      ELSE IF(secder) THEN
        ider=10
      ELSE
        ider=4
      ENDIF
      iders=min0(ider,10)
c
      iclsdr = 0
c
c  entry point for setting numerical second derivatives of
c  Coulomb terms
c
    5 f=1.d1**fl
      t=1.d1**tl
c
      CALL phderc(fl,tl,phi,hst,nosd,notd)
c
c
c  psi and derivatives
      wf= sqrt(1+f)
      psi=2*wf+log(f/((1+wf)*(1+wf)))
      psif=amm*wf
      psff=amm2*f/(2*wf)
c
      zf=x+2*y+anhe0*z
      zf2=zf*zf
      zf3=zf*zf2
c
      ak0=ckh*zf2
c
c  set k*t, in ev and ergs
c
      tkev=ck1*t
      tkergs=ck2*t
c  a0**3*ne/(k*t)
      aa=caa*phi(1)/(tkev*zf3)
c
c  delta mu and derivatives
c
      bmu=tkev+20*ak0
      dmu=aa*bmu
c
      ref=phi(2)
      ret=phi(3)
      dmuf=dmu*ref*amm
      dmut=aa*(bmu*ret-20*ak0)*amm
      dmux=aa*(3*tkev+20*ak0)/zf
c
      IF(secder) THEN
        dmuff=(dmuf*ref+dmu*phi(4)*amm)*amm-psff
        dmuft=(dmut*ref+dmu*phi(5)*amm)*amm
        dmufx=dmux*ref*amm
        dmutt=aa*(20*ak0*(1-2*ret)+bmu*(phi(6)+ret*ret))*amm2
        dmutx=aa*(ret*(3*tkev+20*ak0)-20*ak0)*amm/zf
        dmuxx=aa*(12*tkev+40*ak0)/zf2
      ENDIF
c
      dmu=dmu-psi
      dmuf=dmuf-psif
      idmu=1
c
c  ***************************************
c
c  ionization of heavy elements
c
      IF(ihvz <= 0) THEN
        anuh(1)=anh0
        !
        ! BP 2000-11-14
	  ! TRES dangeureux car, a priori, aucun sens valable en Fortran
	  !
	  anuh(2:10) = 0.
	  ueh(:10) = 0. 
c
      ELSE
c
        CALL hvionac(fl,tl,x,y,z,nosd,notd,anuh,ueh,anuhr,uehr)
c
      ENDIF
c
c  test for complete ionization of h and he
c
      cmplio=(dmu-54.4/tkev)>19
c  test for use of departure coefficient
      IF(idpco >=1) cmplio=cmplio.and.frhi>1.e6
      IF(cmplio) THEN
c
c  complete ionization
c
        xi(1) =1
        xi(11)=0
        xi(21)=1
        DO 10 i=2,22,10
        jj=i+8
        DO 10 j=i,jj
   10   xi(j)=0
c
c  to avoid problems with logics in Coulomb term treatment,
c  initialize iclder
c
        iclder = 0
c
        GOTO 30
c
      ENDIF
c
c  A note on logics of Coulomb iteration: iclder = 0 flags for
c  initial passes iterating for Coulomb correction, and iclder = 1
c  flags for final pass setting derivatives
c  Similarly, when setting second derivatives numeriCALLy, iclsdr
c  flags for modifications of variables
c
c  e h, e he, e he+ and derivatives
c
c  test for initializing Coulomb terms or using previous values
c  Note: for icoulm = 2 switch off effect of Coulomb term on ionization
c  For icoulm = -1 switch off effect of Coulomb term entirely
c  Initialization also depENDs on whether Newton iteration or
c  backsubstitution is used
c
      IF(icoulm==2.OR.icoulm==-1) THEN
        dmuc = 0.   ! BP 2000-11-14    a peu pres ca !!!   CALL zeroc(dmuc,npar*10)
        dmucpr = 0. ! BP 2000-11-14    a peu pres ca !!!   CALL zeroc(dmucpr,npar)
      ELSE
c
c  initialize for Coulomb iteration. Switch off second derivatives
c  for iteration passes.
c
        iclder=0
        secder=.FALSE.
        dmuc = 0.     ! BP 2000-11-14    a peu pres ca !!!   CALL zeroc(dmuc,npar*10)
        IF(ABS(fl-flprev)>0.1.OR.ABS(tl-tlprev)>0.1
     *    .OR.ABS(x-xprev)>0.1) THEN
          dmucpr = 0. ! BP 2000-11-14    a peu pres ca !!!   CALL zeroc(dmucpr,npar)
        ENDIF
      ENDIF
c
      flprev=fl
      tlprev=tl
      xprev=x
      nitdmu=0
c
c  entry point for setting Coulomb derivatives
c
   15 k=-10
c..      WRITE(*,*) 'iclsdr, fl, tl, x', iclsdr, fl, tl, x
      DO 25 ia=1,3
      k=k+10
      ext=ex(ia)/tkev
      eea=dmu-ext+dmuc(ia,1)
      IF(eea+30 <= 0) THEN
c
c  no ionization
c
        DO 16 i=1,10
   16   ea(k+i)=0.e0
        GOTO 25
c
      ENDIF
c
      eea=fxp(eea)
c
c  set statistical weights
c
      IF(ia /= 2) THEN
        eea=eea/2
      ELSE
        eea=2*eea
      ENDIF
c
c  test for departure coefficient for h or he
c
      IF(ia==1.and.idpco>0) THEN
        IF(idpco <= 2) bdcoh=frhi/eea
        eea=bdcoh*eea
      ELSE IF(ia>=2.and.idpco>=2) THEN
        eea=bdcoz*eea
      ENDIF
c
      ea(k+1)=eea
c  first derivatives
      ea(k+2)=(dmuf+dmuc(ia,2))*eea
      ea(k+3)=(amm*ext+dmut+dmuc(ia,3))*eea
      ea(k+4)=(dmux+dmuc(ia,4))*eea
      IF(secder) THEN
c  second derivatives, with Coulomb contribution
        ea(k+5)=(dmuff+dmuc(ia,5))*eea+(dmuf+dmuc(ia,2))*ea(k+2)
        ea(k+6)=(dmuft+dmuc(ia,6))*eea+(dmuf+dmuc(ia,2))*ea(k+3)
        ea(k+7)=(dmufx+dmuc(ia,7))*eea+(dmux+dmuc(ia,4))*ea(k+2)
        ea(k+8)=(dmutt+dmuc(ia,8)-amm2*ext)*eea
     *         +(amm*ext+dmut+dmuc(ia,3))*ea(k+3)
        ea(k+9)=(dmutx+dmuc(ia,9))*eea+(dmux+dmuc(ia,4))*ea(k+3)
        ea(k+10)=(dmuxx+dmuc(ia,10))*eea+(dmux+dmuc(ia,4))*ea(k+4)
      ENDIF
c
   25 CONTINUE
c
c  test for setting up for Coulomb iteration
c
      eahat(:30) = ea(:30)    ! BP 2000-11-14    CALL storec(ea,eahat,30)
c
c  in iteration pass, include dmuc from previous case in ea
c  unless ignoring Coulomb effect on ionization
c
      IF(iclder /= 1.and.icoulm /= -1.and.icoulm /= 2) THEN
        DO 26 i=1,3
        is=1+10*(i-1)
   26   ea(is)=eahat(is)*exp(dmucpr(i))
      ENDIF
c
c  entry point for continuing Coulomb iteration by backsubstitution
c  For subsequent iterations, include Coulomb term
c
   27 IF(iclder /= 1.and.nitdmu > 0) THEN
        DO 28 i=1,3
        is=1+10*(i-1)
   28   ea(is)=eahat(is)*exp(dmuc(i,1))
      ENDIF
c
c  entry point for continuing Coulomb iteration by Newton iteration
c  for diagnostic or iteration purposes, store current ea
c
   29 east(:30) = ea(:30)    ! BP 2000-11-14    CALL storec(ea,east,30)
c
c  set degrees of ionization of H and He
c
      CALL hheion(ea(1), ea(11), ea(21), xi(1), xi(11), xi(21), secder)
c
c  inclusion of h-
c
   30 IF(ihmin==1) THEN
        CALL hmnion(tl, ea(1), ehm, xhm, secder)
      ELSE
        xhm(:10) = 0.   ! BP 2000-11-14     CALL zeroc(xhm,10)
      ENDIF
c
c  ne and derivatives
c
c  combine he fractions for an in xi(20+.), and for ionization
c  energy in xi(10+.)
c
      exhc=exhe+exhep
c
      IF(secder) THEN
        imx=10
      ELSE
        imx=4
      ENDIF
      !
      ! BP 2000-11-14    
	!
!!      DO 35 i=1,21,10
!!   35 CALL storec(xi(i),xii(i),imx)
      !
      xii(01:imx+00) = xi(01:imx+00)
      xii(11:imx+10) = xi(11:imx+10)
	xii(21:imx+20) = xi(21:imx+20)
      xia=xi(1)
      imx=imx+20
c
      DO 36 i=21,imx
      i10=i-10
      xio=exhe*xi(i10)+exhc*xi(i)
      xi(i)=2*xi(i)+xi(i10)
   36 xi(i10)=xio
c
      ider2=min0(ider,10)
      IF(ihmin==1) THEN
c
c  in the case with H- there is a risk of negative (unphysical)
c  Ne. Write error message and reset to values without
c  combine h and h-
c
        dxt1(1)=xi(1)-xhm(1)
        ane(1)=av*(dxt1(1)*x/ah+xi(21)*y/ahe+z*anuh(1))
        IF(ane(1) <= 0) THEN
          WRITE(istdpr,1010) xi(1), xhm(1), anuh(1)
          DO 37 l=1,ider2
   37     dxt1(l)=xi(l)
          ane(1)=av*(dxt1(1)*x/ah+xi(21)*y/ahe+z*anuh(1))
        ELSE
          DO 38 l=2,ider2
   38     dxt1(l)=xi(l)-xhm(l)
        ENDIF
      ELSE
c
c  no H-
c
        DO 39 l=1,ider2
   39   dxt1(l)=xi(l)
c
        ane(1)=av*(dxt1(1)*x/ah+xi(21)*y/ahe+z*anuh(1))
c
      ENDIF
c
c  terms needed for x-derivatives
c
      DO 40 l=1,4
   40 dxt(l)=av*(dxt1(l)/ah-xi(20+l)/ahe)
c
      ii=4
      DO 44 l=1,3
      l1=l+1
      anel=(dxt1(l1)*x/ah+xi(l+21)*y/ahe+z*anuh(l1))*av
      tstl=l==3
      IF(tstl) anel=anel+dxt(1)
      ane(l1)=anel
      IF(secder) THEN
        DO 42 m=l,3
        m1=m+1
        ii=ii+1
        anelm=(dxt1(ii)*x/ah+xi(20+ii)*y/ahe+z*anuh(ii))*av
        IF(tstl) anelm=anelm+dxt(m1)
        IF(m==3) anelm=anelm+dxt(l1)
   42   ane(ii)=anelm
      ENDIF
   44 CONTINUE
c
c  as auxiliary quantities, to avoid over- and underflow, set
c  derivatives of ne divided by ne
c
      DO 46 i=2,iders
   46 aneder(i)=ane(i)/ane(1)
c
c  the density and derivatives (note that rho(2) = dlog rho/dlog f,
c  and so on).
c
      anee=ane(1)
      rho(1)=phi(1)*(crho/anee)
      ii=4
      jj=10
      DO 50 l=1,3
      l1=l+1
      rhol=-aneder(l1)/amm
      tstl=l <= 2
      IF(tstl) rhol=rhol+phi(l1)
      rho(l1)=rhol
      IF(secder) THEN
        DO 48 m=l,3
        ii=ii+1
        m1=m+1
        lm=l+m
        rholm=(aneder(l1)*aneder(m1)-aneder(ii))/amm
        IF(tstl.and.m <= 2) rholm=rholm+amm*phi(2+lm)
        IF(thrder) THEN
          DO 47 n=m,3
          jj=jj+1
          rhd=-2*aneder(l1)*aneder(m1)*aneder(n+1)/amm
          IF(l<3.and.m<3.and.n<3) rhd=rhd+amm2*phi(4+lm+n)
   47     rho(jj)=rhd
        ENDIF
   48   rho(ii)=rholm
      ENDIF
   50 CONTINUE
c
c  test for skipping Coulomb terms entirely
c  Also skip this bit for final pass when setting second derivatives
c
      IF(icoulm==-1.OR.iclsdr==7) GOTO 60
c
c  set Coulomb terms. Prepare for CALLing WD f4 routine
c  Note: we currently skip this in the pass for setting derivatives,
c  unless for full ionization, where Coulomb part has not been passed
c  previously
c  If not setting derivatives in f4der is implemented in future,
c  this may require further thought
c
      IF(iclder /= 1.OR.cmplio.OR.icoulm==2) THEN
c
        anx=x*av/ah
        any=y*av/ahe
c
        sn(1)=(1.d0-xii(1))*anx
        sn(2)=xii(1)*anx
        sn(3)=(1.d0-xii(11)-xii(21))*any
        sn(4)=xii(11)*any
        sn(5)=xii(21)*any
c
c  to avoid rounding error problems in subtractions in sn(1) and sn(3),
c  reset them from the ea-s, or to zero for full ionization
c
        IF(cmplio) THEN
          sn(1)=0
          sn(3)=0
        ELSE
          sn(1)=anx/(1+east(1))
          sn(3)=any/(1+east(11)*(1+east(21)))
        ENDIF
c
c  number of free electrons. For consistency include only those
c  coming from H and He
c
        sn(6)=xii(1)*anx+(xii(11)+2*xii(21))*any
c
        rhl=log10(rho(1))
c        IF(idgeos>=3) WRITE(istdpr,*) 'CALLing f4der at',rhl,tl,sn
c
        CALL f4der(rhl,tl,sn,nspe,f4,e4,p4,h4,d2f4r2,d2f4t2,
     .             d2f4rt,df4,d2f4i,d2f4f,d2f4x,d2f4r,d2f4t,
     .             dp4i,dp4dr,dp4dt,dp4dx,dh4i,dh4dr,dh4dt,dh4dx,npar)
c
        IF(icoulm /= 2) THEN
          dmuc(1,1)=-df4(1)/tkergs
          dmuc(2,1)=-df4(2)/tkergs
          dmuc(3,1)=-df4(3)/tkergs
c
c  test for continuing Coulomb iteration
c
          ddmuc1=dmuc(1,1)-dmucpr(1)
          ddmuc2=dmuc(2,1)-dmucpr(2)
          ddmuc3=dmuc(3,1)-dmucpr(3)
c          IF(idgeos>=2)
c	WRITE(*,*) ' ddmuc1-3:', ddmuc1, ddmuc2, ddmuc3
c        ELSE
          ddmuc1=0
          ddmuc2=0
          ddmuc3=0
        ENDIF
c
      ENDIF
c
c..      WRITE(*,*) 'iclder =',iclder
c
c  test for convergence
c
      IF(.not.cmplio.and.icoulm /= 2.and.iclder==0.and.
     *   (ABS(ddmuc1)>=epsdmu.OR.
     *    ABS(ddmuc2)>=epsdmu.OR.
     *    ABS(ddmuc3)>=epsdmu)) THEN
c
c  test for failure to converge
c
        IF(nitdmu > 20) THEN
          WRITE(istdpr,1100) fl, tl, x, z,
     *      ddmuc1, ddmuc2, ddmuc3
          IF(istdou /= istdpr) WRITE(istdou,1100) fl, tl, x, z,
     *      ddmuc1, ddmuc2, ddmuc3
          STOP
        ENDIF
c
c
        nitdmu=nitdmu+1
c
c  store previous value
c
        DO 51 i=1,npar
   51   dmucpr(i)=dmuc(i,1)
c
c  test for backsubstitution or Newton iteration
c
        IF(iclmit==0) THEN
          GOTO 27
        ELSE
c
c  store derivatives in dmuc
c
c..          WRITE(*,*) 'd2f4f(1,.):',(d2f4f(1,j),j=1,3)
c..          WRITE(*,*) 'd2f4f(2,.):',(d2f4f(2,j),j=1,3)
c..          WRITE(*,*) 'd2f4f(3,.):',(d2f4f(3,j),j=1,3)
          DO 52 i=1,npar
          dmuc(i,2)=-d2f4r(i)/tkergs
          DO 52 j=1,npar
   52     dmuc(i,j+2)=-d2f4f(i,j)/tkergs
c
          IF(idgeos>=4) WRITE(*,*) 'xi(1-3):',xii(1), xii(11), xii(21)
          CALL clmnew(east,eahat,dmuc,ane(1),x,y,ea,npar,nitdmu)
          GOTO 29
        ENDIF
c
c
c  test for needing last pass of ionization section, to set
c  derivatives
c
      ELSE IF(.not.cmplio.and.icoulm /= 2.and.iclder==0) THEN
c
c  set derivatives of delta mu c
c
        DO 54 i=1,npar
        DO 53 j=2,4
        dmuc(i,j)=0
        DO 53 k=1,npar
   53   dmuc(i,j)=dmuc(i,j)-d2f4i(i,k)*xii(10*(k-1)+j)
c
c  add explicit contributions from derivatives of
c  dmuc wrt log rho, log T and X at fixed xi
c  Note: last term in dmuc(i,3) corrects for division by kT.
c
        dmuc(i,2)=dmuc(i,2)-d2f4r(i)*rho(2)
        dmuc(i,3)=dmuc(i,3)-d2f4r(i)*rho(3)-d2f4t(i)+amm*df4(i)
c..        WRITE(*,*) 'dmucx:',i,dmuc(i,4),-d2f4r(i)*rho(4),-d2f4x(i)
        dmuc(i,4)=dmuc(i,4)-d2f4r(i)*rho(4)-d2f4x(i)
c
        DO 54 j=2,4
   54   dmuc(i,j)=dmuc(i,j)/tkergs
c
        iclder=1
c
c  restore secder, unless a later pass will be made to compute
c  second derivatives
c
        secder=.not.nosd
c
        GOTO 15
c
      ELSE
c
c  otherwise store pressure and enthalpy from Coulomb effect,
c  transforming derivatives
c
        pcoul(1)=p4
        hcoul(1)=h4
c
c..        WRITE(*,*) 'setting pcoul(1) =',pcoul(1)
c
      !
      ! BP 2000-11-14
      ! TRES dangeureux car, a priori, aucun sens valable en Fortran
      !
        pcoul(2:10) = 0.     !  CALL zeroc(pcoul(2),9)
        hcoul(2:10) = 0.     !  CALL zeroc(hcoul(2),9)
c
        DO 56 j=2,4
        DO 56 k=1,npar
        pcoul(j)=pcoul(j)+dp4i(k)*xii(10*(k-1)+j)
   56   hcoul(j)=hcoul(j)+dh4i(k)*xii(10*(k-1)+j)
c
c  add explicit contributions from derivatives of
c  dmuc wrt log rho, log T and X at fixed xi
c
        pcoul(2)=pcoul(2)+dp4dr*rho(2)
        pcoul(3)=pcoul(3)+dp4dr*rho(3)+dp4dt
        pcoul(4)=pcoul(4)+dp4dr*rho(4)+dp4dx
        hcoul(2)=hcoul(2)+dh4dr*rho(2)
        hcoul(3)=hcoul(3)+dh4dr*rho(3)+dh4dt
        hcoul(4)=hcoul(4)+dh4dr*rho(4)+dh4dx
c..        WRITE(*,*) 'In eqstf dh4dr, dh4dt =',dh4dr, dh4dt
      ENDIF
c
c  test for setting second derivatives of Coulomb terms numeriCALLy
c
      IF(iclmsd==1.and..not.nosd) THEN
c
        IF(iclsdr==0) THEN
          flc=fl
          tlc=tl
          xc=x
	    !
	    ! BP 2000-11-14      TRES mauvaise programmation : DANGEUREUX
	    !
          ! IF(icoulm /= 2) CALL storec(dmuc,dmucc,4*npar)
          IF(icoulm /= 2) dmucc(:,1:4) = dmuc(:,1:4)
          pcoulc(:4) = pcoul(:4) ! BP 2000-11-14  CALL storec(pcoul,pcoulc,4)
          hcoulc(:4) = hcoul(:4) ! BP 2000-11-14  CALL storec(hcoul,hcoulc,4)
c
          fl=flc-epssdr
          iclsdr=1
          GOTO 5
c
        ELSE
          IF(icoulm /= 2) CALL storec(dmuc,dmucp(1,1,iclsdr),4*npar)
          CALL storec(pcoul,pcoulp(1,iclsdr),4)
          CALL storec(hcoul,hcoulp(1,iclsdr),4)
c
          IF(iclsdr <= 5) THEN
            IF(iclsdr==1) THEN
              fl=flc+epssdr
            ELSE IF(iclsdr==2) THEN
              fl=flc
              tl=tlc-epssdr
            ELSE IF(iclsdr==3) THEN
              tl=tlc+epssdr
            ELSE IF(iclsdr==4) THEN
              tl=tlc
              x=xc-epssdr
              y=1-x-z
            ELSE IF(iclsdr==5) THEN
              x=xc+epssdr
              y=1-x-z
            ENDIF
c
            iclsdr=iclsdr+1
            GOTO 5
c
          ELSE
c
            x=xc
            y=1-x-z
c
            IF(icoulm /= 2) CALL storec(dmucc,dmuc,4*npar)
            CALL storec(pcoulc,pcoul,4)
            CALL storec(hcoulc,hcoul,4)
c
c  set second derivatives
c
            epsdi2=1.d0/(2*epssdr)
            pcoul(5) =(pcoulp(2,2)-pcoulp(2,1))*epsdi2
            hcoul(5) =(hcoulp(2,2)-hcoulp(2,1))*epsdi2
            pcoul(6) =(pcoulp(3,2)-pcoulp(3,1))*epsdi2
            hcoul(6) =(hcoulp(3,2)-hcoulp(3,1))*epsdi2
            pcoul(7) =(pcoulp(4,2)-pcoulp(4,1))*epsdi2
            hcoul(7) =(hcoulp(4,2)-hcoulp(4,1))*epsdi2
            pcoul(8) =(pcoulp(3,4)-pcoulp(3,3))*epsdi2
            hcoul(8) =(hcoulp(3,4)-hcoulp(3,3))*epsdi2
            pcoul(9) =(pcoulp(3,6)-pcoulp(3,5))*epsdi2
            hcoul(9) =(hcoulp(3,6)-hcoulp(3,5))*epsdi2
            pcoul(10)=(pcoulp(4,6)-pcoulp(4,5))*epsdi2
            hcoul(10)=(hcoulp(4,6)-hcoulp(4,5))*epsdi2
c
            IF(icoulm /= 2) THEN
              DO 58 i=1,npar
              dmuc(i,5) =(dmucp(i,2,2)-dmucp(i,2,1))*epsdi2
              dmuc(i,6) =(dmucp(i,3,2)-dmucp(i,3,1))*epsdi2
              dmuc(i,7) =(dmucp(i,4,2)-dmucp(i,4,1))*epsdi2
              dmuc(i,8) =(dmucp(i,3,4)-dmucp(i,3,3))*epsdi2
              dmuc(i,9) =(dmucp(i,3,6)-dmucp(i,3,5))*epsdi2
   58         dmuc(i,10)=(dmucp(i,4,6)-dmucp(i,4,5))*epsdi2
c
c  make final pass of ionization section to set second derivatives
c  correctly
c
              iclder=1
              secder=.true.
              iclsdr=7
              GOTO 15
            ENDIF
c
          ENDIF
c
        ENDIF
c
      ENDIF
 
c
c  ****************************************************
c
c  END skipping Coulomb terms
c
   60 CONTINUE
c
c  start setting total pressure and enthalpy
c
c  as initialization, zero arrays
c
      pt(:20) = 0.    ! CALL zeroc(pt,20)
      ht(:20) = 0.    ! CALL zeroc(ht,20)
c
c
c  delta p, delta h and derivatives.
c
c  test for complete ionization
      IF(cmplio) GOTO 80
c
c  delta ne**2
c
c  note: on 3/1/84 definition of dne was changed to
c        (delta ne **2)/av**2 to avoid overflows in calculation
c        of delta p and delta h on univac (with exponents limited
c        to 38). a corresponding change was made in the definition
c        of ca03 in s/r setcnsc.
c
c
c  reset xi(1) and xi(21) (note that xi(2-10) and xi(22-30) are still
c  consistent).
      xi(1)=-1/(east(1)+1)
      xi(21)=-(2+east(11))/(1+east(11)*(1+east(21)))
c
      DO 65 l=1,4
      dxt1(l)=xi(l+20)/ahe-xi(l)/ah
      dxtl=-xi(l+20)*y/ahe-xi(l)*x/ah-anuh(l)*z
      IF(l==1) dxtl=dxtl+anh0*z
      dxtl2=ane(l)
      IF(l <= 3) GOTO 64
      dxtl=dxtl+dxt1(1)
      dxtl2=dxtl2+avda
   64 dxt(l)=dxtl
   65 dxt2(l)=dxtl2/av
c
c  ne0/av:
c
      ann=ane(1)/av+(x/ah+2*y/ahe+z*anh0)
      xtt=dxt(1)
      xttiav=xtt/av
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
   68 dne(ii)=ane(ii)*xttiav+dxt2(l1)*dxt(m1)+dxt2(m1)*dxt(l1)
     .  +ann*dnlm
   70 CONTINUE
c
c  quantities COMMON to delta p and delta h (dph(15-20) is used as
c  intermediate storage).
c
c  Note: here the constant for delta p and delta H is set to C1 = ca03/2
c
      a03=ca03/zf3/2
c
      !
      ! BP 2000-11-14
      ! TRES dangeureux car, a priori, aucun sens valable en Fortran
      !
      dph(15:20) = 0.    ! CALL zeroc(dph(15),6)
      dxt1(1)=0
      c1=amm*tkev
      dxt1(2)=c1
      dph(18)=amm2*tkev
      dph(19)=3*c1/zf
      dnee=dne(1)
      a03=a03*rho(1)*rho(1)
      nu=2
      nb=20
      k1=1
   75 DO 77 l=2,4
   77 dxt(l-1)=dne(l)+nu*amm*dnee*rho(l)
c
      dxt1(3)=(3*tkev+nb*ak0)/zf
      bmu=tkev+nb*ak0
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
     .  +dnee*(rho(jj)+nu*amm*rho(m1)*rho(l1))))+dxt(l)*dxt1(m)
     .  +dxt(m)*dxt1(l)+dnee*dph(10+jj)
      IF(m==3.and.l==3) dphlm=dphlm+dnee*(12*tkev+2*nb*ak0)/zf2
   78 dph(ii)=a03*dphlm
   79 CONTINUE
      IF(nu==1) GOTO 90
      nu=1
      nb=40
      a03=a03/rho(1)
      k1=11
      GOTO 75
c  complete ionization
   80 dne(:10) = 0.   ! CALL zeroc(dne,10)
      dph(:20) = 0.   ! CALL zeroc(dph,20)
      GOTO 100
   90 CONTINUE
      DO 95 i=1,iders
      pt(i)=pt(i)+dph(i)
   95 ht(i)=ht(i)+dph(i+10)
c
  100 CONTINUE
c
c  *********************************************
c
c  electron pressure and enthalpy
c
      ii=4
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
      tstl=tstl.OR.m==3
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
     .  +anee*hst(4+lm+n)
      pt(jj)=pt(jj)+amm3*pee*(phi(11+l)*phi(11+m)*phi(11+n)
     .  +phi(12+lm)*phi(11+n)+phi(12+l+n)*phi(11+m)
     .  +phi(12+m+n)*phi(11+l)+phi(14+l+m+n))
  105 helm=che*helm
      ht(jj)=ht(jj)+helm
  106 he(jj)=helm
  108 CONTINUE
  110 CONTINUE
c
c  *********************************************
c
c  ionization enthalpy
c
c  (pi is introduced for consistency)
      !
      ! BP 2000-11-14
      ! TRES dangeureux car, a priori, aucun sens valable en Fortran
      !
      pi(:10) = 0.     ! CALL zeroc(pi,10)
      hi(11:20) = 0.   ! CALL zeroc(hi(11),10)
c
      xi(1)=xia
      averg=av*ergev
c  combine h and h-
c	Write (6,*)ider
      DO 111 l=1,ider
  111 dxt1(l)=exh*xi(l)-exhm*xhm(l)
      DO 112 l=1,4
  112 dxt(l)=dxt1(l)/ah-xi(10+l)/ahe
c
      hi1=averg*(dxt1(1)*x/ah+xi(11)*y/ahe+z*ueh(1))
      hi(1)=hi1
      ht(1)=ht(1)+hi1
      ii=4
      DO 115 l=2,4
      dhi=dxt1(l)*x/ah+xi(l+10)*y/ahe+z*ueh(l)
      tstl=l==4
      IF(tstl) dhi=dhi+dxt(1)
      dhi=averg*dhi
      ht(l)=ht(l)+dhi
      hi(l)=dhi
      IF(nosd) GOTO 115
      DO 114 m=l,4
      ii=ii+1
      dhi=dxt1(ii)*x/ah+xi(10+ii)*y/ahe+z*ueh(ii)
      IF(tstl) dhi=dhi+dxt(m)
      IF(m==4) dhi=dhi+dxt(l)
      dhi=averg*dhi
      ht(ii)=ht(ii)+dhi
  114 hi(ii)=dhi
  115 CONTINUE
c
c  *********************************************
c
c  pressure and enthalpy of heavy particles
c
  120 anh=(x/ah+y/ahe+z/az)*av
      rhh=rho(1)
      phh=anh*rhh*tkergs
      ph(1)=phh
      ph(2)=amm*phh*rho(2)
      drht=1+rho(3)
      ph(3)=amm*phh*drht
      ph(4)=amm*phh*rho(4)+rhh*tkergs*avd1
      DO 121 i=1,4
  121 pt(i)=pt(i)+ph(i)
      IF(nosd) GOTO 125
      ph(10)=amm*(phh*rho(10)+rho(4)*(ph(4)+rhh*tkergs*avd1))
      DO 122 k=1,3
      k1=k+1
      ph(k+4)=amm*(ph(k1)*rho(2)+phh*rho(k+4))
      IF(k>1) ph(k+6)=amm*(ph(k1)*drht+phh*rho(k+6))
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
     .  ph(l+1)*rho(k+4)+phh*rho(kl+6))
      IF(k>1) pt(9+kl)=pt(9+kl)+amm*(ph(kl)*drht+ph(k1)*
     .  rho(6+l)+ph(l+1)*rho(6+k)+phh*rho(9+kl))
  124 CONTINUE
      pt(20)=pt(20)+amm*(ph(10)*rho(4)+2*ph(4)*rho(10)
     .  +phh*rho(20)+rhh*tkergs*avd1*(amm*rho(4)*rho(4)+rho(10)))
c
  125 hhh=2.5*anh*tkergs
      hh(1)=hhh
      hh(2)=0
      hh(3)=amm*hhh
      hh(4)=2.5*tkergs*avd1
      hh(5)=0
      hh(6)=0
      hh(7)=0
      hh(8)=amm2*hhh
      hh(9)=amm*hh(4)
      hh(10)=0
      !
      ! BP 2000-11-14
      ! TRES dangeureux car, a priori, aucun sens valable en Fortran
      !
      hh(11:20) = 0.    ! CALL zeroc(hh(11),10)
      hh(17)=amm*hh(8)
      hh(18)=amm*hh(9)
      ht(1)=ht(1)+hhh
      ht(3)=ht(3)+amm*hhh
      hh4=2.5*tkergs*avd1
      ht(4)=ht(4)+hh4
      IF(nosd) GOTO 130
      ht(8)=ht(8)+amm2*hhh
      ht(9)=ht(9)+amm*hh4
      IF(notd) GOTO 130
      ht(17)=ht(17)+amm3*hhh
      ht(18)=ht(18)+amm2*hh4
c
c  *********************************************
c
c  pressure and enthalpy of radiation (included if iprrad /= 0)
c
  130 pr(:10) = 0.   ! CALL zeroc(pr,10)
      hr(:20) = 0.   ! CALL zeroc(hr,20)
      IF(iprrad==0) GOTO 145
      t4=t*t*t*t
      prr=car*t4/3
      pr(1)=prr
      !
      ! BP 2000-11-14
      ! TRES dangeureux car, a priori, aucun sens valable en Fortran
      !
      pr(2:10) = 0.   ! CALL zeroc(pr(2),9)
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
     .  -hrr*rho(jj))
  138 CONTINUE
  140 CONTINUE
      DO 142 i=1,ider
  142 ht(i)=ht(i)+hr(i)
c     IF(notd) GOTO 145
c..      WRITE(istdpr,14201) pe,ph,pr,(dph(i),i=1,10),pt
14201 FORMAT(/' pe:',1p10e12.4/' ph:',1p10e12.4/' pr:',10e12.4/
     *  ' dp:',10e12.4/' pt:',10e12.4/4x,10e12.4)
c
c  pressure and enthalpy from Coulomb effect
c
  145 IF(icoulm /= 1.and.icoulm /= -1) THEN
        DO 150 i=1,ider
        pt(i)=pt(i)+pcoul(i)
  150   ht(i)=ht(i)+hcoul(i)
      ENDIF
c
c  *********************************************
c
c  change to derivatives of log p
c
      ptt=pt(1)
c
c  first divide all derivatives by pt, to avoid over- and
c  underflow
c
      DO 15010 i=2,ider
15010 pt(i)=pt(i)/ptt
c
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
     .  +2*pt(m+1)*pt(l+1)*pt(n+1)))/amm
  155 ii=4
      DO 158 l=2,4
      ptl=pt(l)/amm
      IF(nosd) GOTO 158
      DO 156 m=l,4
      ii=ii+1
  156 pt(ii)=(pt(ii)-pt(l)*pt(m))/amm
  158 pt(l)=ptl
c     IF(.not.notd) WRITE(istdpr,15801) pt
15801 FORMAT(/' pt:',1p10e12.4/4x,10e12.4)
c
c  cp and dad
c
  160 pf=pt(2)
      hf=ht(2)
      dxtt=ht(3)*pf-hf*pt(3)
      lt=6
      DO 165 l=1,3
      lf=l+4
      IF(l>1) lt=6+l
  165 dxt(l+1)=ht(lt)*pf+ht(3)*pt(lf)-ht(lf)*pt(3)-hf*pt(lt)
c
      fcpp=1./(amm*t*pf)
      cpp=dxtt*fcpp
      cp(1)=cpp
      IF(nosd) GOTO 173
      DO 170 l=2,4
      dcp=dxt(l)*fcpp-cpp*pt(l+3)/pf
      IF(l==3) dcp=dcp-cpp*amm
  170 cp(l)=dcp
c
  173 prh=amm*ptt/rhh
      anum=pf*prh-hf
      dad(1)=anum/dxtt
      IF(nosd) GOTO 177
      DO 175 l=2,4
      lf=l+3
  175 dad(l)=(prh*(amm*(pt(l)-rho(l))*pf+pt(lf))-ht(lf)-dad(1)*dxt(l)
     .  )/dxtt
c  further derivatives
  177 rhf=rho(2)
      dxtt=pt(3)*rhf-pf*rho(3)
c
c  Gamma1 and derivatives
c
      dxt(1)=rhf-dad(1)*dxtt
      gm1=pf/dxt(1)
      gmm1(1)=gm1
      IF(secder) THEN
        dxt(2)=rho(5)-dad(2)*dxtt
     *    -dad(1)*(pt(6)*rho(2)+pt(3)*rho(5)-pt(5)*rho(3)-pt(2)*rho(6))
        dxt(3)=rho(6)-dad(3)*dxtt
     *    -dad(1)*(pt(8)*rho(2)+pt(3)*rho(6)-pt(6)*rho(3)-pt(2)*rho(8))
        dxt(4)=rho(7)-dad(4)*dxtt
     *    -dad(1)*(pt(9)*rho(2)+pt(3)*rho(7)-pt(7)*rho(3)-pt(2)*rho(9))
        gmm1(2)=(pt(5)-gm1*dxt(2))/dxt(1)
        gmm1(3)=(pt(6)-gm1*dxt(3))/dxt(1)
        gmm1(4)=(pt(7)-gm1*dxt(4))/dxt(1)
      ELSE
         !
         ! BP 2000-11-14
         ! TRES dangeureux car, a priori, aucun sens valable en Fortran
         !
         gmm1(2:4) = 0.   ! CALL zeroc(gmm1(2),3)
      ENDIF
c
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
      IF(l>2) lt=l+5
  180 dlt(l)=(pt(lt)*rhf+pt(3)*rho(lf)-pt(lf)*rho(3)-pf*rho(lt))/pf
     .  -delta*pt(lf)/pf
  190 xii1(1)=xii(1)
      xii1(2)=xii(11)
      xii1(3)=xii(21)
      xii1(4)=anuh(1)/anh0
      RETURN
 1010 FORMAT(/' **** error in eqstf. With H-, Ne  <=  0.'/
     *         '      x(H+) =',1pe11.3,'  x(H-) =',e11.3,
     *         ' x(heavies) =',e11.3/
     *         '      Ignore H-')
 1100 FORMAT(/' ***** error in s/r eqstf.',
     *  ' Iteration for Coulomb term failed to converge.'/
     *  7x,' log f, log T, X, Z =',4f10.5//
     *  ' last changes in delta mu_c(1-3):'/ 1p3e13.5)
      END SUBROUTINE eqstfc
c
c***************************************************************************
c
	SUBROUTINE eqstpc(pl, tl, xh, yh, zh, nosd, notd, fl, nit)
c
c	Modification de eqstrh par Annie Baglin Juillet 91
c
c  equation of state routine, with indepENDent variables
c    pl = log(p)
c    tl  = log(T)
c
c  Iterates, using eqstfc, to determine log(f), and hence set equation
c  of state variables.
c  log(f) is RETURNed in fl.
c  nit RETURNs number of iterations.
c
c  Accuracy is set in eps below, and should be adjusted for
c  computer.
c
c  If fl  <=  -1.e10 on input, a trial value is set in the routine.
c  Otherwise the value of fl provided is used as trial.
c
c  As for s/r eqstfc, s/r setcnsc must be CALLed before CALLing
c  eqstrh.
c
c  Note: to obtain the trial value of fl, a fit to mue**-1, as a
c  fonction of log T, is used, in the region of partial ionization.
c  this is currently based on solar composition and
c  log T - log rho relation, and hence may give problems under
c  different circumstances.
c
c  Original version: 07/08/87
c
c  This version is set up for use with Coulomb term version
c  of eqstfc. When starting from trial value of fl, make initial
c  iteration for fl without Coulomb terms.

c  Note that present this is the case regardless of input fl

c  Version from 1/5/90.

c  Double precision version.
c  ++++++++++++++++++++++++

c  Dated 9/3/90

c	MODIF : appel a EFF en cas de non convergence P. MOREL oct 1991
c	ajout du COMMON /marche/ avec ok=.true. si c'est bon

      IMPLICIT DOUBLE PRECISION(a-h,o-z)
	IMPLICIT INTEGER(i-n)
      LOGICAL nosd, notd
      DIMENSION flini(2), wt(2)
      COMMON/eqstdc/ xii1(4), ane(10), rho(20), ht(20), p(20), cp(4),
     *  dad(4), dlt(4), gm1, tprh, trhp, rhxp,gmm1(4)
c
c  COMMON defining standard input and output
c
      COMMON/cstdio/ istdou, istdpr
      COMMON/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      COMMON/ccoulm/ epsdmu, icoulm, iclmit, iclmsd, nitdmu, epssdr

	LOGICAL ok
	COMMON /marche/ok 
c
      SAVE
c
      DATA eps /1.d-10/
c
c  store original value of Coulomb case flag
c
      icoulp=icoulm
c
c  test for setting trial
c
c
c  try to use wd setting of trial fl in all cases
c
      icase = 1
      icoulm=-1
c	WRITE(*,*) tl,pl
      CALL inteffc(tl,pl,rhl,flini,wt,icase,iextr)
      fl = wt(1)*flini(1) + wt(2)*flini(2)
c      WRITE(*,*) 'tl, rhl, pl =', tl, rhl,pl
c      WRITE(*,*) 'trial fl set to', fl
      IF(fl <= -1.e10) THEN
       xt=max(0.d0,tl-3.76119)
        xt=xt*xt
        xmufl=3.8*xt/(xt+0.035)-3.83015
        fl=7.829+xmufl+rhl-1.5*tl
c
      ENDIF
c
c  start loop for iteration
c
      nit=0
   10 CALL eqstfc(fl, tl, xh, yh, zh, nosd, notd)
c
c
      IF(p(1) <= 0.d0 .OR. p(2) == 0.d0) THEN
        WRITE(istdpr, 105) fl, tl, pl, p(1)
  	ok=.FALSE.		!modif appel a gong2
        icoulm=icoulp           !yveline 040100
	IF(.true.)RETURN			!modif appel a gong2
        CALL dmpeqs
        STOP
      ENDIF
      pli=log10(p(1))
      dfl=(pl - pli)/p(2)
      nit=nit+1
c
c  limit size of dfl
c
      IF(ABS(dfl) > 0.4d0) dfl=sign(0.4d0,dfl)
c      IF(idgeos>=1) WRITE(istdpr,*) nit,pl, tl, fl, dfl
c
c  test for continuing iteration
c
      IF(ABS(dfl)<eps) THEN
        GOTO 20
      ELSE IF(nit <= 60) THEN
        fl = fl+dfl
        GOTO 10
      ELSE
c
c  diagnostics for failed convergence
c
        WRITE(istdou,110) pl, tl, fl, dfl	!probleme
	ok=.FALSE.
        IF(istdpr /= istdou)THEN		!modif: appel a EFF
	 WRITE(istdpr,110) pl, tl, fl, dfl
	 ok=.FALSE.
	ENDIF
        icoulm=icoulp      !Yveline 040100
	IF(.true.)RETURN
        fl = -1.e11
        GOTO 20
c
      ENDIF
c
c  this is the END
c
   20 CONTINUE
c
c  test for repeat with Coulomb term
c
      IF(icoulm /= icoulp) THEN
        icoulm=icoulp
        nit=0
        GOTO 10
      ENDIF
c..      WRITE(*,*) 'converged fl =',fl
      RETURN
  105 FORMAT(//' **** error in eqstpc. p  <=  0 or p(2) = 0 for'/
     *  ' log f =',f10.5,'  log T =',f10.5,' log p =', f10.5,
     *  ' pi =',1pe13.5,' appel a ETAT_EFF')
  110 FORMAT(//' ***** Iteration failed in s/r eqstrh.'/
     *  7x,'log(p) =',f10.5,'  log(T) =',f10.5,
     *  '  Last log(f) =',f10.5,'  dfl =',1pe13.5,' appel a ETAT_EFF')
      END SUBROUTINE eqstpc
C
C***************************************************************************
C
	SUBROUTINE eqstrh(rhl, tl, xh, yh, zh, nosd, notd, fl, nit)
c
c  equation of state routine, with indepENDent variables
c    rhl = log(rho)
c    tl  = log(T)
c
c  Iterates, using eqstfc, to determine log(f), and hence set equation
c  of state variables.
c  log(f) is RETURNed in fl.
c  nit RETURNs number of iterations.
c
c  Accuracy is set in eps below, and should be adjusted for
c  computer.
c
c  If fl  <=  -1.e10 on input, a trial value is set in the routine.
c  Otherwise the value of fl provided is used as trial.
c
c  As for s/r eqstfc, s/r setcnsc must be CALLed before CALLing
c  eqstrh.
c
c  Note: to obtain the trial value of fl, a fit to mue**-1, as a
c  fonction of log T, is used, in the region of partial ionization.
c  this is currently based on solar composition and
c  log T - log rho relation, and hence may give problems under
c  different circumstances.

c  Original version: 07/08/87

c  This version is set up for use with Coulomb term version
c  of eqstf. When starting from trial value of fl, make initial
c  iteration for fl without Coulomb terms.

c  Note that present this is the case regardless of input fl

c  Version from 1/5/90.

c  Double precision version.
c  ++++++++++++++++++++++++

c  Dated 9/3/90

      IMPLICIT DOUBLE PRECISION(a-h,o-z)
	IMPLICIT INTEGER(i-n)
      LOGICAL :: nosd, notd
      DIMENSION flini(2), wt(2)
      COMMON/eqstdc/ xii(4), ane(10), rho(20), ht(20), p(20), cp(4),
     *  dad(4), dlt(4), gm1, tprh, trhp, rhxp,gmm1(4)
c
c  COMMON defining standard input and output
c
      COMMON/cstdio/ istdou, istdpr
      COMMON/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      COMMON/ccoulm/ epsdmu, icoulm, iclmit, iclmsd, nitdmu, epssdr
c
      SAVE
c
      DATA eps /1.d-10/
c
c  store original value of Coulomb case flag
c
      icoulp=icoulm
c
c  test for setting trial
c
c
c  try to use wd setting of trial fl in all cases
c
      icase = 2
      icoulm=-1
      CALL inteffc(tl,pgl,rhl,flini,wt,icase,iextr)
      fl = wt(1)*flini(1) + wt(2)*flini(2)
c..      WRITE(*,*) 'tl, rhl =', tl, rhl
c..      WRITE(*,*) 'trial fl set to', fl
      IF(fl <= -1.e10) THEN
c..        xt=max(0.d0,tl-3.76119)
c..        xt=xt*xt
c..        xmufl=3.8*xt/(xt+0.035)-3.83015
c..        fl=7.829+xmufl+rhl-1.5*tl
c
      ENDIF
c
c  start loop for iteration
c
      nit=0
   10 CALL eqstfc(fl, tl, xh, yh, zh, nosd, notd)
c
      IF(rho(1) <= 0.d0 .OR. rho(2)  == 0.d0) THEN
        WRITE(istdpr, 105) fl, tl, rhl, rho(1)
        CALL dmpeqs
        STOP
      ENDIF
      rhli=log10(rho(1))
      dfl=(rhl - rhli)/rho(2)
      nit=nit+1
c
c  limit size of dfl
c
      IF(ABS(dfl)>0.4) dfl=sign(0.4d0,dfl)
      IF(idgeos>=1) WRITE(istdpr,*) nit,rhl, tl, fl, dfl
c
c  test for continuing iteration
c
      IF(ABS(dfl)<eps) THEN
        GOTO 20
      ELSE IF(nit <= 60) THEN
        fl = fl+dfl
        GOTO 10
      ELSE
c
c  diagnostics for failed convergence
c
        WRITE(istdou,110) rhl, tl, fl, dfl
        IF(istdpr /= istdou)
     *    WRITE(istdpr,110) rhl, tl, fl, dfl
        fl = -1.e11
        GOTO 20
c
      ENDIF
c
c  this is the END
c
   20 CONTINUE
c
c  test for repeat with Couilomb term
c
      IF(icoulm /= icoulp) THEN
        icoulm=icoulp
        nit=0
        GOTO 10
      ENDIF
c..      WRITE(*,*) 'converged fl =',fl
      RETURN
  105 FORMAT(//' **** error in eqstrh. rho  <=  0 or rho(2) = 0 for'/
     *  ' log f =',f10.5,'  log T =',f10.5,' log rho =', f10.5,
     *  ' rhoi =',1pe13.5)
  110 FORMAT(//' ***** Iteration failed in s/r eqstrh.'/
     *  7x,'log(rho) =',f10.5,'  log(T) =',f10.5,
     *  '  Last log(f) =',f10.5,'  dfl =',1pe13.5)
      END SUBROUTINE eqstrh
C
c***************************************************************************
c
	SUBROUTINE eosder(fl, tl, xh, yh, zh, nosd, notd)
c
c  routine to test EOS derivatives
c
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
	IMPLICIT INTEGER(i-n)
      PARAMETER(npar = 3)
      CHARACTER*6 cvar(3), cvar2(3,4)
      LOGICAL nosd, notd
	DIMENSION plp(10,2), hp(10,2), plc(10), hc(10),
	1	dmucc(3,4),dmucp(3,4,2),xiic(30),xiip(30,2),pcoulc(10),
	2	pcoulp(10,2), hcoulc(10), hcoulp(10,2),ddmuc(3),dddmuc(3),
	3	rhoc(10), rhop(10,2),dxii(3),ddxii(3),
	4	df4c(npar),df4p(npar,2),d2f4ic(npar,npar),d2f4xc(npar),
	5	d2f4rc(npar),d2f4tc(npar),dp4ic(npar),
	6	dh4ic(npar),kjstor(3,4),dxii2(3,4),ddxii2(3,4),
	7	gmm1c(4), gmm1p(4,2)
c
      COMMON/eqscntc/ anh0,anhe0,ihvz,iprrad,ihmin
      COMMON/ln10c/ amm,amm2,amm3
      COMMON/constsc/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
	COMMON/dmudec/ dmu,dmuf,dmut,dmux,dmuff,dmuft,dmufx,dmutt,dmutx,
	1	dmuxx,idmu
      COMMON/dmucdr/ dmuc(npar,10)
      COMMON/df4der/ df4(npar), d2f4i(npar,npar), d2f4f(npar,npar),
     .          d2f4x(npar),d2f4r(npar),d2f4t(npar),
     .          dp4i(npar),dh4i(npar), p4, f4, dp4dr, dp4dt, dp4dx,
     .          dh4dr, dh4dt, dh4dx
      COMMON/eqsoutc/ east(30),xii(30),dne(10),dph(20),he(20),pe(10),
     *  hi(20),pi(10),hh(20),ph(10),hr(20),pr(10),pcoul(10),hcoul(10)
      COMMON/eqstdc/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
      DATA cvar /'log f','log T','X'/
      DATA kjstor /2,3,4,5,6,7,6,8,9,7,9,10/
      DATA cvar2 /'f','T','X','ff','fT','fX','Tf','TT','TX',
     *  'Xf','XT','XX'/
c
c  PRINT contribution of central values to X-derivatives
c
c..      WRITE(*,*) ' contr. to pX:', amm*pt(4), dph(4)/pt(1),
c..     *  pe(4)/pt(1), pi(4)/pt(1), ph(4)/pt(1),
c..     *  pr(4)/pt(1), pcoul(4)/pt(1)
c..      WRITE(*,*) ' contr. to HX:', ht(4)/ht(1), dph(14)/ht(1),
c..     *  he(4)/ht(1), hi(4)/ht(1), hh(4)/ht(1),
c..     *  hr(4)/ht(1), hcoul(4)/ht(1)
c
c  store central values of cluster
c
      CALL eqstfc(fl, tl, xh, yh, zh, nosd, notd)
c
      flc=fl
      tlc=tl
      xhc=xh
      yhc=1-xh-zh
c
      CALL storec(df4,df4c,npar)
      CALL storec(d2f4i,d2f4ic,npar*npar)
      CALL storec(d2f4x,d2f4xc,npar)
      CALL storec(d2f4t,d2f4tc,npar)
      CALL storec(d2f4r,d2f4rc,npar)
      CALL storec(dp4i,dp4ic,npar)
      CALL storec(dh4i,dh4ic,npar)
      dh4drc=dh4dr
      dh4dtc=dh4dt
      dh4dxc=dh4dx
c
      plc(1)=log10(pt(1))
      hc(1)=ht(1)
      rhoc(1)=log10(rho(1))
      CALL storec(pt(2),plc(2),9)
      CALL storec(ht(2),hc(2),9)
      CALL storec(rho(2),rhoc(2),9)
      CALL storec(xii,xiic,30)
      CALL storec(dmuc,dmucc,12)
      CALL storec(pcoul,pcoulc,10)
      CALL storec(hcoul,hcoulc,10)
      CALL storec(gmm1,gmm1c,4)
c
      step=0.00001
      WRITE(*,100) fl, tl, xh, step
c
      DO 80 k=1,3
      WRITE(*,*) ' '
      fl=flc
      tl=tlc
      xh=xhc
      dtl=0
      dxh=0
      stp=step
      DO 35 j=1,2
      IF(k==1) THEN
        fl=flc+stp
      ELSE IF(k==2) THEN
        tl=tlc+stp
        dtl=2*step
      ELSE
        xh=xhc+stp
        yh=yhc-stp
        dxh=2*step
      ENDIF
      CALL eqstfc(fl, tl, xh, yh, zh, nosd, notd)
c
      CALL storec(df4,df4p(1,j),npar)
c
      plp(1,j)=log10(pt(1))
      rhop(1,j)=log10(rho(1))
      CALL storec(pt(2),plp(2,j),9)
      CALL storec(rho(2),rhop(2,j),9)
      CALL storec(ht,hp(1,j),10)
      CALL storec(xii,xiip(1,j),30)
      CALL storec(dmuc,dmucp(1,1,j),12)
      CALL storec(pcoul,pcoulp(1,j),10)
      CALL storec(hcoul,hcoulp(1,j),10)
      CALL storec(gmm1,gmm1p(1,j),4)
   35 stp=-stp
c
c  set and PRINT derivatives
c
      dp_=(plp(1,1)-plp(1,2))/(2*step)
      ddp=dp_/plc(k+1)-1
      drho=(rhop(1,1)-rhop(1,2))/(2*step)
      ddrho=drho/rhoc(k+1)-1
      dh=(hp(1,1)-hp(1,2))/(2*step)
      ddh=dh/hc(k+1)-1
      DO 36 i=1,3
      ddmuc(i)=(dmucp(i,1,1)-dmucp(i,1,2))/(2*step)
      dddmuc(i)=(ddmuc(i)-dmucc(i,k+1))/(1.e-20+ABS(dmucc(i,k+1)))
      ix=1+10*(i-1)
      dxii(i)=(xiip(ix,1)-xiip(ix,2))/(2*step)
      ddxii(i)=(dxii(i)-xiic(ix+k))/(1.e-20+ABS(xiic(ix+k)))
   36 CONTINUE
c
      dpcoul=(pcoulp(1,1)-pcoulp(1,2))/(2*step)
      ddpcou=(dpcoul-pcoulc(k+1))/(ABS(pcoulc(k+1))+1.e-10)
      dhcoul=(hcoulp(1,1)-hcoulp(1,2))/(2*step)
      ddhcou=(dhcoul-hcoulc(k+1))/(ABS(hcoulc(k+1))+1.e-10)
      WRITE(*,120) 'log p',cvar(k),dp_,plc(k+1), ddp
      WRITE(*,120) 'log rho',cvar(k),drho,rhoc(k+1), ddrho
      WRITE(*,120) 'H',cvar(k),dh,hc(k+1), ddh
      WRITE(*,120) 'dmuc(1)',cvar(k),ddmuc(1),dmucc(1,k+1),dddmuc(1)
      WRITE(*,120) 'dmuc(2)',cvar(k),ddmuc(2),dmucc(2,k+1),dddmuc(2)
      WRITE(*,120) 'dmuc(3)',cvar(k),ddmuc(3),dmucc(3,k+1),dddmuc(3)
      WRITE(*,120) 'xii(1)',cvar(k),dxii(1),xiic(1+k),ddxii(1)
      WRITE(*,120) 'xii(2)',cvar(k),dxii(2),xiic(11+k),ddxii(2)
      WRITE(*,120) 'xii(3)',cvar(k),dxii(3),xiic(21+k),ddxii(3)
      WRITE(*,120) 'p Coul',cvar(k),dpcoul,pcoulc(k+1), ddpcou
      WRITE(*,120) 'H Coul',cvar(k),dhcoul,hcoulc(k+1), ddhcou
c
c  second derivatives
c
      DO 40 j=2,4
      kj=kjstor(k,j)
      DO 38 i=1,3
      ix0=10*(i-1)
      ixj=ix0+j
      dxii2(i,j)=(xiip(ixj,1)-xiip(ixj,2))/(2*step)
   38 ddxii2(i,j)=(dxii2(i,j)-xiic(ix0+kj))/(1.e-20+ABS(xiic(ix0+kj)))
c
      WRITE(*,120) 'xii(1)',cvar2(k,j),dxii2(1,j),xiic(kj),ddxii2(1,j)
      WRITE(*,120) 'xii(2)',cvar2(k,j),dxii2(2,j),xiic(10+kj),
     *  ddxii2(2,j)
      WRITE(*,120) 'xii(3)',cvar2(k,j),dxii2(3,j),xiic(20+kj),
     *  ddxii2(3,j)
c
      dpcl2=(pcoulp(j,1)-pcoulp(j,2))/(2*step)
      ddpcl2=(dpcl2-pcoulc(kj))/(ABS(pcoulc(kj))+1.e-10)
      dhcl2=(hcoulp(j,1)-hcoulp(j,2))/(2*step)
      ddhcl2=(dhcl2-hcoulc(kj))/(ABS(hcoulc(kj))+1.e-10)
c
      WRITE(*,120) 'p Coul',cvar2(k,j),dpcl2,pcoulc(kj), ddpcl2
      WRITE(*,120) 'H Coul',cvar2(k,j),dhcl2,hcoulc(kj), ddhcl2
c
   40 CONTINUE
c
c  Gamma1
c
      dgm1=(gmm1p(1,1)-gmm1p(1,2))/(2*step)
      ddgm1=(dgm1-gmm1c(k+1))/(ABS(gmm1c(k+1))+1.e-10)
      WRITE(*,120) 'Gamma1',cvar(k),dgm1,gmm1c(k+1), ddgm1
 
c
c  set change in df4's
c
      WRITE(*,*) dtl, rhop(1,1) - rhop(1,2), dxh
c..     *  ,(xiip(1+(i-1)*10,1) - xiip(1+(i-1)*10,2),i=1,3)
c..      DO 78 i=1,3
c..      sum=0
c..      DO 77 l=1,3
c..   77 sum=sum+d2f4ic(i,l)*(xiip(1+(l-1)*10,1) - xiip(1+(l-1)*10,2))
c..      ddf4=sum+d2f4rc(i)*(rhop(1,1)-rhop(1,2))+d2f4tc(i)*dtl
c..     *  +d2f4xc(i)*dxh
c..   78 WRITE(*,*) i,df4p(i,1)-df4p(i,2),ddf4,sum,
c..     *  d2f4rc(i)*(rhop(1,1)-rhop(1,2)),d2f4tc(i)*dtl,d2f4xc(i)*dxh
c..c
      sump=0
      sumh=0
      DO 79 l=1,3
      sumh=sumh+dh4ic(l)*(xiip(1+(l-1)*10,1) - xiip(1+(l-1)*10,2))
   79 sump=sump+dp4ic(l)*(xiip(1+(l-1)*10,1) - xiip(1+(l-1)*10,2))
c..      dp4=sump+dp4drc*(rhop(1,1)-rhop(1,2))+dp4dtc*dtl
c..     *  +dp4dxc*dxh
      dh4=sumh+dh4drc*(rhop(1,1)-rhop(1,2))+dh4dtc*dtl
     *  +dh4dxc*dxh
c..      WRITE(*,*) ' p4', pcoulp(1,1)-pcoulp(1,2),dp4,sump,
c..     *  dp4drc*(rhop(1,1)-rhop(1,2)),dp4dtc*dtl,
c..     *  dp4dxc*dxh
c..      WRITE(*,*) 'In eosder: dh4drc, dh4dtc =',dh4drc, dh4dtc
      WRITE(*,79901) ' H4', hcoulp(1,1)-hcoulp(1,2),dh4,sumh,
     *  dh4drc*(rhop(1,1)-rhop(1,2)),dh4dtc*dtl,
     *  dh4dxc*dxh
79901 FORMAT(a4,1p6e13.5)
c
   80 CONTINUE
c
c  make sure that original values are RETURNed
c
      fl=flc
      tl=tlc
      xh=xhc
      CALL eqstfc(fl, tl, xh, yh, zh, nosd, notd)
c
      RETURN
  100 FORMAT(/' test derivatives at log f, log T, X =',3f10.5/
     *  ' step =',1pe11.3//
     *  ' derivative, numerical, analytical, rel. error:'/)
  120 FORMAT(a10,' wrt ',a10,':',1p3e13.5)
      END SUBROUTINE eosder
c
c********************************************************************************
c
	SUBROUTINE f4der(rhol,tl,sn,nspe,f4,e4,p4,h4,d2f4r2,d2f4t2,
     .           d2f4rt,df4,d2f4i,d2f4f,d2f4x,d2f4r,d2f4t,
     .           dp4i,dp4dr,dp4dt,dp4dx,dh4i,dh4dr,dh4dt,dh4dx,npar)
c
c     nspe  : number of particles actually used
c     npar  : first DIMENSION of d2f4i,d2f4f (=number of ionization
c             degrees = number of number fractions)
c
c     f4    : Coulomb configurational free energy (Debye-Huckel)
c     e4    : internal energy
c     p4    : pressure
c     h4    : enthalpy
c
c     df4   : derivatives with respect to reaction PARAMETERs ("Saha equations")
c     d2f4i : derivatives of df4 with respect to ionization degrees
c     d2f4f : derivatives of df4 with respect to number fractions
c     d2f4x : derivatives of df4 with respect to X (Y varying accordingly,
c             and at fixed ionization degrees)
c
c     d2f4r : derivatives of df4 with respect to log10 rho
c     d2f4t : derivatives of df4 with respect to log10 t
c     d2f4r2: derivative  of  f4 with respect to log10 rho ** 2
c     d2f4t2: derivative  of  f4 with respect to log10 t   ** 2
c     d2f4rt: derivative  of  f4 with respect to log10 rho and log10 t
c
c     in both d2f4i and d2f4f the element (i,j) denotes the derivative
c     of df4(i) with respect to PARAMETER j (ionization degree or number
c     fraction)
c
c     dp4i  : derivatives of p4 with respect to ionization degrees
c     dp4dr : derivative  of p4 with respect to log10 rho
c     dp4dt : derivative  of p4 with respect to log10 t
c     dp4dx : derivative  of p4 with respect to X (Y varying accordingly
c             and at fixed ionization degrees)

c     dh4i  : derivatives of h4 with respect to ionization degrees
c     dh4dr : derivative  of h4 with respect to log10 rho
c     dh4dt : derivative  of h4 with respect to log10 t
c     dh4dx : derivative  of h4 with respect to X (Y varying accordingly
c             and at fixed ionization degrees)

c=========================================================================

      IMPLICIT REAL*8 (a-h,o-z)
	IMPLICIT INTEGER(i-n)
c
      INTEGER, PARAMETER :: mfe=4, mspes=18
c
      COMMON /freecl/
     .                         e     (mfe   ),        p     (mfe   ),
     .                         f     (mfe   ),               d2fdt2 ,
     .                                d2fdtv ,               d2fdv2 ,
     .                         fscr  (mspes ),        dfdn  (mspes ),
     .                         d2fdnt(mspes ),        d2fdnv(mspes ),
     .                         d2fdn2( mspes ,mspes )
      COMMON/constsc/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
c
      DIMENSION sn(nspe),df4(npar),d2f4i(npar,npar),d2f4f(npar,npar),
     .          d2f4x(npar),d2f4r(npar),d2f4t(npar),
     .          dp4i(npar),dh4i(npar)
c
      DATA umod /2.302585092994046 d0/
c
c...assume ah = 1.008, ahe = 4.0026
      DATA hovhe,heovh /0.251836d0,3.97083d0/
c
c  reset with JC-D values
c
      hovhe = ah/ahe
      heovh = ahe/ah
c
      t   = 10.d0**tl
      vol = 10.d0**(-rhol)
c
      CALL f4n(t,vol,sn,nspe)
c
c     f4    : Coulomb configurational free energy (Debye-Huckel)
c     e4    : internal energy
c     p4    : pressure
c     p4    : enthalpy
c
      f4     =  f(4)
      e4     =  e(4)
      p4     =  p(4)
      h4     =  e(4) + p(4)*vol
c
c     d2f4r2: derivative of  f4 with respect to log10 rho ** 2
c     d2f4t2: derivative of  f4 with respect to log10 t   ** 2
c     d2f4rt: derivative of  f4 with respect to log10 rho and log10 t
c
      d2f4r2 =  d2fdv2*vol*vol*umod*umod
      d2f4rt =  d2fdtv*  t*vol*umod*umod
      d2f4t2 =  d2fdt2*  t*  t*umod*umod
c
c     df4   : derivatives with respect to reaction PARAMETERs ("Saha equations")
c
      df4(1) =  dfdn(2) + dfdn(6)
      df4(2) =  dfdn(4) + dfdn(6)
      df4(3) = -dfdn(4) + dfdn(5) + dfdn(6)
c
c     d2f4r : derivatives of df4 with respect to log10 rho
c
      d2f4r(1) = -(             d2fdnv(2) + d2fdnv(6))*vol*umod
      d2f4r(2) = -(             d2fdnv(4) + d2fdnv(6))*vol*umod
      d2f4r(3) = -(-d2fdnv(4) + d2fdnv(5) + d2fdnv(6))*vol*umod
c
c     d2f4t : derivatives of df4 with respect to log10 t
c
      d2f4t(1) =  (             d2fdnt(2) + d2fdnt(6))*t  *umod
      d2f4t(2) =  (             d2fdnt(4) + d2fdnt(6))*t  *umod
      d2f4t(3) =  (-d2fdnt(4) + d2fdnt(5) + d2fdnt(6))*t  *umod
c
      toth   =  sn(1) + sn(2)
      tothe  =  sn(3) + sn(4) + sn(5)
c
c     d2f4i : derivatives of df4 with respect to ionization degrees
c
      d2f4i(1,1) = toth *(-d2fdn2(2,1)+d2fdn2(2,2)+d2fdn2(2,6)
     .                    -d2fdn2(6,1)+d2fdn2(6,2)+d2fdn2(6,6))
      d2f4i(1,2) = tothe*(-d2fdn2(2,3)+d2fdn2(2,4)+d2fdn2(2,6)
     .                    -d2fdn2(6,3)+d2fdn2(6,4)+d2fdn2(6,6))
      d2f4i(1,3) = tothe*(-d2fdn2(2,3)+d2fdn2(2,5)+d2fdn2(2,6)*2.d0
     .                    -d2fdn2(6,3)+d2fdn2(6,5)+d2fdn2(6,6)*2.d0)
c
      d2f4i(2,1) = toth *(-d2fdn2(4,1)+d2fdn2(4,2)+d2fdn2(4,6)
     .                    -d2fdn2(6,1)+d2fdn2(6,2)+d2fdn2(6,6))
      d2f4i(2,2) = tothe*(-d2fdn2(4,3)+d2fdn2(4,4)+d2fdn2(4,6)
     .                    -d2fdn2(6,3)+d2fdn2(6,4)+d2fdn2(6,6))
      d2f4i(2,3) = tothe*(-d2fdn2(4,3)+d2fdn2(4,5)+d2fdn2(4,6)*2.d0
     .                    -d2fdn2(6,3)+d2fdn2(6,5)+d2fdn2(6,6)*2.d0)
c
      d2f4i(3,1) = toth *( d2fdn2(4,1)-d2fdn2(4,2)-d2fdn2(4,6)
     .                    -d2fdn2(5,1)+d2fdn2(5,2)+d2fdn2(5,6)
     .                    -d2fdn2(6,1)+d2fdn2(6,2)+d2fdn2(6,6))
      d2f4i(3,2) = tothe*( d2fdn2(4,3)-d2fdn2(4,4)-d2fdn2(4,6)
     .                    -d2fdn2(5,3)+d2fdn2(5,4)+d2fdn2(5,6)
     .                    -d2fdn2(6,3)+d2fdn2(6,4)+d2fdn2(6,6))
      d2f4i(3,3) = tothe*( d2fdn2(4,3)-d2fdn2(4,5)-d2fdn2(4,6)*2.d0
     .                    -d2fdn2(5,3)+d2fdn2(5,5)+d2fdn2(5,6)*2.d0
     .                    -d2fdn2(6,3)+d2fdn2(6,5)+d2fdn2(6,6)*2.d0)
c
      fach       =  sn(1)*sn(1)/toth
c
c  Note: temporary fudge to avoid problems when sn(3) or sn(4)
c  are zero. Need correction from WD           10/5/90
c
      IF(sn(3) /= 0.d0) THEN
        ff2      =  sn(4)/sn(3)
      ELSE
        ff2      =  1
      ENDIF
      IF(sn(4) /= 0.d0) THEN
        ff3      =  sn(5)/sn(4)
      ELSE
        ff3      =  1
      ENDIF
      ff22       =  ff2*ff2
      fac2       =  1.d0 + ff2
      fac3       =  1.d0 + ff3
      fache      =  sn(3)*sn(3)/tothe
      dhed2      = -fac3*fache
      dhed3      = -ff2*fache
      dhepd2     =  fache
      dhepd3     = -ff22*fache
      dhe2d2     =  ff3*fache
      dhe2d3     =  ff2*fac2*fache
      ded2       =  (1.d0 + ff3*2.d0)*fache
      ded3        = ff2*(2.d0 + ff2)*fache
c
c     d2f4f : derivatives of df4 with respect to number fractions
c
      d2f4f(1,1) = fach*(-d2fdn2(2,1)+d2fdn2(2,2)+d2fdn2(2,6)
     .                   -d2fdn2(6,1)+d2fdn2(6,2)+d2fdn2(6,6))
      d2f4f(2,1) = fach*(-d2fdn2(4,1)+d2fdn2(4,2)+d2fdn2(4,6)
     .                   -d2fdn2(6,1)+d2fdn2(6,2)+d2fdn2(6,6))
      d2f4f(3,1) = fach*( d2fdn2(4,1)-d2fdn2(4,2)-d2fdn2(4,6)
     .                   -d2fdn2(5,1)+d2fdn2(5,2)+d2fdn2(5,6)
     .                   -d2fdn2(6,1)+d2fdn2(6,2)+d2fdn2(6,6))
c
      d2f4f(1,2) =    ( d2fdn2(2,3)+d2fdn2(6,3))*dhed2
     .              + ( d2fdn2(2,4)+d2fdn2(6,4))*dhepd2
     .              + ( d2fdn2(2,5)+d2fdn2(6,5))*dhe2d2
     .              + ( d2fdn2(2,6)+d2fdn2(6,6))*ded2
      d2f4f(2,2) =    ( d2fdn2(4,3)+d2fdn2(6,3))*dhed2
     .              + ( d2fdn2(4,4)+d2fdn2(6,4))*dhepd2
     .              + ( d2fdn2(4,5)+d2fdn2(6,5))*dhe2d2
     .              + ( d2fdn2(4,6)+d2fdn2(6,6))*ded2
      d2f4f(3,2) =    (-d2fdn2(4,3)+d2fdn2(5,3)+d2fdn2(6,3))*dhed2
     .              + (-d2fdn2(4,4)+d2fdn2(5,4)+d2fdn2(6,4))*dhepd2
     .              + (-d2fdn2(4,5)+d2fdn2(5,5)+d2fdn2(6,5))*dhe2d2
     .              + (-d2fdn2(4,6)+d2fdn2(5,6)+d2fdn2(6,6))*ded2
c
      d2f4f(1,3) =    ( d2fdn2(2,3)+d2fdn2(6,3))*dhed3
     .              + ( d2fdn2(2,4)+d2fdn2(6,4))*dhepd3
     .              + ( d2fdn2(2,5)+d2fdn2(6,5))*dhe2d3
     .              + ( d2fdn2(2,6)+d2fdn2(6,6))*ded3
      d2f4f(2,3) =    ( d2fdn2(4,3)+d2fdn2(6,3))*dhed3
     .              + ( d2fdn2(4,4)+d2fdn2(6,4))*dhepd3
     .              + ( d2fdn2(4,5)+d2fdn2(6,5))*dhe2d3
     .              + ( d2fdn2(4,6)+d2fdn2(6,6))*ded3
      d2f4f(3,3) =    (-d2fdn2(4,3)+d2fdn2(5,3)+d2fdn2(6,3))*dhed3
     .              + (-d2fdn2(4,4)+d2fdn2(5,4)+d2fdn2(6,4))*dhepd3
     .              + (-d2fdn2(4,5)+d2fdn2(5,5)+d2fdn2(6,5))*dhe2d3
     .              + (-d2fdn2(4,6)+d2fdn2(5,6)+d2fdn2(6,6))*ded3
c
c     d2f4x : derivatives of df4 with respect to X (Y varying accordingly)
c
      dnhdx    =  (toth + tothe*heovh)/toth
      dnhedx   = -(toth*hovhe + tothe)/tothe
c
      d2f4x(1) = dnhdx *((d2fdn2(2,1)+d2fdn2(6,1))* sn(1) +
     .                   (d2fdn2(2,2)+d2fdn2(6,2))* sn(2) +
     .                   (d2fdn2(2,6)+d2fdn2(6,6))* sn(2))+
     .           dnhedx*((d2fdn2(2,3)+d2fdn2(6,3))* sn(3) +
     .                   (d2fdn2(2,4)+d2fdn2(6,4))* sn(4) +
     .                   (d2fdn2(2,5)+d2fdn2(6,5))* sn(5) +
     .                   (d2fdn2(2,6)+d2fdn2(6,6))*(sn(4)+2.d0*sn(5)))
c
      d2f4x(2) = dnhdx *((d2fdn2(4,1)+d2fdn2(6,1))* sn(1) +
     .                   (d2fdn2(4,2)+d2fdn2(6,2))* sn(2) +
     .                   (d2fdn2(4,6)+d2fdn2(6,6))* sn(2))+
     .           dnhedx*((d2fdn2(4,3)+d2fdn2(6,3))* sn(3) +
     .                   (d2fdn2(4,4)+d2fdn2(6,4))* sn(4) +
     .                   (d2fdn2(4,5)+d2fdn2(6,5))* sn(5) +
     .                   (d2fdn2(4,6)+d2fdn2(6,6))*(sn(4)+2.d0*sn(5)))
c
      d2f4x(3) = dnhdx *((-d2fdn2(4,1)+d2fdn2(5,1)+d2fdn2(6,1))*
     .                     sn(1) +
     .                   (-d2fdn2(4,2)+d2fdn2(5,2)+d2fdn2(6,2))*
     .                     sn(2) +
     .                   (-d2fdn2(4,6)+d2fdn2(5,6)+d2fdn2(6,6))*
     .                     sn(2))+
     .           dnhedx*((-d2fdn2(4,3)+d2fdn2(5,3)+d2fdn2(6,3))*
     .                     sn(3) +
     .                   (-d2fdn2(4,4)+d2fdn2(5,4)+d2fdn2(6,4))*
     .                     sn(4) +
     .                   (-d2fdn2(4,5)+d2fdn2(5,5)+d2fdn2(6,5))*
     .                     sn(5) +
     .                   (-d2fdn2(4,6)+d2fdn2(5,6)+d2fdn2(6,6))*
     .                    (sn(4)+2.d0*sn(5)))
c
c     dp4i  : derivatives of p4 with respect to ionization degrees
c     dp4dr : derivative  of p4 with respect to log10 rho
c     dp4dt : derivative  of p4 with respect to log10 t
c     dp4dx : derivative  of p4 with respect to X (Y varying accordingly
c             and at fixed ionization degrees)
c
      dp4i(1) = - toth *( -d2fdnv(1)+d2fdnv(2)+d2fdnv(6))
      dp4i(2) = - tothe*( -d2fdnv(3)+d2fdnv(4)+d2fdnv(6))
      dp4i(3) = - tothe*( -d2fdnv(3)+d2fdnv(5)+d2fdnv(6)*2.d0)
c
      dp4dr   =    umod  * vol * d2fdv2
      dp4dt   =  - umod  *   t * d2fdtv
c
      dp4dx   = -dnhdx * (d2fdnv(1)* sn(1) +
     .                    d2fdnv(2)* sn(2) +
     .                    d2fdnv(6)* sn(2))
     .          -dnhedx* (d2fdnv(3)* sn(3) +
     .                    d2fdnv(4)* sn(4) +
     .                    d2fdnv(5)* sn(5) +
     .                    d2fdnv(6)*(sn(4)+2.d0*sn(5)))
c
c     ... analogous for h4 ...
c
      dh4i(1) =   toth *((-dfdn  (1)+dfdn  (2)+dfdn  (6))
     .                  -(-d2fdnt(1)+d2fdnt(2)+d2fdnt(6))*t
     .                  -(-d2fdnv(1)+d2fdnv(2)+d2fdnv(6))*vol  )
c
      dh4i(2) =   tothe*((-dfdn  (3)+dfdn  (4)+dfdn  (6))
     .                  -(-d2fdnt(3)+d2fdnt(4)+d2fdnt(6))*t
     .                  -(-d2fdnv(3)+d2fdnv(4)+d2fdnv(6))*vol  )
c
      dh4i(3) =   tothe*((-dfdn  (3)+dfdn  (5)+dfdn  (6)*2.d0)
     .                  -(-d2fdnt(3)+d2fdnt(5)+d2fdnt(6)*2.d0)*t
     .                  -(-d2fdnv(3)+d2fdnv(5)+d2fdnv(6)*2.d0)*vol  )
c
      dh4dr   =    umod * vol * t * d2fdtv  +  vol * dp4dr
      dh4dt   =  - umod *  t  * t * d2fdt2  +  vol * dp4dt
c
      dh4dx   = dnhdx *((dfdn(1)-t*d2fdnt(1))* sn(1) +
     .                  (dfdn(2)-t*d2fdnt(2))* sn(2) +
     .                  (dfdn(6)-t*d2fdnt(6))* sn(2))+
     .          dnhedx*((dfdn(3)-t*d2fdnt(3))* sn(3) +
     .                  (dfdn(4)-t*d2fdnt(4))* sn(4) +
     .                  (dfdn(5)-t*d2fdnt(5))* sn(5) +
     .                  (dfdn(6)-t*d2fdnt(6))*(sn(4) + 2.d0*sn(5)))+
     .          vol * dp4dx
c
      END SUBROUTINE f4der
c
c*****************************************************************************
c
	SUBROUTINE f4n(t,vol,snn,nspe)
c
c     derivatives with respect to number abundances.
c     CALLs f4 of MHD package and prepares quantities not provided by f4
c     (degeneracy-related stuff).
c
c === all units c.g.s. and degrees Kelvin

      IMPLICIT REAL*8 (a-h,o-z)
	IMPLICIT INTEGER(i-n)

      INTEGER, PARAMETER :: mspes=18, mfe=4, mfd=5

      DIMENSION snn(nspe)

      COMMON /cstmhd/
     .                                camu   ,               cc     ,
     .                                ce     ,               ch     ,
     .                                ck     ,               cme    ,
     .                                cpi    ,               cevw   ,
     .                                cx     ,               cf4    ,
     .                                carad

      COMMON /degmhd/
     .                                eta    ,               detdn  ,
     .                                detdt  ,               detdv  ,
     .                                d2etdn2,               d2etdnt,
     .                                d2etdnv,               d2etdt2,
     .                                d2etdtv,               d2etdv2,
     .                                exeta  ,               thet   ,
     .                                dthet  ,               d2thet ,
     .                                fd    (mfd   )
      COMMON /freecl/
     .                         e     (mfe   ),        p     (mfe   ),
     .                         f     (mfe   ),               d2fdt2 ,
     .                                d2fdtv ,               d2fdv2 ,
     .                         fscr  (mspes ),        dfdn  (mspes ),
     .                         d2fdnt(mspes ),        d2fdnv(mspes ),
     .                         d2fdn2( mspes ,mspes )
      COMMON /mapcl /
     .                                nchem  ,               nspes  ,
     .                                ise
      COMMON /spescl/
     .                         zs    (mspes ),        zsq   (mspes ),
     .                         zmask (mspes ),        sn    (mspes )

	LOGICAL ok
	COMMON /marche/ok

c
      DO 5 is=1,nspes
 5    sn(is) = snn(is)

	CALL neweta(t,vol,sn(ise),ier)
	IF(ier /= 0) THEN
         WRITE(*,*) 'error in s/r f4n. failure in neweta. f4 not CALLed'
	 WRITE(*,*)'appel a ETAT_EFF'
	 ok=.FALSE.	!appel a etat_eff
	 RETURN
	ENDIF
c
      CALL f4mhd(t,vol)
c
      END SUBROUTINE f4n
c
c*****************************************************************************
c
	SUBROUTINE f4mhd(t,vol)
c     ................................
c     (modified from MHD package)
c     ................................

	USE mod_kind

      IMPLICIT REAL*8 (a-h,o-z)
	IMPLICIT INTEGER(i-n)
c
c     free energy of coulomb interactions and derivatives

      INTEGER, PARAMETER :: mspes=18, mfe=4, mfd=5

c  PARAMETERs for tau expansion. xtautr gives transition between
c  expansion and direct expression.
c  With nmax = 15, xtautr = 0.1, the error at the transition point
c  is around 1.e-13 in tau, 1.e-10 in dtau and 1.e-9 in d2tau.
c
      INTEGER, PARAMETER :: nmax=15
      REAL(kind=dp), PARAMETER :: xtautr=0.1d0
c
c-------------------------------------------------------------------
c
      COMMON /cstmhd/
     .                                camu   ,               cc     ,
     .                                ce     ,               ch     ,
     .                                ck     ,               cme    ,
     .                                cpi    ,               cevw   ,
     .                                cx     ,               cf4    ,
     .                                carad
c
      COMMON /degmhd/
     .                                eta    ,               detdn  ,
     .                                detdt  ,               detdv  ,
     .                                d2etdn2,               d2etdnt,
     .                                d2etdnv,               d2etdt2,
     .                                d2etdtv,               d2etdv2,
     .                                exeta  ,               thet   ,
     .                                dthet  ,               d2thet ,
     .                                fd    (mfd   )
      COMMON /freecl/
     .                         e     (mfe   ),        p     (mfe   ),
     .                         f     (mfe   ),               d2fdt2 ,
     .                                d2fdtv ,               d2fdv2 ,
     .                         fscr  (mspes ),        dfdn  (mspes ),
     .                         d2fdnt(mspes ),        d2fdnv(mspes ),
     .                         d2fdn2( mspes ,mspes )
      COMMON /mapcl /
     .                                nchem  ,               nspes  ,
     .                                ise
      COMMON /spescl/
     .                         zs    (mspes ),        zsq   (mspes ),
     .                         zmask (mspes ),        sn    (mspes )
c
c
c  for running on Suns, take out quadruple precision.
c
c  *** Note: we may need to fix this up in expansion later, by
c      including more terms
c
c..      REAL*16 dpx, dp1, dp2, dp3
c
      DIMENSION   dxdn(mspes), d2xdnt(mspes), d2xdnv(mspes),
     .            d2xdn2(mspes, mspes)
      DIMENSION   ctau(nmax), cdtau(nmax), cd2tau(nmax)
c
      EQUIVALENCE (dxdn  , dzdn  ), (d2xdnt, d2zdnt), (d2xdnv, d2zdnv),
     .            (d2xdn2, d2zdn2)
c
      DATA initcf /0/
c
      SAVE initcf, ctau, cdtau, cd2tau
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     sums over charges                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      sumn = dmax1( sdot(nspes - 1, zmask, 1, sn, 1), 1.d-70 )
      zn   = dmax1( sdot(nspes - 1, zs   , 1, sn, 1), 1.d-70 )
      znt  = dmax1( sdot(nspes - 1, zsq  , 1, sn, 1) + sn(ise)*thet,
     .                                              1.d-70 )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     PARAMETER x and its derivatives                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      frat = fd(2) / fd(3)
      x    = cx * frat * (zn/sumn) * sqrt( znt ) / sqrt( vol*t**3 )
c
c-----------------------------------------------------------------------
c     zero everything
c-----------------------------------------------------------------------
c
      DO 2 is = 1, nspes
      fscr  (is) = 0.0d0
      dxdn  (is) = 0.0d0
      d2xdnt(is) = 0.0d0
      d2xdnv(is) = 0.0d0
      DO 1 js = 1, nspes
      d2xdn2(is, js) = 0.0d0
    1 CONTINUE
    2 CONTINUE
c
c-----------------------------------------------------------------------
c     first derivatives                                                c
c-----------------------------------------------------------------------
c
      DO 3 is = 1, nspes - 1
      dxdn(is) = x * ( zs(is)/zn - zmask(is)/sumn+0.5d0*zsq(is)/znt )
    3 CONTINUE
c
      dxdn(ise) = x * ( 0.5d0* (thet + sn(ise)*dthet*detdn)/znt
     .                - 1.5d0* frat * detdn
     .                + 1.0d0 / sn(ise) )
c
      dxdt = x * ( detdt*(thet - 1.5d0*frat + 0.5d0*sn(ise)*dthet/znt)
     .           - 1.5d0/t   )
      dxdv = x * ( detdv*(thet - 1.5d0*frat + 0.5d0*sn(ise)*dthet/znt)
     .           - 0.5d0/vol )
c
c-----------------------------------------------------------------------
c     second derivatives
c-----------------------------------------------------------------------
c
      DO 5 js = 1, nspes - 1
c
      IF(zs(js) /= 0.d0) THEN
            DO 4 is = 1, js
            d2xdn2(is, js) = ( dxdn(is) * dxdn(js) / x
     .                   + x * (       zmask(is) * zmask(js) / sumn**2
     .                         -       zs   (is) * zs   (js) / zn**2
     .                       - 0.5d0 * zsq  (is) * zsq  (js) / znt**2) )
    4       CONTINUE
c
            IF( js > 1 ) CALL scopy( js - 1, d2xdn2( 1,js), 1,
     .                                          d2xdn2(js, 1), mspes )
      d2xdnt(js) = dxdt * dxdn(js) / x
     .               - 0.5d0 * x * zsq(js)*sn(ise) * dthet*detdt/znt**2
      d2xdnv(js) = dxdv * dxdn(js) / x
     .               - 0.5d0 * x * zsq(js)*sn(ise) * dthet*detdv/znt**2
      ENDIF
    5 CONTINUE
c
      DO 6 is = 1, nspes - 1
      d2xdn2(is, ise) = ( dxdn(is) * dxdn(ise) / x
     .   - 0.5d0 * x * zsq(is) * (thet + sn(ise)*dthet*detdn)/znt**2 )
    6 CONTINUE
      CALL scopy( nspes - 1, d2xdn2(1,ise), 1, d2xdn2(ise,1), mspes )
c
      d2xdn2(ise, ise) =
     .    dxdn(ise)**2/x
     .    - x * ( 1.0d0/sn(ise)**2
     .       + 1.5d0 * frat * (d2etdn2 + (thet - 1.5d0*frat)*detdn**2)
     .         - 0.5d0 * ( d2thet * sn(ise) * detdn**2
     .                  + dthet * (2.d0*detdn + sn(ise)*d2etdn2)
     .                - (thet + sn(ise)*dthet*detdn)**2/znt ) / znt  )
c
      d2xdnt(ise) = dxdt*dxdn(ise)/x + x * (
     .  -1.5d0*frat* (d2etdnt + (thet - 1.5d0*frat)*detdn*detdt)
     .  +0.5d0*( d2thet * sn(ise)*detdn*detdt
     .       + dthet * (detdt + sn(ise)*d2etdnt)
     .     - (thet+sn(ise)*dthet*detdn)*sn(ise)*dthet*detdt/znt)/znt )
c
      d2xdnv(ise) = dxdv*dxdn(ise)/x + x * (
     .  -1.5d0*frat* (d2etdnv + (thet - 1.5d0*frat)*detdn*detdv)
     .  +0.5d0*( d2thet * sn(ise)*detdn*detdv
     .       + dthet * (detdv + sn(ise)*d2etdnv)
     .     - (thet+sn(ise)*dthet*detdn)*sn(ise)*dthet*detdv/znt)/znt )
c
      d2xdt2 = dxdt**2/x + x * ( 1.5d0/t**2
     .     + d2etdt2 * ( thet - 1.5d0*frat + 0.5d0*sn(ise)*dthet/znt )
     .       + detdt**2 * ( dthet - 1.5d0*frat*(thet - 1.5d0 * frat)
     .         + 0.5d0*sn(ise)*(d2thet - sn(ise)*dthet**2/znt)/znt ) )
c
      d2xdtv = dxdt*dxdv/x + x * (
     .       d2etdtv * ( thet - 1.5d0*frat + 0.5d0*sn(ise)*dthet/znt )
     .    + detdt * detdv * ( dthet - 1.5d0*frat*(thet - 1.5d0 * frat)
     .         + 0.5d0*sn(ise)*(d2thet - sn(ise)*dthet**2/znt)/znt ) )
c
      d2xdv2 = dxdv**2/x + x * ( 0.5d0/vol**2
     .     + d2etdv2 * ( thet - 1.5d0*frat + 0.5d0*sn(ise)*dthet/znt )
     .     + detdv**2 * ( dthet - 1.5d0*frat*(thet - 1.5d0 * frat)
     .         + 0.5d0*sn(ise)*(d2thet - sn(ise)*dthet**2/znt)/znt ) )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tau and its derivatives                                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c..      IF( x > 1.0d-1 ) THEN
c
c
      IF( x > xtautr ) THEN
          dpx =   x
          dp1 =   3.d0 * ( log(1.d0 + dpx)
     .                   - dpx*(1.d0 - 0.5d0*dpx) ) / dpx**3
          dp2 = - 3.d0 * ( dp1 - 1.d0/(1.d0 + dpx)    ) / dpx
          dp3 =   - ( 4.d0*dp2 + 3.d0/(1.d0 + dpx)**2 ) / dpx
          tau   = dp1
          dtau  = dp2
          d2tau = dp3
      ELSE
c..          tau  = ((((((((x/11. - 1.d0/10.)*x +  1.d0/9.)*x - 1.d0/8.)*x
c..     .                         + 1.d0/7. )*x -  1.d0/6.)*x + 1.d0/5.)*x
c..     .                         - 1.d0/4. )*x +  1.d0/3.)*3.d0
c..          dtau = (((((((8.d0*x/11. - 7.d0/10.)*x +2.d0/3.)*x-5.d0/8.)*x
c..     .                              + 4.d0/7.)*x -1.d0/2.)*x+2.d0/5.)*x
c..     .                              - 1.d0/4. )*3.d0
c..          d2tau= ((((((56.d0*x/11. - 21.d0/5.)*x+10.d0/3.)*x-5.d0/2.)*x
c..     .                             + 12.d0/7.)*x-1.d0)*x+2.d0/5.)*3.d0
c
c  test for setting coefficients
c
        IF(initcf==0) THEN
c
c  set coefficients for expansion
c
          isg=6*mod(nmax,2)-3
          DO n=nmax,3,-1
            ctau(n)  =dfloat(isg)/dfloat(n)
            cdtau(n) =dfloat(isg*(n-3))/dfloat(n)
            cd2tau(n)=dfloat(isg*(n-3)*(n-4))/dfloat(n)
            isg=-isg
          ENDDO
          initcf = 1
        ENDIF
c
c  DO expansion as DO loop, to allow arbitrary high order
c
        tau=ctau(nmax)
        DO n=nmax-1,3,-1
          tau=tau*x+ctau(n)
        ENDDO
c
        dtau=cdtau(nmax)
        DO n=nmax-1,4,-1
          dtau=dtau*x+cdtau(n)
        ENDDO
c
        d2tau=cd2tau(nmax)
        DO n=nmax-1,5,-1
          d2tau=d2tau*x+cd2tau(n)
        ENDDO
      ENDIF
c
c..      WRITE(*,*) 'x, tau, etc =',x,tau, dtau, d2tau
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     free energy, pressure, internal energy                           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c..      WRITE(*,*) 'x, tau, dtau',x, tau,dtau
c     free energy
      f(4) = - cf4 * tau * sqrt( znt**3 / (t * vol) )
c
c     pressure
      p(4) = f(4) * ( 0.5d0/vol-1.5d0* sn(ise) * dthet * detdv / znt
     .                          - dxdv * dtau / tau )
c
c     internal energy
      e(4) = f(4) * ( 1.5d0/t -1.5d0* sn(ise) * dthet * detdt / znt
     .                          - dxdt * dtau / tau ) * t
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     derivatives of f4                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c-----------------------------------------------------------------------
c     first derivatives                                                c
c-----------------------------------------------------------------------
c
      df4dt =   f(4) * ( -0.5d0/ t + 1.5d0*sn(ise) * dthet * detdt / znt
     .                            + dxdt * dtau / tau )
      df4dv = - p(4)
c
      DO 7 is = 1, nspes - 1
      fscr(is) = f(4)*(1.5d0*zsq(is)/znt + dxdn(is)*dtau/tau)
      dfdn(is) = fscr(is)
    7 CONTINUE
      fscr(ise) = f(4)*( 1.5d0* ( thet + sn(ise)*dthet*detdn )/znt
     .                   + dxdn(ise) * dtau/tau )
      dfdn(ise) = fscr(ise)
c
c-----------------------------------------------------------------------
c     second derivatives                                               c
c-----------------------------------------------------------------------
c
      DO 9 js = 1, nspes - 1
c
      IF(zs(js)/= 0.d0) THEN
            DO 8 is = 1, js
            d2fdn2(is, js) = ( fscr(is) * fscr(js) / f(4)
     .            + ( (d2tau/tau - (dtau/tau)**2) * dxdn(is) * dxdn(js)
     .            + d2xdn2(is, js) * dtau / tau
     .            - 1.5d0* zsq(is) * zsq(js) / znt**2 ) * f(4) )
    8       CONTINUE
            IF( js > 1 ) CALL scopy( js - 1, d2fdn2(1 ,js), 1,
     .                                          d2fdn2(js, 1), mspes )
c
            d2fdnt(js) = fscr(js) * df4dt / f(4)
     .              + ( (d2tau/tau - (dtau/tau)**2) * dxdn(js) * dxdt
     .              + d2xdnt(js) * dtau/tau
     .              - 1.5d0*zsq(js)*sn(ise)*dthet*detdt/znt**2)*f(4)
c
            d2fdnv(js) = fscr(js) * df4dv / f(4)
     .              + ( (d2tau/tau - (dtau/tau)**2) * dxdn(js) * dxdv
     .              + d2xdnv(js) * dtau/tau
     .              - 1.5d0*zsq(js)*sn(ise)*dthet*detdv/znt**2)*f(4)
      ENDIF
    9 CONTINUE
c
      DO 10 is = 1, nspes - 1
      d2fdn2(is, ise) = ( fscr(is) * fscr(ise) / f(4)
     .    + f(4)*( (d2tau/tau - (dtau/tau)**2) * dxdn(is) * dxdn(ise)
     .    + d2xdn2(is, ise) * dtau / tau
     .    -1.5d0*zsq(is)*(thet+sn(ise)*dthet*detdn)/znt**2))
   10 CONTINUE
      CALL scopy( nspes - 1, d2fdn2(1,ise), 1, d2fdn2(ise,1), mspes )
c
      d2fdn2(ise, ise) = fscr(ise)**2 / f(4)
     .    + f(4) * ( (d2tau/tau - (dtau/tau)**2) * dxdn(ise)**2
     .    + d2xdn2(ise, ise) * dtau / tau
     .    + 1.5d0*( d2thet * sn(ise) * detdn**2
     .    + dthet * ( 2.d0*detdn + sn(ise) * d2etdn2 )
     .    - (thet + sn(ise)*dthet*detdn)**2/znt) / znt )
c
      d2fdnt(ise) = fscr(ise)*df4dt/f(4) + f(4) *
     . ( (d2tau/tau - (dtau/tau)**2)*dxdt*dxdn(ise)
     . + d2xdnt(ise) * dtau/tau
     . + 1.5d0 *( d2thet*sn(ise)*detdt*detdn
     . +dthet*(detdt + sn(ise)*d2etdnt)
     . -(thet+sn(ise)*dthet*detdn)*sn(ise)*dthet*detdt/znt)/znt )
c
      d2fdnv(ise) = fscr(ise)*df4dv/f(4) + f(4) *
     . ( (d2tau/tau - (dtau/tau)**2)*dxdv*dxdn(ise)
     . + d2xdnv(ise) * dtau/tau
     . + 1.5d0 *( d2thet*sn(ise)*detdv*detdn
     . +dthet*(detdv + sn(ise)*d2etdnv)
     . -(thet+sn(ise)*dthet*detdn)*sn(ise)*dthet*detdv/znt)/znt )
c
      d2fdtv = df4dt*df4dv/f(4) + f(4) *
     . ( (d2tau/tau - (dtau/tau)**2)*dxdt*dxdv
     . + d2xdtv * dtau/tau
     . + 1.5d0 * sn(ise) * ( d2thet*detdt*detdv
     . + dthet*(d2etdtv-sn(ise)*dthet*detdt*detdv/znt))/znt )
c
      d2fdt2 = df4dt**2/f(4) + f(4) *
     . ( 0.5d0 / t**2
     . + (d2tau/tau - (dtau/tau)**2)*dxdt**2
     . + d2xdt2 * dtau/tau
     . + 1.5d0 * sn(ise) * ( d2thet*detdt**2
     . + dthet*(d2etdt2-sn(ise)*dthet*detdt**2/znt))/znt )
c
      d2fdv2 = df4dv**2/f(4) + f(4) *
     . ( 0.5d0 / vol**2
     . + (d2tau/tau - (dtau/tau)**2)*dxdv**2
     . + d2xdv2 * dtau/tau
     . + 1.5d0 * sn(ise) * ( d2thet*detdv**2
     . + dthet*(d2etdv2-sn(ise)*dthet*detdv**2/znt))/znt )
c
      END SUBROUTINE f4mhd
c
c******************************************************************************
c
	SUBROUTINE hheion(eah, eahe, eahep, xih, xihe, xihep, secder)
c
c  given ratios between succesive states of ionization, and
c  derivatives, for H, He and He+ in eah(1-10), eahe(1-10),
c  and eahep(1-10), sets the corresponding degrees of ionization,
c  and derivatives, into xih(1-10), xihe(1-10) and xihep(1-10)
c
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
	IMPLICIT INTEGER(i-n)
      LOGICAL secder
      DIMENSION eah(10), eahe(10), eahep(10), xih(10), xihe(10),
     *  xihep(10), eahs(10), eahes(10), eahepc(10)
c
      IF(secder) THEN
        iders=10
      ELSE
        iders=4
      ENDIF
c
c  set combined ratio for he+
c
      eea=eahep(1)
      IF(eea>0) THEN
        eea1=eahe(1)
        eahepc(1)=eea*eea1
        ii=4
        IF(secder) THEN
          DO 15 l=2,4
          DO 15 m=l,4
          ii=ii+1
   15     eahepc(ii)=eea*eahe(ii)+eahep(l)*eahe(m)+eahe(l)*eahep(m)+
     .      eea1*eahep(ii)
        ENDIF
        DO 20 i=2,4
   20   eahepc(i)=eahep(i)*eea1+eea*eahe(i)
c
      ELSE
c
        DO 22 i=1,iders
   22   eahepc(i)=0
c
      ENDIF
c
c  set x h+, x he+, x he++ and derivatives
c
      dnm=1+eah(1)
      xih(1)=eah(1)/dnm
c
c  hydrogen fraction
c
c  to avoid over- and underflow, divide derivatives
c  of ea by dnm (modification 3/1/84)
c
      DO 25 i=2,iders
   25 eahs(i)=eah(i)/dnm
c
      ii=4
      DO 30 l=1,3
      l1=l+1
      eal=eahs(l1)
      xih(l1)=eal/dnm
      IF(secder) THEN
        DO 28 m=l,3
        m1=m+1
        ii=ii+1
   28   xih(ii)=(eahs(ii)-2*eal*eahs(m1))/dnm
      ENDIF
   30 CONTINUE
c
c  helium fractions
c
      dnm=1+eahe(1)+eahepc(1)
c
c  to avoid over- and underflow, divide derivatives
c  of ea by dnm (modification 3/1/84)
c
      DO 35 i=2,iders
      eahes(i)=eahe(i)/dnm
   35 eahepc(i)=eahepc(i)/dnm
c
      ii=4
      eeahe=eahe(1)
      eeahep=eahepc(1)
      xihe(1)=eeahe/dnm
      xihep(1)=eeahep/dnm
      DO 40 l=2,4
      ealhe=eahes(l)
      ealhep=eahepc(l)
      anmhe=(1+eeahep)*ealhe-eeahe*ealhep
      anmhep=(1+eeahe)*ealhep-eeahep*ealhe
      xihe(l)=anmhe/dnm
      xihep(l)=anmhep/dnm
c
c  second derivatives
c
      IF(secder) THEN
        DO 38 m=l,4
        ii=ii+1
        eamhe=eahes(m)
        eamhep=eahepc(m)
c
c  for xi < 1.e-10 zero second derivatives
c
        IF(xihe(1) <= 1.e-10) THEN
          xihe(ii)=0
        ELSE
          eamhe=eahes(m)
          eamhep=eahepc(m)
          xihe(ii)=ealhe*eamhep-ealhep*eamhe+
     *      ((1+eeahep)*eahes(ii)-eeahe*eahepc(ii)
     *      -2*anmhe*(eamhe+eamhep))/dnm
        ENDIF
c
        IF(xihep(1) <= 1.e-10) THEN
          xihep(ii)=0
        ELSE
          xihep(ii)=ealhep*eamhe-ealhe*eamhep+
     *      ((1+eeahe)*eahepc(ii)-eeahep*eahes(ii)
     *      -2*anmhep*(eamhep+eamhe))/dnm
        ENDIF
   38   CONTINUE
      ENDIF
   40 CONTINUE
c
      END SUBROUTINE hheion
c
c*******************************************************************************
c
	SUBROUTINE hmnion(tl, eah, ehm, xhm, secder)
c
c  given ratio between succesive states of ionization, and
c  derivatives, for H in eah(1-10)
c  sets fraction of H-, and derivatives, into xhm(1-10)
c
c  Note: assumes that fraction of H- is always very small.
c
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
	IMPLICIT INTEGER(i-n)
      LOGICAL secder
      DIMENSION eah(10), ehm(10),xhm(10), eahs(10)
      COMMON/constsc/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
      COMMON/ln10c/ amm,amm2,amm3
	COMMON/dmudec/ dmu,dmuf,dmut,dmux,dmuff,dmuft,dmufx,dmutt,dmutx,
	1	dmuxx,idmu
      COMMON/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
c     IF(idgeos>=3) WRITE(*,*) 'Entering hmnion'
      t=10.d0**tl
      tkev=ck1*t
c
      xhm(:10) = 0.    ! CALL zeroc(xhm,10)
c
      ext=exhm/tkev
      eea=ext-dmu
c
c  test for no h-
c
      IF(eea <= -100) THEN
        IF(idgeos>=3) WRITE(*,*) 'No H-'
        RETURN
      ENDIF
      eea=0.5*exp(eea)
      dnm=1+eah(1)
      xhm(1)=eea/dnm
c
c  to avoid over- and underflow, divide derivatives
c  of ea by dnm (modification 3/1/84)
c
      IF(secder) THEN
        iders=10
      ELSE
       iders=4
      ENDIF
c
      DO 10 i=2,iders
   10 eahs(i)=eah(i)/dnm
c
      ehm(1)=eea
      ehm(2)=-dmuf*eea
      ehm(3)=-(amm*ext+dmut)*eea
      ehm(4)=-dmux*eea
      IF(secder) THEN
        ehm(5)=-dmuff*eea-dmuf*ehm(2)
        ehm(6)=-dmuft*eea-dmuf*ehm(3)
        ehm(7)=-dmufx*eea-dmux*ehm(2)
        ehm(8)=(amm2*ext-dmutt)*eea-(amm*ext+dmut)*ehm(3)
        ehm(9)=-dmutx*eea-dmux*ehm(3)
        ehm(10)=-dmuxx*eea-dmux*ehm(4)
      ENDIF
c
c  derivatives of xh-
c
      ii=4
      DO 20 l=1,3
      l1=l+1
      ehml=ehm(l1)
      eal=eahs(l1)
      xhm(l1)=(ehml-eea*eal)/dnm
      IF(secder) THEN
        DO 15 m=l,3
        m1=m+1
        ii=ii+1
   15   xhm(ii)=(ehm(ii)-(ehml*eahs(m1)+ehm(m1)*eal+eea*eahs(ii))
     *    +2*eea*eal*eahs(m1))/dnm
      ENDIF
   20 CONTINUE
c      IF(idgeos>=3) WRITE(*,*) 'x(H-) =',xhm(1)
      END SUBROUTINE hmnion
c
c*******************************************************************************
c
	SUBROUTINE hvionac(fl,tl,x,y,z,nosd,notd,anu,ue,anur,uer)
c
c  Calculate ionization of heavy elements.
c
c  Modification 17/5/90: include possibility of including full ionization
c     of heavy elements, for ihvz = 4. Note that this involves changing
c     size of arrays in COMMONs /hvredu/ and /xhv/
c
c  Modification 4/6/90: include average degree of ionization for
c     each element in array xhvmn(10) at the beginning of
c     COMMON/hviond/
c
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
	IMPLICIT INTEGER(i-n)
      LOGICAL nosd,notd
      CHARACTER name*4
      DIMENSION anu(*),ue(*),anur(*),uer(*)
      DIMENSION eir(10),dr(10),hr(10),gr(10),xi(29),phi(30),hst(30)
      COMMON/constsc/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
      COMMON/potetcc/ chi(125),am(10),iz(10)
      COMMON/hvname/ name(10)
      COMMON/hvabndc/ ab(10),iab
      COMMON /hvredu/ irab,jfirst,ir(10),chir(125),izr(10),abr(10),
     *  amr(10)
      COMMON /hvcntl/ icount,iWRITE,dptst0,dptst1
      COMMON/hviond/ xhvmn(10),xhv(125)
      COMMON/eqscntc/ anz0,anze0,ihvz,iprrad,ihmin
	COMMON/dmudec/ dmup(10),idmu
      COMMON/eqdpco/ frhi,bdcoh,bdcoz,idpco
      COMMON/ln10c/ amm,amm2,amm3
c
c  COMMON defining standard input and output
c
      COMMON/cstdio/ istdou, istdpr
      COMMON/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
      SAVE
c
c  test for restricting set of elements
c
      IF(icount>0) GOTO 10
      icount=1
      IF(iWRITE==1) THEN
        WRITE(istdpr,200)
        WRITE(istdpr,205)
        lj=0
        DO 3 k=1,10
        izj=iz(k)
        DO 2 i=1,izj
        lj=lj+1
    2   xi(i)=chi(lj)
        im=min0(izj,13)
    3   WRITE(istdpr,210) name(k),ab(k),izj,(xi(i),i=1,im)
      ENDIF
c
c  test for restricting set of elements or levels included.
c  for ihvz = 4 keep everything, but for consistency store
c  full set of inFORMATion in restricted arrays
c
      IF(ihvz==1) THEN
c
c  reset to restricted set of variables
c
        irab=3
        ir(1)=1
        ir(2)=3
        ir(3)=10
c
c  reset abundances, to get reasonable fit to full case
c
        abr(1)=1.02*ab(1)
        abr(2)=1.35*ab(3)
        abr(3)=3.71*ab(10)
c
c  element with lowest ionization potential (here fe)
c
        jfirst=3
c
c  number of elements treated fully
c
        jiz1=2
c
      ELSE
c
c  use full set of elements
c
        irab=10
        DO 76 i=1,iab
        abr(i)=ab(i)
   76   ir(i)=i
c
c  lowest potential is of na
c
        jfirst=5
c
c  number of elements treated fully
c
        IF(ihvz==2) THEN
          jiz1=0
        ELSE IF(ihvz==3) THEN
          jiz1=3
        ELSE IF(ihvz==4) THEN
          jiz1=10
        ELSE
          WRITE(istdpr,110) ihvz
          IF(istdou /= istdpr) WRITE(istdou,110) ihvz
          STOP
        ENDIF
      ENDIF
c
c  set possibly reduced set of ionization DATA
c
   79 j=1
      lj0=0
      ljr=0
      DO 5 i=1,10
c
      izi=iz(i)
      IF(i==ir(j)) THEN
c
c  test for inclusion of all levels
c
        IF(j <= jiz1) THEN
          izrj=izi
        ELSE
          izrj=1
        ENDIF
        izr(j)=izrj
        amr(j)=am(i)
c
c  restore ionization potentials
c
        lj=lj0
        DO 4 k=1,izrj
        lj=lj+1
        ljr=ljr+1
    4   chir(ljr)=chi(lj)
        j=j+1
      ENDIF
c
    5 lj0=lj0+izi
c
c  reset anz0 and anze0
c
      anz0=0
      anze0=0
      DO 77 j=1,irab
      anz0=anz0+abr(j)*izr(j)/amr(j)
   77 anze0=anze0+abr(j)*izr(j)
c
      IF(iWRITE==1) THEN
c
        WRITE(istdpr,220)
        WRITE(istdpr,205)
        lj=0
        DO 7 k=1,irab
        izrj=izr(k)
        DO 6 i=1,izrj
        lj=lj+1
    6   xi(i)=chir(lj)
        im=min0(izrj,13)
    7   WRITE(istdpr,210) name(ir(k)),abr(k),izrj,(xi(i),i=1,im)
        WRITE(istdpr,230) anz0
c
      ENDIF
c
c  change to summed ionization potentials
c
      lj=0
      DO 9 k=1,irab
      sum=0
      izj=izr(k)
      DO 9 i=1,izj
      lj=lj+1
      sum=sum+chir(lj)
    9 chir(lj)=sum

c  set total number of levels for (possibly) reduced set


c  ***************************************

c  END of initialization section

   10 f=1.d1**fl
      t=1.d1**tl
c  test for departure of heavy elements
      IF(idpco < 2) bdcoz=1
c  k*t, in ev
      tk=ck1*t
c
c  test whether phderc has already been CALLed
c
      IF(idmu==1) GOTO 15
c
      CALL phderc(fl,tl,phi,hst,nosd,notd)
c
c
c  psi and derivatives
      wf= sqrt(1+f)
      psi=2*wf+log(f/((1+wf)*(1+wf)))
      psif=amm*wf
      psff=amm2*f/2/wf
c
      zf=x+2*y+anze0*z
      zf2=zf*zf
      zf3=zf*zf2
c
      ak0=ckh*zf2
c  a0**3*ne/(k*t)
      aa=caa*phi(1)/(tk*zf3)
c
c  delta mu and derivatives
c
      bmu=tk+20*ak0
      dmu=aa*bmu
      dmps=dmu-psi
      dmup(1)=dmps
      ref=phi(2)
      ret=phi(3)
      dmuf=dmu*ref*amm
      dmut=aa*(bmu*ret-20*ak0)*amm
      dmux=aa*(3*tk+20*ak0)/zf
      dmup(2)=dmuf-psif
      dmup(3)=dmut
      dmup(4)=dmux
c
      IF(nosd) GOTO 18
      dmup(5)=(dmuf*ref+dmu*phi(4)*amm)*amm-psff
      dmup(6)=(dmut*ref+dmu*phi(5)*amm)*amm
      dmup(7)=dmux*ref*amm
      dmup(8)=aa*(20*ak0*(1-2*ret)+bmu*(phi(6)+ret*ret))*amm2
      dmup(9)=aa*(ret*(3*tk+20*ak0)-20*ak0)*amm/zf
      dmup(10)=aa*(12*tk+40*ak0)/zf2
c
      GOTO 18
   15 dmps=dmup(1)
c
   18 idmx=10
      IF(nosd) idmx=4
      DO 19 i=1,10
      anu(i)=0.e0
      ue(i)=0.e0
      anur(i)=0.e0
   19 uer(i)=0.e0
c
c
      lj=0
      tki=1.d0/tk
      bdcozl=log(bdcoz)
c
c  ***************************************
c
c  start loop over elements
c
      DO 50 j=1,irab
      izoj=iz(ir(j))
      izj=izr(j)
c..      WRITE(*,*) 'j, izj, izoj', j, izj, izoj
c
c  skip elements with zero abundance
c
      IF(abr(j) <= 0) THEN
        anur(j)=0
        uer(j)=0
        xhvmn(j)=0
        lj=lj+izj
        GOTO 50
      ENDIF
c
c  set limit for no ionization
c
      dptstn=dptst1
      IF(j==jfirst) dptstn=dptst0
c
c  set exponents phi in saha equations into xi
c
      lj0=lj
      DO 20 i=1,izj
      lj=lj+1
      xhv(lj)=0
      dxp=i*dmps-chir(lj)*tki+bdcozl
   20 xi(i)=dxp
c
c  -----------------
c
c  test for complete or no ionization
c
      IF(izj-1) 21,21,25
c
c  only one level
c
   21 phm=xi(1)
      IF(phm <= 0.d0)THEN
       imax=0 ; phm=0.d0
c  test for no ionization
       IF(xi(1)<-dptstn) GOTO 50
       GOTO 34
      ENDIF

      imax=1
c  test for complete ionization
      IF(xi(1) > dptst1) GOTO 29
      GOTO 34
c
c  more than one level
c
   25 ii=izj+1
      xi(ii)=0
      imax=1
      phm=xi(1)
c
c  set phm to largest xi
c
      DO 26 i=2,ii
      IF(xi(i) <= phm) GOTO 26
      imax=i
      phm=xi(i)
   26 CONTINUE
c
      IF(imax /= izj) GOTO 30
c
c  test for complete ionization
c
      izj1=izj-1
      dphm=phm-xi(1)
      IF(izj1==1) GOTO 28
      DO 27 i=1,izj1
   27 dphm=min(dphm,phm-xi(i))
   28 IF(dphm <= dptst1) GOTO 34
c
c  complete ionization
c
   29 xhv(lj)=1
c..      WRITE(*,*) 'complete ionization at lj =',lj
      fct1=abr(j)/amr(j)
      anur(j)=fct1*izj
      uer(j)=fct1*chir(lj)
      anu(1)=anu(1)+anur(j)
      ue(1)=ue(1)+uer(j)
      xhvmn(j)=1
      GOTO 50
c
   30 IF(imax /= ii) GOTO 34
c
c  test for no ionization
c
      DO 31 i=1,izj
      IF(xi(i)>-dptstn) GOTO 34
   31 CONTINUE
c
c  no ionization. skip element completely
c
      xhvmn(j)=0
c
      GOTO 50
c
c  ************************************
c
c  general case
c
   34 DO 35 i=1,idmx
      dr(i)=0.e0
      hr(i)=0.e0
   35 gr(i)=0.e0
      IF(phm <= dptst1) dr(1)=omegac(0,izoj)*exp(-phm)
      DO 40 i=1,izj
      lji=lj0+i
      cchi=chir(lji)
      dxp=xi(i)-phm
      IF(dxp<-dptst1) GOTO 40
      dxp=omegac(i,izoj)*exp(dxp)
c..      PRINT *,' j,i,dxp',j,i,dxp
      xhv(lji)=dxp
      eir(1)=1
      DO 36 k=2,4
   36 eir(k)=i*dmup(k)
      eir(3)=eir(3)+cchi*amm*tki
      ii=4
      IF(nosd) GOTO 38
      DO 37 k=2,4
      DO 37 l=k,4
      ii=ii+1
   37 eir(ii)=eir(k)*eir(l)+i*dmup(ii)
      eir(8)=eir(8)-amm2*cchi*tki
c
   38 DO 39 k=1,idmx
      eeir=dxp*eir(k)
      dr(k)=dr(k)+eeir
      hr(k)=hr(k)+i*eeir
   39 gr(k)=gr(k)+cchi*eeir
   40 CONTINUE
c
      dr1i=1.d0/dr(1)
c  scale xhv
      DO 42 i=1,izj
      lji=lj0+i
      xhv(lji)=dr1i*xhv(lji)
c..      WRITE(*,*) 'At i, j, lji =',i,j,lji,'  xhv =',xhv(lji)
   42 CONTINUE
c
      fct1=abr(j)/(amr(j)*dr(1))
      hrdr=hr(1)*dr1i
      grdr=gr(1)*dr1i
      anur(j)=fct1*hr(1)
      uer(j)=fct1*gr(1)
      anu(1)=anu(1)+anur(j)
      ue(1)=ue(1)+uer(j)
c
      xhvmn(j)=hr(1)/(izj*dr(1))
c
c..      PRINT *,' j,anur(j),anu(1)',j,anur(j),anu(1)
c  derivatives
      ii=4
      DO 48 k=2,4
      anu(k)=anu(k)+fct1*(hr(k)-dr(k)*hrdr)
      ue(k)=ue(k)+fct1*(gr(k)-dr(k)*grdr)
      IF(nosd) GOTO 48
      DO 45 l=k,4
      ii=ii+1
      anu(ii)=anu(ii)+fct1*(hr(ii)-(dr(k)*hr(l)+dr(l)*hr(k)+(dr(ii)
     .  -2*dr(k)*dr(l)*dr1i)*hr(1))*dr1i)
   45 ue(ii)=ue(ii)+fct1*(gr(ii)-(dr(k)*gr(l)+dr(l)*gr(k)+(dr(ii)
     .  -2*dr(k)*dr(l)*dr1i)*gr(1))*dr1i)
   48 CONTINUE
   50 CONTINUE
c..      IF(idgeos>=3) WRITE(istdpr,250) (xhv(i),i=1,nlvtot)
c
      RETURN
  110 FORMAT(//' **** error in s/r hvionac. ihvz =',i5,
     *  ' not allowed')
  200 FORMAT(///' Original heavy element DATA.'/)
  205 FORMAT(' Name, abundance, Z, ionization potentials:'/)
  210 FORMAT(1x,a4,f8.5,i4,13f8.2)
  220 FORMAT(///' Heavy element DATA after resetting:'/)
  230 FORMAT(//' now anz0 =',f10.5//)
  250 FORMAT(/' xhv:'/(1p10e12.4))
      END SUBROUTINE hvionac

c*****************************************************************************

	SUBROUTINE inteffc(tl,pgl,rl,flini,wt,icase,iextr)
c>>>>>
c>>>>> finds starting value of flini for eff routines if the argument
c>>>>> is either log10(gas pressure) or log10(density). iextr is
c>>>>> set to 1 if extrapolation outside the table had to be made.
c>>>>>
c>>>>> based on eff-pre-computation with x=.73,y=.25,z=.02
c>>>>>
c>>>>> contains no COMMON statements. PARAMETERs are local.
c>>>>>
c===== icase = 1: input=(tl,pgl),  rl ignored
c===== icase = 2: input=(tl,rl ), pgl ignored
c
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
	IMPLICIT INTEGER(i-n)
      PARAMETER (nm = 11,nml = 5,nmh = 6,ntl = nml*nm,nth = nmh*nm)
c
      DIMENSION tlog(nm),rlog(nm,nm),pglog(nm,nm),flog(nm,nm)
      DIMENSION rll(ntl),rlh(nth),pgll(ntl),pglh(nth),fll(ntl),flh(nth)
      DIMENSION flini(2),wt(2)
c
      EQUIVALENCE(rlog (1,1),rll (1)) ,(rlog (1,nmh),rlh (1))
      EQUIVALENCE(pglog(1,1),pgll(1)) ,(pglog(1,nmh),pglh(1))
      EQUIVALENCE(flog (1,1),fll (1)) ,(flog (1,nmh),flh (1))
c
      DATA tlog/
     .    3.40000d0,   3.80000d0,   4.20000d0,   4.60000d0,   5.00000d0,
     .    5.40000d0,   5.80000d0,   6.20000d0,   6.60000d0,   7.00000d0,
     .    7.40000d0/
c
      DATA rll/
     .  -15.93806d0, -14.59077d0, -13.11393d0, -11.71751d0, -10.38835d0,
     .   -8.97767d0,  -7.59375d0,  -6.16023d0,  -4.74661d0,  -3.34048d0,
     .   -1.93629d0,
     .  -14.39806d0, -13.08988d0, -11.78384d0, -10.46730d0,  -9.15835d0,
     .   -7.85767d0,  -6.54375d0,  -5.24023d0,  -3.92661d0,  -2.62048d0,
     .   -1.31627d0,
     .  -12.79806d0, -11.57078d0, -10.36165d0,  -9.14343d0,  -7.92835d0,
     .   -6.71767d0,  -5.50375d0,  -4.28023d0,  -3.06660d0,  -1.85045d0,
     .   -0.63617d0,
     .  -11.19806d0, -10.06214d0,  -8.95086d0,  -7.82325d0,  -6.69833d0,
     .   -5.57766d0,  -4.45374d0,  -3.33021d0,  -2.20654d0,  -1.08028d0,
     .    0.04407d0,
     .   -9.59806d0,  -8.56014d0,  -7.52659d0,  -6.49503d0,  -5.46791d0,
     .   -4.43760d0,  -3.40366d0,  -2.37003d0,  -1.33614d0,  -0.30961d0,
     .    0.72372d0/
      DATA rlh/
     .   -7.99806d0,  -7.05242d0,  -6.11177d0,  -5.17967d0,  -4.23200d0,
     .   -3.29677d0,  -2.35276d0,  -1.41861d0,  -0.47419d0,   0.46015d0,
     .    1.40286d0,
     .   -6.39806d0,  -5.55147d0,  -4.69828d0,  -3.84877d0,  -3.00142d0,
     .   -2.15675d0,  -1.30499d0,  -0.46287d0,   0.38404d0,   1.23722d0,
     .    2.08971d0,
     .   -4.79807d0,  -4.04144d0,  -3.28452d0,  -2.53132d0,  -1.76802d0,
     .   -1.01603d0,  -0.25685d0,   0.49360d0,   1.25430d0,   2.00696d0,
     .    2.76640d0,
     .   -3.19831d0,  -2.53173d0,  -1.86318d0,  -1.20584d0,  -0.53947d0,
     .    0.12553d0,   0.78489d0,   1.45261d0,   2.11985d0,   2.77988d0,
     .    3.44720d0,
     .   -1.59748d0,  -1.02039d0,  -0.44723d0,   0.12189d0,   0.69032d0,
     .    1.26763d0,   1.83971d0,   2.41156d0,   2.98379d0,   3.55092d0,
     .    4.12641d0,
     .    0.00202d0,   0.48728d0,   0.96517d0,   1.44458d0,   1.92401d0,
     .    2.40470d0,   2.88056d0,   3.36421d0,   3.84159d0,   4.32331d0,
     .    4.80062d0/
c
      DATA pgll/
     .   -4.71634d0,  -2.70110d0,  -0.79422d0,   1.01898d0,   2.74815d0,
     .    4.55883d0,   6.34275d0,   8.17627d0,   9.98989d0,  11.79601d0,
     .   13.60020d0,
     .   -3.17634d0,  -1.30952d0,   0.53582d0,   2.26909d0,   3.97815d0,
     .    5.67883d0,   7.39275d0,   9.09627d0,  10.80989d0,  12.51602d0,
     .   14.22022d0,
     .   -1.57634d0,   0.09093d0,   1.95691d0,   3.59094d0,   5.20815d0,
     .    6.81883d0,   8.43275d0,  10.05627d0,  11.66990d0,  13.28603d0,
     .   14.90028d0,
     .    0.02366d0,   1.56549d0,   3.35740d0,   4.90076d0,   6.43816d0,
     .    7.95883d0,   9.48275d0,  11.00628d0,  12.52992d0,  14.05612d0,
     .   15.58052d0,
     .    1.62366d0,   3.06194d0,   4.76523d0,   6.22483d0,   7.66836d0,
     .    9.09886d0,  10.53279d0,  11.96637d0,  13.40014d0,  14.82664d0,
     .   16.26059d0/
      DATA pglh/
     .    3.22366d0,   4.56932d0,   6.08520d0,   7.53751d0,   8.90120d0,
     .   10.23927d0,  11.58325d0,  12.91715d0,  14.26164d0,  15.59761d0,
     .   16.94023d0,
     .    4.82375d0,   6.07052d0,   7.37716d0,   8.83489d0,  10.11675d0,
     .   11.37463d0,  12.62823d0,  13.87314d0,  15.12352d0,  16.37610d0,
     .   17.63041d0,
     .    6.42746d0,   7.58888d0,   8.76598d0,  10.05827d0,  11.32914d0,
     .   12.52233d0,  13.69271d0,  14.83754d0,  16.00047d0,  17.15703d0,
     .   18.32207d0,
     .    8.15380d0,   9.30678d0,  10.50206d0,  11.73099d0,  12.68725d0,
     .   13.70433d0,  14.76934d0,  15.84437d0,  16.91982d0,  17.98796d0,
     .   19.06506d0,
     .   10.79939d0,  11.94400d0,  13.08449d0,  13.22093d0,  14.17459d0,
     .   15.14145d0,  16.10082d0,  17.05997d0,  18.01959d0,  18.97151d0,
     .   19.93379d0,
     .   13.96914d0,  14.93979d0,  14.51362d0,  15.31529d0,  16.11750d0,
     .   16.92229d0,  17.71977d0,  18.53028d0,  19.33067d0,  20.13707d0,
     .   20.93387d0/
c
      DATA fll/
     .  -15.20875d0, -12.62213d0, -11.68335d0, -10.85416d0, -10.12500d0,
     .   -9.31434d0,  -8.53047d0,  -7.69708d0,  -6.88379d0,  -6.07850d0,
     .   -5.27639d0,
     .  -13.66875d0, -11.40213d0, -10.35335d0,  -9.60416d0,  -8.89500d0,
     .   -8.19434d0,  -7.48047d0,  -6.77708d0,  -6.06379d0,  -5.35850d0,
     .   -4.65639d0,
     .  -12.06875d0, -10.50213d0,  -8.93335d0,  -8.28416d0,  -7.66500d0,
     .   -7.05434d0,  -6.44047d0,  -5.81708d0,  -5.20379d0,  -4.58850d0,
     .   -3.97639d0,
     .  -10.46875d0,  -9.61213d0,  -7.54335d0,  -6.98416d0,  -6.43500d0,
     .   -5.91434d0,  -5.39047d0,  -4.86708d0,  -4.34379d0,  -3.81850d0,
     .   -3.29639d0,
     .   -8.86875d0,  -8.40213d0,  -6.15335d0,  -5.66416d0,  -5.20500d0,
     .   -4.77434d0,  -4.34047d0,  -3.90708d0,  -3.47379d0,  -3.04850d0,
     .   -2.61639d0/
      DATA flh/
     .   -7.26875d0,  -6.92213d0,  -4.97335d0,  -4.35416d0,  -3.97500d0,
     .   -3.63434d0,  -3.29047d0,  -2.95708d0,  -2.61379d0,  -2.27850d0,
     .   -1.93639d0,
     .   -5.66875d0,  -5.42213d0,  -4.11335d0,  -3.09416d0,  -2.77500d0,
     .   -2.50434d0,  -2.25047d0,  -2.00708d0,  -1.75379d0,  -1.49850d0,
     .   -1.24639d0,
     .   -4.06875d0,  -3.91213d0,  -3.30335d0,  -2.10416d0,  -1.64500d0,
     .   -1.40434d0,  -1.20047d0,  -1.03708d0,  -0.87379d0,  -0.71850d0,
     .   -0.55639d0,
     .   -2.46875d0,  -2.40213d0,  -2.26335d0,  -1.25416d0,  -0.26500d0,
     .   -0.17434d0,  -0.11047d0,  -0.03708d0,   0.03621d0,   0.10150d0,
     .    0.17361d0,
     .   -0.85875d0,  -0.88213d0,  -0.90335d0,   1.24584d0,   1.20500d0,
     .    1.17566d0,   1.13953d0,   1.10292d0,   1.06621d0,   1.02150d0,
     .    0.98361d0,
     .    0.92125d0,   0.77787d0,   3.15665d0,   2.99584d0,   2.83500d0,
     .    2.67566d0,   2.50953d0,   2.35292d0,   2.18621d0,   2.02150d0,
     .    1.84361d0/
c
      it=-999
      ip=-999
      ir=-999
      iextr=0
c========== select isotherms for linear inter(extra)polation
      DO 11 i=1,nm
      IF(tl<tlog(i)) GOTO 20
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
      x0 =tlog(it)
      x1 =tlog(it1)
      x2 =x0
      x3 =x1
      x  =  tl
c
c................................. icase = 1 ............................
      IF(icase==2) GOTO 200
c
c========== select pressure points for linear inter(extra)polation
      DO 21 i=1,nm
      IF(pgl<pglog(it,i)) GOTO 30
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
      y0 =pglog(it,ip)
      y1 =pglog(it1,ip)
      y2 =pglog(it,ip1)
      y3 =pglog(it1,ip1)
      z0 =flog(it,ip)
      z1 =flog(it1,ip)
      z2 =flog(it,ip1)
      z3 =flog(it1,ip1)
c
c========== define arguments
      y  = pgl
      GOTO 1000
c................................. icase = 2 ............................
c
c========== select density points for linear inter(extra)polation
 200  DO 31 i=1,nm
      IF(rl<rlog(it,i)) GOTO 40
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
      y0 =rlog(it,ir)
      y1 =rlog(it1,ir)
      y2 =rlog(it,ir1)
      y3 =rlog(it1,ir1)
      z0 =flog(it,ir)
      z1 =flog(it1,ir)
      z2 =flog(it,ir1)
      z3 =flog(it1,ir1)
c
c========== define argument
      y  =  rl
c
c========== CALL bilincear interpolation
c
1000  CONTINUE
c
c....... lower triangle
      CALL  bilinc (x0,y0,z0,x1,y1,z1,x2,y2,z2,x,y,z)
      flini(1) = z
c
c....... upper triangle
      CALL  bilinc (x3,y3,z3,x2,y2,z2,x1,y1,z1,x,y,z)
      flini(2) = z

c....... weights (quite arbitrary)
      wlow = 1./( (x-x0)**2 + (y-y0)**2 + 1.e-5 )
      whig = 1./( (x-x3)**2 + (y-y3)**2 + 1.e-5 )
      wtot = wlow + whig
      wlow = wlow/wtot
      whig = whig/wtot
      wt(1) = wlow
      wt(2) = whig

      RETURN
      END SUBROUTINE inteffc

c******************************************************************************

	FUNCTION ismax(n,a,nstep)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
	IMPLICIT INTEGER(i-n)
      DIMENSION a(*)
c
      ndim=n*nstep
      x=a(1)
      ismax=1
      k=1
      DO 1 i=1,ndim,nstep
      IF(a(i) <= x)GOTO 1
      ismax=k
      x=a(i)
    1 k=k+1
      END FUNCTION ismax
c
c******************************************************************************
c
	FUNCTION isamax(n,a,nstep)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
	IMPLICIT INTEGER(i-n)
      DIMENSION a(*)
c
      ndim=n*nstep
      x=ABS(a(1))
      isamax=1
      k=1
      DO 1 i=1,ndim,nstep
      aa=ABS(a(i))
      IF(aa <= x)GOTO 1
      isamax=k
      x=aa
    1 k=k+1
      RETURN
      END FUNCTION isamax
c
c*****************************************************************************
c
	FUNCTION ismin(n,a,nstep)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
	IMPLICIT INTEGER(i-n)
      DIMENSION a(*)
      ndim=n*nstep
      ismin=1
      k=1
      x=a(1)
      DO 1 i=1,ndim,nstep
      IF(a(i)>=x)GOTO 1
      ismin=k
      x=a(i)
    1 k=k+1
      END FUNCTION ismin
C
C******************************************************************************
C     Bernard PICHON    Janvier 2003
C
      Subroutine LEQ_BP ( A_mat , B_mat , nn , mm , ia , ib , err )
	 
         Integer :: nn, mm, ia, ib
         Double Precision, Dimension(ia,nn) :: A_mat
         Double Precision, Dimension(ib,mm) :: B_mat
         Double Precision :: err
C
         Double Precision, Dimension(SIZE(A_mat)) :: A_vect
         Double Precision, Dimension(SIZE(B_mat)) :: B_vect
C
         If ( ia < nn ) STOP " LEQ_BP : PB 1 "
         If ( ib < nn ) STOP " LEQ_BP : PB 2 "
C
         If ( SIZE(A_mat) > 100 ) STOP " LEQ_BP : ¨PB 3 "
         If ( SIZE(B_mat) > 100 ) STOP " LEQ_BP : ¨PB 4 "
C
         A_vect = RESHAPE( A_mat , SHAPE(A_vect) )
         B_vect = RESHAPE( B_mat , SHAPE(B_vect) )
C
         Call leq(A_vect,B_vect,nn,mm,ia,ib,err)
C
         A_mat = RESHAPE( A_vect , SHAPE(A_mat) )
         B_mat = RESHAPE( B_vect , SHAPE(B_mat) ) 
C
      End Subroutine LEQ_BP
C
      SUBROUTINE leq(a,b,nn,mm,ia,ib,err)
c
c     this routine will find the inverse of a matrix by the method of
c     partial pivoting and gaussian elimination if b is set equal to
c     the identity matrix
c
c     nn - DIMENSION of segment of a to be used
c     mm - number of right hand columns of b to be used
c     ia - the total number of rows in large array a
c     ib - the total number of rows in large array b
c     the matrix equation is    ax=b
c     err = det(a)
c     if mm = 0 leq calculates err = det(a)

c  note on modification on 23/1 1985:

c  previously, the routine contained the statements

c      DO 14 k=i2,m1,ib
c      b(k)=b(k)+b(i1)*r
c   14 i1=i1+ib

c  this caused problems, on some computers, for m = ib = 1.
c  THEN m1 = 1 and the loop was skipped when i2 > 1.
c  this has been changed.

c  Double precision version.
c  ++++++++++++++++++++++++

c  Dated 10/3/90.

c  Note: this DOUBLE PRECISION version of the routine has same name
c  as single precision version, but is distinguished by its file name

      IMPLICIT DOUBLE PRECISION(a-h,o-z)
	IMPLICIT INTEGER(i-n)
      DIMENSION a(100),b(100)
10000 n=nn
      m=mm
      err=0.0d0
      detsc=1.0d0
c
c     treat n  <=  1 separately
      IF(n-1) 280,281,285
c     no equation
  280 WRITE(*,295) n
      RETURN
c     n = 1
  281 err=a(1)
      IF(m <= 0) RETURN
      ai=1.d0/err
      m1=ib*m
      DO 282 j=1,m1,ib
  282 b(j)=b(j)*ai
      RETURN

c     find maximum element in each row and divide the row by it
  285 n1=ia*n
      m1=ib*m
      DO 1 i=1,n
      r= ABS(a(i))
      DO 2 j=  1,n1,ia
      ij=j+i-1
    2 r=max(r, ABS(a(ij)))   
      IF(r== 0.d0)THEN
       WRITE(*,298)i ; RETURN
      ENDIF
      
      DO 3 j=1,n1,ia
      ij=j+i-1
    3 a(ij)=a(ij)/r
      IF(m==0) GOTO 1
      DO 4 j=1,m1,ib
      ij=j+i-1
    4 b(ij)=b(ij)/r
    1 detsc=detsc*r
c
c
c     find maximum element in the i'th column
      n2=n-1
      DO 5 i=1,n2
      ialow=(i-1)*ia+i
      iaup=(i-1)*ia+n
      r= ABS(a(ialow))
      ialow1=ialow+1
      imax=ialow
      DO j=ialow1,iaup
       IF(r-ABS(a(j)) < 0.d0)THEN
        imax=j ; r= ABS(a(j))
       ENDIF
      ENDDO

      IF(imax-ialow)8,8,9
c     replace the i'th row with the row that has the maximum element in
c         the respective column and put the i'th row in its place
    9 im=imax
   72 IF(im-ia)70,70,71
   71 im=im-ia
      GOTO 72
   70 DO 10 j=1,n1,ia
      jj=i+j-1
      ji=im+j-1
      r=a(jj)
      a(jj)=a(ji)
   10 a(ji)=r
c     change sign of determinant
      detsc=-detsc
c
      IF(m==0) GOTO 8
      DO 11 j=1,m1,ib
      jj=i+j-1
      ji=im+j-1
      r=b(jj)
      b(jj)=b(ji)
   11 b(ji)=r
c     multiply the i'th row by (the negative of each i'th column element
c       below the diagonal divided by the diagonal element) and add the
c     resulting row to the respective row of the element used
    8 iaup1=iaup-1

      DO 12 j=ialow,iaup1
      IF(a(ialow) /= 0.d0)GOTO32
      joy=i
      IF(a(ialow1) == 0.d0)THEN
       WRITE(*,299)joy,joy ; RETURN
      ENDIF
      WRITE(*,297)joy,joy
      DO 34 k=1,n1,ia
      jj=joy+k-1
      ji=joy+k
      IF(joy+1-n)35,36,36
   35 WRITE(*,296)
      RETURN
   36 r=a(jj)
      a(jj)=a(ji)
   34 a(ji)=r
c     change sign of determinant
      detsc=-detsc
c
      IF(m==0) GOTO 8
      DO 37 k=1,m1,ib
      jj=joy+k-1
      ji=joy+k
      r=b(jj)
      b(jj)=b(ji)
   37 b(ji)=r
      GOTO 8
   32 j1=j+1
      r=-a(j1)/a(ialow)
      i1=ialow
      DO 13 k=j1,n1,ia
      a(k)=a(k)+a(i1)*r
   13 i1=i1+ia
c
c  loop to reset b has been modified, 25/1/1985.
c
      IF(m==0) GOTO 12
      i1=i
      i2=j-ialow+i+1
      DO 14 k=1,m1,ib
      b(i2)=b(i2)+b(i1)*r
      i1=i1+ib
   14 i2=i2+ib
   12 CONTINUE
c
c
    5 CONTINUE
c
c
c     the matrix is now in triangular form
c     first set err=1.0d0
      err=1.0d0
c     IF(any diagonal element of a is zero x cannot be solved for
      DO i=1,n
      idiag=(i-1)*ia+i
      err=err*a(idiag)
       IF(err== 0.d0)THEN
        WRITE(*,299)i,i ; RETURN
       ENDIF
      ENDDO
      
c     scale determinant
      err=err*detsc

      IF(m==0) RETURN
c     find solution to ax=b
      DO 18 k=1,m
      ka=(n-1)*ia+n
      kb=(k-1)*ib+n
      b(kb)=b(kb)/a(ka)
      DO 19 l=1,n2
      i=n-l
      r=0.0d0
      imax=i+1
      DO 20 j=imax,n
      jj=i+n+1-j
      ja=(jj-1)*ia+i
      jb=(k-1)*ib+jj
   20 r=r+a(ja)*b(jb)
      la=(i-1)*ia+i
      lb=(k-1)*ib+i
   19 b(lb)=(b(lb)-r)/a(la)
   18 CONTINUE
c
c
      RETURN
c
  295 FORMAT(///20h leq CALLed with n =,i4)
  296 FORMAT(///48h the row cannot be changed with the row below it  ,
     .40h because it it the last row in the array   /
     .55h the solution, matrix x, cannot be found by this method  )
  297 FORMAT(5h1  a(,i4,1h,i4,13h) equals zero  /
     .47h   try switching this row with the row below it   ,
     .48h and go back to statement number 8 and try again  )
  298 FORMAT(26h1  all the elements in row,i5,20h  are zero therefore,
     .55h the solution, matrix x, cannot be found by this method )
  299 FORMAT(50h1  the solution, matrix x, xannot be found by this  ,
     .57h method because there is a zero array element in the main ,
     .9h diagonal  / 30x,2ha(,i4,1h,i4,8h) = zero  )
      END SUBROUTINE leq
C
c******************************************************************************
C
	SUBROUTINE ferdir(x,fd)

c	calcul des integrales de Fermi-Dirac from MHD package

c     fd(1) = f    (x) 
c              -1/2 

c     fd(2) = f    (x)
c               1/2

c     fd(3) = f    (x)
c               3/2

c     fd(4) = f'   (x)
c              -1/2

c     fd(5) = f"   (x)
c              -1/2

c	version F95 : P. Morel, Departement J.D. Cassini, O.C.A.

c--------------------------------------------------------------------------

	USE mod_kind
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: x
	REAL (kind=dp), INTENT(out) :: fd(5)
	     
	REAL (kind=dp), DIMENSION(5), SAVE :: p1, p2, p3, p4, p5, p6,
	1	p7, p8, p9, q1, q2, q3, q4, q5, q6, q7, q8, q9
	
	REAL (kind=dp) :: y, p, dr, d2p, q, dq, d2q, root, xsq	

	DATA p1
	1/-1.253314128820d+0,-1.723663557701d+0,-6.559045729258d-1,
	2 -6.342283197682d-2,-1.488383106116d-5/
	DATA q1
	1/+1.000000000000d+0,+2.191780925980d+0,+1.605812955406d+0,
	2 +4.443669527481d-1,+3.624232288112d-2/
	DATA p2
	1/-3.133285305570d-1,-4.161873852293d-1,-1.502208400588d-1,
	2 -1.339579375173d-2,-1.513350700138d-5/
	DATA q2
	1/+1.000000000000d+0,+1.872608675902d+0,+1.145204446578d+0,
	2 +2.570225587573d-1,+1.639902543568d-2/
	DATA p3
	1/-2.349963985406d-1,-2.927373637547d-1,-9.883097588738d-2,
	2 -8.251386379551d-3,-1.874384153223d-5/
	DATA q3
	1/+1.000000000000d+0,+1.608597109146d+0,+8.275289530880d-1,
	2 +1.522322382850d-1,+7.695120475064d-3/
	DATA p4
	1/+1.073812769400d+0,+5.600330366000d+0,+3.688221127000d+0,
	2 +1.174339281600d+0,+2.364193552700d-1/
	DATA q4
	1/+1.000000000000d+0,+4.603184066700d+0,+4.307591067400d-1,
	2 +4.215113214500d-1,+1.183260160100d-2/
	DATA p5
	1/+6.781766266600d-1,+6.331240179100d-1,+2.944796517720d-1,
	2 +8.013207114190d-2,+1.339182129400d-2/
	DATA q5
	1/+1.000000000000d+0,+1.437404003970d-1,+7.086621484500d-2,
	2 +2.345794947350d-3,-1.294499288350d-5/
	DATA p6
	1/+1.153021340200d+0,+1.059155897200d+0,+4.689880309500d-1,
	2 +1.188290878400d-1,+1.943875578700d-2/
	DATA q6
	1/+1.000000000000d+0,+3.734895384100d-2,+2.324845813700d-2,
	2 -1.376677087400d-3, +4.646639278100d-5/
	DATA p7
	1/-8.222559330000d-1, -3.620369345000d+1,-3.015385410000d+3,
	2 -7.049871579000d+4, -5.698145924000d+4/
	DATA q7
	1/+1.000000000000d+0, +3.935689841000d+1,+3.568756266000d+3,
	2 +4.181893625000d+4, +3.385138907000d+5/
	DATA p8
	1/+8.224499762600d-1, +2.004630339300d+1,+1.826809344600d+3,
	2 +1.222653037400d+4, +1.404075009200d+5/
	DATA q8
	1/+1.000000000000d+0, +2.348620765900d+1,+2.201348374300d+3,
	2 +1.144267359600d+4, +1.658471590000d+5/
	DATA p9
	1/+2.467400236840d+0, +2.191675823680d+2,+1.238293790750d+4,
	2 +2.206677249680d+5, +8.494429200340d+5/
	DATA q9
	1/+1.000000000000d+0, +8.911251406190d+1,+5.045756696670d+3,
	2 +9.090759463040d+4, +3.899609156410d+5/
     
c--------------------------------------------------------------------------
 
	IF(x < 1.d0)THEN
	 y= dexp(x)
	 p   = y**2*(   p1(1) + y*(   p1(2) + y*(    p1(3) + y*(   p1(4)
	1                     + y*    p1(5)))))
	 dr  = y**2*(2.*p1(1) + y*(3.*p1(2) + y*( 4.*p1(3) + y*( 5.*p1(4)
	1                     + y* 6.*p1(5)))))
	 d2p = y**2*(4.*p1(1) + y*(9.*p1(2) + y*(16.*p1(3) + y*(25.*p1(4)
	1                     + y*36.*p1(5)))))
	 q   = q1(1) +y*(q1(2) + y*(   q1(3) + y*(   q1(4) + y*    q1(5))))
	 dq  =        y*(q1(2) + y*(2.*q1(3) + y*(3.*q1(4) + y* 4.*q1(5))))
	 d2q =        y*(q1(2) + y*(4.*q1(3) + y*(9.*q1(4) + y*16.*q1(5))))
	 fd(1) = 1.7724 53850 90552 d0*y + p/q
	 fd(2) = y*(0.8862 26925 45276 d0 + y*
	1          (p2(1) + y*(p2(2) + y*(p2(3) + y*(p2(4) + y*p2(5)))))/
	2          (q2(1) + y*(q2(2) + y*(q2(3) + y*(q2(4) + y*q2(5))))))
	 fd(3) = y*(1.3293 40388 17914 d0 + y*
	1          (p3(1) + y*(p3(2) + y*(p3(3) + y*(p3(4) + y*p3(5)))))/
	2          (q3(1) + y*(q3(2) + y*(q3(3) + y*(q3(4) + y*q3(5))))))
	 fd(4) = 1.7724 53850 90552 d0*y + (dr*q - p*dq)/q**2
	 fd(5) = 1.7724 53850 90552 d0*y +
	1               ((d2p*q - p*d2q)*q - 2.*(dr*q -p*dq)*dq)/q**3

	ELSEIF(x < 4.d0)THEN
	 p   =  p4(1)+ x*(p4(2) + x*(   p4(3) + x*(   p4(4) + x*   p4(5))))
	 dr  =            p4(2) + x*(2.*p4(3) + x*(3.*p4(4) + x*4.*p4(5)))
	 d2p =                       2.*p4(3) + 6.*x*(p4(4) + x*2.*p4(5))
	 q   =  q4(1)+ x*(q4(2) + x*(   q4(3) + x*(   q4(4) + x*   q4(5))))
	 dq  =            q4(2) + x*(2.*q4(3) + x*(3.*q4(4) + x*4.*q4(5)))
	 d2q =                       2.*q4(3) + 6.*x*(q4(4) + x*2.*q4(5))
	 fd(1) = p/q
	 fd(2) = (p5(1)  + x*(p5(2) + x*(p5(3) + x*(p5(4) + x*p5(5)))))/
	1        (q5(1)  + x*(q5(2) + x*(q5(3) + x*(q5(4) + x*q5(5)))))
	 fd(3) = (p6(1)  + x*(p6(2) + x*(p6(3) + x*(p6(4) + x*p6(5)))))/
	1        (q6(1)  + x*(q6(2) + x*(q6(3) + x*(q6(4) + x*q6(5)))))
	 fd(4) = (dr *q - p*dq )/q**2
	 fd(5) = (d2p*q - p*d2q)/q**2 - 2.*fd(4)*dq/q
     
	ELSE
	 root  = sqrt(x)
	 xsq   = x**2
	 y     = 1.0d0/xsq
	 p   =        y * (p7(1) + y*(   p7(2) + y*(   p7(3) + y*(    p7(4)
	1                        + y*    p7(5)))))
	 dr  = (-2.*y/x)* (p7(1) + y*(2.*p7(2) + y*(3.*p7(3) + y*( 4.*p7(4)
	1                        + y* 5.*p7(5)))))
	 d2p =  4.*y**2 * (p7(1) + y*(4.*p7(2) + y*(9.*p7(3) + y*(16.*p7(4)
	1                        + y*25.*p7(5))))) - dr/x
	 q   = q7(1) + y * (q7(2) +y*(   q7(3) +y*(   q7(4) +y*    q7(5))))
	 dq  = (-2.*y/x) * (q7(2) +y*(2.*q7(3) +y*(3.*q7(4) +y* 4.*q7(5))))
	 d2q =  4.*y**2  * (q7(2) +y*(4.*q7(3) +y*(9.*q7(4) +y*16.*q7(5))))
	1                -  dq/x
	 fd(1) = root * (2.0d0 + p/q)
	 fd(2) = x*root*(0.66666 66666 66667 d0 + y*
	1        (p8(1) + y*(p8(2) + y*(p8(3)   + y*(p8(4) + y*p8(5)))))/
	2        (q8(1) + y*(q8(2) + y*(q8(3)   + y*(q8(4) + y*q8(5))))))
	 fd(3) = xsq*root*(0.4d0   + y*
	1        (p9(1) + y*(p9(2) + y*(p9(3) + y*(p9(4) + y*p9(5)))))/
	2        (q9(1) + y*(q9(2) + y*(q9(3) + y*(q9(4) + y*q9(5))))))
	 fd(4) =    fd(1)/(2.d0*x) +   root*(dr*q - p*dq)/q**2
	 fd(5) = ( -fd(1)/x + fd(4))/(2.d0*x)
	1           + root* (dr *q - p*dq )/(2.*x*q**2)
	2           + root*((d2p*q - p*d2q) - 2.*dq*(dr*q - p*dq)/q)/q**2

	ENDIF
      
	END SUBROUTINE ferdir
C
c******************************************************************************
C
	SUBROUTINE neweta (t,vol,sne,ier)
c
c     calculate degeneracy PARAMETER and its                           c
c     derivatives. evaluate theta and its derivatives                  c
c     ................................
c     (modified from MHD package)
c     ................................

	IMPLICIT REAL*8 (a-h,o-z)
	IMPLICIT INTEGER(i-n)
c
      PARAMETER (mfd  =  5)
c
      COMMON /degmhd/
     .                                eta    ,               detdn  ,
     .                                detdt  ,               detdv  ,
     .                                d2etdn2,               d2etdnt,
     .                                d2etdnv,               d2etdt2,
     .                                d2etdtv,               d2etdv2,
     .                                exeta  ,               thet   ,
     .                                dthet  ,               d2thet ,
     .                                fd    (mfd   )
      COMMON /cstmhd/
     .                                camu   ,               cc     ,
     .                                ce     ,               ch     ,
     .                                ck     ,               cme    ,
     .                                cpi    ,               cevw   ,
     .                                cx     ,               cf4    ,
     .                                carad
c
      DATA   niter,   conv
     .    /    20, 1.d-10 /
	LOGICAL ok
	COMMON /marche/ok
c
c........... the following factorizing in gse is made to avoid
c........... overflows on simple machines with ranges of about 1.e37
c
      gse  = 2.*cpi*(cme/ch)*(ck/ch)
      gse  = gse**1.5
c
      ceta = sqrt( cpi ) / ( 4.0 * gse )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialization                                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ier = 0
      rhs = ceta * sne   / ( t**1.5 * vol )
ccc      WRITE(*,*) ' rhs,ceta = ',rhs,ceta
c
c........... initialization to a few per cent accuracy
c........... (see ref. in dappen, astron. astrphys., 1980)
c
      IF(rhs<1.d-6) THEN
            eta = dlog(rhs)
      ELSE IF(rhs<100.d0) THEN
            g   =rhs
            g2  =rhs*rhs
            et0 =g+(g2/(4.45+2.257*g+g2))*dexp((1.32934*g)**0.66666667)
            eta = dlog(et0)
      ELSE
            eta = (1.32934*rhs)**0.66666667
      ENDIF
      eta = eta + 0.120782
c
ccc      WRITE(*,*) 'eta,rhs ini: ',eta,rhs
      CALL ferdir( eta, fd )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     newton-raphson iteration
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      iter  = 0
    1 deta  = 2.0d0 * ( rhs - fd(2) ) / fd(1)
      eta = eta + deta
      CALL ferdir( eta, fd )
      IF( ABS(deta)  <=  conv ) GOTO 3
      iter = iter + 1
      IF( iter  <=  niter ) GOTO 1
c
c     failure to converge
      WRITE  ( 6, 2 ) t, vol, sne, eta
    2 FORMAT ( ' nonconvergence of degeneracy PARAMETER ' ,/,
     .         ' t =', 1pg10.2 ,' vol =', g10.2, ' ne =', g10.3,
     .         ' eta   =,' g12.4 )
	ok=.FALSE.
      ier = 1
      RETURN
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     convergence                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
    3 exeta = dexp( - eta )
c
      thet   =  0.5d0 *   fd(1) / fd(2)
      dthet  = thet * ( fd(4) / fd(1) - thet )
      d2thet = thet * ( fd(5) / fd(1) - 3.0d0 * dthet - thet**2 )
c
      detdv   = - 1.0d0 / ( thet * vol )
      detdt   = - 1.5d0 / ( thet * t )
      detdn   =   1.0d0 / ( thet * sne )
      d2etdn2 = - detdn**2 * fd(4) / fd(1)
      d2etdnt = - dthet * detdn * detdt / thet
      d2etdnv = - dthet * detdn * detdv / thet
      d2etdtv = - dthet * detdt * detdv / thet
      d2etdt2 = - ( 1.0d0 / t   + dthet * detdt / thet ) * detdt
      d2etdv2 = - ( 1.0d0 / vol + dthet * detdv / thet ) * detdv
c
      IF(iter>5) WRITE (6,*)' slow convergence in neweta: iter,eta',
     . ',t,vol = ',iter,eta,t,vol
c
ccc      WRITE(*,*) 'eta,fd ',eta,fd
      RETURN
      END SUBROUTINE neweta

c*****************************************************************************
c
	FUNCTION omegac(i,iz)
c  calculates statistical weight of i-th ionization stage of element
c  with number iz
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
	IMPLICIT INTEGER(i-n)
      COMMON/hvomeg/ iom(26),iom1(20)
      COMMON/hvomcl/ iomfll
c
      SAVE
c
      IF(i <= 1.and.iz>=19) GOTO 20
      IF(i==iz) GOTO 15
      omegac=iom(iz-i)
      RETURN
c
c  statistical weight for fully ionized atom.
c  before 5/1/84 was always set to 15, for some bizarre reason.
c  for transition introduce flag iomfll so that iomfll = 0 corresponds
c  to old situation, and iomfll  /=  0 gives correct value omegac = 1.
c
   15 omegac=1
      IF(iomfll==0) omegac=15
      RETURN
   20 omegac=iom1(2*(iz-19)+i+1)
      RETURN
      END FUNCTION omegac

c******************************************************************************

	SUBROUTINE phderc(fl,tl,phi,hst,nosd,notd)
c
c  computes quantities for Eggleton, Faulker & Flannery approximation
c  to partial degeneracy. On input fl and tl are log(F) and log(T).
c  Returns phi(1-30) and hst(1-10). Here phi(1) is defined such that
c  the density is rho = phi(1)*(crho/ane), where ane is the number
c  of electrons per unit mass, phi(11) is defined such that the
c  electron pressure is pe = cpe*phi(11), and hst(1) such that
c  the electron enthalpy per unit mass is He = che*ane*hst(1).
c  The constants crho, cpe and che are given in COMMON/constsc/.
c  phi(21) is related to hst. First, second and third derivatives
c  with respect to log(f) and log(T)
c  of log(phi(1)), log(phi(11)), log(phi(21)) and hst(1) are
c  given in phi(2 - 10), phi(12 - 20), phi(22 - 30) and hst(2 - 10),
c  respectively.

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
	IMPLICIT INTEGER(i-n)
      LOGICAL nosd,notd
      DIMENSION phi(30)
      DIMENSION sig(10),tmn(10),ff(4),gg(4),hst(10)
      COMMON/phdsmsc/ s0,sf,sg,sff,sfg,sgg,sfff,sffg,sfgg,sggg,
     .  cfg,tf,tg,tff,tfg,tgg,tfff,tffg,tfgg,tggg
      COMMON/constsc/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
      COMMON/eqphcs/ c(48),ic
      COMMON/ln10c/ amm,amm2,amm3
c  LOGICAL EQUIVALENCE
      EQUIVALENCE(s0,sig(1)),(cfg,tmn(1))
c     ***********************************

      SAVE

c  number of sums

      imax=10
      IF(notd) imax=6
      IF(nosd) imax=3
c  f,g
      f=1.d1**fl
      t=1.d1**tl
      ts=ct*t
      wf= sqrt(1+f)
      g=ts*wf

c  1/(1+f), 1/(1+g) etc.

      vf=1/(1+f)
      vg=1+g
      fdf=vg*g
      fdf=vf*fdf* sqrt(fdf)
      vg=1/vg
      vfg=vf*vg
      ug=g*vg
      uf=f*vf

      ug2=ug*vg
      ug3=vg*ug2
c  powers of f and g
      ff(1)=f
      gg(1)=1
      DO 10 i=2,ic
      ff(i)=f*ff(i-1)
      gg(i)=g*gg(i-1)
   10 fdf=fdf*vfg

c  test on size of f and g
c  icase is set as flag for size, as follows
c   icase = 1: general case
c   icase = 2: f small
c   icase = 3: g small
c   icase = 4: f and g small

      icase=1
      IF(.not.notd) GOTO 12
      IF(f<1.e-4) icase=icase+1
      IF(g<1.e-4) icase=icase+2

   12 ic2=ic*ic

c  calculate phi* and derivatives

      l=1
      anu=1.5e0
      mio=ic
      an32=2.5e0-ic
      kk=-10

      l0=1
      DO 50 k=1,3
      kk=kk+10
      IF(k-2) 18,16,17
c  reset fdf for k=2 and 3
   16 anu=2.5e0
      fdf=g*fdf
      GOTO 18
   17 mio=mio+1
      fdf=fdf*vf
   18 DO 19 i=1,imax
   19 sig(i)=0.e0
      annu=anu-1

c  the summation

      l=l0

c  select case, based on size of f and g, as determined by icase

      GOTO (20,25,30,35), icase
c  the general case
   20 DO 23 in=1,ic
      annu=annu+1
      DO 23 im=1,ic
c  phimn*(f**(m+1))*(g**n)
      cfg=c(l)*ff(im)*gg(in)

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

c  the sigma tilde (cf (22.2)) are stored in the corresponding sigma.
c  this is o.k. provided that we go backwards.

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
c  g is small
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
c  both f ang g are small
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
      IF(f*g<1.00001e-8) GOTO 42
      sfg=(sfg*s0-sf*sg)/s02
      GOTO 45
c  approximate expression for sfg (may need fixing up, if f = o(1)
c  or g = o(1))
   42 sfg=f*g*(c(l0+5)-c(l0+1)*c(l0+4)/c(l0))/c(l0)

c  phi* and first derivatives

   45 phi(kk+1)=fdf*s0
      pht=an32*ug+sg/s0
      phi(kk+3)=pht
      phi(kk+2)=(pht/2-mio)*uf+sf/s0
      IF(nosd) GOTO 50

c  second derivatives of phi*.

      phtt=an32*ug2+sgg
      phi(kk+6)=phtt
      phi(kk+5)=phtt*uf/2+sfg
      phi(kk+4)=sff+uf*(sfg+vf*(pht/2-mio+f*phtt/4))

      IF(notd) GOTO 50
c  third derivatives
      phttt=an32*ug3*(1-g)+sggg
      phi(kk+10)=phttt
      phi(kk+9)=sfgg+uf*phttt/2
      phfft=sffg+uf*(sfgg+vf*(phtt+f*phttt/2)/2)
      phi(kk+8)=phfft
      phi(kk+7)=sfff+uf*(sffg+phfft/2+vf*(1.5*sfg+f*sfgg/4
     .  +vf*((1-f)*(pht/2-mio)+f*phtt/2)))
   50 l0=l0+ic2

c  h* and its derivatives (pp 23-25)

      DO 55 i=2,imax
      ik=20+i
   55 phi(ik)=phi(ik)-phi(i)

      hs=phi(21)/phi(1)
      wft1=2*g
      hst(1)=hs+wft1

      hf=phi(22)
      ht=phi(23)
      wft2=ts*f/wf
      hst(2)=hs*hf+wft2
      hst(3)=hs*ht+wft1
      IF(nosd) GOTO 58
c  second derivatives
      hff=phi(24)
      hft=phi(25)
      htt=phi(26)
      wft3=uf*(1+f/2)*ts/wf
      hst(4)=hs*(hf*hf+hff)+wft3
      hst(5)=hs*(hf*ht+hft)+wft2
      hst(6)=hs*(ht*ht+htt)+wft1
      IF(notd) GOTO 58
c  third derivatives
      hst(7)=hs*(hf*(hf*hf+3*hff)+phi(27))+uf*vf*(1+f*(2+f)/4)*ts/wf
      hst(8)=hs*(hf*(hf*ht+2*hft)+ht*hff+phi(28))+wft3
      hst(9)=hs*(ht*(ht*hf+2*hft)+hf*htt+phi(29))+wft2
      hst(10)=hs*(ht*(ht*ht+3*htt)+phi(30))+wft1
c  change to derivatives wrt log10 f and log10 t
   58 fct=amm
      DO 60 i=2,imax
      IF(i==4.OR.i==7) fct=fct*amm
   60 hst(i)=hst(i)*fct
   
      RETURN
      
      END SUBROUTINE phderc

c*******************************************************************************

	SUBROUTINE prthvi(tl,xr,ihvz,fin)

c  output results of heavy element ionization calculations

      IMPLICIT DOUBLE PRECISION (a-h, o-z)
	IMPLICIT INTEGER(i-n)
      CHARACTER*(*) fin
      CHARACTER name*4
      DIMENSION ids(4), chipr(30)
      COMMON/potetcc/ chi(125),am(10),iz(10)
      COMMON/hvname/ name(10)
      COMMON/hvabndc/ ab(10),iab
      COMMON /hvredu/ irab,jfirst,ir(10),chir(125),izr(10),abr(10),
     *  amr(10)
      COMMON/hviond/ xhvmn(10),xhv(125)
      COMMON/eqstdc/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)

      INTEGER, SAVE :: init=0

c  test for initializing, printing file headers

      IF(init==0) THEN
        WRITE(20,105) fin
        init=1
        IF(ihvz==1) THEN
          ielem=2
          ids(1)=21
          ids(2)=23
        ELSE
          ielem=4
          DO i=1,ielem
            ids(i)=20+i
          ENDDO
        ENDIF

        lj0=0
        DO i=1,ielem

c  reset original ionization potentials

          chipr(1)=chir(lj0+1)
          IF(izr(i)>1) THEN
            DO j=2,izr(i)
              chipr(j)=chir(lj0+j)-chir(lj0+j-1)
            ENDDO
          ENDIF
          WRITE(ids(i),110) fin,name(ir(i)), ihvz,
     *      (chipr(j),j=1,izr(i))
          WRITE(ids(i),120)
          lj0=lj0+izr(i)
        ENDDO

        WRITE(91,125) fin

      ENDIF

      lj0=0
      WRITE(20,130) xr,tl,xii1(4),gm1,(xii1(i),i=1,3)
      DO i=1,ielem
        WRITE(ids(i),130) xr,tl,xii1(4),gm1,(xhv(lj0+j),j=1,izr(i))
        lj0=lj0+izr(i)
      ENDDO
      WRITE(91,130) xr,tl,xii1(4),gm1,(xhvmn(i),i=1,10)
      RETURN
  105 FORMAT('# Test ionization in model ',a60/'#'/
     *  '#  ionization of H and He'/'#'/
     *  '# r/R, log T, xi(heavy), Gamma_1, xi(H), xi(He), xi(He+):'/'#')
  110 FORMAT('# Test ionization in model ',a60/'#'/
     *  '#  ionization of ',a5/'#'/'# ihvz =',i5/'#'/
     *  '# ionization potentials:'/
     *  '#',13f9.2)
  120 FORMAT('#'/'# r/R, log T, xi(heavy), Gamma_1, x(element):'/'#')
  125 FORMAT('# Test ionization in model ',a60/'#'/
     *  '#  Average ionization of heavy elements'/'#'/
     *  '# r/R, log T, xi(heavy), Gamma_1, xihvmn(1-10):'/'#')
  130 FORMAT(2f10.5,1pe11.3,0pf10.5,1p13e11.3)
      END SUBROUTINE prthvi

c*****************************************************************************

	SUBROUTINE setgm1(tl,xh,xr,n,nn,fin)

c  outputs results for Gamma_1 derivatives to file

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
	IMPLICIT INTEGER(i-n)
      CHARACTER fin*(*)
      PARAMETER(nnmax=1301, iDATA=10, iw=5)
      DIMENSION xrr(nnmax), DATA(iDATA,nnmax), w(iw,nnmax)
      COMMON/constsc/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
      COMMON/eqstdc/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
c
c  set derivatives of Gamma_1 wrt rho, T, X. Other required
c  quantities into DATA.
c
      xrr(n)=xr
      DATA(1,n)=log10(rho(1))
      DATA(2,n)=tl
      DATA(3,n)=xh
      DATA(4,n)=gmm1(1)
      DATA(5,n)=gmm1(2)/rho(2)
      DATA(6,n)=gmm1(3)-gmm1(2)*rho(3)/rho(2)
      DATA(7,n)=gmm1(4)-gmm1(2)*rho(4)/rho(2)
c
c  test for final pass, setting of derivatives and output
c
      IF(n<nn) RETURN
c
      DO 20 i=1,4
   20 CALL derive(xrr,DATA(i,1),w(i,1),nn,iDATA,iw,1,1)
c
      WRITE(30,110) fin
c
      DO 30 n=1,nn
      dgm1rh=w(1,n)*DATA(5,n)
      dgm1t= w(2,n)*DATA(6,n)
      dgm1x= w(3,n)*DATA(7,n)
      dgmtot=dgm1rh+dgm1t+dgm1x
c
   30 WRITE(30,120) xrr(n),(DATA(i,n),i=1,7),
     *              dgm1rh,dgm1t,dgm1x,dgmtot,w(4,n)
c
c  for testing, output also derivatives
c
      open(92,file='ttt.rhotx.derout',status='unknown')
      WRITE(92,140) fin,(xrr(n),(DATA(i,n),i=1,4),(w(i,n),i=1,4),
     *  n =1,nn)
c
      RETURN
  110 FORMAT('# Gamma_1 derivative results.'/'#'/
     *  '# Model file: ',a60/'#'/
     *  '# r/R, log rho, log T, X, Gamma_1, Gamma_1,rho, Gamma_1,T,',
     *  ' Gamma_1,X, (d Gamma_1),rho, (d Gamma_1),T, (d Gamma_1),X,',
     *  ' Sum, d Gamma_1:'/'#')
  120 FORMAT(f12.7,4f10.6,1p3e13.5,5e11.3)
  140 FORMAT('# Derivatives.'/'#'/
     *  '#  Model file: ',a60/'#'/
     *  '#  r/R, log rho, log T, X, Gamma_1, derivatives:'/'#'/
     *  (0p5f12.6,1p4e12.4))
      END SUBROUTINE setgm1
c
c***************************************************************************
c
	SUBROUTINE setf4(nspe)
c
c  Note: COMMONs free, map and spes have been renamed freecl,
c  mapcl and spescl to avoid conflict when CALLing routine
c  from MHD table emulator.
c                                           1/5/90       jcd
c
c  Modified 11/5/90 to allow arbitrary order of expansion for tau
c  in s/r f4mhd (modification necessary because of unavailability
c  of REAL*16).                                          jcd
c
c  Modified 5/6/90 to reset constants with JC-D values for consistency.
c  Note: this requires that s/r setcnsc be CALLed before CALLing setf4.
c  In JC-D usage, setf4 would normally be CALLed from setcnsc.
c
c======================================================================
c
      IMPLICIT REAL*8 (a-h,o-z)
	IMPLICIT INTEGER(i-n)

      INTEGER, PARAMETER :: mspes=18, mfe=4, mfd=  5

      COMMON /cstmhd/
     .                                camu   ,               cc     ,
     .                                ce     ,               ch     ,
     .                                ck     ,               cme    ,
     .                                cpi    ,               cevw   ,
     .                                cx     ,               cf4    ,
     .                                carad
c
      COMMON /degmhd/
     .                                eta    ,               detdn  ,
     .                                detdt  ,               detdv  ,
     .                                d2etdn2,               d2etdnt,
     .                                d2etdnv,               d2etdt2,
     .                                d2etdtv,               d2etdv2,
     .                                exeta  ,               thet   ,
     .                                dthet  ,               d2thet ,
     .                                fd    (mfd   )
      COMMON /freecl/
     .                         e     (mfe   ),        p     (mfe   ),
     .                         f     (mfe   ),               d2fdt2 ,
     .                                d2fdtv ,               d2fdv2 ,
     .                         fscr  (mspes ),        dfdn  (mspes ),
     .                         d2fdnt(mspes ),        d2fdnv(mspes ),
     .                         d2fdn2( mspes ,mspes )
      COMMON /mapcl /
     .                                nchem  ,               nspes  ,
     .                                ise
      COMMON /spescl/
     .                         zs    (mspes ),        zsq   (mspes ),
     .                         zmask (mspes ),        sn    (mspes )
c
c  COMMONs of fundamental constants set by JC-D routine setcnsc.
c  for consistency, replace WD values by these.
c
      COMMON/constsc/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
c
c  COMMON giving fundamental constants
c
      COMMON /fconst/ pi, ebase, amu, ame, clight, planck, boltzm,
     *  cgrav, amsun, echar, ergev, syear, iver
c
      COMMON/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
      IF(idgeos>=1) WRITE(*,*) 'Calling setf4'
c
      camu  =  1.6605655d-24
      cc    =  2.9979246d+10
      ce    =  4.8032426d-10
      ch    =  6.6261764d-27
c
      ck    =  1.3806624d-16
      cme   =  9.1095345d-28
      cpi   =  3.141592654d0
      cevw  =  8.0654653d+03
c
c  resetting with JC-D values
c
      camu  =  amu
      cc    =  clight
      ce    =  echar
      ch    =  planck
      ck    =  boltzm
      cme   =  ame
      cpi   =  pi
c
c  Note: the meaning of cevw is currently unclear. On the other
c  hand, it appears not to be used. Should be fixed up.
 
c
      cf4   = 2.0d0 * sqrt( cpi ) * ce**3 / ( 3.0d0 * sqrt( ck ) )
      cx    = 3.0d0 * cf4 / ck
c
      carad =  7.56567  d-15
c
c  JC-D value
c
      carad =  car
c
      nspes =  nspe
      ise   =  nspe
c
      zmask(1) = 0.d0
      zmask(2) = 1.d0
      zmask(3) = 0.d0
      zmask(4) = 1.d0
      zmask(5) = 1.d0
c
      zs(1) = 0.d0
      zs(2) = 1.d0
      zs(3) = 0.d0
      zs(4) = 1.d0
      zs(5) = 2.d0
c
      zsq(1) = 0.d0
      zsq(2) = 1.d0
      zsq(3) = 0.d0
      zsq(4) = 1.d0
      zsq(5) = 4.d0
c
c................... following values should never be used; if none the less,
c................... they ought to provoke a crash of the program!
      zmask(6) = -1.d200
      zs(6)    = -1.d200
      zsq(6)   = -1.d200
c
      RETURN
      END SUBROUTINE setf4
c
c
c*******************************************************************************
c
C BP	SUBROUTINE seteqs
c
c  dummy SUBROUTINE included for compatibility with dog equation  of
c  state programmes.
c
C BP  RETURN
C BP  END SUBROUTINE seteqs
c

c*****************************************************************************
c
	SUBROUTINE setcnsc
c
c  sets the physical and mathematical constants used in the programme
c
c  numerical constants from CODATA report, with additional values
c  for atomic masses and ionization potentials. this is
c  the GONG set of values.
c
c  corresponds to iver = 2
c
c  modified on 8/12/87 to approach a consistent setting based on
c  minimum set of basic constants. so far keep original values,
c  but PRINT consistent values
c
c  modified 15/3/88 to use consistent numerical constants
c  note that a version number iver has been added in COMMON/fconst/
c  to indicate which set of values have been used.
c
c  also include constants for simple equation of state, with
c  no partial degeneracy, for use for GONG models.
c
c  This version also CALLs WD routine to prepare for inclusion
c  of Coulomb terms.
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 9/3/90
c
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
	IMPLICIT INTEGER(i-n)
      COMMON/constsc/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
c
c  COMMON giving fundamental constants
c
      COMMON /fconst/ pi, ebase, amu, ame, clight, planck, boltzm,
     *  cgrav, amsun, echar, ergev, syear, iver
      COMMON/eqphsm/ rho00, f0
      COMMON/ln10c/ amm,amm2,amm3
      COMMON/eqstdc/ ccc1(94)
      COMMON/eqsoutc/ ccc2(230)
      COMMON/dmudec/ ccc3(10),idmu
      COMMON/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
c
      SAVE
c
c     WRITE(*,100)
c
c  version number for constants
c
c
c  *****************************************************************
c
c  set equation of state version number to 1 for EFF with Coulomb
c  terms (note that actual value may depEND on icoulm, and hence
c  may be reset later, after CALL of eqstfc
c
      ivreos = 1
c
c  *****************************************************************
c
      iver = 2
c
      pi=4.d0*atan(1.d0)
      ebase=exp(1.d0)
c
c  atomic mass unit
      amu=1.6605402d-24
c  electron mass
      ame=9.1093897d-28
c  speed of light
      clight=2.99792458d10
c  planck's constant
      planck=6.6260755d-27
c  boltzman's constant, cgs
      boltzm=1.380658d-16
c  gravitational constant (chosen from planetary orbit DATA, assuming
c  msun = 1.989e33. consistent, but not identical, with CODATA value).
      cgrav=6.67232d-8
c  solar mass
      amsun=1.989d33
c  electron charge, in ESU
      echar=4.8032068d-10
c  number of ergs in 1 ev
      ergev=1.60217733d-12
c  number of seconds in a year
c  (should be replaced with second value)
      syear =3.155692597d7
c
c  start on derived values, largely
c
c     WRITE(*,110)
c  avogadro's number (should be 1/amu)
      av=1.0/amu
c     WRITE(*,120) 'av',avp,av
c  atomic weights of h and he
      ah=1.007825d0
      ahe=4.002603d0
c  average atomic weight for heavy elements.
c  before 4/1/84 was given value 17.8, but was in fact
c  most often reset to the current value in CALLing programme.
      az=16.389d0
c  auxiliary quantities for eggleton fudge calculation
      avda=av*(1/ah-2/ahe)
      avd1=av*(1/ah-1/ahe)
c  boltzmann's constant in ev/deg
      ck1=boltzm/ergev
c     WRITE(*,120) 'ck1',ck1p,ck1
c  the same, in ergs/deg
      ck2=boltzm
c  ionization potentials
      exh=13.595d0
      exhe=24.580d0
      exhep=54.403d0
c  constants for transition to starred variables

      alamc=planck/(ame*clight)
      alamc3=alamc**3
      ct=boltzm/(ame*clight**2)
      crho=8*pi/alamc3
      che=ame*clight**2
      cpe=crho*che
c     WRITE(*,120) 'ct',ctp,ct
c     WRITE(*,120) 'crho',crhop,crho
c     WRITE(*,120) 'cpe',cpep,cpe
c     WRITE(*,120) 'che',chep,che
c  constants for pressure ionization
      ca03p=2.147d-24
      a0=planck/(2*pi*echar)
      a0=(a0/ame)*a0
c..      WRITE(*,*) 'a0 =',a0
c
c  for the moment, set additional fudge factor by hand at this
c  point. need to think about a more suitable way to DO it later
c
      efffac=15
      ca03=efffac*a0**3
      caa=crho*ca03
      ckh=exh
c
c  change from ev to ergs and include factor av**2 to compensate for
c  redefinition of dne (see notes of 3/1/84)
c
      ca03p=av*1.602192d-12*ca03p*av
      ca03=av*ergev*ca03*av
c
c     WRITE(*,120) 'ca03',ca03p,ca03
c     WRITE(*,120) 'caa',caap,caa
c     WRITE(*,120) 'ckh',ckhp,ckh

c  the radiation constant

      car=boltzm/(clight*planck)
c..      WRITE(*,*) 'boltzm/(clight*planck)', carc
      car=8*pi**5*boltzm*car**3/15
c     WRITE(*,120) 'car',carp,car
c  number of ergs in 1 ev
      ergev1=ergev
c  ionization potential for h-
      exhm=0.754d0
c  ln 10
      amm=log(1.d1)
      amm2=amm*amm
      amm3=amm2*amm
c
c constant for phderc, with no partial degeneracy
c
      rho00=sqrt(2*pi)*ebase**2/8
      f0=av/(rho00*crho*ct**1.5d0)
c..     WRITE(*,*) 'f0',f0
c
c  set COMMONs from s/r eqstfc to zero
      ccc1 = 0.    !  CALL zeroc(ccc1,90)
      ccc2 = 0.    !  CALL zeroc(ccc2,210)
      ccc3 = 0.    !  CALL zeroc(ccc3,10)
c
c====================== initialize mhd constants ==================
      nspe=6
      CALL setf4(nspe)
c=====================================================================
c
      RETURN
  100 FORMAT(//2x,75('*')//
     *  ' Set constants from CODATA Report.'/' See Cohen & Taylor,',
     *  ' Rev. Mod. Phys., vol. 59, 1121 (1987)'/
     *  ' ++++ Double precision version'/)
  110 FORMAT(/' variable   old value    current, consistent value'/)
  120 FORMAT(a10,1p2e18.7)
      END SUBROUTINE setcnsc

c************************************************************************

c
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90.
c
c  Note: this DOUBLE PRECISION version of the routine has same name
c  as single precision version, but is distinguished by its file name
c
c Bernard PICHON, 2000-11-14
c
!!!	Double Precision Function SDOT ( n , a , ia , b , ib ) Result ( s )
!!!      Implicit None
!!!      Integer, INTENT(In) :: n , ia, ib
!!!      Double Precision, Dimension(:), INTENT(In) :: a, b
         !
!!!      Integer :: na, nb
	   !
!!!      na = (n-1)*ia + 1
!!!      nb = (n-1)*ib + 1
!!!      s = DOT_PRODUCT( a(1:na:ia) , b(1:nb:ib) )
         !
!!!   End Function SDOT

	FUNCTION sdot(n,a,i1,b,i2)

      IMPLICIT DOUBLE PRECISION(a-h,o-z)
	IMPLICIT INTEGER(i-n)
      DIMENSION a(*),b(*)
      sdot=0.
      j1=1
      j2=1
      DO 1 i=1,n
      sdot=sdot+a(j1)*b(j2)
      j1=j1+i1
      j2=j2+i2
    1 CONTINUE
      RETURN
      END FUNCTION sdot

c**************************************************************************

	FUNCTION ssum(n,a,i1)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
	IMPLICIT INTEGER(i-n)
      DIMENSION a(*)
      ssum=0.
      j1=1
      DO 1 i=1,n
      ssum=ssum+a(j1)
      j1=j1+i1
    1 CONTINUE
      RETURN
      END FUNCTION ssum
c
c******************************************************************************
c
	SUBROUTINE sscal(n,fak,a,i1)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
	IMPLICIT INTEGER(i-n)
      DIMENSION a(*)
      j1=1
      DO 1 i=1,n
      a(j1)=fak*a(j1)
    1 j1=j1+i1
      RETURN
      END SUBROUTINE sscal
c
c******************************************************************************
c
	SUBROUTINE scopy(n,a,na,b,nb)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
	IMPLICIT INTEGER(i-n)
      DIMENSION a(*),b(*)
      ia=1
      ib=1
      DO 1 i=1,n
      b(ib)=a(ia)
      ia=ia+na
      ib=ib+nb
    1 CONTINUE
      RETURN
      END SUBROUTINE scopy

c******************************************************************************

	SUBROUTINE saxpy(n,c,a,na,b,nb)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
	IMPLICIT INTEGER(i-n)
      DIMENSION a(*),b(*)
      ia=1
      ib=1
      DO 1 i=1,n
      b(ib)=c*a(ia)+b(ib)
      ia=ia+na
    1 ib=ib+nb
      RETURN
      END SUBROUTINE saxpy

c*****************************************************************************

	SUBROUTINE storec(a,b,n)

c	stores first n elements of single precision a into single precision b

c	Note: this DOUBLE PRECISION version of the routine has same name
c	as single precision version, but is distinguished by its file name

	IMPLICIT NONE
	
	INTEGER i,n
	
	REAL*8	a(*),b(*)

	DO i=1,n
	 b(i)=a(i)
	ENDDO
	
        RETURN
		
      END SUBROUTINE storec
