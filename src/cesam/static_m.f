
c********************************************************************

	SUBROUTINE static_m(fait,cx,li,y,be,ae,compt,dt,reprend,ip)

c	routine private du module mod_static	

c	Calcul des quantités caractéristiques du problème quasi-statique
c	en lagrangien
c	routine appellée par resout

c	spline collocation avec fonction de répartition

c	avec pression turbulente 7 inconnues 
c	sans pression turbulente 6 inconnues, Ptot=Pgaz

c	pour avoir un algorithme unique pour la diffusion, la composition
c	chimique est toujours tabulée en fonction de mu=(m/Msol)^2/3
c	que ce soit en lagrangien ou en eulérien

c	l'énergie graviphique, TdS/dt=tds est tabulée en fonction
c	de m^2/3 en lagrangien et de m en eulérien,

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c	08 10 96 : introduction de la rotation
c	08 01 96 : rectification de ctem*fac
c	05 05 97 : 2 possibilités de calcul du TdS
c	26 06 97 : remplacement du moment angulaire par vitesse angulaire
c	18 07 97 : correction d'énergie excédentaire de rotation
c	30 07 00 : introduction F95

c entrées
c	fait=1
c	-calcule les résidus be(i) des équations au point xcoll i=1,ne
c	fait=2
c	-calcule le "résidu" be(1) de la condition au point limite li

c	cx : indice du point de collocation
c	li : numéro de la limite
c	y : variables au point de collocation xcoll(cx)
c	compt: compteur du nb. iter. Newton Raphson
c	dt : pas temporel
c	ip : indice de couche pour détermination du facteur de répartition

c sorties
c	be : résidus
c	ae : termes du jacobien pour le point de collocation
c	reprend=.TRUE. : variation relative de U ou Ro trop forte

c	dans cette routine on est toujours en lagrangien
c	en_masse = .TRUE. variables lagrangiennes m23=m**23, r2=r**2

c----------------------------------------------------------------

	USE mod_donnees, ONLY : alpha, cpturb, ctel, ctem, ctep, cter, ctet,
	1 dtmin, d_grav, g, Ipg, Krot, kipp, langue, lsol, msol,
	2 m_ch, m_ptm, m_tds, nchim, ne, nrot, ord_rot, pi, pturb, ro_test,
	3 rsol, secon6, t_inf, nucleo
	USE mod_atm, ONLY : atm
	USE mod_etat, ONLY : etat
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn, no_croiss
	USE mod_variables, ONLY : bp, bp_t, chim, chim_gram, chim_t, inter,
	1 knotc, knotc_t, knotr, knot_tds_t, knot_ptm, knot_t, mc, mct,
	2 mct_t, mc_t, mrot, mrott, mstar, m23_t, n_ch, n_ch_t, n_ptm,
	3 n_qs, n_rot, n_qs_t, n_tds_t, old_ptm, qt_t, q_t, rota, rstar, r2_t,
	4 sortie, tds_t, wrot, xt_ptm, xt_tds_t, x_ptm, x_tds_t
		
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:,0:) :: y	
	REAL (kind=dp), INTENT(in) :: dt		
	INTEGER, INTENT(in) :: cx, compt, fait, li, ip
	
	REAL (kind=dp), INTENT(out), DIMENSION(:,:,0:) :: ae	
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: be	
	LOGICAL, INTENT(out) :: reprend
	
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: yd
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: xchim0
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: bes, dern
	REAL (kind=dp), DIMENSION(nchim) :: depsx, dxchim, dxchim_t, xchim,
	1 xchim_t	
	REAL (kind=dp), DIMENSION(ne) :: dfqs, fqs
	REAL (kind=dp), DIMENSION(nrot) :: dfrot, frot		 
	REAL (kind=dp), DIMENSION(5) :: epsilon
	REAL (kind=dp), DIMENSION(1) :: dfm, dftds, fm, ftds
	 	
	REAL (kind=dp), PARAMETER :: dd=1.d-5 	
	REAL (kind=dp), SAVE ::	cte1, cte2, cte6, cte11, cte13, cte14, cte14t,
	1 cte15, cte21, cte22, cte23, cte24, cte25, dmdl0, dmdr0, dpdl0, dpdr0,
	2 dptdl0, dptdr0, dtdl0, dtdr0, dt_tds, l13, mext0, m13,
	3 pext0, ptext0, ray, text0, unpdd
	
	REAL (kind=dp) :: acc_c, alfa, ay3, ay4, ay5, beta, cp,
	1 dcpm, dcpp, dcpt, dcpx, delta, deltam, deltap, deltat, deltax,
	2 depsm, depsp, depst, dgamlpt, dgamp,
	3 dgampt, dgamt, dgamx, dgamm, dgaml, dgamlp, dgamlpp, dgamr,
	4 dgam1, gradad, dgradadm, dgradadp, dgradadt, dgradadx,
	5 dgradlp, dgradlpt, dgradp, dgradpt, dgradt, dgradx,
	6 dgradm, dgradl, dgradr, dgradlpp, dhpp, dhppt, dhpt, dhpx,
	7 dhpr, dhpm, dkapp, dkapt, dkapx, dlpp, dlppdp, dlppdpt, dmdl,
	8 dmdr, dpn, dpdl, dpdr, dpsitl,
	9 dpsitm, dpsitp, dpsitr, dpsitt, dptdl, dptdr, dp_tm=0.d0, drom, drop, drop_t,
	1 drot, drot_t=0.d0, drox, drox_t=0.d0, dro_tm=0.d0, dstor, dtdl, dtdr,
	2 dtdsl=0.d0, dtdsm=0.d0, dtdsp=0.d0, dtdst=0.d0, dtdst_m=0.d0,
	3 dtetal, dtetam, dtetap, dtetar, dtetat, dtn, dt_tm=0.d0,
	4 dum, dup, dup_t=0.d0, dut, dut_t=0.d0, duv, dux, dux_t=0.d0,
	5 du_tm=0.d0, dvdm, dvdr, dwdm, dyv, dyvdm, dyvdr,
	6 dy0, dy1, dy2, dy3, dy4, dy5, dy71=0.d0, dy710, dy71l=0.d0,
	7 dy71lp=0.d0, dy71lpt=0.d0, dy71m=0.d0, dy71p=0.d0, dy71pt=0.d0,
	8 dy71r=0.d0, dy71t=0.d0, dy72=0.d0, dy720=0.d0, dy72l=0.d0,
	9 dy72lp=0.d0, dy72lpt=0.d0, dy72p=0.d0, dy72pt=0.d0,
	1 dy72m=0.d0, dy72r=0.d0, dy72t=0.d0, eps0, gam, gam0, gam1,
	2 gamma1, grad, grad0 ,gradrad, grad_mu, gravs2, hp, kap,
	3 ln, mext, mk, mn, pext, pgn, prn, psist, psist0, ptext,
	4 p_t=0.d0, p_t0, ro, ro0, ro_t=0.d0, ro_t0, stor, stor0,
	5 tdst=0, tds0, teff, teta, teta0, text, trn, t_t=0.d0, t_t0, u,
	6 u_t=0.d0, u_t0, u0, x00, v, w

	INTEGER, SAVE :: l
	INTEGER :: i, id, j, jd
	
	LOGICAL, SAVE :: deriv=.FALSE., der, init=.TRUE.
	LOGICAL :: radiatif 
	
	CHARACTER (len=7), SAVE, DIMENSION(7) :: variable=(/ 'ln Ptot',
	1 'ln T   ','R**2   ','L**2/3 ','M**23  ','Psi    ','ln Pgaz'/)
	
c------------------------------------------------------------------------

2000	FORMAT(8es10.3)

c	PRINT* ; PRINT*,'entree static_m',fait,compt,cx,ip
c	WRITE(*,2000)dt,xcoll(cx) ; WRITE(*,2000)y(1:ne,0),fac(ip)
c	WRITE(*,2000)y(1:ne,1)	
	
	IF(init)THEN
	 init=.FALSE.
	 cte1=g*msol/rsol**2
	 cte2=2.d0/3.d0*rsol
	 cte11=-3.d0*g/8.d0/pi*(msol/rsol**2)**2	!pour d ln P /d q
	 cte21=cte11*ctep ; cte22=cte11*ctet
	 cte13=msol/rsol*3.d0/4.d0/pi/rsol/rsol		!pour d zeta /d q
	 cte23=cte13*cter
	 cte14=msol/lsol				!pour d lambda/d q
	 cte24=cte14*ctel
	 cte14t=cte14/secon6				! "    " dt /= 0
	 cte15=msol/rsol/4.d0/pi			!pour  la rotation
	 cte25=cte15*ctep
	 cte6=ABS(cpturb)*alpha**2/8.d0		!pour la Pturb
	 der=cpturb < 0.d0	!der=.TRUE. on tient compte de dlnPgaz/dlnPtot
	 dt_tds=10.d0*dtmin	!dt < dt_tds on utilise la tabulation TdS

	 ALLOCATE(xchim0(nchim))	!pour l'atmosphère

	 IF(kipp)THEN
	  WRITE(2,10) ; WRITE(*,10)	  
10	  FORMAT(/,'Approximation de Kippenhahn pour le calcul du TdS',/)
	 ELSE
	  WRITE(*,11) ; WRITE(2,11)	  	  
11	  FORMAT(/,'Calcul exact du TdS')
	 ENDIF	 
	ENDIF		!init
	
c	PRINT*,'fait,cx',fait,cx

c------------------------------------------------------------------

c ensemble des variables

c		y(1,0)=ln Ptot		y(1,1)=(ln Ptot)'
c		y(2,0)=ln T		y(2,1)=(ln T)'	
c		y(3,0)=r**2		y(3,1)=(r**2)'
c		y(4,0)=l**2/3		y(4,1)=(l**2/3)'
c		y(5,0)=m**2/3		y(5,1)=(m**2/3)'
c		y(6,0)=psi=dQ/dq(=cte)	y(6,1)=psi'(=0)
c avec Pturb	y(Ipg,0)=ln Pgaz	y(Ipg,1)=(ln Pgaz)'

c	xcoll(cx) valeur de q au point de collocation q=1, 2, 3, .. , n_qs

c ae(eq,var,der)=ae(ne,ne,0:1)=dérivée de la eq-ieme équation par
c rapport à la der-ieme dérivée de la var-ieme variable

	ae=0.d0 ; be=0.d0 ; reprend=.FALSE.
	SELECT CASE(fait)	
	CASE(1)		!le point courant
	
	 prn=EXP(y(1,0))		!pression Ptot cgs
	 IF(pturb)THEN		!avec pression turbulente
	  pgn=EXP(y(Ipg,0))		!pression Pgaz cgs	 
	 ELSE			!sans pression turbulente
	  pgn=prn
	 ENDIF
	 trn=EXP(y(2,0))	!temperature K
	 ay3=ABS(y(3,0))	!r**2
	 ray=SQRT(ay3)		!rayon
	 ay4=ABS(y(4,0))	!l**2/3
	 l13=SQRT(ay4)		!l**1/3
	 ln=l13**3		!l/lsol
	 ln=SIGN(ln,y(4,0))	!l/lsol	 
	 ay5=ABS(y(5,0))	!m**2/3
	 m13=SQRT(ay5)		!m**1/3
	 mn=m13**3		!masse
c	 PRINT*,'mn=',mn

c l'acc. centrifuge ne doit pas excèder 90% gravité
	 gravs2=cte1*mn/ay3*0.9d0

c vitesse angulaire
	 SELECT CASE(Krot)
	 CASE(0,1,2)
	  w=wrot
	 CASE(3,4,5)
	  CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,
	1 MAX(mrot(1),MIN(ay5,mrot(n_rot))),l,frot,dfrot)	
	  IF(no_croiss)PRINT*,'Pb. en 0 dans static_m'
	  w=frot(1)
	 END SELECT
	 
c accélération centrifuge	 
	 acc_c=cte2*ray*w**2 ; reprend=acc_c > gravs2
	 IF(reprend)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1012)acc_c,gravs2,w,ray,mn ; WRITE(2,1012)acc_c,gravs2,w,ray,mn
1012	   FORMAT('STOP, in static_m, the centrifugal acceleration =',
	1  es10.3,' > 90% gravity = ',es10.3,/,'angular velicity =',es10.3,
	2  ', R/Rsun=',es10.3,', M/Msun=',es10.3)
	  CASE DEFAULT
	   WRITE(*,12)acc_c,gravs2,w,ray,mn ; WRITE(2,12)acc_c,gravs2,w,ray,mn
12	   FORMAT('ARRET, dans static_m, acc. centifuge =',
	1  es10.3,' >  90% gravité = ',es10.3,/,'vitesse angulaire =',es10.3,
	2  ', R/Rsol=',es10.3,', M/Msol=',es10.3)	  
	  END SELECT
	  RETURN		!CALL sortie
	 ENDIF

c si, à cause d'une erreur d'arrondi, ce qui peut arriver au voisinage
c du centre r**2, l**2/3 ou m**23 est negatif, on tente de forcer la
c convergence en utilisant les ABS(**), les dérivées sont inexactes et
c la physique est violée... mais si ça passe les erreurs d'arrondi
c disparaissent, parfois, au cours des iterations
c pour eviter cette disposition supprimer les "c" devant reprend et RETURN
	 IF(MIN(y(3,0),y(4,0),y(5,0)) < 0.d0)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1014)y(3:5,0),SQRT(ABS(y(3,0))),
	1  SQRT(ABS(y(4,0)))**3,SQRT(ABS(y(5,0)))**3,NINT(xcoll(cx))
1014	   FORMAT(/,'static_m, r2=',es10.3,', l23=',es10.3,', m23=',es10.3,/,
	1  'R=',es10.3,', L=',es10.3,', M=',es10.3,/,
	2  'shell=',i4,', one try to converge')	  
	  CASE DEFAULT	  
	   WRITE(*,14)y(3:5,0),SQRT(ABS(y(3,0))),
	1  SQRT(ABS(y(4,0)))**3,SQRT(ABS(y(5,0)))**3,NINT(xcoll(cx))
14	   FORMAT(/,'static_m, r2=',es10.3,', l23=',es10.3,', m23=',es10.3,/,
	1  'R=',es10.3,', L=',es10.3,', M=',es10.3,/,
	2  'couche=',i4,', on tente de converger')	  
	  END SELECT
	 ENDIF
	
c composition chimique au temps t+dt
	 CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1 knotc,.TRUE.,MIN(ay5,mc(n_ch)),l,xchim,dxchim)
	 IF(no_croiss)PRINT*,'Pb. at 1 in static_m'
	 xchim=ABS(xchim)
	 
c	 IF(SUM(xchim*nucleo) > 1.001d0)THEN	
c	  WRITE(*,2000)SUM(xchim*nucleo) ; WRITE(*,2000)xchim ; PAUSE'static_m'
c	 ENDIF	
c	 PRINT*,'cx,nchim/prn,pgn,trn,mn,ray / xchim / y  /dy',cx,nchim
c	 WRITE(*,2000)prn,pgn,trn,mn,ray,xchim(1),dxchim(1)
c	 WRITE(*,2000)xchim ; WRITE(*,2000)dxchim
c	 WRITE(*,2000)y(1:ne,0) ; WRITE(*,2000)y(1:ne,1)
c	 PRINT*,' '
	
	 IF(pturb .AND. der)THEN  !avec pression turbulente 7 inconnues
	  dlpp=y(Ipg,1)/y(1,1)	!dlpp=dln Pgaz/dln Ptot
	  dlppdpt=-dlpp/y(1,1)	!dérivée dlpp /dln Ptot 
	  dlppdp=dlpp/y(Ipg,1)	!dérivée dlpp /dln Pgaz
	 ELSE	!der=.FALSE. on ne tient pas compte de dln Pgaz/dln Ptot
	  dlpp=1.d0
	  dlppdpt=0.d0 ; dlppdp=0.d0
	 ENDIF
	 CALL thermo(prn,pgn,trn,mn,ln,ray,dlpp,xchim,dxchim,
	1 ro,drop,drot,drox,u,dup,dut,dux,
	2 grad,dgradpt,dgradp,dgradt,dgradx,dgradm,dgradl,dgradr,dgradlpp,
	3 gam,dgampt,dgamp,dgamt,dgamx,dgamm,dgaml,dgamr,dgamlpp,		
	4 epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
	5 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	6 gradad,dgradadp,dgradadt,dgradadx,
	7 hp,dhppt,dhpp,dhpt,dhpx,dhpr,dhpm,grad_mu,
	8 gradrad,alfa,beta,gamma1,radiatif)
	 
c	 WRITE(*,2000)mn,prn,pgn,trn,ro,xchim,epsilon(1),depsp,depst,depsx
c	 PRINT*,'couche',cx ; WRITE(*,2000)prn,pgn,trn,mn,ln,ray,ro,grad
c	 WRITE(*,2000)gam,epsilon(1),delta,gradad,hp,gradrad
c	 WRITE(*,2000)xchim ; WRITE(*,2000)epsilon
c	 WRITE(*,2000)prn,pgn,dgradpt,dgradp ; PAUSE'thermo'

c dérivées par rapport à : ln Ptot, ln T, zeta, lambda, mu, ln Pgaz
	 drop=drop*pgn	!der ro /ln Pgaz
	 drot=drot*trn ; drom=drox*dxchim(1)
	 dcpp=dcpp*pgn ; dcpt=dcpt*trn ; dcpm=dcpx*dxchim(1)
	 deltap=deltap*pgn ; deltat=deltat*trn ; deltam=deltax*dxchim(1)
	 dup=dup*pgn ; dut=dut*trn ; dum=dux*dxchim(1)

	 dgradadp=dgradadp*pgn/gradad	!der gradad/ln Pgaz gradad près
	 dgradadt=dgradadt*trn/gradad	!der gradad/ln T à gradad près
	 dgradadm=dgradadx*dxchim(1)/gradad  !dgradad/m**2/3 gradad près

	 dgradpt=dgradpt*prn/grad	!der grad /ln Ptot à grad près
	 dgradp=dgradp*pgn/grad		!der grad /ln Pgaz à grad près
	 dgradt=dgradt*trn/grad		!der grad /ln T
	 dgradr=dgradr/2.d0/ray/grad		!der grad /r**2
	 dgradl=dgradl*3.d0/2.d0*l13/grad	!der grad /l**2/3
	 dgradm=(dgradm*3.d0/2.d0*m13+dgradx*dxchim(1))/grad !d grad/m**2/3
	 dgradlpt=dgradlpp*dlppdpt/grad		!der grad /dln Ptot
	 dgradlp =dgradlpp*dlppdp/grad		!der grad /dln Pgaz

	 dgampt=dgampt*prn				!der gam /ln Ptot
	 dgamp= dgamp*pgn				!der gam /ln Pgaz
	 dgamt= dgamt*trn				!der gam /ln T
	 dgamr=dgamr/2.d0/ray				!der gam /r**2
	 dgaml=dgaml*3.d0/2.d0*l13			!der gam /l**2/3
	 dgamm=dgamm*3.d0/2.d0*m13+dgamx*dxchim(1)	!der gam /m**2/3
	 dgamlpt=dgamlpp*dlppdpt			!der gam /dln Ptot
	 dgamlp =dgamlpp*dlppdp			!der gam /dln Pgaz

	 dhppt=dhppt*prn/hp	!der hp /ln Ptot à hp près
	 dhpp=dhpp*pgn/hp	!der hp /ln Pgaz à hp près
	 dhpt=dhpt*trn/hp	!der hp /ln T à hp près
	 dhpr=dhpr/2.d0/ray/hp	!der hp /r**2 à hp près
	 dhpm=(dhpm*3.d0/2.d0*m13+dhpx*dxchim(1))/hp !der hp/m**2/3 hp près
	 IF(epsilon(1) > 1.d-25)THEN
	  depsp=depsp*pgn/epsilon(1)	!dérivée /ln Pgaz, à 1/epsilon près
	  depst=depst*trn/epsilon(1)
	  depsm=SUM(depsx(1:nchim)*dxchim(1:nchim))/epsilon(1)
	 ELSE
	  depsp=0.d0 ; depst=0.d0 ; depsm=0.d0
	 ENDIF
	 
c sans pression turbulente Ptot=Pgaz
c	IF(.NOT.pturb)THEN
c	 dgradpt=dgradp ; dgradlpt=dgradlp ; dgampt=dgamp
c	 dgamlpt=dgamlp ; dhppt=dhpp
c	ENDIF
	
c on ne tient pas compte de la différence Pgaz, Ptot pour
c l'énergie graviphique TdS=tds et old_ptm sont tabules en fonction
c de m^2/3 en lagrangien et de m en eulérien,
c pour avoir un algorithme unique pour la diffusion, la composition
c chimique est toujours tabulée en fonction de m^2/3
	 IF(dt /= 0.d0)THEN
	  IF(compt == 0 .OR. dt <= dt_tds)THEN	!interp. du TdS au temps t
	   CALL bsp1dn(1,tds_t,x_tds_t,xt_tds_t,n_tds_t,m_tds,
	1  knot_tds_t,.TRUE.,MIN(ay5,x_tds_t(n_tds_t)),l,ftds,dftds)
	   IF(no_croiss)PRINT*,'Pb. at 2 in static_m'
	   tdst=ftds(1) ; dtdst_m=dftds(1)
c	   WRITE(*,2000)tdst,y(5,0),m23_t(n_qs_t)
c	   PAUSE'TdS interpolé'

	  ELSE				!TdS au temps t+dt
	   CALL bsp1dn(1,old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm,.TRUE.,
	1  MIN(ay5,x_ptm(n_ptm)),l,fm,dfm)
	   IF(no_croiss)PRINT*,'Pb. at 3 in static_m'	 	
	   mk=fm(1)		 
	   CALL inter('m23',bp_t,q_t,qt_t,n_qs_t,knot_t,mk,fqs,dfqs,r2_t,m23_t)
	   IF(pturb)THEN		!avec pression turbulente 7 inconnues
	    p_t=EXP(fqs(Ipg)) ; dp_tm=dfqs(Ipg)*p_t	!d p_t / d mu au temps t
	   ELSE			!sans pression turbulente 6 inconnues
	    p_t=EXP(fqs(1)) ; dp_tm=dfqs(1)*p_t	!d p_t / d mu au temps t
	   ENDIF
	   t_t=EXP(fqs(2)) ; dt_tm=dfqs(2)*t_t	!d t_t / d mu au temps t

	   IF(kipp)THEN	!approximation de Kippenhahan
	    dpn=pgn-p_t ; dtn=trn-t_t
	    IF(trn > t_inf)THEN	!contrôle de la variation   
c	     duv=1.d0-ABS((cp*t_t-delta/ro*p_t)/(cp*trn-delta/ro*prn))    
	     duv=MAX(ABS(dpn/pgn),ABS(dtn/trn))
	    ELSE
	     duv=0.d0
	    ENDIF	   
	    reprend=(duv > d_grav) .AND. (compt > 3)
	    IF(reprend)THEN		!TdS varie trop
	     WRITE(usl_static,7)mn/mstar,ray/rstar
7	     FORMAT(/,'Approx. Kipp. trop de variation de TdS en m=',es10.3,
	1    ', R/R*=',es10.3)
	     WRITE(usl_static,13)duv,d_grav
13	     FORMAT('duv=',es10.3,' > d_grav=',es10.3)
	     WRITE(usl_static,8)p_t,pgn,dpn/pgn
8	     FORMAT('P(t)=',es10.3,', P(t+dt)=',es10.3,', dP/P=',es10.3)	      
	     WRITE(usl_static,9)t_t,trn,dtn/trn
9	     FORMAT('T(t)=',es10.3,', T(t+dt)=',es10.3,', dT/T=',es10.3,/)	     
	     RETURN
	    ELSE
	     duv=0.d0
	    ENDIF	!reprend	   
	    tdst=(cp*dtn-delta/ro*dpn)/dt
	   
	   ELSE		!TdS=dU+PdV	  
	    CALL bsp1dn(nchim,chim_t,mc_t,mct_t,n_ch_t,m_ch,
	1   knotc_t,.TRUE.,MIN(mk,mc_t(n_ch_t)),l,xchim_t,dxchim_t)
	    IF(no_croiss)PRINT*,'Pb. at 4 in static_m'
	    CALL chim_gram(xchim_t,dxchim_t)	  
c	    WRITE(*,2000)p_t,t_t,xchim_t(1),xcoll(cx)
	    CALL etat(p_t,t_t,xchim_t,.TRUE.,	!mettre .TRUE. pour d/dX
	1   ro_t,drop_t,drot_t,drox_t,u_t,dup_t,dut_t,dux_t,
	2   delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3   gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	    IF(ro > ro_test .AND. trn > t_inf)THEN    
c	     duv=1.d0-ABS((u_t-prn/ro**2*ro_t)/(u-prn/ro))	    
	     duv=MAX(ABS(u-u_t)/u,ABS(ro-ro_t)/ro)
	    ELSE
	     duv=0.d0
	    ENDIF
	    reprend=(duv > d_grav) .AND. (compt >= 3)
	    IF(reprend)THEN		!TdS varie trop
	     SELECT CASE(langue)
	     CASE('english')
	      WRITE(usl_static,1006)mn/mstar,ray/rstar
1006	      FORMAT(/,'Too large change of dU+PdV at M/M*=',es10.3,', R/R*=',
	1     es10.3)
	     CASE DEFAULT
	      WRITE(usl_static,6)mn/mstar,ray/rstar
6	      FORMAT(/,'Trop de variation de dU+PdV en M/M*=',es10.3', R/R*=',
	1     es10.3)
	     END SELECT
	     WRITE(usl_static,13)duv,d_grav
	     WRITE(usl_static,1)u_t,u,ABS(u-u_t)/u
1	     FORMAT('U(t)=',es10.3,', U(t+dt)=',es10.3,', dU/U=',es10.3)
 	     WRITE(usl_static,2)ro_t,ro,ABS(ro-ro_t)/ro
2	     FORMAT('ro(t)=',es10.3,', ro(t+dt)=',es10.3,', dro/ro=',es10.3)
 	     WRITE(usl_static,3)xchim_t(1:MIN(7,n_ch))
3	     FORMAT('X_t(i=1,7)=  ',7es10.3)
 	     WRITE(usl_static,15)xchim(1:MIN(7,n_ch))
15	     FORMAT('X_t+1(i=1,7)=',7es10.3)
	     WRITE(usl_static,4)p_t,pgn
4	     FORMAT('P(t)=',es10.3,', P(t+dt)=',es10.3)	     
	     WRITE(usl_static,5)t_t,trn
5	     FORMAT('T(t)=',es10.3,', T(t+dt)=',es10.3,/)	          
	     RETURN
	    ENDIF	!reprend
	
	    du_tm =dup_t *dp_tm+dut_t *dt_tm+dux_t *dxchim_t(1)
	    dro_tm=drop_t*dp_tm+drot_t*dt_tm+drox_t*dxchim_t(1)
	    tdst=(u-u_t-prn/ro**2*(ro-ro_t))/dt
	   
c	    PRINT*,'cx/Ptot,Pgaz,T,p_t,t_t,u,u_t,ro,ro_t,tdst',cx
c	    WRITE(*,2000)prn,trn,p_t,t_t,u,u_t,ro,ro_t,tdst

	   ENDIF	!kipp
	  ENDIF		!compt=0
	 ENDIF		!dt /= 0

c la rotation, w ne dépend que de m
c v correspond à l'accélération centrifuge w^2/6 pi R
	 SELECT CASE(Krot)
	 CASE(0,1,2)
	  w=wrot ; dwdm=0.d0
	 CASE(3,4,5)
	  CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,
	1 MAX(mrot(1),MIN(ay5,mrot(n_rot))),l,frot,dfrot)	
	  IF(no_croiss)PRINT*,'Pb. at 5 in static_m'
	  w=frot(1) ; dwdm=dfrot(1)
	 END SELECT
	 IF(w == 0.d0)THEN
	  v=0.d0 ; dvdr=0.d0 ; dvdm=0.d0 ; dyv=0.d0 ; dyvdr=0.d0
	  dyvdm=0.d0
	 ELSE
	  v=SQRT(ay5/ay3)/prn*w**2
	  dvdr=-v/2.d0/ay3 ; dvdm= v/2.d0/ay5+2.d0*v/w*dwdm	
	  dyv=cte25*fac(ip)*v		!dérivée/lnP comme celle de dy1:
	  dyvdr=cte25*fac(ip)*dvdr	!dyv/dlnp=-dyv
	  dyvdm=cte25*fac(ip)*dvdm	  
	 ENDIF	 
	
c la fonction de répartition sur ln Ptot, ln T, l^2/3 et m^2/3

c teta = ctep dksi/dmu
	 dy1=fac(ip)*(ay5/ay3)**2/prn
	 teta=cte21*dy1 
	 dtetap=-teta			!dérivée/ln Ptot
	 dtetar=-teta*2.d0/ay3		!dérivée/r**2
	 dtetam=teta*2.d0/ay5		!dérivée/m**2/3
	 
c teta = teta + rotation
	 teta=teta+dyv
 	 dtetap=dtetap-dyv
	 dtetar=dtetar+dyvdr
	 dtetam=dtetam+dyvdm
	 
c teta = teta + ctet dksi/dmu grad
	 dy2=cte22*dy1*grad ; teta=teta+dy2
	 dtetap=dtetap+dy2*(dgradpt-1.d0)
	 dtetat=dy2*dgradt
	 dtetar=dtetar+dy2*(-2.d0/ay3+dgradr)
	 dtetal=dy2*dgradl
	 dtetam=dtetam+dy2*(2.d0/ay5+dgradm)
	 
c teta = teta + cter d zeta/dmu
	 dy3=cte23*fac(ip)/ro*m13/ray ; teta=teta+dy3
	 dtetap=dtetap-dy3*drop/ro	 
	 dtetat=dtetat-dy3*drot/ro	  
	 dtetar=dtetar-dy3/ay3/2.d0
	 dtetam=dtetam+dy3/ay5/2.d0	 
	 	 	 
c teta = teta + ctel d lambda / d mu	 
	 dy4=cte24*fac(ip)*m13/l13*epsilon(1) ; teta=teta+dy4	 
	 dtetap=dtetap+dy4*depsp			!dérivée/ln P
	 dtetat=dtetat+dy4*depst			!dérivée/ln T
	 dtetal=dtetal-dy4/2.d0/ay4			!dérivée/l**2/3
	 dtetam=dtetam+dy4*(0.5d0/ay5+depsm)		!dérivée/m**2/3	 
	 
c teta = teta + ctem d mu/ d mu 
	 teta=teta+fac(ip)*ctem
	 	 
c psi/teta	 
	 psist=y(6,0)/teta
	 dpsitp=-psist/teta*dtetap	!d (psi/teta) / d ln Ptot
	 dpsitt=-psist/teta*dtetat	!d (psi/teta) / d ln T	  
	 dpsitr=-psist/teta*dtetar	!d (psi/teta) / d r**2
	 dpsitl=-psist/teta*dtetal	!d (psi/teta) / d l**2/3	  
	 dpsitm=-psist/teta*dtetam	!d (psi/teta) / d m**2/3

c équations (Prad est dans Pgas), avec Pturb, y(1,0) est Pgas+Pturb
c sans Pturb, y(1,0) est Pgas

c dérivée/ln P comme celle de dy1: dyv/dlnp=-dyv
	 dyv=cte15*v ; dyvdr=cte15*dvdr ; dyvdm=cte15*dvdm
	 dy1=cte11*(ay5/ay3)**2/prn ; dy0=dy1+dyv
	 be(1)=y(1,1)-dy0*psist	!dln Ptot/dmu=3g/8 pi (mu+1/zeta+1)**2/P
	 ae(1,1,0)=dy0*(psist-dpsitp)		!der./ln Ptot
	 ae(1,1,1)=1.d0				!der./dln Ptot
	 ae(1,2,0)=-dy0*dpsitt			!der./dln T	   
	 ae(1,3,0)=-(dyvdr-dy1*2.d0/ay3)*psist-dy0*dpsitr !dérivée/r**2
	 ae(1,4,0)=-dy0*dpsitl				!der. /l**2/3  
	 ae(1,5,0)=-(dyvdm+dy1*2.d0/ay5)*psist-dy0*dpsitm !der. /m**2/3
	 ae(1,6,0)=-dy0/teta				!dérivée/psi

c dln T/dq=dln Ptot/dq grad	
	 dy2=dy0*psist*grad ; be(2)=y(2,1)-dy2
	 ae(2,1,0)=ae(1,1,0)*grad-dy2*dgradpt	!dérivée/ln Ptot 
	 ae(2,2,0)=ae(1,2,0)*grad-dy2*dgradt	!dérivée/ln T
	 ae(2,2,1)=1.d0				!dérivée/dln T
	 ae(2,3,0)=ae(1,3,0)*grad-dy2*dgradr	!dérivée/r**2
	 ae(2,4,0)=ae(1,4,0)*grad-dy2*dgradl	!dérivée/l**2/3
	 ae(2,5,0)=ae(1,5,0)*grad-dy2*dgradm	!dérivée/m**2/3
	 ae(2,6,0)=ae(1,6,0)*grad		!dérivée/psi		 
	 IF(pturb)THEN				!avec pression turb. 7 inc.
	  ae(2,1,1)=-dy2*dgradlpt		!dérivée/dln Ptot
	  ae(2,Ipg,0)=-dy2*dgradp		!dérivée/ln Pgaz
	  ae(2,Ipg,1)=-dy2*dgradlp		!dérivée/dln Pgaz
	 ENDIF

c dr**2/dmu=-3/4pi /ro SQRT(mu/zeta)	
	 dy3=cte13/ro*m13/ray; be(3)=y(3,1)-dy3*psist
	 IF(pturb)THEN
	  ae(3,1,0)=-dy3*dpsitp			!dérivée/ln Ptot
	  ae(3,Ipg,0)= dy3*drop/ro*psist	!der. /ln Pgaz
	 ELSE					!sans press. turb. 6 inc.
	  ae(3,1,0)=dy3*(drop/ro*psist-dpsitp)!der./ln Ptot 
	 ENDIF
	 ae(3,2,0)=dy3*(psist*drot/ro-dpsitt)	!dérivée/ln T
	 ae(3,3,0)=dy3*(psist/ay3/2.d0-dpsitr)	!dérivée/r**2
	 ae(3,3,1)=1.d0				!dérivée/dr**2
	 ae(3,4,0)=-dy3*dpsitl 			!dérivée/l**2/3  
	 ae(3,5,0)=dy3*(psist*(-0.5d0/ay5+drom/ro)-dpsitm) !dérivée/m**2/3
	 ae(3,6,0)=-dy3/teta				!dérivée/psi

c dl/dmu=epsilon m13/l13 +Tds/dt contribution de epsilon	
	 dy4=cte14*epsilon(1)*m13/l13
	 be(4)=y(4,1)-dy4*psist
	 IF(pturb)THEN		!avec pression turbulente 7 inconnues
	  ae(4,1,0)=-dy4*dpsitp				!dérivée/ln Ptot	 
	  ae(4,Ipg,0)=-dy4*psist*depsp			!dérivée/ln Pgaz
	 ELSE			!sans pression turbulente 6 inconnues
	  ae(4,1,0)=-dy4*(psist*depsp+dpsitp)		!dérivée/ln Ptot
	 ENDIF
	 ae(4,2,0)=-dy4*(psist*depst+dpsitt)		!dérivée/ln T
	 ae(4,3,0)=-dy4*dpsitr				!dérivée/r**2
	 ae(4,4,0)= dy4*(psist/2.d0/ay4-dpsitl)	!der/l**2/3
	 ae(4,4,1)=1.d0				!dérivée/dl**2/3
	 ae(4,5,0)=-dy4*(psist*(0.5d0/ay5+depsm)+dpsitm)!der/m**2/3	
	 ae(4,6,0)=-dy4/teta				!dérivée/psi

c contribution du Tds/dt	
	 IF(dt /= 0.d0)THEN
	  dy5=cte14t*m13/l13 ; tdst=dy5*tdst
	  be(4)=be(4)+tdst*psist		!epsilon-(tdS+d Iw2)/dt 
	  dtdsl=-tdst/2.d0/ay4				!dérivée/l**2/3
	  dtdsm= tdst/2.d0/ay5				!dérivée/m**2/3
	  IF(compt == 0 .OR. dt <= dt_tds)THEN
	   dtdsm= dtdsm+dy5*dtdst_m			!dérivée/m**2/3
	   ae(4,1,0)=ae(4,1,0)+tdst*dpsitp		!dérivée/ln Ptot
	   ae(4,2,0)=ae(4,2,0)+tdst*dpsitt		!dérivée/ln T
	   ae(4,3,0)=ae(4,3,0)+tdst*dpsitr		!dérivée/r**2
	   ae(4,4,0)=ae(4,4,0)+tdst*dpsitl+psist*dtdsl	!dérivée/l**2/3
	   ae(4,5,0)=ae(4,5,0)+tdst*dpsitm+psist*dtdsm	!dérivée/m**2/3
	   ae(4,6,0)=ae(4,6,0)+tdst/teta		!dérivée/psi
	  ELSE
	   IF(kipp)THEN
	    dtdsp=dy5*(dcpp*dtn-delta/ro*((deltap/delta
	1   -drop/ro)*dpn+pgn))/dt			!dérivée/ln Pgaz    
	    dtdst=dy5*(dcpt*dtn+cp*trn
	1   -delta/ro*dpn*(deltat/delta-drot/ro))/dt	!dérivée/ln t
	    dtdsm= dtdsm+dy5*(dcpm*dtn-cp*dt_tm 	!dérivée/m**2/3
	1   -delta/ro*((deltam/delta-drom/ro)*dpn-dp_tm))/dt
	   ELSE
	    dtdsp=dy5*(dup-pgn/ro**2*((ro-ro_t)*
	1   (1.d0-2.d0*drop/ro)+drop))/dt 	!dérivée/ln Pgaz
	    dtdst=dy5*(dut-pgn/ro**2*drot*
	1   (1.d0-2.d0*(ro-ro_t)/ro))/dt	!dérivée/ln t
	    dtdsm= dtdsm+dy5*(dum-du_tm-pgn/ro**2*(-2.d0*drom/ro*(ro-ro_t)
	1   +drom-dro_tm))/dt			!dérivée/m**2/3
	   ENDIF
	   IF(pturb)THEN		!avec pression turbulente 7 inconnues
	    ae(4,1,0)=ae(4,1,0)+tdst*dpsitp		!dérivée/ln Ptot
	    ae(4,Ipg,0)=ae(4,Ipg,0)+psist*dtdsp		!dérivée/ln Pgaz
	   ELSE			!sans pression turbulente 6 inconnues
	    ae(4,1,0)=ae(4,1,0)+tdst*dpsitp+psist*dtdsp
	   ENDIF
	   ae(4,2,0)=ae(4,2,0)+psist*dtdst+tdst*dpsitt	!dérivée/ln T
	   ae(4,3,0)=ae(4,3,0)+tdst*dpsitr		!dérivée/r**2
	   ae(4,4,0)=ae(4,4,0)+tdst*dpsitl+psist*dtdsl	!dérivée/l**2/3
	   ae(4,5,0)=ae(4,5,0)+tdst*dpsitm+psist*dtdsm	!dérivée/m**2/3
	   ae(4,6,0)=ae(4,6,0)+tdst/teta		!dérivée/psi
	  ENDIF	!compt = 0 
	 ENDIF		!dt

c dmu/dq=psi/teta
	 be(5)=y(5,1)-psist
	 ae(5,1,0)=-dpsitp		!dérivée/ln Ptot
	 ae(5,2,0)=-dpsitt		!dérivée/ln T	  
	 ae(5,3,0)=-dpsitr		!dérivée/r**2
	 ae(5,4,0)=-dpsitl		!dérivée/l**2/3  
	 ae(5,5,0)=-dpsitm		!dérivée/m**2/3
	 ae(5,5,1)= 1.d0		!dérivée/dm**2/3)	 
	 ae(5,6,0)=-1.d0/teta		!dérivée/psi

c dpsi/dq=0
	 be(6)=y(6,1) ; ae(6,6,1)=1.d0	!dérivée/dpsi 
	 
c Avec Pturb:
c pour les points de collocation dans les zones radiatives,
c on resout d ln Ptot - d ln Pgas =0
c pour les points de collocation dans les zones convectives,
c on resout -Pturb -Pgaz +Ptot=0

c	 IF(radiatif)WRITE(*,2000)trn,prn,pgn,gradrad,grad,gradad

	 IF(pturb)THEN		!avec pression turbulente 7 inconnues
          IF(radiatif)THEN
           be(Ipg)=y(1,1)-y(Ipg,1)		!d lnPtot = d lnPgaz
           ae(Ipg,1,1)= 1.d0			!dérivée /ln Ptot
           ae(Ipg,Ipg,1)=-1.d0			!dérivée /ln Pgaz

          ELSE
           gam1=gam/(gam+1.d0)
           dgam1=1.d0/(gam+1.d0)**2/gam1

           dy71=cte6*gam1*delta*prn*gradad*dlpp
           dy71pt=dy71*(dgampt*dgam1+1.d0)			!der /ln Ptot
           dy71p= dy71*(dgamp*dgam1+deltap/delta+dgradadp)	!der /ln Pgaz
           dy71t= dy71*(dgamt*dgam1+deltat/delta+dgradadt) 	!der /ln T
           dy71r= dy71*dgamr*dgam1 				!der /r**2
           dy71l= dy71*dgaml*dgam1				!der /l**2/3
           dy71m= dy71*(dgamm*dgam1+deltam/delta+dgradadm)	!der /m**2/3
           dy71lpt=dy71*(dgamlpt*dgam1+dlppdpt/dlpp)  		!der/dln Ptot           
           dy71lp =dy71*(dgamlp *dgam1+dlppdp /dlpp)		!der /dln Pgaz
          
           dy72=cte6*gam1*delta*prn*grad
	   dy72pt=dy72*(dgampt*dgam1+1.d0+dgradpt)		!der /ln Ptot
           dy72p= dy72*(dgamp*dgam1+deltap/delta+dgradp)	!der /ln Pgaz
           dy72t= dy72*(dgamt*dgam1+deltat/delta+dgradt)	!der /ln T
           dy72r= dy72*(dgamr*dgam1+dgradr)			!der /r**2
           dy72l= dy72*(dgaml*dgam1+dgradl)			!der /l**2/3        
           dy72m= dy72*(dgamm*dgam1+deltam/delta+dgradm)	!der /m**2/3
	   dy72lpt=dy72*(dgamlpt*dgam1+dgradlpt)		!der/dln Ptot
	   dy72lp= dy72*(dgamlp *dgam1+dgradlp)			!der/dln Pgaz

	   be(Ipg)=dy71-dy72-pgn+prn 		!-Pturb -Pgaz +Ptot=0
c          WRITE(*,2000)be(Ipg),dy71,dy72,prn,pgn,dy71-dy72,prn-pgn

           ae(Ipg,1,0)=dy71pt-dy72pt+prn	!der /ln Ptot
           ae(Ipg,2,0)=dy71t -dy72t		!der /ln T
           ae(Ipg,3,0)=dy71r -dy72r		!der /r**2
           ae(Ipg,4,0)=dy71l -dy72l		!der /l**2/3         
           ae(Ipg,5,0)=dy71m -dy72m		!der /m**2/3
           ae(Ipg,Ipg,0)=dy71p -dy72p-pgn	!der /ln Pgaz         
           ae(Ipg,1,1)=dy71lpt-dy72lpt		!der /dln Ptot
           ae(Ipg,Ipg,1)=dy71lp -dy72lp		!der /dln Pgaz
          ENDIF	!radiatif
         ENDIF	!Pturb
         
c        PRINT*,'be' ; WRITE(*,2000)be
c        WRITE(*,2000)cte13,ro,ay5,ay3,dy3,psist,y(3,1)
c        WRITE(*,2000)cte14,epsilon(1),m13,l13,dy4,psist,y(4,1)      
c        PRINT*,'ae'
c        DO i=1,ne
c         WRITE(*,2000)ae(i,1:ne,0)
c        ENDDO      
c        DO i=1,ne
c         WRITE(*,2000)ae(i,1:ne,1)
c        ENDDO
c        PAUSE'fin jacob'

c---- test de vérification du jacobien par dérivées numériques-----------

c	 deriv=.TRUE.
c	 deriv=cx <= 2
c	 deriv=cx <= 3
c	 deriv=cx == 1
c	 deriv=cx == 90
c	deriv=cx >= 867 .AND. cx <= 870
c	deriv=cx >= 925 .AND. cx <= 930
c	deriv=cx >= 265 .AND. cx <= 267
c	1	.OR. cx == 50
c	2	.OR. cx == 90
c	2	.OR. cx == 189
c	2	.OR. cx == 190
c	2	.OR. cx >= 188 .AND. cx <= 190
c	2	.OR. cx >= 199 .AND. cx <= 202		
c	3	.OR. (cx > 74 .AND. cx < 76)
c	3	.OR. (cx > 43 .AND. cx < 48)
c	4	.OR. (cx > 88 .AND. cx < 92)
c	5	.OR. cx == 220
c	6	.OR. cx == 103
c	7	.OR. cx. ge. 148
c	8	.OR. cx. ge. 149
c	9	.OR. (cx. ge. 148 .AND. cx <150)
c	1	.OR. cx >= 380
c	2	.OR. ABS(fac(ip)-1.d0) > 1.d-5
c	2	.OR. cx >= 358 .AND. cx <= 362
c	3	.OR.radiatif
c	 deriv=cx == 100
c	 deriv=.NOT.radiatif
c	 deriv= cx >= 299 .AND. cx <= 303

	 IF(.NOT. deriv)RETURN	!tests de mise au point pour dérivées

	 ALLOCATE(bes(ne),dern(ne),yd(ne,0:1)) ; yd=y
	 IF(mn > 0.5d0)THEN
	  unpdd=1.d0-dd
	 ELSE
	  unpdd=1.d0+dd
	 ENDIF
	 IF(radiatif)THEN
	  PRINT*,'RADIATIF'
	 ELSE
	  PRINT*,'CONVECTIF'
	 ENDIF
	 PRINT*,'cx,compt,radiatif,pturb,der,n_qs_t',
	1 cx,compt,pturb,der,n_qs	
	 PRINT*,'ctel,ctem,ctep,cter,ctet,fac,Q'
	 WRITE(*,2000)ctel,ctem,ctep,cter,ctet,fac(ip),
	1 fac(ip)*(yd(4,0)*ctel+yd(1,0)*ctep+yd(5,0)*ctem+yd(3,0)*cter)
	 PRINT*,'cte21,cte22,cte23,cte24,cte25'
	 WRITE(*,2000)cte21,cte22,cte23,cte24,cte25	
	 PRINT*,'mstar,dt,cte6,alpha' ; WRITE(*,2000)mstar,dt,cte6,alpha	
	 PRINT*,'be,xcoll(cx)' ; WRITE(*,2000)be(1:ne),xcoll(cx)
	 PRINT*,'y / dy' ; WRITE(*,2000)yd(1:ne,0) ; WRITE(*,2000)yd(1:ne,1)
	 PRINT*,'mn,prn,pgn,trn,ray,ln,epsilon,ro'
	 WRITE(*,2000)mn,prn,pgn,trn,ray,ln,epsilon(1),ro
	 PRINT*,'xchim(i)' ; WRITE(*,2000)xchim
	 PRINT*,'dxchim(i)' ; WRITE(*,2000)dxchim
	 PRINT*,'despx(i)' ; WRITE(*,2000)depsx
	 PRINT*,'grad,dgradx' ; WRITE(*,2000)grad,dgradx	
	 PRINT*,'dgradpt,dgradt,dgradr,dgradl,dgradm,dgradp,dgradlpt,dgradlp'
	 WRITE(*,2000)dgradpt,dgradt,dgradr,dgradl,dgradm,dgradp,dgradlpt,dgradlp
	 PRINT*,'dgradlpp,gradad,dgradadp,dgradadt,dgradadx,t_t,p_t'
	 WRITE(*,2000)dgradlpp,gradad,dgradadp,dgradadt,dgradadx,t_t,p_t
	 PRINT*,'psist,teta,dpsitp,dpsitt,dpsitr,dpsitl,dpsitm,dpsist'
	 WRITE(*,2000)psist,teta,dpsitp,dpsitt,dpsitr,dpsitl,dpsitm,1.d0/teta
	 PRINT*,'dtetap,dtetat,dtetar,dtetal,dtetam,dp_tm,dt_tm'
	 WRITE(*,2000)dtetap,dtetat,dtetar,dtetal,dtetam,dp_tm,dt_tm
	 PRINT*,'depsp,depst,depsm,tdst,dtdsp,dtdst,dtdsl,dtdsm'
	 WRITE(*,2000)depsp*epsilon(1),depst*epsilon(1),depsm*epsilon(1),tdst,
	1 dtdsp,dtdst,dtdsl,dtdsm
	 PRINT*,'tdst,dlpp,dlppdpt,dlppdp,dy71,dy72'
	 WRITE(*,2000)tdst,dlpp,dlppdpt,dlppdp,dy71,dy72
	 PRINT*,'dy71pt,dy71t,dy71r,dy71l,dy71m,dy71p,dy71lpt,dy71lp'
	 WRITE(*,2000)dy71pt,dy71t,dy71r,dy71l,dy71m,dy71p,dy71lpt,dy71lp
	 PRINT*,'dy72pt,dy72t,dy72r,dy72l,dy72m,dy72p,dy72lpt,dy72lp'
	 WRITE(*,2000)dy72pt,dy72t,dy72r,dy72l,dy72m,dy72p,dy72lpt,dy72lp	
	 PRINT*,'u,dup,dut,dux,dum,u_t,du_tm,dtdst_m'
	 WRITE(*,2000)u,dup,dut,dux,dum,u_t,du_tm,dtdst_m
	 PRINT*,'ro,drop,drot,drox,drom,ro_t,dro_tm'
	 WRITE(*,2000)ro,drop,drot,drox,drom,ro_t,dro_tm
	 PRINT*,'gam,dgamx,dgamlpp' ; WRITE(*,2000)gam,dgamx,dgamlpp	
	 PRINT*,'dgampt,dgamt,dgamr,dgaml,dgamm,dgamp,dgamlpt,dgamlp'
	 WRITE(*,2000)dgampt,dgamt,dgamr,dgaml,dgamm,dgamp,dgamlpt,dgamlp
	 PRINT*,'v,dvdr,dvdm,dyv,dyvdr,dyvdm,w'
	 WRITE(*,2000)v,dvdr,dvdm,dyv,dyvdr,dyvdm,w
	 		
	 psist0=psist ; teta0=teta ; grad0=grad ; eps0=epsilon(1)
	 ro0=ro ; u0=u ; gam0=gam ; x00=xchim(1) ; dy710=dy71 ; dy720=dy72
	 p_t0=p_t ; t_t0=t_t ; ro_t0=ro_t ; u_t0=u_t ; tds0=tdst
	
c	 DO jd=0,1	!pour tester les dérivées/ y(i,0) et y(i,1)
	 DO jd=0,0	!pour tester les dérivées/ y(i,0)	 	
	  DO id=1,ne
	   stor0=yd(id,jd) ; stor=stor0*unpdd
	   IF(ABS(stor) < dd)stor=sign(dd,stor)
	   dstor=stor-stor0 ; yd(id,jd)=stor ; ay3=ABS(yd(3,0))
	   ray=SQRT(ay3)	!rayon
	   ay4=ABS(yd(4,0))
	   l13=SQRT(ay4)	!l**1/3
	   ln=l13**3		!l/lsol
	   ay5=ABS(yd(5,0))
	   m13=SQRT(ay5)	!m**1/3
	   mn=m13**3		!masse/mtot
	   prn=EXP(yd(1,0))	!Ptot cgs
	   IF(pturb)THEN	!avec pression turbulente 7 inconnues
	    pgn=EXP(yd(Ipg,0))	!Pgaz cgs
	   ELSE			!sans pression turbulente 6 inconnues
	    pgn=prn
	   ENDIF
	   trn=EXP(yd(2,0))		!temperature K
c	   PAUSE'avant bsp1dn'
	 
	   CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1  knotc,.TRUE.,MIN(ay5,mc(n_ch)),l,xchim,dxchim)
	   IF(no_croiss)PRINT*,'Pb. at 6 in static_m'	
	   xchim=ABS(xchim)
	 
	   IF(pturb .AND. der)THEN !der=t on tient compte de dln Pgaz/dln Ptot
	    dlpp=yd(Ipg,1)/yd(1,1)		!dlpp=dln Pgaz/dln Ptot
	   ELSE
	    dlpp=1.d0	!der=f on ne tient pas compte de dln Pgaz/dln Ptot
	   ENDIF

	   CALL thermo(prn,pgn,trn,mn,ln,ray,dlpp,xchim,dxchim,
	1  ro,drop,drot,drox,u,dup,dut,dux,
	2  grad,dgradpt,dgradp,dgradt,dgradx,dgradm,dgradl,dgradr,dgradlpp,
	3  gam,dgampt,dgamp,dgamt,dgamx,dgamm,dgaml,dgamr,dgamlpp,		
	4  epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
	5  delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	6  gradad,dgradadp,dgradadt,dgradadx,
	7  hp,dhppt,dhpp,dhpt,dhpx,dhpr,dhpm,grad_mu,
	8  gradrad,alfa,beta,gamma1,radiatif)
	
	   SELECT CASE(Krot)
	   CASE(0,1,2)
	    w=wrot
	   CASE(3,4,5)
	    CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,
	1   ay5,l,frot,dfrot)
	    IF(no_croiss)PRINT*,'Pb. at 7 in static_m'
	    w=frot(1)
	   END SELECT

	   v=SQRT(ay5/ay3)/prn*w**2	!rotation
	   
	   teta=fac(ip)*((ay5/ay3)**2/prn*(cte21+cte22*grad)+cte25*v
	1  +cte23/ro*m13/ray+cte24*m13/l13*epsilon(1)+ctem)
		
	   psist=yd(6,0)/teta  !psi/teta

	   IF(dt /= 0.d0)THEN
c	    PRINT*,dt
	    IF(compt == 0)THEN		!interpolation du TdS au temps t
	     CALL bsp1dn(1,tds_t,x_tds_t,xt_tds_t,n_tds_t,m_tds,
	1    knot_tds_t,.TRUE.,MIN(ay5,x_tds_t(n_tds_t)),l,ftds,dftds)
	     IF(no_croiss)PRINT*,'Pb. at 8 in static_m'
	     tdst=ftds(1) ; dtdst_m=dftds(1)
c	     PRINT*,'tdst,mn,l',tdst,dtdst_m,mn,l

	    ELSE				!TdS au temps t+dt
	     CALL bsp1dn(1,old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm,.TRUE.,
	1    MIN(ay5,x_ptm(n_ptm)),l,fm,dfm)
	     IF(no_croiss)PRINT*,'Pb. at 9 in static_m'	
	     mk=fm(1)	 	  
	     CALL inter('m23',bp_t,q_t,qt_t,n_qs_t,knot_t,mk,fqs,dfqs,
	1    r2_t,m23_t)
	     IF(pturb)THEN		!avec pression turbulente 7 inconnues
	      p_t=EXP(fqs(Ipg))
	     ELSE		!sans pression turbulente 6 inconnues
	      p_t=EXP(fqs(1))	 
	     ENDIF
	     t_t=EXP(fqs(2))
	    
	     IF(kipp)THEN	!approximation de Kippenhahan
	      dpn=pgn-p_t ; dtn=trn-t_t ; tdst=(cp*dtn-delta/ro*dpn)/dt
	     ELSE		!TdS=dU+PdV
	      CALL bsp1dn(nchim,chim_t,mc_t,mct_t,n_ch_t,m_ch,
	1     knotc_t,.TRUE.,MIN(mk,mc_t(n_ch_t)),l,xchim_t,dxchim_t)
	      IF(no_croiss)PRINT*,'Pb. at 10 in static_m'
	      CALL chim_gram(xchim_t,dxchim_t)
	      CALL etat(p_t,t_t,xchim_t,.FALSE.,
	1     ro_t,drop_t,drot_t,drox_t,u_t,dup_t,dut_t,dux_t,
	2     delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3     gradad,dgradadp,dgradadt,dgradadx,	
	4     alfa,beta,gamma1)
	      tdst=(u-u_t-pgn/ro**2*(ro-ro_t))/dt   !+Tds
	     ENDIF		!kipp
	    ENDIF		!compt=0
	   ENDIF		!dt/=0

	   dy1=cte11*(ay5/ay3)**2/prn ; dyv=cte15*v ; dy0=dy1+dyv
	   bes(1)=yd(1,1)-dy0*psist
	   bes(2)=yd(2,1)-dy0*psist*grad  !dln t/dq=dln p/dq grad psi/teta
	   dy3=cte13/ro*m13/ray
	   bes(3)=yd(3,1)-dy3*psist		!dl/dq=epsilon m13 psi/teta
	   bes(4)=yd(4,1)-cte14*epsilon(1)*m13/l13*psist
	   IF(dt /= 0.d0)THEN
	    dy5=cte14t*m13/l13
	    tdst=dy5*tdst
	    bes(4)=bes(4)+tdst*psist		!epsilon-tdS
	   ENDIF				!dt
	   bes(5)=yd(5,1)-psist
	   bes(6)=yd(6,1)
	   IF(pturb)THEN
            IF(radiatif)THEN
             bes(Ipg)=yd(1,1)-yd(Ipg,1)
            ELSE
             gam1=gam/(gam+1.d0) ; dy71=cte6*gam1*delta*prn*gradad*dlpp
             dy72=cte6*gam1*delta*prn*grad ; bes(Ipg)=dy71-dy72-pgn+prn
	    ENDIF		!radiatif
	   ENDIF		!Pturb
 
	   dern=(bes-be)/dstor ; PRINT*
	   SELECT CASE(jd)
	   CASE(0)	 
	    PRINT*,'dérivée / ',variable(id)
	   CASE(1)
	    PRINT*,'dérivée / d',variable(id)
	   END SELECT	 
	   IF(radiatif)THEN
	    PRINT*,'RADIATIF'
	   ELSE
	    PRINT*,'CONVECTIF'
	   ENDIF
c	   WRITE(*,2000)dgradp,dgradp*pgn/grad,pgn,prn
c	   WRITE(*,2000)dgradpt,dgradpt*prn/grad,be(2),bes(2)	   
c	   WRITE(*,2000)teta,teta0,be(2),bes(2)
	   WRITE(*,2000)(dern(j),ae(j,id,jd),j=1,4)
	   WRITE(*,2000)(dern(j),ae(j,id,jd),j=5,ne)	   	  
	   yd(id,jd)=stor0	   	  
	   PRINT*,'grad,teta,eps,psist,ro,ro_t,u,u_t'	 
	   WRITE(*,2000)(grad-grad0)/dstor/grad0,(teta-teta0)/dstor,
	1  (epsilon(1)-eps0)/dstor,(psist-psist0)/dstor,
	2  (ro-ro0)/dstor,(ro_t-ro_t0)/dstor,(u-u0)/dstor,(u_t-u_t0)/dstor
	   PRINT*,'p_t,t_t,X,gam,dy71,dy72,tds'
	   WRITE(*,2000)(tdst-tds0)/dstor,
	1  (p_t-p_t0)/dstor,(t_t-t_t0)/dstor,(xchim(1)-x00)/dstor,
	2  (gam-gam0)/dstor,(dy71-dy710)/dstor,(dy72-dy720)/dstor,
	3  (tdst-tds0)/dstor	   
	  ENDDO	!jd
	  PRINT*,'cx=',cx,ip,fac(ip)
	  PAUSE'test dérivées'
	 ENDDO !id
	 DEALLOCATE(bes,dern,yd)
	 RETURN
	
	CASE(2)		!les limites

c conditions aux limites : résidu be(1) au point limite li
c	 li=1 : limite sur r au centre
c	 li=2 : limite sur l au centre
c	 li=3 : limite sur m au centre
c	 li=4 : limite sur Ptot au raccord
c	 li=5 : limite sur T au raccord
c	 li=6 : limite sur m au raccord
c	 li=Ipg : limite sur Pgaz au raccord avec Pturb

c	 deriv=.TRUE.	   
	 IF(deriv)THEN
	  ALLOCATE(yd(ne,0:1)) ; yd=y ; unpdd=1.d0-dd
	 ENDIF
	 SELECT CASE(li)
	 
c condition sur r au centre	 
	 CASE(1)		!au centre en r
	  be(1)=y(3,0)		!en q=1 r**2=0
	  ae(1,3,0)=1.d0	!dérivée/r**2
c	  PRINT*,'limite au centre',li
c	  WRITE(*,2000)y(1:ne,0),be(1)

c condition sur l au centre
	 CASE(2)		!au centre en l
	  be(1)=y(4,0)		!en q=1 l**2/3=0
	  ae(1,4,0)=1.d0	!dérivée/l**2/3
c	  PRINT*,'limite au centre',li
c	  WRITE(*,2000)y(1:ne,0),be(1)

c condition sur m au centre
	 CASE(3)		!au centre en m
	  be(1)=y(5,0)		!en q=1 m**2/3=0
	  ae(1,5,0)=1.d0	!dérivée/m**2/3
c	   PRINT*,'limite au centre',li
c	   WRITE(*,2000)(y(i),i=1,ne),be(1)
c	   PAUSE'au centre'

c condition sur Ptot au raccord
	 CASE(4)
	  CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1 knotc,.TRUE.,MIN(y(5,0),mc(n_ch)),l,xchim0,dxchim)
	  IF(no_croiss)PRINT*,'Pb. at 11 in static_m'
	  xchim0=ABS(xchim0)	  
	  CALL chim_gram(xchim0,dxchim)
	  ray=SQRT(y(3,0)) ; l13=SQRT(y(4,0)) ; ln=l13**3
	  CALL atm(.FALSE.,ln,ray,xchim0,ptext0,dptdl0,dptdr0,text0,dtdl0,dtdr0,
	1 mext0,dmdl0,dmdr0,pext0,dpdl0,dpdr0,teff)
	  be(1)=y(1,0)-LOG(ptext0)		!condition sur Ptot au raccord
	  ae(1,1,0)=1.d0				!dérivée/ln Ptot
	  ae(1,3,0)=-dptdr0/ptext0/2.d0/ray		!dérivée/r**2
	  ae(1,4,0)=-dptdl0/ptext0*3.d0/2.d0*l13	!dérivée/l**2/3
	  IF(deriv)THEN		!test dérivée
	   PRINT*,'limite pext',li ; WRITE(*,2000)y(1,0),LOG(ptext0)
c	   WRITE(*,2000)y(1:ne,0)
c	   PRINT*,'exterieur lnp,ln ptext,be(1),ray,xchim0'
c	   WRITE(*,2000)y(1,0),LOG(ptext),be(1),ray,ln,xchim0
	   ptext0=ptext ; text0=text ; mext0=mext ; pext0=pext
	   DO i=3,4
	    stor0=yd(i,0) ; stor=stor0*unpdd
	    IF(stor == 0.d0)stor=dd
	    dstor=stor-stor0 ; yd(i,0)=stor
	    ray=SQRT(yd(3,0)) ; ln=SQRT(yd(4,0))**3
	    CALL atm(.FALSE.,ln,ray,xchim0,ptext,dptdl,dptdr,text,dtdl,dtdr,
	1   mext,dmdl,dmdr,pext,dpdl,dpdr,teff)
c	    PRINT*,LOG(ptext),yd(1,0),yd(1,0)-LOG(ptext),be(1)
	    yd(i,0)=stor0
	    PRINT*,'dérivée en li=4 / ',variable(i)
	    WRITE(*,2000)(yd(1,0)-LOG(ptext)-be(1))/dstor,ae(1,i,0)
	    yd(i,0)=stor0	
	   ENDDO	!i
	   PAUSE'test deriv'	  
	  ENDIF		!deriv

c condition sur T au raccord
	 CASE(5)
	  be(1)=y(2,0)-LOG(text0)		!condition sur T au raccord
	  ae(1,2,0)=1.d0			!dérivée/ln t
	  ae(1,3,0)=-dtdr0/text0/2.d0/ray		!dérivée/r**2
	  ae(1,4,0)=-dtdl0/text0*3.d0/2.d0*l13	!dérivée/l**2/3
	  IF(deriv)THEN 	!test de dérivation
	   PRINT*,'limite Text',li ; WRITE(*,2000)y(2,0),LOG(text)
c	   WRITE(*,2000)y(1:ne,0)
c	   PRINT*,'y(2,0),LOG(text),be(1)'
c	   WRITE(*,2000)y(2,0),LOG(text),be(1)
	   DO i=3,4
	    stor0=yd(i,0) ; stor=stor0*unpdd
	    IF(stor == 0.d0)stor=dd
	    dstor=stor-stor0 ; yd(i,0)=stor
	    ray=SQRT(yd(3,0)) ; ln=SQRT(yd(4,0))**3
	    CALL atm(.FALSE.,ln,ray,xchim0,ptext,dptdl,dptdr,text,dtdl,dtdr,
	1   mext,dmdl,dmdr,pext,dpdl,dpdr,teff)
	    PRINT*,'dérivée pour li=5 / ',variable(i)
	    WRITE(*,2000)(yd(2,0)-LOG(text)-be(1))/dstor,ae(1,i,0)
	    yd(i,0)=stor0	     
	   ENDDO	!i
	   PAUSE'test deriv'	   	  
	  ENDIF		!deriv
	  
c condition sur M au raccord	  	  	 
	 CASE(6)
	  be(1)=y(5,0)-mext0**(2.d0/3.d0)	!en q=n_qs, m=mext(r,l)
	  ae(1,3,0)=-dmdr0/3.d0/ray/m13		!dérivée/r**2
	  ae(1,4,0)=-dmdl0*l13/m13		!dérivée/l**2/3
	  ae(1,5,0)=1.d0			!dérivée/m**2/3
c	  WRITE(*,2000)y(5,0)
	  IF(deriv)THEN 	!test de dérivation
	   PRINT*,'limite mext',li
	   WRITE(*,2000)y(5,0),mext0**(2.d0/3.d0)
	   WRITE(*,2000)y(1:ne,0)
	   PRINT*,'limite y(5),mext',li
	   PRINT*,'y(5,0),mext,be(1)'
	   WRITE(*,2000)y(5,0),mext0,be(1)
	   DO i=3,4
	    stor0=yd(i,0) ; stor=stor0*unpdd
	    IF(stor == 0.d0)stor=dd
	    dstor=stor-stor0 ; yd(i,0)=stor
	    ray=SQRT(yd(3,0)) ; ln=SQRT(yd(4,0))**3
	    CALL atm(.FALSE.,ln,ray,xchim0,ptext,dptdl,dptdr,text,dtdl,dtdr,
	1   mext,dmdl,dmdr,pext,dpdl,dpdr,teff)
	    PRINT*,'dérivée pour li=6 / ',variable(i)
	    WRITE(*,2000)(yd(5,0)-mext**(2.d0/3.d0)-be(1))/dstor,ae(1,i,0)
	    yd(i,0)=stor0	    
	   ENDDO	!i
	   PAUSE'test deriv'	  
	  ENDIF		!deriv	  

c condition sur Pgaz au raccord	  	  	  
	 CASE(7)
	  be(1)=y(Ipg,0)-LOG(pext0)			!en q=n
	  ae(1,3,0)=-dpdr0/pext0/2.d0/ray		!dérivée/r**2
	  ae(1,4,0)=-dpdl0/pext0*3.d0/2.d0*l13	!dérivée/l**2/3
	  ae(1,Ipg,0)=1.d0			!dérivée/ln Pgaz
	  IF(deriv)THEN 	!test de dérivation
	   PRINT*,'limite pext, uint',li
	   WRITE(*,2000)y(Ipg,0),pext0 ; WRITE(*,2000)y(1:ne,0)
	   PRINT*,'limite y(Ipg)',li ; PRINT*,'y(Ipg,0),pext,be(1)'
	   WRITE(*,2000)y(Ipg,0),pext0,be(1)
	   DO i=3,4
	    stor0=yd(i,0) ; stor=stor0*unpdd
	    IF(stor == 0.d0)stor=dd
	    dstor=stor-stor0 ; yd(i,0)=stor
	    ray=SQRT(yd(3,0)) ; ln=SQRT(yd(4,0))**3	   
	    CALL atm(.FALSE.,ln,ray,xchim,ptext,dptdl,dptdr,text,dtdl,dtdr,
	1   mext,dmdl,dmdr,pext,dpdl,dpdr,teff)
	    PRINT*,'dérivée pour li=7 / ',variable(i)
	    WRITE(*,2000)(yd(Ipg,0)-LOG(pext)-be(1))/dstor,ae(1,i,0)
	    yd(i,0)=stor0	    
	   ENDDO	!i
	   PAUSE'test deriv'	   	  
	  ENDIF		!deriv

	 CASE DEFAULT
	  WRITE(*,"('static_m erreur sur li =',i3,', ne=',i3)")li,ne
	  PRINT*,'ARRET'
	  STOP		
	 END SELECT
	 IF(deriv)DEALLOCATE(yd)	 
c	 PRINT*,'li=',li ; PAUSE'fin des limites'
	
	CASE DEFAULT
	 PRINT*,' static_m erreur sur fait = 1, ou 2, fait=',fait
	 PRINT*,'ARRET' ; STOP	
	END SELECT		!fait

	RETURN

	END SUBROUTINE static_m
