
c********************************************************************

      SUBROUTINE static_r(fait,cx,li,y,be,ae,compt,dt,reprend,ip)

c     routine private du module mod_static

c     Calcul des quantités caractéristiques du problème quasi-statique
c     en eulérien
c     routine appellée par resout

c     spline collocation avec fonction de répartition

c     avec pression turbulente 7 inconnues 
c     sans pression turbulente 6 inconnues, Ptot=Pgaz

c     pour avoir un algorithme unique pour la diffusion, la composition
c     chimique est toujours tabulée en fonction de mu=(m/Msol)^2/3
c     que ce soit en lagrangien ou en eulérien

c     l'énergie graviphique, TdS/dt=tds est tabulée en fonction
c     de m23=m^2/3 en lagrangien et de m23=m en eulérien,

c     Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c     CESAM2k

c    08 10 96 : introduction de la rotation
c    08 01 96 : rectification de ctem*fac
c    05 05 97 : 2 possibilités de calcul du TdS
c    26 06 97 : remplacement du moment angulaire par la vitesse angulaire
c    18 07 97 : correction d'énergie excédentaire de rotation
c    30 07 00 : introduction F95

c entrées
c     fait=1
c     -calcule les résidus be(i) des équations au point xcoll i=1,ne
c     fait=2
c     -calcule le "résidu" be(1) de la condition au point limite li

c     cx : indice du point de collocation
c     li : numéro de la limite
c     y : variables au point de collocation xcoll(cx)
c     compt: compteur du nb. iter. Newton Raphson
c     dt : pas temporel
c     ip : indice de la couche pour détermination du facteur de répartition

c sorties
c     be : résidus
c     ae : termes du jacobien pour le point de collocation
c     reprend=.TRUE. : variation relative de U ou Ro trop forte

c     dans cette routine on est toujours en eulérien
c     en_masse = .FALSE. variables eulériennes m23=m, r2=r

c----------------------------------------------------------------

      USE mod_donnees, ONLY : alpha, cpturb, ctel, ctep, ctem,
     1 cter, ctet, dtmin, d_grav, g, Ipg, Krot, kipp, langue, lsol,
     2 msol, m_ch, m_ptm, m_rot, m_tds, nchim, ne, nrot, ord_rot, pi, pturb,
     3 ro_test, rsol, secon6, t_inf
      USE mod_atm, ONLY : atm
      USE mod_etat, ONLY : etat
      USE mod_kind
      USE mod_numerique, ONLY : bsp1dn, no_croiss
      USE mod_variables, ONLY : bp, bp_t, chim, chim_gram, chim_t, inter,
     1 knotc, knotc_t, knotr, knot_tds_t, knot_ptm, knot_t, mc, mct, mct_t,
     2 mc_t, mrot, mrott, mstar, m23_t, n_ch, n_ch_t, n_ptm, n_qs_t,
     3 n_rot, n_tds_t, old_ptm, qt_t, q_t, rota, r2_t, sortie, tds_t, wrot,
     4 xt_ptm, xt_tds_t, x_ptm, x_tds_t

      IMPLICIT NONE
      
      REAL (kind=dp), INTENT(in), DIMENSION(:,0:) :: y    
      REAL (kind=dp), INTENT(in) :: dt        
      INTEGER, INTENT(in) :: cx, compt, fait, li, ip
      
      REAL (kind=dp), INTENT(out), DIMENSION(:,:,0:) :: ae    
      REAL (kind=dp), INTENT(out), DIMENSION(:) :: be 
      LOGICAL, INTENT(out) :: reprend
      REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: yd
      REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: bes, dern
      REAL (kind=dp), DIMENSION(nchim) :: depsx, dxchim,
     1 dxchim_t, xchim, xchim_t
      REAL (kind=dp), DIMENSION(ne) :: dfqs, fqs
      REAL (kind=dp), DIMENSION(nrot) :: dfrot, frot  
      REAL (kind=dp), DIMENSION(5) :: epsilon
      REAL (kind=dp), DIMENSION(1) :: dfm, dftds, fm, ftds
             
      REAL (kind=dp), PARAMETER :: dd=1.d-5   
      REAL (kind=dp), SAVE :: cte2, cte6, cte11, cte13, cte14, cte14t,
     1 cte15, cte21, cte22, cte23, cte24, cte25, dmdl, dmdr, dpdl,
     2 dpdr, dptdl, dptdr, dtdl, dtdr, dt_tds, mext, pext, unpdd
      
      REAL (kind=dp) :: acc_c, alfa, am13, am23, ay3, ay5, beta, cp, dcpm,
     1 dcpp, dcpt, dcpx, delta, deltam, deltap, deltat, deltax, depsm,
     2 depsp,
     3 depst, dgamlpt, dgamp, dgampt, dgamt, dgamx, dgamm, dgaml, dgamlp,
     4 dgamlpp, dgamr, dgam1, gradad, dgradadm, dgradadp, dgradadt,
     5 dgradadx, dgradlp, dgradlpt, dgradp, dgradpt, dgradt, dgradx,
     6 dgradm, dgradl, dgradr, dgradlpp, dhpm, dhpp, dhppt, dhpr,
     7 dhpt, dhpx, dkapp, dkapt,
     8 dkapx, dlpp, dlppdp, dlppdpt, dpn, dpsitl, dpsitm, dpsitp, dpsitpt,
     9 dpsitr, dpsitt, drom, drop, drop_t, drot, drot_t, drox, drox_t,
     1 dro_tm=0.d0, dstor, dtdsm=0.d0, dtdsp=0.d0, dtdsr, dtdst=0.d0,
     2 dtdst_m=0.d0, dtetal, dtetam, dtetap, dtetapt, dtetar, dtetat, dtn,
     3 dp_tm=0.d0, dt_tm=0.d0, dum, dup, dup_t, dut, dut_t, duv, dux,
     4 dux_t, du_tm=0.d0, dvdm, dvdp, dvdr, dvdt, dwdm,
     5 dyv, dyvdm, dyvdp, dyvdr, dyvdt, dy0, dy1, dy2, dy3, dy4,
     6 dy5, dy71=0.d0, dy710, dy71l=0.d0, dy71lp=0.d0, dy71lpt=0.d0,
     7 dy71m=0.d0, dy71p=0.d0, dy71pt=0.d0, dy71r=0.d0, dy71t=0.d0,
     8 dy72=0.d0, dy720, dy72l=0.d0, dy72lp=0.d0, dy72lpt=0.d0,
     9 dy72p=0.d0, dy72pt=0.d0, dy72m=0.d0, dy72r=0.d0, dy72t=0.d0,
     1 eps0, gam, gam0, gam1, gamma1, grad, grad0 , gradrad, grad_mu, gravs2,
     2 hp, kap, mext0, mk, pext0, pgn, prn, psist, psist0,
     3 ptext, ptext0, p_t=0.d0, p_t0, ro, ro0, ro_t=0.d0, ro_t0,
     4 stor, stor0, tdst=0.d0, tds0, teff, teta, teta0, text, text0,
     5 trn, t_t=0.d0, t_t0, u, u_t=0.d0, u_t0, u0, x00, v, w

      INTEGER :: i, id, j, jd, l=1
      
      LOGICAL, SAVE :: init=.TRUE., der
      LOGICAL :: radiatif, deriv=.FALSE.
      
      CHARACTER (len=7), SAVE, DIMENSION(7) :: variable=(/ 'ln Ptot',
     1 'ln T   ','R      ','L      ','M      ','Psi    ', 'ln Pgaz'/)
      
c------------------------------------------------------------------------

2000	FORMAT(10es10.3)

c	PRINT* ; PRINT*,'entrée static_r',fait,compt,cx,ip
c	WRITE(*,2000)dt,xcoll(cx) ; WRITE(*,2000)y(1:ne,0),fac(ip)
c	WRITE(*,2000)y(1:ne,1)  
      
	IF(init)THEN
	 init=.FALSE.
     
c le formalisme en eulérien n'est pas disponible										
	IF(.TRUE.)THEN
c	IF(.FALSE.)THEN
	 	 
	 SELECT CASE(langue)
	 CASE('english')
          WRITE(2,1021) ; WRITE(*,1021) ; CALL sortie  
1021      FORMAT(/,'STOP', /,
	1 'The formalism with eulerian variables is not available, sorry')	
	 CASE DEFAULT	
          WRITE(2,21) ; WRITE(*,21) ; CALL sortie   
21         FORMAT(/,'ARRET',/,
	1 'Le formalisme avec les variables eulériennes n''est pas disponible')
	 END SELECT
	 
	ENDIF      
     
	 cte2=2.d0/3.d0*rsol
	 cte11=-g*msol/rsol		!pour d ln P /d q
	 cte21=cte11*ctep ; cte22=cte11*ctet    
	 cte13=4.d0*pi*rsol**3/msol	!pour d mu /d q
	 cte23=cte13*ctem     
	 cte14=4.d0*pi*rsol**3/lsol     !pour d lambda/d q
	 cte24=cte14*ctel
	 cte14t=cte14/secon6            ! "    " dt .ne. 0   
	 cte15=2.d0/3.d0*rsol           !pour  la rotation
	 cte25=cte15*ctep    
	 cte6=ABS(cpturb)*alpha**2/8.d0     !pour la Pturb
	 der=cpturb < 0.d0  !der=.TRUE. on tient compte de dlnPgaz/dlnPtot
	 dt_tds=3.d0*dtmin	!dt < dt_tds on utilise la tabulation TdS

	 IF(kipp)THEN
	  SELECT CASE(langue)
	   CASE('english')
	   WRITE(2,1010) ; WRITE(*,1010)   
1010	   FORMAT(/,'Kippenhahn''s approximation for TdS',/)	
	  CASE DEFAULT	
	   WRITE(2,10) ; WRITE(*,10)   
10	   FORMAT(/,'Approximation de Kippenhahn pour le calcul du TdS',/)
	  END SELECT
	 ELSE
	  SELECT CASE(langue)
	   CASE('english')
	   WRITE(*,1011) ; WRITE(2,1011)       
1011	   FORMAT(/,'Full calculation of TdS',/)       
	  CASE DEFAULT       
	   WRITE(*,11) ; WRITE(2,11)       
11	   FORMAT(/,'Calcul exact du TdS',/)
	  END SELECT
	 ENDIF
	ENDIF       !init
      
c	PRINT*,'fait,cx',fait,cx

c------------------------------------------------------------------

c ensemble des variables
c         y(1,0)=ln Ptot  		y(1,1)=(ln Ptot)'
c         y(2,0)=ln T    	 	y(2,1)=(ln T)'  
c         y(3,0)=r        		y(3,1)=r'
c         y(4,0)=l        		y(4,1)=l'
c         y(5,0)=m        		y(5,1)=m'
c         y(6,0)=psi=dQ/dq(=cte)	y(6,1)=psi'(=0)
c avec Pturb  y(Ipg,0)=ln Pgaz      	y(Ipg,1)=(ln Pgaz)'

c xcoll(cx) valeur de q au point de collocation q=1, 2, 3, .. , n

c ae(eq,var,der)=ae(ne,ne,0:1)=dérivée de la eq-ieme équation par
c rapport a la der-ieme dérivée de la var-ieme variable
	ae=0.d0 ; be=0.d0 ; reprend=.FALSE.

	SELECT CASE(fait)   
	CASE(1)     !le point courant      
	 prn=EXP(y(1,0))    !pression Ptot cgs
	 IF(pturb)THEN      !avec pression turbulente 7 inconnues
	  pgn=EXP(y(Ipg,0))   !pression Pgaz cgs   
	 ELSE           !sans pression turbulente 6 inconnues
	  pgn=prn
	 ENDIF
	 trn=EXP(y(2,0))        !temperature K
	 ay3=ABS(y(3,0))        !rayon
	 ay5=ABS(y(5,0))        !masse
	 am23=ay5**(2.d0/3.d0)  !m**2/3 pour int. comp. chim. 
	 am13=SQRT(am23)        !m**1/3 pour int. comp. chim. 

c l'acc. centrifuge ne doit pas excéder	90% gravité
	 gravs2=-cte11*ay5/ay3*0.9d0

c vitesse angulaire
	 SELECT CASE(Krot)
	 CASE(0,1,2)
	  w=wrot
	 CASE(3,4,5)
	  CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,
	1 MAX(mrot(1),MIN(am23,mrot(n_rot))),l,frot,dfrot)	
	  IF(no_croiss)PRINT*,'Pb. en 0 dans static_r'
	  w=frot(1)
	 END SELECT
	 
c accélération centrifuge	 
	 acc_c=cte2*ay3*w**2 ; reprend=acc_c > gravs2
	 IF(reprend)THEN
	  SELECT CASE(langue)
	  CASE('english')
	   WRITE(*,1012)acc_c,gravs2,w,ay3,ay5
	   WRITE(2,1012)acc_c,gravs2,w,ay3,ay5
1012	   FORMAT('STOP, in static_r, the centrifugal acceleration =',
	1  es10.3,' > 90% gravity = ',es10.3,/,'angular velicity =',es10.3,
	2  ', R/Rsun=',es10.3,', M/Msun=',es10.3)	   
	  CASE DEFAULT
	   WRITE(*,12)acc_c,gravs2,w,ay3,ay5 ; WRITE(2,12)acc_c,gravs2,w,ay3,ay5
12	   FORMAT('ARRET, dans static_r, accélération centifuge =',
	1  es10.3,' > 90% gravité = ',es10.3,/,'vitesse angulaire =',es10.3
	2  ', R/Rsun=',es10.3,', M/Msun=',es10.3)	  
	   END SELECT
	   CALL sortie
	  ENDIF

c si, à cause d'une erreur d'arrondi, ce qui peut arriver au voisinage
c du centre r, l ou m est negatif, on tente de forcer la
c convergence en utilisant les ABS(**), les dérivées sont inexactes et
c la physique est violee... mais si ca passe les erreurs d'arrondi
c disparaissent, parfois, au cours des iterations
c pour éviter cette disposition supprimer les "c" devant reprend et RETURN

	  IF(MIN(y(3,0),y(4,0),y(5,0)) < 0.d0)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1014)y(3:5,0),NINT(xcoll(cx))
1014	    FORMAT('static_r, r=',es10.3,', l=',es10.3,', m=',es10.3,/,
	1  'shell=',i4,', one try to converge')	  
	  CASE DEFAULT	  
	   WRITE(*,14)y(3:5,0),NINT(xcoll(cx))
14	   FORMAT('static_r, r=',es10.3,', l=',es10.3,', m=',es10.3,/,
	1  'couche=',i4,', on tente de converger')	  
	  END SELECT
	 ENDIF
     
c composition chimique au temps t+dt
       CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
     1 knotc,.TRUE.,MIN(am23,mc(n_ch)),l,xchim,dxchim)
       IF(no_croiss)PRINT*,'Pb. at 1 in static_r'
       xchim=ABS(xchim) ; dxchim=dxchim*2.d0/3.d0/am13     !dérivée/ m
      
c      PRINT*,'cx,nchim/prn,pgn,trn / xchim / y  /dy',cx,nchim
c      WRITE(*,2000)prn,pgn,trn,xchim(1),dxchim(1)
c      WRITE(*,2000)xchim ; WRITE(*,2000)dxchim
c      WRITE(*,2000)y(1:ne,0) ; WRITE(*,2000)y(1:ne,1)
c      PRINT*,' '
      
       IF(pturb .AND. der)THEN    !avec pression turbulente 7 inconnues
        dlpp=y(Ipg,1)/y(1,1)    !dlpp=dln Pgaz/dln Ptot
        dlppdpt=-dlpp/y(1,1)  !dérivée dlpp /dln Ptot 
        dlppdp=dlpp/y(Ipg,1)    !dérivée dlpp /dln Pgaz
       ELSE
        dlpp=1.d0 !der=.FALSE. on ne tient pas compte de dln Pgaz/dln Ptot
        dlppdpt=0.d0 ; dlppdp=0.d0
       ENDIF
       CALL thermo(prn,pgn,trn,ay5,ABS(y(4,0)),ay3,dlpp,xchim,dxchim,
     1 ro,drop,drot,drox,u,dup,dut,dux,
     2 grad,dgradpt,dgradp,dgradt,dgradx,dgradm,dgradl,dgradr,dgradlpp,
     3 gam,dgampt,dgamp,dgamt,dgamx,dgamm,dgaml,dgamr,dgamlpp,     
     4 epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
     5 delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
     6 gradad,dgradadp,dgradadt,dgradadx,hp,dhppt,dhpp,dhpt,dhpx,
     7 dhpr,dhpm,grad_mu,gradrad,alfa,beta,gamma1,radiatif)
       
c      WRITE(*,2000)y(5,0),prn,pgn,trn,ro,xchim,epsilon(1),depsp,depst,
c	1 depsx
c      PRINT*,'couche',cx
c      WRITE(*,2000)prn,pgn,trn,y(5,0),y(4,0),y(3,0),ro,grad
c      WRITE(*,2000)gam,epsilon(1),delta,gradad,hp,gradrad     
c      WRITE(*,2000)xchim ; WRITE(*,2000)epsilon ; PAUSE'thermo'

c dérivées par rapport a : ln Ptot, ln T, r, l, m, ln Pgaz
c sans Pturb les dérivées / Ptot et / Pgaz sont identiques         
       drop=drop*pgn             !der ro /ln Pgaz
       drot=drot*trn ; drom=drox*dxchim(1)

       dcpp=dcpp*pgn ; dcpt=dcpt*trn ; dcpm=dcpx*dxchim(1)
       
       deltap=deltap*pgn ; deltat=deltat*trn ; deltam=deltax*dxchim(1)
       
       dup=dup*pgn ; dut=dut*trn ; dum=dux*dxchim(1)
       
       dgradadp=dgradadp*pgn/gradad      !der gradad/ln Pgaz gradad près
       dgradadt=dgradadt*trn/gradad      !der gradad/ln T a gradad près
       dgradadm=dgradadx*dxchim(1)/gradad    !der gradad/m**2/3 a gradad près

       dgradpt=dgradpt*prn/grad      !der grad /ln Ptot a grad près
       dgradp=dgradp*pgn/grad        !der grad /ln Pgaz   
       dgradt=dgradt*trn/grad        !der grad /ln T
       dgradr=dgradr/grad            !der grad /r
       dgradl=dgradl/grad            !der grad /l
       dgradm=(dgradm+dgradx*dxchim(1))/grad     !der grad/m
       dgradlpt=dgradlpp*dlppdpt/grad        !der grad /dln Ptot  
       dgradlp =dgradlpp*dlppdp/grad         !der grad /dln Pgaz
               
       dgampt=dgampt*prn             !der gam /ln Ptot
       dgamp= dgamp*pgn              !der gam /ln Pgaz 
       dgamt= dgamt*trn              !der gam /ln T
       dgamm=dgamm+dgamx*dxchim(1)           !der gam /m   
       dgamlpt=dgamlpp*dlppdpt           !der gam /dln Ptot
       dgamlp =dgamlpp*dlppdp            !der gam /dln Pgaz
       
       dhppt=dhppt*prn/hp    !der hp /ln Ptot a hp près
       dhpp=dhpp*pgn/hp  !der hp /ln Pgaz a hp près
       dhpt=dhpt*trn/hp  !der hp /ln T a hp près
       dhpr=dhpr/hp      !der hp /r a hp près    
       dhpm=(dhpm+dhpx*dxchim(1))/hp !der hp/m a hp près
       
       IF(epsilon(1) > 1.d-25)THEN
        depsp=depsp*pgn/epsilon(1)   !dérivée /ln Pgaz, a 1/epsilon près
        depst=depst*trn/epsilon(1)
        depsm=SUM(depsx*dxchim)/epsilon(1)
       ELSE
        depsp=0.d0 ; depst=0.d0 ; depsm=0.d0
       ENDIF
   
c on ne tient pas compte de la différence Pgaz, Ptot pour
c l'énergie graviphique TdS=tds et old_ptm sont tabulés en fonction
c de m^2/3 en lagrangien et de m en eulérien,
c pour avoir un algorithme unique pour la diffusion, la composition
c chimique est toujours tabulée en fonction de m^2/3
       IF(dt /= 0.d0)THEN
        IF(compt == 0 .OR. dt <= dt_tds)THEN	!interp. du TdS au temps t
         CALL bsp1dn(1,tds_t,x_tds_t,xt_tds_t,n_tds_t,m_tds,
     1   knot_tds_t,.TRUE.,MIN(ay5,x_tds_t(n_tds_t)),l,ftds,dftds)
         IF(no_croiss)PRINT*,'Pb. at 2 in static_r'
         tdst=ftds(1) ; dtdst_m=dftds(1)
c        WRITE(*,2000)tdst,m_tds,y(5,0),x_tdst_t(n_tds_t)
c        PAUSE'TdS interpole'

        ELSE              !TdS au temps t+dt
         CALL bsp1dn(1,old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm,.TRUE.,
     1   MIN(ay5,x_ptm(n_ptm)),l,fm,dfm)
         IF(no_croiss)PRINT*,'Pb. at 3 in static_r'
         mk=fm(1)       
         CALL inter('m23',bp_t,q_t,qt_t,n_qs_t,knot_t,mk,fqs,dfqs,r2_t,m23_t)
         IF(pturb)THEN        !avec pression turbulente 7 inconnues
          p_t=EXP(fqs(Ipg)) ; dp_tm=dfqs(Ipg)*p_t   !d p_t / d mu au temps t
         ELSE         !sans pression turbulente 6 inconnues
          p_t=EXP(fqs(1)) ; dp_tm=dfqs(1)*p_t   !d p_t / d mu au temps t
         ENDIF
         t_t=EXP(fqs(2)) ; dt_tm=dfqs(2)*t_t    !d t_t / d mu au temps t

         IF(kipp)THEN !approximation de Kippenhahn
          dpn=pgn-p_t ; dtn=trn-t_t
          IF(trn > t_inf)THEN !controle de la variation
	   duv=1.d0-ABS((cp*t_t-delta/ro*p_t)/(cp*trn-delta/ro*prn))	  
c          duv=MAX(ABS(dpn/pgn),ABS(dtn/trn))
          ELSE
           duv=0.d0
          ENDIF     
          reprend=(duv > d_grav) .AND. (compt > 3)
          IF(reprend)THEN       !TdS varie trop
	   WRITE(*,7)ay5/mstar
7	   FORMAT('Approx. Kipp. trop de var. de TdS en m/Mstar=',es10.3)
	   WRITE(*,13)duv,d_grav
13	   FORMAT('duv=',es10.3,' > d_grav=',es10.3)
	   WRITE(*,8)p_t,pgn,dpn/pgn
8	   FORMAT('P(t)=',es10.3,', P(t+dt)=',es10.3,', dP/P=',es10.3)   
	   WRITE(*,9)t_t,trn,dtn/trn
9	   FORMAT('T(t)=',es10.3,', T(t+dt)=',es10.3,', dT/T=',es10.3)  
	   RETURN
          ELSE
           duv=0.d0
          ENDIF   !reprend        
          tdst=(cp*dtn-delta/ro*dpn)/dt
         
         ELSE     !TdS=dU+PdV
          mk=mk**(2.d0/3.d0)    
          CALL bsp1dn(nchim,chim_t,mc_t,mct_t,n_ch_t,m_ch,
     1    knotc_t,.TRUE.,MIN(mk,mc_t(n_ch_t)),l,xchim_t,dxchim_t)
          IF(no_croiss)PRINT*,'Pb. at 4 in static_r'
          xchim_t=ABS(xchim_t) ; dxchim_t=dxchim_t*2.d0/3.d0/am13
          CALL chim_gram(xchim_t,dxchim_t)      
c         WRITE(*,2000)p_t,t_t,xchim_t(1),xcoll(cx)
          CALL etat(p_t,t_t,xchim_t,.TRUE.,    !mettre .TRUE. pour d/dX
     1    ro_t,drop_t,drot_t,drox_t,u_t,dup_t,dut_t,dux_t,
     2    delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
     3    gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
          IF(ro > ro_test .AND. trn > t_inf)THEN
	   duv=1.d0-ABS((u_t-prn/ro**2*ro_t)/(u-prn/ro))
c          duv=MAX(ABS(u-u_t)/u,ABS(ro-ro_t)/ro)
          ELSE
           duv=0.d0
          ENDIF
          reprend=(duv > d_grav) .AND. (compt >= 3)
          IF(reprend)THEN       !TdS varie trop	  
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(*,1006)ay5/mstar
1006	    FORMAT('Too large change of dU+PdV at m/Mstar=',es10.3)
	   CASE DEFAULT
	    WRITE(*,6)ay5/mstar
6	    FORMAT('Trop de variation de dU+PdV en m/Mstar=',es10.3)
	   END SELECT
	   WRITE(*,13)duv,d_grav
	   WRITE(*,1)u_t,u,ABS(u-u_t)/u
1	   FORMAT('U(t)=',es10.3,', U(t+dt)=',es10.3,', dU/U=',es10.3)
 	   WRITE(*,2)ro_t,ro,ABS(ro-ro_t)/ro
2	   FORMAT('ro(t)=',es10.3,', ro(t+dt)=',es10.3,', dro/ro=',es10.3)
 	   WRITE(*,3)xchim_t(1),xchim(1)
3	   FORMAT('X(t)=',es10.3,', X(t+dt)=',es10.3)
	   WRITE(*,4)p_t,pgn
4	   FORMAT('P(t)=',es10.3,', P(t+dt)=',es10.3)	     
	   WRITE(*,5)t_t,trn
5	   FORMAT('T(t)=',es10.3,', T(t+dt)=',es10.3)	          
           RETURN
          ENDIF   !reprend
      
          du_tm =dup_t *dp_tm+dut_t *dt_tm+dux_t *dxchim_t(1)
          dro_tm=drop_t*dp_tm+drot_t*dt_tm+drox_t*dxchim_t(1)
          tdst=(u-u_t-prn/ro**2*(ro-ro_t))/dt
         
c         PRINT*,'cx/Ptot,Pgaz,T,p_t,t_t,u,u_t,ro,ro_t,tdst',cx
c         WRITE(*,2000)prn,trn,p_t,t_t,u,u_t,ro,ro_t,tdst

         ENDIF        !kipp
        ENDIF     !compt=0
       ENDIF      !dt /= 0

c vitesse angulaire
       SELECT CASE(Krot)
       CASE(0,1,2)
        w=wrot ; dwdm=0.d0
       CASE(3,4,5)
        CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,knotr,.TRUE.,
     1  MAX(mrot(1),MIN(am23,mrot(n_rot))),l,frot,dfrot)
        IF(no_croiss)PRINT*,'Pb. at 5 in static_r'
        w=frot(1) ; dwdm=dfrot(1)*2.d0/3.d0/am13
       END SELECT

c v correspond à l'accélération centrifuge omega^2 / 6 pi R
c      WRITE(*,2000)y(5,0),w,y(3,0)
       IF(w == 0.d0)THEN
        v=0.d0 ; dvdr=0.d0 ; dvdm=0.d0 ; dyv=0.d0 ; dyvdr=0.d0 ; dyvdm=0.d0
        dvdp=0.d0 ; dyvdt=0.d0 ; dyvdp=0.d0
        dvdt=0.d0
       ELSE
c       PAUSE'vérifier les algorithmes + dérivées si V /= 0'
        v=ro*ay3/prn*w**2
        dvdp=v*drop/ro    !dérivée/ln Pgaz
        dvdt=v*drot/ro    !dérivée/ln T
        dvdr=v/ay3        !dérivée/r
        dvdm=v*(drom/ro+2.d0/w*dwdm)
        dyv=cte25*fac(ip)*v   !dérivée/lnP comme celle de dy1: dyv/dlnp=-dyv
        dyvdp=cte25*fac(ip)*dvdp
        dyvdt=cte25*fac(ip)*dvdt ; dyvdr=cte25*fac(ip)*dvdr
        dyvdm=cte25*fac(ip)*dvdm
       ENDIF
      
c la fonction de répartition
c en eulérien, il y a un dtetapt pour Ptot et un dtetap pour Pgaz
c      car ro dépend de Pgaz = Pgaz + Prad

c teta = ctep dksi / dr
       dy1=fac(ip)*ay5/ay3**2/prn*ro
       teta=cte21*dy1
       IF(pturb)THEN
        dtetapt=-teta			!dérivée/ln Ptot
       ELSE
        dtetapt=teta*(drop/ro-1.d0)	!dérivée/ln Ptot = ln Pgaz
       ENDIF
       dtetat=teta*drot/ro		!dérivée/ln T
       dtetar=-teta*2.d0/ay3		!dérivée/r
       dtetam=teta*(drom/ro+1.d0/ay5)	!dérivée/m      
       dtetap=teta*drop/ro		!dérivée/ln Pgaz
       
c teta = teta + rotation       
       teta=teta+dyv
       dtetap=dtetap-dyvdp
       dtetat=dtetat+dyvdt
       dtetar=dtetar+dyvdr
       dtetam=dtetam+dyvdm       
       
c teta = teta + ctet dksi / dr grad      
       dy2=cte22*dy1*grad ; teta=teta+dy2
       IF(pturb)THEN
        dtetapt=dtetapt-dy2			!dérivée/ln Ptot
       ELSE
        dtetapt=dtetapt+dy2*(dgradpt+drop/ro-1.d0)!dérivée/ln Ptot=ln Pgaz
       ENDIF
       dtetat=dtetat+dy2*(dgradt+drot/ro)
       dtetar=dtetar+dy2*(-2.d0/ay3+dgradr)
       dtetal=dy2*dgradl      
       dtetam=dtetam+dy2*(1.d0/ay5+dgradm)
       dtetap=dtetap+dy2*dgradp
    
c teta = teta + ctem d m/ d r
       dy3=fac(ip)*cte23*ro*ay3**2 ; teta=teta+dy3
       IF(.NOT.pturb)dtetapt=dtetapt+dy3*drop/ro !dérivée/ln Ptot = lnPgaz
       dtetap=dtetap+dy3*drop/ro	!dérivée/ln Pgaz  
       dtetat=dtetat+dy3*drot/ro	!dérivée/ln T
       dtetar=dtetar+dy3*2.d0/ay3	!dérivée/r
       dtetam=dtetam+dy3*drom/ro	!dérivée/m

c teta = teta + ctel d l / d r
       dy4=cte24*fac(ip)*ro*ay3**2*epsilon(1) ; teta=teta+dy4
       IF(.NOT.pturb)dtetapt=dtetapt+dy4*(drop/ro+depsp)!dérivée/ln P      
       dtetap=dtetap+dy4*(drop/ro+depsp)	!dérivée/ln Pgaz 
       dtetat=dtetat+dy4*(drot/ro+depst)	!dérivée/ln T
       dtetar=dtetar+dy4*2.d0/ay3		!dérivée/r
       dtetam=dtetam+dy4*(drom/ro+depsm)	!dérivée/m

c teta = teta + cter d zeta/ d zeta       
       teta=teta+fac(ip)*cter          
            
c psi/teta
       psist=y(6,0)/teta		!psi/teta
       dpsitpt=-psist/teta*dtetapt   	!d (psi/teta) / d ln Ptot   
       dpsitt=-psist/teta*dtetat 	!d (psi/teta) / d ln T
       dpsitr=-psist/teta*dtetar 	!d (psi/teta) / d r
       dpsitl=-psist/teta*dtetal	! (psi/teta) / d l      
       dpsitm=-psist/teta*dtetam 	!d (psi/teta) / d m
       dpsitp=-psist/teta*dtetap 	!d (psi/teta) / d ln Pgas    

c équations
       dyv=cte15*v        !dérivée/ln P comme celle de dy1: dyv/dlnp=-dyv
       dyvdp=cte15*dvdp ; dyvdt=cte15*dvdt 
       dyvdr=cte15*dvdr ; dyvdm=cte15*dvdm

       dy1=cte11*ay5/ay3**2/prn*ro ; dy0=dy1+dyv
       be(1)=y(1,1)-dy0*psist 				!dln Ptot/dq
       IF(pturb)THEN
        ae(1,1,0)=dy0*(psist-dpsitpt)      		!der./ln Ptot
        ae(1,Ipg,0)=-dy0*dpsitp-(dyvdp+dy1*drop/ro)*psist  !dérivée/ln Pgaz	
       ELSE
        ae(1,1,0)=dy0*(psist*(1.d0-drop/ro)-dpsitpt)!d/ln Ptot = ln Pgaz
       ENDIF
       ae(1,1,1)=1.d0                      !der./dln Ptot
       ae(1,2,0)=-dy0*dpsitt-(dy1*drot/ro+dyvdt)*psist     !dérivée/ln T
       ae(1,3,0)=-dy0*dpsitr-(-dy1*2.d0/ay3+dyvdr)*psist   !dérivée/r 
       ae(1,4,0)=-dy0*dpsitl        
       ae(1,5,0)=-dy0*dpsitm-(dy1*(drom/ro+1.d0/ay5)+dyvdm)*psist !dér/m
       ae(1,6,0)=-dy0/teta                     !dérivée/psi
   
       dy2=dy0*psist*grad
       be(2)=y(2,1)-dy2           !dln T/dq=dln Ptot/dq grad
       IF(pturb)THEN             !avec pression turb. 7 inc.       
        ae(2,1,0)=ae(1,1,0)*grad-dy2*dgradpt    !dérivée/ln Ptot
        ae(2,1,1)=-dy2*dgradlpt      !dérivée/dln Ptot
        ae(2,Ipg,0)=ae(1,Ipg,0)*grad-dy2*dgradp  !dérivée/ln Pgaz
        ae(2,Ipg,1)=-dy2*dgradlp       !dérivée/dln Pgaz     
       ELSE                  !sans pression turb. 6 inc.
        ae(2,1,0)=ae(1,1,0)*grad-dy2*dgradpt   !dérivée/ln Ptot=lnPgaz
       ENDIF 
       ae(2,2,0)=ae(1,2,0)*grad-dy2*dgradt	!dérivée/ln T
       ae(2,2,1)=1.d0                		!dérivée/dln T
       ae(2,3,0)=ae(1,3,0)*grad-dy2*dgradr	!dérivée/r
       ae(2,4,0)=ae(1,4,0)*grad-dy2*dgradl	!dérivée/l
       ae(2,5,0)=ae(1,5,0)*grad-dy2*dgradm	!dérivée/m
       ae(2,6,0)=ae(1,6,0)*grad          	!dérivée/psi           
      
       dy3=cte13*ro*ay3**2
       be(3)=y(5,1)-dy3*psist     !d m /dq = 4 pi r**2 ro psi/teta
       IF(pturb)THEN             !avec press. turb. 7 inc. 
        ae(3,1,0)=-dy3*dpsitpt			!dérivée/ln Ptot       
        ae(3,Ipg,0)=-dy3*(psist*drop/ro+dpsitp)    !der. /ln Pgaz
       ELSE                  !sans press. turb. 6 inc.
        ae(3,1,0)=-dy3*(psist*drop/ro+dpsitpt)   !der./ln Ptot 
       ENDIF     
       ae(3,2,0)=-dy3*(psist*drot/ro+dpsitt)	!dérivée/ln T
       ae(3,3,0)=-dy3*(psist*2.d0/ay3+dpsitr)	!dérivée/r
       ae(3,4,0)=-dy3*dpsitl      
       ae(3,5,0)=-dy3*(psist*drom/ro+dpsitm)	!dérivée/m
       ae(3,5,1)=1.d0                		!dérivée/dm
       ae(3,6,0)=-dy3/teta               	!dérivée/psi
      
       dy4=cte14*ro*ay3**2*epsilon(1) !dl/dq=dm/dq epsilon
       be(4)=y(4,1)-dy4*psist         !contribution de epsilon
       IF(pturb)THEN     !avec pression turbulente 7 inconnues
        ae(4,1,0)=-dy4*dpsitpt
        ae(4,Ipg,0)=-dy4*(psist*(drop/ro+depsp)+dpsitp)!dérivée/ln Pgaz
       ELSE          !sans pression turbulente 6 inconnues
        ae(4,1,0)=-dy4*(psist*(drop/ro+depsp)+dpsitpt) !dérivée/ln Ptot
       ENDIF
       ae(4,2,0)=-dy4*(psist*(drot/ro+depst)+dpsitt)	!dérivée/ln T
       ae(4,3,0)=-dy4*(psist*2.d0/ay3+dpsitr)		!dérivée/r
       ae(4,4,0)=-dy4*dpsitl				!dérivée/l       
       ae(4,4,1)=1.d0					!dérivée/dl					!der/dl
       ae(4,5,0)=-dy4*(psist*(drom/ro+depsm)+dpsitm)	!der/m
       ae(4,6,0)=-dy4/teta								!dérivée/psi
     
       IF(dt /= 0.d0)THEN     !contribution du Tds/dt
        dy5=cte14t*ro*ay3**2 ; tdst=dy5*tdst
        be(4)=be(4)+tdst*psist        !epsilon-(tdS+d Iw2)/dt 
        dtdsp=tdst*drop/ro ; dtdst=tdst*drot/ro
        dtdsr=tdst*2.d0/ay3 ; dtdsm=tdst*drom/ro
        IF(compt == 0 .OR. dt <= dt_tds)THEN
         dtdsm=dtdsm+dy5*dtdst_m     !dérivée/m
        ELSE
         IF(kipp)THEN
          dtdsp=dtdsp+dy5*(dcpp*dtn-delta/ro*((deltap/delta
     1    -drop/ro)*dpn+pgn))/dt              !dérivée/ln Pgaz
          dtdst=dtdst+dy5*(dcpt*dtn+cp*trn        !dérivée/ln T
     1    -delta/ro*dpn*(deltat/delta-drot/ro))/dt
          dtdsm=dtdsm+dy5*((dcpm*dtn-cp*dt_tm     !dérivée/m
     1    -delta/ro*((deltam/delta-drom/ro)*dpn-dp_tm))/dt)
         ELSE
          dtdsp=dtdsp             !dérivée/ln Pgaz
     1    +dy5*(dup-pgn/ro**2*((ro-ro_t)*(1.d0-2.d0*drop/ro)+drop))/dt 
          dtdst=dtdst
     1    +dy5*(dut-pgn/ro**2*drot*(2.d0*ro_t/ro-1.d0))/dt !dérivée/ln T
          dtdsm=dtdsm+dy5*(dum-du_tm-pgn/ro**2*(-2.d0*drom*(ro-ro_t)/ro
     1    +drom-dro_tm))/dt               !dérivée/m
         ENDIF   !kipp
        ENDIF    !compt == 0 
        IF(pturb)THEN        !avec pression turbulente 7 inconnues
	 ae(4,1,0)=ae(4,1,0)+tdst*dpsitpt     !dérivée/ln Ptot
         ae(4,Ipg,0)=ae(4,Ipg,0)+dtdsp*psist !dérivée/ln Pgaz
        ELSE         !sans pression turbulente 6 inconnues
         ae(4,1,0)=ae(4,1,0)+tdst*dpsitpt+dtdsp*psist !dérivée/ln Ptot=Pgaz
        ENDIF
        ae(4,2,0)=ae(4,2,0)+dtdst*psist+tdst*dpsitt	!dérivée/ln T
        ae(4,3,0)=ae(4,3,0)+dtdsr*psist+tdst*dpsitr	!dérivée/r
        ae(4,4,0)=ae(4,4,0)+tdst*dpsitl			!dérivée/l	
        ae(4,5,0)=ae(4,5,0)+dtdsm*psist+tdst*dpsitm	!dérivée/m
        ae(4,6,0)=ae(4,6,0)+tdst/teta        		!dérivée/psi
       ENDIF      !dt

       be(5)=y(3,1)-psist     !d zeta/dq=psi/teta
       IF(pturb)THEN     !avec pression turbulente 7 inconnues
        ae(5,Ipg,0)=-dpsitp        !dérivée/ln Pgaz
        ae(5,1,0)=-dpsitpt	!dérivée/ln Ptot	
       ELSE           !sans pression turbulente 6 inconnues
        ae(5,1,0)=-dpsitpt !dérivée/ln Ptot = ln Pgaz
       ENDIF       
       ae(5,2,0)=-dpsitt     	!dérivée/ln T     
       ae(5,3,0)=-dpsitr	!dérivée/r
       ae(5,3,1)=1.d0		!dérivée/dr
       ae(5,4,0)=-dpsitl	!dérivée/l         
       ae(5,5,0)=-dpsitm	!dérivée/m
       ae(5,6,0)=-1.d0/teta	!dérivée/psi

       be(6)=y(6,1)           !dpsi/dq=0
       ae(6,6,1)=1.d0         !dérivée/dpsi
      
c Avec Pturb:
c pour les points de collocation dans les zones radiatives,
c on resout d ln Ptot - d ln Pgaz =0
c pour les points de collocation dans les zones convectives,
c on resout -Pturb -Pgaz +Ptot=0
       IF(pturb)THEN      !avec pression turbulente 7 inconnues
	IF(radiatif)THEN
         be(Ipg)=y(1,1)-y(Ipg,1)        !d lnPtot = d lnPgaz
         ae(Ipg,1,1)= 1.d0    !dérivée /ln Ptot
         ae(Ipg,Ipg,1)=-1.d0    !dérivée /ln Pgaz

        ELSE
         gam1=gam/(gam+1.d0)
         dgam1=1.d0/(gam+1.d0)**2/gam1

         dy71=cte6*gam1*delta*prn*gradad*dlpp
         dy71pt=dy71*(dgampt*dgam1+1.d0)            !der /ln Ptot
         dy71p= dy71*(dgamp*dgam1+deltap/delta+dgradadp)    !der /ln Pgaz
         dy71t= dy71*(dgamt*dgam1+deltat/delta+dgradadt)    !der /ln T
         dy71r= dy71*dgamr*dgam1                !der /r
         dy71l= dy71*dgaml*dgam1                !der /l
         dy71m= dy71*(dgamm*dgam1+deltam/delta+dgradadm)    !der /m
         dy71lpt=dy71*(dgamlpt*dgam1+dlppdpt/dlpp)          !der/dln Ptot           
         dy71lp =dy71*(dgamlp *dgam1+dlppdp /dlpp)      !der /dln Pgaz
          
         dy72=cte6*gam1*delta*prn*grad
         dy72pt=dy72*(dgampt*dgam1+1.d0+dgradpt)      !der /ln Ptot
         dy72p= dy72*(dgamp*dgam1+deltap/delta+dgradp)  !der /ln Pgaz
         dy72t= dy72*(dgamt*dgam1+deltat/delta+dgradt)  !der /ln T
         dy72r= dy72*(dgamr*dgam1+dgradr)       !der /r
         dy72l= dy72*(dgaml*dgam1+dgradl)       !der /l       
         dy72m= dy72*(dgamm*dgam1+deltam/delta+dgradm)  !der /m
         dy72lpt=dy72*(dgamlpt*dgam1+dgradlpt)    !der/dln Ptot
         dy72lp= dy72*(dgamlp *dgam1+dgradlp) !der/dln Pgaz

         be(Ipg)=dy71-dy72-pgn+prn
c        WRITE(*,2000)be(Ipg),dy71,dy72,prn,pgn,dy71-dy72,prn-pgn

         ae(Ipg,1,0)=dy71pt-dy72pt+prn       !der /ln Ptot
         ae(Ipg,2,0)=dy71t -dy72t        !der /ln T
         ae(Ipg,3,0)=dy71r -dy72r        !der /r
         ae(Ipg,4,0)=dy71l -dy72l        !der /l    
         ae(Ipg,5,0)=dy71m -dy72m        !der /m
         ae(Ipg,Ipg,0)=dy71p -dy72p-pgn        !der /ln Pgaz         
         ae(Ipg,1,1)=dy71lpt-dy72lpt     !der /dln Ptot
         ae(Ipg,Ipg,1)=dy71lp -dy72lp      !der /dln Pgaz
        ENDIF       !radiatif
       ENDIF    !Pturb
         
c      PRINT*,'be' ; WRITE(*,2000)be
c      WRITE(*,2000)cte13,ro,ay5,ay3,dy3,psist,y(3,1)
c      WRITE(*,2000)cte14,epsilon(1),am13,l13,dy4,psist,y(4,1)      
c      PRINT*,'ae'
c      DO i=1,ne
c       WRITE(*,2000)ae(i,1:ne,0)
c      ENDDO      
c      DO i=1,ne
c       WRITE(*,2000)ae(i,1:ne,1)
c      ENDDO
c      PAUSE'fin du jacobien'
 
c----------------------pour tester les dérivées analytiques---------------

c	deriv=.FALSE.
c	deriv=.TRUE.
c	deriv=cx == 1
c	deriv=cx <= 2
c	deriv=cx <= 3
c	deriv=cx <= 4
c	deriv=cx <= 5
c	deriv=cx == 5
c	1    .OR. cx == 6
c	1    .OR. cx == 50
c	2    .OR. cx == 90
c	2    .OR. cx == 189
c	2    .OR. cx == 190
c	2    .OR. cx >= 188 .AND. cx <= 190
c	2    .OR. cx >= 199 .AND. cx <= 202      
c	3    .OR. (cx > 74 .AND. cx < 76)
c	3    .OR. (cx > 43 .AND. cx < 48)
c	4    .OR. (cx > 88 .AND. cx < 92)
c	5    .OR. cx == 220
c	6    .OR. cx == 103
c	7    .OR. cx. ge. 148
c	8    .OR. cx. ge. 149
c	9    .OR. (cx. ge. 148 .AND. cx <150)
c	1    .OR. cx >= 380
c	2    .OR. ABS(fac(ip)-1.d0) > 1.d-5
c	deriv=cx == 100
c	deriv=.NOT.radiatif
c	deriv= cx >= 299 .AND. cx <= 303

       IF(.NOT. deriv)RETURN  !tests de mise au point pour dérivées
       ALLOCATE(bes(ne),dern(ne),yd(ne,0:1)) ; yd=y    
       IF(ay5 > 0.5d0)THEN
        unpdd=1.d0-dd
       ELSE
        unpdd=1.d0+dd
       ENDIF
       IF(radiatif)THEN
	PRINT*,'RADIATIF'
       ELSE
	PRINT*,'CONVECTIF'
       ENDIF 
       PRINT*,'cx,compt,radiatif,pturb,der',cx,compt,radiatif,pturb,der   
       PRINT*,'ctep,cter,fac,Q,mstar,dt,cte6,alpha'
       WRITE(*,2000)ctep,cter,fac(ip),fac(ip)*(yd(1,0)*ctep+ay3*cter),
     1 mstar,dt,cte6,alpha
       PRINT*,'be,xcoll(cx)'
       WRITE(*,2000)be(1:ne),xcoll(cx)
       PRINT*,'y / dy'
       WRITE(*,2000)yd(1:ne,0)
       WRITE(*,2000)yd(1:ne,1)
       PRINT*,'prn,pgn,trn,ln,epsilon,ro'
       WRITE(*,2000)prn,pgn,trn,y(4,0),epsilon(1),ro
       PRINT*,'xchim(i)'
       WRITE(*,2000)xchim
       PRINT*,'dxchim(i)'
       WRITE(*,2000)dxchim
       PRINT*,'despx(i)'
       WRITE(*,2000)depsx
       PRINT*,'grad,dgradx'   
       WRITE(*,2000)grad,dgradx   
       PRINT*,'dgradpt,dgradt,dgradr,dgradl,dgradm,dgradp,dgradlpt,dgradlp'
       WRITE(*,2000)dgradpt,dgradt,dgradr,dgradl,dgradm,dgradp,dgradlpt,dgradlp
       PRINT*,'dgradlpp,gradad,dgradadp,dgradadt,dgradadx,t_t,p_t'
       WRITE(*,2000)dgradlpp,gradad,dgradadp,dgradadt,dgradadx,t_t,p_t
       PRINT*,'psist,teta,dpsitp,dpsitt,dpsitpt,dpsitr,dpsitm'
       WRITE(*,2000)psist,teta,dpsitp,dpsitt,dpsitpt,dpsitr,dpsitm
       PRINT*,'dtetap,dtetapt,dtetat,dtetar,dtetam,tdst,dp_tm,dt_tm'
       WRITE(*,2000)dtetap,dtetapt,dtetat,dtetar,dtetam,tdst,dp_tm,dt_tm
       PRINT*,'depsp,depst,depsm,tdst,dtdsp,dtdst,dtdsm'
       WRITE(*,2000)depsp*epsilon(1),depst*epsilon(1),depsm*epsilon(1),tdst,
     1 dtdsp,dtdst,dtdsm
       PRINT*,'tdst'
       WRITE(*,2000)tdst
       PRINT*,'dlpp,dlppdpt,dlppdp,dy71,dy72'
       WRITE(*,2000)dlpp,dlppdpt,dlppdp,dy71,dy72
       PRINT*,'dy71pt,dy71t,dy71r,dy71l,dy71m,dy71p,dy71lpt,dy71lp'
       WRITE(*,2000)dy71pt,dy71t,dy71r,dy71l,dy71m,dy71p,dy71lpt,dy71lp
       PRINT*,'dy72pt,dy72t,dy72r,dy72l,dy72m,dy72p,dy72lpt,dy72lp'
       WRITE(*,2000)dy72pt,dy72t,dy72r,dy72l,dy72m,dy72p,dy72lpt,dy72lp   
       PRINT*,'u,dup,dut,dux,dum,u_t,du_tm,dtdst_m'
       WRITE(*,2000)u,dup,dut,dux,dum,u_t,du_tm,dtdst_m
       PRINT*,'ro,drop,drot,drox,drom,ro_t,dro_tm'
       WRITE(*,2000)ro,drop,drot,drox,drom,ro_t,dro_tm
       PRINT*,'gam,dgamx,dgamlpp' 
       WRITE(*,2000)gam,dgamx,dgamlpp 
       PRINT*,'dgampt,dgamt,dgamr,dgaml,dgamm,dgamp,dgamlpt,dgamlp'
       WRITE(*,2000)dgampt,dgamt,dgamr,dgaml,dgamm,dgamp,dgamlpt,dgamlp
              
       psist0=psist ; teta0=teta ; grad0=grad ; eps0=epsilon(1)
       ro0=ro ; u0=u ; gam0=gam ; x00=xchim(1) ; dy710=dy71 ; dy720=dy72
       p_t0=p_t ; t_t0=t_t ; ro_t0=ro_t ; u_t0=u_t
       tds0=tdst
      
c       DO jd=0,1  !pour tester les dérivées/ y(i,0) et y(i,1)
       DO jd=0,0  !pour tester les dérivées/ y(i,0)       
        DO id=1,ne
         stor0=yd(id,jd) ; stor=stor0*unpdd
         IF(ABS(stor) < dd)stor=sign(dd,stor)
         dstor=stor-stor0 ; yd(id,jd)=stor
         ay3=ABS(yd(3,0)) !r
         ay5=ABS(yd(5,0)) !m
         am23=ay5**(2.d0/3.d0)!m**2/3
         am13=SQRT(am23)  !m**1/3
         prn=EXP(yd(1,0)) !Ptot cgs
         IF(pturb)THEN    !avec pression turbulente 7 inconnues
          pgn=EXP(yd(Ipg,0))    !Pgaz cgs
         ELSE         !sans pression turbulente 6 inconnues
          pgn=prn
         ENDIF
         trn=EXP(yd(2,0))     !temperature K
c        PAUSE'avant bsp1dn1'
       
         CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
     1    knotc,.TRUE.,MIN(am23,mc(n_ch)),l,xchim,dxchim)
         IF(no_croiss)PRINT*,'Pb. at 6 in static_r'
         xchim=ABS(xchim) ; dxchim=dxchim*2.d0/3.d0/am13
         
         IF(pturb .AND. der)THEN !der=t on tient compte de dln Pgaz/dln Ptot
          dlpp=yd(Ipg,1)/yd(1,1)        !dlpp=dln Pgaz/dln Ptot
         ELSE
          dlpp=1.d0   !der=f on ne tient pas compte de dln Pgaz/dln Ptot
         ENDIF

         CALL thermo(prn,pgn,trn,ay5,ABS(yd(4,0)),ay3,dlpp,xchim,dxchim,
     1   ro,drop,drot,drox,u,dup,dut,dux,
     2   grad,dgradpt,dgradp,dgradt,dgradx,dgradm,dgradl,dgradr,dgradlpp,
     3   gam,dgampt,dgamp,dgamt,dgamx,dgamm,dgaml,dgamr,dgamlpp,     
     4   epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
     5   delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
     6   gradad,dgradadp,dgradadt,dgradadx,
     7   hp,dhppt,dhpp,dhpt,dhpx,dhpr,dhpm,grad_mu,
     8   gradrad,alfa,beta,gamma1,radiatif)

c la rotation
 	 SELECT CASE(Krot)
	 CASE(0,1,2)
	  w=wrot ; dwdm=0.d0
	 CASE(3,4,5)
	  CALL bsp1dn(nrot,rota,mrot,mrot,n_rot,ord_rot,knotr,.TRUE.,
	1   am23,l,frot,dfrot)
	  IF(no_croiss)PRINT*,'Pb. at 7 in static_r'
	  w=frot(1) ; dwdm=dfrot(1)*2.d0/3.d0/am13
	 END SELECT

         v=ro*ay3/prn*w**2	!rotation

         teta=fac(ip)*(ro*(ay5/ay3**2/prn*(cte21+cte22*grad)+cte25*v
     1   +ay3**2*(cte23+cte24*epsilon(1)))+cter)
         psist=yd(6,0)/teta	!psi/teta

         IF(dt /= 0.d0)THEN
c         PRINT*,dt
          IF(compt == 0)THEN      !interpolation du TdS au temps t
           CALL bsp1dn(1,tds_t,x_tds_t,xt_tds_t,n_tds_t,m_tds,
     1     knot_tds_t,.TRUE.,MIN(ay5,x_tds_t(n_tds_t)),l,ftds,dftds)
           IF(no_croiss)PRINT*,'Pb. at 8 in static_r'
           tdst=ftds(1) ; dtdst_m=dftds(1)
c          PRINT*,'tdst,l',tdst,dtdst_m,l

          ELSE                !TdS au temps t+dt
           CALL bsp1dn(1,old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm,.TRUE.,
     1     MIN(ay5,x_ptm(n_ptm)),l,fm,dfm)
           IF(no_croiss)PRINT*,'Pb. at 9 in static_r'
           mk=fm(1)          
           CALL inter('m23',bp_t,q_t,qt_t,n_qs_t,knot_t,mk,fqs,dfqs,r2_t,m23_t)
           IF(pturb)THEN      !avec pression turbulente 7 inconnues
            p_t=EXP(fqs(Ipg))
           ELSE           !sans pression turbulente 6 inconnues
            p_t=EXP(fqs(1))  
           ENDIF
           t_t=EXP(fqs(2))
          
           IF(kipp)THEN   !approximation de Kippenhahan
            dpn=pgn-p_t ; dtn=trn-t_t ; tdst=(cp*dtn-delta/ro*dpn)/dt
           ELSE       !TdS=dU+PdV
            mk=mk**(2.d0/3.d0)
            CALL bsp1dn(nchim,chim_t,mc_t,mct_t,n_ch_t,m_ch,
     1      knotc_t,.TRUE.,MIN(mk,mc_t(n_ch_t)),l,xchim_t,dxchim_t)
            IF(no_croiss)PRINT*,'Pb. at 10 in static_r'
            xchim_t=ABS(xchim_t) ; dxchim_t=dxchim_t*2.d0/3.d0/am13
            CALL chim_gram(xchim_t,dxchim_t)
            CALL etat(p_t,t_t,xchim_t,.FALSE.,
     1      ro_t,drop_t,drot_t,drox_t,u_t,dup_t,dut_t,dux_t,
     2      delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
     3      gradad,dgradadp,dgradadt,dgradadx,  
     4      alfa,beta,gamma1)
            tdst=(u-u_t-pgn/ro**2*(ro-ro_t))/dt   !+Tds
           ENDIF      !kipp
          ENDIF       !compt=0
         ENDIF        !dt/=0

         dy1=cte11*ay5/ay3**2/prn*ro
         dyv=cte15*v              !rotation    
         dy0=dy1+dyv
         bes(1)=yd(1,1)-dy0*psist
         bes(2)=yd(2,1)-dy0*psist*grad  !dln T/dq=dln p/dq grad psi/teta
         dy3=cte13*ro*ay3**2
         bes(3)=yd(5,1)-dy3*psist     !dl/dq=epsilon am13 psi/teta
         bes(4)=yd(4,1)-cte14*ro*ay3**2*epsilon(1)*psist
         IF(dt /= 0.d0)THEN
          dy5=cte14t*ro*ay3**2 ; tdst=dy5*tdst
          bes(4)=bes(4)+tdst*psist        !epsilon-tdS
         ENDIF                !dt

         bes(5)=yd(3,1)-psist
	 bes(6)=yd(6,1)
       
         IF(pturb)THEN
          IF(radiatif)THEN
           bes(Ipg)=yd(1,1)-yd(Ipg,1)
          ELSE
           gam1=gam/(gam+1.d0) ; dy71=cte6*gam1*delta*prn*gradad*dlpp
           dy72=cte6*gam1*delta*prn*grad ; bes(Ipg)=dy71-dy72-pgn+prn
          ENDIF       !radiatif
         ENDIF        !Pturb
 
         dern=(bes-be)/dstor ; PRINT*
         SELECT CASE(jd)
         CASE(0)   
          PRINT*,'dérivée / ',variable(id)
         CASE(1)
          PRINT*,'dérivée / d',variable(id)
         END SELECT       
         WRITE(*,2000)(dern(j),ae(j,id,jd),j=1,4)
         WRITE(*,2000)(dern(j),ae(j,id,jd),j=5,ne)          
         yd(id,jd)=stor0    
         PRINT*,'psist,grad,teta,eps,ro,ro_t,u,u_t,tdst,p_t,t_t,X,gam,dy71,dy72'   
         WRITE(*,2000)(psist-psist0)/dstor,(grad-grad0)/dstor/grad0,
	1 (teta-teta0)/dstor,(epsilon(1)-eps0)/dstor,
	2 (ro-ro0)/dstor,(ro_t-ro_t0)/dstor,(u-u0)/dstor,
	3 (u_t-u_t0)/dstor,(tdst-tds0)/dstor,
	4 (p_t-p_t0)/dstor,(t_t-t_t0)/dstor,(xchim(1)-x00)/dstor,
	5 (gam-gam0)/dstor,(dy71-dy710)/dstor,(dy72-dy720)/dstor
c il y a dpsitp=d psist / d lnP et dpsitpt=d psist / d lnPtot	
	
        ENDDO !jd
        PRINT*,'cx=',cx,ip,fac(ip) ; PAUSE'fin dérivée'
       ENDDO !id
       DEALLOCATE(bes,dern,yd)
       RETURN
      
      CASE(2)     !les limites

c conditions aux limites : résidu be(1) au point limite li
c      li=1 : limite sur r au centre
c      li=2 : limite sur l au centre
c      li=3 : limite sur m au centre
c      li=4 : limite sur Ptot au raccord
c      li=5 : limite sur T au raccord
c      li=6 : limite sur m au raccord
c      li=7 : limite sur Pgaz au raccord avec Pturb

c deriv=.FALSE. !; deriv=.TRUE.    
       IF(deriv)THEN
        ALLOCATE(yd(ne,0:1)) ; yd=y ; unpdd=1.d0-dd
       ENDIF
       SELECT CASE(li)
       CASE(1)        !au centre en r
        be(1)=y(3,0)      !en q=1 r=0
        ae(1,3,0)=1.d0    !dérivée/r
c       PRINT*,'limite au centre',li
c       WRITE(*,2000)y(1:ne,0),be(1)

       CASE(2)        !au centre en l
        be(1)=y(4,0)      !en q=1 l=0
        ae(1,4,0)=1.d0    !dérivée/l
c       PRINT*,'limite au centre',li
c       WRITE(*,2000)y(1:ne,0),be(1)

       CASE(3)        !au centre en m
        be(1)=y(5,0)      !en q=1 m=0
        ae(1,5,0)=1.d0    !dérivée/m
c       PRINT*,'limite au centre',li
c       WRITE(*,2000)(y(i),i=1,ne),be(1)
c       PAUSE'au centre'

       CASE(4)        !Ptot au raccord
        am23=y(5,0)**(2.d0/3.d0)
        CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
     1  knotc,.TRUE.,MIN(am23,mc(n_ch)),l,xchim,dxchim)
        IF(no_croiss)PRINT*,'Pb. at 11 in static_r'
        CALL chim_gram(xchim,dxchim)
        CALL atm(.FALSE.,y(4,0),y(3,0),xchim,ptext,dptdl,dptdr,
     1	text,dtdl,dtdr,mext,dmdl,dmdr,pext,dpdl,dpdr,teff)
        be(1)=y(1,0)-LOG(ptext)		!condition sur Ptot au raccord
        ae(1,1,0)=1.d0			!dérivée/ln Ptot
        ae(1,3,0)=-dptdr/ptext		!dérivée/r
        ae(1,4,0)=-dptdl/ptext		!dérivée/l
        IF(deriv)THEN			!test dérivée
         PRINT*,'limite pext',li ; WRITE(*,2000)y(1,0),LOG(ptext)  
c        WRITE(*,2000)y(1:ne,0)
c        PRINT*,'exterieur lnp,ln ptext,be(1),ray,xchim'
c        WRITE(*,2000)y(1,0),LOG(ptext),be(1),ray,ln,xchim
         ptext0=ptext ; text0=text ; mext0=mext ; pext0=pext
         DO i=3,4
          stor0=yd(i,0) ; stor=stor0*unpdd
          IF(stor == 0.d0)stor=dd
          dstor=stor-stor0 ; yd(i,0)=stor
          CALL atm(.FALSE.,yd(4,0),yd(3,0),xchim,ptext,dptdl,dptdr,
     1    text,dtdl,dtdr,mext,dmdl,dmdr,pext,dpdl,dpdr,teff)
c         PRINT*,LOG(ptext),yd(1,0),yd(1,0)-LOG(ptext),be(1)
          yd(i,0)=stor0
          PRINT*,'dérivée en li=4 / ',variable(i)
          WRITE(*,2000)(yd(1,0)-LOG(ptext)-be(1))/dstor,ae(1,i,0)
          yd(i,0)=stor0   
         ENDDO    !i
         ptext=ptext0 ; text=text0 ; mext=mext0 ; pext=pext0
         PAUSE'test deriv'   
        ENDIF     !deriv

       CASE(5)      !T au raccord
        be(1)=y(2,0)-LOG(text)		!condition sur T au raccord
        ae(1,2,0)=1.d0			!dérivée/ln T
        ae(1,3,0)=-dtdr/text		!dérivée/r
        ae(1,4,0)=-dtdl/text		!dérivée/l
        IF(deriv)THEN     !test de derivation
         PRINT*,'limite Text',li ; WRITE(*,2000)y(2,0),LOG(text)
c        WRITE(*,2000)y(1:ne,0)
c        PRINT*,'y(2,0),LOG(text),be(1)'
c        WRITE(*,2000)y(2,0),LOG(text),be(1)
         pext0=pext ; text0=text ; mext0=mext ; pext0=pext
         DO i=3,4
          stor0=yd(i,0) ; stor=stor0*unpdd
          IF(stor == 0.d0)stor=dd
          dstor=stor-stor0 ; yd(i,0)=stor
          CALL atm(.FALSE.,yd(4,0),yd(3,0),xchim,ptext,dptdl,dptdr,
     1    text,dtdl,dtdr,mext,dmdl,dmdr,pext,dpdl,dpdr,teff)
           PRINT*,'dérivée pour li=5 / ',variable(i)
           WRITE(*,2000)(yd(2,0)-LOG(text)-be(1))/dstor,ae(1,i,0)
           yd(i,0)=stor0       
         ENDDO    !i
         ptext=ptext0 ; text=text0 ; mext=mext0 ; pext=pext0
         PAUSE'test deriv'       
        ENDIF     !deriv
           
       CASE(6)    !condition sur M au raccord
        be(1)=y(5,0)-mext	!en q=n_qs, m=mext(r,l)
        ae(1,3,0)=-dmdr		!dérivée/r
        ae(1,4,0)=-dmdl		!dérivée/l
        ae(1,5,0)=1.d0		!dérivée/m
c       WRITE(*,2000)y(5,0)
        IF(deriv)THEN     !test de derivation
         PRINT*,'limite mext',li ; WRITE(*,2000)y(5,0),mext
c        WRITE(*,2000)y(1:ne,0)
c        PRINT*,'limite y(5),mext',li
c        PRINT*,'y(5,0),mext,be(1)'
c        WRITE(*,2000)y(5,0),mext,be(1)
         ptext0=ptext ; text0=text ; mext0=mext ; pext0=pext
         DO i=3,4
          stor0=yd(i,0) ; stor=stor0*unpdd
          IF(stor == 0.d0)stor=dd
          dstor=stor-stor0 ; yd(i,0)=stor
          CALL atm(.FALSE.,yd(4,0),yd(3,0),xchim,ptext,dptdl,dptdr,
     1    text,dtdl,dtdr,mext,dmdl,dmdr,pext,dpdl,dpdr,teff)
          PRINT*,'dérivée pour li=6 / ',variable(i)
          WRITE(*,2000)(yd(5,0)-mext-be(1))/dstor,ae(1,i,0)
          yd(i,0)=stor0       
         ENDDO    !i    
         ptext=ptext0 ; text=text0 ; mext=mext0 ; pext=pext0
         PAUSE'test deriv'   
        ENDIF     !deriv
	
       CASE(7)    !condition sur Pgaz au raccord
        be(1)=y(Ipg,0)-LOG(pext)        !en q=n
        ae(1,3,0)=-dpdr/pext		!dérivée/r
        ae(1,4,0)=-dpdl/pext		!dérivée/l
        ae(1,Ipg,0)=1.d0		!dérivée/ln Pgaz
        IF(deriv)THEN     !test de derivation
         PRINT*,'limite pext, uint',li ; WRITE(*,2000)y(Ipg,0),pext
c        WRITE(*,2000)y(1:ne,0)
c        PRINT*,'limite y(Ipg)',li
c        PRINT*,'y(Ipg,0),pext,be(1)'
c        WRITE(*,2000)y(Ipg,0),pext,be(1)
         ptext0=ptext ; text0=text ; mext0=mext ; pext0=pext
         DO i=3,4
          stor0=yd(i,0) ; stor=stor0*unpdd
          IF(stor == 0.d0)stor=dd
          dstor=stor-stor0 ; yd(i,0)=stor
          CALL atm(.FALSE.,yd(4,0),yd(3,0),xchim,ptext,dptdl,dptdr,
     1    text,dtdl,dtdr,mext,dmdl,dmdr,pext,dpdl,dpdr,teff)
          PRINT*,'dérivée pour li=7 / ',variable(i)
          WRITE(*,2000)(yd(Ipg,0)-LOG(pext)-be(1))/dstor,ae(1,i,0)
          yd(i,0)=stor0       
         ENDDO    !i    
         ptext=ptext0 ; text=text0 ; mext=mext0 ; pext=pext0
         PAUSE'test deriv'       
        ENDIF     !deriv

       CASE DEFAULT       
        WRITE(*,"(' static_r erreur sur li =',i3,', ne=',i3)")li,ne
        PRINT*,'ARRET'	
        STOP      
       END SELECT
       IF(deriv)DEALLOCATE(yd)
       
c      PRINT*,'li=',li ; PAUSE'à la fin'
      
      CASE DEFAULT   
       PRINT*,' static_r erreur sur fait = 1, ou 2, fait=',fait
       PRINT*,'ARRET' ; STOP  
      END SELECT       !fait
      
      RETURN
         
      END SUBROUTINE static_r
