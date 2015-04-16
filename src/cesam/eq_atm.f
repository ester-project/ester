
c************************************************************************

	SUBROUTINE eq_atm(fait,li,xchim,cx,y,be,ae,r_rac,l_rac,xcoll_atm)

c routine private	du module mod_atm

c équations de l'atmosphère pour la méthode de collocation
c le rayon, la masse sont mesures en tau23 fixe au point n23_atm
c correction pour assurer la continuite du gradient au point
c de raccord avec l'enveloppe
	 
c avec pression turbulente 8 inconnues 
c sans pression turbulente 7 inconnues, Ptot=Pgaz

c modifs :
c 19 11 99 : suppression de nh1, nhe1, nhe2, lamb
c 04 07 01 : mise en place de rotation uniforme avec conservation du moment
c angulaire
c 16 08 01 : F95

c Auteur: P. Morel, Département J.D. Cassini, O.C.A., CESAM2k

c entrées:
c	fait=1 : point courant	
c	     2 : limite
c	mstar: masse totale, temps t+dt, avec perte de masse
c	wrot : vitesse angulaire solide
c	li : indice de la limite
c	xchim : composition chimique / gramme
c	cx : indice du point de collocation
c	r_rac :rayon au fond de l'atmosphère
c	l_rac : luminosité=cte
c	y(ne_atm,0), y(ne_atm,1) : variables, dérivées

c sorties
c	be, ae : équations et dérivées

c-----------------------------------------------------------------------

	USE mod_donnees, ONLY : alpha, aradia, clight, cpturb, g, Krot,
	1 lim_ro, lsol, msol, pi, pturb, rsol
	USE mod_etat, ONLY : etat	
	USE mod_kind
	USE mod_opa, ONLY : opa
	USE mod_variables, ONLY: dim_rot, mstar, rota, wrot
	
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:,0:) :: y
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim, xcoll_atm
	REAL (kind=dp), INTENT(in) :: r_rac, l_rac
	INTEGER, INTENT(in) :: fait, cx, li		
	REAL (kind=dp), INTENT(out), DIMENSION(:,:,0:) :: ae	
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: be

	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: yb
	REAL (kind=dp), SAVE, ALLOCATABLE, DIMENSION(:) :: Lxchim
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: bes, dern
	
	REAL (kind=dp), PARAMETER :: dd=0.00001d0, unpdd=1.d0+dd
	
	REAL (kind=dp), SAVE ::	cte1, cte2, cte3, cte4, cte5, cte6	
	
	REAL (kind=dp) ::ro, drop, drot, drox, tau, u, dup, dut, dux,
	1 grav, teff, dtsdtau, dtsdteff, dtsdg, dtbr23, trn, kap0, w, ro0,
	2 kap, dkapt, tb, dkapp, dkapro, dkapx, dgravre, dteffre,
	3 f_tau, df_tau, d2f_tau, vt, acc_c,
	4 dy1, dy2, dy3, dy5, prn, stor, stor0,
	5 gradrad, dstor, dlpp, dgradlpt, dgradlp,
	6 dgradlpp, dlppdp, pgn, dgradpt, hp, dhppt, dhpp, dhpt, dhpr,
	7 dhpm, dy81, dy81pt, dy81t, dy81r, dy81rs, dy81m, dy81tau,
	8 dy81p, dy81lpt, dy81lp, dy82, dy82pt, dy82t, dy82r, dy82rs,
	9 dy82m, dy82tau, dy82p, dy82lpt, dy82lp, gam, dgampt, dgamp,
	1 dgamt, dgamr, dgamrs, dgamm, dgamtau, dgamlpp, dgamlpt, dgamlp,
	2 gam1, dgam1, gam0, dlppdpt, delfi, ddelfi, ro_ext, dro_grav,
	3 dro_teff, grad, dgradp, dgradt, dgradr, dgradrs, dgradm,
	4 dgradtau, grad0, gradad, dgradadp, dgradadt, dgradadx,
	5 alfa, beta, gamma1, delta, deltap, deltat, deltax, cp,
	6 dcpp, dcpt, dcpx
	
	INTEGER :: i, j
		
	CHARACTER (len=7), PARAMETER, DIMENSION(9) ::
	1 variable=(/'ln Ptot','ln T   ','R      ','Rstar  ',
	2 'M      ','lntau23','ln tau ','ln rho ', 'ln Pgaz'/)
	
	LOGICAL, SAVE :: der, init=.TRUE.	
	LOGICAL :: deriv, radiatif
	
c-----------------------------------------------------------------------

2000	FORMAT(8es10.3)

c	PRINT*,'entrée eqatm, tau_min,tau_max,n_atm,n23_atm,ne_atm',
c	1 tau_min,tau_max,n_atm,n23_atm,ne_atm
	
	IF(init)THEN
	 init=.FALSE.
	 cte1=g*msol/rsol**2 ; cte2=lsol/pi/rsol**2/aradia/clight
	 cte3=4.d0*pi*rsol**2/msol ; cte4=2.d0/3.d0*rsol
	 cte5=-1.d0/rsol ; cte6=abs(cpturb)*alpha**2/8.d0 !pour la Pturb
	 der=cpturb < 0.d0	!der=.TRUE. avec dln Pgaz/dln Ptot
	 ALLOCATE(Lxchim(SIZE(xchim,1)))
	ENDIF		

	ae=0.d0 ; be=0.d0
	
	SELECT CASE(fait)	
	CASE(1)
	
c------------------	équations au point courant	------------------

c		y(1,0)=ln Ptot		y(1,1)=(ln Ptot)'
c		y(2,0)=ln T		y(2,1)=(ln T)'
c		y(3,0)=R/Rsol		y(3,1)=(R/Rsol)'
c		y(4,0)=R23/Rsol		y(4,1)=(R23/Rsol)'
c		y(5,0)=M/Msol		y(5,1)=(M/Msol)'
c		y(6,0)=ln tau23		y(6,1)=(ln tau23)'
c		y(7,0)=ln tau		y(7,1)=(ln tau)'
c avec Pturb	y(Ipgt,0)=ln Pgaz	y(Ipgt,1)=(ln Pgaz)'

c	 ae(eq,var,der)=dérivée de la eq-ième
c	 équation par rapport à la der-ième dérivée de la var-ième variable

	 prn=EXP(y(1,0))	!pression totale cgs
	 trn=EXP(y(2,0))	!temperature cgs
	 tau=EXP(y(7,0))	!profondeur optique
	 IF(pturb)THEN		!avec pression turbulente 8 inconnues
	  pgn=EXP(y(Ipgt,0))	!pression gaz+rad cgs	 
	 ELSE			!sans pression turbulente 7 inconnues
	  pgn=prn
	 ENDIF
c	 WRITE(*,2000)y(:,0) ; WRITE(*,2000)y(:,1) ; PAUSE'y dans eq_atm'
		
c température effective	
	 teff=(cte2*l_rac/y(4,0)**2)**0.25d0
	 
c gravité 
	 grav=cte1*mstar/y(4,0)**2
	 
c vitesse angulaire
	 SELECT CASE(Krot)
	 CASE(0,1,2)
	  w=wrot
	 CASE(3,4,5)
	  w=rota(1,dim_rot)
	 END SELECT
	 w=w**2
	 
c accélération centrifuge	 
	 acc_c=cte4*y(4,0)*w
	 
c gravité effective	 
	 grav=grav-acc_c

c la fonction g(tau) et dérivées		
	 CALL tdetau(tau,teff,grav,tb,dtsdtau,dtsdteff,dtsdg,
	1 ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)	

c	 WRITE(*,2000)tau,l_rac,y(4,0),teff,tb,f_tau,df_tau,d2f_tau
c	 PAUSE'après t(tau)'

c les gradients, der permet d'éliminer dln Pgaz/dln Ptot
	 IF(pturb .and. der)THEN
	  dlpp=y(Ipgt,1)/y(1,1)		!dlpp=dln Pgaz/dln Ptot
	  dlppdpt=-dlpp/y(1,1)		!der dlpp /dln Ptot
	  dlppdp=dlpp/y(Ipgt,1)		!der dlpp /dln Pgaz
	 ELSE
	  dlpp=1.d0 ; dlppdpt=0.d0 ; dlppdp=0.d0	 	 
	 ENDIF
	
	 Lxchim=xchim
	 CALL thermo_atm(prn,pgn,trn,Lxchim,y(5,0),l_rac,y(3,0),dlpp,
	1 tau,df_tau,d2f_tau,y(4,0),ro,drop,drot,kap,dkapp,dkapt,gradad,
	2 dgradadp,dgradadt,grad,dgradpt,dgradp,dgradt,dgradr,dgradrs,
	3 dgradm,dgradtau,dgradlpp,
	4 gam,dgampt,dgamp,dgamt,dgamr,dgamrs,dgamm,dgamtau,dgamlpp,	
	5 hp,dhppt,dhpp,dhpt,dhpr,dhpm,delta,deltap,deltat,
	6 cp,gradrad,alfa,beta,gamma1,radiatif,.TRUE.)
	
	 drop=drop*pgn/ro	!der ro /ln Pgaz à ro près		
	 drot=drot*trn/ro	!der ro /ln T à ro près
		
	 dkapp=dkapp*pgn/kap	!der kap/ln Pgaz à kap près
	 dkapt=dkapt*trn/kap	!der kap/ln t à kap près	

	 dgradadp=dgradadp*pgn/gradad	!der gradad /ln Pgas à gradad près
	 dgradadt=dgradadt*trn/gradad	!der gradad /ln t à gradad près
	
	 dgradpt=dgradpt*prn/grad	!der grad /ln Ptot à grad près
	 dgradp= dgradp *pgn/grad	!der grad /ln Pgaz à grad près	
	 dgradt=  dgradt*trn/grad	!der grad /ln t à grad près 
	 dgradr=dgradr/grad		!der grad /r à grad près 
	 dgradrs=dgradrs/grad		!der grad /rstar à grad près 	
	 dgradm=dgradm/grad		!der grad /m  à grad près	
	 dgradtau=dgradtau*tau/grad	!der grad /ln tau  à grad près
	 dgradlpt=dgradlpp*dlppdpt/grad	!der grad /dln Ptot à grad près 
	 dgradlp= dgradlpp*dlppdp/grad	!der grad /dln Pgaz à grad près
	
	 dgampt=dgampt*prn		!der gam /ln Ptot
	 dgamp=dgamp*pgn		!der gam /ln Pgaz 	
	 dgamt=dgamt*trn		!der gam /ln T 
	 dgamtau=dgamtau*tau		!der gam /ln tau 
	 dgamlpt=dgamlpp*dlppdpt	!der gam /dln Ptot 
	 dgamlp= dgamlpp*dlppdp		!der gam /dln Pgaz
	
	 dhppt=dhppt*prn/hp	!der hp /ln Ptot à hp près
	 dhpp= dhpp *pgn/hp	!der hp /ln Pgaz à hp près	
	 dhpt= dhpt *trn/hp	!der hp /ln T à hp près
	 dhpr=dhpr/hp		!der hp /r à hp près	
	 dhpm=dhpm/hp		!der hp /m à hp près
	
	 deltap=deltap*pgn/delta	!der delta /ln Pgas à delta près
	 deltat=deltat*trn/delta	!der delta /ln T à delta près

c Definition de la fonction fi
         IF(xcoll_atm(cx) <= n23_atm) THEN
          delfi=(y(6,0)-ltauf)*delfim  
          ddelfi=delfim   !dérivée fi /y(6,0)
         ELSE
          delfi=(y(6,0)-ltaue)*delfip
          ddelfi=delfip   !dérivée fi /y(6,0)
         ENDIF
	 ddelfi=ddelfi/delfi		!a delfi près

c équations
	 dy1=(cte1*y(5,0)/y(3,0)**2-cte4*w*y(3,0))*delfi/kap
	1 *EXP(y(7,0)-y(1,0))
	 be(1)=y(1,1)-dy1		!d log p/d tau=gm/kr**2/p
	 ae(1,1,0)=dy1			!dérivée /ln Ptot
	 ae(1,1,1)=1.d0			!dérivée /dln Ptot
	 ae(1,2,0)=dy1*dkapt		!dérivée /ln T
	 ae(1,3,0)=(2.*cte1*y(5,0)/y(3,0)**3+cte4*w)*
	1 delfi/kap*EXP(y(7,0)-y(1,0))		!dérivée /r
	 ae(1,5,0)=-cte1/y(3,0)**2*delfi/kap*EXP(y(7,0)-y(1,0))	!dérivée /m
         ae(1,6,0)=-ddelfi*dy1		!derive/ln tau23
         ae(1,7,0)=-dy1			!dérivée /ln tau
	 IF(pturb)THEN		!avec pression turbulente 8 inconnues
 	  ae(1,Ipgt,0)=dy1*dkapp		!dérivée /ln Pgaz 
	 ELSE			!sans pression turbulente 7 inconnues
	  ae(1,1,0)=ae(1,1,0)+dy1*dkapp		!dérivée /ln Ptot
	 ENDIF

	 dy2=dy1*grad
	 be(2)=y(2,1)-dy2
	 ae(2,1,0)=ae(1,1,0)*grad-dy2*dgradpt	!dérivée /ln Ptot
	 ae(2,2,0)=ae(1,2,0)*grad-dy2*dgradt	!dérivée /ln T
	 ae(2,2,1)=1.d0				!dérivée /dln T
	 ae(2,3,0)=ae(1,3,0)*grad-dy2*dgradr	!dérivée /r
	 ae(2,4,0)=-dy2*dgradrs			!dérivée /r23
	 ae(2,5,0)=ae(1,5,0)*grad-dy2*dgradm	!dérivée /m
         ae(2,6,0)=ae(1,6,0)*grad   		!d/ln tau23	
	 ae(2,7,0)=ae(1,7,0)*grad-dy2*dgradtau	!dérivée /lntau
	 IF(pturb)THEN		!avec pression turbulente 8 inconnues
	  ae(2,1,1)=-dy2*dgradlpt		!der/dln Ptot 	        
          ae(2,Ipgt,0)=ae(1,Ipgt,0)*grad-dy2*dgradp	!dérivée /ln Pgaz
	  ae(2,Ipgt,1)=-dy2*dgradlp		!dérivée /dln Pgaz	
	 ELSE			!sans pression turbulente 7 inconnues
	  ae(2,1,0)=ae(2,1,0)-dy2*dgradp	!dérivée /ln Ptot 
	 ENDIF
	
	 dy3=cte5/kap/ro*tau*delfi
	 be(3)=y(3,1)-dy3
	 ae(3,2,0)=dy3*(drot+dkapt)	!dérivée /ln t
	 ae(3,3,1)=1.d0			!dérivée /d(R/Rsol)
	 ae(3,6,0)=-dy3*ddelfi		!dérivée /ln tau23
	 ae(3,7,0)=-dy3			!dérivée /ln tau
	 IF(pturb)THEN		!avec pression turbulente 8 inconnues
	  ae(3,Ipgt,0)=dy3*(drop+dkapp)	!dérivée /ln Pgaz  
	 ELSE			!sans pression turbulente 7 inconnues
	  ae(3,1,0)=dy3*(drop+dkapp)	!dérivée /ln Ptot  	 
	 ENDIF

         be(4)=y(4,1)
         ae(4,4,1)=1.d0			!dérivée /d(r23/Rsol)

	 dy5=-cte3*y(3,0)**2*tau*delfi/kap
	 be(5)=y(5,1)-dy5
	 ae(5,2,0)=dy5*dkapt		!dérivée /ln t
	 ae(5,3,0)=-2.d0*dy5/y(3,0)	!dérivée /R/Rsol
	 ae(5,5,1)=1.d0			!dérivée /d(M/Msol)
	 ae(5,6,0)=-dy5*ddelfi		!dérivée /ln tau23
	 ae(5,7,0)=-dy5			!dérivée /ln tau
 	 IF(pturb)THEN		!avec pression turbulente 8 inconnues
          ae(5,Ipgt,0)=dy5*dkapp		!dérivée /ln Pgaz
	 ELSE			!sans pression turbulente 7 inconnues
          ae(5,1,0)=dy5*dkapp		!dérivée /ln Ptot	 
	 ENDIF
        
         be(6)=y(6,1)
	 ae(6,6,1)=1.d0			!dérivée /d ln(tau23)

         be(7)=y(7,1)-delfi
         ae(7,6,0)=-ddelfi*delfi	!dérivée /ln tau23
         ae(7,7,1)=1.d0			!dérivée /d ln(tau)
    
	 IF(pturb)THEN		!avec pression turbulente 8 inconnues 
          IF(radiatif)THEN
           be(Ipgt)=y(1,1)-y(Ipgt,1)
           ae(Ipgt,1,1)=1.d0	!dérivée /ln Ptot
           ae(Ipgt,Ipgt,1)=-1.d0	!dérivée /ln Pgaz       
	  ELSE
           gam1=gam/(gam+1.d0)
           dgam1=1.d0/(gam+1.d0)**2/gam1
                 
           dy81=cte6*gam1*delta*prn*gradad*dlpp
           dy81pt=dy81*(1.d0+dgampt*dgam1)		!der /ln Ptot        
           dy81p =dy81*(dgradadp+deltap+dgamp*dgam1) 	!der /ln Pgaz
           dy81t=dy81* (dgradadt+deltat+dgamt*dgam1) 	!der /ln T
           dy81r=dy81 *dgamr *dgam1			!der /r
           dy81rs=dy81*dgamrs*dgam1			!der /rstar
           dy81m=dy81 *dgamm *dgam1			!der /m
           dy81tau=dy81*dgamtau*dgam1			!der /ln tau
           dy81lpt=dy81*(dgamlpt*dgam1+dlppdpt/dlpp)	!der /dln Ptot
           dy81lp =dy81*( dgamlp*dgam1+dlppdp /dlpp)	!der /dln Pgaz
                   
           dy82=cte6*gam1*delta*prn*grad
           dy82pt=dy82*(dgampt*dgam1+1.d0+dgradpt)	!der /ln Ptot
           dy82p= dy82*(deltap+dgamp*dgam1+dgradp)	!der /ln Pgaz
           dy82t= dy82*(deltat+dgamt*dgam1+dgradt)	!der /ln T
           dy82r=dy82*(dgamr*dgam1+dgradr)		!der /r
           dy82rs=dy82*(dgamrs*dgam1+dgradrs)		!der /rs
           dy82m=dy82*(dgamm*dgam1+dgradm)		!der /m	
           dy82tau=dy82*(dgamtau*dgam1+dgradtau)	!der /tau
	   dy82lpt=dy82*(dgamlpt*dgam1+dgradlpt)	!der/dln Ptot
	   dy82lp= dy82*(dgamlp *dgam1+dgradlp)		!der/dln Pgaz

           be(Ipgt)=dy81-dy82-pgn+prn
c          WRITE(*,2000)be(Ipgt),dy81,dy82,prn,pgn,dy81-dy82,prn-pgn
          
           ae(Ipgt,1,0)=dy81pt-dy82pt+prn	!der /ln Ptot
           ae(Ipgt,1,1)=dy81lpt-dy82lpt	!der /dln Ptot
           ae(Ipgt,2,0)=dy81t-dy82t	!der /ln T         
           ae(Ipgt,3,0)=dy81r-dy82r	!der /R/Rsol
           ae(Ipgt,4,0)=dy81rs-dy82rs	!der /Rstar/Rsol        
           ae(Ipgt,5,0)=dy81m-dy82m	!der /M/Msol
           ae(Ipgt,7,0)=dy81tau-dy82tau	!der /ln tau
           ae(Ipgt,Ipgt,0)=dy81p-dy82p-pgn	!der /ln Pgaz
           ae(Ipgt,Ipgt,1)=dy81lp-dy82lp 	!der /dln Pgaz 
          ENDIF		!radiatif
         ENDIF		!pturb
c	 PRINT*,cx ; WRITE(*,2000)ae
c	 PRINT*,'be' ; WRITE(*,2000)be ; PAUSE	 
        
c	 deriv=cx <= 2 .OR. cx == 23          
c	 deriv=cx <= 3
c	1 .OR. (cx >= n23_atm-1 .and. cx <= n23_atm+1)
c	2 .OR. (cx >= n_atm-1 .and. cx <= n_atm)
c	 deriv=cx >= n_atm-2
c	 deriv=.TRUE.
c	 deriv=.FALSE.
c	 deriv=cx >= n23_atm-3 .and. cx <= n23_atm
c	 deriv=.not.radiatif
c	 deriv=cx == 12
c	 deriv=radiatif

	 deriv=.FALSE.
	 IF(deriv)THEN	!tests de mise au point pour dérivées
	 
	  ALLOCATE(bes(ne_atm),dern(ne_atm))
	  PRINT*,'rad,radiatif,cx,n23_atm',rad,radiatif,cx,n23_atm
	  PRINT*,'prn,pgn,trn,tau,dlpp,teff,l_rac,xchim(1)'
	  WRITE(*,2000)prn,pgn,trn,tau,dlpp,teff,l_rac,xchim(1)
	  PRINT*,'grav,m,r,rstar'
	  WRITE(*,2000)grav,y(5,0),y(3,0),y(4,0)
	  PRINT*,'y/dy'
	  WRITE(*,2000)(y(i,0),i=1,ne_atm)
	  WRITE(*,2000)(y(i,1),i=1,ne_atm)
	  PRINT*,'grad,dgradpt,dgradp,dgradt,dgradr,dgradrs,dgradm,
	1 dgradtau'
	  WRITE(*,2000)grad,dgradpt,dgradp,dgradt,dgradr,dgradrs,
	1 dgradm,dgradtau
	  PRINT*,'dlpp,dgradlpp dgradlpt dgradlp'
	  WRITE(*,2000)dlpp,dgradlpp,dgradlpt,dgradlp	
c	  PRINT*,'prn,grav,dy1,dy2,dy3,dy5,kap,tau'
c	  WRITE(*,2000)prn,grav,dy1,dy2,dy3,dy5,kap,tau
c	  PRINT*,'case_rot,iz/delfi,w',case_rot,iz 
c	  WRITE(*,2000)delfi,w,cte1*y(5,0)/y(3,0)**2,cte4*w*y(3,0)
c	  WRITE(*,2000)(be(i),i=1,ne_atm)
c         PRINT*,'cte2',cte1,cte2,l_rac,lsol,pi,rsol,aradia,clight,ne_atm
c	  PRINT*,'dkapp,dkapt,dtbr23,drop,drot,dtsdtau',xcoll_atm(cx),cx
c	  WRITE(*,2000)dkapp,dkapt,dtbr23,drop,drot,dtsdtau*tau/tb
	  PRINT*,'dkapp,dkapt,drop,drot'
	  WRITE(*,2000)dkapp,dkapt,drop,drot
	  PRINT*,'gam,dgampt,dgamt,dgamr,dgamrs,dgamm,dgamtau,dgamp,
	1 dgamlpt,dgamlp'	
	  WRITE(*,2000)gam,dgampt,dgamt,dgamr,dgamrs,dgamm,dgamtau,dgamp,
	1 dgamlpt,dgamlp	
	  PRINT*,' '
	  kap0=kap ; ro0=ro ; grad0=grad ; gam0=gam ; i=1
	  ALLOCATE(yb(size(y,1),0:1)) ; yb=y
	  DO WHILE(i <= 2*ne_atm)
	   IF(i <= ne_atm)THEN
            stor0=yb(i,0) ; stor=stor0*unpdd	
            IF(stor == 0.d0)stor=dd ; dstor=stor-stor0 ; yb(i,0)=stor
	   ELSE
            stor0=yb(i-ne_atm,1) ; stor=stor0*unpdd
	    IF(stor == 0.d0)stor=dd
	    dstor=stor-stor0 ; yb(i-ne_atm,1)=stor
	   ENDIF

	   prn=EXP(yb(1,0))	!pression totale cgs
	   trn=EXP(yb(2,0))	!temperature cgs
	   tau=EXP(yb(7,0))	!profondeur optique
	   IF(pturb)THEN	!avec pression turbulente 8 inconnues	 
	    pgn=EXP(yb(Ipgt,0))	 !pression gaz+rad cgs
	   ELSE			!sans pression turbulente 7 inconnues
	    pgn=prn
	   ENDIF
	 	 
c	   Omega

	   SELECT CASE(Krot)
	   CASE(0,1,2)
	    w=wrot
	   CASE(3,4,5)
	    w=rota(1,dim_rot)	 
	   END SELECT
	   w=w**2 
	
c	   la fonction f(tau)	
	
	   teff=(cte2*l_rac/yb(4,0)**2)**0.25d0
	   grav=cte1*mstar/yb(4,0)**2
	
	   CALL tdetau(tau,teff,grav,tb,dtsdtau,dtsdteff,dtsdg,
	1    ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)
	
	  IF(pturb .and. der)THEN
	   dlpp=yb(Ipgt,1)/yb(1,1)	!dlnPgaz/dlnPtot		
	  ELSE
	   dlpp=1.d0
	  ENDIF
	
	  CALL thermo_atm(prn,pgn,trn,Lxchim,yb(5,0),l_rac,yb(3,0),dlpp,
	1  tau,df_tau,d2f_tau,yb(4,0),	
	2  ro,drop,drot,kap,dkapp,dkapt,gradad,dgradadp,dgradadt,
	3  grad,dgradpt,dgradp,dgradt,dgradr,dgradrs,dgradm,dgradtau,dgradlpp,
	4  gam,dgampt,dgamp,dgamt,dgamr,dgamrs,dgamm,dgamtau,dgamlpp,	
	5  hp,dhppt,dhpp,dhpt,dhpr,dhpm,delta,deltap,deltat,
	6  cp,gradrad,alfa,beta,gamma1,radiatif,.FALSE.)
	
c	   WRITE(*,2000)dgradpt,dgradp,dgradt,dgradr,dgradrs,dgradm,dgradtau,dgradlpp
c	   WRITE(*,2000)dgampt,dgamp,dgamt,dgamr,dgamrs,dgamm,dgamtau,dgamlpp
c	   PAUSE'après les grad'

c	   Définition de la fonction fi

           IF(xcoll_atm(cx) <= n23_atm)THEN
            delfi=(yb(6,0)-ltauf)*delfim  
           ELSE
            delfi=(yb(6,0)-ltaue)*delfip
           ENDIF
	 
	   dy1=(cte1*yb(5,0)/yb(3,0)**2-cte4*w*yb(3,0))*delfi/kap*
	1  EXP(yb(7,0)-yb(1,0))
	   bes(1)=yb(1,1)-dy1

c	   dy2=yb(1,1)*grad ; dy2=dy1*grad ; bes(2)=yb(2,1)-dy2 
	 
           dy3=cte5/kap/ro*tau*delfi ; bes(3)=yb(3,1)-dy3 ; bes(4)=yb(4,1)

           dy5=-cte3*yb(3,0)**2*tau*delfi/kap
           bes(5)=yb(5,1)-dy5
	   bes(6)=yb(6,1)
	   bes(7)=yb(7,1)-delfi
         
	   IF(pturb)THEN		!avec pression turbulente 8 inconnues
            IF(radiatif)THEN
             bes(Ipgt)=yb(1,1)-yb(Ipgt,1)
            ELSE         
             gam1=gam/(gam+1.d0) ; dy81=cte6*gam1*delta*prn*gradad*dlpp
             dy82=cte6*gam1*delta*prn*grad ; bes(Ipgt)=dy81-dy82-pgn+prn
            ENDIF
           ENDIF
  
	   DO j=1,ne_atm
	    dern(j)=(bes(j)-be(j))/dstor
	   ENDDO	!j
	   PRINT*	 
	   IF(i <= ne_atm)THEN
	    PRINT*,'dérivée / ',variable(i)
	    WRITE(*,2000)(dern(j),ae(j,i,0),j=1,ne_atm)
	    yb(i,0)=stor0	  	  
	   ELSE
	    PRINT*,'dérivée /dérivée de ',variable(i-ne_atm)
	    WRITE(*,2000)(dern(j),ae(j,i-ne_atm,1),j=1,ne_atm)
	    yb(i-ne_atm,1)=stor0	  	  
	   ENDIF	 
	   PRINT*,'dgrad,dro,dkap,dgam'
	   WRITE(*,2000)(grad-grad0)/dstor/grad,(ro-ro0)/dstor/ro,
	1  (kap-kap0)/dstor/kap,(gam-gam0)/dstor
	   PRINT*
	   i=i+1
	  ENDDO	!i
c	  PAUSE'limite 1'
	  deriv=.FALSE.
	  DEALLOCATE(bes,dern,yb)
	 ENDIF

c-------------	limites	---------------------------

	CASE(2)
	 SELECT CASE(li)	 
	 CASE(1)	!tau=tau_max, q=1 R/Rsol=y(3)=Rb/Rsol=r_rac
          be(1)=y(3,0)-r_rac ; ae(1,3,0)=1.d0             !dérivée /R/Rsol
c         PRINT*,'case 1' ; WRITE(*,2000)y(3,0),r_rac,be(1)

	 CASE(2) 	!tau=tau*, q=n23_atm, ln T=y(2)=ln T(tau,Teff,..)
	  tau=EXP(y(7,0))		!profondeur optique
	  teff=(cte2*l_rac/y(4,0)**2)**.25
	  grav=cte1*mstar/y(4,0)**2
	  CALL tdetau(tau,teff,grav,tb,dtsdtau,dtsdteff,dtsdg,
	1   ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)
  	  be(1)=y(2,0)-log(tb)
	  dteffre=-teff/y(4,0)/2.d0			!dérivée teff / r23
	  dgravre=-2.d0*grav/y(4,0)			!dérivée grav / r23
	  dtbr23=dtsdteff*dteffre+dtsdg*dgravre		!d tb / r23
	  ae(1,2,0)=1.d0				!dérivée /ln T
	  ae(1,4,0)=-dtbr23/tb				!dérivée /r23/Rsol
	  ae(1,7,0)=-dtsdtau*tau/tb			!dérivée /ln tau
c         PRINT*,'case 2' ; WRITE(*,2000)y(2,0),tb,be(1)

	  deriv=.FALSE.
c	  deriv=.TRUE.
	  IF(deriv)THEN
	   ALLOCATE(yb(size(y,1),0:1)) ; yb=y
	   deriv=.FALSE. ; PRINT*,'limite',li
           DO i=1,ne_atm
            stor0=yb(i,0) ; stor=stor0*unpdd	
            IF(stor == 0.d0)stor=dd ; dstor=stor-stor0 ; yb(i,0)=stor
	    tau=EXP(yb(7,0))		!profondeur optique
	    teff=(cte2*l_rac/yb(4,0)**2)**0.25d0
	    grav=cte1*mstar/yb(4,0)**2
	    CALL tdetau(tau,teff,grav,tb,dtsdtau,dtsdteff,dtsdg,
	1     ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)  
	    vt=yb(2,0)-log(tb) ; yb(i,0)=stor0
	    PRINT*,'dérivée',i ; WRITE(*,2000)(vt-be(1))/dstor,ae(1,i,0)
	   ENDDO
	   PAUSE'limite 2' ; DEALLOCATE(yb)
	  ENDIF

	 CASE(3) 	!tau=tau*, q=n23_atm, T=Teff
	  teff=(cte2*l_rac/y(4,0)**2)**0.25d0
	  be(1)=y(2,0)-log(teff)
	  dteffre=-teff/y(4,0)/2.d0	!dérivée teff / r23
	  ae(1,2,0)=1.d0 		!dérivée / ln T
	  ae(1,4,0)=-dteffre/teff 	!dérivée / r23/Rsol
c	  PRINT*,'case 3' ; WRITE(*,2000)y(2,0), teff, be(1)	  
	  

	  deriv=.FALSE.
c	  deriv=.TRUE.
	  IF(deriv)THEN
	   ALLOCATE(yb(size(y,1),0:1)) ; yb=y
	   deriv=.FALSE. ; PRINT*,'limite',li
           DO i=1,ne_atm
            stor0=yb(i,0) ; stor=stor0*unpdd	
            IF(stor == 0.d0)stor=dd ; dstor=stor-stor0 ; yb(i,0)=stor
	    teff=(cte2*l_rac/yb(4,0)**2)**0.25d0 ; vt=yb(2,0)-log(teff)
	    yb(i,0)=stor0 ; PRINT*,'dérivée',i
	    WRITE(*,2000)(vt-be(1))/dstor,ae(1,i,0)
	   ENDDO
	   PAUSE'limite 3'
	   DEALLOCATE(yb)
	  ENDIF

	 CASE(4) 	!tau=tau*
	  be(1)=y(6,0)-y(7,0)
	  ae(1,6,0)=1.d0		!dérivée /log tau23
	  ae(1,7,0)=-1.d0		!dérivée /log tau
c	  PRINT*,'case 4' ; WRITE(*,2000)be(1),y(6,0), y(7,0) 	  

	 CASE(5)	!tau=tau*, q=n23_atm R/Rsol=y(3)=y(4)=R23/Rsol
	  be(1)=y(3,0)-y(4,0)
          ae(1,3,0)=1.d0           !dérivée /R/Rsol
          ae(1,4,0)=-1.d0          !dérivée /R23/Rsol
c	  PRINT*,'case 5' ; WRITE(*,2000)be(1),y(3,0),y(4,0)	  

	 CASE(6)	!tau=tau*, q=n23_atm M=y(5)=M*/Msol
	  be(1)=y(5,0)-mstar
          ae(1,5,0)=1.d0           !dérivée /M/Mstar
c	  PRINT*,'case 6' ; WRITE(*,2000)be(1),y(5,0),mstar	  

	 CASE(7)	!en tau=tau_min, q=N, P=g/kappa tau_min
	  Lxchim=xchim
	  IF(lim_ro)THEN		!limite sur la densité
	   teff=(cte2*l_rac/y(4,0)**2)**0.25d0
	   dteffre=-teff/y(4,0)/2.d0			!dérivée teff / r23
	   grav=cte1*mstar/y(4,0)**2
	   dgravre=-2.d0*grav/y(4,0)			!dérivée grav / r23
	   CALL tdetau(tau,teff,grav,tb,dtsdtau,dtsdteff,dtsdg,
	1    ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)
	   IF(pturb)THEN		!avec pression turbulente 8 inconnues	
            pgn=EXP(y(Ipgt,0))
	   ELSE			!sans pression turbulente 7 inconnues
            pgn=EXP(y(1,0))	 
	   ENDIF
           trn=EXP(y(2,0))
           CALL etat(pgn,trn,Lxchim,.FALSE.,
	1    ro,drop,drot,drox,u,dup,dut,dux,
	2    delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3    gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
         
	   be(1)=ro-ro_ext			!ro(n)=ro_ext
	   IF(pturb)THEN	!avec pression turbulente 8 inconnues
	    ae(1,Ipgt,0)=drop*pgn	!dérivée /log Pgas
	   ELSE			!sans pression turbulente 7 inconnues
	    ae(1,1,0)=drop*pgn	!dérivée /log Ptot	 
	   ENDIF
           ae(1,2,0)=drot*trn	!dérivée /log t
           ae(1,4,0)=-dro_teff*dteffre-dro_grav*dgravre	!dérivée / r23
c	   PRINT*,'case 7' ; WRITE(*,2000)be(1),ro,ro_ext	

	   deriv=.FALSE.
c	   deriv=.TRUE.
	   IF(deriv)THEN
	    ALLOCATE(yb(size(y,1),0:1)) ; yb=y
	    deriv=.FALSE. ; PRINT*,'limite',li
            DO i=1,ne_atm
             stor0=yb(i,0) ; stor=stor0*unpdd ; IF(stor == 0.d0)stor=dd
	     dstor=stor-stor0 ; yb(i,0)=stor
	     teff=(cte2*l_rac/yb(4,0)**2)**0.25d0 ; grav=cte1*mstar/yb(4,0)**2
	     CALL tdetau(tau,teff,grav,tb,dtsdtau,dtsdteff,dtsdg,
	1      ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)
	     IF(pturb)THEN	!avec pression turbulente 8 inconnues
              pgn=EXP(yb(Ipgt,0))
	     ELSE		!sans pression turbulente 7 inconnues
              pgn=EXP(yb(1,0))	 
	     ENDIF
             trn=EXP(yb(2,0))
             CALL etat(pgn,trn,Lxchim,.FALSE.,
	1    ro,drop,drot,drox,u,dup,dut,dux,
	2    delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3    gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	     vt=ro-ro_ext                         !ro(n)=ro_ext

             yb(i,0)=stor0
	     PRINT*,'dérivée',i ; WRITE(*,2000)(vt-be(1))/dstor,ae(1,i,0)
	    ENDDO
c	    PAUSE'limite 7-1'
	    DEALLOCATE(yb)
	   ENDIF
	  
	  ELSE			!limite sur la gravite
	   grav=cte1*mstar/y(4,0)**2
	   dgravre=-2.d0*grav/y(4,0)		!dérivée grav / r23
           prn=EXP(y(1,0)) ; trn=EXP(y(2,0))
	   IF(pturb)THEN	!avec pression turbulente 8 inconnues
            pgn=EXP(y(Ipgt,0))
	   ELSE			!sans pression turbulente 7 inconnues
	    pgn=prn
	   ENDIF
	   CALL etat(pgn,trn,Lxchim,.FALSE.,
	1  ro,drop,drot,drox,u,dup,dut,dux,
	2  delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3  gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	   CALL opa(Lxchim,trn,ro,kap,dkapt,dkapro,dkapx)
	   dkapp=drop*dkapro*pgn/kap	!d kap/d ln p à kap près
	   dkapt=(drot*dkapro+dkapt)*trn/kap	!d kap/d ln t à kap près
	 
	   dy1=grav/kap*tau_min
	   be(1)=prn-dy1
	   ae(1,1,0)=prn				!dérivée / ln Ptot
	   ae(1,2,0)= dy1*dkapt				!dérivée / ln T
	   ae(1,4,0)=-dy1/grav*dgravre			!dérivée / log R*
	   IF(pturb)THEN	!avec pression turbulente 8 inconnues
	    ae(1,Ipgt,0)=dy1*dkapp			!dérivée / ln Pgaz
	   ELSE			!sans pression turbulente 7 inconnues
	    ae(1,1,0)=ae(1,1,0)+dy1*dkapp		!dérivée / ln Ptot
	   ENDIF
	  	  
c	   PRINT*,'prn,dy1,ro' ; WRITE(*,2000)prn,dy1,ro	 
	 
	   deriv=.FALSE.
c	   deriv=.TRUE.
	   IF(deriv)THEN
	    ALLOCATE(yb(size(y,1),0:1)) ; yb=y ; deriv=.FALSE.
	    PRINT*,'limite',li
            DO i=1,ne_atm
             stor0=yb(i,0) ; stor=stor0*unpdd ; IF(stor == 0.d0)stor=dd
	     dstor=stor-stor0 ; yb(i,0)=stor ; grav=cte1*mstar/yb(4,0)**2
             prn=EXP(yb(1,0)) ; trn=EXP(yb(2,0))
	     IF(pturb)THEN	!avec pression turbulente 8 inconnues
              pgn=EXP(yb(Ipgt,0))
	     ELSE			!sans pression turbulente
	      pgn=prn
	     ENDIF
             CALL etat(pgn,trn,Lxchim,.FALSE.,
	1    ro,drop,drot,drox,u,dup,dut,dux,
	2    delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
	3    gradad,dgradadp,dgradadt,dgradadx,alfa,beta,gamma1)
	     CALL opa(Lxchim,trn,ro,kap,dkapt,dkapro,dkapx)
	     dy1=grav/kap*tau_min ; vt=prn-dy1 ; yb(i,0)=stor0
	     PRINT*,'dérivée',i ; WRITE(*,2000)(vt-be(1))/dstor,ae(1,i,0)
	    ENDDO
	    PAUSE'limite 7-2'
	    DEALLOCATE(yb)
	   ENDIF
	  ENDIF		!lim_ro	  

	 CASE(8)	!en tau=tau_min, q=N, ln Pgas=ln Ptot
	  be(1)=y(1,0)-y(Ipgt,0) ; ae(1,1,0)=1.d0 ; ae(1,Ipgt,0)=-1.d0
	  
	 CASE DEFAULT
	  PRINT*,'dans eqatm erreur sur li=1, 2, ..., Ipgt, li=',li
	  PRINT*,'arrêt' ; STOP
	 END SELECT		!li
	 
	CASE DEFAULT
	 PRINT*,'dans eqatm erreur sur fait=1, 2 ou 3, fait=',fait
	 PRINT*,'arrêt' ; STOP
	
	END SELECT		!fait

	RETURN
	
	END SUBROUTINE eq_atm
