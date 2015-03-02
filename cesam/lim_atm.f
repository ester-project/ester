
c****************************************************************

	SUBROUTINE lim_atm(list,l_rac,r_rac,xchim,
	1 pt_rac,dptsdl,dptsdr,t_rac,dtsdl,dtsdr,
	2 m_rac,dmsdl,dmsdr,p_rac,dpsdl,dpsdr,t_eff)

c	routine public du module mod_atm
c	restitution de l'atmosphère

c	modifs :
c	04 07 01 : mise en place de rotation uniforme avec conservation
c	du moment angulaire
c	16 08 01 : F95

c	pour la restitution de l'atmosphère on utilise des splines
c	du même ordre que pour le modèle quasi statique m=m_qs et
c	pour les équa. dfiff. r_atm=1

c	Auteur: P. Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

cYLD     YLD : Modifiée pour faire varier le tau_max dans le cas d'une loi
c	T(tau) non radiative : septembre 2003

c entrées :
c	list=.TRUE. : calcul réduit pour une liste
c	r_rac : rayon au raccord
c	l_rac : luminosité
c	xchim composition chimique par gramme
c	mstar: masse avec perte de masse
c	wrot : vitesse angulaire solide

c sorties :
c	pt_rac : Ptot au raccord avec l'enveloppe
c	t_rac : température au raccord avec l'enveloppe
c	m_rac : masse au raccord avec l'enveloppe
c	p_rac : Pgaz au raccord avec l'enveloppe
c	d*r, d*l : dérivées /r /l
c	t_eff : température effective
c	rstar : rayon de l'étoile
c	pt_atm, t_atm, r_atm, m_atm, p_atm, dlpp_atm : Ptot, temperature,
c	rayon, masse, Pgaz, dln Pgaz/dln Ptot 
c	tau : profondeur optique	
c	avec pression turbulente 8 inconnues 
c	sans pression turbulente 7 inconnues, Ptot=Pgaz

c-----------------------------------------------------------------------

	USE mod_donnees, ONLY : aradia, clight, g, granr, langue,
	1 lsol, lim_ro, msol, m_qs, nom_fich2, n_atm, pi, pturb,
	2 rep_atm, rsol, r_qs, tau_max
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn, newspl, noedif, no_croiss, pause
	USE mod_variables, ONLY : mstar, rstar
      
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim
	REAL (kind=dp), INTENT(in) :: l_rac, r_rac
	LOGICAL, INTENT(in) :: list
	REAL (kind=dp), INTENT(out) :: dmsdl, dmsdr, dpsdl, dpsdr, dptsdl,
	1 dptsdr, dtsdl, dtsdr, m_rac, p_rac, pt_rac, t_eff, t_rac
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: b, bd
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: dfxdx, fx
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: ltau, mltau,
	1 mltaut, to, tot
	REAL (kind=dp), PARAMETER :: unpdx=1.00001d0
	REAL (kind=dp), SAVE :: cte1, cte2
	REAL (kind=dp) :: bid, delfi, df_tau, dl, dr, dro_grav, dro_teff,
	1 dtsdg, dtsdtau, dtsdteff, d2f_tau, f_tau, grav, ld, ltau23, rd,
	2 ro_ext, tau23, teff

	INTEGER, SAVE :: dim_atm, dim_atm_rep, knot_atm, ms_atm, ord_atm
	INTEGER :: i, j, knot0, ll, npt, n_atm_rep

	LOGICAL, SAVE ::  ini=.TRUE., ok
	
	CHARACTER (len=80) :: fich_atm

c---------------------------------------------------------------------

2000	FORMAT(8es10.3)

c	PRINT*,'lim_ext,tau_max,n_atm,rep_atm',tau_max,n_atm,rep_atm
c	WRITE(*,2000)l_rac,r_rac ; PAUSE'l_rac'

	IF(ini)THEN
	 ini=.FALSE.
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1006) ; WRITE(2,1006)
1006	  FORMAT(/,'------------ATMOSPHERE (initialisation)-----------',/)	 
	 CASE DEFAULT	 
	  WRITE(*,6) ; WRITE(2,6)
6	  FORMAT(/,'------------ATMOSPHERE (initialisation)-----------',/)
	 END SELECT
	 
c les constantes ``splines'' pour l'intégration
c	 par souci de cohérence ordre des splines est le même que celui de
c	 l'équilibre quasi-statique

c	 ne_atm: nombre de variables
c	 ms_atm: ordre des splines
c	 r_qs: ordre des équa. diff.
c	 ord_atm: ordre de l'interpolation
c	 ncoll_atm: nombre de points de collocation
c	 dim_atm: dimension de l'espace des splines
c	 n23_atm: indice pour lequel T=Teff

	 cte1=g*msol/rsol**2 ; cte2=lsol/pi/rsol**2/aradia/clight
	 teff=(cte2*l_rac/r_rac**2)**0.25d0 ; grav=cte1*mstar/r_rac**2
	 
	 IF(pturb)THEN		!avec pression turbulente 8 inconnues
	  ne_atm=8 ; Ipgt=8	!nombre d'inconnues
	 ELSE			!sans pression turbulente 7 inconnues
	  ne_atm=7 ; Ipgt=-100
	 ENDIF
	 ms_atm=m_qs ; ord_atm=ms_atm+r_qs
	 dim_atm=(n_atm-1)*ms_atm+r_qs ; n23_atm=NINT(3.d0/4.d0*n_atm)
	 
c	 reprise d'une atmosphère, test de cohérence
c	 si incohérence on utilise un modèle interne

c	 PRINT*,rep_atm ; PAUSE'rep_atm1'

	 IF(rep_atm)THEN 	 
	  fich_atm=TRIM(nom_fich2)//'_B.atm' 
	  INQUIRE(file=TRIM(fich_atm),exist=ok) 
c	  PRINT*,ok ; PAUSE'ici'
          IF(ok)THEN
	   OPEN(unit=31,form='unformatted',status='unknown',
	1    file=TRIM(nom_fich2)//'_B.atm')
	   READ(31)dim_atm_rep,n_atm_rep,knot_atm
	   IF(dim_atm_rep /= dim_atm .OR. n_atm_rep /= n_atm)THEN
	    SELECT CASE(langue)
	    CASE('english')
	     WRITE(*,1009)dim_atm_rep,dim_atm,n_atm_rep,n_atm
1009	     FORMAT('unconsistent model of atmosphere ',/,'dim_atm_rep=',
	1    i3,', dim_atm=',i3,'; n_atm_rep=',i3,', n_atm=',i3,:,
	2    'use of an internal model')	    
	    CASE DEFAULT	   
	     WRITE(*,9)dim_atm_rep,dim_atm,n_atm_rep,n_atm
9	     FORMAT('atmosphère reprise incohérente',/,'dim_atm_rep=',i3,
	1    ', dim_atm=',i3,'; n23_atm_rep=',i3,', n23_atm=',i3,:,
	2    'on utilise un modèle interne')
	    END SELECT	
	    rep_atm=.FALSE. ; CLOSE(UNIT=31)
	   ELSE	   
	    ALLOCATE(bp_atm(ne_atm,dim_atm),bp_atm_t(ne_atm,dim_atm),
	1   x_atm(n_atm),x_atm_t(n_atm),xt_atm(knot_atm),
	2   xt_atm_t(knot_atm))
	    REWIND(UNIT=31) 	  
	    READ(31)dim_atm,n_atm,knot_atm,bp_atm,x_atm,xt_atm
	    CLOSE(UNIT=31)	    
	    SELECT CASE(langue)
	    CASE('english')
             WRITE(*,1007)TRIM(fich_atm) ; WRITE(2,1007)TRIM(fich_atm)
1007	     FORMAT('the atmosphere is initialized with the model : ',a)	    
	    CASE DEFAULT	    	   
             WRITE(*,7)TRIM(fich_atm) ; WRITE(2,7)TRIM(fich_atm)
7	     FORMAT('initialisation de l''atmosphère avec le modèle : ',a)
	    END SELECT
	   ENDIF
	  ELSE
	   rep_atm=.FALSE.
	  ENDIF
	 ENDIF
	  
c reprise de l'atmosphère ou initialisation avec un modèle interne

c 	 PRINT*,rep_atm ; PAUSE'rep_atm2'
	 IF(rep_atm)THEN

c	  appel fictif à T(tau) pour init. des delfim, delfip, etc.. 

	  CALL tdetau(1.d0,6.d3,5.d4,bid,dtsdtau,dtsdteff,dtsdg,
	1  ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)

c	  avec une loi T(tau) non purement radiative tau_max est imposé

	  IF(.NOT.rad)THEN
cYLD	   tau_max=MAX(tau_max,20.d0)
	   WRITE(2,1)tau_max,tau_min,n_atm,n23_atm
	   WRITE(*,1)tau_max,tau_min,n_atm,n23_atm    
	  ELSE
	   WRITE(2,2)tau_max,tau_min,n_atm,n23_atm
	   WRITE(*,2)tau_max,tau_min,n_atm,n23_atm    
	  ENDIF
	  IF(lim_ro)THEN
	   WRITE(2,3) ; WRITE(*,3)
	  ELSE
	   WRITE(2,4) ; WRITE(*,4)
	  ENDIF

c recherche de tau23 où T=Teff

	  teff=(cte2*l_rac/r_rac**2)**0.25d0 ; grav=cte1*mstar/r_rac**2
	  CALL taueff(teff,grav,tau23)

c	  changement de variable tau --> psi linéaire par morceaux 

	  ltauf=log(tau_max) ; ltaue=log(tau_min) ; ltau23=log(tau23)
	  delfim=1.d0/(n23_atm-1) ; delfip=1.d0/(n23_atm-n_atm)

	  ALLOCATE(fx(ne_atm),dfxdx(ne_atm))

	 ELSE
	  SELECT CASE(langue)
	  CASE('english')
           WRITE(*,1008)teff ; WRITE(2,1008)teff
1008	   FORMAT('init. atm. with internal model for Teff=',es10.3)	  
	  CASE DEFAULT	 
           WRITE(*,8)teff ; WRITE(2,8)teff
8	   FORMAT('init. atm. avec modèle interne pour Teff=',es10.3)
	  END SELECT

c	  valeurs initiales des pressions et profondeurs optiques

	  npt=60 ; ALLOCATE(bd(2,npt),to(npt))

	  IF(teff < 7000.d0 ) THEN
	   bd(1,:)=(/
	1  3.501D+05,3.426D+05,3.353D+05,3.284D+05,3.218D+05,3.154D+05,
	1  3.093D+05,3.035D+05,2.979D+05,2.926D+05,2.874D+05,2.824D+05,
	2  2.776D+05,2.729D+05,2.683D+05,2.638D+05,2.595D+05,2.553D+05,
	3  2.511D+05,2.471D+05,2.432D+05,2.393D+05,2.355D+05,2.318D+05,
	4  2.280D+05,2.242D+05,2.204D+05,2.165D+05,2.125D+05,2.084D+05,
	5  2.042D+05,1.998D+05,1.953D+05,1.906D+05,1.855D+05,1.802D+05,
	6  1.746D+05,1.687D+05,1.624D+05,1.560D+05,1.493D+05,1.424D+05,
	7  1.354D+05,1.285D+05,1.218D+05,9.094D+04,6.713D+04,4.924D+04,
	8  3.597D+04,2.616D+04,1.895D+04,1.370D+04,9.919D+03,7.187D+03,
	9  5.195D+03,3.739D+03,2.671D+03,1.877D+03,1.269D+03,7.785D+02/)      
	   to=(/
	1  1.000D+02,8.917D+01,7.951D+01,7.090D+01,6.323D+01,5.638D+01,
	1  5.027D+01,4.483D+01,3.998D+01,3.565D+01,3.179D+01,2.834D+01,
	2  2.527D+01,2.254D+01,2.010D+01,1.792D+01,1.598D+01,1.425D+01,
	3  1.271D+01,1.133D+01,1.010D+01,9.009D+00,8.034D+00,7.164D+00,
	4  6.388D+00,5.696D+00,5.079D+00,4.529D+00,4.039D+00,3.602D+00,
	5  3.212D+00,2.864D+00,2.554D+00,2.277D+00,2.031D+00,1.811D+00,
	6  1.615D+00,1.440D+00,1.284D+00,1.145D+00,1.021D+00,9.103D-01,
	7  8.117D-01,7.238D-01,6.454D-01,3.596D-01,2.004D-01,1.117D-01,
	8  6.221D-02,3.467D-02,1.932D-02,1.076D-02,5.997D-03,3.341D-03,
	9  1.862D-03,1.037D-03,5.781D-04,3.221D-04,1.795D-04,1.000D-04/)
          ELSE     
	   bd(1,:)=(/
	1  7.512D+04,6.740D+04,6.037D+04,5.397D+04,4.817D+04,4.295D+04,
	1  3.828D+04,3.411D+04,3.039D+04,2.706D+04,2.408D+04,2.141D+04,
	2  1.903D+04,1.691D+04,1.503D+04,1.338D+04,1.193D+04,1.067D+04,
	3  9.577D+03,8.638D+03,7.836D+03,7.158D+03,6.593D+03,6.126D+03,
	4  5.743D+03,5.430D+03,5.177D+03,4.971D+03,4.804D+03,4.668D+03,
	5  4.555D+03,4.462D+03,4.382D+03,4.313D+03,4.252D+03,4.196D+03,
	6  4.144D+03,4.093D+03,4.044D+03,3.995D+03,3.944D+03,3.893D+03,
	7  3.840D+03,3.785D+03,3.728D+03,3.401D+03,3.039D+03,2.692D+03,
	8  2.407D+03,2.198D+03,2.059D+03,1.972D+03,1.920D+03,1.889D+03,
	9  1.872D+03,1.862D+03,1.856D+03,1.853D+03,1.852D+03,1.851D+03/)
	   to=(/
	1  1.000D+02,8.917D+01,7.951D+01,7.090D+01,6.323D+01,5.638D+01,
	1  5.027D+01,4.483D+01,3.998D+01,3.565D+01,3.179D+01,2.834D+01,
	2  2.527D+01,2.254D+01,2.010D+01,1.792D+01,1.598D+01,1.425D+01,
	3  1.271D+01,1.133D+01,1.010D+01,9.009D+00,8.034D+00,7.164D+00,
	4  6.388D+00,5.696D+00,5.079D+00,4.529D+00,4.039D+00,3.602D+00,
 	5  3.212D+00,2.864D+00,2.554D+00,2.277D+00,2.031D+00,1.811D+00,
 	6  1.615D+00,1.440D+00,1.284D+00,1.145D+00,1.021D+00,9.103D-01,
	7  8.117D-01,7.238D-01,6.454D-01,3.596D-01,2.004D-01,1.117D-01,
	8  6.221D-02,3.467D-02,1.932D-02,1.076D-02,5.997D-03,3.341D-03,
	9  1.862D-03,1.037D-03,5.781D-04,3.221D-04,1.795D-04,1.000D-04/)
	  ENDIF
c	  bd(2,:) températures aux points tau=to

	  DO i=1,npt
	   CALL tdetau(to(i),teff,grav,bd(2,i),dtsdtau,dtsdteff,dtsdg,
	1  ro_ext,dro_grav,dro_teff,f_tau,df_tau,d2f_tau)
	  ENDDO
c	  PRINT*,'teff, grav, températures initiales'
c	  WRITE(*,2000)teff,grav ; WRITE(*,2000)bd(2,:) ; PAUSE'teff'

c	  retournement et en ln

	  bd(:,1:npt:+1)=bd(:,npt:1:-1) ; to(1:npt:+1)=to(npt:1:-1)
	  bd=log(bd) ; to=log(to)

c	  avec une loi T(tau) non purement radiative tau_max est imposé

	  IF(.not.rad)THEN
cYLD	   tau_max=MAX(tau_max,20.d0)
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(2,1001)tau_max,tau_min,n_atm,n23_atm
	    WRITE(*,1001)tau_max,tau_min,n_atm,n23_atm    
1001	    FORMAT('atmosphere restaured with connection of the gradient',/,
	1   'tau at the bottom of the atmosphere tau_max=',es10.3,/,
	2   'imposed with a law T(tau) no purely radiative',/,
	3   'radius and mass fixed at tau* as T(tau*,Teff)=Teff',/,
	4   'tau exterior',es10.3,', number of shells : ',i3,
	5   ', R* at the shell ',i2)	   
	   CASE DEFAULT	   
	    WRITE(2,1)tau_max,tau_min,n_atm,n23_atm
	    WRITE(*,1)tau_max,tau_min,n_atm,n23_atm    
1	    FORMAT('atmosphère restituée avec raccord du gradient',/,
	1   'tau au fond de l''atmosphère tau_max=',es10.3,/,
	2   'imposé avec une loi T(tau) non purement radiative',/,
	3   'rayon et masse pris a tau* tel que T(tau*,Teff)=Teff',/,
	4   'tau extérieur',es10.3,', nombre de couches: ',i3,
	5   ', R* à la couche ',i2)
	   END SELECT		
	  ELSE
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(2,1002)tau_max,tau_min,n_atm,n23_atm
	    WRITE(*,1002)tau_max,tau_min,n_atm,n23_atm    
1002	    FORMAT('atmosphere restaured with connection of the gradient',/,
	1   'tau at the bottom of the atmosphere tau_max=',es10.3,/,
	2   'radius and mass fixed at tau* as T(tau*,Teff)=Teff',/,
	3   'tau exterior',es10.3,', number of shells : ',i3,
	4   ', R* at the shell # ',i2)	   
	   CASE DEFAULT	   
	    WRITE(2,2)tau_max,tau_min,n_atm,n23_atm
	    WRITE(*,2)tau_max,tau_min,n_atm,n23_atm    
2	    FORMAT('atmosphère restituée avec raccord du gradient',/,
	1   'tau au fond de l''atmosphère tau_max=',es10.3,/,
	2   'rayon et masse pris a tau* tel que T(tau*,Teff)=Teff',/,
	3   'tau exterieur',es10.3,', nombre de couches: ',i3,
	4   ', R* à la couche # ',i2)
	   END SELECT
	  ENDIF
	  IF(lim_ro)THEN
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(2,1003) ; WRITE(*,1003)
1003	    FORMAT('external limit on the density')	   
	   CASE DEFAULT	   
	    WRITE(2,3) ; WRITE(*,3)
3	    FORMAT('limite externe portant sur la densité')
	   END SELECT
	  ELSE
	   SELECT CASE(langue)
	   CASE('english')
	    WRITE(2,1004) ; WRITE(*,1004)
1004	    FORMAT('external limit on the gravity')	   
	   CASE DEFAULT	   
	    WRITE(2,4) ; WRITE(*,4)
4	    FORMAT('limite externe portant sur la gravité')
	   END SELECT
	  ENDIF

c	  recherche de tau23 ou T=Teff

	  teff=(cte2*l_rac/r_rac**2)**0.25d0 ; grav=cte1*mstar/r_rac**2
	  CALL taueff(teff,grav,tau23)
	 
c	  WRITE(*,2000)teff,grav,tau23,l_rac,r_rac
c	  PRINT*,'initialisation lim_atm'

c	  changement de variable tau --> psi linéaire par morceaux 

	  ltauf=log(tau_max) ; ltaue=log(tau_min) ; ltau23=log(tau23)
	  delfim=1.d0/(n23_atm-1) ; delfip=1.d0/(n23_atm-n_atm)

	  ALLOCATE(ltau(n_atm))
	  DO i=1,n_atm
	   IF(i <= n23_atm)THEN
	    delfi=(ltau23-ltauf)*delfim ; ltau(i)=ltauf+delfi*(i-1)
	   ELSE
	    delfi=(ltau23-ltaue)*delfip ; ltau(i)=ltaue+delfi*(i-n_atm)
	   ENDIF
	  ENDDO

c	  PRINT*,'n_atm,n23_atm,m_qs',n_atm,n23_atm,m_qs
c	  PRINT*,'tau_max,tau23,tau_min' ; WRITE(*,2000)tau_max,tau23,tau_min
c	  PRINT*,'ln tau' ; WRITE(*,2000)ltau ; PRINT*,'tau'
c	  WRITE(*,2000)EXP(ltau) ; PRINT*,'les tau lim_atm'

c	  projection des pressions et temperatures sur une base de
c	  B-splines d'ordre 2 en ln to

	  ALLOCATE(tot(npt+2),fx(ne_atm),dfxdx(ne_atm))
	  CALL bsp1dn(2,bd,to,tot,npt,2,knot0,.FALSE.,to(1),ll,fx,dfxdx)
 	  
c	  PRINT*,'projection de bd'

c	  formation du vecteur solution initiale sur les ln to

	  ALLOCATE(b(ne_atm,n_atm))
	  DO i=1,n_atm
	   CALL bsp1dn(2,bd,to,tot,npt,2,knot0,.TRUE.,ltau(i),ll,fx,dfxdx)
	   IF(no_croiss)PRINT*,'Pb. at 1 in lim_atm'
	   b(1,i)=fx(1) ; b(2,i)=fx(2) ; b(3,i)=r_rac; b(4,i)=r_rac
	   b(5,i)=mstar ; b(6,i)=ltau23 ; b(7,i)=ltau(i)
	   IF(pturb)b(Ipgt,i)=fx(1)
	  ENDDO
c	  PRINT*,'les B'

c	  projection du vecteur solution initiale sur une base de
c	  B-splines d'ordre 2 en ltau
c	  ltau étant décroissant on utilise mltau=-ltau

	  ALLOCATE(mltau(n_atm),mltaut(n_atm+2)) ; mltau=-ltau

	  CALL bsp1dn(ne_atm,b,mltau,mltaut,n_atm,2,knot0,.FALSE.,mltau(1),
	1  ll,fx,dfxdx)
	
	  DO i=1,ne_atm
	   CALL bsp1dn(ne_atm,b,mltau,mltaut,n_atm,2,knot0,.FALSE.,mltau(i),
	1    ll,fx,dfxdx)
c	     WRITE(*,2000)mltau(i),fx
	  ENDDO
c	  PRINT*,'b sur ltau'

c	  abscisses d'intégration
	 
	  ALLOCATE(x_atm(n_atm),x_atm_t(n_atm),xt_atm(dim_atm+ord_atm),
	1   xt_atm_t(dim_atm+ord_atm))
	  x_atm=(/ (i, i=1,n_atm) /)

c	  projection du vecteur solution initiale sur une base de
c	  B-splines d'ordre 2 en x_atm (mltaut:VT)

	  CALL bsp1dn(ne_atm,b,x_atm,mltaut,n_atm,2,knot0,.FALSE.,x_atm(1),
	1   ll,fx,dfxdx)
c	  PRINT*,'b sur x_atm'

c	  vecteur nodal en ltau pour l'intégration

	  CALL noedif(x_atm,n_atm,ms_atm,r_qs,xt_atm,knot_atm)
          IF(no_croiss)THEN
           PRINT*,'Arrêt 5 dans lim_atm' ; STOP
	  ENDIF

c	  on place la solution initiale dans la base de x_atm, b-->bp_atm

	  ALLOCATE(bp_atm(ne_atm,dim_atm),bp_atm_t(ne_atm,dim_atm))
	  CALL newspl(ne_atm,mltau,mltaut,knot0,2,xt_atm,knot_atm,ord_atm,
	1  b,bp_atm)
c	  DO i=1,ne_atm
c	   CALL bsp1dn(ne_atm,bp_atm,x_atm,xt_atm,n_atm,ord_atm,knot_atm,
c	1  .TRUE.,x_atm(i),ll,fx,dfxdx)
c	   WRITE(*,2000)x_atm(i),fx
c	  ENDDO
c	  PAUSE'bp_atm'

c	  suppression des tableaux de travail

	  DEALLOCATE(bd,b,to,tot,ltau,mltau,mltaut)
	 
	 ENDIF

c	 allocations à n_atm des tableaux pt_atm etc..

	 ALLOCATE(pt_atm(n_atm),t_atm(n_atm),m_atm(n_atm),tau(n_atm),
	1 p_atm(n_atm),r_atm(n_atm),dlpp_atm(n_atm))

	ENDIF

c	intégration

	IF(.NOT.ini)THEN
	 PRINT*
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1011)l_rac,r_rac,teff   	   
1011	  FORMAT('-----Restauration of the atmosphère (begin)-----------'
	1 ,/,'Lrac=',es10.3,', Rrac=',es10.3,', Teff='es10.3)	 
	 CASE DEFAULT
	  WRITE(*,11)l_rac,r_rac,teff   	   
11	  FORMAT('-----Restitution de l''atmosphère (début)-----------',
	1 /,'Lrac=',es10.3,', Rrac=',es10.3,', Teff='es10.3)
	 END SELECT
	ENDIF
	CALL coll_atm(r_rac,l_rac,xchim,ord_atm,knot_atm,dim_atm,ms_atm)
	rstar=bp_atm(4,1)           !bp(4,i) est R_star, tous i
	t_eff=(cte2*l_rac/rstar**2)**0.25d0 ; grav=cte1/rstar**2
	pt_rac=EXP(bp_atm(1,1)) ; t_rac=EXP(bp_atm(2,1))
	m_rac=bp_atm(5,1)
	IF(pturb)THEN   !avec pression turbulente
	 p_rac=EXP(bp_atm(Ipgt,1))
	ELSE        !sans pression turbulente
	 p_rac=pt_rac
	ENDIF
c	PRINT*,'pt_rac,t_rac,m_rac,p_rac'
c	WRITE(*,2000)pt_rac,t_rac,m_rac,p_rac ; PAUSE'après coll_atm'

	IF(list)THEN
c	IF(.TRUE.)THEN
	 DO i=1,n_atm
	  CALL bsp1dn(ne_atm,bp_atm,x_atm,xt_atm,n_atm,ord_atm,knot_atm,
	1 .TRUE.,x_atm(i),ll,fx,dfxdx)
	  IF(no_croiss)PRINT*,'Pb. at 2 in lim_atm'	
c	  WRITE(*,2000)(EXP(fx(j)),j=1,2),(fx(j),j=3,5),(EXP(fx(j)),
c	  j=6,ne_atm)
	  pt_atm(i)=EXP(fx(1)) ; t_atm(i)=EXP(fx(2))
	  r_atm(i)=fx(3) ; m_atm(i)=fx(5) ; tau(i)=EXP(fx(7))
	  IF(pturb)THEN     !avec pression turbulente
	   p_atm(i)=EXP(fx(Ipgt)) ; dlpp_atm(i)=dfxdx(Ipgt)/dfxdx(1)
	  ELSE          !sans pression turbulente
	   p_atm(i)=pt_atm(i) ; dlpp_atm(i)=1.d0
	  ENDIF
c	  WRITE(*,2000)tau(i),pt_atm(i),p_atm(i),t_atm(i),r_atm(i),
c	  m_atm(i)
	 ENDDO
c	 PAUSE'lim_atm'

	 OPEN(unit=31,form='unformatted',status='unknown',
	1 file=TRIM(nom_fich2)//'_B.atm')
	 WRITE(31)dim_atm,n_atm,knot_atm,bp_atm,x_atm,xt_atm	
	 CLOSE(unit=31)

	ELSE        !dérivées / r et l
	 SELECT CASE(langue)
	 CASE('english')
	  PRINT*,'numerical partial derivatives / r et l'	 
	 CASE DEFAULT	   
	  PRINT*,'dérivées partielles numériques / r et l'
	 END SELECT	 
	 rd=r_rac*unpdx ;  dr=rd-r_rac
	 CALL coll_atm(rd,l_rac,xchim,ord_atm,knot_atm,dim_atm,ms_atm)     
	 dptsdr=(EXP(bp_atm(1,1))-pt_rac)/dr
	 dtsdr=(EXP(bp_atm(2,1))-t_rac)/dr ; dmsdr=(bp_atm(5,1)-m_rac)/dr
	 IF(pturb)THEN      !avec pression turbulente
	  dpsdr=(EXP(bp_atm(Ipgt,1))-p_rac)/dr
	 ELSE           !sans pression turbulente
	  dpsdr=dptsdr
	 ENDIF

	 ld=l_rac ; ld=ld*unpdx ; dl=ld-l_rac
	 CALL coll_atm(r_rac,ld,xchim,ord_atm,knot_atm,dim_atm,ms_atm)         
	 dptsdl=(EXP(bp_atm(1,1))-pt_rac)/dl
	 dtsdl=(EXP(bp_atm(2,1))-t_rac)/dl ; dmsdl=(bp_atm(5,1)-m_rac)/dl
	 IF(pturb)THEN      !avec pression turbulente
	  dpsdl=(EXP(bp_atm(Ipgt,1))-p_rac)/dl
	 ELSE           !sans pression turbulente 7 inconnues
	  dpsdl=dptsdl
	 ENDIF
	ENDIF

c	IF(.TRUE.)THEN
	IF(.FALSE.)THEN
	 DO i=1,n_atm
	  CALL bsp1dn(ne_atm,bp_atm,x_atm,xt_atm,n_atm,ord_atm,knot_atm,
	1 .TRUE.,x_atm(i),ll,fx,dfxdx)
	  IF(no_croiss)PRINT*,'Pb. at 3 in lim_atm'	
	  WRITE(*,2000)(EXP(fx(j)),j=1,2),(fx(j),j=3,5),(EXP(fx(j)),
	1 j=6,ne_atm)
	 ENDDO
	 CALL pause('solution atm')
	ENDIF
	SELECT CASE(langue)
	CASE('english')
	 WRITE(*,1005)
1005	 FORMAT('------- Restauration of the atmosphere (end) ------')	
	CASE DEFAULT
	 WRITE(*,5)
5	 FORMAT('------- Restitution de l''atmosphère (fin) ------')
 	END SELECT
		
	RETURN

	END SUBROUTINE lim_atm
