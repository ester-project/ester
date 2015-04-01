
c*************************************************************************

	SUBROUTINE rq_reac(comp,t,ro,r,drt,dro,drx,q,dqt,dqo,dqx,mue,dmuex)

c routine private du module mod_nuc

c calcul des R et Q des réactions thermonucléaires
c et de l'effet d'écran, on utilise la tabulation tabul_nuc

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c références :
c	effet d'écran faible ou fort de Mitler ApJ 212, 513, 1977
c	mitler=.TRUE. : écrantage de Mitler, sinon écran faible

c entrées :
c	comp : abondances en mole
c	t : température
c	ro : densité

c sorties :
c	r(i) : taux de réaction g /s /mole/mole (/mole)
c	q(i): énergie ( erg /sec /gr )
c	mue : 1./masse moléculaire moyenne par électron

c--------------------------------------------------------------------

	USE mod_donnees, ONLY : amu, echarg, eve, kbol, hpl, me, mitler,
	1 nchim, nucleo, pi, t_inf, zi
	USE mod_numerique, ONLY : bsp1dn	
	USE mod_kind
	USE mod_variables, ONLY : sortie

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: comp
	REAL (kind=dp), INTENT(in) :: t, ro
	REAL (kind=dp), INTENT(out), DIMENSION(:,:) :: drx, dqx
	REAL (kind=dp), INTENT(out), DIMENSION(:) :: r, drt, dro, q,
	1 dqt, dqo, dmuex	
	REAL (kind=dp), INTENT(out) :: mue

	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: decrx, dzetax
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: datx, dbidx,
	1 dc1x, decro, decrt, dfx, df12x, dkx, dnex, dr, dzstarx, ecran,
	2 lambda, mz1z2, tx, z1z2
		
	REAL (kind=dp), SAVE, DIMENSION(2) :: zeta, dzetat, dzetaro	
	REAL (kind=dp), SAVE :: cte1, cte2, cte3, cte4, cte5, cte6, cte7,
	1 cte8, cte9, cte10
	REAL (kind=dp) :: zstar, ne, c1, f12, at, k, bid, bid1, df,
	1 e0, eps, beta, dnero, dc1ro, df12t, df12ro, datt, datro, dkt,
	2 dkro, dbidt, dbidro, dft, dfro
	
	INTEGER, SAVE :: l=1
	INTEGER :: i, j

	LOGICAL, SAVE :: init=.TRUE.
	
c-------------------------------------------------------------------

2000	FORMAT(8es10.3)
2001	FORMAT(10es8.1)
2002	FORMAT(20i4)

	IF(init)THEN	!en cas de modifications assurer la cohérence
	 init=.FALSE.		!avec reac_nuc

c	 PRINT*,nchim ; PAUSE'rq_reac entrée'

c tabulation des réactions nucléaires
	 CALL tabul_nuc	 

c écritures diverses
	 IF(mitler)THEN
	  WRITE(2,4) ; WRITE(*,4)
4         FORMAT('écrantage de Mitler')
	 ELSE
	  WRITE(2,5) ; WRITE(*,5)
5         FORMAT('écrantage faible')
	 ENDIF
	 	 	 
c	 WRITE(*,2000)q0 ; PAUSE'les q'
	 
	 WRITE(2,1)t_inf ; WRITE(*,1)t_inf
1	 FORMAT('Température d''amorçage des réactions nucléaires:',
	1 es10.3)
 	 WRITE(2,2)t_sup ; WRITE(*,2)t_sup
2	 FORMAT('Température maximale pour les réactions nucléaires:',
	1 es10.3,/)
	
c quelqes constantes
	 cte1=2.d0/3.d0*pi*echarg ; cte2=4.d0*pi*echarg**2/kbol
	 cte3=hpl**3/4.d0/pi/SQRT(2.d0*me*kbol)**3*3.d0/2.d0
	 cte4=3.d0/4.d0/pi ; cte5=-18.d0/5.d0
	 cte6=1.2204d0*(1.d-6)**(2.d0/3.d0)*eve*1.d3	!3', Lang 4-44 A en MeV
	 cte7=5.d0/3.d0 ; cte8=0.001037d0
	 cte9=7.d0*5.47d2/9.6d2 ; cte10=5.d0/8.d0
	 
c quantités à conserver
	 ALLOCATE(z1z2(nreac),lambda(nreac),mz1z2(nreac))
	 
	 DO i=1,nreac
	  z1z2(i)=izz(i,1)*izz(i,2)
	  mz1z2(i)=cte6*(z1z2(i)**2*ar(i))**(1.d0/3.d0)	!Lang 4-44	 
	  lambda(i)=z1z2(i)*echarg**2		!eq. 30'
	 ENDDO
	 
c autres allocations 
	 ALLOCATE(ecran(nreac),dr(nreac),decrt(nreac),
	1 decro(nreac),dzetax(2,nchim),decrx(nreac,nchim),dzstarx(nchim),
	2 dc1x(nchim),dnex(nchim),df12x(nchim),datx(nchim),dkx(nchim),
	3 dbidx(nchim),dfx(nchim),tx(nchim))
	 
c	 PAUSE'sortie init'
	 
	 RETURN	!le premier appel est un appel d'initialisation
     
	ENDIF

	IF(t < t_inf)THEN
	 r=0.d0 ; q=0.d0 ; drt=0.d0 ; dro=0.d0
	 q=0.d0 ; dqt=0.d0 ; dqo=0.d0
	 RETURN
	ELSEIF(t > t_sup)THEN
	 WRITE(*,3)t,t_sup
3	 FORMAT('T=',es10.3,' > T_sup='es10.3,
	1', hors de la tabulation des réactions nucléaires')
	 WRITE(2,3)t,t_sup ; CALL sortie	!STOP	 	 
	ENDIF

c effet d'écran
	IF(mitler)THEN	!écran Mitler

	 mue=0.d0	!1/poids moléculaire par électron
	 zstar=0.d0
	 DO i=1,nchim	!on néglige les e- des éléments à l'équilibre
	  mue=mue+zi(i)*comp(i) ; dmuex(i)=zi(i)
	  zstar=zstar+zi(i)**2*comp(i)	!eq. 6"
	  dzstarx(i)=zi(i)**2
	 ENDDO
	 zstar=zstar/mue		!eq. 6"
	 dnero=mue/amu	 
	 ne=ro*dnero	!nb. électrons, mue=1/poids mol. par électron
	 c1=cte1*ne				!eq. 17'
	 dc1ro=cte1*dnero
	 DO i=1,nchim
	  dzstarx(i)=dzstarx(i)/mue-zstar/mue*dmuex(i)
	  dnex(i)=ne/mue*dmuex(i) ; dc1x(i)=cte1*dnex(i)
	 ENDDO
	 f12=cte3*ne/SQRT(t)**3		!Clayton eq. 2.59 X 3/2
	 df12t=-3.d0/2.d0*f12/t ; df12ro=f12/ne*dnero
	 DO i=1,nchim
	  df12x(i)=f12*dnex(i)/ne
	 ENDDO	 
	 at=4.d0/9.d0*f12**(4.d0/3.d0) ; datt= 4.d0/3.d0*at/f12*df12t
	 datro=4.d0/3.d0*at/f12*df12ro
	 DO i=1,nchim
	  datx(i)=4.d0/3.d0*at/f12*df12x(i)
	 ENDDO	 	 
	 at=1.d0/SQRT(at+1.d0)			!eq. 22"
	 datt=  -datt/2.d0*at**3		!eq. 22"
	 datro=-datro/2.d0*at**3		!eq. 22"
	 DO i=1,nchim
	  datx(i)=-datx(i)/2.d0/at**3
	 ENDDO	 	 
	 k=cte2*ne/t*(zstar+at) ; dkt= k*(datt /(zstar+at)-1.d0/t)
	 dkro=k*(datro/(zstar+at)+dnero/ne)
	 DO i=1,nchim
	  dkx(i)=k*((dzstarx(i)+datx(i))/(zstar+at)+dnex(i)/ne)
	 ENDDO	 	 	 
	 k=SQRT(k)				!eq. 22'
	 dkt=dkt/2.d0/k ; dkro=dkro/2.d0/k	
	 DO i=1,nchim
	  dkx(i)=dkx(i)/2.d0/k
	 ENDDO	 	 
	  
c	 PRINT*,'zstar,mue,ne,c1,f12,at,k'	
c	 WRITE(*,2000)zstar,mue,ne,c1,f12,at,k
	
	 DO i=1,nreac
	  IF(z1z2(i) /= 0.d0)THEN	!pour passer les captures e-
	   DO j=1,2
	    zeta(j)=cte4*k**3/ne*izz(i,j)	!eq. 24
	    dzetat(j)= zeta(j) *3./k*dkt
	    dzetaro(j)=zeta(j)*(3./k*dkro-dnero/ne)
	    DO l=1,nchim
	     dzetax(j,l)=zeta(j)*(3./k*dkx(l)-dnex(l)/ne)
	    ENDDO	    
	   ENDDO
	   bid=1.d0
	   bid1=(zeta(1)+zeta(2)+1.d0)**cte7		!eq. 25
	   bid=bid+bid1
	   dbidt= bid1*cte7/(zeta(1)+zeta(2)+1.d0)*(dzetat(1) +dzetat(2))
	   dbidro=bid1*cte7/(zeta(1)+zeta(2)+1.d0)*(dzetaro(1)+dzetaro(2))
	   DO l=1,nchim
	    dbidx(l)=bid1*cte7/(zeta(1)+zeta(2)+1.d0)*(dzetax(1,l)
	1   +dzetax(2,l))
	   ENDDO	   
	   bid1=-(zeta(1)+1.d0)**cte7
	   bid=bid+bid1	   
	   dbidt= dbidt +bid1*cte7/(zeta(1)+1.d0)*dzetat(1)
	   dbidro=dbidro+bid1*cte7/(zeta(1)+1.d0)*dzetaro(1)
	   DO l=1,nchim	
	    dbidx(l)=dbidx(l)+bid1*cte7/(zeta(1)+1.d0)*dzetax(1,l)	   
	   ENDDO	   	   
	   bid1=-(zeta(2)+1.d0)**cte7
	   bid=bid+bid1		   
	   dbidt= dbidt +bid1*cte7/(zeta(2)+1.d0)*dzetat(2)
	   dbidro=dbidro+bid1*cte7/(zeta(2)+1.d0)*dzetaro(2)
	   DO l=1,nchim	
	    dbidx(l)=dbidx(l)+bid1*cte7/(zeta(2)+1.d0)*dzetax(2,l)	   
	   ENDDO	   	   
	   df=cte5*c1**2/k**5*bid		!eq. 25
	   dft= df*(-5.d0*dkt /k+dbidt /bid)   
	   dfro=df*(-5.d0*dkro/k+dbidro/bid+2.d0*dc1ro/c1)
	   DO l=1,nchim	
	    dfx(l)=df*(-5.d0*dkx(l)/k+dbidx(l)/bid+2.d0*dc1x(l)/c1)   
	   ENDDO
	   	      	   
	   ecran(i)=ABS(df) ; decrt(i)=SIGN(dft,df) 
	   decro(i)=SIGN(dfro,df)
	   DO l=1,nchim	
	    decrx(i,l)=SIGN(dfx(l),df)	   
	   ENDDO	   	   
	   
	   IF(.FALSE.)THEN	!suppression du terme d'ordre superieur	   
	    e0=t**(2.d0/3.d0)*mz1z2(i)	!eq. 3' i.e. Lang 4-44 E0 Gamov
	    bid1=(zeta(1)+zeta(2)+1.d0)**cte7-cte7*(zeta(1)+zeta(2))-1.d0!eq.32'
	    beta=cte8*df/lambda(i)/k*bid1/bid  !eq. 32'
	    eps=beta*(lambda(i)*k/(e0+df))**3	!eq. 38' 39 + remarque eq. 44
	    ecran(i)=ecran(i)-cte10*e0*eps*(1.d0+cte9*eps)	!eq. 58
	   ENDIF
	   ecran(i)=ecran(i)/kbol/t			!eq. 58	   
	   decrt(i)=decrt(i)/kbol/t-ecran(i)/t	
	   decro(i)=decro(i)/kbol/t
	   DO l=1,nchim
	    decrx(i,l)=decrx(i,l)/kbol/t	   
	   ENDDO		   
	   		      
	   ecran(i)=EXP(ecran(i))			!eq. 58
	   decrt(i)=ecran(i)*decrt(i) ; decro(i)=ecran(i)*decro(i)
	   DO l=1,nchim
	    decrx(i,l)=ecran(i)*decrx(i,l)	   
	   ENDDO		   
	   	   
c	   WRITE(*,2000)df,cte10*e0*eps*(1.d0+cte9*eps)
c	   PRINT*,'zeta(1),zeta(2),bid,df,e0,bid1,beta,eps'
c	   WRITE(*,2000)zeta(1),zeta(2),bid,df,e0,bid1,beta,eps
c	   PRINT*,'t,ecran(i),izz(i,1),izz(i,2)',izz(i,1),izz(i,2)
c	   WRITE(*,2000)t,ecran(i),EXP(.188*SQRT(ro*(comp(1)+3.)/2./
c	1  (t*1.d-6)**3))
c	   WRITE(*,2000)EXP(.188*SQRT(ro*(comp(1)+3.)/2./(t*1.d-6)**3))	  
c	   PAUSE
	  ELSE
	   ecran(i)=1.d0 !si z1 ou z2=0 par exemple pour Be7(e-,nu g)Li7
	   decrt(i)=0.d0 ; decro(i)=0.d0 ; decrx(i,:)=0.d0
	  ENDIF   
	 ENDDO
	 	  
	ELSE		!écran faible	
	 zstar=0.d0 ; mue=0.d0	
	 DO i=1,nchim
	  mue=mue+zi(i)*ABS(comp(i)) ; dmuex(i)=zi(i)
	  zstar=zstar+zi(i)*(zi(i)+1.d0)*ABS(comp(i))
	  dzstarx(i)=zi(i)*(zi(i)+1.d0)
	 ENDDO	 
	 DO i=1,nreac	  
	  IF(izz(i,1) > 0)THEN
	   ecran(i)=1.88d8*z1z2(i)*SQRT(zstar*ro/t**3)	!Clayton 4-221
	   decrt(i)=-3.d0/2.d0*ecran(i)/t ; decro(i)=ecran(i)/ro/2.d0	   
	   DO j=1,nchim
	    decrx(i,j)=ecran(i)/2.d0/zstar*dzstarx(j)
	   ENDDO	   	   
	   ecran(i)=EXP(ecran(i)) ; decrt(i)=ecran(i)*decrt(i) 
	   decro(i)=ecran(i)*decro(i)  
	   DO j=1,nchim
	    decrx(i,j)=ecran(i)*decrx(i,j)
	   ENDDO
	   
c photodissociations ou capture électroniques	   	   	   
	  ELSE
	   ecran(i)=1.d0 !si z1 ou z2=0 par exemple pour Be7(e-,nu g)Li7
	   decrt(i)=0.d0 ; decro(i)=0.d0 ; decrx(i,1:nchim)=0.d0
	  ENDIF
	 ENDDO
	ENDIF		!mitler
			
c interpolation en ln t
c	PRINT*,nreac,n_temp,m_temp,knot_temp,l
c	WRITE(*,2000)temp(1),ttemp(1),LOG(t)
c	WRITE(*,2000)taux_reac(1,1)
	CALL bsp1dn(nreac,taux_reac,temp,ttemp,n_temp,m_temp,knot_temp,
	1 .TRUE.,LOG(t),l,r,dr)
	
c limitation de la réaction 6 : Be7(e-,nu g)Li7 (Caughlan & Fowler)
	IF(t < 1.d6)THEN
	 bid=LOG(1.57d-7*2.d0/ro/(1.d0+comp(1)*nucleo(1)))	
	 r(6)=MIN(r(6),bid)
	ENDIF
	
c Rij et Q
	DO i=1,nreac
	 IF(i == i3al)THEN	 
	  r(i)=EXP(r(i))*ro**2*ecran(i)		!réaction 3 alpha
	  drt(i)=r(i)*(dr(i)/t+decrt(i)/ecran(i))
	  dro(i)=r(i)*(2.d0/ro+decro(i)/ecran(i))
	 ELSE	 
	  r(i)=EXP(r(i))*ro*ecran(i)
	  drt(i)=r(i)*(dr(i)/t+decrt(i)/ecran(i))
	  dro(i)=r(i)*(1.d0/ro+decro(i)/ecran(i))
	 ENDIF	
	 DO j=1,nchim
	  drx(i,j)=r(i)/ecran(i)*decrx(i,j)
	 ENDDO
	 q(i)=r(i)*q0(i) ; dqt(i)=drt(i)*q0(i) ; dqo(i)=dro(i)*q0(i)
	 DO j=1,nchim
	  dqx(i,j)=drx(i,j)*q0(i)
	 ENDDO
	ENDDO	!i
	
	RETURN

	END SUBROUTINE rq_reac
