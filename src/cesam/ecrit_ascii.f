
	SUBROUTINE ecrit_ascii
	
c subroutine subordonnée de cesam.f introduite par INCLUDE
c créé et écrit les fichiers ACII de type output

c ATTENTION
c ne pas mettre en 1-ère ligne c***********  c'est refusé par le compilateur???
c en_masse=.FALSE. (en rayon) n'est pas toujours prévu

c	Auteur: P.Morel, Département Cassiopée, O.C.A.
c	CESAM2k

c-----------------------------------------------------------------------	

	USE mod_kind
	USE mod_numerique, ONLY : delete_doubles, shell
	USE mod_variables, ONLY : tot_conv, tot_rad
	
	IMPLICIT NONE
	
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: compgo
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: cpo, dcapdro,
	1 dcapdto, deltao, depsdro, depsdto, epsilono, gammao, grado, gradado,
	2 gradrado, grad_mjo, grad_muo, kapo, lo, mo, muo, mueo, po, pto, q_out,
	3 ryo, roo, temp, to, vaissalao, wo

	real (kind=dp), allocatable, dimension(:) :: epsilonnuc, epsilongrav !YLD

	REAL (kind=dp), DIMENSION(nchim) :: dxchimg, xchimg, z_bar
	REAL (kind=dp) :: ad, add, alfa, beta, degene, dlpp, hp, nel, qad, qo,
	1 u

	INTEGER :: n_out	

	LOGICAL, SAVE :: init=.TRUE.
	LOGICAL :: convec	
	
c-----------------------------------------------------------------------	
	
2000	FORMAT(8es10.3)
	
	IF(init)THEN
	 iglob=19		!nb de global
	 ALLOCATE(glob(iglob)) ; init=.FALSE.
	ENDIF
	
c	PAUSE'ecrit_ascii'

c abscisses q_out des points ajouter autour des limites ZR/ZC
c ajout de l0 points entre lim_ZC et lim_ZC+1
c ajout de l0-1 points entre lim_ZC+1 et lim_ZC+2
c ajout de l0-2 points entre lim_ZC+2 et lim_ZC+3
c..............................
c ajout de 1 point entre lim_ZC+l0-1 et lim_ZC+l0
c et par symétrie entre lim_ZC et lim_ZC-l0
c avec l0=3 on ajoute 4, 3, 2, 1 points dans chaque intervalle à partir de
c q(jlim(i)), i=1,..,lim

c ajout d'un point à q0 du centre

c ajout de l0 points dans la couche 1
c ajout de l0-1 points dans la couche 2 
c ajout de l0-2 points dans la couche 3
c..............................
c ajout de 1 point dans la couche l0 

c ajout total : 1 (q0), l0(l0+1)/2 (au début), l0(l0+1)*lim (autour des limites)

c on réduit l'extension des ajouts s'il y a risque de débordement
c au centre et à la surface
	IF(tot_conv .OR. tot_rad)THEN
	 IF(q0 > 0.d0)THEN
	 
c un point à q0 du centre 
	  n_out=n_qs+1 ; ALLOCATE(q_out(n_out))
	  q_out(1:n_qs)=q(1:n_qs) ; q_out(n_out)=1.d0+q0

	 ELSE
	  n_out=n_qs ; ALLOCATE(q_out(n_out))
	  q_out(1:n_qs)=q(1:n_qs)
	 ENDIF
	ELSE
	
c un point à q0 du centre	
	 IF(q0 > 0.d0)THEN
	  n_out=n_qs+(l0+1)*l0/2*(2*lim+1)+1
	  ALLOCATE(q_out(n_out))	
	  q_out(1:n_qs)=q(1:n_qs) ; q_out(n_out)=1.d0+q0
c	  WRITE(*,2000)q_out(1:n_qs)
	 ELSE
	  n_out=n_qs+(l0+1)*l0/2*(2*lim+1) 
	  ALLOCATE(q_out(n_out))	
	  q_out(1:n_qs)=q(1:n_qs) !; WRITE(*,2000)q_out(1:n_qs)
	 ENDIF

c de part et d'autre des limites
	 j=n_qs
c	 PRINT*,'n_qs,n_out,l0,lim' ; PRINT*,n_qs,n_out,l0,lim
c	 PRINT*,'q(jlim)' ; WRITE(*,2000)q(jlim(1:lim))
	 DO k=1,lim
	  qo=q(jlim(k))
	  DO lq=0,l0-1
	   add=REAL(l0+1-lq,dp)
	   DO i=1,l0-lq
	    ad=REAL(i,dp)/add+REAL(lq,dp)
	    qad=qo+ad
	    IF(qad < n_qs)THEN
	     j=j+1 ; q_out(j)=qad
c	     PRINT*,'k,lq,i,j' ; PRINT*,k,lq,i,j
c	     PRINT*,'qo,add,ad,qad' ; WRITE(*,2000)qo,add,ad,qad	    
	    ENDIF
	    qad=qo-ad 
	    IF(qad > 1)THEN
	     j=j+1 ; q_out(j)=qad
c	     PRINT*,'qo,ad,qad' ; WRITE(*,2000)qo,ad,qad ; PAUSE'-'
	    ENDIF
	   ENDDO
	  ENDDO	
	 ENDDO
	 
c dans les l0 premières couches	 
	 qo=q(1)
	 DO lq=0,l0-1
	  add=REAL(l0+1-lq,dp)
	  DO i=1,l0-lq
	   ad=REAL(i,dp)/add+REAL(lq,dp)
	   qad=qo+ad ; j=j+1 ; q_out(j)=qad
c	   PRINT*,'k,lq,i,j' ; PRINT*,k,lq,i,j
c	   PRINT*,'qo,add,ad,qad' ; WRITE(*,2000)qo,add,ad,qad ; PAUSE'0'    
	  ENDDO
	 ENDDO	
	 
c s'il y a un point de plus au centre 
	 IF(q0 > 0.d0)j=j+1

c il ne faut pas trop de points	
	 IF(j > n_out)THEN
c ERREUR : trop de points	
	  PRINT*,'STOP, erreur dans ecrit_ascii : j =',j,' > ',n_out,'= n_out'
	  STOP
	 
c moins de points à cause de la proximité des limites réallocation de q_out	 
	 ELSEIF(j < n_out)THEN
	  IF(q0 > 0.d0)q_out(j)=q_out(n_out)
	  ALLOCATE(temp(j)) ; temp(1:j)=q_out(1:j) ; DEALLOCATE(q_out)
	  ALLOCATE(q_out(j)) ; q_out=temp ; DEALLOCATE(temp) ; n_out=j
	 ENDIF	
	ENDIF
	
c tri par ordre croissant
	CALL shell(n_out,q_out)
c	PRINT*,n_out,n_qs,lim,jlim(1:lim) ; WRITE(*,2000)q_out ; PAUSE'q_out'

c suppression d'abscisses égales
	CALL delete_doubles(n_out,q_out)
c	PRINT*,n_out ; PAUSE'doubles'
			    
c formation des global	    	    
	glob(1)=mstar*msol		!Masse a l'intérieur de R*
	glob(2)=rstar*rsol		!R* rayon a R(tau*)
	glob(3)=l(n_qs)*lsol		!Luminosité
	glob(4)=z0			!Z initial
	glob(5)=x0			!X initial
	glob(6)=alpha		!paramètre de longueur de melange
	glob(7)=9.d0/4.d0	!paramètre de convection arbitraire
	glob(8)=1.d0/162.d0	!paramètre de convection arbitraire
	glob(9)=compg(1,n_qs)		!X dans ZC
	IF(ihe4 > 1)THEN	!Y dans ZC
	 glob(10)=compg(ihe4,n_qs)+compg(ihe4-1,n_qs)
	ELSE
	 glob(10)=1.d0-glob(9)-z0
	ENDIF
	glob(11)=d2p		!dérivée seconde au centre
	glob(12)=d2ro		!dérivée seconde au centre
	glob(13)=age		!durée de l'évolution
	glob(14)=wrot 	!vitesse de rotation globale
	glob(15)=w_rot	!vitesse de rotation initiale
	glob(16)=g	!constante de la gravitation utilisée
	glob(17)=msol
	glob(18)=rsol
	glob(19)=lsol
	  
c	WRITE(*,2000)glob ; PAUSE'glob'

c formation de la partie atmosphère	   
!	nadd=max(0,n_atm-1) ; itot=n_out+nadd ; ivar=46
	nadd=max(0,n_atm-1) ; itot=n_out+nadd ; ivar=48  ! YLD ajout eps_nuc, eps_grav
	ALLOCATE(var(ivar+nchim,itot))
c	PRINT*,'q0,itot,n_out,nadd',q0,itot,n_out,nadd ; PAUSE'itot'	   
	   
	IF(nadd > 0)THEN	!on ajoute l'atmosphère sauf la couche 1
	 DO i=1,nadd
          ipn=n_atm-i+1
          var(1,i)=r_atm(ipn)*rsol      	!Rayon
          var(2,i)=LOG(m_atm(ipn)/mstar)	!Ln M/Mtot
          var(3,i)=t_atm(ipn)           	!Température
          var(4,i)=pt_atm(ipn)          	!Pression
          var(5,i)=ro_atm(ipn)          	!Densité
          var(6,i)=grad_mj_a(ipn)       	!Gradient
          var(7,i)=grad_atm(ipn)        	!Gradient
	  var(8,i)=l(n_qs)*lsol         	!Luminosité
          var(9,i)=k_atm(ipn)           	!Opacité
          var(10,i)=0.d0                	!Energie nuc+grav
          var(11,i)=gamma_atm(ipn)      	!Grand Gamma1
          var(12,i)=grada_atm(ipn)      	!Gradient adiabatique
          var(13,i)=delta_atm(ipn)      	!Delta
          var(14,i)=cp_atm(ipn)         	!Cp
          var(15,i)=mue_atm(ipn)        	!Mue
          var(16,i)=mu_atm(ipn)         	!Mu      
          var(17,i)=vais_atm(ipn)       	!Vaissala
          var(18,i)=w(n_qs)             	!Omega
          var(19,i)=ldcapdt_a(ipn)*t_atm(ipn)/k_atm(ipn) !dln kap/dln T
          var(20,i)=ldcapdr_a(ipn)*ro_atm(ipn)/k_atm(ipn)!dln kap/dln ro
          var(21,i)=0.d0                	!d epsilon(nuc) / d ln T
          var(22,i)=0.d0                	!d epsilon(nuc) / d ln ro
	  IF(pturb)THEN
           var(23,i)=pt_atm(ipn)/p_atm(ipn)	!Ptot/Pgaz
	  ELSE
	   var(23,i)=0.d0			!grad_mu dans atmosphère
	  ENDIF
          var(24,i)=gradr_atm(ipn)		!Gradient radiatif       
          DO j=1,nchim
           var(ivar+j,i)=xchim1g(j)		!Abondances / gramme
          ENDDO
         var(47:48,i)=0.d0    ! YLD   energy nuc ; energy grav    
         ENDDO
	ENDIF		!atmosphère	
c	PRINT*,'n_atm,nadd,n_out,n_qs,itot',n_atm,nadd,n_out,n_qs,itot
c	PRINT*,var(1,1),var(1,nadd)	

c les variables de la structure interne	
	ALLOCATE(cpo(n_out),compgo(nchim,n_out),dcapdro(n_out),dcapdto(n_out),
	1 deltao(n_out),depsdro(n_out),depsdto(n_out),epsilono(n_out),
	2 gammao(n_out),grado(n_out),gradado(n_out),gradrado(n_out),
	3 grad_mjo(n_out),grad_muo(n_out),kapo(n_out),lo(n_out),mo(n_out),
	4 muo(n_out),mueo(n_out),po(n_out),pto(n_out),roo(n_out),ryo(n_out),
	5 to(n_out), vaissalao(n_out),wo(n_out))

      allocate(epsilonnuc(n_out), epsilongrav(n_out)) ! YLD
c calcul des variables aux abscisses q_out
c utilisation des valeurs déjà calculées dans cesam en dehors des ajouts 
	j=1
	DO i=1,n_out                  !extraction p,t,r,l,m
	 IF(q_out(i) == q(j))THEN
c points déjà calculés
	  ryo(i)=r(j)	  	  
	  mo(i)=m(j)
	  to(i)=t(j)
	  po(i)=p(j)
	  pto(i)=pt(j)
	  roo(i)=ro(j)
	  grad_mjo(i)=grad_mj(j)
	  grad_muo(i)=grad_mu(j)
	  grado(i)=grad(j)
	  lo(i)=l(j)
	  kapo(i)=kap(j)
	  epsilono(i)=epsilon(1,j) 

	      epsilonnuc(i)=epsilon(2,j)+epsilon(3,j)+epsilon(4,j)  !YLD nuclear energy
	      epsilongrav(i)=-epsilon(5,j)                          ! YLD gravitation energy

	  gammao(i)=gamma(j)
	  gradado(i)=gradad(j)
	  deltao(i)=delta(j)
	  cpo(i)=cp(j)
	  mueo(i)=mue(j)
	  muo(i)=mu(j)
	  vaissalao(i)=vaissala(j)
	  wo(i)=w(j)
	  dcapdto(i)=dcapdt(j)
	  dcapdro(i)=dcapdr(j)
	  depsdto(i)=depsdt(j)
	  depsdro(i)=depsdr(j)
	  gradrado(i)=gradrad(j)
	  compgo(:,i)=compg(:,j)   
	  j=MIN(n_qs,j+1)	  	
	 ELSE
	 
c points ajoutés
c	  PRINT*,'point ajouté,i,j ',i,j
	  CALL bsp1ddn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q_out(i),lq,fd)
	  IF(no_croiss)PRINT*,'Pb. at 1 in ecrit_ascii'
	  pto(i)=EXP(fd(1,0))    	!Ptot        
	  IF(pturb)THEN         	!avec pression turbulente 8 inconnues
	   po(i)=EXP(fd(Ipg,0))	!Pgas
	   dlpp=fd(Ipg,1)/fd(1,1)	!dln Pgaz/dln Ptot
	  ELSE                  !sans pression turbulente 7 inconnues
	   po(i)=pto(i)
	   dlpp=1.d0		!dln Pgaz/dln Ptot
	  ENDIF      
	  to(i)=EXP(fd(2,0))     !variable ln T
	  IF(en_masse)THEN
	   ryo(i)=SQRT(ABS(fd(3,0)))		!rayon/Rsol
	   lo(i)=SQRT(ABS(fd(4,0)))**3		!l/Lsol
	   mo(i)=SQRT(ABS(fd(5,0)))**3		!m/Msol
	   IF(fd(3,0) > 0.d0)THEN
	    grav=cte1*mo(i)/fd(3,0)		!fd(3,0)=r^2         
	   ELSE
	    grav=0.d0
	   ENDIF 
	  ELSE       
	   ryo(i)=ABS(fd(3,0))		!rayon/Rsol
	   lo(i)=fd(4,0)	     	!l/Lsol
	   mo(i)=ABS(fd(5,0))    	!m/Msol
	   fd(5,0)=mo(i)**(2.d0/3.d0)
	   IF(fd(3,0) > 0.d0)THEN
	    grav=cte1*mo(i)/fd(3,0)**2
	   ELSE
	    grav=0.d0
	   ENDIF
	  ENDIF	   
	  	       
c le gradient réel pour MJo
	  grad_mjo(i)=fd(2,1)/fd(1,1)
c	  WRITE(*,2000)p(i),t(i),r(i),l(i),m(i),grad_mj(i)
c	  WRITE(*,2000)r(i),r2(i),m(i),m23(i)
	 
c la composition chimique
	  CALL bsp1dn(nchim,chim,mc,mct,n_ch,m_ch,
	1 knotc,.TRUE.,MAX(mc(1),MIN(fd(5,0),mc(n_ch))),lq,xchim,dxchim)
	  IF(no_croiss)PRINT*,'Pb. at 17 in ecrit_ascii'
	  IF(fd(5,0) > 0.d0)THEN
	   dxchimg=dxchim*2.d0/3.d0/SQRT(fd(5,0))
	  ELSE
	   dxchimg=0.d0
	  ENDIF		  	
	  xchimg=xchim ; CALL chim_gram(xchimg,dxchimg)
          compgo(:,i)=xchimg(:)

c la rotation
	  SELECT CASE(Krot)
	  CASE(0,1,2)
	   wo(i)=wrot
	  CASE(3,4,5)
	   CALL bsp1dn(nrot,rota,mrot,mrott,n_rot,ord_rot,
	1  knotr,.TRUE.,MAX(mrot(1),MIN(fd(5,0),mrot(n_rot))),lq,frot,dfrot)
	   IF(no_croiss)PRINT*,'Pb. at 18 in cesam'	
	   wo(i)=frot(1)
	  END SELECT

c la thermo +	
	  CALL thermo(pto(i),po(i),to(i),mo(i),lo(i),ryo(i),dlpp,xchim,dxchim,
	1 roo(i),drop,drot,drox,u,dup,dut,dux,
	2 grado(i),dgradpt,dgradp,dgradt,dgradx,dgradm,dgradl,dgradr,
	3 dgradlpp,gam,dgampt,dgamp,dgamt,dgamx,dgamm,dgaml,dgamr,
	4 dgamlpp,epsilo,depsp,depst,depsx,kapo(i),dkapp,dkapt,dkapx,
	5 deltao(i),deltap,deltat,deltax,cpo(i),dcpp,dcpt,dcpx,
	6 gradado(i),dgradadp,dgradadt,dgradadx,
	7 hp,dhppt,dhpp,dhpt,dhpx,dhpr,dhpm,grad_muo(i),
	8 gradrado(i),alfa,beta,gammao(i),convec)	

c dérivées de l'énergie nucléaire
	  depsdro(i)=depsp/drop ; depsdto(i)=depst-depsdro(i)*drot
	  
c ajout du TdS	  
	  CALL bsp1dn(1,tds,x_tds,xt_tds,n_tds,m_tds,knot_tds,.TRUE.,
	1 MAX(x_tds(1),MIN(fd(5,0),x_tds(n_tds))),lq,ftds,dftds)
	  IF(no_croiss)PRINT*,'Pb. at 19 in ecrit_ascii'
	  epsilono(i)=epsilo(1)-ftds(1)/secon6
	      epsilonnuc(i)=epsilo(1)  !YLD nuclear energy
	      epsilongrav(i)=-ftds(1)/secon6 ! YLD gravitation energy
	  	   
c dérivées de l'opacité
	  dcapdro(i)=dkapp/drop ; dcapdto(i)=dkapt-dcapdro(i)*drot
   
c estimation des taux d'ionisation et calcul des poids
c moléculaires moyens, si T < 4500K on suppose les éléments
c totalement recombinés (saha ne converge pas)
	  IF(to(i) < 4.5d3)THEN
	   z_bar=0.d0 ; mueo(i)=1.d30
	  ELSE
	   IF(mu_saha)THEN	   	   
	    CALL saha(xchim,to(i),roo(i),ioni,z_bar,nel,degene)
	    mueo(i)=1.d0/DOT_PRODUCT(z_bar,xchim)
	   ELSE
	    z_bar=zi
	   ENDIF
	   mueo(i)=1.d0/DOT_PRODUCT(z_bar,xchim)   	    
	  ENDIF
	  muo(i)=1.d0/DOT_PRODUCT((1.d0+z_bar),xchim)
	   
c calcul d'une des formes de Vaissala  (Kippenhahan & Weigert p. 42)
c 1/gamma1 dlnP/dln r - dln ro/dln r
	  IF(hp <= 0.d0)THEN            !au centre
	   vaissalao(i)=0.d0
	   
c avec mu0 ou utilisation de la tabulation vth	   
c fvth(1)=lnP, fvth(2)=lnT, fvth(3)=r**2, fvth(4)=l**2/3, fvth(5)=ln ro, 
c fvth(6)= cp, fvth(7)=delta, fvth(8)=gamma1, fvth(9)=ln µ, fvth(10)=ln kap
	  ELSEIF(new_bv)THEN
c	   vaissalao(i)=ryo(i)*rsol/hp*(alfa*(dlpp-1.d0)
c	1  +deltao(i)*(gradado(i)-grad_mjo(i)
c	2  +beta/(4.d0-3.d0*beta)*grad_muo(i)))

	   CALL bsp1dn(nvth,vth,mc,mct,n_ch,m_ch,
	1  knotc,.TRUE.,fd(5,0),lq,fvth,dfvth)
	   IF(no_croiss)PRINT*,'Pb. at 20 in ecrit_ascii'		
	   vaissalao(i)=(dfvth(1)/fvth(8)-dfvth(5))*fvth(3)/dfvth(3)*2.d0

c ancienne formulation
c r / Hp delta ( grad_ad - grad ) - 4 pi r^3 dro/dX  dX/dm
c le 4 pi r^3 dro/dX  dX/dm est une approx. de phi grad_mu		
	  ELSE	  
	   vaissalao(i)=ryo(i)*rsol/hp*(alfa*(dlpp-1.d0)
	1  +deltao(i)*(gradado(i)-grad_mjo(i)))-cte3*ryo(i)**3*drox*dxchimg(1)
	  ENDIF
	 ENDIF
	ENDDO

c formation de la partie structure interne
	DO i=n_out,1,-1
	 ipn=n_out-i+nadd+1
	 var(1,ipn)=ryo(i)*rsol			!Rayon
	 IF(mo(i) > 0.d0)THEN
	  var(2,ipn)=LOG(mo(i)/mstar)		!Ln M/Mtot
	 ELSE
	  var(2,ipn)=-1.d38		        !Ln M/Mtot
	 ENDIF
	 var(3,ipn)=to(i)			!Température
	 var(4,ipn)=pto(i)			!Pression
	 var(5,ipn)=roo(i)			!Densité
	 var(6,ipn)=grad_mjo(i) 		!Gradient
	 var(7,ipn)=grado(i)			!Gradient	    
	 var(8,ipn)=lo(i)*lsol			!Luminosité
	 var(9,ipn)=kapo(i)			!Opacité
	 var(10,ipn)=epsilono(i)		!Energie nuc+grav
	 var(11,ipn)=gammao(i)			!Grand Gamma1
	 var(12,ipn)=gradado(i)			!Gradient adiabatique
	 var(13,ipn)=deltao(i)			!Delta
	 var(14,ipn)=cpo(i)			!cp
	 var(15,ipn)=mueo(i)			!Mue
	 var(16,ipn)=muo(i)			!Mu	    
	 var(17,ipn)=vaissalao(i)		!Vaissala
	 var(18,ipn)=wo(i)			!Omega
	 var(19,ipn)=dcapdto(i)*to(i)/kapo(i)	!dln kappa / dln T
	 var(20,ipn)=dcapdro(i)*roo(i)/kapo(i)	!dln kappa / dln ro
	 var(21,ipn)=depsdto(i)*to(i)		!d epsilon(nuc) / d ln T
	 var(22,ipn)=depsdro(i)*roo(i)		!d epsilon(nuc) / d ln ro
	 IF(pturb)THEN
	  var(23,ipn)=pto(i)/po(i)		!Ptot/Pgaz
	 ELSE
	  var(23,ipn)=grad_muo(i)		!grad_mu
	 ENDIF
	 var(24,ipn)=gradrado(i)		!Gradient radiatif     
	 DO j=1,nchim
	  var(ivar+j,ipn)=compgo(j,i)		!Abondances / gramme
	 ENDDO
	   var(47,ipn)=epsilonnuc(i) !YLD nuclear energy
	   var(48,ipn)=epsilongrav(i) !YLD gravitation energy
	ENDDO
	
	DEALLOCATE(cpo,compgo,dcapdro,dcapdto,deltao,depsdro,depsdto,
	1 epsilono,gammao,grado,gradado,gradrado,grad_mjo,grad_muo,kapo,lo,mo,
	2 muo,mueo,po,pto,q_out,roo,ryo,to,vaissalao,wo)	

      deallocate(epsilonnuc, epsilongrav) !YLD
	
c on complète les données pour les sorties ASCII	   
	CALL add_ascii(var,glob,itot,ivar)
	
c écritures du fichier ASCII
	CALL output(var,glob,itot,ivar)

c le nombre de couches pouvant changer	
	DEALLOCATE(var) ; IF(sort)DEALLOCATE(glob)
	    	    
	RETURN
	
	END SUBROUTINE ecrit_ascii	   
