
c***********************************************************************

	SUBROUTINE coll_atm(r_rac,l_rac,xchim,ord_atm,knot_atm,dim_atm,
	1 ms_atm)

c	routine PRIVATE du module mod_atm

c	résolution du système d'équations différentielles non linéaires
c	de l'atmosphère par développement sur B-splines
c	par itération Newton avec dérivées analytiques

c	modifs :
c	04 07 01 : mise en place de rotation uniforme avec conservation
c		   du moment angulaire
c	16 08 01 : F95

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c entrées
c	ray : rayon		!au point tau_f
c	lum : luminosite	!au point tau_f
c	xchim : composition chimique par gramme
c	wrot : vitesse angulaire solide

c sorties
c	bp_atm, x_atm, xt_atm, knot_atm: solution spline

c	le réseau de points de raccord x, possède p points
c	la solution, les coefficients des splines d'ordre m+r sur le
c	réseau de points de table xt, sont dans b

c	les initialisations, les résidus et les dérivées sont 
c	calculés dans eqatm

c	avec pression turbulente 8 inconnues 
c	sans pression turbulente 7 inconnues, Ptot=Pgaz

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : g, Krot, langue, msol, n_atm, pturb, rsol,
	1 r_qs
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn, bvald, coll, gauss_band, linf
	USE mod_variables, ONLY : dim_rot, if_conv, mstar, rota, rstar,
	1 r_zc_conv, sortie, wrot
		
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in), DIMENSION(:) :: xchim
	REAL (kind=dp), INTENT(in) :: r_rac, l_rac
	INTEGER, INTENT(in) :: dim_atm, ms_atm, ord_atm
	INTEGER, INTENT(inout) :: knot_atm
		 	
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: ae, dl, ds
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: a, b, y
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: d
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: be, dfxdx, fx,
	1 xcoll_atm, xl_atm
	
	REAL (kind=dp), SAVE :: cte1, preci_atm, to_max, to_min, rdx
	REAL (kind=dp) :: acc_c, corr, er, err, gravs2, w	

	INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: indpc
	INTEGER, SAVE :: l=1, ncoll_atm, compt_max, nc_atm, nl_atm 
	INTEGER :: compt, couche, cx, eq, i, id, ind, indice, j, ligne,
	1 spl, var, vare

	LOGICAL, SAVE :: init=.TRUE.
	LOGICAL :: inversible

c---------------------------------------------------------------------
	
2000	FORMAT(8es10.3)

c	PRINT*,'entrée dans coll_atm,ms_atm,n_atm',ms_atm,n_atm
c	PAUSE'coll_atm'

	IF(init)THEN
	 init=.FALSE.
	 cte1=g*msol/rsol**2
	 to_min=0.2d0 ; to_max=0.7d0		!bornes pour tau*
	 to_min=log(to_min) ; to_max=log(to_max)
	 
c	 création de la table des points de collocation
c	 ncoll_atm: nombre de points de collocation

	 ncoll_atm=(n_atm-1)*ms_atm
	 ALLOCATE(xcoll_atm(ncoll_atm))
	 CALL coll(x_atm,n_atm,ms_atm,xcoll_atm)
c	 PRINT*,'xcoll_atm' ; WRITE(*,2000)xcoll_atm
	 
c	 allocations diverses
c	 nc_atm : longueur du bloc du jacobien
c	 nl_atm : nombre de ligne du jacobein

	 nc_atm=ne_atm*ord_atm ; nl_atm=ne_atm*(ncoll_atm+r_qs)
	 
c	 PRINT*,ne_atm,r_qs,nl_atm,nc_atm,ord_atm ; PAUSE

c	 allocations et définitions diverses
	
	 ALLOCATE(ae(ne_atm,ne_atm,0:r_qs),be(ne_atm),a(nl_atm,nc_atm),
	1 b(1,nl_atm),indpc(nl_atm),y(ne_atm,0:r_qs),xl_atm(ne_atm),
	2 d(0:r_qs,ord_atm),ds(ord_atm,0:r_qs,ord_atm),
	3 dl(ne_atm,0:r_qs,ord_atm),fx(ne_atm),dfxdx(ne_atm))

c	 limites

	 xl_atm(1)=1	     !en tau=tau_max, q=1 R/Rsol=y(3)=Rb/Rsol
	 xl_atm(2)=n23_atm   !en tau=tau*, q=n23 ln T =y(2)=ln T(tau,Teff)
	 xl_atm(3)=n23_atm   !en tau=tau*, q=n23 ln T =y(2)=ln Teff
	 xl_atm(4)=n23_atm   !en tau=tau*, q=n23 ln tau =y(7)=ln tau*=y(6)
	 xl_atm(5)=n23_atm   !en tau=tau*, q=n23 R/Rsol=y(3)=y(4)=R23/Rsol
	 xl_atm(6)=n23_atm   !en tau=tau*, q=n23 M=y(5)=M*
	 xl_atm(7)=n_atm     !en tau=tau_min, q=N P = g / kappa tau_min
	 IF(pturb)xl_atm(Ipgt)=n_atm      !en tau=tau_min  ln Ptot=ln Pgaz

c	 les points de raccord x_atm sont les entiers 1, 2,..n_atm,
c	 entre deux points de raccord les abscisses des ord_atm
c	 points de collocation ne diffèrent que par un nombre entier 	
c	 les valeurs des ord_atm B-splines non id. nulles aux ms_atm
c	 points de collocation sur [i,i+1] peuvent être tabulées
c	 une fois pour toutes
	
	 DO cx=1,ms_atm
	  CALL linf(xcoll_atm(cx),xt_atm,knot_atm,l)
	  CALL bvald(xcoll_atm(cx),xt_atm,ord_atm,l,r_qs,d)
	  ds(cx,:,:)=d
c	  PRINT*,cx ; WRITE(*,2000)ds(cx,:,:)
	 ENDDO
	 
c	 idem pour les limites

	 DO i=1,ne_atm
	  CALL linf(xl_atm(i),xt_atm,knot_atm,l)
	  CALL bvald(xl_atm(i),xt_atm,ord_atm,l,r_qs-1,d) ; dl(i,:,:)=d
	 ENDDO

c	 d est désormais inutile
	 
	 DEALLOCATE(d) 
	 
	 preci_atm=1.d-7
c	 preci_atm=1.d-2
	 rdx=0.6d0	!facteur d'amortissement Newton-Raphson
c	 rdx=0.1d0	!facteur d'amortissement Newton-Raphson

	 compt_max=20		!nombre maximal d'itérations

	ENDIF
	
c l'acc. centrifuge ne doit pas excéder	90% de la gravité
	gravs2=cte1*mstar/rstar**2*0.9d0
	
c vitesse angulaire
	SELECT CASE(Krot)
	CASE(0,1,2)
	 w=wrot
	CASE(3)
	 w=rota(1,dim_rot)	 
	END SELECT
	 
c accélération centrifuge	 
	acc_c=rsol*rstar*w**2
	IF(acc_c > gravs2)THEN
	 SELECT CASE(langue)
	 CASE('english')
	   WRITE(*,1004)acc_c,gravs2,w ; WRITE(2,1004)acc_c,gravs2,w
1004	   FORMAT('STOP, in the atmophere, the centrifugal acceleration =',
	1    es10.3,' > 90% gravity = ',es10.3,', angular velicity =',
	2    es10.3)	   
	  CASE DEFAULT
	   WRITE(*,4)acc_c,gravs2,w ; WRITE(2,4)acc_c,gravs2,w
4	   FORMAT('ARRET, dans atmosphère, accélération centrifuge =',
	1    es10.3,' >  90% gravité = ',es10.3,', vitesse angulaire =',
	2    es10.3)	  
	  END SELECT
	  CALL sortie
	 ENDIF

c	DO i=1,n_atm
c	 CALL bsp1dn(ne_atm,bp_atm,x_atm,xt_atm,n_atm,ord_atm,knot_atm,
c	1 .TRUE.,x_atm(i),l,fx,dfxdx)
c	 WRITE(*,2000)EXP(fx(1:2)),fx(3:5),EXP(fx(6:7))
c	ENDDO
c	WRITE(*,2000)xchim(1:8)
c	WRITE(*,2000)rota(1,dim_rot)
c	PAUSE'entrée coll_atm'	
	
c	intégration, boucle infinie de NR

	compt=0
	B1: DO

c	on place la limite de la ZC externe au rayon bolométrique R*
c	R* = bp(4,i),tous i	
	 rstar=bp_atm(4,1)
	 r_zc_conv(if_conv)=rstar
	 
c	 pour les points de collocation

	 ligne=0 ; l=ord_atm ; a=0.d0 ; b=0.d0
	 DO i=1,ncoll_atm,ms_atm 	 
	  DO j=1,ms_atm
	   cx=i+j-1
	   CALL bsp1dn(ne_atm,bp_atm,x_atm,xt_atm,n_atm,ord_atm,knot_atm,
	1  .TRUE.,xcoll_atm(cx),l,fx,dfxdx)
	   y(:,0)=fx ; y(:,1)=dfxdx ; ind=ne_atm*(l-ord_atm)+1
c	   PRINT*,i,j,cx,l,ind ; WRITE(*,2000)xcoll_atm(cx)
c	   WRITE(*,2000)y(:,0) ; WRITE(*,2000)y(:,1)
c	   PAUSE'y dans coll_atm'

	   CALL eq_atm(1,1,xchim,cx,y,be,ae,r_rac,l_rac,xcoll_atm)
	   
	   DO eq=1,ne_atm
	    ligne=ligne+1 ; b(1,ligne)=be(eq) ; indpc(ligne)=ind
	    DO spl=1,ord_atm
	     DO var=1,ne_atm
	      indice=ne_atm*(spl-1)+var
	      DO id=0,r_qs
	       a(ligne,indice)=a(ligne,indice)+ae(eq,var,id)*ds(j,id,spl)
	      ENDDO	!id
	     ENDDO	!var
	    ENDDO	!spl
	   ENDDO	!eq	   	   
	  ENDDO		!j	  
	 ENDDO		!i
c	 PAUSE'avant lim'
	  
c	 pour les limites, comme r_qs=1, il y a ne_atm limites
	   
	 DO i=1,ne_atm
	  CALL bsp1dn(ne_atm,bp_atm,x_atm,xt_atm,n_atm,ord_atm,knot_atm,
	1 .TRUE.,xl_atm(i),l,fx,dfxdx)	
	  y(:,0)=fx ; y(:,1)=dfxdx
	  
	  CALL eq_atm(2,i,xchim,cx,y,be,ae,r_rac,l_rac,xcoll_atm)
	  
	  ligne=ligne+1 ; b(1,ligne)=be(1)
	  indpc(ligne)=ne_atm*(l-ord_atm)+1
	  DO spl=1,ord_atm
	   DO var=1,ne_atm
	    indice=ne_atm*(spl-1)+var
	    a(ligne,indice)=a(ligne,indice)+ae(1,var,0)*dl(i,0,spl)
	   ENDDO	!var
	  ENDDO		!spl
c	  PRINT*,ligne ; WRITE(*,2000)a(ligne,:),b(1,ligne) ; PAUSE'a'
	 ENDDO		!i
c	 PAUSE'apres lim'
	 
c	 PRINT*,ligne,nl_atm,nc_atm ; PRINT*,indpc
c	 WRITE(*,2000)b(1,:) ; PAUSE'b avant'

	 CALL gauss_band(a,b,indpc,nl_atm,nl_atm,nc_atm,1,inversible)
	 IF(.not.inversible)THEN
	  PRINT*,'jacobien non inversible dans coll_atm' ; STOP
	 ENDIF
c	 WRITE(*,2000)b(1,:) ; PAUSE'b après'	 
	 
c	 PAUSE'après gauss_band'	 
	 
c	 limitation des corrections

	 corr=1.d0
	 DO spl=1,dim_atm
	  DO var=1,ne_atm
	   indice=ne_atm*(spl-1)+var
c	   PRINT*,indice,corr
	   DO WHILE(abs(bp_atm(var,spl)) > preci_atm .and.
	1   corr*abs(b(1,indice)) > rdx*abs(bp_atm(var,spl)))
	    corr=corr/2.d0
	   ENDDO	!while
	  ENDDO	!var
	 ENDDO	!spl

c	 corrections et b --> bp

	 err=0.d0
	 DO spl=1,dim_atm
	  DO var=1,ne_atm
	   er=0.d0 ; indice=ne_atm*(spl-1)+var
	   IF(abs(bp_atm(var,spl)) > preci_atm)er=abs(b(1,indice)/
	1    bp_atm(var,spl))
	   bp_atm(var,spl)=bp_atm(var,spl)-b(1,indice)*corr
	   IF(var == 6)bp_atm(var,spl)=MAX(to_min,MIN(to_max,
	1    bp_atm(var,spl)))
	   err=max(er,err)
	   IF(er == err)THEN
	    vare=var ; couche=spl/ms_atm+1
	   ENDIF
	  ENDDO	!var
	 ENDDO	!spl
	 
c	 DO i=1,n_atm
c	  CALL bsp1dn(ne_atm,bp_atm,x_atm,xt_atm,n_atm,ord_atm,knot_atm,
c	1 .TRUE.,x_atm(i),l,fx,dfxdx)
c	  WRITE(*,2000)EXP(fx(1:2)),fx(3:5),EXP(fx(6:7))
c	 ENDDO
c	 PAUSE'solution colatm1'

	 compt=compt+1
	 WRITE(*,1)compt,err,nom_atm(vare),couche,corr
1	 FORMAT('atmosphère iter.',i3,' err. max.',es8.1,
	1 ', var:',a,', couche',i3,', corr',es8.1)
	 
c	 test de sortie

	 IF(err < preci_atm)THEN
	  EXIT B1
	 ELSEIF(compt < compt_max)THEN
	  CYCLE B1
	 ELSEIF(err < 1.d-2)THEN
	  WRITE(*,2) ; WRITE(2,2)
2	  FORMAT('on force la convergence pour l''atmosphère')
	  EXIT B1	 
	 ELSE
	  WRITE(*,3)compt ; WRITE(2,3)compt
3	  FORMAT('pas de conv. pour l''atmosphère après ',i3,' itérations',
	1 /,'on force la convergence')
	  EXIT B1
	 ENDIF
	ENDDO B1
	
	RETURN

	END SUBROUTINE coll_atm
