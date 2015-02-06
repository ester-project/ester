
c***********************************************************************

	SUBROUTINE diffus(ok,dt,convd,convf,nzc,mc_tmp,nc_tmp)
	
c	routine PRIVATE du module mod_evol	
	
c	gestion de la résolution des systèmes des équations de diffusion
c	du moment cinétique et des éléments chimiques
c	méthode des éléments finis Petrov-Galerkin

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c entrées

c	dt: pas temporel
c	convd,convf,nzc: limites des ZC, nombre de ZM
c	mc_tmp,nc_tmp: abscisses temporaires pour l'intégration, et nombre 

c sorties
c	ok=.TRUE.: il y a eu convergence

c NOTATIONS (hélas incohérentes)
c	n_ch : nombre VARIABLE de points, élément de mod_variables
c	nch : nombre FIXE de fonctions, élément de mod_donnees
c	m_ch : ordre FIXE des splines, élément de mod_donnees 
c	mch(n_ch) : abscisses VARIABLES, élément de mod_variables

c----------------------------------------------------------------------

	USE mod_donnees, ONLY : Krot, langue, mdot, m_ch, m_rot,
	1 nchim, nrot, rsol, secon6
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn, bsp1ddn, bsp_dis, genere_bases,
	1 no_croiss, newspl, noeud
	USE mod_variables, ONLY : chim, dim_ch, dim_rot,
	1 knotc, knotr, mc, mct, mrot, mrott, n_ch, n_qs, n_rot,
	2 rota, rstar, sortie, tot_conv, tot_rad
       
	IMPLICIT NONE

	REAL (kind=dp),  INTENT(in), DIMENSION(:) :: mc_tmp
	REAL (kind=dp), INTENT(in) :: dt
	INTEGER, INTENT(in), DIMENSION(0:) :: convf
	INTEGER, INTENT(in), DIMENSION(:) :: convd		
	INTEGER, INTENT(in) :: nc_tmp, nzc
	LOGICAL, INTENT(out) :: ok
   
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: chim_tmp, fmixd,
	1 rota_tmp
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: mct_tmp, mrott1	
	REAL (kind=dp), DIMENSION(nmix) :: dgmix, gmix
	REAL (kind=dp), DIMENSION(nrot) :: dfrot, frot
	REAL (kind=dp) :: nu
	INTEGER, ALLOCATABLE, DIMENSION(:) :: idis, mdis
	INTEGER, SAVE :: l=1	
	INTEGER :: i, izc, j, knotc_tmp, ndis
	INTEGER :: knotr1
	
	LOGICAL, ALLOCATABLE, DIMENSION(:) :: ldis      
	LOGICAL, SAVE :: init=.TRUE.
	
c-------------------------------------------------------------------------

2000	FORMAT(8es10.3)

c	PRINT*,tot_rad,nzc ; PRINT*,'convd',convd(1:nzc+1)
c	PRINT*,'convf',convf(0:nzc)
	
c définitions et initialisation 

	IF(init)THEN
         init=.FALSE.

	 SELECT CASE(langue)
	 CASE('english')
          WRITE(*,1001) ; WRITE(2,1001)
1001      FORMAT(/,'Diffusion of chemicals using the finite element method')	 
          WRITE(*,1004)m_ch ; WRITE(2,1004)m_ch     
1004      FORMAT('order of B-splines:',i3)
	 CASE DEFAULT 
          WRITE(*,1) ; WRITE(2,1)
1         FORMAT(/,'Diffusion des éléments chimiques par éléments finis')
          WRITE(*,4)m_ch ; WRITE(2,4)m_ch     
4         FORMAT('ordre des B-splines:',i3)
	 END SELECT

         IF(mdot /= 0.d0)THEN
	  SELECT CASE(langue)
	  CASE('english')
           WRITE(*,1005)mdot ; WRITE(2,1005)mdot
1005       FORMAT('rate of external mass change = ',es10.3,' Msun/year')
	  CASE DEFAULT	  
           WRITE(*,5)mdot ; WRITE(2,5)mdot
5          FORMAT('taux de variation de masse = ',es10.3,' Msol/an')
	  END SELECT		
         ELSE
	  SELECT CASE(langue)
	  CASE('english')
           WRITE(*,1006) ; WRITE(2,1006)
1006       FORMAT('without mass change')	  
	  CASE DEFAULT
           WRITE(*,6) ; WRITE(2,6)
6          FORMAT('sans variation externe de masse')
	  END SELECT
         ENDIF
	 
	 IF(Krot == 3)THEN
	  ALLOCATE(xlim_rot(nlim_rot))
	  SELECT CASE(langue)
	  CASE('english')
           WRITE(*,1007) ; WRITE(2,1007)
1007       FORMAT(/,'Diffusion of ang. mom. with the finite element method')	 
           WRITE(*,1008)m_rot ; WRITE(2,1008)m_rot     
1008       FORMAT('order of B-splines:',i3)
	  CASE DEFAULT 
           WRITE(*,7) ; WRITE(2,7)
7          FORMAT(/,'Diffusion du moment cinétique par éléments finis')
           WRITE(*,8)m_rot ; WRITE(2,4)m_rot     
8          FORMAT('ordre des B-splines:',i3)
	  END SELECT  
	 ENDIF
	 
	ENDIF       !initialisation

c	le vecteur de mélange 1 dans ZR, -1 dans ZC	
c	x_mix, x_mixt, fmix, n_mix, knot_mix, m_mix

	IF(ALLOCATED(fmix))DEALLOCATE(fmix,x_mix,x_mixt)	
	IF(tot_rad .OR. tot_conv)THEN
	 n_mix=nc_tmp
         ALLOCATE(fmix(nmix,n_mix),idis(0),ldis(0),x_mix(n_mix),
	1 x_mixt(n_mix+m_mix))
         x_mix=mc_tmp(1:n_mix) ; ndis=0
         IF(tot_rad)THEN
          fmix=1.d0
         ELSE
          fmix=-1.d0
         ENDIF	
         CALL bsp1dn(nmix,fmix,x_mix,x_mixt,n_mix,m_mix,knot_mix,.FALSE.,
	1  x_mix(1),l,gmix,dgmix)

	ELSE

c        il y a des ZR / ZC, 
c	 conventions pour les indices des délimitations
c	 convf(0)=1 si convd(1) /= 1, sinon convf(0)=10000
c	 convd(nzc+1)=nc_tmp si convf(nzc) /= nc_tmp, sinon
c	 convd(nzc+1)=-10000
c	 modèle tot.conv. nzc=1, convd(1)=1, convf(1)=n_ch,
c	 lconv(1)=.FALSE.
c	 modèle tot.rad.  nzc=0, convd(1)=n_ch, convf(0)=1
c	 PRINT*,nzc ; PRINT*,'convd',convd ; PRINT*,'convf',convf
c	 CALL pause('nzc')
c	la ZC d'indice i=1,nzc commence en convd(i) et finit en convf(i)-eps
c	la ZR d'indice i=1,nzc commence en convf(i) et finit en convd(i)-eps

c estimation de ndis : nombre de discontinuités
	 ndis=2*nzc ; n_mix=nc_tmp
	 IF(convd(1) == 1)ndis=ndis-1
	 IF(convf(nzc) == n_mix)ndis=ndis-1	 	
	 knot_mix=n_mix+ndis+m_mix	 	  
	 ALLOCATE(fmix(nmix,n_mix+ndis),fmixd(nmix,ndis),idis(0:ndis+1),
	1  ldis(ndis),x_mix(n_mix),x_mixt(knot_mix)) 
	 x_mix=mc_tmp(1:n_mix) ; idis=-100 ; idis(0)=1
	 idis(ndis+1)=n_mix ; j=0
	 DO izc=1,nzc
	  IF(convd(izc) /= 1)THEN
	   j=j+1 ; idis(j)=convd(izc) ; ldis(j)=.TRUE.
	  ENDIF
	  IF(convf(izc) /= n_mix)THEN 
	   j=j+1 ; idis(j)=convf(izc) ; ldis(j)=.FALSE.
	  ENDIF
	 ENDDO
	 fmix(1,:)=1.d0 ; fmixd(1,:)=1.d0
	 DO izc=1,nzc
	  DO j=convd(izc)+1,convf(izc)
	   fmix(1,j)=-1.d0
	  ENDDO
	 ENDDO
	 IF(convd(1) == 1)fmix(1,1)=-1.d0
	 DO j=1,ndis
	  IF(ldis(j))fmixd(1,j)=-1.d0
	 ENDDO
         CALL bsp_dis(nmix,x_mix,fmix,ndis,idis,fmixd,eps,n_mix,m_mix,
	1  x_mixt,knot_mix)
         IF(no_croiss)THEN
          PRINT*,'Arrêt 2 dans diffus' ; CALL sortie
         ENDIF 
	 DEALLOCATE(fmixd)
	ENDIF	
	
c on s'assure que la zone externe est convective
c	nu=(x_mix(n_mix)+x_mix(n_mix-1))/2.d0
c	CALL bsp1dn(nmix,fmix,x_mix,x_mixt,n_mix,m_mix,knot_mix,.TRUE.,
c	1 nu,l,frot,dfrot)	
c	IF(frot(1) > 0.d0)THEN
c	 SELECT CASE(langue)
c	 CASE('english')
c	  WRITE(*,1002) ; WRITE(2,1002)
c1002	  FORMAT('ARRET, the external part is radiative')	 
c	 CASE DEFAULT
c	  WRITE(*,2) ; WRITE(2,2)
c2	  FORMAT('ARRET, la partie externe est radiative')
c	 END SELECT
	 
c	 PRINT*,idis ; PRINT*,ndis,nzc ; PRINT*,convd(1:nzc) ; PRINT*,convf(1:nzc)
c	 PRINT*,x_mix(convd(1:nzc)) ; PRINT*,x_mix(convf(1:nzc))
	 
c	 j=1
c	 DO i=1,n_mix
c	  nu=x_mix(i)
c	  CALL bsp1dn(nmix,fmix,x_mix,x_mixt,n_mix,m_mix,knot_mix,.TRUE.,
c	1   nu,l,frot,dfrot)	
c	  WRITE(*,2000)nu,frot(1)
c	  IF(i == idis(j))THEN
c	   PRINT*,'limite ZR/ZC j=',j,idis(j),i	   	  
c	   nu=un_eps*x_mix(i)+eps*x_mix(i+1)	   
c	   PRINT*,x_mix(i-1),x_mix(i),nu,x_mix(i+1)	   
c	   j=MIN(j+1,ndis)
c	   CALL bsp1dn(nmix,fmix,x_mix,x_mixt,n_mix,m_mix,knot_mix,.TRUE.,
c	1    nu,l,frot,dfrot)	
c	   WRITE(*,2000)nu,frot(1)
c	  ENDIF
c	 ENDDO	 
c	 CALL sortie
c	ENDIF	

c test	 
c	IF(.TRUE.)THEN
	IF(.FALSE.)THEN
	 PRINT*,'test limites'	
	 j=1
	 DO i=1,n_mix
	  nu=x_mix(i)
	  CALL bsp1dn(nmix,fmix,x_mix,x_mixt,n_mix,m_mix,knot_mix,.TRUE.,
	1   nu,l,frot,dfrot)	
	  WRITE(*,2000)nu,frot(1)
	  IF(i == idis(j))THEN
	   PRINT*,'limite ZR/ZC j=',j,idis(j),i	   	  
	   nu=un_eps*x_mix(i)+eps*x_mix(i+1)
	   
	   PRINT*,x_mix(i-1),x_mix(i),nu,x_mix(i+1)
	   
	   j=MIN(j+1,ndis)
	   CALL bsp1dn(nmix,fmix,x_mix,x_mixt,n_mix,m_mix,knot_mix,.TRUE.,
	1    nu,l,frot,dfrot)	
	   WRITE(*,2000)nu,frot(1)
	  ENDIF
	 ENDDO
	 PAUSE'test limites'
	ENDIF
	
c pour l'initialisation des itérations NR, on utilise la comp.chim.
c	se trouvant dans chim. Au 1-ier appel à evol (ie. pour
c	un nouveau pas temporel) c'est la comp.chim. de chim_t=chim sur
c	mc_tmp qui sert pour l'initialisation primaire. Au cours des
c	itérations NR on utilise celle obtenue à l'itération NR précédente.
c	Le nombre, la position des couches, les limites ZR/ZC,
c	la masse totale changent d'une itération NR
c	quasi-statique à la suivante. Pour le calcul par éléménts finis
c	de chim on utilise une base avec discontinuité de
c	la dérivée 1-ière aux limites ZR/ZC.

c mc --> mc_tmp	
	ALLOCATE(chim_tmp(nchim,SIZE(chim,2)),mct_tmp(SIZE(mct)))	
	chim_tmp=chim ; mct_tmp=mct ; knotc_tmp=knotc ; n_ch=nc_tmp

	DEALLOCATE(mc,mct) ; ALLOCATE(mc(n_ch))
	mc=mc_tmp(1:n_ch) ; mc(1)=0.d0
	 
c génération de la base (dérivée 1-ère discontinue)
	ALLOCATE(mdis(ndis)) ; mdis=m_ch-1
	CALL genere_bases(idis,m_ch,mdis,n_ch,ndis,mc,knotc,
	1 mct)

c bricolage sans importance, il s'agit d'initialisations,
c on modifie la position des limites au temps t 
c si le rayon externe a augmenté 
	IF(mct_tmp(knotc_tmp) < mc(n_ch))THEN
	 mct_tmp(knotc_tmp-m_ch+1:knotc_tmp)=mc(n_ch)
	ENDIF
	
c spline sur le nouveau vecteur nodal mct	
	dim_ch=knotc-m_ch
	DEALLOCATE(chim) ; ALLOCATE(chim(nchim,dim_ch))	 	 
	CALL newspl(nchim,mc,mct_tmp,knotc_tmp,m_ch,mct,knotc,m_ch,
	1 chim_tmp,chim)
	IF(no_croiss)PRINT*,'Pb. en 1 dans diffus'	
	DEALLOCATE(chim_tmp,mct_tmp)

c pour le moment cinétique on utilise une base différente de mct
c rota(n_rota,mrot,mrott,m_rot,knotr) pour avoit la liberté d'utiliser un ordre
c de spline différent m_rot différent de m_ch, toutefois mrot=mc
	IF(Krot == 3)THEN
	 ALLOCATE(mrott1(knotr),rota_tmp(nrot,SIZE(rota,2)))
	 rota_tmp=rota ; mrott1=mrott ; knotr1=knotr ; n_rot=n_ch

	 DEALLOCATE(mrot,mrott) ; ALLOCATE(mrot(n_rot))
	 mrot=mc ; nu_min=eps*mrot(2)

c génération de la base
	 mdis=m_rot-1
	 CALL genere_bases(idis,m_rot,mdis,n_rot,ndis,mrot,knotr,
	1  mrott)
		
c spline sur le nouveau vecteur nodal mrott	
	 dim_rot=knotr-m_rot
	 DEALLOCATE(rota) ; ALLOCATE(rota(nrot,dim_rot))	 	 

c bricolage sans importance, il s'agit d'initialisations,
c on modifie la position des limites au temps t 
c si le rayon externe a augmenté 
	 IF(mrott1(knotr1) < mrot(n_rot))THEN
	  mrott1(knotr1-m_rot+1:knotr1)=mrot(n_rot)
	 ENDIF
	 
c idem au centre  
	 mrott1(1:m_rot)=mrott(1:m_rot)

c nouvelle spline pour initialisation	 
	 CALL newspl(nrot,mrot,mrott1,knotr1,m_rot,mrott,knotr,m_rot,
	1  rota_tmp,rota)
	 IF(no_croiss)PRINT*,'Pb. en 2 dans diffus'	
	 DEALLOCATE(mrott1)

c écriture
c	IF(.TRUE.)THEN
	IF(.FALSE.)THEN
	 DEALLOCATE(rota_tmp) ; ALLOCATE(rota_tmp(5,0:1))
	 DO i=1,n_rot
	  CALL bsp1ddn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,
	1   mrot(i),l,rota_tmp)
	  CALL bsp1dn(nmix,fmix,x_mix,x_mixt,n_mix,m_mix,knot_mix,.TRUE.,
	1   mrot(i),l,frot,dfrot)
c	  PRINT*,i
	  WRITE(*,2000)rota_tmp(1:5,0),frot(1),REAL(i,dp)
c	  WRITE(*,2000)rota_tmp(1:5,1)
	  IF(MOD(i,300) == 0)PAUSE'rota'
	 ENDDO
	 PAUSE'fin rota initialisation'
	ENDIF
	DEALLOCATE(mdis)
	
c	nettoyage et abscisses des limites	
	 DEALLOCATE(rota_tmp) 	 
	 xlim_rot(1)=mrot(1) ; xlim_rot(2)=mrot(n_rot) 	 	 
	ENDIF

c	 diffusion du moment cinétique		 
	IF(Krot == 3)THEN
	 PRINT*
	 SELECT CASE(langue)	 
	 CASE('english')
	  PRINT*,'---- diffusion of angular momentun (begin) ----'
	 CASE DEFAULT
	  PRINT*,'---- diffusion du moment cinétique (début) ----'
	 END SELECT

c	tabulation des coefficient de diffusion
	 CALL coeff_rota(dt,idis,ndis)

c résolution des équations  
	 CALL resout_rota(dt,ok) ; IF(.NOT.ok)RETURN

c écriture

c	IF(.TRUE.)THEN
	IF(.FALSE.)THEN	
	 ALLOCATE(rota_tmp(nrot,0:m_rot-1))
	 DO i=1,n_rot
	  CALL bsp1ddn(nrot,rota,mrot,mrott,n_rot,m_rot,knotr,.TRUE.,
	1   mrot(i),l,rota_tmp)
c	    WRITE(*,2000)rota_tmp(1:2,0),rota_tmp(3,0:1),rota_tmp(5,0)
	  WRITE(*,2000)rota_tmp(:,0)
	 ENDDO
	 DEALLOCATE(rota_tmp) 
	 PAUSE'solution'
	ENDIF
  
	 SELECT CASE(langue)	 
	 CASE('english')
	  PRINT*,'---- diffusion of angular momentun (end) ----'
	 CASE DEFAULT
	  PRINT*,'---- diffusion du moment cinétique (fin) ----'
	 END SELECT  
	ENDIF
	
c diffusion de la composition chimique
	PRINT*
	SELECT CASE(langue)
	CASE('english')
	 PRINT*,'---- diffusion of chemical species (begin) ----'
	CASE DEFAULT
	 PRINT*,'---- diffusion des éléments chimiques (début) ----'
	END SELECT
			
	CALL resout_chim(dt,ok) ; IF(.NOT.ok)RETURN

	SELECT CASE(langue)
	CASE('english')
	 PRINT*,'---- diffusion of chemical species (end) ----'
	CASE DEFAULT
	 PRINT*,'---- diffusion des éléments chimiques (fin) ----'
	END SELECT
	
c nettoyage
	DEALLOCATE(fmix,idis,x_mix,x_mixt)
		
	RETURN

	END SUBROUTINE diffus
