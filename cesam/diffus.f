
c***********************************************************************

	SUBROUTINE diffus(ok,dt,mc_tmp,nc_tmp)
	
croutine PRIVATE du module mod_evol	
	
c gestion de la résolution des systèmes des équations de diffusion
c du moment cinétique et des éléments chimiques
c méthode des éléments finis Petrov-Galerkin

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

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

	USE mod_donnees, ONLY : langue, mdot, m_ch, nchim, rsol
	USE mod_kind
	USE mod_numerique, ONLY : bsp1ddn, no_croiss, noeud, newspl_gal
	USE mod_variables, ONLY : chim, dim_ch,
	1 knotc, mc, mct, mstar, n_ch, rstar, sortie
       
	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(:) :: mc_tmp
	REAL (kind=dp), INTENT(in) :: dt
	INTEGER, INTENT(in) :: nc_tmp
	LOGICAL, INTENT(out) :: ok
   
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:,:) :: chim_tmp
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: mct_tmp

	INTEGER :: knotc_tmp
	
	LOGICAL, SAVE :: init=.TRUE.
	
c-------------------------------------------------------------------------

2000	FORMAT(8es10.3)

c	PRINT*,nzc ; PRINT*,'convd',convd(1:nzc+1)
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
	ENDIF       !initialisation
	
c pour l'initialisation des itérations NR, on utilise la comp.chim.
c	se trouvant dans chim. Au 1-ier appel à evol (ie. pour
c	un nouveau pas temporel) c'est la comp.chim. de chim_t=chim sur
c	mc_tmp qui sert pour l'initialisation primaire. Au cours des
c	itérations NR on utilise celle obtenue à l'itération NR précédente.
c	Le nombre, la position des couches, les limites ZR/ZC,
c	la masse totale changent d'une itération NR
c	quasi-statique à la suivante. Pour le calcul par éléments finis
c	de chim on utilise :
c	  --une base continue avec discontinuité de la dérivée 1-ière aux
c	    limites ZR/ZC, avantages : plus propre moins de calculs,
c	    inconvénient localisation des limites.
c	  --une base avec discontinuité de la dérivée 1-ière en tout point de
c	    raccord, avantage : pas de localisation des mimites, inconvénient
c	    le double de calculs

c avec rotation on fixe la grille à celle de la 1-ère itération
c sans rotation on adapte la grille

c mc_tmp --> mc	
	ALLOCATE(chim_tmp(nchim,SIZE(chim,2)),mct_tmp(SIZE(mct)))	
	chim_tmp=chim ; mct_tmp=mct ; knotc_tmp=knotc ; n_ch=nc_tmp
	
	DEALLOCATE(mc) ; ALLOCATE(mc(n_ch))
	mc=mc_tmp(1:n_ch) ; mc(1)=0.d0
	 
c génération de la base pour la composition chimique			
c base avec dérivée 1ère discontinue aux limites ZR/ZC
	CALL base_chim
	
c spline sur le nouveau vecteur nodal mct	
	dim_ch=knotc-m_ch
	DEALLOCATE(chim) ; ALLOCATE(chim(nchim,dim_ch))
c	CALL newspl(nchim,mc,mct_tmp,knotc_tmp,m_ch,mct,knotc,m_ch,
c	1 chim_tmp,chim)
	CALL newspl_gal(nchim,mc,mct_tmp,knotc_tmp,m_ch,chim_tmp,
	1 mct,knotc,m_ch,chim)
		
	IF(no_croiss)PRINT*,'Pb. en 1 dans diffus'	
	DEALLOCATE(chim_tmp,mct_tmp)
	
c diffusion de la composition chimique
	PRINT*
	SELECT CASE(langue)
	CASE('english')
	 WRITE(usl_evol,*)'---- diffusion of chemical species (begin) ----'
	CASE DEFAULT
	 WRITE(usl_evol,*)'---- diffusion des éléments chimiques (début) ----'
	END SELECT
			
	CALL resout_chim(dt,ok) ; IF(.NOT.ok)RETURN
	
c	PRINT*,chim(1,1:20) ; PAUSE'diffus1'
c	chim=REAL(chim,kind=sp)
c	PRINT*,chim(1,1:20) ; PAUSE'diffus2'

	SELECT CASE(langue)
	CASE('english')
	 WRITE(usl_evol,*)'---- diffusion of chemical species (end) ----'
	CASE DEFAULT
	 WRITE(usl_evol,*)'---- diffusion des éléments chimiques (fin) ----'
	END SELECT

c	PAUSE'diffus solution'
			
	RETURN

	CONTAINS
	 INCLUDE 'base_chim.f'

	END SUBROUTINE diffus
