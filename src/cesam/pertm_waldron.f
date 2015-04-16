
c**********************************************************************

	SUBROUTINE pertm_waldron(dt)

c routine private du module mod_static
c perte de masse empirique de Waldron
c cite dans l'article A.Ap 229 (1990) 469-474

c routine d'interpolation m(t+dt)--->m(t) en tenant compte 	
c de la perte de masse (mdot>0 : gain de masse,
c mdot < 0 : perte de masse)
c la perte de masse est concentrée dans la couche n_ptm-1 n_ptm
c on assure la stricte croissance des masses et on calcule la nouvelle
c masse totale


c Auteur: B. Pichon, Département J.D. Cassini, O.C.A., CESAM2k

c entrées
c	bp,q,n,qt,knot,chim,mc,nc,mct,knotc : modèle au temps t
c	dt : pas temporel
c	mstar_t: mstar au temps age

c entrée/sortie
c	mstar: mstar au temps age+dt

c sortie
c	old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm : interpo. de
c	l'ancienne masse
c	en fonction de la nouvelle (en m**2/3) normalise (Mstar x Msol)
c	old_ptm et x_ptm sont identiques ce n'est pas le cas si
c	on tient compte de la perte de masse due a E=mc**2

c----------------------------------------------------------------

	USE mod_donnees, ONLY : en_masse, mdot, m_ptm, nchim, ne, ord_qs
	USE mod_kind
	USE mod_nuc, ONLY : l_planet, planetoides
	USE mod_numerique, ONLY : bsp1dn, no_croiss
	USE mod_variables, ONLY : bp, knot, knot_ptm, mstar, mstar_t,
	1 n_ptm, n_qs, old_ptm, q, qt, rstar, xt_ptm, x_ptm
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: dt

	REAL (kind=dp), DIMENSION(ne) :: df, f
	REAL (kind=dp), DIMENSION(n_qs) :: tmp1, tmp2	
	
	REAL (kind=dp) :: dm, log_mdot, luminosite, m_planet	
		
	INTEGER :: i, l

c--------------------------------------------------------------------------
	
2000	FORMAT(8es10.3)
	 	
c extraction des masses au temps t+dt
c m^2/3 en lagrangien, m en eulérien on a n_ptm=n_qs
	DO i=1,n_qs
	 CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),l,f,df)
	 tmp1(i)=f(5) 		!; WRITE(*,2000)f
	 IF(i == n_qs)THEN
	  IF(en_masse)THEN
	   luminosite=SQRT(f(4))**3
	  ELSE
	   luminosite=f(4)
	  ENDIF
	 ENDIF	 
	ENDDO
	
c on s'assure de la stricte croissance
	n_ptm=1 ; tmp2(n_ptm)=tmp1(1)
	DO i=2,n_qs
	 IF(tmp1(i) > tmp2(n_ptm))THEN
	  n_ptm=n_ptm+1 ; tmp2(n_ptm)=tmp1(i)
	 ENDIF
	ENDDO
	  
c allocations
	IF(ALLOCATED(old_ptm))DEALLOCATE(old_ptm,x_ptm,xt_ptm)
	ALLOCATE(old_ptm(1,n_ptm),x_ptm(n_ptm),xt_ptm(n_ptm+m_ptm))
	x_ptm=tmp2(1:n_ptm) ; old_ptm(1,1:n_ptm)=tmp2(1:n_ptm)

c tabulation de l'ancienne répartition en fonction de la nouvelle
	CALL bsp1dn(1,old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm,.FALSE.,
	1 x_ptm(1),l,f,df)
	
c calcul de log_mdot (rstar est en Rsun) puis de mdot (en Msun/an)
	log_mdot=1.07d0*LOG10(luminosite)+1.77d0*LOG10(rstar)-14.3d0
	mdot = 10.d0**log_mdot
	
c dm : perte/gain de masse, selon signe de mdot en Msol/an
	dm=mdot*1.d6*dt
	
c contribution des planétoïdes, m_planet en Msol/My
	IF(l_planet)THEN
	 CALL planetoides(m_planet=m_planet)
	 dm=dm+m_planet*dt
	ENDIF
	
c calcul de mstar(t+dt), en entrée mstar_t(Msol) est la masse au temps t
c mstar contient l'atmosphère alors que m(n_ptm) ne la contient pas	
	mstar=mstar_t+dm	

	RETURN

	END SUBROUTINE pertm_waldron
