
c**********************************************************************

      SUBROUTINE pertm_waldron(dt)

c	routine private du module mod_static

c perte de masse empirique de Waldron
c cite dans l'article A.Ap 229 (1990) 469-474

c	routine d'interpolation m(t+dt)--->m(t) en tenant compte 	
c	de la perte de masse (mdot>0 : gain de masse,
c	mdot < 0 : perte de masse)
c	la perte de masse est concentrée dans la couche n_ptm-1 n_ptm

c	utilisation par sbsp1dn (m**2/3 ---> m**2/3 ancien)

c	Auteur: B. Pichon, Département J.D. Cassini, O.C.A.,

c	25 08 97 : mise en place des variables eulériennes

c	en_masse = .true. variables lagrangiennes m23=m**23, r2=r**2
c	en_masse = .false. variables eulériennes m23=m, r2=r

c	CESAM2k

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

	USE mod_donnees, ONLY : en_masse, mdot, m_ptm, ne, ord_qs
	USE mod_kind
	USE mod_numerique, ONLY : bsp1dn, no_croiss
	USE mod_variables, ONLY : bp, knot, knot_ptm, mstar, mstar_t,
	1 n_ptm, n_qs, old_ptm, q, qt, rstar, xt_ptm, x_ptm
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: dt
	REAL (kind=dp), ALLOCATABLE, SAVE, DIMENSION(:) :: df, f
	REAL (kind=dp), ALLOCATABLE, DIMENSION(:) :: tmp1, tmp2
	
	REAL (kind=dp) :: log_mdot, luminosite	
		
	INTEGER :: i, l
	
	LOGICAL, SAVE :: init=.TRUE.

c--------------------------------------------------------------------------
	
2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE. ; ALLOCATE(df(ne),f(ne))
	ENDIF
	 	
c extraction des masses au temps t+dt
c	m^2/3 en lagrangien, m en eulérien on a n_ptm=n_qs

	ALLOCATE(tmp1(n_qs),tmp2(n_qs))	
	DO i=1,n_qs
	 CALL bsp1dn(ne,bp,q,qt,n_qs,ord_qs,knot,.TRUE.,q(i),l,f,df)
	 tmp1(i)=f(5) 		!; WRITE(*,2000)f	 
	ENDDO
	
	n_ptm=1 ; tmp2(n_ptm)=tmp1(1)
	DO i=2,n_qs
	 IF(tmp1(i) > tmp2(n_ptm))THEN
	  n_ptm=n_ptm+1 ; tmp2(n_ptm)=tmp1(i)
	 ENDIF
	ENDDO
	  
c	allocations

	IF(ALLOCATED(old_ptm))THEN
	 DEALLOCATE(old_ptm,x_ptm,xt_ptm)
	ENDIF
	ALLOCATE(old_ptm(1,n_ptm),x_ptm(n_ptm),xt_ptm(n_ptm+m_ptm))
	x_ptm=tmp2(1:n_ptm) ; old_ptm(1,1:n_ptm)=tmp2(1:n_ptm)
	
	DEALLOCATE(tmp1,tmp2)
	
c	calcul de mstar(t+dt)
C     d'abord calcul de log_mdot (rstar est en Rsun)
C      puis de mdot (en Msun/an)

	IF(en_masse)THEN
	 luminosite=SQRT(bp(4,n_qs))**3
	ELSE
	 luminosite=bp(4,n_qs)
	ENDIF
	log_mdot=1.07d0*LOG10(luminosite)+1.77d0*LOG10(rstar)-14.3d0
	mdot = 10.d0**log_mdot

	mstar=mstar_t+mdot*1.d6*dt
	
c	tabulation

	CALL bsp1dn(1,old_ptm,x_ptm,xt_ptm,n_ptm,m_ptm,knot_ptm,.FALSE.,
	1 x_ptm(1),l,f,df)
        IF(no_croiss)THEN
         PRINT*,'Arrêt 1 dans pertm_waldron' ; STOP
        ENDIF
	
c	WRITE(*,2000)mtot-mstar ; PRINT*,mtot,mstar ; PAUSE

	RETURN

	END SUBROUTINE pertm_waldron
