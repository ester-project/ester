	
c************************************************************************

         SUBROUTINE pertw_ptm(r,mw_dot)

c calcul de la perte de moment cinétique due à la perte de masse
c ne concerne que la ZC ext

c routine private du module mod_evol

c Auteur: P.Morel, Département Cassiopée, O.C.A., CESAM2k

c----------------------------------------------------------------

	USE mod_donnees, ONLY : langue, mdot, p_pertw, rsol, secon6
	USE mod_kind
	USE mod_nuc, ONLY : mzc_ext
	USE mod_variables, ONLY : dim_rot, mstar, rota
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: r	
	REAL (kind=dp), INTENT(out) :: mw_dot
	
	REAL (kind=dp),	SAVE :: coef
		
	LOGICAL, SAVE :: init=.TRUE.	
		
c---------------------------------------------------------------- 

	IF(init)THEN
	 init=.FALSE.
	 
c il ne peut y avoir que des pertes de moment cinétique
	 lw_perte = mdot <= 0.d0
	 IF(.NOT.lw_perte)RETURN
	 	 
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1000)p_pertw ; WRITE(2,1000)p_pertw
1000	  FORMAT('The ang. momentum change due to the mass loss',/,
	1 'loss of ang. momentum=',es10.3,', in the external ZC.')
	 CASE DEFAULT	 
	  WRITE(*,1)p_pertw ; WRITE(2,1)p_pertw
1	  FORMAT('Variation de moment cinétique dû à la perte de masse',/,
	1 'perte de moment cinétique :',es10.3,', dans la ZC externe')
	 END SELECT
	 
c 1.d6/secon6 nombre de sec/an	 
	 coef=mdot*1.d6/secon6*rsol**2
	ENDIF

c perte de moment cinétique Msol/Myr
	mw_dot=coef/mzc_ext*r**2*rota(1,dim_rot)
	
	RETURN
	
	END SUBROUTINE pertw_ptm
