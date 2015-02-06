	
c************************************************************************

         SUBROUTINE pertw_ptm(nu,r,mw_dot)

c	routine du calcul de la perte/gain de moment cinétique lié à la
c	perte de masse

c	routine private du module mod_evol

c	Auteur: P.Morel, Département Cassiopée, O.C.A.
c	CESAM2k

c----------------------------------------------------------------

	USE mod_donnees, ONLY : langue, mdot, pw_extend, p_pertw,
	1 rsol, secon6
	USE mod_kind
	USE mod_variables, ONLY : mstar
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: nu, r	
	REAL (kind=dp), INTENT(out) :: mw_dot
	
	REAL (kind=dp),	SAVE :: coef
		
	LOGICAL, SAVE :: init=.TRUE.	
		
c---------------------------------------------------------------- 

c	il ne peut y avoir que des apports de moment cinétique

	IF(mdot <= 0.d0)THEN
	 mw_dot=0.d0 ; RETURN
	ENDIF

	IF(init)THEN
	 init=.FALSE.	 
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1000)p_pertw,pw_extend ; WRITE(2,1000)
1000	  FORMAT('The ang. momentum change due to the mass loss',/,
	1   'input of ang. momentum=',es10.3,', in m/Mtot :',es10.3)
	 CASE DEFAULT	 
	  WRITE(*,1)p_pertw,pw_extend ; WRITE(2,1)p_pertw,pw_extend
1	  FORMAT('Variation de moment cinétique dû à la perte de masse',/,
	1   'input de moment cinétique :',es10.3,', répartit sur m/Mtot :',
	2   es10.3)
	 END SELECT
	 mlim_w=(mstar*(1.d0-pw_extend))**(2.d0/3.d0)
	 coef=p_pertw*mdot/secon6/mlim_w
	ENDIF

c	en deça de mlim_w il n'y a pas d'apport de moment cinétique
	IF(nu < mlim_w)THEN
	 mw_dot=0.d0
	ELSE
	 mw_dot=coef*(r*rsol)**2
	ENDIF
	
	RETURN
	
	END SUBROUTINE pertw_ptm
