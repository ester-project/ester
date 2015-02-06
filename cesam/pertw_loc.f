	
c************************************************************************

         SUBROUTINE pertw_loc(nu,omega,r,mw_dot)

c	routine du calcul de la perte/gain de moment cinétique local

c	routine private du module mod_evol

c	Auteur: P.Morel, Département Cassiopée, O.C.A.
c	CESAM2k

c----------------------------------------------------------------

	USE mod_donnees, ONLY : langue, pw_extend, p_pertw, rsol, secon6
	USE mod_kind
	USE mod_variables, ONLY : mstar
		
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: nu, omega, r	
	REAL (kind=dp), INTENT(out) :: mw_dot
	
	REAL (kind=dp), SAVE :: mlim_w
	
	LOGICAL, SAVE :: init=.TRUE.	
		
c---------------------------------------------------------------- 

	IF(init)THEN
	 init=.FALSE.

	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1001)p_pertw,pw_extend ; WRITE(2,1001)p_pertw,pw_extend
1001	  FORMAT('The ang. momentum loss obeys to a local law, gamma=',
	1  es10.3,/,'localized in the external mass fraction=',es10.3,/)
	 CASE DEFAULT	 
	  WRITE(*,1)p_pertw,pw_extend ; WRITE(2,1)p_pertw,pw_extend
1	  FORMAT('Variation du moment cinétique local, gamma=',es10.3,/,
	1   'fraction de masse de la partie externe concernée=',es10.3,/)
	 END SELECT
	 mlim_w=(mstar*(1.d0-pw_extend))**(2.d0/3.d0)
	ENDIF
	
c	en deça de mlim_w il n'y a pas d'apport de moment cinétique		 
	IF(nu < mlim_w)THEN
	 mw_dot=0.d0
	ELSE 
	 mw_dot=p_pertw*(r*rsol)**2*omega**2
	ENDIF
			
	RETURN
	
	END SUBROUTINE pertw_loc
