	
c************************************************************************

         SUBROUTINE pertw_sch(nu,omega,r,mw_dot)

c	routine du calcul de la perte/gain de moment cinétique selon
c	Schumanisch

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
	  WRITE(*,1000)p_pertw ; WRITE(2,1000)p_pertw
1000	  FORMAT('The ang. momentum loss obeys to a Schumanisch law',
	1   'a=',es10.3)
	 CASE DEFAULT	 
	  WRITE(*,1)p_pertw ; WRITE(2,1)p_pertw
1	  FORMAT('Perte de moment cinétique suivant la loi de Schumanisch',
	1   'a=',es10.3)
	 END SELECT
	 mlim_w=(mstar*(1.d0-pw_extend))**(2.d0/3.d0)	 
	ENDIF
	
c	en deça de mlim_w il n'y a pas d'apport de moment cinétique		 
	IF(nu < mlim_w)THEN
	 mw_dot=0.d0
	ELSE
	 mw_dot=p_pertw*(r*rsol)**2*omega**3
	ENDIF
	
	RETURN
	
	END SUBROUTINE pertw_sch
