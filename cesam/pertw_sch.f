	
c************************************************************************

         SUBROUTINE pertw_sch(r,mw_dot)

c routine du calcul de la perte/gain de moment cinétique selon Schumanisch

c routine private du module mod_evol

c Auteur: P.Morel, Département Cassiopée, O.C.A. CESAM2k

c----------------------------------------------------------------

	USE mod_donnees, ONLY : langue, p_pertw, rsol
	USE mod_kind
	USE mod_variables, ONLY : dim_rot, rota
				
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: r	
	REAL (kind=dp), INTENT(out) :: mw_dot
	
	REAL (kind=dp), SAVE :: coeff
	
	LOGICAL, SAVE :: init=.TRUE.	
		
c---------------------------------------------------------------- 

	IF(init)THEN
	 init=.FALSE.
	 lw_perte=.TRUE. 
	 SELECT CASE(langue)
	 CASE('english')
	  WRITE(*,1000)p_pertw ; WRITE(2,1000)p_pertw
1000	  FORMAT('The ang. momentum loss obeys to a Schumanisch law, a=',es10.3)
	 CASE DEFAULT	 
	  WRITE(*,1)p_pertw ; WRITE(2,1)p_pertw
1	  FORMAT('Perte de moment cinétique, loi de Schumanisch, a=',es10.3)
	 END SELECT	 
	 coeff=	p_pertw*rsol**2	 
	ENDIF
	
	mw_dot=coeff*r**2*rota(1,dim_rot)**3
	
	RETURN
	
	END SUBROUTINE pertw_sch
