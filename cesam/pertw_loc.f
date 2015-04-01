	
c************************************************************************

         SUBROUTINE pertw_loc(r,mw_dot)

c	routine du calcul de la perte/gain de moment cinétique local

c	routine private du module mod_evol

c	Auteur: P.Morel, Département Cassiopée, O.C.A.
c	CESAM2k

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
	  WRITE(*,1001)p_pertw ; WRITE(2,1001)p_pertw
1001	  FORMAT('The ang. momentum loss obeys to a local law, gamma=',
	1 es10.3,/,'localized in the external CZ',/)
	 CASE DEFAULT	 
	  WRITE(*,1)p_pertw ; WRITE(2,1)p_pertw
1	  FORMAT('Variation du moment cinétique local, gamma=',es10.3,/,
	1 'localisée dans la ZC ext',/)
	 END SELECT
	  
	 coeff=p_pertw*rsol**2 
	 
	ENDIF
	
	mw_dot=coeff*r**2*rota(1,dim_rot)**2
			
	RETURN
	
	END SUBROUTINE pertw_loc
