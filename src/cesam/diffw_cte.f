
c*************************************************************

	SUBROUTINE diffw_cte(deff,dh,dv)

c routine private du module mod_evol

c calcul des coefficients de diffusion turbulente pour la rotation
c Dh, Dv, Deff, sont constants


c sorties :
c	Deff, Dh, Dv : notations évidentes

c	Auteur: P. Morel, Département Cassiopée, O.C.A., CESAM2k

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : langue, nom_diffw
	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp), INTENT(out) :: deff, dh, dv
	
	LOGICAL, SAVE :: init=.TRUE.
				
c------------------------------------------------------------------------
	
2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE.
	 SELECT CASE(langue)	 
	 CASE('english')
	  WRITE(*,1001)nom_diffw ; WRITE(2,1001)nom_diffw
1001	  FORMAT('Fixed coefficients for the diffusion of angular momentum',a)
	 CASE DEFAULT
	  WRITE(*,1)nom_diffw ; WRITE(2,1)nom_diffw
1	  FORMAT('Coefficients de diffusion du moment cinétique constants',a)
	 END SELECT
	ENDIF
	
c coefficients
	deff=3.d2 ; dh=1.d6 ; dv=2.5d2
	 
	RETURN
	
	END SUBROUTINE diffw_cte
