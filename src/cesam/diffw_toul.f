
c*************************************************************

	SUBROUTINE diffw_toul(r,y,deff,dh,dv)

c routine private du module mod_evol

c calcul des coefficients de diffusion turbulente pour la rotation
c Dh, Dv, Deff, suivant Castro et al. A&A 463,755, 2007

c entrées :
c	dlnro : d ln ro / dnu
c	grand_k : K
c	nu : m^2/3
c	nu_cin : viscosité cinématique
c	y(:,0:1) : variables propres à la rotation et dérivées / nu

c sorties :
c	Deff, Dh, Dv : notations évidentes

c	Auteur: P. Morel, Département Cassiopée, O.C.A., CESAM2k

c---------------------------------------------------------------------

	USE mod_donnees, ONLY : langue, nom_diffw, rsol 
	USE mod_kind

	IMPLICIT NONE

	REAL (kind=dp), INTENT(in), DIMENSION(2,0:1) :: y	
	REAL (kind=dp), INTENT(in) :: r
	REAL (kind=dp), INTENT(out) :: deff, dh, dv
	
	REAL (kind=dp), PARAMETER :: alfa_turb=6.d0
	
	REAL (kind=dp), SAVE :: ch, cv 
	
	LOGICAL, SAVE :: init=.TRUE.
			
c------------------------------------------------------------------------
	
2000	FORMAT(8es10.3)

c initialisations
	IF(init)THEN
	 init=.FALSE.
	 ch=5.d4 
	 cv=alfa_turb-1.d0/30.d0/ch
	 	 
	 SELECT CASE(langue)	 
	 CASE('english')
	  WRITE(*,1001)nom_diffw ; WRITE(2,1001)nom_diffw
1001	  FORMAT('Diffusion of angular momentum of Toulouse',a,/,
	1 'Castro et al. A&A 463,755, 2007')
	 CASE DEFAULT
	  WRITE(*,1)nom_diffw ; WRITE(2,1)nom_diffw
1	  FORMAT('Coefficients de diffusion du moment cinétique toul',a,/,
	1 'Castro et al. A&A 463,755, 2007')
	 END SELECT
	 
	ENDIF

c Dh
	dh=rsol*ch*r*ABS(y(2,0))

c Dv
	dv=rsol*cv*r*ABS(y(2,0))

c Deff
	deff=rsol*alfa_turb*r*ABS(y(2,0))

	RETURN
	
	END SUBROUTINE diffw_toul
