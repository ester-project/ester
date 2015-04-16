
c*****************************************************************

	 REAL (kind=dp) FUNCTION tacho_ext(r)
	 
c routine subordonnée de la routine générique difft.f	 

c calcul du coefficient de diffusion turbulente dans la tachocline de la ZC
c externe selon Castro & al. A&A 463,755, 2007
c cette routine n'est appelée que dans les ZR

c entrées:
c	r : abscisse R/Rsol

c sortie
c	tacho : coefficient de diffusion

c------------------------------------------------------------------

	USE mod_donnees, ONLY : rsol
	USE mod_kind
	USE mod_variables, ONLY : lim, r_zc
	
	IMPLICIT NONE
	
	REAL (kind=dp), INTENT(in) :: r
	
c épaisseur de la tachocline solaire en Rsol
	REAL (kind=dp), SAVE :: dbzc, delta

	REAL (kind=dp) :: dzc
	
	LOGICAL, SAVE :: init=.TRUE.
	
c------------------------------------------------------------------------
	
2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE.
c	 dbzc=8.d5	!valeur des Toulousains, trop fort Li7=7d-11 mi pms	 
	 dbzc=5.d4	!valeur pour modèle solaire CESAM2k calibré avec PMS
	 delta=0.048d10/rsol	![Li7]=1
	ENDIF
	
c dzc distance à la ZC externe
	IF(r < r_zc(nzc))THEN
	dzc=(r_zc(nzc)-r)/delta
	
c coeff de diffusion turbulente	
	tacho_ext=dbzc/2.d0**dzc
	
	ELSE
	 tacho_ext=1.d0
	ENDIF	
		
	RETURN
	
	END FUNCTION tacho_ext
