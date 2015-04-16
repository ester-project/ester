
c*****************************************************************

	 REAL (kind=dp) FUNCTION tacho_full(r)
	 
c routine subordonnée de la routine générique difft.f	 

c calcul du coefficient de diffusion turbulente dans la tachocline
c de la plus proche ZC, selon Castro & al. A&A 463,755, 2007
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
	
	INTEGER :: i
	
	LOGICAL, SAVE :: init=.TRUE.
	
c------------------------------------------------------------------------
	
2000	FORMAT(8es10.3)

	IF(init)THEN
	 init=.FALSE.
c	 dbzc=8.d5	!valeur des Toulousains, trop fort Li7=7d-11 mi pms	 
	 dbzc=5.d4	!valeur pour modèle solaire CESAM2k calibré avec PMS
	 delta=0.048d10/rsol	![Li7]=1
	ENDIF
	
c dzc distance de la plus proche limite ZR/ZC
	dzc=HUGE(1.d0)
	DO i=1,lim	
	 dzc=MIN(dzc, ABS(r-r_zc(i)))	
	ENDDO
	
c en fraction de l'épaisseur de la tachocline solaire	
	dzc=dzc/delta
	
c coeff de diffusion turbulente	
	tacho_full=dbzc/2.d0**dzc
		
	RETURN
	
	END FUNCTION tacho_full
