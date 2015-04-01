
c******************************************************************

	SUBROUTINE tab_vth

c subroutine PUBLIC du module mod_evol

c formation de l'interpolation de variables thermodynamiques
c sur la base mc, mct, knotc
c vth(1)=lnP, vth(2)=lnT, vth(3)= r**2, vth(4)=ln ro, vth(5)=ln µ

c Auteur: P.Morel, Département Cassiopée, O.C.A., CESAM2k

c-----------------------------------------------------------------------

	USE mod_donnees, ONLY :  m_ch, nvth
	USE mod_kind	
	USE mod_numerique, ONLY : bsp_gal
	USE mod_variables, ONLY : dim_ch, knotc, mct, vth
	
	IMPLICIT NONE
	
c----------------------------------------------------------------------

2000	FORMAT(8es10.3)

	IF(ALLOCATED(vth))DEALLOCATE(vth)		
	ALLOCATE(vth(nvth,dim_ch))

c formation des coefficients
	CALL bsp_gal(nvth,vth,m_ch,mct,knotc,coeff_vth)
	
	RETURN

	END SUBROUTINE tab_vth
