
c******************************************************************

	MODULE mod_etat

c	MODULE regroupant les routines de CESAM2k
c	concernant la physique de l'équation d'état
c	la routine etat est la routine générique pour le calcul des EOS

c La signification des variables fait l'objet du paragraphe F5 de la notice
c	 de CESAM2k

c 	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c-----------------------------------------------------------------

	PRIVATE
c	PUBLIC :: coeff_vth, df_rotx, etat, mu_mol, saha
	PUBLIC :: df_rotx, etat, mu_mol, saha
	
	
	CONTAINS

c-----------------------------------------------------------------------
	 
c	 INCLUDE 'coeff_vth.f'
	 INCLUDE 'df_rotx.f'
	 INCLUDE 'etat.f'	 
	 INCLUDE 'etat_ceff.f'
	 INCLUDE 'etat_eff.f'
	 INCLUDE 'etat_gong1.f'		
	 INCLUDE 'etat_gong2.f'	
	 INCLUDE 'etat_mhd.f'
	 INCLUDE 'etat_opal.f'
	 INCLUDE 'etat_opalX.f'	 
	 INCLUDE 'etat_opalZ.f'
	 INCLUDE 'mu_mol.f'
	 INCLUDE 'saha.f'

	END MODULE mod_etat
