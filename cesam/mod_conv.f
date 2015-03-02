
c*******************************************************************

	MODULE mod_conv
      
c	Module regroupant les routines relatives à la convection
c 	Auteur : P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k	
	
c fonctions PRIVATE:

c	conv_a0 : convection MLT, avec longueur de melange tendant
c       vers 0 aux limites ZR/ZC
c	conv_cgm_reza : calcul du gradient convectif selon CGM
c	conv_cm : convection suivant Canuto & Mazitelli avec l = alpha Hp
c	conv_cml : calcul du gradient convectif selon Canuto Mazitelli
c	avec longueur de melange = + courte distance des bords de la ZC
c	conv_cm_reza : convec. Canuto & Mazitelli adapté par S.Reza
c	conv_jmj : convection MLT avec l = alpha Hp

c fonction PUBLIC:
c	conv : routine générique de convection

	PRIVATE     
	PUBLIC :: conv

	CONTAINS

c------------------------------------------------------------------

	INCLUDE 'conv.f'       
	INCLUDE 'conv_a0.f' 
	INCLUDE 'conv_cgm_reza.f'
c ----------------------------------------------------------------------------
	INCLUDE 'conv_cgm_reza_r.f'
c ----------------------------------------------------------------------------
	INCLUDE 'conv_cm.f'
	INCLUDE 'conv_cml.f'
	INCLUDE 'conv_cm_reza.f'
	INCLUDE 'conv_jmj.f'

	END MODULE mod_conv
