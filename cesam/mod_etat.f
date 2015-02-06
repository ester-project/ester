
c******************************************************************

	MODULE mod_etat

c	MODULE regroupant les routines de CESAM2k
c	concernant la physique de l'équation d'état
c	la routine etat est la routine générique pour le calcul des EOS
c	qui sont différenciées par leur nom: nom_etat est lu dans lit_nl
c	et mis dans mod_donnees

c fonctions private
c	etat_ceff: EOS de Eggleton et al., avec corrections coulombiennes,
c	appel à etat_eff en cas de Pb
c	etat_eff: EOS de Eggleton et al., appel à etat_gong2 en cas de Pb
c	etat_gong1: EOS de GONG1, H et He4 sont seuls considérés, et
c	supposés tot. ionisés
c	etat_gong2: EOS de GONG2, H et He4 sont seuls considérés
c	etat_mhd: EOS Mihalas, Hummer & Dappen, appel à etat_eff si Pb
c	etat_opal: EOS de OPAL, appel à etat_eff en cas de Pb

c fonction public:

c	df_rotx : transformation des dérivées partielles P,T,X --> ro,T,X
c	etat : fonction générique de calcul de l'EOS
c	saha : équation de Saha

c 	Auteurs: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c	fonctions private:

	PRIVATE

c	function public:

	PUBLIC :: df_rotx, etat, saha

	CONTAINS

c-----------------------------------------------------------------------
	 
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
	 INCLUDE 'etat_opal5Z.f'
	 INCLUDE 'saha.f'

	END MODULE mod_etat
