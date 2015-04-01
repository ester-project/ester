
c******************************************************************

	MODULE mod_numerique

c	module regroupant les outils numériques et utilitaires pour CESAM2k

c La signification des variables est décrite au paragraphe F2 de la notice
c	 de CESAM2k

c	Auteurs: P.Morel, B.Pichon Département J.D. Cassini, O.C.A.
c	CESAM2k

c--------------------------------------------------------------------

	USE mod_kind
	
	LOGICAL, SAVE, PUBLIC :: no_croiss=.FALSE.

	PUBLIC :: arb_rom, boite, box, bsp1ddn, bsp1dn, bsp_dis, bsp_gal,
	1 bvald, bval1, coll, delete_doubles, fermi_dirac,
	2 gauss_band, intgauss, left_right, linf, matinv, min_max,
	3 max_local, neville, newton, noedif, noein, noeud, newspl,
	4 newspl_gal, pause, polyder, shell, sum_n, zoning

	CONTAINS

c--------------------------------------------------------------------

	INCLUDE 'arb_rom.f'
	INCLUDE 'boite.f'
	INCLUDE 'box.f'		
	INCLUDE 'bsp1dn.f'
	INCLUDE 'bsp1ddn.f'
	INCLUDE 'bsp_dis.f'
	INCLUDE 'bsp_gal.f'	
	INCLUDE 'bval0.f'
	INCLUDE 'bval1.f'
	INCLUDE 'bvald.f'
	INCLUDE 'coll.f'
	INCLUDE 'colpnt.f'
	INCLUDE 'delete_doubles.f'
	INCLUDE 'difdiv.f'    
	INCLUDE 'fermi_dirac.f'     
	INCLUDE 'gauss_band.f'
	INCLUDE 'horner.f'
	INCLUDE 'intgauss.f'
	INCLUDE 'left_right.f'
	INCLUDE 'linf.f'
	INCLUDE 'matinv.f'
	INCLUDE 'max_local.f'
	INCLUDE 'min_max.f'
	INCLUDE 'neville.f'  
	INCLUDE 'newspl.f'
	INCLUDE 'newspl_gal.f'
	INCLUDE 'newton.f'
	INCLUDE 'noedif.f'
	INCLUDE 'noein.f'
	INCLUDE 'noeud.f'	
	INCLUDE 'noeu_dis.f'
	INCLUDE 'pause.f'       
	INCLUDE 'polyder.f'
	INCLUDE 'schu58_n.f' 
	INCLUDE 'shell.f'
	INCLUDE 'sum_n.f'
	INCLUDE 'zoning.f'

	END MODULE mod_numerique
