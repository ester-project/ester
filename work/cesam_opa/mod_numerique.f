
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

	INCLUDE 'numeric/arb_rom.f'
	INCLUDE 'numeric/boite.f'
	INCLUDE 'numeric/box.f'		
	INCLUDE 'numeric/bsp1dn.f'
	INCLUDE 'numeric/bsp1ddn.f'
	INCLUDE 'numeric/bsp_dis.f'
	INCLUDE 'numeric/bsp_gal.f'	
	INCLUDE 'numeric/bval0.f'
	INCLUDE 'numeric/bval1.f'
	INCLUDE 'numeric/bvald.f'
	INCLUDE 'numeric/coll.f'
	INCLUDE 'numeric/colpnt.f'
	INCLUDE 'numeric/delete_doubles.f'
	INCLUDE 'numeric/difdiv.f'    
	INCLUDE 'numeric/fermi_dirac.f'     
	INCLUDE 'numeric/gauss_band.f'
	INCLUDE 'numeric/horner.f'
	INCLUDE 'numeric/intgauss.f'
	INCLUDE 'numeric/left_right.f'
	INCLUDE 'numeric/linf.f'
	INCLUDE 'numeric/matinv.f'
	INCLUDE 'numeric/max_local.f'
	INCLUDE 'numeric/min_max.f'
	INCLUDE 'numeric/neville.f'  
	INCLUDE 'numeric/newspl.f'
	INCLUDE 'numeric/newspl_gal.f'
	INCLUDE 'numeric/newton.f'
	INCLUDE 'numeric/noedif.f'
	INCLUDE 'numeric/noein.f'
	INCLUDE 'numeric/noeud.f'	
	INCLUDE 'numeric/noeu_dis.f'
	INCLUDE 'numeric/pause.f'       
	INCLUDE 'numeric/polyder.f'
	INCLUDE 'numeric/schu58_n.f' 
	INCLUDE 'numeric/shell.f'
	INCLUDE 'numeric/sum_n.f'
	INCLUDE 'numeric/zoning.f'

	END MODULE mod_numerique
