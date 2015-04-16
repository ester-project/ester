
c*****************************************************************

	SUBROUTINE ctes_94
	
c	routine private du module mod_donnees
	
c	routine d'initialisation des constantes physiques

c	constantes physiques selon  CRC Handbook of
c	chemistry and physics 75th edition 1994 D.R. Lide et al. CRC press
c	Boka Raton : Ann Arbor : London : Tokyo

c	masses atomiques en amu	de la compilation NACRE
c	valeurs solaires de Guenter & al

c	la routine print_ctes permet de les écrire
c	CESAM2k

c	Auteur: P.Morel, Département J.D. Cassini : O.C.A.

c-------------------------------------------------------------------

	USE mod_kind
		
	IMPLICIT NONE
	
c-------------------------------------------------------------------

c	origine des principales données
	
	source='Handbook chem. phys. 94 + NACRE'
	
c	données

	amu=1.6605402d-24		!masse atom. unité, Avogadro=1/amu
	aradia=7.565912199839849d-15	!constante de la radiation
	clight=2.99792458d10		!célérité de la lumière
	echarg=4.8032068d-10		!charge de l'électron
	eve=1.60217733d-12		!électron volt
	gmsol=1.32712438d26		!produit G Msol
	hpl=6.6260755d-27		!planck	
	kbol=1.380658d-16		!boltzman
	granr=kbol/amu			!constante des gaz parfaits	
	ln10=LOG(10.d0)
	lbol0=3.055d33			!point 0 des Mbol	
	me=9.1093897d-28		!masse électron			
	pi=acos(-1.d0)
	secon6=365.25d0*86400.d6       !nb. de s. dans 1Myr (norme de 1967)
	sigma=aradia*clight/4.d0	!constante de Stefan
	zsx_sol=0.0245
	fesh_sol=7.5d0-12.00d0		![Fe/H] de GN 93
	
	lsol=3.846d33			!luminosité solaire		
	msol=1.98919d33			!masse solaire	
	rsol=6.9599d10			!rayon solaire
	
	g=gmsol/msol			!gravité
	
	an=1.008665d0		!masse atomique en amu du neutron
	ap=1.00727647d0		!masse atomique en amu du proton
	ah=1.007825d0		!masse atomique en amu de l'hydrogène
	ah2=2.0141018d0		!masse atomique en amu du deutérium
	ahe3=3.0160293d0	!masse atomique en amu de l'hélium 3
	ahe4=4.0026033d0	!masse atomique en amu de l'hélium 4
	ali6=6.015121d0		!masse atomique en amu du lithium 6	
	ali7=7.0160040d0	!masse atomique en amu du lithium 7	
	abe7=7.0169292d0	!masse atomique en amu du béryllium 7
	abe9=9.0121821d0	!masse atomique en amu du béryllium 9
	ab11=11.0093055d0       !masse atomique en amu du bore 11	
	ac12=12.d0		!masse atomique en amu du carbone 12
	ac13=13.0033548d0	!masse atomique en amu du carbone 13
	an13=13.0057386d0	!masse atomique en amu de l'azote 13
	an14=14.003074d0	!masse atomique en amu de l'azote 14
	an15=15.001089d0	!masse atomique en amu de l'azote 15
	ao16=15.9949146d0	!masse atomique en amu de l'oxygène 16
	ao17=16.9991315d0	!masse atomique en amu de l'oxygène 17
	ao18=17.9991604d0	!masse atomique en amu de l'oxygène 18
	af18=18.0009377d0	!masse atomique en amu du fluor 18
	af19=18.9984032d0	!masse atomique en amu du fluor 19
	ane20=19.9924402d0	!masse atomique en amu du néon 20
	ane21=20.9938467d0	!masse atomique en amu du néon 21
	ane22=21.9913855d0	!masse atomique en amu du néon 22
	ana23=22.9897697d0	!masse atomique en amu du sodium 23
	amg23=22.9941249d0	!masse atomique en amu du magnésium 23
	amg24=23.9850419d0	!masse atomique en amu du magnésium 24
	amg25=24.985837d0	!masse atomique en amu du magnésium 25
	amg26=25.982593d0	!masse atomique en amu du magnésium 26
	afe56=55.847d0		!masse atomique en amu du fer 56

	RETURN
			
c----------------------------------------------------------------------
	
	END SUBROUTINE ctes_94
	
