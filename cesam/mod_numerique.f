
c******************************************************************

	MODULE mod_numerique

c	module regroupant les outils numériques et utilitaires pour CESAM2k

c	Auteur: P.Morel, B.Pichon Département J.D. Cassini, O.C.A.
c	CESAM2k

c variable public:
c	no_croiss=.TRUE. : une suite des abscisses n'est pas
c	strictement croisaante

c fonctions private:
c	bval0 : calcul des B-splines non id. nulle en un point
c	de leur support
c	colpnt : calcul de l'abscisses du i-ieme point de collocation
c	difdiv : algorithme des différences divisées (pour interpolation
c	polynomiale formule de Newton)
c	horner : algorithme de Horner: calcul de la valeur d'un polynôme et
c	de ses dérivées en un point
c	noeu_dis : formation du vecteur nodal avec discontinuités
c	schu58_n : interpolation pour n fonctions par B-splines (algorithme
c	5-8 p.194 de Schumaker)

c fonctions public:
c	arb_rom : transformation numérotation arabe ==> romaine
c	boite : dessin d'une boite centrée
c	box : dessin d'une boite assymétrique
c	bsp1ddn : interpolation 1D de n fonctions par B-spline avec dérivées 
c	bsp1dn : interpolation 1D de n fonctions par B-spline
c	bsp_dis : calcul des coefficients des B-splines pour interpolation
c	avec discontinuités, contient noeu_dis
c	bvald : calcul des B-splines non id. nulle en un point
c	de leur support avec dérivées
c	bval1 : calcul des B-splines non id. nulle en un point
c	de leur support avec dérivées premières	 
c	coll : détermimation des points de collocation pour
c	intégration d'éq. diff. avec B-splines
c	fermi_dirac: calcul des intégrales de Fermi Dirac
c	gauss_band : résolution système linéaire, pivot partiel,
c	adapté au cas des matrices bande
c	genere_bases : formation d'une base avec discontinuités
c	et de la base continue associée
c	intgauss : initialisation des poids et abscisses pour intégration
c	de Gauss
c	linf : recherche de l'encadrement d'un nombre dans un tableau
c	ordonné de façon croissante
c	matinv : inversion de matrice
c	neville : interpolation de Lagrange par algorithme de Neville
c	newspl : changement de base de B-spline pour n fonctions
c	newton : intepolation polynomiale, formule de Newton
c	noedif : formation du vecteur nodal pour intégration
c	d'éq. diff. avec B-splines
c	noein : formation du vecteur nodal pour interpolation B-splines
c	noeud : formation du vecteur nodal a partir du vecteur de multiplicité
c	pause : pause avec commentaire
c	polyder : valeurs et dérivees d'un polynôme, algorithme de Horner  
c	shell : routine de tri
c	sum_n : intégration de n fonctions par B-spline
c	zoning : détermination des abcisses assurant des incréments constants
c	pour une fonction monotone tabulée

c--------------------------------------------------------------------

	USE mod_kind
	
	LOGICAL, SAVE, PUBLIC :: no_croiss=.FALSE.

	PUBLIC :: arb_rom, boite, box, bsp1ddn, bsp1dn, bsp_dis, bvald,
	1 bval1, coll, fermi_dirac, gauss_band, genere_bases, intgauss,
	2 linf, matinv, neville, newton, noedif, noein, noeud, newspl,
	3 pause, polyder, shell, sum_n, zoning

	CONTAINS

c--------------------------------------------------------------------

	INCLUDE 'arb_rom.f'
	INCLUDE 'boite.f'
	INCLUDE 'box.f'
	INCLUDE 'bsp1dn.f'
	INCLUDE 'bsp1ddn.f'
	INCLUDE 'bsp_dis.f'
	INCLUDE 'bval0.f'
	INCLUDE 'bval1.f'
	INCLUDE 'bvald.f'
	INCLUDE 'coll.f'
	INCLUDE 'colpnt.f'
	INCLUDE 'difdiv.f'    
	INCLUDE 'fermi_dirac.f'     
	INCLUDE 'gauss_band.f'
	INCLUDE 'genere_bases.f'
	INCLUDE 'horner.f'
	INCLUDE 'intgauss.f'
	INCLUDE 'linf.f'
	INCLUDE 'matinv.f'
	INCLUDE 'neville.f'  
	INCLUDE 'newspl.f'
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
