
c******************************************************************

	MODULE mod_kind
	
c	module de définition des types
c	sp: simple précision, dp: double précision

c	Auteur: P.Morel + B. Pichon
c	Laboratoire JD. Cassini, OCA

c-----------------------------------------------------------------
	
	INTEGER, PARAMETER, public :: dp=kind(1.d0), sp=kind(1.)
	
	END module mod_kind
