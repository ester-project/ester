
c******************************************************************

	MODULE mod_nuc

c	module regroupant les routines de CESAM2k
c	concernant les réactions thermonucléaires
c	et les routines propres d'exploitation
c	les paramètres propres, valeurs max. des dimensions de tableaux

c	le calcul des taux des réactions est effectué par la routine
c	générique nuc
c	les réseaux de  réactions sont différenciés par leur nom: nom_nuc
c	lu par lit_nl, fonction publique du module mod_donnees

c paramètre private:
c	nelem_ini : nombre d'éléments chimiques initiaux

c variables privates:
c	ab : mixture, initialisé dans abon_ini
c	abon_rela : abondance relative des éléments dans Z, initialisé
c	dans abon_ini
c	m : masses atomiques des éléments, initialisé dans abon_ini 
c	c : charges des éléments, initialisé dans abon_ini

c	be7sbe9 : rapport isotopique Be7/Be9, initialisé dans abon_ini 
c	be7sz : rapport isotopique Be7/Z, initialisé dans abon_ini
c	c13sc12 : rapport isotopique C13/C12, initialisé dans abon_ini
c	h2sh1 : rapport isotopique H2/H1, initialisé dans abon_ini 
c	he3she4 : rapport isotopique He3/He4, initialisé dans abon_ini
c	he3she4z : rapport isotopique He3/He4Z, initialisé dans abon_ini
c	li6sli7 : rapport isotopique Li6/Li7, initialisé dans abon_ini
c	mg25smg24 : rapport isotopique Mg25/Mg24, initialisé dans abon_ini
c	mg26smg25 : rapport isotopique Mg26/Mg25, initialisé dans abon_ini
c	ne22sne20 : rapport isotopique Ne22/Ne20, initialisé dans abon_ini
c	nom_abon : type d'abondance initiale , initialisé dans lit_nl
c	n15sn14 : rapport isotopique N15/N14, initialisé dans abon_ini 
c	o17so16 : rapport isotopique O17/O16, initialisé dans abon_ini
c	izz : charges des noyaux utilisés, initialisé dans tabul_nuc 
c	nreac : nombre de réactions thermonucléaires, initialisé dans
c	tabul_nuc 	
c	elem : noms des éléments, initialisé dans abon_ini

c variables public:
c	t_sup : température maximale des tabulations des réactions
c	nucléaires, initialisé dans tabul_nuc

c fonctions private:
c	abon_ini : détermination des abondances initiales	
c	iben : routine de réactions nucléaires fictive, détermination de
c	la constante de contraction PMS
c	pp1 : cycle PP1 simplifié
c	ppcno10Fe : cycle PPCNO 10 éléments + Fe56
c	ppcno10K : cycle PPCNO 10 éléments + K
c	ppcno10 : cycle PPCNO 10 éléments
c	ppcno11 : cycle PPCNO 11 éléments
c	ppcno12Be : cycle PPCNO 12 éléments + Be9
c	ppcno12Li : cycle PPCNO 12 éléments + Li6
c	ppcno12 : cycle PPCNO 12 éléments
c	ppcno3a9 : cycle PPCNO + 3alpha 9 éléments
c	ppcno3ac10 : cycle PPCNO + 3alpha + carbone 10 éléments
c	ppcno9 : cycle PPCNO 9 éléments 
c	rq_reac : interpolation des taux de réaction et effet d'écran
c	tabul_nuc : tabulation des réactions thermonucléaires 
c	taux_nuc : calcul des taux des réactions thermonucléaires 	

c fonction public:
c	nuc : routine générique de réactions thermonucléaires

c	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c--------------------------------------------------------------------

	USE mod_kind
	
	IMPLICIT NONE
	
	INTEGER, PRIVATE, PARAMETER :: nelem_ini=28, niso_tot=28,
	1 nreac_tot=45
	REAL (kind=dp), SAVE, PRIVATE, DIMENSION(nelem_ini) :: ab,
	1 abon_rela, m, c
	REAL (kind=dp), SAVE, PRIVATE ::  be7sbe9, be7sz, c13sc12, h2sh1,
	1 he3she4, he3she4z, li6sli7, mg25smg24, mg26smg25, ne22sne20,
	2 nom_abon, n15sn14, o17so16
	REAL (kind=dp), SAVE, PUBLIC :: t_sup	
	INTEGER, SAVE, PRIVATE, ALLOCATABLE, DIMENSION(:,:) :: izz	
	INTEGER, SAVE, PRIVATE :: nreac
	CHARACTER (len=2), PRIVATE, DIMENSION(nelem_ini) :: elem
			
	PRIVATE
!	PUBLIC :: nuc, vent
	PUBLIC :: nuc, vent, rq_reac, tabul_nuc

	CONTAINS

c------------------------------------------------------------------------

	INCLUDE 'abon_ini.f'
	INCLUDE 'iben.f'
	INCLUDE 'nuc.f'
c	INCLUDE 'neutrinos.f'				!YLD
	INCLUDE 'pp1.f'
	INCLUDE 'pp3.f'	
	INCLUDE 'ppcno10BeBFe.f'	
	INCLUDE 'ppcno10Fe.f'
	INCLUDE 'ppcno10K.f'
	INCLUDE 'ppcno10.f'	
	INCLUDE 'ppcno11.f'
	INCLUDE 'ppcno12Be.f'
	INCLUDE 'ppcno12BeBFe.f'
	INCLUDE 'ppcno12Li.f'
	INCLUDE 'ppcno12.f'
	INCLUDE 'ppcno3a9.f'
	INCLUDE 'ppcno3ac10.f'
	INCLUDE 'ppcno9.f'
	INCLUDE 'ppcno9Fe.f'	
	INCLUDE 'rq_reac.f'
	INCLUDE 'tabul_nuc.f'
	INCLUDE 'taux_nuc.f'
	INCLUDE 'vent.f'	 	 	 	 	 	
	
	END MODULE mod_nuc
