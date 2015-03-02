
c*****************************************************************

	MODULE mod_atm

c	Module regroupant les routines de restitution de l'atmosphère

c 	Auteur: P. Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k	

c variables private:
c	delfim, delfip : Delta phi -, Delta phi +
c	ltaue, ltauf : ln tau_max, ln tau_min
c	tau_min : épaisseur optique de la couche la plus externe
c	rad=.true. : la loi T(tau) est purement radiative
c	ne_atm : nombre d'équations pour la restitution de l'atmosphère
c	n23_atm : indice de la couche où T(tau*)=Teff

c variables public:
c	dlpp_atm : tableau des gradients d ln T/ d ln P
c	m_atm : tableau des masses, abscisses lagrangiennes  
c	p_atm : tableau des pressions gazeuses
c	pt_atm : tableau des pressions totales
c	r_atm : tableau des rayons abscisses euleriennes 
c	tau : tableau des profondeurs optiques
c	t_atm : tableau des températures

c fonction PUBLIC:
c	atm : fonction générique du calcul d'atmosphère
c	tdetau : fonction générique de loi T(tau)
c	thermo_atm : fonctions thermo. pour l'atmosphère

c fonctions PRIVATE:

c	coll_atm : intégration des équations de restitution de
c	l'atmosphère
c	edding : loi T(tau) d'Eddington
c	eq_atm : calcul des coefficients des équations de restitution de
c	l'atmosphère
c	hopf : loi T(tau) de Hopf
c	k5750, k5777 : loi T(tau) solaires de Kurucz Teff=5750, 5777
c	lim_atm : routine de gestion de restitution de l'atmosphère
c	lim_gong1 : routine de restitution 1 couche de l'atmosphère
c	cas GONG1
c	lim_tau1 : routine de restitution 1 couche de l'atmosphère
c	roger00, roger02, roger05, roger10a : loi T(tau) de Kurucz,
c	[Fe/H]=0, -0.2, -0.5, -1 enhanced
c	taueff : détermination de tau* tel que T(tau*)=Teff
c	trho : interpolation de T et ro pour lois t(tau) roger*

c-----------------------------------------------------------------

	USE mod_kind
	
	IMPLICIT NONE

	REAL (kind=dp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:):: bp_atm,
	1 bp_atm_t	
	REAL (kind=dp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:) :: dlpp_atm,
	1 m_atm, p_atm, pt_atm, r_atm, tau, t_atm, x_atm, x_atm_t, xt_atm,
	2 xt_atm_t
	REAL (kind=dp), SAVE, PRIVATE :: delfim, delfip, ltaue, ltauf,
	1 tau_min
	
	INTEGER, SAVE, PRIVATE :: Ipgt, ne_atm, n23_atm

	LOGICAL, SAVE, PRIVATE :: rad
	
	CHARACTER (len=3), PARAMETER, PRIVATE, DIMENSION(8) :: nom_atm=
	1 (/ ' Pt', ' T ', ' R ', ' R*', ' M ', 'to*', ' to', ' Pg' /)
		 
	PRIVATE
	PUBLIC :: atm, tdetau, thermo_atm
	
	CONTAINS	

c----------------------------------------------------------------

	INCLUDE 'atm.f'
	INCLUDE 'coll_atm.f'
    INCLUDE 'edding.f'
	INCLUDE 'eq_atm.f'
	INCLUDE 'hopf.f'
	INCLUDE 'k5750.f'
	INCLUDE 'k5777.f'
	INCLUDE 'kurucz.f'
	INCLUDE 'lim_atm.f'
	INCLUDE 'lim_gong1.f'		
	INCLUDE 'lim_tau1.f'
	INCLUDE 'roger00.f'
	INCLUDE 'rogerYL.f'
	INCLUDE 'roger02.f'
	INCLUDE 'roger05.f'
	INCLUDE 'roger10a.f'
	INCLUDE 'taueff.f'
	INCLUDE 'tdetau.f'	
	INCLUDE 'thermo_atm.f'
	INCLUDE 'trho.f'				
	INCLUDE 'trhoYL.f'				

	END MODULE mod_atm
