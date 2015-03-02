
c******************************************************************

	MODULE mod_static

c	Module regroupant les routines permettant la résolution 
c	des équations de l'équilibre quasi-statique etla gestion de
c	l'évolution

c variables private:
c	fac : facteur de répartition, initialisé dans resout
c	xcoll : table des points de collocation pour l'équilibre
c	quasi-statique, initialisé dans coll

c fonctions private:
c	coll_qs : résolution des équations de l'équilibre quasi-statique
c	dgrad : calcul du gradient de température
c	lim_zc : répartition des couches pour l'équilibre quasi-statique,
c	détermination des limites ZR/ZC
c	pertm_ext : détermination de la perte de masse
c	pertm_msol : détermination de la perte de masse cas M >= Msol
c	pertm_tot : détermination de la perte de masse tenant compte de
c	E=mc^2  
c	static_m : détermination des coefficients des équations de
c	l'équilibre quasi-static, cas lagrangien
c	static_r : détermination des coefficients des équations de
c	l'équilibre quasi-static, cas eulérien
c	update : on passe les variables du temps t+dt à t
	
c fonctions public:
c	resout : gestion de modèles initiaux, évolution temporelle de la
c	composition chimique ==> équilibre quasi-statique ==> évolution
c	temporelle de la composition chimique...
c	thermo : calcul des principales fonctions thermodynamiques, du
c	gradient convectif

c 	Auteur: P.Morel, Département J.D. Cassini, O.C.A.
c	CESAM2k

c--------------------------------------------------------------------

	USE mod_kind
	
	IMPLICIT NONE
	
	REAL (kind=dp), SAVE, PRIVATE, ALLOCATABLE, DIMENSION(:) :: fac,
	1 xcoll
	
	CHARACTER (len=3), PARAMETER, PRIVATE, DIMENSION(8) :: nom_qs=
	1(/ ' Pt', ' T ', ' R ', ' L ', ' M ', 'psi', ' ro', ' Pg' /)
	
	PRIVATE
	PUBLIC :: pertm, resout, thermo
	
	CONTAINS
	
c--------------------------------------------------------------

	INCLUDE 'coll_qs.f'
	INCLUDE 'dgrad.f'
	INCLUDE 'lim_zc.f'
	INCLUDE 'pertm.f'
	INCLUDE 'pertm_ext.f'
	INCLUDE 'pertm_msol.f'
	INCLUDE 'pertm_tot.f'
	INCLUDE 'pertm_waldron.f'		
	INCLUDE 'resout.f' 
	INCLUDE 'static.f'
	INCLUDE 'static_m.f'
	INCLUDE 'static_r.f'
	INCLUDE 'thermo.f'
	INCLUDE 'update.f'

	END MODULE mod_static
