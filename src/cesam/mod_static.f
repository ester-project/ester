
c******************************************************************

	MODULE mod_static

c Module regroupant les routines permettant la résolution 
c des équations de l'équilibre quasi-statique et la gestion de l'évolution

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c--------------------------------------------------------------------

	USE mod_kind
	
	IMPLICIT NONE
		
	CHARACTER (len=3), PARAMETER, PRIVATE, DIMENSION(8) :: nom_qs=
	1(/ ' Pt', ' T ', ' R ', ' L ', ' M ', 'psi', ' ro', ' Pg' /)
	
	REAL (kind=dp), SAVE, PRIVATE, ALLOCATABLE, DIMENSION(:) :: fac, xcoll

	INTEGER, SAVE, PUBLIC, DIMENSION(7) :: iter_qs	
	INTEGER, SAVE, PRIVATE :: iter_stag, usl_static
		
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
