
c******************************************************************

	MODULE mod_nuc

c module regroupant les routines concernant les réactions thermonucléaires

c le calcul des taux des réactions est effectué par la routine	générique nuc
c les réseaux de  réactions sont différenciés par leur nom: nom_nuc
c lu par lit_nl, fonction publique du module mod_donnees

c La signification des variables est décrite au paragraphe F9 de la notice
c	 de CESAM2k

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c--------------------------------------------------------------------

	USE mod_kind

	IMPLICIT NONE

	INTEGER, PUBLIC, PARAMETER :: m_temp=4, niso_tot=32, nreac_tot=64
	INTEGER, PRIVATE, PARAMETER :: nelem_ini=28	

	REAL (kind=dp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: taux_reac
	REAL (kind=dp),SAVE, PUBLIC , ALLOCATABLE, DIMENSION(:) :: ar, q0, temp,
	1 ttemp	
	REAL (kind=dp), SAVE, PRIVATE, DIMENSION(nelem_ini) :: ab, abon_rela,
	1 m, c
	REAL (kind=dp), SAVE, PRIVATE :: be7sbe9, be7sz, c13sc12, h2sh1,
	1 he3she4, he3she4z, li6sli7, mg25smg24, mg26smg25, ne22sne20,
	2 n15sn14, o17so16, o18so16
	REAL (kind=dp), SAVE, PUBLIC :: age_deb, age_fin, dt_planet, t_sup
	REAL (kind=dp), PUBLIC :: mzc_ext, nuzc_ext

	INTEGER, SAVE, PRIVATE, ALLOCATABLE, DIMENSION(:,:) :: izz

c i3al : indice de la réaction 3alpha			
	INTEGER, SAVE, PRIVATE :: i3al=0
	INTEGER, SAVE, PUBLIC ::  knot_temp, nreac, n_temp


	LOGICAL, SAVE, PUBLIC :: l_planet, l_vent

	CHARACTER (len=2), SAVE, PRIVATE, DIMENSION(nelem_ini) :: elem

	PRIVATE
	PUBLIC :: abon_ini, nuc, planetoides, taux_nuc, vent

	CONTAINS

c------------------------------------------------------------------------

	INCLUDE 'abon_ini.f'
	INCLUDE 'iben.f'
	INCLUDE 'nuc.f'
	INCLUDE 'planetoides.f'
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
	INCLUDE 'ppcno3a12Ne.f'
	INCLUDE 'ppcno3a9.f'
	INCLUDE 'ppcno3aco.f'
	INCLUDE 'ppcno9.f'
	INCLUDE 'ppcno9Fe.f'
	INCLUDE 'rq_reac.f'
	INCLUDE 'tabul_nuc.f'
	INCLUDE 'taux_nuc.f'
	INCLUDE 'vent.f'

	END MODULE mod_nuc
