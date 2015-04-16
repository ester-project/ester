
c******************************************************************

	MODULE mod_evol

c module contenant les routines concernant l'initialisation
c et l'évolution temporelle de la composition chimique avec
c et sans diffusion à l'exclusion des routines relatives au nucléaire
c ainsi que celles relatives à la diffusion du moment cinétique

c La signification des variables est décrite au paragraphe F11 de la notice
c de CESAM2k

c Auteur: P.Morel, Département J.D. Cassini, O.C.A., CESAM2k

c NOTATIONS (hélas incohérentes) pour les développements en B-splines
c n_ch : nombre VARIABLE de points élément de mod_variables
c nch : nombre FIXE de fonctions élément de mod_donnees
c m_ch : ordre FIXE des splines élément de mod_donnees 
c mch(n_ch) : abscisses VARIABLES élément de mod_variables

c eps : taux des écarts à droite pour les discontinuités
c mini valeur MIN du coefficient de diffusion Dv
c mzc_ext : masse de la ZC externe / mstar
c nuzc_ext : m^2/3 de la limite de la ZC externe

c-------------------------------------------------------------------------

	USE mod_donnees, ONLY : pnzc
	USE mod_kind
	USE mod_nuc, ONLY : nuc

	REAL (kind=dp), PRIVATE, PARAMETER :: eps=1.d-6, nu_min=1.d-6,
	1 t_gab=1.d6, un_eps=1.d0-eps

	REAL (kind=dp), PRIVATE, SAVE, ALLOCATABLE, DIMENSION(:) :: xcoll_rot
	REAL (kind=dp), PRIVATE, SAVE, DIMENSION(2*pnzc) :: x_mix=-100.d0
	REAL (kind=dp), PRIVATE, DIMENSION(pnzc,2) :: omega_tau=-100.d0
	
	INTEGER, PRIVATE, PARAMETER :: ncoeff=30

	INTEGER, PRIVATE, SAVE, DIMENSION(0:2*pnzc) :: idis
	INTEGER, PRIVATE, SAVE, DIMENSION(pnzc+1) :: convd
	INTEGER, PRIVATE, SAVE, DIMENSION(0:pnzc) :: convf
	INTEGER, PRIVATE, SAVE :: ndis, nzc, ncoll_rot, n_mix=0, usl_evol

	LOGICAL, PRIVATE, SAVE, DIMENSION(2*pnzc) :: mix=.TRUE.
	LOGICAL, PRIVATE, SAVE :: lw_perte=.TRUE.

	PRIVATE
	PUBLIC :: coeff_rota, coeff_vth, ecrit_rota, evol, initialise_rota,
	1 rk_imps, tab_vth

	CONTAINS

c********************************************************************

	INCLUDE 'alecian1.f'
	INCLUDE 'coeff_rota.f'
	INCLUDE 'coeff_rota3.f'
	INCLUDE 'coeff_rota4.f'
	INCLUDE 'coeff_vth.f'
	INCLUDE 'collision.f'
	INCLUDE 'coulomb.f'
	INCLUDE 'diffm.f'
	INCLUDE 'diffm_br.f'
	INCLUDE 'diffm_mp.f'
	INCLUDE 'difft.f'
	INCLUDE 'difft_gab.f'
	INCLUDE 'difft_nu.f'
	INCLUDE 'difft_smc.f'
	INCLUDE 'difft_sun.f'
	INCLUDE 'diffus.f'
	INCLUDE 'diffw.f'
	INCLUDE 'diffw_cte.f'
	INCLUDE 'diffw_mpz.f'
	INCLUDE 'diffw_p03.f'
	INCLUDE 'diffw_toul.f'
	INCLUDE 'ecrit_rota.f'
	INCLUDE 'eq_diff_chim.f'
	INCLUDE 'eq_diff_rota3.f'
	INCLUDE 'eq_diff_rota4.f'
	INCLUDE 'eq_ini_rota4.f'
	INCLUDE 'evol.f'
	INCLUDE 'f_rad.f'
	INCLUDE 'initialise_rota.f'
	INCLUDE 'initialise_rota4.f'
	INCLUDE 'lmix.f'
	INCLUDE 'pertw.f'
	INCLUDE 'pertw_loc.f'
	INCLUDE 'pertw_ptm.f'
	INCLUDE 'pertw_sch.f'
	INCLUDE 'resout_chim.f'
	INCLUDE 'resout_rota.f'
	INCLUDE 'rk_imps.f'
	INCLUDE 'tab_vth.f'

	END MODULE mod_evol
