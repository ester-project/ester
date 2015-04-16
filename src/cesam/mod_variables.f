
c************************************************************************

      MODULE mod_variables

c ce MODULE regroupe les paramètres et les variables de CESAM2k

c La signification des variables est décrite dans l'annexe E4 de la notice
c de CESAM2k

c Auteur: P.Morel Département J.D. Cassini O.C.A.

c-------------------------------------------------------------------------

	USE mod_donnees, ONLY: nrot, pnzc
	USE mod_kind

	IMPLICIT NONE

c variables public:	
	REAL (kind=dp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: bp,
	1 bp_t, chim, chim_t, old_ptm, rota, rota_t, tds, tds_t, vth, vth_t
	REAL (kind=dp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mc,
	1 mct, mc_fixe, mc_t, mct_t, mrot, mrott, mrott_t, mrot_t, m23,
	2 m23_t, q, qt, q_t, qt_t, r2, r2_t, xl, xt_ptm, xt_tds, xt_tds_t,
	3 x_ptm, x_tds, x_tds_t, x_planet
     	REAL (kind=dp), SAVE, PUBLIC, DIMENSION(0:2*pnzc+1) :: r_zc_conv
	REAL (kind=dp), SAVE, PUBLIC, DIMENSION(2*pnzc) :: m_zc, m_zc_t,
	1 r_zc, r_zc_t, r_ov, r_ov_t

	REAL (kind=dp), SAVE, PUBLIC :: age, c_iben, mstar, mstar_t,
	1 mw_tot=0.d0, mw_tot_t, psi0, rstar, wrot, wrot_t

	INTEGER, SAVE, PUBLIC, DIMENSION(2*pnzc) :: jlim, jlim_t
	INTEGER, SAVE, PUBLIC :: dim_ch, dim_qs, dim_rot, id_conv, if_conv,
	1 knot, knotc, knotc_t, knot_ptm, knotr=0, knotr_t, knot_t,
	2 knot_tds, knot_tds_t, lim, lim_t, model_num=-1, nb_modeles=-1,
	3 nc_fixe, n_ch, n_ch_t=-100, n_ptm=-100, n_qs, n_qs_t=-100, n_rot=0,
	4 n_rot_t, n_tds, n_tds_t

	LOGICAL, SAVE, PUBLIC, DIMENSION(0:2*pnzc) :: lconv=.FALSE.,
	1 lconv_t=.FALSE.
	LOGICAL, SAVE, PUBLIC :: lhe_stop=.FALSE., lt_stop=.FALSE.,
	1 lx_stop=.FALSE., tot_conv, tot_rad

c fonctions PUBLIC: 
c chim_gram : transformation des abondances /mole ==> /gramme 
c inter : recherche des variables à m^23 ou r^2 (m ou r) cte

	PRIVATE
	PUBLIC :: chim_gram, inter, sortie

	CONTAINS
	
c----------------------------------------------------------------------	
	
	INCLUDE 'chim_gram.f'
	INCLUDE 'inter.f'
	INCLUDE 'sortie.f'
	
	END MODULE mod_variables
