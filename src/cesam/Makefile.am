export HDF5_CXX = $(HDF5_BACKEND)

AM_CXXFLAGS = -I. -I$(top_srcdir)/src/include
AM_LDFLAGS  = -lgfortran

lib_LTLIBRARIES	= libcesam.la

namelistdir = $(datadir)/ester/tables/cesam_nuc
namelist_DATA = test.don

libcesam_la_SOURCES = mod_evol.f mod_kind.f mod_nuc.f mod_donnees.f \
					  mod_etat.f mod_opa.f mod_bp_for_alecian.f		\
					  mod_numerique.f mod_variables.f				\
					  mod_static.f mod_atm.f mod_conv.f

mod_nuc.lo:				mod_donnees.lo mod_variables.lo				\
						abon_ini.f iben.f nuc.f planetoides.f		\
						pp1.f pp3.f ppcno10BeBFe.f ppcno10Fe.f		\
						ppcno10K.f ppcno10.f ppcno11.f ppcno12Be.f	\
						ppcno12BeBFe.f ppcno12Li.f ppcno12.f 		\
						ppcno3a12Ne.f ppcno3a9.f ppcno3aco.f 		\
						ppcno9.f ppcno9Fe.f rq_reac.f tabul_nuc.f	\
						taux_nuc.f vent.f

mod_evol.lo:			mod_kind.lo mod_nuc.lo mod_etat.lo mod_opa.lo		\
						mod_bp_for_alecian.lo mod_conv.lo					\
						alecian1.f coeff_rota.f coeff_rota3.f				\
						coeff_rota4.f coeff_vth.f collision.f coulomb.f		\
						diffm.f diffm_br.f diffm_mp.f difft.f difft_gab.f	\
						difft_nu.f difft_smc.f difft_sun.f diffus.f 		\
						diffw.f diffw_cte.f diffw_mpz.f diffw_p03.f 		\
						diffw_toul.f ecrit_rota.f eq_diff_chim.f			\
						eq_diff_rota3.f eq_diff_rota4.f eq_ini_rota4.f		\
						evol.f f_rad.f initialise_rota.f initialise_rota4.f	\
						lmix.f pertw.f pertw_loc.f pertw_ptm.f pertw_sch.f	\
						resout_chim.f resout_rota.f rk_imps.f tab_vth.f

mod_donnees.lo:			mod_numerique.lo journal.f lit_nl.f ini_ctes.f \
						print_ctes.f

mod_variables.lo:		mod_donnees.lo chim_gram.f inter.f sortie.f

mod_etat.lo:			mod_donnees.lo
mod_conv.lo:			mod_donnees.lo
mod_opa.lo:				mod_donnees.lo
mod_numerique.lo:		mod_kind.lo
mod_bp_for_alecian.lo:	mod_donnees.lo mod_variables.lo
mod_static.lo:			mod_atm.lo mod_evol.lo
mod_atm.lo:				mod_conv.lo

EXTRA_DIST= abon_ini.f alecian1.f arb_rom.f atm.f boite.f box.f bsp1ddn.f \
			bsp1dn.f bsp_dis.f bval0.f bval1.f bvald.f chim_gram.f \
			coeff_rota.f coll_atm.f coll.f collision.f coll_qs.f \
			colpnt.f conv_a0.f conv_cgm_reza.f conv_cgm_reza_r.f conv_cm.f \
			conv_cml.f conv_cm_reza0206.f conv_cm_reza.f conv.f conv_jmj.f \
			coulomb.f ctes_85.f ctes_94.f df_rotx.f dgrad.f difdiv.f \
			diffm_br.f diffm.f diffm_mp.f difft.f difft_nu.f difft_ventura.f \
			diffus.f diffw.f diffw_mpz.f edding.f eq_atm.f eq_diff_chim.f \
			eq_diff_rota.f etat_ceff.f etat_eff.f etat.f etat_gong1.f \
			etat_gong2.f etat_mhd.f etat_opal5Z.f etat_opal.f etat_opalX.f \
			etat_opalZ.f evol.f fermi_dirac.f f_rad.f from_alecian.f \
			gauss_band.f genere_bases.f hopf.f horner.f iben.f ini_ctes.f \
			inter.f intgauss.f k5750.f k5777.f kappa_cond.f kurucz.f \
			lim_atm.f lim_gong1.f lim_tau1.f lim_zc.f linf.f lit_nl.f \
			matinv.f mod_atm.f mod_bp_for_alecian.f mod_conv.f mod_donnees.f \
			mod_etat.f mod_evol.f mod_kind.f mod_nuc.f mod_numerique.f \
			mod_opa.f mod_static.f mod_variables.f neville.f newspl.f \
			newton.f noedif.f noein.f noeud.f noeu_dis.f nuc.f opa.f \
			opa_gong.f opa_houdek9.f opa_int_zsx.f opa_opal2.f opa_opalCO.f \
			opa_yveline.f pause.f pertm_ext.f pertm.f pertm_msol.f \
			pertm_tot.f pertm_waldron.f pertw.f pertw_loc.f pertw_ptm.f \
			pertw_sch.f polyder.f pp1.f pp3.f ppcno10BeBFe.f ppcno10.f \
			ppcno10Fe.f ppcno10K.f ppcno11.f ppcno12BeBFe.f ppcno12Be.f \
			ppcno12.f ppcno12Li.f ppcno3a9.f ppcno3ac10.f ppcno9.f \
			ppcno9Fe.f print_ctes.f resout_chim.f resout.f resout_rota.f \
			rk_imps.f roger00.f roger02.f roger05.f roger10a.f rogerYL.f \
			rq_reac.f saha.f schu58_n.f shell.f sortie.f static.f static_m.f \
			static_r.f sum_n.f tabul_nuc.f taueff.f taux_nuc.f tdetau.f \
			thermo_atm.f thermo.f trho.f trhoYL.f update.f v2_cgm.f vent.f \
			z14xcotrin21.f zoning.f journal 

clean-local:
	    rm -f *.mod
