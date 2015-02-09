export HDF5_CXX=$(HDF5_BACKEND)

AM_CXXFLAGS=-I. -I$(top_srcdir)/src/include
AM_LDFLAGS=-lgfortran

lib_LTLIBRARIES	= libcesam.la

library_includedir  	= $(includedir)
library_include_HEADERS = cesam.h

libcesam_la_SOURCES = mod_evol.f mod_kind.f mod_nuc.f mod_donnees.f \
					  mod_etat.f mod_opa.f mod_bp_for_alecian.f \
					  mod_numerique.f mod_variables.f \
					  interface.cpp interface_f.F90

mod_nuc.lo:				mod_donnees.lo mod_variables.lo
mod_evol.lo:			mod_kind.lo mod_nuc.lo mod_etat.lo mod_opa.lo \
						mod_bp_for_alecian.lo
mod_donnees.lo:			mod_numerique.lo
mod_variables.lo:		mod_donnees.lo
mod_etat.lo:			mod_donnees.lo
mod_opa.lo:				mod_donnees.lo
mod_numerique.lo:		mod_kind.lo
mod_bp_for_alecian.lo:	mod_donnees.lo mod_variables.lo
interface_f.lo:			mod_evol.lo

EXTRA_DIST= abon_ini.f alecian1.f arb_rom.f boite.f box.f bsp1ddn.f bsp1dn.f \
			bsp_dis.f bval0.f bval1.f bvald.f chim_gram.f coeff_rota.f coll.f \
			collision.f colpnt.f coulomb.f ctes_85.f ctes_94.f df_rotx.f \
			difdiv.f diffm_br.f diffm.f diffm_mp.f difft.f difft_nu.f \
			difft_ventura.f diffus.f diffw.f diffw_mpz.f eq_diff_chim.f \
			eq_diff_rota.f etat_ceff.f etat_eff.f etat.f etat_gong1.f \
			etat_gong2.f etat_mhd.f etat_opal5Z.f etat_opal.f etat_opalX.f \
			etat_opalZ.f evol.f fermi_dirac.f f_rad.f from_alecian.f \
			gauss_band.f genere_bases.f horner.f iben.f ini_ctes.f \
			inter.f intgauss.f kappa_cond.f linf.f lit_nl.f matinv.f \
			mod_bp_for_alecian.f mod_donnees.f mod_etat.f mod_evol.f \
			mod_kind.f mod_nuc.f mod_numerique.f mod_opa.f mod_variables.f \
			neville.f newspl.f newton.f noedif.f noein.f noeud.f \
			noeu_dis.f nuc.f opa.f opa_gong.f opa_houdek9.f opa_int_zsx.f \
			opa_opal2.f opa_opalCO.f opa_yveline.f pause.f pertw.f pertw_loc.f \
			pertw_ptm.f pertw_sch.f polyder.f pp1.f pp3.f ppcno10BeBFe.f \
			ppcno10.f ppcno10Fe.f ppcno10K.f ppcno11.f ppcno12BeBFe.f \
			ppcno12Be.f ppcno12.f ppcno12Li.f ppcno3a9.f ppcno3ac10.f \
			ppcno9.f ppcno9Fe.f print_ctes.f resout_chim.f resout_rota.f \
			rk_imps.f rq_reac.f saha.f schu58_n.f shell.f sortie.f sum_n.f \
			tabul_nuc.f taux_nuc.f v2_cgm.f vent.f z14xcotrin21.f zoning.f \
			journal

clean-local:
	    rm -f *.mod
