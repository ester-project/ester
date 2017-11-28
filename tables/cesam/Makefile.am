lib_LTLIBRARIES	= libcesam.la

tablesdir	= $(datadir)/ester/tables/cesam
tables_DATA	= opa/opa_yveline.data

libcesam_la_SOURCES = mod_kind.f		\
					  mod_numerique.f	\
					  mod_donnees.f		\
					  sortie.f			\
					  mod_variables.f	\
					  mod_nuc.f			\
					  mod_opa.f

mod_numerique.lo	: mod_kind.lo
mod_donnees.lo		: mod_kind.lo mod_numerique.lo
sortie.lo			: mod_donnees.lo
mod_variables.lo	: mod_kind.lo mod_donnees.lo mod_numerique.lo
mod_nuc.lo			: mod_kind.lo mod_donnees.lo mod_variables.lo mod_numerique.lo
mod_opa.lo			: mod_kind.lo mod_donnees.lo mod_variables.lo mod_numerique.lo

EXTRA_DIST	= \
			  nuc/abon_ini.f nuc/iben.f nuc/nuc.f nuc/pp1.f nuc/pp3.f \
			  nuc/ppcno10BeBFe.f nuc/ppcno10.f nuc/ppcno10Fe.f \
			  nuc/ppcno10K.f nuc/ppcno11.f \
			  nuc/ppcno12BeBFe.f nuc/ppcno12Be.f nuc/ppcno12.f \
			  nuc/ppcno12Li.f nuc/ppcno3a12Ne.f nuc/ppcno3a9.f \
			  nuc/ppcno3aco.f nuc/ppcno9.f nuc/ppcno9Fe.f \
			  nuc/rq_reac.f nuc/tabul_nuc.f nuc/taux_nuc.f \
			  \
			  opa/kappa_cond.f opa/opa_compton.f opa/opa.f opa/opa_gong.f \
			  opa/opa_int_zsx.f opa/opa_opal2.f opa/opa_opalCO.f \
			  opa/opa_yveline.f opa/opa_yveline_lisse.f \
			  \
			  numeric/arb_rom.f numeric/boite.f numeric/box.f \
			  numeric/bsp1ddn.f numeric/bsp1dn.f numeric/bsp_dis.f \
			  numeric/bsp_gal.f numeric/bval0.f \
			  numeric/bval1.f numeric/bvald.f numeric/coll.f \
			  numeric/colpnt.f numeric/delete_doubles.f \
			  numeric/difdiv.f numeric/fermi_dirac.f \
			  numeric/gauss_band.f numeric/horner.f numeric/intgauss.f \
			  numeric/left_right.f numeric/linf.f numeric/matinv.f \
			  numeric/max_local.f numeric/min_max.f numeric/neville.f \
			  numeric/newspl.f numeric/newspl_gal.f numeric/newton.f \
			  numeric/noedif.f numeric/noein.f numeric/noeud.f \
			  numeric/noeu_dis.f numeric/pause.f numeric/polyder.f \
			  numeric/schu58_n.f numeric/shell.f numeric/sum_n.f \
			  numeric/zoning.f \
			  \
			  z14xcotrin21.f journal.f ini_ctes.f lit_nl.f print_ctes.f \
			  chim_gram.f inter.f sortie.f

clean-local:
	rm -f *.mod
