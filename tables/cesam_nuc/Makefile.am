export HDF5_CXX=$(HDF5_BACKEND)

lib_LTLIBRARIES	= libnuc.la

libnuc_la_SOURCES = mod_kind.f		\
					mod_numerique.f	\
					mod_donnees.f	\
					sortie.f		\
					mod_variables.f	\
					mod_nuc.f

mod_numerique.lo	: mod_kind.lo
mod_donnees.lo		: mod_kind.lo mod_numerique.lo
sortie.lo			: mod_donnees.lo
mod_variables.lo	: mod_kind.lo mod_donnees.lo mod_numerique.lo
mod_nuc.lo			: mod_kind.lo mod_donnees.lo mod_variables.lo mod_numerique.lo

mod_nuc.o: abon_ini.f iben.f nuc.f pp1.f pp3.f ppcno10BeBFe.f ppcno10Fe.f \
           ppcno10K.f ppcno10.f ppcno11.f ppcno12Be.f ppcno12BeBFe.f \
           ppcno12Li.f ppcno12.f ppcno3a12Ne.f ppcno3a9.f ppcno3aco.f \
           ppcno9.f ppcno9Fe.f rq_reac.f tabul_nuc.f taux_nuc.f

mod_donnees.o: ini_ctes.f lit_nl.f print_ctes.f

EXTRA_DIST	= arb_rom.f abon_ini.f arb_rom.f bsp1ddn.f bsp1dn.f			\
			  bsp_dis.f bsp_gal.f bval0.f bval1.f bvald.f 				\
			  chim_gram.f coll.f colpnt.f delete_doubles.f difdiv.f 	\
			  fermi_dirac.f gauss_band.f horner.f iben.f ini_ctes.f 	\
			  inter.f intgauss.f journal.f left_right.f linf.f 			\
			  lit_nl.f matinv.f max_local.f min_max.f 					\
			  neville.f newspl.f newspl_gal.f newton.f noedif.f noein.f \
			  noeud.f noeu_dis.f nuc.f pause.f polyder.f 				\
			  pp1.f pp3.f ppcno10BeBFe.f ppcno10.f ppcno10Fe.f 			\
			  ppcno10K.f ppcno11.f ppcno12BeBFe.f ppcno12Be.f ppcno12.f \
			  ppcno12Li.f ppcno3a12Ne.f ppcno3a9.f ppcno3aco.f ppcno9.f \
			  ppcno9Fe.f print_ctes.f rq_reac.f schu58_n.f shell.f 		\
			  sum_n.f tabul_nuc.f taux_nuc.f zoning.f

tablesdir		= $(datadir)/ester/tables/cesam_nuc
tables_DATA	= test.don

clean-local:
	    rm -f *.mod
