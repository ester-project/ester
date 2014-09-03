BUILT_SOURCES	= main/read_config.cpp

AM_CPPFLAGS = -I$(top_srcdir)/include
AM_CPPFLAGS +=-I$(srcdir)/include
AM_CPPFLAGS +=-I../tables/cesam_nuc
AM_CPPFLAGS	+= -DESTER_DATADIR='"$(datadir)"'

noinst_LTLIBRARIES	= libnuc.la libopint.la
libnuc_la_SOURCES	=
libnuc_la_LIBADD	= ../tables/cesam_nuc/libnuc.la

libopint_la_SOURCES	=
libopint_la_LIBADD	= ../tables/houdek/v9/libopint.la


bin_PROGRAMS		= star1d star2d gen_output version star_evol

star1d_LDFLAGS		= -lester -lnuc_cesam -lnuc -lopint
star2d_LDFLAGS		= -lester -lnuc_cesam -lnuc -lopint
gen_output_LDFLAGS	= -lester -lnuc_cesam -lnuc -lopint
version_LDFLAGS		= -lester -lnuc_cesam -lnuc -lopint
star_evol_LDFLAGS	= -lester -lnuc_cesam -lnuc -lopint


star1d_SOURCES		= main/star1d.cpp main/read_config.cpp
star2d_SOURCES		= main/star2d.cpp main/read_config.cpp
gen_output_SOURCES	= main/gen_output.cpp
version_SOURCES		= main/version.cpp
star_evol_SOURCES	= main/star_evol.cpp main/read_config.cpp

star1d		: libnuc.la libopint.la libester.la
star2d		: libnuc.la libopint.la libester.la
gen_output	: libnuc.la libopint.la libester.la
version 	: libnuc.la libopint.la libester.la
star_evol	: libnuc.la libopint.la libester.la

EXTRA_DIST			= main/ester_info.cpp					\
					  include/cblas.h include/lapack.h		\
					  $(top_srcdir)/config/1d_default.par	\
					  $(top_srcdir)/config/2d_default.par	\
					  $(top_srcdir)/config/star.cfg


bin_SCRIPTS	= main/ester

lib_LTLIBRARIES			= libester.la

libester_la_SOURCES		= global/global.cpp									\
                          graphics/graphics.cpp								\
                          mapping/mapping.cpp mapping/remapper.cpp 			\
                          matrix/mat_math.cpp matrix/matrix_block_diag.cpp	\
                          matrix/matrix.cpp matrix/matrix_map.cpp			\
                          matrix/mat_spec.cpp 								\
                          numdiff/diff_gl.cpp numdiff/diff_leg.cpp 			\
                          parser/iofile.cpp parser/parser.cpp 				\
                          physics/atm_onelayer.cpp physics/composition.cpp	\
                          physics/eos_ideal.cpp physics/eos_idealrad.cpp	\
                          physics/eos_opal.cpp physics/nuc_cesam.cpp		\
                          physics/nuc_simple.cpp physics/opa_houdek.cpp		\
                          physics/opa_kramer.cpp physics/opa_opal.cpp		\
                          physics/physics.cpp 								\
                          physics/nuc_cesam_for.F physics/EOS5_xtrin.F		\
                          physics/xztrin21.F physics/ZFS_interp_EOS5.F		\
                          solver/RKF_solver.cpp solver/SDIRK_solver.cpp		\
                          solver/solver.cpp solver/solver_full.cpp 			\
                          star/convection.cpp star/star1d_class.cpp			\
                          star/star1d_extra.cpp star/star1d_solvers.cpp 	\
                          star/star2d_class.cpp star/star2d_extra.cpp		\
                          star/star2d_solvers.cpp star/star_evol.cpp		\
                          star/star_map.cpp star/star_phys.cpp 				\
                          symbolic/rational.cpp symbolic/symbolic.cpp		\
                          symbolic/sym.cpp symbolic/sym_expr.cpp			\
                          symbolic/sym_tens.cpp symbolic/sym_vec.cpp

library_includedir  	= $(includedir)/
library_include_HEADERS = $(top_srcdir)/include/constants.h \
						  $(top_srcdir)/include/ester.h		\
						  $(top_srcdir)/include/mapping.h	\
						  $(top_srcdir)/include/numdiff.h	\
						  $(top_srcdir)/include/physics.h	\
						  $(top_srcdir)/include/star.h		\
						  $(top_srcdir)/include/version.h	\
						  $(top_srcdir)/include/debug.h		\
						  $(top_srcdir)/include/graphics.h	\
						  $(top_srcdir)/include/matrix.h	\
						  $(top_srcdir)/include/parser.h	\
						  $(top_srcdir)/include/solver.h	\
						  $(top_srcdir)/include/symbolic.h	\
						  $(top_srcdir)/include/read_config.h

defaultconfigdir	= $(datadir)/ester
defaultconfig_DATA	= $(top_srcdir)/config/1d_default.par \
					  $(top_srcdir)/config/2d_default.par \
					  $(top_srcdir)/config/star.cfg


