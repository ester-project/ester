SUBDIRS	= physics

BUILT_SOURCES	= main/read_config.cpp

AM_CPPFLAGS = -I$(srcdir)/include
AM_CPPFLAGS +=-I$(top_srcdir)/tables/cesam_nuc
AM_CPPFLAGS	+= -DESTER_DATADIR='"$(datadir)"'

noinst_LTLIBRARIES		= libnuc.la libopint.la libnuc_cesam.la
libnuc_la_SOURCES		=
libnuc_la_LIBADD		= ../tables/cesam_nuc/libnuc.la

libopint_la_SOURCES		=
libopint_la_LIBADD		= ../tables/houdek/v9/lib/libopint.la

libnuc_cesam_la_SOURCES	=
libnuc_cesam_la_LIBADD	= physics/libnuc_cesam.la


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

EXTRA_DIST	= main/ester_info.cpp					\
			  $(top_srcdir)/config/1d_default.par	\
			  $(top_srcdir)/config/2d_default.par	\
			  $(top_srcdir)/config/star.cfg			\
			  $(srcdir)/include/version.h


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
library_include_HEADERS = $(srcdir)/include/constants.h \
						  $(srcdir)/include/ester.h		\
						  $(srcdir)/include/mapping.h	\
						  $(srcdir)/include/numdiff.h	\
						  $(srcdir)/include/physics.h	\
						  $(srcdir)/include/star.h		\
						  $(srcdir)/include/debug.h		\
						  $(srcdir)/include/graphics.h	\
						  $(srcdir)/include/matrix.h	\
						  $(srcdir)/include/parser.h	\
						  $(srcdir)/include/solver.h	\
						  $(srcdir)/include/symbolic.h	\
						  $(srcdir)/include/read_config.h

defaultconfigdir	= $(datadir)/ester
defaultconfig_DATA	= $(top_srcdir)/config/1d_default.par \
					  $(top_srcdir)/config/2d_default.par \
					  $(top_srcdir)/config/star.cfg


