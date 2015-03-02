SUBDIRS	= cesam physics
AM_CXXFLAGS=-Wall

export HDF5_CXX=$(HDF5_BACKEND)

BUILT_SOURCES	= main/read_config.cpp

AM_CPPFLAGS = -I$(srcdir)/include -I$(top_srcdir)/cesam
AM_CPPFLAGS +=-I$(top_srcdir)/tables/cesam_nuc -I$(srcdir)/physics

noinst_LTLIBRARIES		= libopint.la libcesam.la

libopint_la_SOURCES		= dummy.cpp
libopint_la_LIBADD		= ../tables/houdek/v9/lib/libopint.la

libcesam_la_SOURCES	= dummy.cpp
libcesam_la_LIBADD	= cesam/libcesam.la

bin_PROGRAMS		= star1d star2d gen_output version star_evol ester_info

star1d_LDFLAGS		= -lester -lopint -lcesam
star2d_LDFLAGS		= -lester -lopint -lcesam
gen_output_LDFLAGS	= -lester -lopint -lcesam
version_LDFLAGS		= -lester -lopint -lcesam
star_evol_LDFLAGS	= -lester -lopint -lcesam
ester_info_LDFLAGS	= -lester -lopint -lcesam

star1d_SOURCES		= main/star1d.cpp main/read_config.cpp
star2d_SOURCES		= main/star2d.cpp main/read_config.cpp
star_evol_SOURCES	= main/star_evol.cpp main/read_config.cpp
gen_output_SOURCES	= main/gen_output.cpp
version_SOURCES		= main/version.cpp
ester_info_SOURCES	= main/ester-info.cpp

star1d		: libopint.la libester.la
star2d		: libopint.la libester.la
gen_output	: libopint.la libester.la
version 	: libopint.la libester.la
star_evol	: libopint.la libester.la
ester_info	: libopint.la libester.la

EXTRA_DIST	= $(srcdir)/main/ester.in				\
			  $(srcdir)/include/lapack.h			\
			  $(srcdir)/include/utils.h


bin_SCRIPTS	= main/ester
CLEANFILES  = $(bin_SCRIPTS)

do_subst = sed -e 's,[@]datadir[@],$(datadir),g' \
               -e 's,[@]prefix[@],$(prefix),g'   \
               -e 's,[@]PYTHON[@],$(PYTHON),g'

main/ester: main/ester.in Makefile
	$(do_subst) < $< > $@

lib_LTLIBRARIES			= libester.la
libester_la_SOURCES		= global/global.cpp									\
                          graphics/graphics.cpp								\
                          mapping/mapping.cpp mapping/remapper.cpp 			\
                          matrix/mat_math.cpp matrix/matrix_block_diag.cpp	\
                          matrix/matrix.cpp matrix/matrix_map.cpp			\
                          matrix/mat_spec.cpp 								\
                          numdiff/diff_gl.cpp numdiff/diff_leg.cpp 			\
                          utils/iofile.cpp utils/parser.cpp 				\
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
libester_la_LIBADD	= physics/libester-physics.la

library_includedir  	= $(includedir)/
library_include_HEADERS = $(srcdir)/include/constants.h \
						  $(srcdir)/include/ester.h		\
						  $(srcdir)/include/mapping.h	\
						  $(srcdir)/include/numdiff.h	\
						  $(srcdir)/include/physics.h	\
						  $(srcdir)/include/star.h		\
						  $(srcdir)/include/graphics.h	\
						  $(srcdir)/include/matrix.h	\
						  $(srcdir)/include/parser.h	\
						  $(srcdir)/include/solver.h	\
						  $(srcdir)/include/symbolic.h	\
						  $(srcdir)/include/read_config.h

