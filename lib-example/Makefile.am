EXTRA_DIST = src/poly.h

AM_CXXFLAGS=-I$(top_srcdir)/src/include

AM_LIBS = -lm

noinst_bindir = $(abs_builddir)/bin
noinst_bin_PROGRAMS = polytrope1d polytrope2d polytrope2d_final
polytrope1d_SOURCES = src/polytrope1d.cpp src/polytrope.cpp
polytrope1d_LDADD = $(top_builddir)/src/libester.la

polytrope2d_SOURCES = src/polytrope2d.cpp src/polytrope.cpp
polytrope2d_LDADD = $(top_builddir)/src/libester.la

polytrope2d_final_SOURCES = src/polytrope2d_final.cpp src/polytrope.cpp
polytrope2d_final_LDADD = $(top_builddir)/src/libester.la
