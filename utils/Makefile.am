AM_CXXFLAGS=-Wall

BUILT_SOURCES	= star1dR.cpp star2dR.cpp read_config.cpp

AM_CPPFLAGS = -I$(top_srcdir)/src/include -I$(srcdir)/graphics
AM_CPPFLAGS +=-I$(top_srcdir)/tables/cesam -I$(srcdir)/starR
AM_CPPFLAGS +=-I$(top_srcdir)/src/graphics

EXTRA_DIST = $(srcdir)/starR/starR.h

bin_PROGRAMS		= star1dR star2dR app_luminosity

star1dR_LDADD		= $(top_builddir)/src/libester.la
star2dR_LDADD		= $(top_builddir)/src/libester.la


star1dR_SOURCES		= star1dR.cpp read_config.cpp starR/star1dR_class.cpp
star2dR_SOURCES		= star2dR.cpp read_config.cpp starR/star2dR_class.cpp
app_luminosity_SOURCES		= app_luminosity.cpp
app_luminosity_LDADD		= $(top_builddir)/src/libester.la

do_subst = sed -e 's,[@]datadir[@],$(datadir),g' \
               -e 's,[@]prefix[@],$(prefix),g'   \
               -e 's,[@]PYTHON[@],$(PYTHON),g' \
               -e 's,star1d,star1dR,g' \
               -e 's,star2d,star2dR,g' \
               -e '1a\#include \"starR.h\"'

star1dR.cpp: $(top_srcdir)/src/main/star1d.cpp Makefile
	$(do_subst) < $< > $@

star2dR.cpp: $(top_srcdir)/src/main/star2d.cpp Makefile
	$(do_subst) < $< > $@

read_config.cpp: $(top_srcdir)/src/main/read_config.cpp Makefile
	cp $< $@

CLEANFILES  = star1dR.cpp star2dR.cpp  read_config.cpp
