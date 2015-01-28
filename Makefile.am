ACLOCAL_AMFLAGS = -I m4

export HDF5_CXX=$(HDF5_BACKEND)

SUBDIRS = . tables src test/models

all-local: ester-config.h
clean-local:
	rm -f ester-config.h

bin_SCRIPTS	= ester_visu.py
CLEANFILES  = $(bin_SCRIPTS)

do_subst = sed -e 's,[@]datadir[@],$(datadir),g' \
               -e 's,[@]prefix[@],$(prefix),g'   \
               -e 's,[@]PYTHON[@],$(PYTHON),g'

ester_visu.py: python/ester_visu.in Makefile
	$(do_subst) < $< > $@

ester-config.h: config.h Makefile
	$(do_subst) < $< > $@

helpdir		= @datadir@/ester/doc/help/
help_DATA	= doc/help/1d	\
			  doc/help/2d	\
			  doc/help/evol \
			  doc/help/help \
			  doc/help/info \
			  doc/help/output

EXTRA_DIST = $(top_srcdir)/README \
			 $(top_srcdir)/NEWS \
			 $(top_srcdir)/INSTALL \
			 $(top_srcdir)/LICENSE \
			 $(top_srcdir)/python/ester_visu.in \
			 $(top_srcdir)/doc/help/1d \
			 $(top_srcdir)/doc/help/2d \
			 $(top_srcdir)/doc/help/evol \
			 $(top_srcdir)/doc/help/help \
			 $(top_srcdir)/doc/help/info \
			 $(top_srcdir)/doc/help/output \
			 $(top_srcdir)/doc/doxygen.conf \
			 $(top_srcdir)/m4 \
			 $(top_srcdir)/autoconf/ester_get_version.sh

doxygen:
	cd $(top_srcdir) && doxygen doc/doxygen.conf

