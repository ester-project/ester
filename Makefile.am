ACLOCAL_AMFLAGS = -I m4

SUBDIRS = tables src test/models

bin_SCRIPTS	= ester-read-hdf5
CLEANFILES  = $(bin_SCRIPTS)

do_subst = sed -e 's,[@]datadir[@],$(datadir),g' \
               -e 's,[@]prefix[@],$(prefix),g'   \
               -e 's,[@]PYTHON[@],$(PYTHON),g'

ester-read-hdf5: python/read-hdf5.in Makefile
	        $(do_subst) < $< > $@

helpdir		= @datadir@/ester/doc/help/
help_DATA	= doc/help/1d	\
			  doc/help/2d	\
			  doc/help/evol \
			  doc/help/help \
			  doc/help/info \
			  doc/help/output

EXTRA_DIST	= $(top_srcdir)/python/read-hdf5.in \
			  $(top_srcdir)/doc/help/1d			\
			  $(top_srcdir)/doc/help/2d			\
			  $(top_srcdir)/doc/help/evol		\
			  $(top_srcdir)/doc/help/help		\
			  $(top_srcdir)/doc/help/info		\
			  $(top_srcdir)/doc/help/output		\
			  $(top_srcdir)/doc/doxygen.conf

doxygen:
	cd $(top_srcdir) && doxygen doc/doxygen.conf

