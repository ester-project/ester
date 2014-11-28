ACLOCAL_AMFLAGS = -I m4

SUBDIRS = tables src test/models

helpdir		= @datadir@/ester/doc/help/
help_DATA	= doc/help/1d	\
			  doc/help/2d	\
			  doc/help/evol \
			  doc/help/help \
			  doc/help/info \
			  doc/help/output

EXTRA_DIST	= $(top_srcdir)/doc/help/1d		\
			  $(top_srcdir)/doc/help/2d		\
			  $(top_srcdir)/doc/help/evol	\
			  $(top_srcdir)/doc/help/help	\
			  $(top_srcdir)/doc/help/info	\
			  $(top_srcdir)/doc/help/output	\
			  $(top_srcdir)/doc/doxygen.conf
			  $(top_srcdir)/m4

doxygen:
	cd $(top_srcdir) && doxygen doc/doxygen.conf

