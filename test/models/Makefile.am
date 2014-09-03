EXTRA_DIST = 1d.par 		\
			 2d.par 		\
			 out1			\
			 out2			\
			 out3			\
			 out4			\
			 star.out		\
			 template_1d	\
			 template_2d	\
			 $(top_srcdir)/test/models/test_models.in

do_subst = sed -e 's,[@]datadir[@],$(datadir),g'	\
		   -e 's,[@]prefix[@],$(prefix),g'			\
		   -e 's,[@]PYTHON[@],$(PYTHON),g'

bin_SCRIPTS = test_models
CLEANFILES = $(bin_SCRIPTS)

test_models: $(top_srcdir)/test/models/test_models.in Makefile
	$(do_subst) < $< > $@
	chmod +x $@

checksdir	= $(datadir)/ester/test/models
checks_DATA	= 1d.par 2d.par out1 out2 out3 out4 star.out template_1d template_2d

