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

if USE_HDF5
hdf5=1
else
hdf5=0
endif
do_subst = sed -e 's,[@]datadir[@],$(datadir),g'	\
		   -e 's,[@]prefix[@],$(prefix),g'			\
		   -e 's,[@]use_hdf5[@],$(hdf5),g'			\
		   -e 's,[@]PYTHON[@],$(PYTHON),g'

bin_SCRIPTS = ester-check-models
CLEANFILES = $(bin_SCRIPTS)

ester-check-models: $(top_srcdir)/test/models/test_models.in Makefile
	$(do_subst) < $< > $@
	chmod +x $@

checksdir	= $(datadir)/ester/test/models
checks_DATA	= 1d.par 2d.par out1 out2 out3 out4 star.out template_1d template_2d

