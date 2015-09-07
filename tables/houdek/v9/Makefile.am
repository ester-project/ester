EXTRA_DIST = opal95/opal95e.bin alex94/alex94.bin \
			 opal95/opal95e+alex94-pd.bin \
			 opal95/ival95.dat

SUBDIRS = lib

houdekdir	= $(datadir)/ester/tables/houdek/v9
houdek_DATA = OPINTPATH_AX \
			  opal95/opal95e.bin \
			  alex94/alex94.bin \
			  opal95/opal95e+alex94-pd.bin \
			  opal95/ival95.dat

OPINTPATH_AX:
	echo "$(datadir)/ester/tables/houdek/v9/opal95e.bin"			>  $@
	echo "$(datadir)/ester/tables/houdek/v9/alex94.bin" 			>> $@
	echo "$(datadir)/ester/tables/houdek/v9/opal95e+alex94-pd.bin"	>> $@
	echo "$(datadir)/ester/tables/houdek/v9/ival95.dat" 			>> $@

clean-local:
	rm -f OPINTPATH_AX
