tables: opal houdek cesam_nuc

opal:
	tar zxf opal/eos_tables.tar.gz --directory opal

houdek:
	make FC=$(FC) -C houdek/v9
	echo "$(CURDIR)/houdek/v9/opal95/opal95e.bin">$(CURDIR)/houdek/v9/OPINTPATH_AX
	echo "$(CURDIR)/houdek/v9/alex94/alex94.bin">>$(CURDIR)/houdek/v9/OPINTPATH_AX
	echo "$(CURDIR)/houdek/v9/opal95/opal95e+alex94-pd.bin">>$(CURDIR)/houdek/v9/OPINTPATH_AX
	echo "$(CURDIR)/houdek/v9/opal95/ival95.dat">>$(CURDIR)/houdek/v9/OPINTPATH_AX

cesam_nuc:
	make FC=$(FC) FCFLAGS="$(FCFLAGS)" -C cesam_nuc

clean: opal_clean houdek_clean cesam_nuc_clean

opal_clean:
	rm -f -r opal/eos_tables
	
houdek_clean:
	make -C houdek/v9 clean
	
cesam_nuc_clean:
	make -C cesam_nuc clean
	
.PHONY: tables opal houdek clean opal_clean houdek_clean cesam_nuc cesam_nuc_clean
