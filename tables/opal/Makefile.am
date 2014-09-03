EXTRA_DIST = eos_tables.tar.gz GN93hz

tablesdir		= $(datadir)/ester/tables/opal
eostablesdir	= $(datadir)/ester/tables/opal/eos_tables
tables_DATA	= GN93hz
eostables_DATA	= eos_tables/EOS5_00z0x		\
				  eos_tables/EOS5_00z10x	\
				  eos_tables/EOS5_00z2x		\
				  eos_tables/EOS5_00z4x		\
				  eos_tables/EOS5_00z6x		\
				  eos_tables/EOS5_00z8x		\
				  eos_tables/EOS5_02z0x		\
				  eos_tables/EOS5_02z2x		\
				  eos_tables/EOS5_02z4x		\
				  eos_tables/EOS5_02z6x		\
				  eos_tables/EOS5_02z8x		\
				  eos_tables/EOS5_04z0x		\
				  eos_tables/EOS5_04z2x		\
				  eos_tables/EOS5_04z4x		\
				  eos_tables/EOS5_04z6x		\
				  eos_tables/EOS5_04z8x

eos_tables/EOS5_00z0x	: eos_tables
eos_tables/EOS5_00z10x	: eos_tables
eos_tables/EOS5_00z2x 	: eos_tables
eos_tables/EOS5_00z4x 	: eos_tables
eos_tables/EOS5_00z6x 	: eos_tables
eos_tables/EOS5_00z8x 	: eos_tables
eos_tables/EOS5_02z0x 	: eos_tables
eos_tables/EOS5_02z2x 	: eos_tables
eos_tables/EOS5_02z4x 	: eos_tables
eos_tables/EOS5_02z6x 	: eos_tables
eos_tables/EOS5_02z8x 	: eos_tables
eos_tables/EOS5_04z0x 	: eos_tables
eos_tables/EOS5_04z2x 	: eos_tables
eos_tables/EOS5_04z4x 	: eos_tables
eos_tables/EOS5_04z6x 	: eos_tables
eos_tables/EOS5_04z8x 	: eos_tables

eos_tables:
	tar zxf $(top_srcdir)/tables/opal/eos_tables.tar.gz

clean-local:
	    rm -rf eos_tables

