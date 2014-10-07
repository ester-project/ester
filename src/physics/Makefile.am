AM_CPPFLAGS = -I$(top_srcdir)/include
AM_CPPFLAGS +=-I$(srcdir)/include
AM_CPPFLAGS +=-I../../tables/cesam_nuc
AM_CPPFLAGS	+= -DESTER_DATADIR='"$(datadir)"'

lib_LTLIBRARIES			= libnuc_cesam.la

libnuc_cesam_la_SOURCES = nuc_cesam_for.F EOS5_xtrin.F		\
                          xztrin21.F ZFS_interp_EOS5.F 		\
                          atm_onelayer.cpp composition.cpp	\
                          eos_ideal.cpp eos_idealrad.cpp	\
                          eos_opal.cpp nuc_cesam.cpp		\
                          nuc_simple.cpp opa_houdek.cpp		\
                          opa_kramer.cpp opa_opal.cpp		\
                          physics.cpp
