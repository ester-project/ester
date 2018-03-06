AM_CPPFLAGS  = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/graphics
AM_CPPFLAGS += -I../../tables/cesam
AM_CXXFLAGS=-Wall

lib_LTLIBRARIES			= libnuc_cesam.la

libnuc_cesam_la_SOURCES = nuc_cesam_for.F EOS5_xtrin.F		\
                          xztrin21.F ZFS_interp_EOS5.F 		\
                          atm_onelayer.cpp composition.cpp	\
                          eos_ideal.cpp eos_idealrad.cpp	\
                          eos_opal.cpp nuc_cesam.cpp		\
                          nuc_simple.cpp opa_houdek.cpp		\
                          opa_kramer.cpp opa_opal.cpp		\
                          physics.cpp opa_cesam.cpp 		\
						  opa_cesam_for.F90 eos_freeeos.cpp
