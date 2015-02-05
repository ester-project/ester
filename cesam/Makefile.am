export HDF5_CXX=$(HDF5_BACKEND)

lib_LTLIBRARIES	= libcesam.la

EXTRA_DIST = rk_imps.f abon_ini.f

libcesam_la_SOURCES = mod_evol.f mod_kind.f mod_nuc.f mod_donnees.f \
					  interface.cpp

mod_nuc.lo:			mod_donnees.lo mod_variables.lo
mod_evol.lo:		mod_kind.lo mod_nuc.lo
mod_donnees.lo:		mod_numerique.lo
mod_variables.lo:	mod_donnees.lo

clean-local:
	    rm -f *.mod
