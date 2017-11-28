CFLAGS       = -c -w 
FFLAGS       = -c -w -ffixed-line-length-132
INCLUDE	     = /home/rieutord/Ester/src/include

test_opa: test_opa.f
mod_opa:  kappa_cond.f opa.f opa_compton.f opa_gong.f \
          opa_int_zsx.f opa_opalCO.f opa_opal2.f \
          opa_yveline.f opa_yveline_lisse.f z14xcotrin21.f

ascii2bin_opa.o: ascii2bin_opa.f
a2b: mod_kind.o mod_numerique.o ascii2bin_opa.o
	gfortran mod_numerique.o mod_kind.o ascii2bin_opa.o -o a2b -lpgplot

test_opa : mod_kind.o mod_numerique.o mod_donnees.o mod_variables.o \
           mod_opa.o test_opa.o
	gfortran mod_donnees.o  mod_numerique.o mod_variables.o \
           mod_kind.o mod_opa.o test_opa.o -o test_opa -lpgplot

.f.o:
	gfortran $(FFLAGS) $*.f
.cpp.o:
	g++ $(CFLAGS) -I$(INCLUDE) $*.cpp
