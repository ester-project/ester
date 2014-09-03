UNAME:=$(shell uname -s)

FCFLAGS:=$(FCFLAGS) -I$(ESTER_ROOT)/tables/cesam_nuc
LDFLAGS:=$(LDFLAGS) -L$(ESTER_ROOT)/tables/houdek/v9/lib -L$(ESTER_ROOT)/tables/cesam_nuc
LIBS:=-lopint -lnuc $(LIBS) 

CCFLAGS:=$(CCFLAGS) -I$(ESTER_ROOT)/include 
LDFLAGS:=$(LDFLAGS) -L$(ESTER_ROOT)/lib -Wl,-rpath,$(ESTER_ROOT)/lib

F90 = $(FC)
FFLAGS90     = -O -c -w $(FCFLAGS)

AR           = ar
RANLIB       = ranlib

LIBRARY     = libnuc.a

SRCS90      = mod_kind.f mod_numerique.f mod_donnees.f sortie.f mod_variables.f mod_nuc.f

OBJS        = mod_kind.o mod_numerique.o mod_donnees.o sortie.o mod_variables.o mod_nuc.o

all:
	$(F90) $(FFLAGS90) $(SRCS90)
	$(AR) cr $(LIBRARY) $(OBJS)
	$(RANLIB) $(LIBRARY)
	rm -f $(OBJS)

clean:
	rm -f $(LIBRARY)
	rm -f *.mod

mod_nuc.o : abon_ini.f iben.f nuc.f pp1.f pp3.f ppcno10BeBFe.f ppcno10Fe.f ppcno10K.f ppcno10.f ppcno11.f ppcno12Be.f ppcno12BeBFe.f ppcno12Li.f ppcno12.f ppcno3a12Ne.f ppcno3a9.f ppcno3aco.f ppcno9.f ppcno9Fe.f rq_reac.f tabul_nuc.f taux_nuc.f
