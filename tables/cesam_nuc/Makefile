.IGNORE:

PROGRAM        	= test_tabul_reac

#----------------------------------------------------------------------------
# Next lines are automatically produced by configure.
# Edit Makefile.inc then launch configure
# configure customize here.  DO NOT CHANGE THIS LINE ----------------
AR           = ar
CLEANOBJECTS = *.o *.mod $(DIR_MOD)/*.mod
CPP          = /lib/cpp -P
DEST         = ..
DIR_MOD      = ./lib/$(OS)
DIR_MOD_LIB  = ../../lib/$(OS)
F77          = ifort
F90          = ifort
FFLAGS       = -O -c -w -traceback -extend-source
FFLAGS90     = -O -c -w -traceback -real-size 64 -extend-source
FLAG_MOD     = -I$(DIR_MOD)
FLAG_MOD_LIB = -I$(DIR_MOD_LIB)
INSTALL      = mv
LD           = ifort
LDFLAGS      = -Bdynamic
LIBPGPLOT    = /usr/lib/libpgplot.so
LIBBLAS      = 
LIBLAPACK    =
OS           = Linux
RANLIB       = ranlib
SIZE         = size
# end of configure customizations. DO NOT CHANGE THIS LINE ----------
#----------------------------------------------------------------------------
# Libraries

MYLIBDIRPREFIX = ./lib
MYLIBDIR       = $(MYLIBDIRPREFIX)/$(OS)
LIBNUMERIC     = $(MYLIBDIR)/libnumeric.a

LOWLIBS      =  $(LIBNUMERIC) 

MYLIBS       =  $(LOWLIBS)
LIBS         =  $(LIBPGPLOT)

#----------------------------------------------------------------------------
OBJS	=   	mod_donnees.o sortie.o mod_variables.o mod_nuc.o


#----------------------------------------------------------------------------
MAKE	= 	make MACHINE="$(MACHINE)" RANLIB="$(RANLIB)" \
		OS="$(OS)" F77="$(F77)" FFLAGS="$(FFLAGS)" \
		F90="$(F90)" FFLAGS90="$(FFLAGS90)" AR="$(AR)" \
		INSTALL="$(INSTALL)" CPP="$(CPP)" FLAG_MOD_LIB="$(FLAG_MOD_LIB)"

MAKECLEAN = 	make clean

#----------------------------------------------------------------------------
.SUFFIXES: .f90 .F

.F.f90:
		$(CPP) -D$(OS) $*.F >$*.f90
.F.o:
		$(CPP) -D$(OS) $*.F >$*.f90
		$(F90) $(FFLAGS90) $(FLAG_MOD) $*.f90

.f.o:
		$(F90) $(FFLAGS90) $(FLAG_MOD) $*.f
.f90.o:
		$(F90) $(FFLAGS90) $(FLAG_MOD) $*.f90

#----------------------------------------------------------------------------
all:		$(PROGRAM) 

$(PROGRAM):     low_libs $(OBJS) $(PROGRAM).o
		$(LD) $(LDFLAGS) $(PROGRAM).o $(OBJS)  \
		$(MYLIBS) $(LIBS) -o $(PROGRAM)
		$(SIZE) $(PROGRAM)

install:	$(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@$(INSTALL) $(PROGRAM) $(DEST)

low_libs:
	@test -d $(MYLIBDIR) || mkdir $(MYLIBDIR) ;\
	test $(MYLIBDIRPREFIX) = "./lib" && MYLIBDEST="../$(OS)" ;\
	test $(MYLIBDIRPREFIX) = "./lib" || MYLIBDEST=$(MYLIBDIR) ;\
	for lib in $(LOWLIBS) ; do \
	DIR=`echo "$$lib"|sed 's/[^\/]*\/\([a-zA-Z0-9_]*\)\.[^.]*$$/\1/'`;\
	LIB=`basename "$$lib"`;\
	(cd $$DIR ; $(MAKE) DEST="$$MYLIBDEST" $$LIB);\
	done

clean:
	@rm -f $(PROGRAM) degrill test_format.f90 $(MYLIBS) $(DIR_MOD)/*.*
	@rm -rf $(CLEANOBJECTS)
	for lib in $(MYLIBS) ; do \
	DIR=`echo "$$lib"|sed 's/[^\/]*\/\([a-zA-Z0-9_]*\)\.[^.]*$$/\1/'`;\
	LIB=`basename "$$lib"`;\
	(cd $$DIR ; $(MAKECLEAN) );\
	done

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# dependencies:
# DO NOT DELETE THIS LINE - used by make depend
mod_nuc.o : abon_ini.f iben.f nuc.f pp1.f pp3.f ppcno10BeBFe.f ppcno10Fe.f ppcno10K.f ppcno10.f ppcno11.f ppcno12Be.f ppcno12BeBFe.f ppcno12Li.f ppcno12.f ppcno3a12Ne.f ppcno3a9.f ppcno3aco.f ppcno9.f ppcno9Fe.f rq_reac.f tabul_nuc.f taux_nuc.f

