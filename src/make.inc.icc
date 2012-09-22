ESTER_ROOT=$(CURDIR)/..

USE_PGPLOT=1

BLAS_LIB=-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
CBLAS_LIB=
LAPACK_LIB=

CC=icpc
CCFLAGS=-O2 -fp-model precise -ip -fPIC  -I/usr/include/x86_64-linux-gnu
FC=ifort
FCFLAGS=-O2 -r8 -nowarn -fp-model precise -fpp -fPIC
LD=icpc
LDFLAGS=
LIBS=$(LAPACK_LIB) $(CBLAS_LIB) $(BLAS_LIB) -liomp5 -lpthread -lifcore -lifport -lgfortran

ifeq ($(USE_PGPLOT),1) 
	LIBS:=-lcpgplot -lpgplot -lX11 -lpng12 $(LIBS)
endif

BLAS_H=$(ESTER_ROOT)/src/include/cblas.h
LAPACK_H=$(ESTER_ROOT)/src/include/lapack.h

