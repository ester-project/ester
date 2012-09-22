ESTER_ROOT=$(CURDIR)/..

USE_PGPLOT=1

BLAS_LIB=$(HOME)/usr/atlas/lib/libptf77blas.a $(HOME)/usr/atlas/lib/libatlas.a
CBLAS_LIB=$(HOME)/usr/atlas/lib/libptcblas.a
LAPACK_LIB=$(HOME)/usr/atlas/lib/liblapack.a

CC=g++
CCFLAGS=-O3 -fPIC -pthread
#CCFLAGS=-O0 -fPIC -g  
FC=gfortran
FCFLAGS=-O3 -fPIC -fno-align-commons -fdefault-real-8 -pthread
#FCFLAGS=-O0 -fPIC -fno-align-commons -fdefault-real-8 -g
LD=g++
LDFLAGS=-pthread
#LDFLAGS=-pg
LIBS=$(LAPACK_LIB) $(CBLAS_LIB) $(BLAS_LIB) -lgfortran

ifeq ($(USE_PGPLOT),1) 
	LIBS:=-lcpgplot -lpgplot -lX11 -lpng12 $(LIBS)
endif

BLAS_H=$(ESTER_ROOT)/src/include/cblas.h
LAPACK_H=$(ESTER_ROOT)/src/include/lapack.h

