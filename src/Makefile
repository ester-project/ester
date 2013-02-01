
include make.inc

UNAME:=$(shell uname -s)

LDFLAGS:=$(LDFLAGS) -L$(ESTER_ROOT)/tables/houdek/v9/lib
LIBS:=-lopint $(LIBS) 

VARS=-DCBLAS=\<$(BLAS_H)\> -DLAPACK=\<$(LAPACK_H)\> -DESTER_ROOT=\"$(ESTER_ROOT)\"
ifeq ($(USE_PGPLOT),1) 
	VARS:=$(VARS) -DUSE_PGPLOT=1
endif
CCFLAGS:=$(CCFLAGS) -I$(ESTER_ROOT)/include 
LDFLAGS:=$(LDFLAGS) -L$(ESTER_ROOT)/lib -Wl,-rpath,$(ESTER_ROOT)/lib

export

DIRS = global matrix numdiff graphics mapping parser solver symbolic physics star
ifeq ($(UNAME),Darwin)
	LIB=$(ESTER_ROOT)/lib/libester.dylib
else
	LIB=$(ESTER_ROOT)/lib/libester.so
endif
ALL_LIBS = $(foreach dir,$(DIRS),$(dir)/lib$(dir).a)

CLEANDIRS = $(DIRS:%=clean-%)
DISTCLEANDIRS = $(DIRS:%=distclean-%)

.PHONY: all dirs $(DIRS) clean $(CLEANDIRS) distclean $(DISTCLEANDIRS) main

all: main $(LIB) $(ESTER_ROOT)/bin/ester_build

tables: 
	make -C ../tables	

clean: $(CLEANDIRS)
	make -C main clean

distclean: clean $(DISTCLEANDIRS)
	rm -f $(ESTER_ROOT)/bin/*
	rm -f $(LIB)
	make -C ../tables clean

main: $(LIB)
	$(MAKE) -C main

$(LIB): $(ALL_LIBS)
ifeq ($(UNAME),Darwin)
	$(LD) $(LDFLAGS) -dynamiclib -o $@ \
		$(foreach lib,$(ALL_LIBS),-force_load $(lib)) $(LIBS)
else
	$(LD) $(LDFLAGS) -shared -o $@ -Wl,--whole-archive $(ALL_LIBS) \
		-Wl,--no-whole-archive
endif

$(ALL_LIBS): $(DIRS)

$(DIRS):
	$(MAKE) -C $@ DIR=$@
	
$(CLEANDIRS):
	$(MAKE) -C $(@:clean-%=%) clean DIR=$(@:clean-%=%)

$(DISTCLEANDIRS):
	$(MAKE) -C $(@:distclean-%=%) distclean DIR=$(@:distclean-%=%)


.PHONY: ester_build

$(ESTER_ROOT)/bin/ester_build: $(LIB)
		echo '#!/bin/bash'> $@
		echo "file=\$$1" >> $@
		echo "shift" >> $@
		echo "$(CC) \$$file $(CCFLAGS) $(LDFLAGS) -lester $(LIBS)" '$(VARS)' \$$@>> $@
		chmod +x $@


global:
matrix: 
numdiff: global matrix
graphics: matrix
mapping: global matrix numdiff graphics solver
parser: matrix
solver: matrix
symbolic: matrix mapping solver
physics: global matrix numdiff
star: global matrix numdiff mapping parser physics solver graphics symbolic

dist: 
	tar -C $(ESTER_ROOT)/.. -cf ester.tar ester/src --exclude *.o --exclude *.a --exclude ester.tar.gz
	tar -C $(ESTER_ROOT)/.. -rf ester.tar ester/bin --exclude ester/bin/*
	tar -C $(ESTER_ROOT)/.. -rf ester.tar ester/config 
	tar -C $(ESTER_ROOT)/.. -rf ester.tar ester/include  --exclude ester/include/*
	tar -C $(ESTER_ROOT)/.. -rf ester.tar ester/lib --exclude ester/lib/*
	tar -C $(ESTER_ROOT)/.. -rf ester.tar ester/python
	tar -C $(ESTER_ROOT)/.. -rf ester.tar ester/tables
	tar -C $(ESTER_ROOT)/.. -rf ester.tar ester/utils
	gzip ester.tar
	
test:
	make -C ../test/models

