SRC=eos_ideal.cpp eos_idealrad.cpp eos_opal.cpp \
	nuc_simple.cpp nuc_cesam.cpp \
	opa_houdek.cpp opa_kramer.cpp opa_opal.cpp \
	atm_onelayer.cpp \
	composition.cpp \
	physics.cpp
SRC_F=EOS5_xtrin.F xztrin21.F ZFS_interp_EOS5.F nuc_cesam_for.F
SRC_H=physics.h
DEPS=constants.h matrix.h numdiff.h

LIB=lib$(DIR).a
OBJS=$(SRC:.cpp=.o)
OBJS:=$(OBJS) $(SRC_F:.F=.o)
HEADERS=$(addprefix $(ESTER_ROOT)/include/,$(SRC_H))
DEPS:=$(addprefix $(ESTER_ROOT)/include/,$(DEPS))

.PHONY: clean distclean

$(LIB): $(HEADERS) $(OBJS)
	ar rvu $(LIB) $(OBJS)
	ranlib $(LIB)

%.o: %.cpp $(HEADERS) $(DEPS)
	$(CC) -c $< $(CCFLAGS) $(VARS)
	
%.o: %.F 
	$(FC) -c $< $(FCFLAGS) $(VARS)

$(ESTER_ROOT)/include/%.h: %.h $(SRC_H) $(DEPS)
	cp $< $(ESTER_ROOT)/include

clean: 
	rm -f *.o
	rm -f $(LIB)

distclean:
	rm -f $(HEADERS)
	
