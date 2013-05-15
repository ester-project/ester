SRC=star1d_class.cpp star1d_solvers.cpp star_phys.cpp star1d_extra.cpp \
	star2d_class.cpp star2d_solvers.cpp star2d_extra.cpp star_map.cpp \
	convection.cpp star_evol.cpp
SRC_H=star.h
DEPS=constants.h version.h debug.h matrix.h numdiff.h mapping.h parser.h physics.h solver.h graphics.h symbolic.h

LIB=lib$(DIR).a
OBJS=$(SRC:.cpp=.o)
HEADERS=$(addprefix $(ESTER_ROOT)/include/,$(SRC_H))
DEPS:=$(addprefix $(ESTER_ROOT)/include/,$(DEPS))

.PHONY: clean distclean

$(LIB): $(HEADERS) $(OBJS)
	ar rvu $(LIB) $(OBJS)
	ranlib $(LIB)

%.o: %.cpp $(HEADERS) $(DEPS)
	$(CC) -c $< $(CCFLAGS) $(VARS)

$(ESTER_ROOT)/include/%.h: %.h $(SRC_H) $(DEPS)
	cp $< $(ESTER_ROOT)/include

clean: 
	rm -f *.o
	rm -f $(LIB)

distclean:
	rm -f $(HEADERS)
	
