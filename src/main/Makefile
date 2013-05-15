
MAIN = star1d star2d gen_output version ester_info star_evol ester
MAIN_EXES = $(addprefix $(ESTER_ROOT)/bin/,$(MAIN))

all: $(MAIN_EXES)

$(ESTER_ROOT)/bin/star2d: star2d.cpp read_config.cpp $(LIB)
	$(CC) -c star2d.cpp $(CCFLAGS) $(VARS)
	$(LD) -o $@ star2d.o $(LDFLAGS) -lester $(LIBS)

$(ESTER_ROOT)/bin/star1d: star1d.cpp read_config.cpp $(LIB)
	$(CC) -c star1d.cpp $(CCFLAGS) $(VARS)
	$(LD) -o $@ star1d.o $(LDFLAGS) -lester $(LIBS)

$(ESTER_ROOT)/bin/star_evol: star_evol.cpp read_config.cpp $(LIB)
	$(CC) -c star_evol.cpp $(CCFLAGS) $(VARS)
	$(LD) -o $@ star_evol.o $(LDFLAGS) -lester $(LIBS)

$(ESTER_ROOT)/bin/gen_output: gen_output.cpp $(LIB)
	$(CC) -c gen_output.cpp $(CCFLAGS) $(VARS)
	$(LD) -o $@ gen_output.o $(LDFLAGS) -lester $(LIBS)

$(ESTER_ROOT)/bin/version: version.cpp $(LIB)
	$(CC) -c version.cpp $(CCFLAGS) $(VARS)
	$(LD) -o $@ version.o $(LDFLAGS) -lester $(LIBS)
	
$(ESTER_ROOT)/bin/ester_info: ester_info.cpp $(LIB)
	$(CC) -c ester_info.cpp $(CCFLAGS) $(VARS)
	$(LD) -o $@ ester_info.o $(LDFLAGS) -lester $(LIBS)

$(ESTER_ROOT)/bin/ester: ester.py
	sed -e 's|ESTER_ROOT|\"$(ESTER_ROOT)\"|' ester.py > $(ESTER_ROOT)/bin/ester
	chmod +x $(ESTER_ROOT)/bin/ester

clean:
	rm -f *.o
