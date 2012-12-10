
MAIN = star1d star2d gen_output
MAIN_EXES = $(addprefix $(ESTER_ROOT)/bin/,$(MAIN))

all: $(MAIN_EXES)

$(ESTER_ROOT)/bin/star2d: star2d.cpp read_config.cpp $(LIB)
	$(CC) -c star2d.cpp $(CCFLAGS) $(VARS)
	$(LD) -o $@ star2d.o $(LDFLAGS) -lester $(LIBS)

$(ESTER_ROOT)/bin/star1d: star1d.cpp read_config.cpp $(LIB)
	$(CC) -c star1d.cpp $(CCFLAGS) $(VARS)
	$(LD) -o $@ star1d.o $(LDFLAGS) -lester $(LIBS)

$(ESTER_ROOT)/bin/gen_output: gen_output.cpp $(LIB)
	$(CC) -c gen_output.cpp $(CCFLAGS) $(VARS)
	$(LD) -o $@ gen_output.o $(LDFLAGS) -lester $(LIBS)

clean:
	rm -f *.o
