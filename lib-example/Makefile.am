all: polytrope1d polytrope2d

AM_LDFLAGS = -Wl,-rpath=$(prefix)/lib
LIBS += -lester -lopint -lnuc -lcesam -lnuc_cesam -lm

%:src/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(AM_LDFLAGS) -o $@ $^ $(LIBS)
