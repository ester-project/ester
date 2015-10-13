EXTRA_DIST = src/polytrop1d.cpp src/polytrop2d.cpp

noinst_bindir = .
noinst_bin_PROGRAMS = polytrope1d polytrope2d polytrope2d_final
polytrope1d_SOURCES = src/polytrope1d.cpp
polytrope2d_SOURCES = src/polytrope2d.cpp
polytrope2d_final_SOURCES = src/polytrope2d_final.cpp

# AM_LDFLAGS = -Wl,-rpath=$(prefix)/lib
polytrope1d_CXXFLAGS = -I$(prefix)/include
polytrope2d_CXXFLAGS = -I$(prefix)/include
polytrope2d_final_CXXFLAGS = -I$(prefix)/include

polytrope1d_LDFLAGS = -L$(prefix)/lib -lester -lopint -lnuc -lcesam -lnuc_cesam -lm -Wl,-rpath=$(prefix)/lib
polytrope2d_LDFLAGS = -L$(prefix)/lib -lester -lopint -lnuc -lcesam -lnuc_cesam -lm -Wl,-rpath=$(prefix)/lib
polytrope2d_final_LDFLAGS = -L$(prefix)/lib -lester -lopint -lnuc -lcesam -lnuc_cesam -lm -Wl,-rpath=$(prefix)/lib

# %:src/%.cpp
# 	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(AM_LDFLAGS) -o $@ $^ $(LIBS)
