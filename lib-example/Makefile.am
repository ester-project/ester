EXTRA_DIST = src/polytrop1d.cpp src/polytrop2d.cpp

AM_CXXFLAGS=-I$(prefix)/include -Wall
AM_LDFLAGS = -L$(prefix)/lib -Wl,-rpath=$(prefix)/lib \
			 -lester -lnuc_cesam -lnuc -lopint -lm

noinst_bindir = .
noinst_bin_PROGRAMS = polytrope1d polytrope2d polytrope2d_final
polytrope1d_SOURCES = src/polytrope1d.cpp src/polytrope.cpp
polytrope2d_SOURCES = src/polytrope2d.cpp src/polytrope.cpp
polytrope2d_final_SOURCES = src/polytrope2d_final.cpp src/polytrope.cpp

