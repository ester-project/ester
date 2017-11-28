FFLAGS= -c -w -fdefault-real-8 -fdefault-double-8 -fd-lines-as-comments


INCLUDE	     = /home/rieutord/Ester/src/include

opa : xztrin21.o opa.o $(patsubst %.f,%.o,$(wildcard Houdek/*.f))
	g++ $^ -o opa -lgfortran

.cpp.o:
	g++ -c $(FFLAGS) -I$(INCLUDE) $*.cpp 

%.o: Houdek/%.f
	
	gfortran -c -o $@ $^


	
