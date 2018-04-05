FFLAGS= -w -fdefault-real-8 -fdefault-double-8 -fd-lines-as-comments


INCLUDE	     = /home/rieutord/Ester/src/include
LIB	     = /home/rieutord/Ester/local/lib/

scan : scan.o
	g++ $^ -L$(LIB) -o scan -lphysics -lester -lfreeeos

.cpp.o:
	g++ -c $(FFLAGS) -I$(INCLUDE) $*.cpp 
