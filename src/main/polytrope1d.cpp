#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "utils.h"
#include "matplotlib.h"
#include "star.h"
#include "read_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iomanip>
#include <signal.h>

int killed=0;

void sig_handler(int sig) {
    char yn;
    if (sig == SIGINT) {
        printf("\nFinish iteration and save model ([y]/n)?");
        if (scanf(" %c",&yn) == 1) {
            if (yn=='y') {
                killed = 1;
                return;
            }
        }
        else {
            killed = 1;
        }
    }
    exit(sig);
}

int main(int argc,char *argv[]) {

    int nit,last_it;
    tiempo t;

    configuration config;
    config.read_config(argc, argv);

    signal(SIGINT, sig_handler);

    t.start();


            double n = 1.5;
	    int nr=90;
            matrix h = solve_poly1d(n, 1e-10, nr, 1e-10);

	    printf("solve_poly1d finished\n");
            mapping map;
            map.set_ndomains(1);
            map.set_npts(nr);
            map.gl.set_xif(0., 1.);
            map.set_nt(1);
            map.init();
            h = map.gl.eval(h, map.r);

            matrix p = pow(h, n+1);
            matrix rho = pow(h, n);
            matrix T = p/rho;

    FILE *ficflux = fopen("polytrope.txt", "w");
    fprintf(ficflux, "index = %e\n", n);
    for(int k=0;k<nr;k++) fprintf(ficflux, "%d %e %e\n", k,map.r(k), T(k));
    fclose(ficflux);

    if(config.verbose) {
	//printf("Virial test: %e Energy test: %e\n",A.virial(),A.energy_test());
    }

    t.stop();
    if(config.verbose) printf("%2.2f seconds\n",t.value());

    return 0;
}
