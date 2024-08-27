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
    double err;
    tiempo t;
    // double t_plot;

    configuration config(argc, argv);

    signal(SIGINT, sig_handler);

    t.start();

    star1d A;
    solver *op;

    if(A.init(config.input_file,config.param_file,argc,argv)) {
        ester_critical("Could not initialize star");
        return 1;
    }

    nit=0;

    matrix tt(config.maxit+1,1),error(config.maxit+1,1);
    matrix_map error_map;
    error_map["Phi"] = zeros(config.maxit+1, 1);
    error_map["log_p"] = zeros(config.maxit+1, 1);
    error_map["log_T"] = zeros(config.maxit+1, 1);
    error_map["log_pc"] = zeros(config.maxit+1, 1);
    error_map["log_Tc"] = zeros(config.maxit+1, 1);
    error_map["Ri"] = zeros(config.maxit+1, 1);

    // t_plot=0;
    last_it = (nit>=config.maxit) || killed;
    op=A.init_solver();
    if(config.verbose>2) op->verbose=1;
    A.config.newton_dmax=config.newton_dmax;
    if(config.verbose>1) A.config.verbose=1;

    // If no input file, ignore core convection until the model starts to converge
    int core_convec_set=A.core_convec;
    int env_convec_set=A.env_convec;
    if(*config.input_file==0) {
        A.core_convec=0;
        A.env_convec=0;
    }
    err=1;

    int last_plot_it = -100;

    if (config.noplot == false) {
        plt::init();
        plt::figure(1, 10, 4);
        A.plot(error_map.block(0, nit-1, 0 ,0));
        plt::show();
    }

    while(!last_it) {

        try {
            if(err<0.1&&!*config.input_file) {
                A.core_convec=core_convec_set;
                A.env_convec=env_convec_set;
            }
            nit++;
            //A.check_jacobian(op,"log_T");exit(0);

            err=A.solve(op, error_map, nit-1);

            tt(nit-1)=t.value();
            error(nit-1)=err;
            last_it=(err<config.tol&&nit>=config.minit)||nit>=config.maxit || killed;
            if(config.verbose) {
                printf("it=%d err=%e\n",nit,err);
            }
            if (config.noplot == false && (nit+1 - last_plot_it > config.plot_interval || last_it)) {
                last_plot_it = nit;
                A.plot(error_map.block(0, nit-1, 0 ,0));
            }
        }
        catch (runtime_exception) {
            debugger d(argc, argv, A);
            d.exec();
            return 1;
        }

    }

    if(config.verbose) {
        printf("Mass=%3.3f Msun  Radius=%3.3f Rsun  Luminosity=%3.3f Lsun  Teff=%1.1f K\n",
                A.M/M_SUN,A.R/R_SUN,A.luminosity()/L_SUN,A.Teff()(0));
        printf("X=%3.3f (Xc/X=%3.3f) Z=%3.3f\n",A.X0,A.Xc,A.Z0);
        printf("rhoc=%e Tc=%e pc=%e\n",A.rhoc,A.Tc,A.pc);
        A.test_virial=A.virial(); A.test_energy=A.energy_test();
	printf("Virial test: %e Energy test: %e\n",A.test_virial,A.test_energy);

        if(A.conv) printf("r_cz=%3.3f Rsun\n",*(A.map.gl.xif+A.conv)*A.R/R_SUN);
    }
    delete op;
    A.write(config.output_file);

    // if(config.verbose) {
    //     delete fig;
    // }

    t.stop();
    if(config.verbose)
        printf("%2.2f seconds\n",t.value());

    plt::show(true);

    return 0;
}
