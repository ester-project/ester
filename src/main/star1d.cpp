#include "ester-config.h"
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
	int last_plot_it=0;
	configuration config(argc,argv);
	
	signal(SIGINT, sig_handler);
	t.start();
	
	star1d A;
	solver *op;
	
	if(A.init(config.input_file,config.param_file,argc,argv)) {

        ester_err("Could not initialize star");
        return 1;
                                                                  }
        //A.config.verbose=config.verbose;
	printf("verbose= %d\n",A.config.verbose);

	nit=0;
        A.time+=A.dtime;
	
	matrix tt(config.maxit+1,1),error(config.maxit+1,1);
        matrix_map error_map;
        error_map["Phi"] = zeros(config.maxit+1, 1);
        error_map["log_p"] = zeros(config.maxit+1, 1);
        error_map["log_T"] = zeros(config.maxit+1, 1);
        error_map["log_pc"] = zeros(config.maxit+1, 1);
        error_map["log_Tc"] = zeros(config.maxit+1, 1);
        error_map["Ri"] = zeros(config.maxit+1, 1);

	last_it=nit>=config.maxit; // last_it=0 normally
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
	A.n_essai=0;
	A.delta=0;  // steady solution used in solve_Xh
	A.global_err=1;
	A.glit=0;
	A.details=1;
	A.config.input_file=*config.input_file;
	printf("check config.input_file= %d\n",A.config.input_file);
    //int last_plot_it = -100;

    if (config.noplot == false) {
        plt::figure(1, 10, 4);
        A.plot(error_map.block(0, nit-1, 0 ,0));
    }

	while(!last_it) {
        try {
		if(A.global_err<0.01&&!*config.input_file) { // global_err<0.1 and no input file
			A.core_convec=core_convec_set;
			A.env_convec=env_convec_set;
		}
		nit++;
		A.glit++;
		//A.check_jacobian(op,"log_T");exit(0);
		A.global_err=A.solve(op, error_map, nit-1);
		
		tt(nit-1)=t.value();
		error(nit-1)=A.global_err;
		last_it=(A.global_err<config.tol&&nit>=config.minit)||nit>=config.maxit || killed;
		if(config.verbose) {
		  printf("it=%d err=%e\n",nit,A.global_err);
		}
                if (config.noplot == false && (nit - last_plot_it > 4 || last_it)) {
                     last_plot_it = nit;
                     A.plot(error_map.block(0, nit-1, 0 ,0));
                }
	}
        catch (runtime_exception) {
            debugger d(argc, argv, A);
            d.exec();
            return 1;
        }

  	if (nit > 200) last_it=1;
	} // End of the while loop

	if(config.verbose) {
		printf("Mass=%3.3f Msun  Radius=%3.3f Rsun  Luminosity=%3.3f Lsun  Teff=%1.1f K\n",
				A.M/M_SUN,A.R/R_SUN,A.luminosity()/L_SUN,A.Teff()(0));
		printf("X=%3.3f (Xc/X=%3.3f) Z=%3.3f\n",A.X0,A.Xc,A.Z0);
		printf("rhoc=%e Tc=%e pc=%e\n",A.rhoc,A.Tc,A.pc);
		int jc=0;
		for (int n=0;n<A.nd_core;n++) jc+=A.map.gl.npts[n];
		if(A.nd_core != 0) printf("r_cz=%3.3f Rsun\n",(A.r(jc))*A.R/R_SUN);
	}
	delete op;
	A.write(config.output_file,config.output_mode);
	//A.plot(error.block(0, nit-1, 0 ,0));
	
	t.stop();
	if(config.verbose) 
		printf("%2.2f seconds\n",t.value());	
	plt::show(true);
	return 0;
}
