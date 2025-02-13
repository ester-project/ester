#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "utils.h"
#include "star.h"
#include "read_config.h"
#include "matplotlib.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <signal.h>

int killed=0;

void sig_handler(int sig);

int main(int argc,char *argv[]) {

	int nit,last_it;
	double err=1;
	double t_plot;
	tiempo t;
	// figure *fig = NULL;

    // Create config object
    configuration config;
    // Parse configuration file star.cfg and command line simulation parameters (not star ones)
	config.read_config(argc, argv);

	signal(SIGINT,sig_handler);

    plt::figure(1, 10, 4);
	
	t.start();
		
	// if(config.verbose) {
	// 	fig=new figure(config.plot_device);
	// 	fig->subplot(2,2);
	// }
	
	star2d A;
	solver *op;
	
	if(A.init(config.input_file, config.param_file, argc, argv)) {
        ester_err("Could not initialize star");
        return 1;
    }
	
	nit=0;
	
	matrix tt(config.maxit+1,1),error(config.maxit+1,1);
    matrix_map error_map;
    error_map["Phi"] = zeros(config.maxit+1, 1);
    error_map["p"] = zeros(config.maxit+1, 1);
    error_map["T"] = zeros(config.maxit+1, 1);
    error_map["log_pc"] = zeros(config.maxit+1, 1);
    error_map["log_Tc"] = zeros(config.maxit+1, 1);
    error_map["Ri"] = zeros(config.maxit+1, 1);

	t_plot=0;
	last_it=nit>=config.maxit;
	op=A.init_solver();
	if(config.verbose>2) op->verbose=1;
	A.config.newton_dmax=config.newton_dmax;
	if(config.verbose>1) A.config.verbose=1;
	
	// If no input file, ignore core convection until the model starts to converge
	int core_convec_set=A.core_convec;
	if(*config.input_file==0) {
		A.core_convec=0;
	}
	
    while(!last_it) {
        if (A.config.dump_iter) {
            char *filename = NULL;
            if (asprintf(&filename, "%s-iter%d.hdf5",
                        config.output_file,
                        nit) != -1) {
                A.hdf5_write(filename);
                free(filename);
            }
        }

        nit++;

		if(err<0.1&&!*config.input_file) {
            A.core_convec=core_convec_set;
		}
		
		err=A.solve(op, error_map, nit-1);

		tt(nit-1)=t.value();
		error(nit-1)=err;
		last_it=(err<config.tol&&nit>=config.minit)||nit>=config.maxit;
		if(killed) last_it=1;
		if(config.verbose) {
			printf("it=%d err=%e (%2.2fs)\n",nit,err,t.value());
			printf("\tOmega=%e (%2.2f%%) eps=%.2e M=%f\n",A.Omega,A.Omega/A.Omegac*100,1-1./A.map.leg.eval_00(A.r.row(-1),PI/2)(0),A.m*A.rhoc*A.R*A.R*A.R/M_SUN);

			if(tt(nit-1)-t_plot>config.plot_interval||last_it) {
#if 0
				fig->semilogy(error.block(0,nit-1,0,0));
				fig->label("Iteration number","Relative error","");
				fig->colorbar();
				A.spectrum(fig,A.rho);
				fig->label("Density (normalized spectrum)","","");
				fig->colorbar();
				A.drawi(fig,A.w,100,64);
				fig->label("Differential rotation","","");
				A.drawci(fig,A.G,100,64,15,11);
				fig->label("Meridional circulation","","");
				t_plot=tt(nit-1);
#else
                static bool plotOpen = false;
                if (plotOpen == false) {
                    plt::show();
                    plotOpen = true;
                }
                A.plot(error_map.block(0, nit-1, 0, 0));
#endif
			}

		}
		
	}

	if(config.verbose) {
		printf("\nMass=%.3f Msun  Luminosity=%.3f Lsun\n",
				A.M/M_SUN,A.luminosity()/L_SUN);
		double Re=A.units.r*A.map.leg.eval_00(A.r.row(-1),PI/2)(0);
		printf("Radius(p)=%.3f Rsun    Radius(e)=%.3f Rsun\n",
			A.R/R_SUN,Re/R_SUN);
		printf("Teff(p)=%.2f K      Teff(e)=%.2f K\n",
			A.map.leg.eval_00(A.Teff(),0)(0),A.map.leg.eval_00(A.Teff(),PI/2)(0));
		printf("log(geff)(p)=%.2f       log(geff)(e)=%.2f\n",
			A.map.leg.eval_00(log10(A.gsup()),0)(0),A.map.leg.eval_00(log10(A.gsup()),PI/2)(0));
		double wp=A.map.leg.eval_00(A.w.row(-1),0)(0)*A.units.Omega;
		double we=A.map.leg.eval_00(A.w.row(-1),PI/2)(0)*A.units.Omega;
		double Pp=2*PI/wp/3600./24.;
		double Pe=2*PI/we/3600./24.;
		printf("P_rot(p)=%.3f days     P_rot(e)=%.3f days      v_eq=%.2f km/s\n",
			Pp,Pe,we*Re/1e5);
		printf("P_rot(c)=%.3f days\n",2*PI/(A.w(0,0)*A.units.Omega)/3600./24.);
		printf("X=%3.4f (Xc/X=%3.4f) Z=%3.4f\n",A.X0,A.Xc,A.Z0);
		printf("rhoc=%e Tc=%e pc=%e\n",A.rhoc,A.Tc,A.pc);
		if(A.conv) printf("R. conv. core (p)=%3.3f Rsun\n",*(A.map.gl.xif+A.conv)*A.R/R_SUN);
                A.test_virial=A.virial(); A.test_energy=A.energy_test();
		printf("Virial test: %e Energy test: %e\n",A.virial(),A.energy_test());
        printf("\n");
	}

	A.write(config.output_file);

	// if(config.verbose) {
	// 	delete fig;
	// }

	delete op;
	t.stop();
	if(config.verbose) 
		printf("%2.2f seconds\n",t.value());

    plt::show(true);

	return 0;
}

void sig_handler(int sig) {

	char yn;

	if(sig==SIGINT) {
		printf("\nFinish iteration and save model ([y]/n)?");
		if (scanf(" %c",&yn) == 1) {
            if(yn=='y') {
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



