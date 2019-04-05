#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "utils.h"
#include "ester.h"
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include "matplotlib.h"

#include "read_config.h"

void sig_handler(int sig) {
    exit(sig);
}

int main(int argc,char *argv[]) {

	signal(SIGINT, sig_handler);

	configuration config(argc,argv);
	cmdline_parser cmd;

	double step=1, age=0;
	double maxStep = 10, minStep = 1e-4;
	double min_max_err = 1e-2, max_max_err = 5e-2;
	
	char *arg,*val;
	cmd.open(argc,argv);
	while(int err_code=cmd.get(arg,val)) {
		if(err_code==-1) exit(1);
		err_code=0;
		if(!strcmp(arg,"step")) {
			if(val==NULL) err_code=2;
			else step=atof(val);
		} else if(!strcmp(arg,"age")) {
			if(val==NULL) err_code=2;
			else age=atof(val);
		} else if(!strcmp(arg,"max_step")) {
			if(val==NULL) err_code=2;
			else maxStep=atof(val);
		} else if(!strcmp(arg,"min_step")) {
			if(val==NULL) err_code=2;
			else minStep=atof(val);
		} else if(!strcmp(arg,"min_error")) { // Newton error
			if(val==NULL) err_code=2;
			else min_max_err=atof(val);
		} else if(!strcmp(arg,"max_error")) { // Newton error
			if(val==NULL) err_code=2;
			else max_max_err=atof(val);
		} else err_code=1;
		/*if(err_code==1) {
			fprintf(stderr,"Unknown parameter %s\n",arg);
			exit(1);
		}*/
		if(err_code==2) {
			fprintf(stderr,"Argument to %s missing\n",arg);
			exit(1);
		}
		if(err_code==0) cmd.ack(arg,val);
	}
	cmd.close();

	if(*config.input_file==0) {
		fprintf(stderr,"Must specify an input file\n");
		exit(1);
	}
	if(*config.output_file==0) {
		strcpy(config.input_file,config.output_file);
	}

	star_evol A;
	if(A.init(config.input_file,config.param_file,argc,argv)) {
		ester_err("Could not initialize star");
		return 1;
	}

	/*if(A.read(config.input_file)) {
		star1d A1d;
		if(A1d.read(config.input_file)) {
			fprintf(stderr,"Error reading input file %s\n",config.input_file);
			exit(1);
		}
		A=A1d;
	}*/

	//figure *fig;
	solver *op;
	SDIRK_solver *rk;
	/*
	if(config.verbose) {
		fig=new figure(config.plot_device);
		fig->subplot(2,2);
	}
	*/
	op = A.init_solver();
	rk = A.init_time_solver();
	if(config.verbose>2) op->verbose=1;
	A.config.newton_dmax=config.newton_dmax;
	if(config.verbose>1) A.config.verbose=1;
	
	int n=0;
	char outfile[256];
	int state, stage = 0;
	rk->set_step(step);
	star_evol A0 = A;
	double max_err = 0;

	while ((state = rk->solve(A.age, age)) != RK_END) {
		A.init_step(rk);
		int last_it=0,nit=0;
		double err;

		printf("Age = %f Myr (step = %e)\n", A.age, step);

		printf("RK stage %d / %d:\n", ++stage, rk->number_of_stages());
		while(!last_it) {
			nit++;
			err=A.solve(op);
			last_it=(err<config.tol&&nit>=config.minit)||nit>=config.maxit;
			if(config.verbose) {
				printf("\tit=%d err=%e\n",nit,err);
				//printf("\t\tOmega=%e (%2.2f%%) eps=%.4f M=%f\n",A.Omega,A.Omega/A.Omegac*100,1-1./A.map.leg.eval_00(A.r.row(-1),PI/2)(0),A.m*A.rhoc*A.R*A.R*A.R/M_SUN);
			}
			max_err = err > max_err ? err : max_err;
			if (err > max_max_err && step > minStep) break;
		}
		if (nit >= config.maxit && step == minStep) {
			printf("No convergence\n");
			exit(1);
		}
		if ((max_err > max_max_err || nit >= config.maxit) && step > minStep) {
			step /= 2;
			if (step < minStep) step = minStep;
			A = A0;
			A.reset_time_solver(rk);
			printf("Decreasing step (%e)\nRestarting iteration...\n", step);
			rk->set_step(step);
			max_err = 0;
			stage = 0;
			continue;
		}

		if(state == RK_STEP) { // Final step of the 3 RK steps (save the results)
			stage = 0;
			if (max_err < min_max_err) {
				step *= 2;
				if (step > maxStep) {
					step = maxStep;
				}
				else {
					printf("Increasing step (%e)\n", step);
					rk->set_step(step);
				}
			}
			max_err = 0;
			static matrix logL(log10(A.luminosity()/L_SUN)*ones(1,1));
			static matrix logR(log10(A.R/R_SUN)*ones(1,1));
			static matrix logTeff(log10(A.Teff()(-1))*ones(1,1));
			static matrix Wbk(A.Omega_bk*ones(1,1));
			static matrix RR(A.R*ones(1,1));
			static matrix tt(A.age*ones(1,1));
			if(n) {
				logL=logL.concatenate(log10(A.luminosity()/L_SUN)*ones(1,1));
				logR=logR.concatenate(log10(A.R/R_SUN)*ones(1,1));
				logTeff=logTeff.concatenate(log10(A.Teff()(-1))*ones(1,1));
				Wbk=Wbk.concatenate(A.Omega_bk*ones(1,1));
				RR=RR.concatenate(A.R*ones(1,1));
				tt=tt.concatenate(A.age*ones(1,1));
			}
		
			if(config.verbose) {
				printf("Mcore/M: %e\n", A.Mcore()/A.M);
				printf("Mean Rcore/R: %e (conv=%d)\n", (A.Rcore(), A.map.It)(0)/(A.r.row(-1), A.map.It)(0)/A.R, A.conv);
				printf("Xcore: %e\n", A.comp["H"](0,0));

				matrix theta = vector(0, M_PI/2., 64);
				matrix r = A.map.leg.eval_00(A.r, theta);
				matrix w = A.map.leg.eval_00(A.w, theta);
				matrix vr = A.map.leg.eval_00(A.vr, theta);
				matrix X = A.map.leg.eval_00(A.comp.X(), theta);

				matrix cost = cos(theta);
				matrix sint = sin(theta);

				matrix x = r*sint;
				matrix y = r*cost;

				plt::clf();

				plt::subplot(231);
				plt::pcolormesh(x, y, w);
				plt::colorbar();
				plt::axis("scaled");
				plt::title("$\\Omega$");

				plt::subplot(232);
				plt::pcolormesh(x, y, vr);
				plt::clim(-max(abs(vr)), max(abs(vr)));
				plt::colorbar();
				plt::axis("scaled");
				plt::title("$V_r$");

				plt::subplot(233);
				plt::pcolormesh(x, y, X);
				plt::clim(0, A.X0);
				plt::colorbar();
				plt::axis("scaled");
				plt::title("X");


				plt::subplot(234);
				plt::plot(-logTeff, logL);
				plt::xlabel("$-\\log T_\\mathrm{eff, p}$");
				plt::ylabel("$\\log L$");

				plt::subplot(235);
				plt::plot(tt, Wbk);
				plt::xlabel("Age (Myr)");
				plt::title("$\\Omega / \\Omega_c$");

				plt::subplot(236);
				plt::plot(tt, RR/R_SUN);
				plt::xlabel("Age (Myr)");
				plt::title("$R_p / R_\\odot$");

				plt::draw();
				plt::pause();

/*
				fig->subplot(2,4);
				fig->colorbar();
				A.drawi(fig,A.w,100,64);
				fig->label("Differential rotation","","");
				A.drawi(fig,A.vr,100,64);
				fig->label("Meridional circulation","","");
				fig->caxis(0, A.X0);
				fig->colorbar();
				A.drawi(fig,A.comp.X(),100,64);
				fig->label("X","","");
				A.spectrum(fig, A.vt, 11);
				fig->label("Meridional circulation","","");
				fig->plot(-logTeff,logL);
				fig->label("-log(Teff(pol))","log(L/Lsun)","");
				fig->plot(tt,Wbk);
				fig->label("age","Omega_bk","");
				fig->plot(tt,RR/R_SUN);
				fig->label("age","Radius","");
				A.spectrum(fig, A.p);
				fig->label("Pressure","","");
*/
			}
			char *filename = strdup(config.output_file);
			char *ext = strstr(filename, ".");
			if (ext != NULL) {
				*ext = '\0';
				ext++;
				sprintf(outfile,"%s_%04d.%s",filename,n,ext);
			}
			else {
				sprintf(outfile,"%s_%04d",filename,n);
			}

			free(filename);
			printf("Writing %s\n", outfile);
			A.write(outfile,config.output_mode);
			n++;
		}
		A.finish_step(rk, state);
		if(state == RK_STEP) A0 = A;
	}
}




