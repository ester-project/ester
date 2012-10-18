#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<signal.h>
#include"star.h"

figure *fig;

int killed=0;
class configuration {
public:
	int minit,maxit;
	double tol,newton_dmax,spec_diff;
	int verbose;
	char plot_device[64];
	double plot_interval;
	char input_file[256];
	char param_file[256];
	char output_file[256];
	char output_mode;
	configuration(int argc,char *argv[]);
	~configuration(){};
	void missing_argument(const char *arg);
	int check_arg(const char *arg,const char *val);
};

void sig_handler(int sig);

int main(int argc,char *argv[]) {

	int nit,last_it;
	double err;
	double t_plot,t_output;
	configuration config(argc,argv);
	tiempo t;
	
	signal(SIGINT,sig_handler);
	
	t.start();
		
	if(config.verbose) {
		fig=new figure(config.plot_device);
		fig->subplot(2,2);
	}
	
	star2d A;
	solver *op;
	
	if(!A.init(config.input_file,config.param_file,argc,argv)) return 1;
	
	nit=0;
	
	matrix tt(config.maxit+1,1),error(config.maxit+1,1);
	t_output=0;
	t_plot=0;
	last_it=nit>=config.maxit;
	op=A.init_solver();
	if(config.verbose>2) op->verbose=1;
	A.config.newton_dmax=config.newton_dmax;
	if(config.verbose>1) A.config.verbose=1;
	while(!last_it) {
		nit++;
		
		//A.check_jacobian(op,"gsup");exit(0);
		err=A.solve(op);

		tt(nit-1)=t.value();
		error(nit-1)=err;
		last_it=(err<config.tol&&nit>=config.minit)||nit>=config.maxit;
		if(killed) last_it=1;
		if(config.verbose) {
			printf("it=%d err=%e (%2.2fs)\n",nit,err,t.value());
			printf("\tOmega=%e (%2.2f%%) eps=%.4f M=%f\n",A.Omega,A.Omega/A.Omegac*100,A.map.eps(A.map.gl.ndomains()-1),A.m*A.rhoc*A.R*A.R*A.R/M_SUN);
			t_output=tt(nit-1);

			if(tt(nit-1)-t_plot>config.plot_interval||last_it) {
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
				
/*
				static figure fig2("/XSERVE");
				fig2.subplot(2,2);
				fig2.colorbar();
				A.drawi(&fig2,A.psi,100,100);
				fig2.colorbar();
				A.drawi(&fig2,A.vt,100,100,11);
				fig2.colorbar();
				A.drawi(&fig2,A.vr,100,100);
				fig2.colorbar();
				A.drawi(&fig2,A.G,100,100,11);
				static figure fig3("/XSERVE");
				fig3.subplot(2,2);
				fig3.colorbar();
				A.spectrum(&fig3,A.psi);
				fig3.colorbar();
				A.spectrum(&fig3,A.vt,11);
				fig3.colorbar();
				A.spectrum(&fig3,A.vr);
				fig3.colorbar();
				A.spectrum(&fig3,A.G,11);
*/
				//A.boundary_layer();				
				
			}

		}
		
	}

	if(config.verbose) {
		printf("\nMass=%.3f Msun  Luminosity=%.3f Lsun\n",
				A.M/M_SUN,A.luminosity()/L_SUN);
		double Re=A.units.r*A.map.leg.eval_00(A.r.row(A.nr()-1),PI/2)(0);
		printf("Radius(p)=%.3f Rsun    Radius(e)=%.3f Rsun\n",
			A.R/R_SUN,Re/R_SUN);
		printf("Teff(p)=%.2f K      Teff(e)=%.2f K\n",
			A.map.leg.eval_00(A.Teff(),0)(0),A.map.leg.eval_00(A.Teff(),PI/2)(0));
		printf("log(geff)(p)=%.2f       log(geff)(e)=%.2f\n",
			A.map.leg.eval_00(log10(A.gsup()),0)(0),A.map.leg.eval_00(log10(A.gsup()),PI/2)(0));
		double wp=A.map.leg.eval_00(A.w.row(A.nr()-1),0)(0)*A.units.Omega;
		double we=A.map.leg.eval_00(A.w.row(A.nr()-1),PI/2)(0)*A.units.Omega;
		double Pp=2*PI/wp/3600./24.;
		double Pe=2*PI/we/3600./24.;
		printf("P_rot(p)=%.3f days     P_rot(e)=%.3f days      v_eq=%.2f km/s\n",
			Pp,Pe,we*Re/1e5);
		printf("P_rot(c)=%.3f days\n",2*PI/(A.w(0,0)*A.units.Omega)/3600./24.);
		printf("X=%3.3f (Xc/X=%3.3f) Z=%3.3f\n",A.X,A.Xc,A.Z);
		printf("rhoc=%e Tc=%e pc=%e\n",A.rhoc,A.Tc,A.pc);
		if(A.conv) printf("R. conv. core (p)=%3.3f Rsun\n",*(A.map.gl.xif+A.conv)*A.R/R_SUN);
		printf("Virial test: %e Energy test: %e\n",A.virial(),A.energy_test());
	}

	A.write(config.output_file,config.output_mode);



	delete fig;

	op->destroy();
	t.stop();
	if(config.verbose) 
		printf("%2.2f seconds\n\n",t.value());	
	
	return 0;
}

void sig_handler(int sig) {

	char yn;

	if(sig==SIGINT) {
		printf("\nFinish iteration and save model (y/n)?");
		scanf(" %c",&yn);
		if(yn=='y') {
			killed=1;
			return;
		}
	}
	exit(sig);
}

configuration::configuration(int argc,char *argv[]) {
	
	int i,i_arg,k;
	char *arg,*val;
	char file[256],line[256];
	cmdline_parser cmd;
	file_parser fp;
	
	verbose=1;
	strcpy(plot_device,"/NULL");
	plot_interval=10;
	strcpy(output_file,"star.out");
	*input_file=0;
	*param_file=0;
	output_mode='b';
	minit=1;
	maxit=200;
	tol=1e-8;
	newton_dmax=0.5;
	
	sprintf(file,"%s/config/star.cfg",ESTER_ROOT);
	if(!fp.open(file)) 
		printf("Can't open configuration file %s\n",file);
	else {
		while(k=fp.get(arg,val)) {			
			if(i=check_arg(arg,val)) {
				printf("Sintax error in configuration file %s, line %d\n",file,k);
				if(i==2) missing_argument(arg);
				if(i==1) {
					printf("Unknown parameter %s\n",arg);
					exit(1);
				}
			}
		}
		fp.close();
	}
	
	cmd.open(argc,argv);
	while(int err_code=cmd.get(arg,val)) {
		if(err_code==-1) exit(1);
		err_code=check_arg(arg,val);
		if(err_code==2) missing_argument(arg);
		if(err_code==0) cmd.ack(arg,val);
	}
	cmd.close();

}

int configuration::check_arg(const char *arg,const char *val) {

	int err=0;

	if(!strcmp(arg,"v0"))
		verbose=0;
	else if(!strcmp(arg,"v1"))
		verbose=1;
	else if(!strcmp(arg,"v2"))
		verbose=2;
	else if(!strcmp(arg,"v3"))
		verbose=3;
	else if(!strcmp(arg,"v4"))
		verbose=4;
	else if(!strcmp(arg,"verbose")) {
		if(val==NULL) return 2;
		verbose=atoi(val);
		verbose=verbose>4?4:verbose;
		verbose=verbose<0?0:verbose;
	} 
	else if(!strcmp(arg,"o")||!strcmp(arg,"output_file")) {
		if(val==NULL) return 2;
		strcpy(output_file,val);
	}
	else if(!strcmp(arg,"i")||!strcmp(arg,"input_file")) {
		if(val==NULL) return 2;
		strcpy(input_file,val);
	}  
	else if(!strcmp(arg,"p")||!strcmp(arg,"param_file")) {
		if(val==NULL) return 2;
		strcpy(param_file,val);
	}  
	else if(!strcmp(arg,"ascii")) 
		output_mode='t';
	else if(!strcmp(arg,"binary"))
		output_mode='b';
	else if(!strcmp(arg,"output_mode")) {
		if(val==NULL) return 2;
		if(val[0]!='b'&&val[0]!='t') 
			printf("Ignoring unknown output_mode %s\n",val);
		else output_mode=val[0];
	}
	else if(!strcmp(arg,"plot_interval")) {
		if(val==NULL) return 2;
		plot_interval=atof(val);
	}
	else if(!strcmp(arg,"plot_device")) {
		if(val==NULL) return 2;
		strcpy(plot_device,val);
	}
	else if(!strcmp(arg,"noplot")) {
		strcpy(plot_device,"/NULL");
	}
	else if(!strcmp(arg,"maxit")) {
		if(val==NULL) return 2;
		maxit=atoi(val);
	}
	else if(!strcmp(arg,"minit")) {
		if(val==NULL) return 2;
		minit=atoi(val);
	}
	else if(!strcmp(arg,"tol")) {
		if(val==NULL) return 2;
		tol=atof(val);
	}
	else if(!strcmp(arg,"newton_dmax")) {
		if(val==NULL) return 2;
		newton_dmax=atof(val);
	}
	else err=1;

	return err;

}

void configuration::missing_argument(const char *arg) {
	
	fprintf(stderr,"Error: Argument to '%s' missing\n",arg);
	exit(1);
}



