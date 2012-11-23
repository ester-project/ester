#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include"star.h"

figure *fig;

class configuration {
public:
	int minit,maxit;
	double tol,newton_dmax;
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

int main(int argc,char *argv[]) {

	int nit,last_it;
	double err;
	tiempo t;
	double t_plot,t_output;
	configuration config(argc,argv);
	
	t.start();
	if(config.verbose) {
		fig=new figure(config.plot_device);
		fig->subplot(2,1);
	}

	
	star1d A;
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
	
	int conv_set=A.conv;
	double Xc_set=A.Xc;
	// If no input file, ignore core convection until the model starts to converge
	if(*config.input_file==0) {
		if(A.conv) A.Xc=1;
		A.conv=0;
	}
	err=1;
	while(!last_it) {
		if(err<0.1&&!*config.input_file) {
			A.conv=conv_set;
			A.Xc=Xc_set;
		}
		
		nit++;
		//A.check_jacobian(op,"log_T");exit(0);
		err=A.solve(op);
		
		tt(nit-1)=t.value();
		error(nit-1)=err;
		last_it=(err<config.tol&&nit>=config.minit)||nit>=config.maxit;
		if(config.verbose) {
			printf("it=%d err=%e\n",nit,err);
			t_output=tt(nit-1);
			
			if(tt(nit-1)-t_plot>config.plot_interval||last_it) {
				fig->semilogy(error.block(0,nit-1,0,0));
				fig->label("Iteration number","Relative error","");
				A.spectrum(fig,A.rho);
				fig->label("Density (normalized spectrum)","","");
				t_plot=tt(nit-1);
			}

		}
	}
	if(config.verbose) {
		printf("Mass=%3.3f Msun  Radius=%3.3f Rsun  Luminosity=%3.3f Lsun  Teff=%1.1f K\n",
				A.M/M_SUN,A.R/R_SUN,A.luminosity()/L_SUN,A.Teff()(0));
		printf("X=%3.3f (Xc/X=%3.3f) Z=%3.3f\n",A.X,A.Xc,A.Z);
		printf("rhoc=%e Tc=%e pc=%e\n",A.rhoc,A.Tc,A.pc);
		if(A.conv) printf("r_cz=%3.3f Rsun\n",*(A.map.gl.xif+A.conv)*A.R/R_SUN);
	}
	op->destroy();
	A.write(config.output_file,config.output_mode);
	
	delete fig;

	t.stop();
	if(config.verbose) 
		printf("%2.2f seconds\n",t.value());	
	return 0;
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
	else if(!strcmp(arg,"verbose")) {
		if(val==NULL) return 2;
		verbose=atoi(val);
		verbose=verbose>3?3:verbose;
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
	else if(!strcmp(arg,"noplot")) {
		strcpy(plot_device,"/NULL");
	}
	else if(!strcmp(arg,"plot_device")) {
		if(val==NULL) return 2;
		strcpy(plot_device,val);
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



