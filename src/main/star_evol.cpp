#include "ester-config.h"
#include "utils.h"
#include "ester.h"
#include <string.h>
#include <stdlib.h>

#include "read_config.h"

int main(int argc,char *argv[]) {

	configuration config(argc,argv);
	cmdline_parser cmd;
	
	double dXc=0.05,Xcmin=0.05;
	
	char *arg,*val;
	cmd.open(argc,argv);
	while(int err_code=cmd.get(arg,val)) {
		if(err_code==-1) exit(1);
		err_code=0;
		if(!strcmp(arg,"dXc")) {
			if(val==NULL) err_code=2;
			else dXc=atof(val);
		} else if(!strcmp(arg,"Xcmin")) {
			if(val==NULL) err_code=2;
			else Xcmin=atof(val);
		} else err_code=1;
		if(err_code==1) {
			ester_err("Unknown parameter %s",arg);
			exit(1);
		}
		if(err_code==2) {
			ester_err("Argument to %s missing",arg);
			exit(1);
		}
		cmd.ack(arg,val);
	}
	cmd.close();
	
	if(*config.input_file==0) {
		ester_err("Must specify an input file");
		exit(1);
	}
	if(*config.output_file==0) {
		strcpy(config.input_file,config.output_file);
	}

	star_evol A;

	if(A.read(config.input_file)) {
		star1d A1d;
		if(A1d.read(config.input_file)) {
			ester_err("Error reading input file %s", config.input_file);
			exit(1);
		}
		A=A1d;
	}
	
	figure *fig;
	solver *op;
	
	if(config.verbose) {
		fig=new figure(config.plot_device);
		fig->subplot(2,2);
	}
	
	op=A.init_solver();
	if(config.verbose>2) op->verbose=1;
	A.config.newton_dmax=config.newton_dmax;
	if(config.verbose>1) A.config.verbose=1;
	
	double Xc=A.Xc;
	int n=0;
	char outfile[256];
	
	while((Xc >= Xcmin && Xc <= 1)) {
		printf("Xc=%f\n",Xc);
		int last_it=0,nit=0;
		double err;
		A.Xc=Xc;
		while(!last_it) {
			nit++;
            A.converged = false;
			err=A.solve(op);
			last_it=(err<config.tol&&nit>=config.minit)||nit>=config.maxit;
			if(config.verbose) {
				printf("\tit=%d err=%e\n",nit,err);
				printf("\t\tOmega=%e (%2.2f%%) eps=%.4f M=%f\n",A.Omega,A.Omega/A.Omegac*100,1-1./A.map.leg.eval_00(A.r.row(-1),PI/2)(0),A.m*A.rhoc*A.R*A.R*A.R/M_SUN);

			}
		}
        A.converged = true;
		
		static matrix logL(log10(A.luminosity()/L_SUN)*ones(1,1));
		static matrix logR(log10(A.R/R_SUN)*ones(1,1));
		static matrix logTeff(log10(A.Teff()(-1))*ones(1,1));
		static matrix Wbk(A.Omega_bk*ones(1,1));
		static matrix XX(A.X0*A.Xc*ones(1,1));
		static matrix Lzc(A.Lzcore()*ones(1,1));
		if(n) {
			logL=logL.concatenate(log10(A.luminosity()/L_SUN)*ones(1,1));
			logR=logR.concatenate(log10(A.R/R_SUN)*ones(1,1));
			logTeff=logTeff.concatenate(log10(A.Teff()(-1))*ones(1,1));
			Wbk=Wbk.concatenate(A.Omega_bk*ones(1,1));
			XX=XX.concatenate(A.X0*A.Xc*ones(1,1));
			Lzc=Lzc.concatenate(A.Lzcore()*ones(1,1));
		}
		
		if(config.verbose) {
			fig->subplot(2,2);
			fig->colorbar();
			A.draw(fig,A.w);
			fig->label("Differential rotation","","");
			fig->plot(-logTeff,logL);
			fig->label("-log(Teff(pol))","log(L/Lsun)","");
			fig->axis(1.05*max(XX),0,min(Wbk)*0.95,max(Wbk)*1.05);
			fig->plot(XX,Wbk);
			fig->label("X(core)","Omega_bk","");
			fig->axis(1.05*max(XX),0,min(Lzc/A.Lz())*0.95,max(Lzc/A.Lz())*1.05);
			fig->plot(XX,Lzc/A.Lz());
			fig->label("X(core)","Lz_core/Lz","");
		}
		sprintf(outfile,"%s_%04d",config.output_file,n);
        {
            char *filename = NULL;
            asprintf(&filename, "%s-%04d.hdf5",
                    config.output_file,
                    n);
            A.hdf5_write(filename);
            free(filename);
        }	A.write(outfile,config.output_mode);
		Xc = A.comp["H"](0)/A.comp["H"](-1);
        A.fill(); // this updated the chemical composition
		n++;
	}
}

