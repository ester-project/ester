#include"star.h"
#include<stdlib.h>
#include<string.h>

void star1d::atmosphere() {

	if(!strcmp(atm_name,"simple")) {
		atm_simple();
	} else if(!strcmp(atm_name,"test")) {
    	atm_test();
    } else {
    	printf("Unknown atmosphere type: %s\n",atm_name);
    	exit(1);
    }
}
void star1d::solve_atm(solver *op) {

	if(!strcmp(atm_name,"simple")) {
		solve_atm_simple(op);
	} else if(!strcmp(atm_name,"test")) {
    	solve_atm_test(op);
    } else {
    	printf("Unknown atmosphere type: %s\n",atm_name);
    	exit(1);
    }
}

void star1d::atm_simple() {

	ps=surff*(D.row(nr()-1),phi)(0)/opa.k(nr()-1)/rhoc/R;
	Ts=surff*p(nr()-1)/ps/Tc/Tc/Tc/R*Frad(nr()-1)*Lambda/SIG_SB;
	Ts=pow(Ts,1./4.);
	
}

void star1d::solve_atm_simple(solver *op) {

	double qq;
	matrix q;

	qq=1.;
	op->bc_top1_add_d(ndomains()-1,"ps","ps",qq*ones(1,1));
	qq=-1./(D,phi)(nr()-1);
	op->bc_top1_add_l(ndomains()-1,"ps","phi",qq*ones(1,1),D.block(ndomains()-1).row(gl.npts[ndomains()-1]-1));
	qq=1.;
	op->bc_top1_add_d(ndomains()-1,"ps","rhoc",qq*ones(1,1));
	qq=1.;
	op->bc_top1_add_d(ndomains()-1,"ps","R",qq*ones(1,1));
	
	q=-1./eos.chi_rho*(1+opa.dlnxi_lnrho);
	q=q.row(nr()-1);
	op->bc_top1_add_d(ndomains()-1,"ps","p",q);
	q=-1./eos.chi_rho*(1+opa.dlnxi_lnrho);
	q=q.row(nr()-1);
	op->bc_top1_add_d(ndomains()-1,"ps","pc",q);
	q=1.*(3+eos.d+eos.d*opa.dlnxi_lnrho-opa.dlnxi_lnT);
	q=q.row(nr()-1);
	op->bc_top1_add_d(ndomains()-1,"ps","T",q);
	q=1.*(3+eos.d+eos.d*opa.dlnxi_lnrho-opa.dlnxi_lnT);
	q=q.row(nr()-1);
	op->bc_top1_add_d(ndomains()-1,"ps","Tc",q);
	
	op->bc_top1_add_d(ndomains()-1,"ps","dx",ones(1,1)/(gl.xif[ndomains()]-gl.xif[ndomains()-1]));
	
	op->set_rhs("ps",zeros(1,1));
	
	qq=4.;
	op->bc_top1_add_d(ndomains()-1,"Ts","Ts",qq*ones(1,1));
	
	qq=1.;
	op->bc_top1_add_d(ndomains()-1,"Ts","R",qq*ones(1,1));
	qq=1.;
	op->bc_top1_add_d(ndomains()-1,"Ts","ps",qq*ones(1,1));
	qq=-1.;
	op->bc_top1_add_d(ndomains()-1,"Ts","Lambda",qq*ones(1,1));
	
	q=3*ones(1,1);
	op->bc_top1_add_d(ndomains()-1,"Ts","Tc",q);
	q=-ones(1,1);
	op->bc_top1_add_d(ndomains()-1,"Ts","p",q);
	
	q=-ones(1,1)/Frad(nr()-1);
	op->bc_top1_add_d(ndomains()-1,"Ts","Frad",q);
	
	op->set_rhs("Ts",zeros(1,1));
	
}

// Michel's test model

void star1d::atm_test() {

	double Lt=1/0.4;
	
	ps=1.2e-5;
	//ps=1e-4;
	Ts=Lambda*Lt*Frad(nr()-1)/opa.xi(nr()-1);
}

void star1d::solve_atm_test(solver *op) {

	double qq;
	matrix q;

	qq=1.;
	op->bc_top1_add_d(ndomains()-1,"ps","ps",qq*ones(1,1));
	
	op->set_rhs("ps",zeros(1,1));
	
	qq=1.;
	op->bc_top1_add_d(ndomains()-1,"Ts","Ts",qq*ones(1,1));
	q=-ones(1,1)/Frad(nr()-1);
	op->bc_top1_add_d(ndomains()-1,"Ts","Frad",q);
	qq=-1.;
	op->bc_top1_add_d(ndomains()-1,"Ts","Lambda",qq*ones(1,1));
	
	q=opa.dlnxi_lnrho/eos.chi_rho;
	q=q(nr()-1)*ones(1,1);
	op->bc_top1_add_d(ndomains()-1,"Ts","p",q);
	op->bc_top1_add_d(ndomains()-1,"Ts","pc",q);
	q=opa.dlnxi_lnT-opa.dlnxi_lnrho*eos.d;
	q=q(nr()-1)*ones(1,1);
	op->bc_top1_add_d(ndomains()-1,"Ts","T",q);
	op->bc_top1_add_d(ndomains()-1,"Ts","Tc",q);
	
	op->set_rhs("Ts",zeros(1,1));
	
}
