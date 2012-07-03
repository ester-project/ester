#include"star.h"
#include<stdlib.h>
#include<string.h>

void star2d::atmosphere() {

	if(!strcmp(atm_name,"simple")) {
		atm_simple();
	} else if(!strcmp(atm_name,"test")) {
    	atm_test();
    } else {
    	printf("Unknown atmosphere type: %s\n",atm_name);
    	exit(1);
    }
}
void star2d::solve_atm(solver *op) {

	if(!strcmp(atm_name,"simple")) {
		solve_atm_simple(op);
	} else if(!strcmp(atm_name,"test")) {
    	solve_atm_test(op);
    } else {
    	printf("Unknown atmosphere type: %s\n",atm_name);
    	exit(1);
    }
}

void star2d::atm_simple() {

	ps=surff*gsup_/opa.k.row(nr()-1)/rhoc/R;
	Ts=surff*p.row(nr()-1)/ps/SIG_SB/Tc/Tc/Tc/R*Frad.row(nr()-1)*Lambda;
	Ts=pow(Ts,1./4.);
}

void star2d::solve_atm_simple(solver *op) {

	matrix q;
	double qq;
	
	q=ones(1,nth());
	op->bc_top1_add_d(ndomains()-1,"ps","ps",q);
	
	op->bc_top1_add_d(ndomains()-1,"ps","gsup",-1/gsup_);
	
	qq=1.;
	op->bc_top1_add_d(ndomains()-1,"ps","rhoc",qq*ones(1,nth()));
	qq=1.;
	op->bc_top1_add_d(ndomains()-1,"ps","R",qq*ones(1,nth()));
	
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
	
	op->set_rhs("ps",zeros(1,nth()));
	
	q=4.*ones(1,nth());
	op->bc_top1_add_d(ndomains()-1,"Ts","Ts",q);
	
	qq=1.;
	op->bc_top1_add_d(ndomains()-1,"Ts","R",qq*ones(1,nth()));
	qq=-1.;
	op->bc_top1_add_d(ndomains()-1,"Ts","Lambda",qq*ones(1,nth()));
	q=ones(1,nth());
	op->bc_top1_add_d(ndomains()-1,"Ts","ps",q);
	
	q=3*ones(1,nth());
	op->bc_top1_add_d(ndomains()-1,"Ts","Tc",q);
	q=-ones(1,nth());
	op->bc_top1_add_d(ndomains()-1,"Ts","p",q);
	
	q=-1/Frad.row(nr()-1);
	op->bc_top1_add_d(ndomains()-1,"Ts","Frad",q);
	
	
	op->set_rhs("Ts",zeros(1,nth()));
	
}

// Michel's test model

void star2d::atm_test() {

	double Lt=1/0.4;

	ps=1.2e-5*ones(1,nth());
	Ts=Lambda*Lt*Frad.row(nr()-1)/opa.xi.row(nr()-1)*map.leg.eval_00(r.row(nr()-1),PI/2)(0);
	
}

void star2d::solve_atm_test(solver *op) {

	matrix q,TT;
	
	q=ones(1,nth());
	op->bc_top1_add_d(ndomains()-1,"ps","ps",q);	
	op->set_rhs("ps",zeros(1,nth()));
	
	
	q=ones(1,nth());
	op->bc_top1_add_d(ndomains()-1,"Ts","Ts",q);
	
	q=-1/Frad.row(nr()-1);
	op->bc_top1_add_d(ndomains()-1,"Ts","Frad",q);
	
	q=-ones(1,nth());
	op->bc_top1_add_d(ndomains()-1,"Ts","Lambda",q);
	
	map.leg.eval_00(th,PI/2,TT);
	op->bc_top1_add_ri(ndomains()-1,"Ts","Ri",-ones(1,nth()),TT,1./r.row(nr()-1));
	
	q=opa.dlnxi_lnrho/eos.chi_rho;
	q=q.row(nr()-1);
	op->bc_top1_add_d(ndomains()-1,"Ts","p",q);
	op->bc_top1_add_d(ndomains()-1,"Ts","pc",q);
	q=opa.dlnxi_lnT-opa.dlnxi_lnrho*eos.d;
	q=q.row(nr()-1);
	op->bc_top1_add_d(ndomains()-1,"Ts","T",q);
	op->bc_top1_add_d(ndomains()-1,"Ts","Tc",q);
	
	
	op->set_rhs("Ts",zeros(1,nth()));
	
}

