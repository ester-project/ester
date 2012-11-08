#include"star.h"
#include<stdlib.h>
#include<string.h>

void star2d::atmosphere() {

	if(!strcmp(atm_name,"simple")) {
		atm_simple();
    } else {
    	printf("Unknown atmosphere type: %s\n",atm_name);
    	exit(1);
    }
}
void star2d::solve_atm(solver *op) {

	if(!strcmp(atm_name,"simple")) {
		solve_atm_simple(op);
    } else {
    	printf("Unknown atmosphere type: %s\n",atm_name);
    	exit(1);
    }
}
// Using tau_phot=1 and T^4=tau*Teff^4
void star2d::atm_simple() {

	ps=surff*gsup()/opa.k.row(-1)/pc;
	
	matrix q;
	q=surff*p.row(-1)/ps;
	Ts=pow(q,0.25)*Teff()/Tc;
}

void star2d::solve_atm_simple(solver *op) {

	matrix q;
	double qq;
	int n=ndomains-1;
	
	op->bc_top1_add_d(n,"ps","ps",1/ps);
	op->bc_top1_add_d(n,"ps","gsup",-1/gsup());
	op->bc_top1_add_d(n,"ps","opa.k",1/opa.k.row(-1));
	op->bc_top1_add_d(n,"ps","log_pc",ones(1,nth));
	op->set_rhs("ps",zeros(1,nth));
	
	op->bc_top1_add_d(n,"Ts","Ts",1/Ts);
	op->bc_top1_add_d(n,"Ts","p",-0.25/p.row(-1));
	op->bc_top1_add_d(n,"Ts","ps",0.25/ps);
	op->bc_top1_add_d(n,"Ts","Teff",-1/Teff());
	op->bc_top1_add_d(n,"Ts","log_Tc",ones(1,nth));
	op->set_rhs("Ts",zeros(1,nth));
	
}

/* // Using tau_phot=2/3 and T^4=(3/4*tau+1/2)*Teff^4
void star2d::atm_simple() {

	ps=2./3.*surff*gsup()/opa.k.row(-1)/pc;
	
	matrix q;
	q=0.5*(surff*p.row(-1)/ps+1);
	Ts=pow(q,0.25)*Teff()/Tc;
}

void star2d::solve_atm_simple(solver *op) {

	matrix q;
	double qq;
	int n=ndomains-1;
	
	op->bc_top1_add_d(n,"ps","ps",1/ps);
	op->bc_top1_add_d(n,"ps","gsup",-1/gsup());
	op->bc_top1_add_d(n,"ps","opa.k",1/opa.k.row(-1));
	op->bc_top1_add_d(n,"ps","log_pc",ones(1,nth));
	op->set_rhs("ps",zeros(1,nth));
	
	op->bc_top1_add_d(n,"Ts","Ts",1/Ts);
	op->bc_top1_add_d(n,"Ts","p",-0.25*surff/(surff*p.row(-1)+ps));
	op->bc_top1_add_d(n,"Ts","ps",0.25*surff*p.row(-1)/ps/(surff*p.row(-1)+ps));
	op->bc_top1_add_d(n,"Ts","Teff",-1/Teff());
	op->bc_top1_add_d(n,"Ts","log_Tc",ones(1,nth));
	op->set_rhs("Ts",zeros(1,nth));
	
}
*/
