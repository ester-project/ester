#include "ester-config.h"
#include "physics.h"
#include "constants.h"

int eos_idealrad(const matrix &X,double Z,const matrix &T,const matrix &p,
		matrix &rho,eos_struct &eos) {

	matrix mu,b;
	
	mu=4/(3+5*X-Z);
	eos.prad=A_RAD/3*pow(T,4);
	rho=(p-eos.prad)/T/K_BOL*mu*HYDROGEN_MASS;
	b=(p-eos.prad)/p;
	eos.G1=b+2./3.*(4-3*b)*(4-3*b)/(8-7*b);
	eos.cp=1.5*K_BOL/mu/HYDROGEN_MASS*((8-7*b)/b+2./3.*(4-3*b)*(4-3*b)/b/b);;
	eos.del_ad=1./(4-1.5*b*b/(4-3*b));
	eos.G3_1=2./3.*(4-3*b)/(8-7*b);
	eos.cv=b*eos.cp/eos.G1;
	eos.d=(4-3*b)/b;
	eos.chi_T=4-3*b;
	eos.chi_rho=b;
	eos.s=K_BOL/mu/HYDROGEN_MASS*log(pow(T,1.5)/rho)+4.*A_RAD/3.*T*T*T/rho;
	
	return 0;
}


