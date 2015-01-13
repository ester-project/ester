#include "ester-config.h"
#include"physics.h"
#include"constants.h"

int eos_ideal(const matrix &X,double Z,const matrix &T,const matrix &p,
		matrix &rho,eos_struct &eos) {
		
	matrix mu,b;
	
	mu=4/(3+5*X-Z);
	eos.prad=zeros(T.nrows(),T.ncols());
	rho=p/T/K_BOL*mu*HYDROGEN_MASS;
	b=ones(T.nrows(),T.ncols());
	eos.G1=5./3.*b;
	eos.cp=2.5*K_BOL/mu/HYDROGEN_MASS;
	eos.del_ad=2./5.*b;
	eos.G3_1=2./3.*b;
	eos.cv=1.5*K_BOL/mu/HYDROGEN_MASS;
	eos.d=b;
	eos.chi_T=b;
	eos.chi_rho=b;
	eos.s=K_BOL/mu/HYDROGEN_MASS*log(pow(T,1.5)/rho);
	
	return 0;
}


