#include"physics.h"
#include"constants.h"

int eos_test(const matrix &X,double Z,const matrix &T,const matrix &p,
		matrix &rho,eos_struct &eos) {
		
	matrix mu,b;
	double RM=8.314385e7;
	
	mu=4./(3+5*X-Z);
	eos.prad=zeros(T.nrows(),T.ncols());
	rho=p/T*mu/RM;
	b=ones(T.nrows(),T.ncols());
	eos.G1=5./3.*b;
	eos.cp=2.5*RM/mu;
	eos.del_ad=2./5.*b;
	eos.G3_1=2./3.*b;
	eos.cv=1.5*RM/mu;
	eos.d=b;
	eos.chi_T=b;
	eos.chi_rho=b;
	
	return 0;
}


