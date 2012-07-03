#include<math.h>
#include"physics.h"

int nuc_test(const matrix &X,double Z,const matrix &T,const matrix &rho,nuc_struct &nuc) {

	matrix T9;
	
	T9=T/1e9;

	nuc.pp=8.24e4*X*X*pow(T9,-2./3.)*exp(-3.6*pow(T9,-1./3.));
	nuc.cno=0*nuc.pp;

	nuc.pp*=rho;
	nuc.cno*=rho;
	nuc.eps=nuc.pp+nuc.cno;
	nuc.dlneps_lnrho=ones(T.nrows(),T.ncols());
	nuc.dlneps_lnT=(-2./3.+1.2*pow(T9,-1./3.));
	
	return 0;
}
