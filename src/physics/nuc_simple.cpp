#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include<cmath>
#include"physics.h"

int nuc_simple(const composition_map &comp,const matrix &T,const matrix &rho,nuc_struct &nuc) {

	matrix T9;
	double C,N,XCNO;
	double Z=comp.Z()(-1);
	matrix X(comp.X());
	
	T9=T/1e9;
	C=0.173285;
	N=0.053152;
	// O=0.482273;
	XCNO=Z*(C+N);

// nuc.pp from the CESAM code (cited in Espinosa Lara Rieutord 2007)
// Actually used in routine pp1.f of CESAM code citing GONG

	nuc.pp=8.2398e4*X*X*pow(T9,-2./3.)*exp(-3.6*pow(T9,-1./3.));
// nuc.cno from Kippenhahn & Weigert (1991)
	nuc.cno=8.67e25*XCNO*X*pow(T9,-2./3.)*exp(-15.228*pow(T9,-1./3.))
		*(1+0.027*pow(T9,1./3.)-0.778*pow(T9,2./3.)-0.149*T9);

	nuc.pp*=rho;
	nuc.cno*=rho;
	nuc.eps=nuc.pp+nuc.cno;
	nuc.dlneps_lnrho=ones(T.nrows(),T.ncols());
	nuc.dlneps_lnT=nuc.pp/nuc.eps*(-2./3.+1.2*pow(T9,-1./3.))+
		nuc.cno/nuc.eps*(-2./3.+15.228/3.*pow(T9,-1./3.)+
		(0.027/3.*pow(T9,1./3.)-2.*0.778/3.*pow(T9,2./3.)-0.149*T9)/(1+0.027*pow(T9,1./3.)-0.778*pow(T9,2./3.)-0.149*T9));
	
	return 0;
}
