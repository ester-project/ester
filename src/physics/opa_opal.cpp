#include<cmath>
#include"matrix.h"
#include"constants.h"
#include"physics.h"

extern"C" {
	void opac_(int *izi,int *mzin,double *X,double *t6,double *r);
	void opacgn93_(double *Z,double *X,double *t6,double *r);
	extern struct{
		double opact,dopact,dopacr,dopactd;
	} e_;
}

double opa_opal_i(double X,double Z,double T,double rho,double &dlogkt,double &dlogkr) {
    
    double t6,r;
    
    t6=T*1e-6;
    r=rho/t6/t6/t6;
    
    opacgn93_(&Z,&X,&t6,&r);
    dlogkt=e_.dopact;
    dlogkr=e_.dopacr;
    
	return e_.opact;
}

int opa_opal(const matrix &X,double Z,const matrix &T,const matrix &rho,
		opa_struct &opa) {

	int i,N,error=0;
	matrix dlnkT,dlnkrho;

	opa.k.dim(T.nrows(),T.ncols());
	dlnkT.dim(T.nrows(),T.ncols());
	dlnkrho.dim(T.nrows(),T.ncols());
	N=T.nrows()*T.ncols();
	double dlogkt,dlogkr;
	
	for(i=0;i<N;i++) {
		opa.k(i)=opa_opal_i(X(i),Z,T(i),rho(i),dlogkt,dlogkr);
		if(opa.k(i)==-99) error=1;
		dlnkT(i)=dlogkt;
		dlnkrho(i)=dlogkr;
	}
	opa.k=pow(10,opa.k);
	dlnkT-=3*dlnkrho;
	opa.xi=16*SIG_SB*pow(T,3)/(3*opa.k*rho);
	opa.dlnxi_lnrho=-1-dlnkrho;
    opa.dlnxi_lnT=3-dlnkT;
	if(error) printf("Values outside OPAL opacity table\n");
	
	return error;
		
}
		
