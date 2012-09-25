#include<math.h>
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

double opa_opal_i(double X,double Z,double T,double rho,double *dlogkt,double *dlogkr) {
    
    double t6,r;
    
    t6=T*1e-6;
    r=rho/t6/t6/t6;
    
    opacgn93_(&Z,&X,&t6,&r);
    *dlogkt=e_.dopact;
    *dlogkr=e_.dopacr;
    
	return e_.opact;
}

int opa_opal(const matrix &X,double Z,const matrix &T,const matrix &rho,
		opa_struct &opa) {

	int i,N,error=0;
	double *pX,*pT,*prho,*pk,*pdkt,*pdkr;
	matrix dlnkT,dlnkrho;

	opa.k.dim(T.nrows(),T.ncols());
	dlnkT.dim(T.nrows(),T.ncols());
	dlnkrho.dim(T.nrows(),T.ncols());
	N=T.nrows()*T.ncols();
	pX=X.data();
	pT=T.data();
	prho=rho.data();
	pk=opa.k.data();
	pdkt=dlnkT.data();
	pdkr=dlnkrho.data();
	
	for(i=0;i<N;i++) {
		pk[i]=opa_opal_i(pX[i],Z,pT[i],prho[i],pdkt+i,pdkr+i);
		if(pk[i]==-99) error=1;
	}
	opa.k=pow(10,opa.k);
	dlnkT-=3*dlnkrho;
	opa.xi=16*SIG_SB*pow(T,3)/(3*opa.k*rho);
	opa.dlnxi_lnrho=-1-dlnkrho;
    opa.dlnxi_lnT=3-dlnkT;
	if(error) printf("Values outside OPAL opacity table\n");
	
	return error;
		
}
		
