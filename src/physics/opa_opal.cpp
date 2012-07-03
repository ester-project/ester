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
    int mzin,izi,iz1e4;
    
    iz1e4=round(Z*1e4);
    switch(iz1e4) {
    	case 0:
    		mzin=0;
    		break;
    	case 1:
    		mzin=2;
    		break;
    	case 3:
    		mzin=3;
    		break;
    	case 10:
    		mzin=4;
    		break;
    	case 20:
    		mzin=5;
    		break;
    	case 40:
    		mzin=6;
    		break;
    	case 100:
    		mzin=7;
    		break;
    	case 200:
    		mzin=8;
    		break;
    	case 300:
    		mzin=9;
    		break;
    	case 400:
    		mzin=10;
    		break;
    	case 600:
    		mzin=11;
    		break;
    	case 800:
    		mzin=12;
    		break;
    	case 1000:
    		mzin=13;
    		break;
    	default:
    		mzin=0;
    }

    izi=0;
    
    t6=T*1e-6;
    r=rho/t6/t6/t6;
    if(mzin) opac_(&izi,&mzin,&X,&t6,&r);
    	else opacgn93_(&Z,&X,&t6,&r);
    
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
		
