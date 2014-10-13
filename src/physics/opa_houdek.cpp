#include<cmath>
#include<string.h>
#include"matrix.h"
#include"constants.h"
#include"physics.h"

extern"C" {
	void maceps_(double *eps);
	void opinit_(double *eps,int *iorder,char *tabnam,int *imode,int tabnamlen);
	void opints_(double *x,double *z,double *tlg,double *rlg,double *opalg,
			double *opr,double *opt,double *opx,double *opz,int *iexp,int *ier);
}

double opa_houdek_i(double X,double Z,double T,double rho,double &dlogkt,double &dlogkr) {
    
    double logT,t6,logR,logk,dlogkx,dlogkz;
    int iexp,ier;
        
    logT=log10(T);
    t6=T*1e-6;
    logR=log10(rho/t6/t6/t6);
    
    opints_(&X,&Z,&logT,&logR,&logk,&dlogkr,&dlogkt,&dlogkx,&dlogkz,&iexp,&ier);
    
	return logk;
}

int opa_houdek(const matrix &X,double Z,const matrix &T,const matrix &rho,
		opa_struct &opa) {

	int i,N,error=0;
	matrix dlnkT,dlnkrho;
    static int init=0;
	char tabnam[80];
	int iorder,imode;
	double eps;

	if(!init) {
    	iorder=4;
		imode=2;
		maceps_(&eps);
		sprintf(tabnam, "%s/ester/tables/houdek/v9/OPINTPATH_AX", ESTER_DATADIR);
		for(i=strlen(tabnam);i<80;i++) tabnam[i]=' ';
		opinit_(&eps,&iorder,tabnam,&imode,80);	
    	init=1;
    }

	opa.k.dim(T.nrows(),T.ncols());
	dlnkT.dim(T.nrows(),T.ncols());
	dlnkrho.dim(T.nrows(),T.ncols());
	N=T.nrows()*T.ncols();
	double dlogkt,dlogkr;
	
	for(i=0;i<N;i++) {
		opa.k(i)=opa_houdek_i(X(i),Z,T(i),rho(i),dlogkt,dlogkr);
		dlnkT(i)=dlogkt;
		dlnkrho(i)=dlogkr;
	}
	opa.k=pow(10,opa.k);
	opa.xi=16*SIG_SB*pow(T,3)/(3*opa.k*rho);
	opa.dlnxi_lnrho=-1-dlnkrho;
    opa.dlnxi_lnT=3-dlnkT;
	
	return error;
		
}
		
