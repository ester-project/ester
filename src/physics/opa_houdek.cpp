#include<math.h>
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

double opa_houdek_i(double X,double Z,double T,double rho,double *dlogkt,double *dlogkr) {
    
    double logT,t6,logR,logk,dlogkx,dlogkz;
    int iexp,ier;
        
    logT=log10(T);
    t6=T*1e-6;
    logR=log10(rho/t6/t6/t6);
    
    opints_(&X,&Z,&logT,&logR,&logk,dlogkr,dlogkt,&dlogkx,&dlogkz,&iexp,&ier);
    
	return logk;
}

int opa_houdek(const matrix &X,double Z,const matrix &T,const matrix &rho,
		opa_struct &opa) {

	int i,N,error=0;
	double *pX,*pT,*prho,*pk,*pdkt,*pdkr;
	matrix dlnkT,dlnkrho;
    static int init=0;
	char tabnam[81];
	int iorder,imode;
	double eps;

	if(!init) {
    	iorder=4;
		imode=2;
		maceps_(&eps);
		sprintf(tabnam,"%s/tables/houdek/v9/OPINTPATH_AX",ESTER_ROOT);
		for(i=strlen(tabnam);i<80;i++) tabnam[i]=' ';
		tabnam[80]='\0';
		opinit_(&eps,&iorder,tabnam,&imode,80);	
    	init=1;
    }

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
		pk[i]=opa_houdek_i(pX[i],Z,pT[i],prho[i],pdkt+i,pdkr+i);
	}
	opa.k=pow(10,opa.k);
	opa.xi=16*SIG_SB*pow(T,3)/(3*opa.k*rho);
	opa.dlnxi_lnrho=-1-dlnkrho;
    opa.dlnxi_lnT=3-dlnkT;
	
	return error;
		
}
		
