#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "utils.h"
#include "physics.h"
#include "constants.h"
#include <string.h>
#include <cmath>

int atm_onelayer(const matrix &X,double Z,const matrix &g,const matrix &Teff,
		const char *eos_name,const char *opa_name,atm_struct &atm) {

	int n=g.nrows();
	int m=g.ncols();

	atm.ps.dim(n,m);atm.Ts.dim(n,m);
	atm.dlnps_lng.dim(n,m);atm.dlnps_lnTeff.dim(n,m);
	atm.dlnTs_lng.dim(n,m);atm.dlnTs_lnTeff.dim(n,m);
	
	atm.Ts=Teff;
	atm.dlnTs_lnTeff=ones(n,m);
	atm.dlnTs_lng=zeros(n,m);
	
	for(int i=0;i<n;i++) {
		for(int j=0;j<m;j++) {
			double logps,logg;
			matrix p,T,rho;
			eos_struct eos;
			opa_struct opa;
	
			logps=4.2; //Initial value
			logg=log10(g(i,j));
			T=Teff(i,j)*ones(1,1);
			strcpy(eos.name,eos_name);
			strcpy(opa.name,opa_name);
			int fin=0;
			int nit=0;
			while(fin<2) {
				nit++;
				if(nit>100) {
					ester_critical("(atm_onelayer) No convergence");
					return 1; // Useless because ester_critical already exit
				}
				double F,dF,dlogps;
				p=pow(10,logps)*ones(1,1);
				if(eos_calc(X(i,j)*ones(1,1),Z,T,p,rho,eos)) return 1;
				if(opa_calc(X(i,j)*ones(1,1),Z,T,rho,opa)) return 1;
				F=logps+log10(opa.k)(0)-logg-log10(2./3.);
				dF=1.-(1+opa.dlnxi_lnrho(0))/eos.chi_rho(0);
				dlogps=-F/dF;
				if (fabs(dlogps)<1e-8) fin++;
				logps+=dlogps;
			}
			atm.ps(i,j)=pow(10,logps);
			double kp,kT;
			kp=-(1+opa.dlnxi_lnrho(0))/eos.chi_rho(0);
			kT=3-opa.dlnxi_lnT(0)+(1+opa.dlnxi_lnrho(0))*eos.d(0);
			atm.dlnps_lng(i,j)=1./(1+kp);
			atm.dlnps_lnTeff(i,j)=-kT/(1+kp);
		}
	}
	
	return 0;
}
