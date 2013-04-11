#include<cmath>
#include"matrix.h"
#include"constants.h"
#include"physics.h"

int opa_kramer(const matrix &T,const matrix &rho,opa_struct &opa) {

	matrix ke,ki,dlnk_lnrho,dlnk_lnT;
	ki=7.1548412e13*pow(rho,0.138316)*pow(T,-1.97541);
	ke=1.6236784e-33*pow(rho,0.407895)*pow(T,9.28289);
    opa.k=1/(1/ki+0/ke);
    dlnk_lnrho=opa.k/ki*0.138316+0*opa.k/ke*0.407895;
    dlnk_lnT=-opa.k/ki*1.97541+0*opa.k/ke*9.28289;      
    
    opa.xi=16*SIG_SB*pow(T,3)/(3*opa.k*rho);
    opa.dlnxi_lnrho=-1-dlnk_lnrho;
    opa.dlnxi_lnT=3-dlnk_lnT;
    
    return 0;
}

