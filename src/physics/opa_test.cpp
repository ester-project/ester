#include<math.h>
#include"matrix.h"
#include"constants.h"
#include"physics.h"

int opa_test(const matrix &T,const matrix &rho,opa_struct &opa) {

	matrix ki,dlnk_lnrho,dlnk_lnT;
	double beta,nu;
	beta=1.975;
	nu=0.138;
	ki=7.1548412e13*pow(rho,nu)*pow(T,-beta);
    opa.k=ki;
    dlnk_lnrho=opa.k/ki*nu;
    dlnk_lnT=-opa.k/ki*beta;    
    opa.xi=16*SIG_SB*pow(T,3)/(3*opa.k*rho);
    opa.dlnxi_lnrho=-1-dlnk_lnrho;
    opa.dlnxi_lnT=3-dlnk_lnT;
    
    return 0;
}

