#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "utils.h"
#include "matrix.h"
#include "constants.h"
#include "physics.h"
#include <cmath>

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

        double T6_debug = T(i)*1e-6;
        double R_opal   = rho(i)/(T6_debug*T6_debug*T6_debug);
        double logT     = log10(T(i));
        double logR     = log10(R_opal);

        if (i >= 315 && i <= 325) {
            printf(
                "OPAL_DIAG "
                "i=%4d  "
                "T=% .16e  "
                "rho=% .16e  "
                "T6=% .16e  "
                "R_OPAL=% .16e  "
                "logT=% .16e  "
                "logR=% .16e  "
                "opact=% .16e\n",
                i,
                T(i),
                rho(i),
                T6_debug,
                R_opal,
                logT,
                logR,
                opa.k(i)
            );
        }

        if (opa.k(i) == -99) {

            printf("\nValues outside OPAL opacity table\n");
            printf("  i        = %d\n", i);
            printf("  X        = %.16e\n", X(i));
            printf("  Z        = %.16e\n", Z);
            printf("  T        = %.16e K\n", T(i));
            printf("  rho      = %.16e g/cm^3\n", rho(i));
            printf("  T6       = %.16e\n", T6_debug);
            printf("  R_OPAL   = %.16e\n", R_opal);
            printf("  logT     = %.16e\n", logT);
            printf("  logR     = %.16e\n", logR);

            print_stack();
            exit(EXIT_FAILURE);
        }

        if (error) {

            double T6_debug = T(i)*1e-6;
            double R_debug  = rho(i)/(T6_debug*T6_debug*T6_debug);

            printf("Values outside OPAL opacity table\n");
            printf("  X   = %e\n", X(i));
            printf("  Z   = %e\n", Z);
            printf("  T   = %e\n", T(i));
            printf("  rho = %e\n", rho(i));
            printf(" i = %i\n",i);
            printf("  R_OPAL   = %e\n", R_debug);
            print_stack();
            exit(EXIT_FAILURE);
        }
	}
	opa.k=pow(10,opa.k);
	dlnkT-=3*dlnkrho;
	opa.xi=16*SIG_SB*pow(T,3)/(3*opa.k*rho);
	opa.dlnxi_lnrho=-1-dlnkrho;
    opa.dlnxi_lnT=3-dlnkT;

	return error;

}

