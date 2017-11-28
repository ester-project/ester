#include "ester-config.h"
#include <cmath>
#include <string.h>
#include "matrix.h"
#include "constants.h"
#include "physics.h"

extern"C" {
    int init_cesam_opa_();
    void opa_cesam_(double *, double *,
            double *, double *, double *, double *, double *);
}

int opa_cesam(const matrix& X, double Z, const matrix& T, const matrix& rho,
		opa_struct& opa) {

    int error = 0;
    static bool init = false;

	opa.k.dim(T.nrows(), T.ncols());
    opa.xi.dim(T.nrows(), T.ncols());
    opa.dlnxi_lnrho.dim(T.nrows(), T.ncols());
    opa.dlnxi_lnT.dim(T.nrows(), T.ncols());

    if (!init) {
        init_cesam_opa_();
    }

    for (int i=0; i<X.nrows(); i++) {
        for (int j=0; j<X.ncols(); j++) {
            double x[3];
            double t, ro, kap, dkapt, dkapro, dkapx;

            x[0] = X(i, j);
            x[1] = 1.0 - X(i, j) - Z;
            x[2] = Z;
            // x[3] = 3.425E-03;
            // x[4] = 4.128E-05;
            // x[5] = 1.059E-03;
            // x[6] = 4.167E-06;
            // x[7] = 9.640E-03;
            // x[8] = 3.903E-06;
            // x[9] = 5.827E-03;

            t = T(i, j);
            ro = rho(i, j);

            opa_cesam_(x, &t, &ro, &kap, &dkapt, &dkapro, &dkapx);

            opa.k(i, j) = kap;
            opa.xi(i, j) = 16*SIG_SB*pow(t,3)/(3*kap*ro);
            opa.dlnxi_lnrho(i, j) = -1.0 - ro/kap * dkapro;
            opa.dlnxi_lnT(i, j) = 3.0 - t/kap * dkapt;

        }
    }

	return error;

}
