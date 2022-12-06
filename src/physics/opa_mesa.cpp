#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "utils.h"
#include "matrix.h"
#include "constants.h"
#include "physics.h"
#include <cmath>
#include <iostream>
using namespace std;

extern"C" {
    int init_mesa_opa_();
    void opa_mesa_(double *, double *,
            double *, double *, double *, double *, double *);
}

int opa_mesa(const matrix& X, double Z, const matrix& T, const matrix& rho,
		opa_struct& opa) {

    int error = 0;
    static bool init = false;

	opa.k.dim(T.nrows(), T.ncols());
    opa.xi.dim(T.nrows(), T.ncols());
    opa.dlnxi_lnrho.dim(T.nrows(), T.ncols());
    opa.dlnxi_lnT.dim(T.nrows(), T.ncols());

    if (!init) {
        init_mesa_opa_();
    }

    for (int i=0; i<X.nrows(); i++) {
        for (int j=0; j<X.ncols(); j++) {
            double x[10];
            double t, ro, kap, dlnkap_dlnT, dlnkap_dlnd, dkapx;
//             cout << "i= : " << i << endl;
//             cout << "j= : " << j << endl;        
//              x[0] = X(i, j);
//              x[1] = 1.0 - X(i, j) - Z;
//              x[2] = Z*(1.71243418737847e-01 + 2.06388003380057e-03);
//              x[3] = Z*(5.29501630871695e-02 + 2.08372414940812e-04);
//              x[4] = Z*(4.82006487350336e-01 + 1.95126448826986e-04);
//              
//              
//              
             
             
             x[0] = 0.0000000;
             x[1] = X(i, j);
             x[2] = 3.15247417638132e-04*(1.-X(i, j)-Z);
             x[3] = 1.0 - X(i, j) - Z-x[2];
             x[4] = Z*1.71243418737847e-01;
             x[5] = Z*2.06388003380057e-03;
             x[6] = Z*5.29501630871695e-02;
             x[7] = Z*2.08372414940812e-04;
             x[8] = Z*4.82006487350336e-01;
             x[9] = Z*1.95126448826986e-04;

             
//              x[2] = Z*(1.71243418737847e-01 + 2.06388003380057e-03);
//              x[3] = Z*(5.29501630871695e-02 + 2.08372414940812e-04);
//              x[4] = Z*(4.82006487350336e-01 + 1.95126448826986e-04);
             
             
             
//              x[5] = 1.0 - X(i, j) - Z;
//              x[6] = 1.0 - X(i, j) - Z;
             
//             x[2] = Z;
            // x[3] = 3.425E-03;
            // x[4] = 4.128E-05;
            // x[5] = 1.059E-03;
            // x[6] = 4.167E-06;
            // x[7] = 9.640E-03;
            // x[8] = 3.903E-06;
            // x[9] = 5.827E-03;

            t = T(i, j);
            ro = rho(i, j);

            opa_mesa_(x, &t, &ro, &kap, &dlnkap_dlnT, &dlnkap_dlnd, &dkapx);
//             cout << "after opa_mesa_ " << endl;
//             cout << "opa_mesa  kap =" << kap << endl;
//             cout << "opa_mesa  ro =" << ro << endl;
//             cout << "opa_mesa  SIG_SB =" << SIG_SB << endl;
//             cout << "opa_mesa  t =" << t << endl;
//             cout << "opa_mesa  dkapro =" << dkapro << endl;

            opa.k(i, j) = kap;
            opa.xi(i, j) = 16*SIG_SB*pow(t,3)/(3*kap*ro);
            opa.dlnxi_lnrho(i, j) = -1.0 - dlnkap_dlnd;
//             opa.dlnxi_lnrho(i, j) = -1.0 - ro/kap * 0.1;
            opa.dlnxi_lnT(i, j) = 3.0 - dlnkap_dlnT; 
//             opa.dlnxi_lnT(i, j) = 3.0 - t/kap * 0.1; 
            
//             cout << "opa_mesa  opa.k(i, j) =" << opa.k(i, j) << endl;
//             cout << "opa_mesa  opa.xi(i, j) =" << opa.xi(i, j) << endl;
//             cout << "opa_mesa  opa.dlnxi_lnrho(i, j) =" << opa.dlnxi_lnrho(i, j) << endl;
//             cout << "opa_mesa  opa.dlnxi_lnT(i, j) =" << opa.dlnxi_lnT(i, j) << endl;
//             cout << "opa_mesa  SIG_SB =" << SIG_SB << endl;
//             cout << "opa_mesa  t =" << t << endl;
//             cout << "opa_mesa  dkapro =" << dkapro << endl;
            

        }
    }

	return error;

}

