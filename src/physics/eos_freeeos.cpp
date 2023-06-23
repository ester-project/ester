#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "utils.h"
#include "physics.h"
#include "constants.h"

#include <iostream>
#include <cmath>

extern "C" {
    void free_eos_(int *ifoption, int *ifmodified, int *ifion, int *kif_in,
            double *eps, int *neps, double *match_variable, double *tl, double*
            fl, double *t, double *rho, double *rl, double *p, double *pl,
            double *cf, double *cp, double *qf, double *qp, double *sf, double
            *st, double *grada, double *rtp, double *rmue, double *fh2, double
            *fhe2, double *fhe3, double *xmu1, double *xmu3, double *eta, double
            *gamma1, double *gamma2, double *gamma3, double *h2rat, double
            *h2plusrat, double *lambda, double *gamma_e, double *degeneracy,
            double *pressure, double *density, double *energy, double *enthalpy,
            double *entropy, int *iteration_count);
}

int eos_freeeos(const matrix &X, double Z, const matrix &T, const matrix &p,
        matrix &rho, eos_struct &eos) {

    double t;
    int ifoption = 1;
    int ifmodified = 2;
    int ifion = 0;
    int kif_in = 1;
    int neps = 20;
    double *eps = new double[neps];
    double tl;
    double fl;
    double rl;
    double pl;
    double cf;
    double cp;
    double qf;
    double qp;
    double sf;
    double st;
    double grada;
    double rtp;

    double rmue;
    double fh2;
    double fhe2;
    double fhe3;
    double xmu1;
    double xmu3;
    double eta;

    double gamma1;
    double gamma2;
    double gamma3;
    double h2rat;
    double h2plusrat;
    double lambda;
    double gamma_e;

    double degeneracy[3];
    double pressure[3];
    double density[3];
    double energy[3];
    double enthalpy[3];
    double entropy[3];

    int iteration_count;

    int N=T.nrows()*T.ncols();

    rho.dim(T.nrows(), T.ncols());
    eos.s.dim(T.nrows(), T.ncols());
    eos.G1.dim(T.nrows(), T.ncols());
    eos.del_ad.dim(T.nrows(), T.ncols());
    eos.G3_1.dim(T.nrows(), T.ncols());
    eos.d.dim(T.nrows(), T.ncols());
    eos.cp.dim(T.nrows(), T.ncols());
    eos.cv.dim(T.nrows(), T.ncols());
    eos.chi_rho.dim(T.nrows(), T.ncols());
    eos.chi_T.dim(T.nrows(), T.ncols());
    eos.prad.dim(T.nrows(), T.ncols()); // added MR june 2023

    for (int i=0; i<N; i++) {

        double_map comp = initial_composition(X(i), Z);

        eps[0] = comp["H"] / 1.008e0;                       // H
        eps[1] = (comp["He3"] + comp["He4"]) / 4.0026e0;   // He3 + He4
        eps[2] = (comp["C12"] + comp["C13"]) / 12.0111e0; // C12 + C13
        eps[3] = (comp["N14"] + comp["N15"]) / 14.0067e0; // N14 + N15
        eps[4] = (comp["O16"] + comp["O17"]) / 15.9994e0; // O16 + O17
        eps[5] = 0.0; // Ne
        eps[6] = 0.0; // Na
        eps[7] = 0.0; // Mg
        eps[8] = 0.0; // AL
        eps[9] = 0.0; // Si
        eps[10] = .0; // P
        eps[11] = .0; // S
        eps[12] = .0; // Cl
        eps[13] = .0; // A
        eps[14] = .0; // Ca
        eps[15] = .0; // Ti
        eps[16] = .0; // Cr
        eps[17] = .0; // Mn
        eps[18] = .0; // Fe
        eps[19] = .0; // Ni

        double pi = p(i);
        double match_variable = log(pi);
        double rhoi;
        t = T(i);
        tl = log(t);
        pl = log(pi);
        free_eos_(&ifoption, &ifmodified,
                &ifion, &kif_in, eps, &neps, &match_variable, &tl, &fl,
                &t, &rhoi, &rl, &pi, &pl, &cf, &cp, &qf, &qp, &sf, &st, &grada, &rtp,
                &rmue, &fh2, &fhe2, &fhe3, &xmu1, &xmu3, &eta,
                &gamma1, &gamma2, &gamma3, &h2rat, &h2plusrat, &lambda, &gamma_e,
                degeneracy, pressure, density, energy, enthalpy, entropy,
                &iteration_count);
        rho(i) = rhoi;
        if (iteration_count < 0) {
            ester_err(
                    "Values outside freeEOS eos table:\n"
                    "  X = %e\n"
                    "  Z = %e\n"
                    "  T = %e\n"
                    "  p = %e", X(i), Z, t, p(i));
        }
        eos.s(i) = entropy[0];
        eos.G1(i) = gamma1;
        eos.del_ad(i) = grada;
        eos.G3_1(i) = gamma1*(gamma2-1.0)/gamma2;
        eos.d(i) = -density[2];          // -d(lnRho)/d(lnT)
        eos.cp(i) = cp;
    	// eos.cv(i)=1e6*(*(eeos_.eos+4));
        eos.cv(i) = energy[2] * (1.0/t); // dE/dT (energy[2] is dE/dlnT)
        eos.chi_rho(i) = 1.0/density[1];   // dlogP/dlogRho
        eos.chi_T(i) = -density[2] / density[1];     // dlogP/dlogT
        eos.prad(i)=A_RAD/3*pow(t,4); // added MR june 2023
    }
    delete[] eps;

    return 0;
}
