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


        // Only isotops that exists in non negligible quantities have been considered
        // according to their wikipedia page on the 3rd of July 2023
        eps[0] = comp["H"] / 1.008e0;                     // H
        eps[1] = comp["He3"] / 3. + comp["He4"] / 4.;     // He3 + He4
        eps[2] = comp["C12"] / 12. + comp["C13"] / 13.;   // C12 + C13
        eps[3] = comp["N14"] / 14. + comp["N15"] / 15.;   // N14 + N15
        eps[4] = comp["O16"] / 16. + comp["O17"] / 17.;   // O16 + O17
        eps[5] = comp["Ne20"] / 20. + comp["Ne22"] / 22.; // Ne20 + Ne22
        eps[6] = comp["Na23"] / 23.; // Na23
        eps[7] = comp["Mg24"] / 24 + comp["Mg25"] / 25 + comp["Mg26"] / 26; // Mg24 + Mg25 + Mg26
        eps[8] = comp["Al27"] / 27.0; // Al27
        eps[9] = comp["Si28"] / 28.0; // Si28
        eps[10] = comp["P31"] / 31.0; // P31
        eps[11] = comp["S32"] / 32.0; // S32
        eps[12] = comp["Cl35"] / 35.0 + comp["Cl37"] / 37.0; // Cl35 + Cl37
        eps[13] = comp["A40"] / 40.0; // A40
        eps[14] = comp["Ca40"] / 40.0; // Ca40
        eps[15] = comp["Ti"] / 47.867; // Ti //there are too many isotops
        eps[16] = comp["Cr"] / 52.0; // Cr //there are too many isotops
        eps[17] = comp["Mn55"] / 55.0; // Mn
        eps[18] = comp["Fe"] / 55.845; // Fe //there are too many isotops
        eps[19] = comp["Ni"] / 58.693; // Ni //there are too many isotops

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
