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

    AbundanceMap& abundance_map = global_abundance_map; // using atomic weight values

    // Loading initial composition from text files before the loop. 
    // There is still room to improve as this is called on every iteration
    // Most efficient would be to load this only once at the beginning of the run - PH --> has been done - MG 

    //CompositionData initial_comp = parse_composition_data();
	CompositionData initial_comp = global_abundance_map.comp_data; // already been made in initilisation

    for (int i=0; i<N; i++) {

        double_map comp = update_initial_composition(initial_comp, X(i), Z);

        eps[0] = comp["H"]/abundance_map.A_weights["H"] ; // H1 + H2
        eps[1] = (comp["He3"] + comp["He4"]) / abundance_map.A_weights["He"];   // He3 + He4
        eps[2] = (comp["C12"] + comp["C13"]) / abundance_map.A_weights["C"]; // C12 + C13
        eps[3] = (comp["N14"] + comp["N15"]) / abundance_map.A_weights["N"]; // N14 + N15
        eps[4] = (comp["O16"] + comp["O17"] + comp["O18"]) / abundance_map.A_weights["O"]; // O16 + O17 + O18 
        eps[5] = (comp["Ne20"] + comp["Ne21"] + comp["Ne22"]) / abundance_map.A_weights["Ne"]; // Ne
        eps[6] = comp["Na23"] / abundance_map.A_weights["Na"]; // Na
        eps[7] = (comp["Mg24"] + comp["Mg25"] + comp["Mg26"]) / abundance_map.A_weights["Mg"]; // Mg
        eps[8] = comp["Al27"] / abundance_map.A_weights["Al"]; // AL
        eps[9] = (comp["Si28"] + comp["Si29"] + comp["Si30"]) / abundance_map.A_weights["Si"]; // Si
        eps[10] = comp["P31"]/abundance_map.A_weights["P"]; // P
        eps[11] = (comp["S32"] + comp["S33"] + comp["S34"] + comp["S36"]) / abundance_map.A_weights["S"]; // S
        eps[12] = (comp["Cl35"] + comp["Cl37"]) / abundance_map.A_weights["Cl"]; // Cl
        eps[13] = (comp["Ar36"] + comp["Ar38"] + comp["Ar40"]) / abundance_map.A_weights["Ar"]; // Ar
        eps[14] = (comp["Ca40"] + comp["Ca42"] + comp["Ca43"] + comp["Ca44"] + comp["Ca46"] + comp["Ca48"]) / abundance_map.A_weights["Ca"]; // Ca
        eps[15] = (comp["Ti46"] + comp["Ti47"] + comp["Ti48"] + comp["Ti49"] + comp["Ti50"]) / abundance_map.A_weights["Ti"]; // Ti
        eps[16] = (comp["Cr50"] + comp["Cr52"] + comp["Cr53"] + comp["Cr54"]) / abundance_map.A_weights["Cr"]; // Cr
        eps[17] = comp["Mn55"] / abundance_map.A_weights["Mn"]; // Mn
        eps[18] = (comp["Fe54"] + comp["Fe56"] + comp["Fe57"] + comp["Fe58"]) / abundance_map.A_weights["Fe"]; // Fe
        eps[19] = (comp["Ni58"] + comp["Ni60"] + comp["Ni61"] + comp["Ni62"] + comp["Ni64"]) / abundance_map.A_weights["Ni"]; // Ni	

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
    }
    delete[] eps;

    return 0;
}
