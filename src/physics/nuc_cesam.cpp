#include "ester-config.h"
#include "physics.h"

static bool init = false;

extern "C" {
    void nuc_cesam_init_();
    void nuc_cesam_set_time_step_(int *nyears);
    void nuc_cesam_init_abon_(double *X, double *Z, double *comp);
    void nuc_cesam_eps_(double *t, double *ro, double *comp,
            double *epsilon, double *et, double *ero, double *ex);
    void nuc_cesam_dcomp_(double *t, double *ro, double *comp,
            double *dcomp, double *jac, int *dp);
}

void nuc_cesam_init() {
    int ny = 1000;  // time step: 1000 years
    nuc_cesam_init_();
    double XX = 0.7, ZZ = 0.02;
    matrix ab(10);
    nuc_cesam_init_abon_(&XX, &ZZ, ab.data());
    init = true;
    nuc_cesam_set_time_step_(&ny);
}

matrix nuc_cesam_abon(double_map comp) {

    matrix ab(10);

    ab(0) = comp["H"];
    ab(1) = comp["He3"];
    ab(2) = comp["He4"];
    ab(3) = comp["C12"];
    ab(4) = comp["C13"];
    ab(5) = comp["N14"];
    ab(6) = comp["N15"];
    ab(7) = comp["O16"];
    ab(8) = comp["O17"];
    ab(9) = 1. - comp["H"] - comp["He3"] - comp["He4"]
        - comp["C12"] - comp["C13"] - comp["N14"]
        - comp["N15"] - comp["O16"] - comp["O17"];
    return ab;
}

double_map nuc_cesam_abon(matrix ab) {

    double_map comp;

    comp["H"] = ab(0);
    comp["He3"] = ab(1);
    comp["He4"] = ab(2);
    comp["C12"] = ab(3);
    comp["C13"] = ab(4);
    comp["N14"] = ab(5);
    comp["N15"] = ab(6);
    comp["O16"] = ab(7);
    comp["O17"] = ab(8);

    return comp;
}

void nuc_cesam_init_jac(composition_map &comp) {
    composition_map::iterator it;

    for (it = comp.begin(); it != comp.end(); it++)
        comp.jac[it->first] = 0 * comp;
}

void nuc_cesam_jac(composition_map &comp, const matrix &J, int i, int j) {

    std::map < std::string, matrix_map > &jac = comp.jac;

    jac["H"](i, j) = nuc_cesam_abon(J.row(0));
    jac["He3"](i, j) = nuc_cesam_abon(J.row(1));
    jac["He4"](i, j) = nuc_cesam_abon(J.row(2));
    jac["C12"](i, j) = nuc_cesam_abon(J.row(3));
    jac["C13"](i, j) = nuc_cesam_abon(J.row(4));
    jac["N14"](i, j) = nuc_cesam_abon(J.row(5));
    jac["N15"](i, j) = nuc_cesam_abon(J.row(6));
    jac["O16"](i, j) = nuc_cesam_abon(J.row(7));
    jac["O17"](i, j) = nuc_cesam_abon(J.row(8));
}

int nuc_cesam(const composition_map &comp, const matrix &T, const matrix &rho,
        nuc_struct &nuc) {

    if (!init)
        nuc_cesam_init();

    nuc.eps = zeros(T.nrows(), T.ncols());
    nuc.pp = zeros(T.nrows(), T.ncols());
    nuc.cno = zeros(T.nrows(), T.ncols());
    nuc.dlneps_lnrho = zeros(T.nrows(), T.ncols());
    nuc.dlneps_lnT = zeros(T.nrows(), T.ncols());
    matrix epsilon(4), ex(10);
    double et, ero, t, ro;
    for (int j = 0; j < T.ncols(); j++) {
        for (int i = 0; i < T.nrows(); i++) {
            matrix ab;
            ab = nuc_cesam_abon(comp(i, j));
            t = T(i, j);
            ro = rho(i, j);
            nuc_cesam_eps_(&t, &ro, ab.data(), epsilon.data(),
                    &et, &ero, ex.data());
            nuc.eps(i, j) = epsilon(0);
            nuc.pp(i, j) = epsilon(1);
            nuc.cno(i, j) = epsilon(2);
            if (epsilon(0) == 0) {
                nuc.dlneps_lnrho(i, j) = 0;
                nuc.dlneps_lnT(i, j) = 0;
            }
            else {
                nuc.dlneps_lnrho(i, j) = ero / epsilon(0) * rho(i, j);
                nuc.dlneps_lnT(i, j) = et / epsilon(0) * T(i, j);
            }
        }
    }
    return 0;
}

int nuc_cesam_dcomp(composition_map &comp, const matrix &T, const matrix &rho,
        nuc_struct &nuc) {

    if (!init)
        nuc_cesam_init();

    comp.dt = 0 * comp;
    nuc_cesam_init_jac(comp);
    double t, ro;
    for (int j = 0; j < T.ncols(); j++) {
        for (int i = 0; i < T.nrows(); i++) {
            int doprint = 0;
            matrix ab, dab(10), jac(10, 10);
            ab = nuc_cesam_abon(comp(i, j));
            t = T(i, j);
            ro = rho(i, j);
            if (j == 0 && i == 0) {
                doprint = 1;
            }
            nuc_cesam_dcomp_(&t, &ro, ab.data(), dab.data(), jac.data(),
                    &doprint);

            comp.dt(i, j) = nuc_cesam_abon(dab);
            nuc_cesam_jac(comp, jac, i, j);

            nuc_cesam_dcomp_(&t, &ro, ab.data(), dab.data(), jac.data(),
                    &doprint);
        }
    }
    // As part of the mass is transformed in energy (and neutrinos), the
    // temporal derivatives given by cesam doesn't verify Sum(dXi/dt)=0. We need
    // to make a correction in this derivatives to ensure that sum(Xi)=1 during
    // the evolution.
    // Let Xi be the normalized abundances and Xi' the non-normalized ones:
    //            Xi=Xi'/sum(Xi')
    // then
    //                      dXi/dt=1/sum(Xi')*dXi/dt-Xi'/sum(Xi')^2*sum(dXi'/dt)
    // dXi'/dt are the derivatives given by cesam. If at a given step Xi'=Xi
    // then the correction at this step reads:
    //                      dXi/dt=dXi'/dt-Xi*sum(dXi'/dt)
    // where the new values verify sum(dXi/dt)=0
    // 
    // Note that, after this correction, even for elements that doesn't
    // participate in the reaction we will have dXi/dt!=0. This is a consequence
    // of the change in the total mass caused by the reaction.

    matrix dXtot;
    dXtot = comp.dt.sum();
    composition_map::iterator it, it2;
    std::map < std::string, matrix_map > Jnew;
    matrix_map Jsum;
    for (it = comp.begin(); it != comp.end(); it++) {
        Jsum[it->first] = zeros(T.nrows(), T.ncols());
        for (it2 = comp.begin(); it2 != comp.end(); it2++) {
            Jsum[it->first] += comp.jac[it2->first][it->first];
        }
    }
    for (it = comp.begin(); it != comp.end(); it++) {
        for (it2 = comp.begin(); it2 != comp.end(); it2++) {
            Jnew[it->first][it2->first] =
                comp.jac[it->first][it2->first] -
                comp[it->first] * Jsum[it2->first];
        }
        Jnew[it->first][it->first] -= dXtot;
    }
    comp.jac = Jnew;
    comp.dt -= comp * dXtot;

    return 0;
}

