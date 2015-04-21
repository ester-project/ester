#include "ester-config.h"
#include "cesam.h"
#include "utils.h"

#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <unistd.h>

extern "C" {
    void update_comp_(double *t0, double *t1, double *rho0, double *rho1,
            double *comp, double *r, int *nconv, double *dt);
    void init_evol_(int *order);
}

namespace cesam {

void init_evol(int order) {
    // We don't want too much output from cesam:
    // redirect stdout (fd = 1) to the file `cesam.log'
    redirect_stdout("cesam.log");

    nuc_cesam_init();
    init_evol_(&order);

    // And restore previous stdout
    restore_stdout();
}

compo::compo(star2d *s, int nchim): nchim(nchim) {
    compo_matrix = NULL;
    chims = new char*[nchim];
    if (nchim < 9) {
        ester_err("nchim < 9");
        exit(EXIT_FAILURE);
    }
    chims[0] = strdup("H");
    chims[1] = strdup("He3");
    chims[2] = strdup("He4");
    chims[3] = strdup("C12");
    chims[4] = strdup("C13");
    chims[5] = strdup("N14");
    chims[6] = strdup("N15");
    chims[7] = strdup("O16");
    chims[8] = strdup("O17");
    chims[9] = strdup("other");
}

double *compo::to_cesam(const composition_map &c) {
    nr = c["H"].nrows();
    compo_matrix = new double[nchim*nr];
    if (compo_matrix == NULL) {
        ester_err("out of memory");
        exit(EXIT_FAILURE);
    }

    for (int i=0; i<nr; i++) {
        double rem = 1.;
        for (int ic=0; ic<nchim-1; ic++) {
            rem -= c[chims[ic]](i);
            compo_matrix[i*nchim+ic] = c[chims[ic]](i);
        }
        compo_matrix[i*nchim+nchim] = rem;
    }

    return compo_matrix;
}

matrix_map *compo::from_cesam(double *c) {
    matrix_map *cm = new matrix_map();

    for (int ic=0; ic<nchim-1; ic++) {
        (*cm)[chims[ic]].dim(nr, 1);
        for (int i=0; i<nr; i++) {
            (*cm)[chims[ic]](i) = compo_matrix[i*nchim+ic];
        }
    }

    return cm;
}

compo::~compo() {
    if (compo_matrix == NULL)
        return;
    delete[] compo_matrix;
    for (int i=0; i<nchim; i++) {
        free(chims[i]);
    }
}

// subroutine update_comp(t0, t1, rho0, rho1, comp, r, nr, nchim, &
//                        nconv, dt)
int update_comp(double *t0, double *t1, double *rho0, double *rho1,
        double *comp, double *r, int nchim, int nconv, double dt) {

    redirect_stdout("cesam.log");
    update_comp_(t0, t1, rho0, rho1, comp, r, &nconv, &dt);
    restore_stdout();
    return 0;
}

}
