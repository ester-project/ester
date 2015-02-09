#include "ester-config.h"
#include "cesam.h"
#include "utils.h"
#include <iostream>
#include <cstdlib>
#include <stdio.h>


extern "C" {
    void update_comp_(double *T, double *rho, double *comp,
            double *r, int *nr, int *nchim);
}

namespace cesam {

compo::compo(star2d *s, int nchim): nr(s->nr), nchim(nchim) {
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
    compo_matrix = new double[nchim*nr];
    if (compo_matrix == NULL) {
        ester_err("out of memory");
        exit(EXIT_FAILURE);
    }

    matrix remainings = ones(nr, 1);
    int offset = 0;
    for (int i=0; i<nchim-1; i++, offset+=nr) {
        remainings -= c[chims[i]];
        memcpy(compo_matrix+offset,
                c[chims[i]].data(),
                nr * sizeof(*compo_matrix));
    }
    memcpy(compo_matrix+offset,
            remainings.data(),
            nr * sizeof(*compo_matrix));

    return compo_matrix;
}

matrix_map *compo::from_cesam(double *c) {
    matrix_map *cm = new matrix_map();

    int offset = 0;
    for (int i=0; i<nchim-1; i++, offset+=nr) {
        (*cm)[chims[i]].dim(nr, 1);
        memcpy((*cm)[chims[i]].data(),
                c+offset,
                nr * sizeof(*c));
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

int update_comp(double *T, double *rho, double *comp, double *r,
        int nr, int nchim) {

    update_comp_(T, rho, comp, r, &nr, &nchim);
    return 0;
}

}
