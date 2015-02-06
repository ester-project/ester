#include "ester-config.h"
#include "cesam.h"
#include <iostream>
#include <cstdlib>
#include <stdio.h>


extern "C" {
    void update_comp_(double *T, double *rho, double *comp,
            double *r, int *nr, int *nchim);

}

namespace cesam {

int update_comp(double *T, double *rho, double *comp, double *r,
        int nr, int nchim) {

    update_comp_(T, rho, comp, r, &nr, &nchim);
    return 0;
}

}
