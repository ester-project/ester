#include "ester-config.h"
#include "cesam.h"
#include <iostream>
#include <cstdlib>


extern "C" {
    void update_comp_(double *T, double *rho, double **comp, double *r,
            double *th, int *nr, int *nth);

}

namespace cesam {

int update_comp(double *T, double *rho, double **comp, double *r, double *th,
        int nr, int nth) {
    update_comp_(T, rho, comp, r, th, &nr, &nth);
    return 0;
}

}
