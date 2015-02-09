#ifndef _CESAM_H_
#define _CESAM_H_

#include "star.h"

namespace cesam {

class compo {
    private:
        double *compo_matrix;
        int nr, nchim;
        char **chims;

    public:
        compo(star2d *s, int nchim);
        ~compo();
        double *to_cesam(const composition_map &c);
        matrix_map *from_cesam(double *c);
};

int update_comp(double *T, double *rho, double *comp, double *r,
        int nr, int nchim);
}

#endif //_CESAM_H_
