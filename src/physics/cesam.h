#ifndef _CESAM_H_
#define _CESAM_H_

#include "star.h"

extern int nuc_cesam_init();

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

int update_comp(double *t0, double *t1, double *rho0, double *rho1,
        double *comp, double *r, int nchim, int nconv, double dt);
void init_evol(int order);
}

#endif //_CESAM_H_
