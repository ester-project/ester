#ifndef POLY_H
#define POLY_H

class solution {
    public:
        matrix *phi;
        double lambda;
        solution(matrix *phi, double lambda) {
            this->phi = phi;
            this->lambda = lambda;
        }
};

solution *solve_poly1d(double n, double tol, int nr);
solution *solve_poly2d(double n, double tol, int nr, int nt, int nex, double omega,
        matrix *init_phi = NULL, double lambda = 0);

#endif
