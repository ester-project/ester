#ifndef POLY_H
#define POLY_H

class solution {
    public:
        matrix *phi;
        matrix *phi_ex;
        double lambda;
        matrix *Rs;
        inline solution(matrix *phi, matrix *phi_ex, matrix *Rs, double lambda) {
            this->phi = phi;
            this->phi_ex = phi_ex;
            this->Rs = Rs;
            this->lambda = lambda;
        }
};

solution *solve_poly1d(double n, double tol, int nr);
solution *solve_poly2d(double n, double tol, int nr, int nt, int nex, double omega,
        solution *init_guess = NULL);

#endif
