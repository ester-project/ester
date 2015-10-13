#include "ester.h"
#include "poly.h"
#include <cstdlib>

/* Solve the rotating polytrope in 2d.
 *
 *	Usage:
 *			polytrope2d polytropic_index omega
 *
 *
 */

int main(int argc, char *argv[]) {

	double n = atof(argv[1]); // Polytropic index (first command-line argument)
	double tol = 1e-12; // Required tolerance
	int nr = 50; // # of radial points
	int nt = 32; // # of points in theta
	int nex = 25; // # of points in external domain
	double omega = atof(argv[2]); // Angular velocity (second command-line argument)


    nr = 10;
    nt = 5;

    printf("Solving 1D Polytropic problem:\n");
    solution *sol_1d = solve_poly1d(n, tol, nr);

    // Build initial solution for 2d
    matrix *guess = new matrix(nr, nt);
    for (int it=0; it<nt; it++) {
        guess->setcol(it, *sol_1d->phi);
    }

    printf("\n\nSolving 2D Polytropic problem with guess:\n");
    solve_poly2d(n, tol, nr, nt, nex, omega, guess, sol_1d->lambda);

    delete(sol_1d);
    delete(guess);

    return 0;
}
