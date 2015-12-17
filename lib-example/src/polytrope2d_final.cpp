#include "ester.h"
#include "poly.h"
#include <cstdlib>

/* Solve the rotating polytrope in 2d.
 *
 *	Usage:
 *			polytrope2d_final polytropic_index omega nr nt
 *
 *
 */

int main(int argc, char *argv[]) {

	double n = atof(argv[1]); // Polytropic index (first command-line argument)
	double tol = 1e-12; // Required tolerance
	int nr = atof(argv[3]); // # of radial points
	int nt = atof(argv[4]); // # of points in theta
	int nex = 25; // # of points in external domain
	double omega = atof(argv[2]); // Angular velocity (second command-line argument)


    printf("Solving 1D Polytropic problem:\n");
    solution *sol_1d = solve_poly1d(n, tol, nr);

    // Build initial solution for 2d
    matrix *guess = new matrix(nr, nt);
    for (int it=0; it<nt; it++) {
        guess->setcol(it, *sol_1d->phi);
    }
    sol_1d->phi = guess;

    printf("\n\nSolving 2D Polytropic problem with guess:\n");
    if (omega <= 0.4) solve_poly2d(n, tol, nr, nt, nex, omega, sol_1d);
    if (omega > 0.4) {
      solution *sol_2d = solve_poly2d(n, tol, nr, nt, nex, 0.4, sol_1d);
      printf("\n\nSolving 2D Polytropic problem with (2D) guess:\n");
      solve_poly2d(n, tol, nr, nt, nex, omega, sol_2d);
    }

    delete(sol_1d);
    delete(guess);

    return 0;
}
