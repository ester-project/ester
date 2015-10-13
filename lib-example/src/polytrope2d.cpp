#include "ester.h"
#include "poly.h"
#include <cstdlib>

int main(int argc, char *argv[]) {

    double n = atof(argv[1]); // Polytropic index (first command-line argument)
    double tol = 1e-12; // Required tolerance
    int nr = atof(argv[3]); // # of radial points
    int nt = atof(argv[4]); // # of points in theta
    int nex = 25; // # of points in external domain
    double omega = atof(argv[2]); // Angular velocity (second command-line argument)

    solve_poly2d(n, tol, nr, nt, nex, omega);

    return 0;
}
