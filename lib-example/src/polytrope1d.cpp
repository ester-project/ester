#include "ester.h"
#include "poly.h"
#include <cstdlib>

int main(int argc, char *argv[]) {

    double n = 1.5; // Polytropic index
    double tol = 1e-12; // Required tolerance
    int nr = 50; // # of points

    solve_poly1d(n, tol, nr);

    return 0;
}
