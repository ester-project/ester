#include "ester.h"
// #include "matplotlib.h"

#include <cstdlib>

matrix solve_poly1d(double n, double tol, int nr, double hsurf) {
    // Create a mapping
    mapping map;
    map.set_ndomains(1);
    map.set_npts(nr);
    map.gl.set_xif(0., 1.); // This changes the values of zeta of the limits of the domain directly (do it before map.init())
    map.set_nt(1); // 1d
    map.init();
    /* Instead of changing the values of zeta directly, we can also do (after map.init()):

       map.R.setrow(0, 0.*ones(1, 1)); // or map.R(0) = 0; (only because nth = 1)
       map.R.setrow(1, 1.*ones(1, 1)); // or map.R(1) = 1;
       map.remap();

       This is more useful for 2d maps.
       Note that in this case this is not really necessary to change the interval as the default is already (0, 1)
       */

    // Create a symbolic object for the equation with 3 variables: Phi, Lambda and Phi0
    // Phi0 and Lambda are scalar variables but it doesn't matter for the symbolic object
    symbolic S;
    S.set_map(map);
    sym symPhi = S.regvar("Phi");
    sym symLambda = S.regvar("Lambda");
    sym symPhi0 = S.regvar("Phi0");

    sym eq = lap(symPhi) - pow(1 - symLambda * (symPhi-symPhi0), n); // This is our equation

    // Create numerical variables for the solution with the initial guesses
    matrix Phi = map.r*map.r;
    double Lambda = 1.;
    double Phi0 = 0.;

    double error = 1.;
    int it = 0;

    // Create a solver for the three variables. The size of each variable is determined by the
    // solver object based on the size of the equations that we will define for them.
    // The "Phi" equation will have nr x 1 points while for "Phi0" and "Lambda" we will introduce
    // only boundary conditions, so the resulting equations will have 1 x 1 size
    solver op;
    op.init(1, 3, "full");
    op.regvar("Phi");
    op.regvar("Lambda");
    op.regvar("Phi0");
    op.set_nr(map.npts);

    while(error>tol && it<10000) {
        // Put the current values of variables in the symbolic object
        S.set_value("Phi", Phi);
        S.set_value("Lambda", Lambda*ones(1, 1)); // Note that the assigned value should be of type matrix
        S.set_value("Phi0", Phi0*ones(1, 1));

        op.reset(); // Delete the equations of the previous iteration

        // Define the equation for "Phi", here we will use the symbolic object to automatically calculate
        // the required terms
        eq.add(&op, "Phi", "Phi"); // Add the jacobian of eq with respect to "Phi"
        eq.add(&op, "Phi", "Lambda"); // Add the jacobian of eq with respect to "Lambda"
        eq.add(&op, "Phi", "Phi0"); // Add the jacobian of eq with respect to "Phi0"

        // Add the boundary conditions
        op.bc_bot2_add_l(0, "Phi", "Phi", ones(1, 1), map.D.block(0).row(0));
        op.bc_top1_add_l(0, "Phi", "Phi", ones(1, 1), map.D.block(0).row(-1));
        op.bc_top1_add_d(0, "Phi", "Phi", ones(1, 1));

        // RHS for "Phi"
        matrix rhs = -eq.eval();
        rhs(0) = -(map.D, Phi)(0);
        rhs(-1) = -(map.D, Phi)(-1) - Phi(-1);
        op.set_rhs("Phi", rhs);

        // Equation for "Phi0":    dPhi(0) - dPhi0 = -( Phi(0) - Phi0 )
        // We add it as a boundary condition at the bottom of the domain. This way, it uses the value of "Phi" at r = 0
        op.bc_bot2_add_d(0, "Phi0", "Phi", ones(1, 1));
        op.bc_bot2_add_d(0, "Phi0", "Phi0", -ones(1, 1));
        op.set_rhs("Phi0", -(Phi(0) - Phi0) * ones(1, 1));


        // Equation for "Lambda", recall that Lambda = 1./(Phi(1)-Phi(0)), so we will use the equation
        // 		Lambda * (dPhi(1) - dPhi0) + dLambda * ( Phi(1)-Phi0 ) = -( Lambda*(Phi(1)-Phi0) - 1)
        // We add it as a boundary condition at the top of the domain. This way, it uses the value of "Phi" at r = 1
        // ("Phi0" is defined in the whole domain)
        op.bc_top1_add_d(0, "Lambda", "Phi", Lambda*ones(1, 1));
        op.bc_top1_add_d(0, "Lambda", "Phi0", -Lambda*ones(1, 1));
        op.bc_top1_add_d(0, "Lambda", "Lambda", (Phi(-1)-Phi0)*ones(1, 1));
        op.set_rhs("Lambda", -(Lambda*(Phi(-1)-Phi0) -1) * ones(1, 1));

        op.solve(); // Solve the equations

        matrix dPhi = op.get_var("Phi");
        error = max(abs(dPhi));  // Calculate the error (absolute)

        double relax = 1.;
        if(error>0.01) relax = 0.2; // Relax if error too large

        // Update variables
        Phi += relax*dPhi;
        Phi0 += relax*op.get_var("Phi0")(0);
        Lambda += relax*op.get_var("Lambda")(0);

        it++;
    }

    if(error>tol) {
        ester_err("No converge\n");
    }

    matrix h = 1.0 - (Phi - Phi(0))*Lambda;
    double ri = .5;
    double hi = map.gl.eval(h, ri)(0);
    matrix dh = (map.D, h);
    double dhi = map.gl.eval(dh, ri)(0);

    int nit = 0;
    int maxit = 100;
    // Newton method to find ri such that: h(ri) = hsurf
    while (fabs(hi - hsurf) > 1e-12 && ++nit < maxit) {
        ri += (hsurf - hi)/dhi;
        hi = map.gl.eval(h, ri)(0);
        dhi = map.gl.eval(dh, ri)(0);
    }
    if (nit >= maxit) {
        ester_err("No convergence finding polytrope surface such that h(surface) = %e", hsurf);
    }

    matrix h2 = map.gl.eval(h, map.r*ri);

    // printf("h(r=%e) = %e\n", ri, hi);

    return h2;
}

