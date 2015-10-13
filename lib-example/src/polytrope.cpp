#include "ester.h"
#include "poly.h"
#include <cstdlib>

solution *solve_poly1d(double n, double tol, int nr) {
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
        printf("Error: %e\n", error);

        double relax = 1.;
        if(error>0.01) relax = 0.2; // Relax if error too large

        // Update variables
        Phi += relax*dPhi;
        Phi0 += relax*op.get_var("Phi0")(0);
        Lambda += relax*op.get_var("Lambda")(0);

        it++;
    }

    if(error>tol) {
        printf("No converge\n");
        return NULL;
    }

    figure fig("/XSERVE");
    fig.subplot(2, 1);

    fig.plot(map.r, Phi);
    fig.label("r", "Phi", "");

    S.set_value("Phi", Phi);
    S.set_value("Lambda", Lambda*ones(1, 1));
    S.set_value("Phi0", Phi0*ones(1, 1));
    fig.semilogy(map.r.block(1, -1, 0, 0), abs(eq.eval()).block(1, -1, 0, 0));
    fig.label("r", "Residual", "");

    printf("\nLambda = %f\n", Lambda);
    printf("Phi(0) = %f\n", Phi(0));
    printf("Phi(1) = %f\n", Phi(-1));
    printf("Boundary conditions:\n");
    printf("dPhi/dr(0) = %e\n", (map.D, Phi)(0));
    printf("dPhi/dr(1) + Phi(1) = %e\n", (map.D, Phi)(-1)+Phi(-1));

    return new solution(new matrix(Phi), Lambda);
}


solution *solve_poly2d(double n, double tol, int nr, int nt, int nex, double omega,
        matrix *init_phi, double lambda) {

	/**********  Mapping *************/

	/* Let's define a mapping for our problem using one internal domain and a one external
	 * domain. The initial internal domain will be spherical, extending from r=0 to r=1.
	 * Later, we will modify this mapping in each iteration in order to adapt it to the
	 * shape of the solution
	 */

	mapping map;
	map.set_ndomains(1); // # of sub-domains in the internal domain
	map.set_npts(nr); // # of points in the radial direction
	map.gl.set_xif(0., 1.); // Define values of zeta on the boundaries of the domains (only the internal ones)
	map.set_nex(nex); // # of radial points in the external domain
	map.set_nt(nt); // # of points in theta
	map.init(); // Initialize mapping

	/* The mapping is initialized to the default one, that is spherical boundaries with
	 * radii equal to the corresponding values of zeta (set with set_xif). If we want
	 * to change the shape of the domains we have to modify the variable map.R. Each row
	 * of map.R corresponds to a boundary, starting with the innermost boundary (in this
	 * case r=0). After that, we have to call map.remap()
	 */





	/********** Symbolic equations *************/

	/* We will now create a symbolic object to hold the equations that will be used in
	 * the calculation. The symbolic object is used as a helper for writing the equations
	 * in the solver, it can calculate automatically the jacobian of a given equation
	 * with respect to a variable, that simplifies a lot the task of adding this terms
	 * to the solver
	 *
	 * First, we create the object and add some variables to it
	 */

	symbolic S;
	sym symPhi = S.regvar("Phi");
	sym symPhiex = S.regvar("Phiex");
	sym symLambda = S.regvar("Lambda");
	sym symPhi0 = S.regvar("Phi0");

	// Now, we will define three equations.

	// The differential equation for the internal domain

	sym h = 1 - symLambda * (symPhi - symPhi0)
			+ 0.5 * omega * omega * S.r * S.r * sin(S.theta) * sin(S.theta);
	sym eq = lap(symPhi) - pow(sqrt(h * h), n);

	// The differential equation for the external domain

	sym eq_ex = lap(symPhiex);







	/************* Create solver **************/

	// Let's create a solver for 2 domains and 9 variables.

	solver op;
	op.init(2, 9, "full");

	// Now we have to add some variables to the solver.
	// We use two variables for Phi, one for the internal domain (Phi) and another
	// one for the external domain (Phiex)

	op.regvar("Phi");
	op.regvar("Phiex");

	// We add also scalar variables

	op.regvar("Lambda");
	op.regvar("Phi0");

	// The mapping will change during the calculation, so we should add also the variable
	// "r" to the solver

	op.regvar_dep("r");

	/* The variable "r" is not changed directly. In order to adapt the mapping we change
	 * only the boundaries between domains, and apply a fixed relation to obtain "r".
	 * In each domain, we have a mapping given by:
	 * 			r = f(eta,deta,R,dR,zeta,theta)
	 * zeta and theta are fixed for each point, so we don't take them into account, the
	 * other variables are defined in each domain as:
	 * 		- eta : Value of zeta of the internal boundary (scalar)
	 * 		- deta : Difference between the value of zeta of the external and internal
	 * 				boundaries (scalar)
	 * 		- R : Shape of the internal boundary (function of theta)
	 * 		- dR : Difference between the shape of the external and internal boundaries
	 * 				(function of theta)
	 * Note that r in each domain can depend only on variables defined in the same domain.
	 * That is the reason why we have added "r" as a "dependent variable", that means
	 * that we can use it to write equations in the solver, and the solver will substitute
	 * it by the real variables eta, deta, R, dR. Dependent variables doesn't increase
	 * the size of the jacobian matrix to be solved.
	 * We will now add the real variables.
	 */

	op.regvar("eta");
	op.regvar("deta");
	op.regvar("R");
	op.regvar("dR");

	// Instead of always writing the equations independently in each domain, sometimes
	// we want to be able to write them once and let the solver distribute them in the
	// different domains. For that, the solver has to now the size of each domain
	// (normally only the internal ones)

	op.set_nr(map.npts); // This doesn't take into account the external domain

	/* Each of the variables that we have defined "lives" in one or several domains and
	 * has a given size. We don't have to tell this directly to the solver. It will be
	 * inferred from the equations when we write them. But we have to have in mind that
	 * regular equations can include only variables defined in the same domain, while
	 * boundary conditions can mix variables from two consecutive domains.
	 * The structure of the variables will be as follows:
	 *
	 * 				Internal domain			External domain
	 * 				---------------			---------------
	 * 		Phi:		nr x nt					N/A
	 * 		Phiex:		N/A						nex x nt
	 * 		Lambda:		1 x 1					N/A
	 * 		Phi0:		1 x 1					N/A
	 * 		r:			nr x nt					nex x nt
	 * 		eta:		1 x 1					1 x 1
	 * 		deta: 		1 x 1					N/A
	 * 		R:			1 x nt					1 x nt
	 * 		dR:			1 x nt					N/A
	 *
	 */





	/************** Initial values ***************/

	// Create numerical variables for the solution with the initial guesses

	matrix Phi = map.r * map.r - 1;
	matrix Phiex = zeros(nex, nt);
	double Lambda = 1.;
	double Phi0 = 0.;

	/* Let's define two transformation matrices that we will use later:
	 Tpole: Given a variable f(theta) (1 x nt),  (f,Tpole) will give
	 	 	 the value of f in theta=0
	 Teq: Given a variable f(theta) (1 x nt),  (f,Teq) will give
	 	 	 the value of f in theta=pi/2
	 Tmean: Given a variable f(theta) (1 x nt),  (f,Tmean) will give
	 	 	 the mean value of f (l=0 component of
	 	 	 the expansion in Legendre polynomials)
	 */

	matrix Tmean = map.leg.P_00.col(0);
	matrix Tpole;
	map.leg.eval_00(ones(1, nt), 0, Tpole);
	matrix Teq;
	map.leg.eval_00(ones(1, nt), PI/2, Teq);


    // set initial guess if provided
    if (init_phi != NULL) {
        if (Phi.ncols() != init_phi->ncols() || Phi.nrows() != init_phi->nrows()) {
            printf("Initial solution doesn't have the correct size\n");
            return NULL;
        }
        Phi = *init_phi;
        matrix linspace = vector_t(1.0, 0.0, nex);
        Phiex = (vector_t(1.0, 0.0, nex), Phi.row(-1));
        Phi0 = Phi(0);
        Lambda = lambda;
    }


	// More initialization

	double error = 1.;
	int it = 0;

	/*********** Main loop (Newton's iteration) ***********/

	while (error > tol && it < 10000) {

		// Write the current values in the symbolic object (this is needed for the numeric
		// evaluation of the symbolic equations

		S.set_value("Phi", Phi);
		S.set_value("Lambda", Lambda * ones(1, 1)); // Note that the assigned value should be of type matrix
		S.set_value("Phi0", Phi0 * ones(1, 1));
		S.set_value("Phiex", Phiex);
		S.set_map(map);

		// Delete the equations of the previous iteration

		op.reset();

		/***** Equation for Phi *****/

		/* Here we use the symbolic equation eq that we have defined before.
		 * We have to write (only in the internal domain):
		 *
		 * 	D(eq)            D(eq)               D(eq)           D(eq)
		 *  ------D(Phi) + ---------D(Lambda) + -------D(Phi0) + -----D(r) = -eq
		 *  D(Phi)         D(Lambda)            D(Phi0)          D(r)
		 *
		 *  with the boundary conditions:
		 *
		 *  d(D(Phi))   -dPhi
		 *  --------- = -----      at the internal boundary (r=0)
		 *   dzeta      dzeta
		 *
		 *  and
		 *
		 *  D(Phi) - D(Phiex) = -(Phi -Phiex) at the external boundary (continuity)
		 *
		 *  Here, "d" means regular differentiation and "D" designates the corrections
		 *  to be calculated using Newton's method.
		 *
		 */

		eq.add(&op, "Phi", "Phi"); // Add the jacobian of eq with respect to "Phi"
		eq.add(&op, "Phi", "Lambda"); // Add the jacobian of eq with respect to "Lambda"
		eq.add(&op, "Phi", "Phi0"); // Add the jacobian of eq with respect to "Phi0"
		eq.add(&op, "Phi", "r"); // Add the jacobian of eq with respect to "r"

		// Internal boundary condition

		op.bc_bot2_add_l(0, "Phi", "Phi", ones(1, nt), map.D.block(0).row(0));

		// External boundary condition

		op.bc_top1_add_d(0, "Phi", "Phi", ones(1, nt));
		op.bc_top2_add_d(0, "Phi", "Phiex", -ones(1, nt));

		// Right-hand side

		matrix rhs = -eq.eval();

		// RHS for boundary conditions

		rhs.setrow(0, -(map.D, Phi).row(0));
		rhs.setrow(-1, -(Phi.row(-1) - Phiex.row(0)));

		// Write the RHS in the solver

		op.set_rhs("Phi", rhs);



		/*** Equation for Phiex ***/

		/* We use the symbolic equation eq_ex for writing (now in the external domain)
		 *
		 * 	D(eq_ex)           D(eq_ex)
		 *  --------D(Phiex) + --------D(r) = -eq_ex
		 *  D(Phiex)              D(r)
		 *
		 *  with boundary conditions:
		 *
		 *  - Continuity of the derivative (dPhi/dzeta)/rz at zeta=1. Note
		 *     that rz is NOT continuous at the boundary
		 *
		 *	D(Phiex) = -Phiex    at the external boundary (infinity)
		 *
		 */


		// We will write an equation for the external domain so we have to write the
		// values of r of the external domain in the symbolic object, otherwise the
		// internal values of r would be used in the numerical evaluation of symbolic
		// equations

		S.set_map(map.ex); // Change to external mapping

		// Now we add the equation

		eq_ex.add_ex(&op, 1, "Phiex", "Phiex");
		eq_ex.add_ex(&op, 1, "Phiex", "r");

		// and the boundary conditions

		op.bc_bot2_add_l(1, "Phiex", "Phiex", 1./map.ex.rz.row(0), map.ex.D.row(0));
		op.bc_bot1_add_l(1, "Phiex", "Phi", -1./map.rz.row(-1), map.D.row(-1));
		op.bc_bot2_add_l(1, "Phiex", "r", -1./map.ex.rz.row(0)/map.ex.rz.row(0)*
				(map.ex.D.row(0),Phiex),map.ex.D.row(0));
		op.bc_bot1_add_l(1, "Phiex", "r", 1./map.rz.row(-1)/map.rz.row(-1)*
						(map.D.row(-1),Phi),map.D.row(-1));

		op.bc_top1_add_d(1, "Phiex", "Phiex", ones(1, nt));

		// Finally the right-hand side

		rhs = -eq_ex.eval();

		rhs.setrow(0, -((map.ex.D.row(0), Phiex)/map.ex.rz.row(0)
				- (map.D.row(-1), Phi)/map.rz.row(-1)));
		rhs.setrow(-1,-Phiex.row(-1));

		op.set_rhs("Phiex", rhs);

		// We are done with this equation, so we can restore the internal mapping
		// in the symbolic object

		S.set_map(map);



		/**** Equation for Phi0 *****/

		/* Phi0 is just the value of Phi in r=0. We use it as a separate variable because
		 * we need to write an equation for Lambda that mixes the values of Phi(0) and
		 * Phi(pole) that correspond to different points (it is not a local equation).
		 * We will write the equation as a boundary condition at r=0 as we need the
		 * value of Phi in this point.
		 *
		 *   (D(Phi),Tmean) - D(Phi0) = -(Phi0 - (Phi,Tmean))    at r=0
		 *
		 *  Here we use the transformation Tmean to get the mean value of Phi(0,theta)
		 */

		// Variable to hold the RHS

		op.bc_bot2_add_r(0, "Phi0", "Phi", ones(1, 1), Tmean);
		op.bc_bot2_add_d(0, "Phi0", "Phi0", -ones(1, 1));

		rhs=-((Phi.row(0), Tmean)(0) - Phi0)*ones(1,1);

		op.set_rhs("Phi0", rhs);


		/**** Equation for Lambda ****/

		/* Lambda is defined as
		 *                    1
		 *   Lambda = ----------------
		 *            Phi(pole) - Phi(0)
		 *
		 * We will rewrite is as
		 *
		 *   Lambda * (Phi(pole) - Phi(0)) - 1 = 0
		 *
		 * And applying Newton's method we have the final equation that
		 * we will write as a boundary condition for the external boundary
		 *
		 *   Lambda*(D(Phi),Tpole) - Lambda*D(Phi0) + ((Phi,Tpole)-Phi0)*D(Lambda)
		 *                            = -(Lambda*((Phi,Tpole)-Phi0)-1)
		 *                                            at the external boundary
		 *
		 * Here we have used the transformation Tpole defined before.
		 */


		op.bc_top1_add_r(0, "Lambda", "Phi", Lambda * ones(1, 1), Tpole);
		op.bc_top1_add_d(0, "Lambda", "Phi0", -Lambda * ones(1, 1));

		op.bc_top1_add_d(0, "Lambda", "Lambda", ((Phi.row(-1),Tpole)(0) - Phi0) * ones(1, 1));

		rhs = -(Lambda * ((Phi.row(-1),Tpole)(0) - Phi0) - 1)*ones(1,1);

		op.set_rhs("Lambda", rhs);



		/******** Adaptive grid **********/

		/***** Equation for r *****/

		/* Now, we have to define the relation between the variable "r" and the real
		 * variables of the mapping eta, deta, R and dR.
		 * We use for that the jacobian of the mapping
		 *
		 * 	D(r) = map.J[0] * D(eta) + map.J[1] * D(deta) + map.J[2] * D(R) + map.J[3] * D(dR)
		 *
		 * Note that, because r is defined as a dependent variable, only terms of type
		 * "d" are allowed
		 */

		// For the internal domain

		op.add_d(0, "r", "eta", map.J[0]);
		op.add_d(0, "r", "deta", map.J[1]);
		op.add_d(0, "r", "R", map.J[2]);
		op.add_d(0, "r", "dR", map.J[3]);

		// and for the external domain

		op.add_d(1, "r", "eta", map.ex.J[0]);
		op.add_d(1, "r", "R", map.ex.J[2]);


		/***** Equation for eta *****/

		/* eta should be equal to the value of R at theta=0
		 */

		op.add_d(0, "eta", "eta", ones(1, 1)); // Internal domain
		op.add_d(1, "eta", "eta", ones(1, 1)); // External domain
		op.add_r(1,"eta","R",-ones(1,1),Tpole);

		op.set_rhs("eta", zeros(2, 1));


		/***** Equation for deta *****/

		/* By definition
		 *   D(deta) + D(eta[0]) - D(eta[1]) = 0
		 *
		 * where eta[0] is the value of eta in the internal domain and eta[1] its
		 * value in the external domain.
		 * As we are mixing variables in the internal and external domains, we write
		 * it as a boundary condition
		 * Note that deta is defined only in the internal domain
		 */

		op.bc_top1_add_d(0, "deta", "deta", ones(1, 1));
		op.bc_top1_add_d(0, "deta", "eta", ones(1, 1)); // top1 refers to the internal domain
		op.bc_top2_add_d(0, "deta", "eta", -ones(1, 1)); // top2 refers to the external domain

		op.set_rhs("deta", zeros(1, 1));

		/****** Equation for dR *******/

		/* Again, by definition
		 *   D(dR) + D(R[0]) - D(R[1]) = 0
		 */

		op.bc_top1_add_d(0, "dR", "dR", ones(1, nt));
		op.bc_top1_add_d(0, "dR", "R", ones(1, nt));
		op.bc_top2_add_d(0, "dR", "R", -ones(1, nt));

		op.set_rhs("dR", zeros(1, nt));


		/********* Equation for R *********/

		/* Here, we will have different equations for the internal and external domains.
		 * The internal domain is simple (R=0), but in the external domain we have to
		 * write our condition for the surface, that we will used to adapt the shape of
		 * the polytrope.
		 */

		// We initialize a variable to hold the RHS

		rhs = zeros(2, nt);

		/* Let's start with the internal domain. In the internal domain R is always 0,
		 * so we will just write
		 *
		 * 	D(R) = 0
		 */

		op.add_d(0, "R", "R", ones(1, nt));
		rhs.setrow(0, zeros(1, nt));

		/* Now comes the tricky part. We will write two different conditions in the
		 * external domain. The first one is written only in one point and is used
		 * to fix the equatorial radius (1). The second one is the condition for the surface
		 * that will be written in the remaining points of the boundary.
		 */

		// We define a variable to select the point where we will write the first
		// condition

		matrix qeq = zeros(1, nt);
		qeq(0) = 1;

		// Fix the equatorial radius

		op.bc_bot2_add_r(1, "R", "R", qeq, (Teq, ones(1, nt)));
		rhs.setrow(1,-qeq*0);

		/* Now we'll write the condition for the surface
		 *
		 *   1 - Lambda * (Phi(1,theta) - Phi0) + 0.5 * omega^2 * r^2 * sin(theta)^2 = 0
		 *
		 * We have to write it only on the points where 1-qpole = 1
		 *
		 */

		op.bc_bot1_add_d(1,"R","Lambda",-(1-qeq)*(Phi.row(-1)-Phi0));
		op.bc_bot1_add_d(1,"R","Phi0",(1-qeq)*Lambda);
		op.bc_bot2_add_d(1,"R","Phiex",-(1-qeq)*Lambda);
		op.bc_bot2_add_d(1,"R","R",(1-qeq)*(omega*omega*map.r.row(-1)*sin(map.th)*sin(map.th)));

		rhs.setrow(1,rhs.row(1)-(1-qeq)*
				(1-Lambda*(Phi.row(-1)-Phi0)+0.5*omega*omega*map.r.row(-1)*map.r.row(-1)*sin(map.th)*sin(map.th)));

		op.set_rhs("R",rhs);

		/********* Get solution ********/

		op.solve(); // Solve the equations

		matrix dPhi = op.get_var("Phi");
		error = max(abs(dPhi));  // Calculate the error (absolute)
		//printf("Error: %e\n", error);
		printf("Error: %d\t%e\n",it,error);

		double relax = 1.;
		if (error > 0.1)
			relax = 0.2; // Relax if error too large

		// Update variables
		Phi += relax * dPhi;
		Phiex += relax * op.get_var("Phiex");
		Phi0 += relax * op.get_var("Phi0")(0);
		Lambda += relax * op.get_var("Lambda")(0);

		// Update mapping
		map.R.setrow(1, map.R.row(1) + relax * op.get_var("R").row(1));
		map.remap();

		it++;

	}

	if (error > tol) {
		printf("No convergence\n");
		return NULL;
	}

	figure fig("/XSERVE");
	fig.subplot(2, 1);

	fig.colorbar();
	map.draw(&fig, Phi);
	fig.label("", "", "phi");

	S.set_value("Phi", Phi);
	S.set_value("Lambda", Lambda * ones(1, 1));
	S.set_value("Phi0", Phi0 * ones(1, 1));
	fig.semilogy(map.r.block(1, -1, 0, 0), abs(eq.eval()).block(1, -1, 0, 0));
	fig.label("r", "Residual", "");

	printf("Lambda = %f\n",Lambda);
	float Re = map.leg.eval_00(map.r.row(-1),PI/2)(0);
	float Rp = map.leg.eval_00(map.r.row(-1),0)(0);
	printf("eps = %f \n",(Re-Rp)/Re);

	return new solution(new matrix(Phi), Lambda);
}
