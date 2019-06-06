#include "ester.h"
// #include "matplotlib.h"

#include <cstdlib>

matrix solve_ester_cesam(double M, double tol, int nr) {
    // Create a mapping

    symbolic::spherical = true; // a definir initialement

    printf("parametres: tol %e, nr %d\n",tol,nr);
    mapping map;
    map.set_ndomains(1);
    map.set_npts(nr);
    map.gl.set_xif(0., 1.); // This changes the values of zeta of the limits of the domain directly (do it before map.init())
    map.set_nt(1); // 1d
    map.init();


// The structure differential equations

// Additional equation for nabla

// Create numerical variables for the solution with the initial guesses
    matrix T = 1.-0.5*map.r*map.r;
    matrix P = T;

    double error = 1.;
    double M = 5*M_SUN;
    int it = 0;

    // Create a solver for the three variables. The size of each variable is determined by the
    // solver object based on the size of the equations that we will define for them.
    // The "Phi" equation will have nr x 1 points while for "Phi0" and "Lambda" we will introduce
    // only boundary conditions, so the resulting equations will have 1 x 1 size
    star1d A;
    solver *op;
    op=A.init_solver()
    op->regvar("lnP");
    op->regvar("lnT");
    op->regvar("x");
    op->regvar("lam");
    op->regvar("rhoc");
    op->regvar("Tc");
    op->regvar("R");
    op->regvar("Lum");
    op->regvar("xi");
    op->regvar("eps");
    op->regvar("rho");
    op->regvar("nabla");

    op.set_nr(map.npts);

    while(error>tol && it<10000) {
        op.reset(); // Delete the equations of the previous iteration







        op.solve(); // Solve the equations

        matrix dPhi = op.get_var("Phi");
        error = max(abs(dPhi));  // Calculate the error (absolute)

        double relax = 1.;
        if(error>0.01) relax = 0.2; // Relax if error too large

        // Update variables
        Phi += relax*dPhi;
        Phi0 += relax*op.get_var("Phi0")(0);
        Lambda += relax*op.get_var("Lambda")(0);
	//printf("it = %d, Lambda =  %e\n",it,Lambda);

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
	//printf("err %e Lambda = %e \n",fabs(hi - hsurf),Lambda);
        ri += (hsurf - hi)/dhi;
        hi = map.gl.eval(h, ri)(0);
        dhi = map.gl.eval(dh, ri)(0);
    }
    if (nit >= maxit) {
        ester_err("No convergence finding polytrope surface such that h(surface) = %e", hsurf);
    }

    matrix h2 = map.gl.eval(h, map.r*ri);

    return h2;
}
//--------------------------------------------------------------------
void solve_lnp(solver *op) {
// set the pressure equation
    symbolic S;
    sym sym_mu=S.r;
    sym sym_lnP = S.regvar("lnP");
    sym sym_P = exp(sym_lnP);
    sym sym_x = S.regvar("x");
    sym sym_R = S.regconst("R");
    sym sym_Pc = exp(S.regconst("lnPc"));
    sym sym_kappa_s = S.regvar("kappa_s");
    sym _fac1 = S.regconst("fac1")
    sym _fac2 = S.regconst("fac2")

    sym eqP = Dz(sym_lnP)+_fac1/pow(sym_R,4)/sym_Pc/sym_P*pow(sym_mu/sym_x,2);
    sym bcP = sym_P - _fac2/sym_Pc/sym_kappa_s/sym_R/sym_R;

    S.set_map(map);
    S.set_value("lnPc",lnPc);
    S.set_value("R",R);
    S.set_value("lnP",lnP);
    S.set_value("mu", mu);
    S.set_value("fac1", 3*GRAV*M*M/8/PI);
    S.set_value("fac2", 2*GRAV*M/3);
// We now make the insertion into the jacobian
    eqP.add(op, "lnP", "lnP");
    eqP.add(op, "lnP", "mu"); 
    eqP.add(op, "lnP", "x"); 
    eqP.add(op, "lnP", "R");
    eqP.add(op, "lnP", "Pc");
    rhs=-eqP.eval();

// Central BC
    op->bc_bot2_add_d(0,"lnP","lnP",ones(1,1));
    rhs(0)=-lnP(0);

// Surface BC
    bcP.bc_top1_add(op, "lnP", "Pc");
    bcP.bc_top1_add(op, "lnP", "R");
    bcP.bc_top1_add(op, "lnP", "kappa_s");
    //rhs(-1)=-(P(-1)-2*GRAV*M/3/Pc/kappa(-1)/R/R);
    rhs(-1)=-bcP.eval()(-1); // est-ce possible ?
    op->set_rhs("lnP",rhs);
}
//-------------------------------------------------------------------
void solve_lnT(solver *op) {
// Set the temperature equation

    matrix rhs;
    symbolic S;
    sym sym_mu=S.r;
    sym sym_lnP = S.regvar("lnP");
    sym sym_P = exp(sym_lnP);
    sym sym_lnT = S.regvar("lnT");
    sym sym_T = exp(sym_lnT); //    relation (1)
    sym sym_nabla = S.regvar("nabla");
    sym sym_lam = S.regvar("lam");

// Some dependent variables
    sym sym_rho = S.regvar("rho");
    sym sym_xi = S.regvar("xi");
// The global scalars
    sym sym_Pc = S.regconst("Pc");
    sym sym_rhoc = S.regconst("rhoc");
    sym sym_Tc = S.regconst("Tc");
    sym sym_Lum = S.regconst("Lum");

    sym eqT = Dz(sym_lnT) - sym_nabla*Dz(sym_lnP);
    sym eq_nabla = sym_nabla - sym_Lum/4/PI/GRAV/M*sym_P*sym_Pc
         /sym_rho/sym_rhoc/sym_Tc/sym_T/sym_xi*pow(sym_lam/sym_mu,1.5);
    sym bcT = sym_T - pow(sym_lum/4/PI/SIG_SB,0.25)/sym_Tc/sqrt(sym_R);

    S.set_map(map);
    S.set_value("lnP",lnP);
    S.set_value("mu", mu);
    S.set_value("Pc",Pc);
    S.set_value("rho",rho);
    S.set_value("rhoc",rhoc);
    S.set_value("lnT",lnT);
// est-ce qu'il faut un set_value pour T où est-ce que (1) suffit ?
// meme question pour P.
    S.set_value("Tc",Tc);
    S.set_value("xi",opa.xi);
    S.set_value("lam",lam);
    S.set_value("nabla",nabla);
// We now make the insertion into the jacobian
    eqT.add(op, "lnT", "lnT");
    eqT.add(op, "lnT", "lnP"); 
    eqT.add(op, "lnT", "nabla");
    rhs=-eqT.eval();
    op->set_rhs("lnT",rhs);

    eq_nabla.add(op, "nabla", "nabla");
    eq_nabla.add(op, "nabla", "lnP");
    eq_nabla.add(op, "nabla", "Pc");
    eq_nabla.add(op, "nabla", "lnT");
    eq_nabla.add(op, "nabla", "Tc");
    eq_nabla.add(op, "nabla", "rho");
    eq_nabla.add(op, "nabla", "xi");
    eq_nabla.add(op, "nabla", "lam");
    eq_nabla.add(op, "nabla", "mu");
    eq_nabla.add(op, "nabla", "lum");
    eq_nabla.add(op, "nabla", "rhoc");
    rhs=zeros(nr,1);
    op->set_rhs("nabla",rhs);

// Central BC
    op->bc_bot2_add_d(0,"lnT","lnT",ones(1,1));
    rhs(0)=-lnT(0);

// Surface BC
    bcT.bc_top1_add(op, "lnT", "Lum");
    bcT.bc_top1_add(op, "lnT", "Tc");
    bcT.bc_top1_add(op, "lnT", "R");
    rhs(-1)=-bcT.eval()(-1); // est-ce possible ?
    op->set_rhs("lnT",rhs);
}
//-------------------------------------------------------------------
void solve_x(solver *op) {
// Set the x-variable equation
    matrix rhs;
    symbolic S;
    sym sym_mu=S.r;
    sym sym_x = S.regvar("x");
    sym sym_rho = S.regvar("rho");
    sym sym_rhoc = S.regconst("rhoc");
    sym sym_R = S.regconst("R");

    sym eqx = Dz(sym_x) - 3*M/8/PI/pow(sym_R,3)/sym_rhoc/sym_rho
                            *sqrt(sym_mu/sym_x);

    S.set_map(map);
    S.set_value("R",R);
    S.set_value("rho",rho);
    S.set_value("rhoc",rhoc);
    S.set_value("mu", mu);
    S.set_value("x", x);
// We now make the insertion into the jacobian
    eqx.add(op, "x", "x"); 
    eqx.add(op, "x", "mu");
    eqx.add(op, "x", "rho");
    eqx.add(op, "x", "rhoc");
    eqx.add(op, "x", "R");
    rhs=-eqx.eval();
    op->set_rhs("x",rhs);

// Central BC
    op->bc_bot2_add_d(0,"x","x",ones(1,1));
    rhs(0)=-x(0);

// Surface BC
    op->bc_top1_add(0, "x", "x",ones(1,1));
    rhs(-1)=1d0-x(-1);

    op->set_rhs("x",rhs);
}
//-------------------------------------------------------------------
void solve_lam(solver *op) {
// Set the luminosity-variable equation

    symbolic S;
    sym sym_mu=S.r;
    sym sym_eps = S.regvar("eps");
    sym sym_lam = S.regvar("lam");
    sym sym_Lum = S.regconst("Lum");
    sym eqL = Dz(sym_lam) - M/sym_Lum*sym_eps*sqrt(sym_mu/sym_lam);
// We now make the insertion into the jacobian
    eqL.add(op, "lam", "lam");
    eqL.add(op, "lam", "mu");
    eqL.add(op, "lam", "eps");
    eqL.add(op, "lam", "Lum");
    rhs=-eqL.eval();
    op->set_rhs("lam",rhs);

// Central BC
    op->bc_bot2_add_d(0,"lam","lam",ones(1,1));
    rhs(0)=-lam(0);

// Surface BC
    op->bc_top1_add(0, "lam", "lam",ones(1,1));
    rhs(-1)=1d0-lam(-1);
}

void solve_definitions(solver *op) {

        op->add_d("T","lnT",T);
        op->add_d("P","lnP",P);

        op->add_d("rho","lnP",rho/eos.chi_rho);
        op->add_d("rho","lnT",-rho*eos.d);
        op->add_d("rho","lnPc",rho/eos.chi_rho);
        op->add_d("rho","lnTc",-rho*eos.d);
        op->add_d("rho","ln_rhoc",-rho);

// Expression of \delta\khi (thermal conductivity)      
        op->add_d("xi","rho",opa.dlnxi_lnrho*opa.xi/rho);
        op->add_d("xi","ln_rhoc",opa.dlnxi_lnrho*opa.xi);
        op->add_d("xi","lnT",opa.dlnxi_lnT*opa.xi);
        op->add_d("xi","lnTc",opa.dlnxi_lnT*opa.xi);
        op->add_d("xi","ln_xic",-opa.xi);

// Expression of \delta\kappa (opacity)         
        op->add_d("opa.k","lnT",3*opa.k);
        op->add_d("opa.k","lnTc",3*opa.k);
        op->add_d("opa.k","rho",-opa.k/rho);
        op->add_d("opa.k","ln_rhoc",-opa.k);
        op->add_d("opa.k","xi",-opa.k/opa.xi);
        op->add_d("opa.k","ln_xic",-opa.k);

        op->add_d("nuc.eps","rho",nuc.dlneps_lnrho*nuc.eps/rho);
        op->add_d("nuc.eps","ln_rhoc",nuc.dlneps_lnrho*nuc.eps);
        op->add_d("nuc.eps","lnT",nuc.dlneps_lnT*nuc.eps);
        op->add_d("nuc.eps","lnTc",nuc.dlneps_lnT*nuc.eps);
        op->add_d("nuc.eps","ln_epsc",-nuc.eps);
}
