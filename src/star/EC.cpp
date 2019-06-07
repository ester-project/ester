#include "constants.h"
#include "mapping.h"
#include "matrix.h"
#include "numdiff.h"
#include "parser.h"
#include "physics.h"
#include "solver.h"
#include "symbolic.h"

#include <cstdlib>

int main(int argc,char *argv[]) {

    symbolic::spherical = true; // a definir initialement

    int nr=100;
    int ndomains=1;
    mapping map;
    map.set_ndomains(ndomains);
    map.set_npts(nr);
    map.gl.set_xif(0., 1.);
    map.set_nt(1); // 1d
    map.init();



// the initial guesses
    matrix lnT = log(1.-0.5*map.r*map.r);
    matrix lnP = lnT;
    matrix x = map.r;
    matrix lam = map.r;
    double M = 5*M_SUN;

// Initialize the solver
    solver op;
    int nvar=16;
    op.init(ndomains,nvar,"full")
    op.set_nr(map.npts);

    op.regvar("lnP");
    op.regvar("lnT");
    op.regvar("x");
    op.regvar("lam");
    op.regvar("lnPc");
    op.regvar("lnTc");
    op.regvar("lnR");
    op.regvar("lnLum");
    op.regvar_dep("P");
    op.regvar_dep("T");
    op.regvar_dep("rho");
    op.regvar_dep("lnrho");
    op.regvar_dep("lnrhoc");
    op.regvar_dep("lnxi");
    op.regvar_dep("lnxic");
    op.regvar_dep("eps");
    op.regvar_dep("nabla");


    int it = 0;
    double error = 1.;
    while(error>tol && it<10000) {

        op.reset(); // Delete the equations of the previous iteration

	solve_definitions(op);
	solve_lnP(op);
	solve_lnT(op);
	solve_x(op);
	solve_lam(op);

        op.solve(); // Solve the equations

// get the variations
        matrix dlnP = op.get_var("lnP");
        matrix dlnT = op.get_var("lnT");
        matrix dx = op.get_var("x");
        matrix dlam = op.get_var("lam");
// Variations on the global scalar constant
// ont-elles besoin d'etre des matrix ????
        matrix dlnPc = op.get_var("lnPc");
        matrix dlnTc = op.get_var("lnTc");
        matrix dR = op.get_var("R");
        matrix dL = op.get_var("Lum");

// get the amplitude of the variation on the lnT field
	error=max(abs(dlnT));

        double relax = 1.;
        if(error>0.01) relax = 0.2; // Relax if error too large

// Update main variables
        lnT += relax*dlnT;
        lnP += relax*dlnP;
        x += relax*dx;
        lam += relax*dlam;
        R += relax*op.get_var("R"); // matrix pas matrix ??
        lnTc += relax*op.get_var("lnTc");
        lnPc += relax*op.get_var("lnPc");
        Lum += relax*op.get_var("Lum");
	printf("it = %d  lnTc = %e\n",it,lnTc);

// update all the dependent fields: xi,eps,rho, etc...
	fields();

        it++;
    }

    if(error>tol) {
        ester_err("No converge\n");
    }

}
//--------------------------------------------------------------------
void solve_lnp(solver *op) {
// First the non-dimensional factor:
// Mass and Radius in solar units, Central pressure in 1e17 cgs
    double P_ref=1e17;
    double kappa_ref=0.1;
    double fac1=3*GRAV*M_SUN*M_SUN/8/PI/pow(R_SUN,4)/P_ref;
// fac1=0.0134698
    double fac2=2*GRAV*M_SUN/3/P_ref/kappa_ref/R_SUN/R_SUN;
// fac2=1.82952e-12

// set the pressure equation
    symbolic S;
    sym sym_mu=S.r;
    sym sym_lnP = S.regvar("lnP");
    sym sym_P = exp(sym_lnP);  // valid syntax ??
    sym sym_x = S.regvar("x");
    sym sym_lnR = log(S.regconst("R")); // valid syntax ??
    sym sym_Pc = exp(S.regconst("lnPc"));  // valid syntax ??
    sym log_kappa = log(S.regvar("opa.k"));
    sym _fac1 = S.regconst("fac1");
    sym log_fac2 = log(S.regconst("fac2"));
    sym log_M = log(S.regconst("M"));

    sym eqP = Dz(sym_lnP)+_fac1/pow(sym_R,4)/sym_Pc/sym_P*pow(sym_mu/sym_x,2);
    sym bcP = sym_lnP - log_fac2-log_M+sym_lnPc+log_kappa+2*sym_lnR;

    S.set_map(map);
    S.set_value("lnPc",log(Pc));
    S.set_value("R",R);
    S.set_value("lnP",lnP);
    S.set_value("mu", mu);
    S.set_value("fac1", fac1);
    S.set_value("fac2", fac2);
    S.set_value("M", M);
// We now make the insertion into the jacobian
    eqP.add(op, "lnP", "lnP");
    eqP.add(op, "lnP", "mu"); 
    eqP.add(op, "lnP", "x"); 
    eqP.add(op, "lnP", "R");
    eqP.add(op, "lnP", "lnPc");
    rhs=-eqP.eval();

// Central BC
    op->bc_bot2_add_d(0,"lnP","lnP",ones(1,1));
    rhs(0)=-lnP(0);

// Surface BC
    bcP.bc_top1_add(op, "lnP", "lnPc");
    bcP.bc_top1_add(op, "lnP", "R");
    bcP.bc_top1_add(op, "lnP", "opa.k");
    //rhs(-1)=-(P(-1)-2*GRAV*M/3/Pc/kappa(-1)/R/R);
    rhs(-1)=-bcP.eval()(-1); // est-ce possible ?
    op->set_rhs("lnP",rhs);
}
//-------------------------------------------------------------------
void solve_lnT(solver *op) {
// Set the temperature equation
// Luminosity, Mass and Radius in solar units
    double P_ref=1e17;
    double kappa_ref=0.1;
    double T_ref=1e6;
    double rho_ref=1e1;
    double xi_ref=1e14;
    double log_fac=log(L_SUN/4/PI/GRAV/M_SUN*P_ref/T_ref/rho_ref/xi_ref);
// fac=230.167
    double log_fac_bc=log(pow(L_SUN/4/PI/SIG_SB/R_SUN/R_SUN,0.25)/T_ref);
// fac_bc=5.77716e-3

    matrix rhs;
    symbolic S;
    sym sym_mu=S.r;
    sym sym_lam = S.regvar("lam");
    sym log_P = S.regvar("lnP");
    sym log_T = S.regvar("lnT");
    sym log_rho = S.regvar("lnrho");
    sym log_xi = S.regvar("lnxi");
    sym sym_nabla = S.regvar("nabla");
    sym log_nabla = log(S.regvar("nabla"));

// The global scalars
    sym log_fac = S.regconst("log_fac");
    sym log_facbc = S.regconst("log_facbc");
    sym log_Pc = S.regconst("lnPc");
    sym log_rhoc = S.regconst("lnrhoc");
    sym log_Tc = S.regconst("lnTc");
    sym log_L = log(S.regconst("Lum"));
    sym log_R = log(S.regconst("R"));
    sym log_M = log(S.regconst("M"));

    sym eqT = Dz(sym_lnT) - sym_nabla*Dz(sym_lnP);

    sym eq_nabla = log_nabla - log_fac-log_Lum-log_M-log_Pc+log_rhoc
     +log_Tc-log_P+log_xi+log_rho+log_T-1.5*log(sym_lam/sym_mu);

    sym bcT = sym_lnT - log_facbc-0.25*log_L+0.5*log_R+log_Tc;

    S.set_map(map);
    S.set_value("lnP",lnP);
    S.set_value("lnT",lnT);
    S.set_value("mu", mu);
    S.set_value("lam",lam);
    S.set_value("rho",rho);
    S.set_value("nabla",nabla);
    S.set_value("lnrhoc",log(rhoc));
    S.set_value("lnTc",log(Tc));
    S.set_value("lnPc",log(Pc));
    S.set_value("lnxi",log(opa.xi));
    S.set_value("M",M);
    S.set_value("log_fac",log_fac);
    S.set_value("log_fac_bc",log_fac_bc);
// We now make the insertion into the jacobian
    eqT.add(op, "lnT", "lnT");
    eqT.add(op, "lnT", "lnP"); 
    eqT.add(op, "lnT", "nabla");
    rhs=-eqT.eval();
    op->set_rhs("lnT",rhs);

    eq_nabla.add(op, "nabla", "nabla");
    eq_nabla.add(op, "nabla", "lnP");
    eq_nabla.add(op, "nabla", "lnPc");
    eq_nabla.add(op, "nabla", "lnT");
    eq_nabla.add(op, "nabla", "lnTc");
    eq_nabla.add(op, "nabla", "lnrho");
    eq_nabla.add(op, "nabla", "lnrhoc");
    eq_nabla.add(op, "nabla", "lnxi");
    eq_nabla.add(op, "nabla", "lam");
    eq_nabla.add(op, "nabla", "mu");
    eq_nabla.add(op, "nabla", "lum");
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
    double rho_ref=1e1;
    double fac=3*M_SUN/8/PI/pow(R_SUN,3)/rho_ref;
// fac=0.07057
    matrix rhs;
    symbolic S;
    sym sym_mu=S.r;
    sym sym_x = S.regvar("x");
    sym sym_rho = S.regvar("rho");
    sym sym_rhoc = S.regconst("rhoc");
    sym sym_R = S.regconst("R");

    sym eqx = Dz(sym_x) - fac*M/pow(sym_R,3)/sym_rhoc/sym_rho
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
    double eps_ref=1e3;
    double fac=M_SUN/L_SUN*eps_ref;
// fac=518.049

    symbolic S;
    sym sym_mu=S.r;
    sym sym_eps = S.regvar("eps");
    sym sym_lam = S.regvar("lam");
    sym sym_Lum = S.regconst("Lum");
    sym eqL = Dz(sym_lam) - fac*M/sym_Lum*sym_eps*sqrt(sym_mu/sym_lam);

    S.set_map(map);
    S.set_value("R",R);
    S.set_value("lnrho",log(rho));
    S.set_value("lnrhoc",log(rhoc));
    S.set_value("mu", mu);
    S.set_value("lam", lam);
// We now make the insertion into the jacobian
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
        op->add_d("rho","lnrho",rho);

        op->add_d("lnrho","lnP",1./eos.chi_rho);
        op->add_d("lnrho","lnT",-eos.d);
        op->add_d("lnrho","lnPc",1./eos.chi_rho);
        op->add_d("lnrho","lnTc",-eos.d);
        op->add_d("lnrho","ln_rhoc",-ones(nr,1));

// Expression of \delta\khi (thermal conductivity)      
        op->add_d("lnxi","lnrho",opa.dlnxi_lnrho);
        op->add_d("lnxi","ln_rhoc",opa.dlnxi_lnrho);
        op->add_d("lnxi","lnT",opa.dlnxi_lnT);
        op->add_d("lnxi","lnTc",opa.dlnxi_lnT);
        op->add_d("lnxi","ln_xic",-ones(nr,1));

// Expression of \delta\kappa (opacity)         
        op->add_d("opa.k","lnT",3*opa.k);
        op->add_d("opa.k","lnTc",3*opa.k);
        op->add_d("opa.k","lnrho",-opa.k);
        op->add_d("opa.k","ln_rhoc",-opa.k);
        op->add_d("opa.k","lnxi",-opa.k);
        op->add_d("opa.k","ln_xic",-opa.k);

        op->add_d("nuc.eps","lnrho",nuc.dlneps_lnrho*nuc.eps);
        op->add_d("nuc.eps","ln_rhoc",nuc.dlneps_lnrho*nuc.eps);
        op->add_d("nuc.eps","lnT",nuc.dlneps_lnT*nuc.eps);
        op->add_d("nuc.eps","lnTc",nuc.dlneps_lnT*nuc.eps);
        op->add_d("nuc.eps","ln_epsc",-nuc.eps);
}

void fields(){

        Y0=1.-X0-Z0;
        init_comp();
        eq_state();
	opacity();
	nuclear();
	atmosphere();

}

void opacity() {
        int error;

        error=opa_calc(comp.X(),Z0,Tc*T,rhoc*rho,opa);
        xic=opa.xi(0);
        opa.xi=opa.xi/xic;

        if(error) exit(1);
}

void nuclear() {
        int error;

        error=nuc_calc(comp,T*Tc,rho*rhoc,nuc);
        epsc=nuc.eps(0);
        nuc.eps=nuc.eps/epsc;

        if(error) exit(1);
}

void eq_state() {
   int error;

   matrix rhoc_m(1,1);
   eos_calc(comp.X()(0,0)*ones(1,1),Z0,ones(1,1)*Tc,ones(1,1)*pc,rhoc_m,eos);
   rhoc=rhoc_m(0);

   error=eos_calc(comp.X(),Z0,T*Tc,p*pc,rho,eos);

   rho=rho/rhoc;

   if(error) exit(1);
}

void atmosphere() {
        int error;
        if(!strcmp(atm.name,"simple")) {
                ps=surff*gsup()/opa.k.row(-1)/pc;
                matrix q;
                q=surff*p.row(-1)/ps;
                Ts=pow(q,0.25)*Teff()/Tc;
                return;
        }

        error=atm_calc(comp.X(),Z0,gsup(),Teff(),eos.name,opa.name,atm);

        if(error) exit(1);

        matrix q;
        ps=surff*atm.ps/pc;
        q=surff*p.row(-1)/ps;//q=ones(1,nth);
        Ts=pow(q,0.25)*atm.Ts/Tc;
}
