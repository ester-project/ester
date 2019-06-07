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
// set the pressure equation
    symbolic S;
    sym sym_mu=S.r;
    sym sym_lnP = S.regvar("lnP");
    sym sym_P = exp(sym_lnP);  // valid syntax ??
    sym sym_x = S.regvar("x");
    sym sym_R = S.regconst("R");
    sym sym_Pc = exp(S.regconst("lnPc"));  // valid syntax ??
    sym sym_kappa_s = S.regvar("opa.k");
    sym _fac1 = S.regconst("fac1")
    sym _fac2 = S.regconst("fac2")

    sym eqP = Dz(sym_lnP)+_fac1/pow(sym_R,4)/sym_Pc/sym_P*pow(sym_mu/sym_x,2);
    sym bcP = sym_P - _fac2/sym_Pc/sym_kappa_s/sym_R/sym_R;

    S.set_map(map);
    S.set_value("lnPc",log(Pc));
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
