#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include"star.h"
#include<stdlib.h>
#include<sys/time.h>
#include<string.h>
#include"symbolic.h"

#include <omp.h>
#define TNEW

/// \brief Initialize star's chemical composition, equation of state, opacity,
/// nuclear reaction, mass definition, pi_c, Lambda, velocity, units, atmosphere,
/// flatness, scaled keplerian angular velocity
void star2d::fill() {
	Y0=1.-X0-Z0;
	//init_comp();

	eq_state();
	opacity();
	nuclear();

	m=2*PI*(map.gl.I,rho*r*r*map.rz,map.leg.I_00)(0);

	R=pow(M/m/rhoc,1./3.);

	pi_c=(4*PI*GRAV*rhoc*rhoc*R*R)/pc;
	Lambda=rhoc*R*R/Tc;

	calc_units();

	atmosphere();

	//Omegac=sqrt(map.leg.eval_00(((Dex,phiex)/rex/map.ex.rz).row(0),PI/2)(0));
	double eps=1.-1./map.leg.eval_00(r.row(-1),PI/2)(0);
	Omegac=sqrt(pi_c*m/4/PI*(1-eps)*(1-eps)*(1-eps));

	calc_vangle();

}

void star2d::init_comp() {
// Update the object comp

	comp=initial_composition(X0,Z0)*ones(nr,nth);

	if(!conv) return;

// Count the number of point in the core:
    int n = 0;
    for (int i=0; i<conv; i++) {
        n += map.gl.npts[i];
    }

    if(stratified_comp == 0) {
		printf("Calling initial_composition with Xc = %e.", Xc*X0);
        comp.setblock(0,n-1,0,-1,initial_composition(Xc*X0,Z0)*ones(n,nth));
    }
	else if (stratified_comp == 2)
	{
		printf("Taking X profile from model.");	
	}
    else {
        comp.setblock(0, n-1, 0, -1,
                initial_composition(Xc*X0, Z0)*ones(n, nth));
        int m = 0;
        int l = n;
        double a = (1.-exp((1.-(nr/n))))/(X0*(1.-Xc));
        for(int i=conv+1; i<ndomains; i++) {
            m = map.gl.npts[i];
            comp.setblock(l, l+m-1, 0, -1,
                    initial_composition(((Xc*X0)+((1./a)*(1.-exp((1.-(l/n)))))),Z0)*ones(m,nth));
            l += m;
        }
    }
}

solver *star2d::init_solver(int nvar_add) {
	int nvar;
	solver *op;

	nvar=37;
	op=new solver;
	op->init(ndomains+1,nvar+nvar_add,"full");

	op->maxit_ref=10;op->use_cgs=1;op->maxit_cgs=20;op->debug=0;
	op->rel_tol=1e-12;op->abs_tol=1e-20;
	register_variables(op);

	return op;
}

void star2d::register_variables(solver *op) {
	int i,var_nr[ndomains+1];

	for(i=0;i<ndomains;i++)
		var_nr[i]=map.gl.npts[i];
	var_nr[ndomains]=nex;
	op->set_nr(var_nr);

	op->regvar("Phi");
	op->regvar_dep("p");
	op->regvar("log_p");
	op->regvar("pi_c");
	op->regvar_dep("T");
	op->regvar("log_T");
	op->regvar("Lambda");
	op->regvar("eta");
	op->regvar("deta");
	op->regvar("Ri");
	op->regvar("dRi");
	op->regvar_dep("r");
	op->regvar_dep("rz");
	op->regvar("Omega");
	op->regvar_dep("log_rhoc");
	op->regvar("log_pc");
	op->regvar("log_Tc");
	op->regvar("log_R");
	op->regvar("m");
	op->regvar("ps");
	op->regvar("Ts");
	op->regvar("Teff");
	op->regvar("gsup");
	op->regvar("w");
	op->regvar("vr");
	op->regvar("vt");
	op->regvar_dep("rho");
	op->regvar_dep("opa.xi");
	op->regvar_dep("opa.k");
	op->regvar_dep("nuc.eps");
	op->regvar_dep("s");
	op->regvar("gamma");
	op->regvar("Omega0");
	op->regvar_dep("log_Re");
	op->regvar("sin_vangle");
	op->regvar_dep("cos_vangle");


}

double star2d::solve(solver *op) {
    matrix_map error_map;
    return solve(op, error_map, 0);
}

void star2d::write_eqs(solver *op) {

#ifdef MKL
	int mkl_threads = mkl_get_max_threads();
	mkl_set_num_threads(1);   // Faster with 1 thread !?
#endif

#ifdef THREADS
	int num_threads = 12;
	std::thread t[num_threads];
	t[0] = std::thread(&star2d::solve_definitions, this, op);
	t[1] = std::thread(&star2d::solve_poisson, this, op);
	t[2] = std::thread(&star2d::solve_mov, this, op);
	t[3] = std::thread(&star2d::solve_temp, this, op);
	t[4] = std::thread(&star2d::solve_dim, this, op);
	t[5] = std::thread(&star2d::solve_map, this, op);
	t[6] = std::thread(&star2d::solve_Omega, this, op);
	t[7] = std::thread(&star2d::solve_atm, this, op);
	t[8] = std::thread(&star2d::solve_gsup, this, op);
	t[9] = std::thread(&star2d::solve_Teff, this, op);
	t[10] = std::thread(&star2d::solve_cont, this, op);
	t[11] = std::thread(&star2d::solve_vangle, this, op);
	for(int i=0; i< num_threads; i++)
		t[i].join();
#else
	solve_definitions(op);
	solve_poisson(op);
	solve_mov(op);
	solve_temp(op);
	solve_dim(op);
	solve_map(op);
	solve_Omega(op);
	solve_atm(op);
	solve_gsup(op);
	solve_Teff(op);
	solve_cont(op);
	solve_vangle(op);
#endif
#ifdef MKL
	mkl_set_num_threads(mkl_threads);
#endif
}

double star2d::update_solution(solver *op, double &h, matrix_map& error_map, int nit) {

	double err,err2,dmax;

// h : relaxation parameter for Newton solver: useful for the first
// iterations
	dmax=config.newton_dmax;

	matrix dphi,dphiex,dp,dT,dpc,dTc;
	dphi=op->get_var("Phi").block(0,nr-1,0,-1);
	dphiex=op->get_var("Phi").block(nr,nr+nex-1,0,-1);
	err=max(abs(dphi/phi));
    error_map["Phi"](nit) = err;

	//printf("err(phi)=%e\n",err);
	dp=op->get_var("p");
	err2=max(abs(dp/p));err=err2>err?err2:err;
    error_map["p"](nit) = err2;
	while(exist(abs(h*dp/p)>dmax)) h/=2;

	//printf("err(p)=%e\n",err2);
	dT=op->get_var("T");
	err2=max(abs(dT/T));err=err2>err?err2:err;
    error_map["T"](nit) = err2;
	while(exist(abs(h*dT/T)>dmax)) h/=2;

	//printf("err(T)=%e\n",err2);
	dpc=op->get_var("log_pc");
	err2=fabs(dpc(0));err=err2>err?err2:err;
    error_map["log_pc"](nit) = err2;
	while(fabs(h*dpc(0))>dmax) h/=2;

	//printf("err(pc)=%e\n",err2);
	dTc=op->get_var("log_Tc");
	err2=fabs(dTc(0));err=err2>err?err2:err;
    error_map["log_Tc"](nit) = err2;
	while(fabs(h*dTc(0))>dmax) h/=2;
	//printf("err(Tc)=%e\n",err2);

	phi+=h*dphi;
	phiex+=h*dphiex;
	p+=h*dp;
	T+=h*dT;
	pc*=exp(h*dpc(0));
	Tc*=exp(h*dTc(0));
	Omega=Omega+h*op->get_var("Omega")(0);
	w+=h*op->get_var("w");
	vr+=h*op->get_var("vr");
	vt+=h*op->get_var("vt");

	matrix dRi;
	dRi=op->get_var("Ri");
    error_map["Ri"](nit) = max(abs(dRi));
	update_map(h*dRi);
	err2=max(abs(dRi));err=err2>err?err2:err;

	return err;
}

/// \brief Performs one step of the Newton algorithm to compute the star's
/// internal structure.
double star2d::solve(solver *op, matrix_map& error_map, int nit) {
	int info[5];

	if(Omega==0 && Omega_bk != 0) {
		Omega=Omega_bk*Omegac;
		w=Omega*ones(nr,nth);
	}

	check_map();

// Clear up the preceding equations (clear up is not necessary
// if the system is linear) reset does not clear up the LU factorization
// which can be reused:
	op->reset();

	if(config.verbose) {printf("Writing equations...");fflush(stdout);}
	write_eqs(op);
	if(config.verbose) printf("Done\n");

// Solving the system:
	op->solve(info);

// Some output verbose ----------------------
	if (config.verbose) {
		if(info[2]) {
			printf("CGS Iteration: ");
			if(info[4]>0)
				printf("Converged after %d iterations\n",info[4]);
			else
				printf("Not converged (Error %d)\n",info[4]);
		}
		if(info[0]) 
			printf("Solved using LU factorization\n");
		if(info[1]) {
			printf("CGS Refinement: ");
			if(info[3]>0)
				printf("Converged after %d iterations\n",info[3]);
			else
				printf("Not converged (Error %d)\n",info[3]);
		}
	}
// End  output verbose ----------------------

	double h = 1.;
	double err = update_solution(op, h, error_map, nit);

	matrix rho0=rho;

	fill();

	double err2=max(abs(rho-rho0));err=err2>err?err2:err;

	return err;

}

// Special treatment for updating R_i because we need
// a different relaxation parameter "h".

void star2d::update_map(matrix dR) {
	double h=1,dmax=config.newton_dmax;

	matrix R0;
	R0=map.R;
	dR.setrow(0,zeros(1,nth));
	while(exist(abs(h*dR)>dmax*R0)) h/=2;
	map.R+=h*dR;
	while(map.remap()) {
		h/=2;
		map.R=R0+h*dR;
	}

}

/// \brief insert the definitions depending on opacity and eos tables into the solver,
/// and the definitions used by the mapping (eta,deta,Ri,dRi,...), and the entropy
void star2d::solve_definitions(solver *op) {

// EOS written rho(P,T). eos.chi_rho is from OPAL,
// and eos.d=chi_t/chi_rho
	op->add_d("rho","p",rho/eos.chi_rho/p);
	op->add_d("rho","T",-rho*eos.d/T);
	op->add_d("rho","log_pc",rho/eos.chi_rho);
	op->add_d("rho","log_Tc",-rho*eos.d);
	op->add_d("rho","log_rhoc",-rho);

// Constitutive relation for thermal radiative conductivity with xi=xi(rho,T)
	op->add_d("opa.xi","rho",opa.dlnxi_lnrho*opa.xi/rho);
	op->add_d("opa.xi","log_rhoc",opa.dlnxi_lnrho*opa.xi);
	op->add_d("opa.xi","T",opa.dlnxi_lnT*opa.xi/T);
	op->add_d("opa.xi","log_Tc",opa.dlnxi_lnT*opa.xi);

// Constitutive relation for opacity (k) k=16sigma*T^3/rho/xi
	op->add_d("opa.k","T",3*opa.k/T);
	op->add_d("opa.k","log_Tc",3*opa.k);
	op->add_d("opa.k","rho",-opa.k/rho);
	op->add_d("opa.k","log_rhoc",-opa.k);
	op->add_d("opa.k","opa.xi",-opa.k/opa.xi);

// Constitutive relation for nuclear heat generation; formal dependence nuc.eps=nuc.eps(rho,T)
	op->add_d("nuc.eps","rho",nuc.dlneps_lnrho*nuc.eps/rho);
	op->add_d("nuc.eps","log_rhoc",nuc.dlneps_lnrho*nuc.eps);
	op->add_d("nuc.eps","T",nuc.dlneps_lnT*nuc.eps/T);
	op->add_d("nuc.eps","log_Tc",nuc.dlneps_lnT*nuc.eps);

// Here we write the variation of r as a function of the variations of the mapping
// parameters according to equation 6.10 of the manual
	op->add_d("r","eta",map.J[0]);
	op->add_d("r","deta",map.J[1]);
	op->add_d("r","Ri",map.J[2]);
	op->add_d("r","dRi",map.J[3]);

// same as above but for the last vacuum domain
	op->add_d(ndomains,"r","eta",map.ex.J[0]);
	op->add_d(ndomains,"r","Ri",map.ex.J[2]);

// same as for r but for rz
	op->add_d("rz","eta",(D,map.J[0]));
	op->add_d("rz","deta",(D,map.J[1]));
	op->add_d("rz","Ri",(D,map.J[2]));
	op->add_d("rz","dRi",(D,map.J[3]));

// last domain for rz
	op->add_d(ndomains,"rz","eta",(Dex,map.ex.J[0]));
	op->add_d(ndomains,"rz","Ri",(Dex,map.ex.J[2]));

// From the variation of entropy with respect to the variation of pressure and
// temperature but valid only for homogeneus composition !!!!!!
	op->add_d("s","T",eos.cp/T);
	op->add_d("s","log_Tc",eos.cp);
	op->add_d("s","p",-eos.cp*eos.del_ad/p);
	op->add_d("s","log_pc",-eos.cp*eos.del_ad);

	/*eos_struct eos2;
	matrix rho2;
	strcpy(eos2.name, eos.name);
	double dlogT = 1e-4;
	matrix T2 = T*(1+dlogT);
	eos_calc(comp.X(), Z0, T2*Tc, p*pc, rho2, eos2);
	matrix dcpdlnT = (eos2.cp - eos.cp)/dlogT;
	matrix ddeladdlnT = (eos2.del_ad - eos.del_ad)/dlogT;
	double dlogp = 1e-4;
	matrix p2 = p*(1+dlogp);
	eos_calc(comp.X(), Z0, T*Tc, p2*pc, rho2, eos2);
	matrix dcpdlnp = (eos2.cp - eos.cp)/dlogp;
	matrix ddeladdlnp = (eos2.del_ad - eos.del_ad)/dlogp;
	op->add_d("eos.cp", "log_p", dcpdlnp);
	op->add_d("eos.cp", "log_pc", dcpdlnp);
	op->add_d("eos.cp", "log_T", dcpdlnT);
	op->add_d("eos.cp", "log_Tc", dcpdlnT);
	op->add_d("eos.del_ad", "log_p", ddeladdlnp);
	op->add_d("eos.del_ad", "log_pc", ddeladdlnp);
	op->add_d("eos.del_ad", "log_T", ddeladdlnT);
	op->add_d("eos.del_ad", "log_Tc", ddeladdlnT);*/

}

void star2d::calc_vangle() {
	static symbolic S;
	static bool sym_inited = false;
	static sym sin_vangle;
	if (!sym_inited) {
		sym p = S.regvar("p");
		sym_vec thvec = S.r*grad(S.theta);
		sin_vangle = (grad(p), thvec)/sqrt((grad(p), grad(p)));
		sym_inited = true;
	}
	S.set_value("p", p);
	S.set_map(map);
	vangle = asin(sin_vangle.eval());
	vangle.setrow(0, zeros(1, nth));
}

///\brief vangle is the angle between e_r and the
/// normal of an isobar. Positive in the northern hemisphere
void star2d::solve_vangle(solver *op) {
	static symbolic S;
	static bool sym_inited = false;
	static sym eq;
	if (!sym_inited) {
		sym p = S.regvar("p");
		sym sin_vangle = S.regvar("sin_vangle");
		sym_vec thvec = S.r*grad(S.theta);
		eq = sin_vangle - (grad(p), thvec)/sqrt((grad(p), grad(p)));
		sym_inited = true;
	}
	S.set_value("p", p);
	S.set_value("sin_vangle", sin(vangle), 11);
	S.set_map(map);

	eq.add(op, "sin_vangle", "sin_vangle");
	eq.add(op, "sin_vangle", "p");
	eq.add(op, "sin_vangle", "r");
	op->bc_bot2_add_d(0, "sin_vangle", "sin_vangle", ones(1, nth));
	op->set_rhs("sin_vangle", zeros(nr, nth));

	op->add_d("cos_vangle", "sin_vangle", -sin(vangle)/cos(vangle));
}

/// \brief Writes Poisson equation and interface conditions into the solver.
void star2d::solve_poisson(solver *op) {
	matrix q,rhs1,rhs2,rhs;
	int n,j0;
	matrix &rz=map.rz;

	symbolic S;
	sym lap_phi;
	{
        sym phi;
        phi=S.regvar("Phi");
        lap_phi=lap(phi);
	}

	S.set_value("Phi",phi);
	S.set_map(map);

	lap_phi.add(op,"Phi","Phi");
// The equation is named Phi and depends on variable tagged "Phi"
	lap_phi.add(op,"Phi","r");
// The equation is named Phi and depends on variable tagged "r"

	rhs1=-lap_phi.eval()+pi_c*rho; // Expression of the RHS inside the star

	//rho
	op->add_d("Phi","rho",-pi_c*ones(nr,nth));  // Coefficient of delta rho

	//pi_c
	op->add_d("Phi","pi_c",-rho);  // Coefficient of delta pi_c

// phiex:
// Outside the star we define phiex but this is still equation "Phi"
	S.set_value("Phi",phiex);
	S.set_map(map.ex);
// add exterior domain with variable tagged "Phi":
	lap_phi.add_ex(op,ndomains,"Phi","Phi");
// add exterior domain with variable tagged "r":
	lap_phi.add_ex(op,ndomains,"Phi","r");

	rhs2=-lap_phi.eval();     // Expression of RHS outside the star (rho=0)

	rhs=zeros(nr+nex,nth);          // set the RHS vector
	rhs.setblock(0,nr-1,0,-1,rhs1);
	rhs.setblock(nr,nr+nex-1,0,-1,rhs2);

	j0=0;

// Loop on the domains from n=0 to n=ndomains
// n=0 is the first domain
// n=ndomains-1 is the last domain inside the star
// n=ndomain is the domain outside the star (the remaining Universe !)
// Each domain has a top and a bottom. On the top we impose the continuity
// of phi, on the bottom we impose the continuity of (1/rz)(dphi/dzeta)
// 1 designate the "Left condition" associated with the (n-1)th domain,
// and 2 designate the "Right condition" associated with the nth domain.

// set the interface conditions in the Jacobian matrix
// i.e. continuity of Phi and (1/rz)*dphi/dzeta

	for(n=0;n<ndomains+1;n++) {
	  if(n==0) {
// First start with the bottom conditions
// In the first domain we demand that dphi/dzeta
// vanishes at the center (right BC):
		op->bc_bot2_add_l(n,"Phi","Phi",ones(1,nth),D.block(0).row(0));
		rhs.setrow(0,-(D,phi).row(0));
	   } else {
	if(n<ndomains)
// for all stellar domains, at the bottom (1/rz)d(delta phi)/dzeta
// is continuuous (Right block tagged 2))
                op->bc_bot2_add_l(n,"Phi","Phi",1/rz.row(j0),D.block(n).row(0));
// in the last stellar domain continuity demanded with the "ex" fields
// (Right block)
	else op->bc_bot2_add_l(n,"Phi","Phi",1/map.ex.rz.row(0),Dex.row(0));

// for all stellar domains, at the bottom (1/rz)d(delta phi)/dzeta
// is continuuous (left block assoc. with domain n-1)
	op->bc_bot1_add_l(n,"Phi","Phi",-1/rz.row(j0-1),D.block(n-1).row(-1));

// Functional derivative of "1/rz times dphi/dzeta",s
// complements the interface cond. (Right cond. block n)
	if(n<ndomains) op->bc_bot2_add_d(n,"Phi","rz",-1/rz.row(j0)/rz.row(j0)*(D,phi).row(j0));
	else op->bc_bot2_add_d(n,"Phi","rz",-1/map.ex.rz.row(0)/map.ex.rz.row(0)*(Dex,phiex).row(0));

// For all stellar domain Left cond. in block n-1 of the complement of IC
  op->bc_bot1_add_d(n,"Phi","rz",1/rz.row(j0-1)/rz.row(j0-1)*(D,phi).row(j0-1));

// set the interface conditions in the RHS
	if(n<ndomains)
      rhs.setrow(j0,-(D,phi).row(j0)/rz.row(j0)+(D,phi).row(j0-1)/rz.row(j0-1));
	else
      rhs.setrow(j0,-(Dex,phiex).row(0)/map.ex.rz.row(0)+(D,phi).row(j0-1)/rz.row(j0-1));
		} // end of block n#0

// Continuity of Phi imposed on top of the domains
// left condition for Jacobian (in domain n-1):
	op->bc_top1_add_d(n,"Phi","Phi",ones(1,nth));

// Right condition for Jacobian (in domain n)
	if(n<ndomains) op->bc_top2_add_d(n,"Phi","Phi",-ones(1,nth));

// Prepare the RHS for the continuity of phi
   if(n<ndomains) rhs.setrow(j0+map.gl.npts[n]-1,-phi.row(j0+map.gl.npts[n]-1));
   else rhs.setrow(nr+nex-1,-phiex.row(nex-1)); // left terms of the last domain
   if(n<ndomains-1) rhs.setrow(j0+map.gl.npts[n]-1,rhs.row(j0+map.gl.npts[n]-1)
	+phi.row(j0+map.gl.npts[n])); // add the right term
// add the right term for the last stellar domain but here the
// right term is phiex:
   else if(n==ndomains-1) rhs.setrow(j0+map.gl.npts[n]-1,
    rhs.row(j0+map.gl.npts[n]-1)+phiex.row(0));
   if(n<ndomains) j0+=map.gl.npts[n];
// Note the Right BC for phi is phi(infty)=0,
// so no right term in the rhs of the n==ndomain case
	}
	op->set_rhs("Phi",rhs);
}

///\brief Solve the continuity equation
void star2d::solve_cont(solver *op) {

	static symbolic S;
	static sym eq, flux;
	static bool sym_inited = false;
	if (!sym_inited) {
		sym rho = S.regvar("rho");
		sym vr = S.regvar("vr");
		sym vt = S.regvar("vt");

		sym_vec phivec(COVARIANT);
		phivec(0) = 0*S.one; phivec(1) = 0*S.one; phivec(2) = S.r*sin(S.theta);
		sym_vec rvec = grad(S.r);
		sym_vec thvec = cross(phivec, rvec);
		sym_vec nvec = grad(S.zeta);
		nvec = nvec / sqrt((nvec, nvec));

		sym_vec V = vr*rvec + vt*thvec;
		eq = div(rho*V);
		flux = rho * (nvec, V);
		sym_inited = true;
	}

	S.set_value("vr", vr);
	S.set_value("vt", vt, 11);
	S.set_value("rho", rho);
	S.set_map(map);

	eq.add(op, "vr", "vr");
	eq.add(op, "vr", "vt");
	eq.add(op, "vr", "rho");
	eq.add(op, "vr", "r");

	matrix rhs = -eq.eval();
	matrix flux_val = flux.eval();

	int j0 = 0;
	for (int n = 0; n < ndomains; n++) {
		int j1 = j0 + map.npts[n] - 1;

		if (n == 0) {
			op->bc_bot2_add_d(n, "vr", "vr", ones(1, nth));
			rhs.setrow(0, -vr.row(0));
		}

		if (n == conv-1) {
			flux.bc_top1_add(op, n, "vr", "vr");
			flux.bc_top1_add(op, n, "vr", "vt");
			flux.bc_top1_add(op, n, "vr", "rho");
			flux.bc_top1_add(op, n, "vr", "r");
			flux.bc_top2_add(op, n, "vr", "vr", -ones(1, nth));
			flux.bc_top2_add(op, n, "vr", "vt", -ones(1, nth));
			flux.bc_top2_add(op, n, "vr", "rho", -ones(1, nth));
			flux.bc_top2_add(op, n, "vr", "r", -ones(1, nth));
			rhs.setrow(j1, -flux_val.row(j1) + flux_val.row(j1+1));
		}
		else if (n == ndomains - 1) {
			flux.bc_top1_add(op, n, "vr", "vr");
			flux.bc_top1_add(op, n, "vr", "vt");
			flux.bc_top1_add(op, n, "vr", "rho");
			flux.bc_top1_add(op, n, "vr", "r");
			rhs.setrow(j1, -flux_val.row(j1));
		}
		else {
			op->bc_top1_add_d(n, "vr", "vr", ones(1, nth));
			op->bc_top2_add_d(n, "vr", "vr", -ones(1, nth));
			rhs.setrow(j1, -vr.row(j1) + vr.row(j1+1));

		}
		j0 += map.npts[n];
	}

	op->set_rhs("vr", rhs);


}

void star2d::solve_mov(solver *op) {
	static bool eqinit = false;
	static symbolic S;
	static sym eq_phi, bc_w, ic_w, ic_t, ic_dt, ic_visc, bc_t, eqmov_r, eqmov_t;
	if(!eqinit) {
		eqinit = true;

		sym p = S.regvar("p");  // Pressure
		sym vt = S.regvar("vt");
		sym vr = S.regvar("vr");
		sym w = S.regvar("w");  // Differential Rotation
		sym rho = S.regvar("rho"); // Density
		sym phi = S.regvar("Phi"); // Gravitational Potential
		sym reynolds_v = S.regconst("reynolds_v"); // Re_vertical
		sym reynolds_h = S.regconst("reynolds_h"); // Re_Horiz.
		//sym Re = exp(S.regconst("log_Re"));
		sym Re = exp(S.regvar("log_Re"));
		sym visc_ratio = S.regconst("visc_ratio"); // nu_v/nu_h
		sym sin_vangle = S.regvar("sin_vangle");   // sin(v)
		sym cos_vangle = S.regvar("cos_vangle"); // v= angle(e_r,-gradP)

		sym_vec phivec(COVARIANT);
		phivec(0) = 0*S.one; phivec(1) = 0*S.one; phivec(2) = S.r*sin(S.theta);
		sym_vec rvec = grad(S.r);
		sym_vec thvec = cross(phivec, rvec); // rvec=e_r
		sym_vec nvec = grad(S.zeta);
		nvec = nvec / sqrt((nvec, nvec)); // nvec = E^zeta normalized
		sym_vec tvec = cross(phivec, nvec);

		sym_vec V, svec;
		sym s;
		s = S.r*sin(S.theta);
		svec = grad(s);
		V = vr*rvec + vt*thvec;
// define v_vec perp isobar h_vec// isobar
		sym_vec v_vec = cos_vangle*rvec - sin_vangle*thvec;
		sym_vec h_vec = sin_vangle*rvec + cos_vangle*thvec;

		// The momentum equation  is eqmov=0; it gives the pressure
		sym_vec eqmov = grad(p)/rho + grad(phi) - s*w*w*svec;

// nu_mc = tensorial kin. viscosity
// Re=R*R/visc_h/MYR = diffusion time in Myrs
// reynolds_v = true Reynolds number
		sym_tens nu_mc = tensor(v_vec, v_vec)*Re/reynolds_v + tensor(h_vec, h_vec)*Re/reynolds_h + tensor(phivec, phivec)*Re/reynolds_h;
		sym_tens SS = (grad(V), nu_mc.T()) + (nu_mc, grad(V).T()) - rational(2, 3)*S.g*(nu_mc%grad(V));;
		//sym_tens SS = stress(V);
		sym_vec visc = (div(rho*SS))/rho; // viscous force/rho

		sym_vec eqmov_visc = covariant(eqmov - visc);
		eqmov_r = eqmov_visc(0);
		eqmov_t = eqmov_visc(1);

// The angular momentum equation

		sym_tens nu = tensor(v_vec, v_vec)*visc_ratio/Re + tensor(h_vec, h_vec)/Re;
		eq_phi = rho*(V, grad(s*s*w)) - div(s*s*rho*(nu, grad(w)));

		ic_t = (tvec, V);  // tangential comp of V (iso zeta)
		ic_visc = ((SS, nvec), tvec); // tangential comp of stress on isozeta
		ic_dt = S.Dz(vt)/S.rz; // dv_theta/dr 

		ic_w = S.Dz(w)/S.rz; // dOmega/dr

		bc_t = ic_visc;  // associated boundary condition
		bc_w = (nvec, grad(w)); // stress-free for Omega
	}

	int j0;

	S.set_value("p",p);
	S.set_value("vt",vt,11); // 11= symmetry  of v_theta
	S.set_value("vr", vr);
	S.set_value("w",w);
	S.set_value("rho",rho);
	S.set_value("Phi",phi);
	S.set_value("reynolds_v", reynolds_v*ones(1,1));
	S.set_value("reynolds_h", reynolds_h*ones(1,1));
    matrix log_Re = log(R*R/visc_h/MYR)*ones(nr, nth);
	double log_Re_conv = log(R*R/1e8/MYR); // Set the viscosity in the convective core to 1e11 cm^2/s. We assume the vertical and horizontal viscosity are equal.
    // Viscosity in the core which is enhanced!!
    //if (conv) {
    //	int nc = 0;
    //	for (int n = 0; n < conv; n++) nc += map.npts[n];
    //	log_Re.setblock(0, nc-1, 0, -1, ones(nc, nth) * log_Re_conv);
    //}
	S.set_value("log_Re", log_Re);

	//S.set_value("log_Re", log(R*R/visc_h/MYR)*ones(1,1));
	S.set_value("visc_ratio", visc_v/visc_h*ones(1,1));
	S.set_value("sin_vangle", sin(vangle), 11);// 11= symmetry  of v_theta
	S.set_value("cos_vangle", cos(vangle));

	S.set_map(map);

	op->add_d("p","log_p",p);
	op->add_d("log_Re", "log_R", 2*ones(nr, nth)); // comfort equation

// First component of the equation of motion tagged "log_p"
// we explicit the dependences
	eqmov_r.add(op,"log_p","p");
	eqmov_r.add(op,"log_p","w");
	eqmov_r.add(op,"log_p","rho");
	eqmov_r.add(op,"log_p","Phi");
	eqmov_r.add(op,"log_p","r");
	eqmov_r.add(op,"log_p","vr");
	eqmov_r.add(op,"log_p","vt");
	eqmov_r.add(op,"log_p","log_Re");
	eqmov_r.add(op,"log_p","sin_vangle");
	eqmov_r.add(op,"log_p","cos_vangle");
	op->set_rhs("log_p",-eqmov_r.eval());


	matrix rhs;

	// log_p - Pressure
	rhs=op->get_rhs("log_p");

	j0 = 0;
	for(int n=0; n<ndomains; n++) {
		//int j1 = j0 + map.npts[n] - 1;
		if (n == 0) {
			op->bc_bot2_add_d(0,"log_p","p",ones(1,nth));
			rhs.setrow(0,-p.row(0)+1); // central pressure is unity
		}
		else {
			op->bc_bot2_add_d(n,"log_p","p",ones(1,nth));
			op->bc_bot1_add_d(n,"log_p","p",-ones(1,nth));
			rhs.setrow(j0,-p.row(j0)+p.row(j0-1));
		}

		if (n == ndomains-1) {
		}
		else {
		}

		j0 += map.npts[n];
	}
	op->set_rhs("log_p",rhs);

	//pi_c -  non-dimensional parameter
	rhs=zeros(ndomains,1);
	j0=0;
	for(int n=0;n<ndomains;n++) {
		if(n<ndomains-1) {
			op->bc_top1_add_d(n,"pi_c","pi_c",ones(1,1));
			op->bc_top2_add_d(n,"pi_c","pi_c",-ones(1,1));
		} else {
			matrix TT;
			map.leg.eval_00(th,0,TT);
			op->bc_top1_add_r(n,"pi_c","ps",-ones(1,1),TT);
			op->bc_top1_add_r(n,"pi_c","p",ones(1,1),TT);
			rhs(ndomains-1)=(ps-p.row(-1),TT)(0);
		}

		j0+=map.gl.npts[n];
	}
	op->set_rhs("pi_c",rhs);

	if(Omega==0) {
		op->add_d("w","w",ones(nr,nth));
		op->set_rhs("w",zeros(nr,nth));

		op->add_d("vt","vt",ones(nr,nth));
		op->set_rhs("vt",zeros(nr,nth));

		for(int n=0; n<ndomains; n++) op->add_d(n, "Omega0", "Omega0", ones(1,1));
		op->set_rhs("Omega0",zeros(ndomains,1));
		return;
	}


	// Second component of the equation of motion, dependences...
	eqmov_t.add(op,"vt","p");
	eqmov_t.add(op,"vt","w");
	eqmov_t.add(op,"vt","rho");
	eqmov_t.add(op,"vt","Phi");
	eqmov_t.add(op,"vt","r");
	eqmov_t.add(op,"vt","vr");
	eqmov_t.add(op,"vt","vt");
	eqmov_t.add(op,"vt","log_Re");
	eqmov_t.add(op,"vt","sin_vangle");
	eqmov_t.add(op,"vt","cos_vangle");
	op->set_rhs("vt",-eqmov_t.eval());

	// Equation of angular momentum, dependences...
	eq_phi.add(op,"w","w");
	eq_phi.add(op,"w","vr");
	eq_phi.add(op,"w","vt");
	eq_phi.add(op,"w","rho");
	eq_phi.add(op,"w","r");
	eq_phi.add(op,"w","sin_vangle");
	eq_phi.add(op,"w","cos_vangle");
	eq_phi.add(op,"w","log_Re");
	op->set_rhs("w",-eq_phi.eval());

	// w - Differential Rotation

	rhs = op->get_rhs("w");
	matrix ic_w_val = ic_w.eval();
	j0 = 0;
	for(int n=0; n<ndomains; n++) {
		int j1 = j0 + map.npts[n] - 1;
		if (n == 0) {
			op->bc_bot2_add_d(n, "w", "w", ones(1, nth));
			op->bc_bot2_add_d(n, "w", "Omega0", -ones(1, nth));
			rhs.setrow(j0, zeros(1, nth));
// Omega0 is an auxiliary variable to impose at equator both Omega=Omega(input) and
// stress-free condition
// So there is an extra-scalar equation for Omega0
// cf Lambda equation at center.

// w(theta=pi/2,r=1)=w_eq imposed from input
		}
		else {
// continuity of Omega at interface
			op->bc_bot2_add_d(n, "w", "w", ones(1, nth));
			op->bc_bot1_add_d(n, "w", "w", -ones(1, nth));
			rhs.setrow(j0, -w.row(j0) + w.row(j0-1));
		}
		if (n == ndomains - 1) {
			bc_w.bc_top1_add(op, n, "w", "w");
			bc_w.bc_top1_add(op, n, "w", "r");
			rhs.setrow(j1, -bc_w.eval().row(-1));
		}
		else {
// continuity of dOmega/dr at interface
			ic_w.bc_top1_add(op, n, "w", "w");
			ic_w.bc_top1_add(op, n, "w", "r");
			ic_w.bc_top2_add(op, n, "w", "w", -ones(1, nth));
			ic_w.bc_top2_add(op, n, "w", "r", -ones(1, nth));
			rhs.setrow(j1, -ic_w_val.row(j1) + ic_w_val.row(j1+1));
		}
		j0 += map.npts[n];
	}
	op->set_rhs("w",rhs);

	rhs = zeros(ndomains, 1);
	for(int n=0; n<ndomains; n++) {
		if (n == ndomains-1) {
			matrix TT;
			map.leg.eval_00(th, PI/2, TT);
			op->bc_top1_add_r(n, "Omega0", "w", ones(1,1), TT);
			op->bc_top1_add_d(n, "Omega0", "Omega", -ones(1,1));
			rhs(n) = -(w.row(-1), TT)(0) + Omega;
		}
		else {
			op->bc_top1_add_d(n, "Omega0", "Omega0", ones(1,1));
			op->bc_top2_add_d(n, "Omega0", "Omega0", -ones(1,1));
		}
	}
	op->set_rhs("Omega0", rhs);

	//Vtheta - Meridional circulation
	rhs=op->get_rhs("vt");
	matrix ic_t_val = ic_t.eval();
	matrix ic_dt_val = ic_dt.eval();
	matrix ic_visc_val = ic_visc.eval();

	j0 = 0;
	for(int n=0; n<ndomains; n++) {
		int j1 = j0 + map.npts[n] - 1;
		if (n == 0) {
			op->bc_bot2_add_d(n ,"vt", "vt", ones(1, nth));
			rhs.setrow(j0, -vt.row(j0));
		}
		else if (n == conv){
		// impose continuity of tangential velocity
			ic_t.bc_bot2_add(op, n, "vt", "vt");
			ic_t.bc_bot2_add(op, n, "vt", "vr");
			ic_t.bc_bot2_add(op, n, "vt", "r");
			ic_t.bc_bot1_add(op, n, "vt", "vt", -ones(1, nth));
			ic_t.bc_bot1_add(op, n, "vt", "vr", -ones(1, nth));
			ic_t.bc_bot1_add(op, n, "vt", "r", -ones(1, nth));
			rhs.setrow(j0, -ic_t_val.row(j0) + ic_t_val.row(j0-1));
		}
		else {
			op->bc_bot2_add_d(n ,"vt", "vt", ones(1, nth));
			op->bc_bot1_add_d(n ,"vt", "vt", -ones(1, nth));
			rhs.setrow(j0, -vt.row(j0) + vt.row(j0-1));
		}

		if (n == ndomains-1) {
			bc_t.bc_top1_add(op, n, "vt", "p"); // ? may be not necessary
			bc_t.bc_top1_add(op, n, "vt", "rho");
			bc_t.bc_top1_add(op, n, "vt", "Phi"); // ?
			bc_t.bc_top1_add(op, n, "vt", "w");
			bc_t.bc_top1_add(op, n, "vt", "r");
			bc_t.bc_top1_add(op, n, "vt", "vr");
			bc_t.bc_top1_add(op, n, "vt", "vt");
			rhs.setrow(j1, -bc_t.eval().row(j1));
		}
		else if (n == conv-1) {
			ic_visc.bc_top1_add(op, n, "vt", "vr");
			ic_visc.bc_top1_add(op, n, "vt", "vt");
			ic_visc.bc_top1_add(op, n, "vt", "r");
			ic_visc.bc_top2_add(op, n, "vt", "vr", -ones(1, nth));
			ic_visc.bc_top2_add(op, n, "vt", "vt", -ones(1, nth));
			ic_visc.bc_top2_add(op, n, "vt", "r", -ones(1, nth));
			rhs.setrow(j1, -ic_visc_val.row(j1) + ic_visc_val.row(j1+1));
		}
		else {
			ic_dt.bc_top1_add(op, n, "vt", "vt");
			ic_dt.bc_top1_add(op, n, "vt", "r");
			ic_dt.bc_top2_add(op, n, "vt", "vt", -ones(1, nth));
			ic_dt.bc_top2_add(op, n, "vt", "r", -ones(1, nth));
			rhs.setrow(j1, -ic_dt_val.row(j1) + ic_dt_val.row(j1+1));
		}

		j0 += map.npts[n];
	}
	op->set_rhs("vt",rhs);

}


void star2d::solve_temp(solver *op) {

	double xi_conv = 1.e13; // PARAMETER dimensional 

	static symbolic S;
	static sym eq, icflux;
	static bool sym_inited = false;
	if (!sym_inited){
		sym T = S.regvar("T");
		sym p = S.regvar("p");
		sym xi = S.regvar("opa.xi");
		sym eps = S.regvar("nuc.eps");
		sym Lambda = S.regvar("Lambda");
		sym rho = S.regvar("rho");
		sym vr = S.regvar("vr");
		sym vt = S.regvar("vt");
		sym lnTc = S.regconst("log_Tc");
		sym xi_conv = S.regconst("xi_conv");
		sym cp = S.regvar("eos.cp");
		sym del_ad = S.regvar("eos.del_ad");
		sym MYR = S.regconst("MYR");

		sym_vec phivec(COVARIANT);
		phivec(0) = 0*S.one; phivec(1) = 0*S.one; phivec(2) = S.r*sin(S.theta);
		sym_vec rvec = grad(S.r);
		sym_vec thvec = cross(phivec, rvec);
		sym_vec nvec = grad(S.zeta);
		nvec = nvec / sqrt((nvec, nvec));

		sym_vec grad_s = cp * (grad(T)/T - del_ad * grad(p)/p);
		sym_vec flux = -xi * grad(T) - xi_conv * grad_s;
		icflux = (nvec, flux)/Lambda;

		sym_vec rhoV = rho*(vr*rvec + vt*thvec);

		sym advec = exp(lnTc) / MYR * T * (rhoV, grad_s);

		eq = div(flux)/Lambda - rho*eps + advec;
		sym_inited = true;

	}

	matrix xi_conv_mat = zeros(nr, nth);
	int j0 = 0;
	for (int n=0; n<conv; n++) {
		int j1 = j0 + map.npts[n] - 1;
		xi_conv_mat.setblock(j0, j1, 0, -1, xi_conv*ones(map.npts[n], nth));
		j0 += map.npts[n];
	}

	S.set_value("T", T);
	S.set_value("p", p);
	S.set_value("opa.xi", opa.xi);
	S.set_value("nuc.eps", nuc.eps);
	S.set_value("Lambda", Lambda*ones(1,1));
	S.set_value("rho", rho);
	S.set_value("xi_conv", xi_conv_mat);
	S.set_value("eos.cp", eos.cp);
	S.set_value("eos.del_ad", eos.del_ad);
	S.set_value("vt", vt, 11);
	S.set_value("vr", vr);
	S.set_value("log_Tc", log(Tc)*ones(1,1));
	S.set_value("MYR", MYR*ones(1,1));
	S.set_map(map);

	op->add_d("T","log_T",T);

	eq.add(op, "log_T", "T");
	eq.add(op, "log_T", "p");
	eq.add(op, "log_T", "opa.xi");
	eq.add(op, "log_T", "nuc.eps");
	eq.add(op, "log_T", "Lambda");
	eq.add(op, "log_T", "rho");
	eq.add(op, "log_T", "vr");
	eq.add(op, "log_T", "vt");
	eq.add(op, "log_T", "log_Tc");
	eq.add(op, "log_T", "r");


	//eq.add(op, "log_T", "eos.cp");
	//eq.add(op, "log_T", "eos.del_ad");

	matrix rhs = -eq.eval();
	matrix icflux_val = icflux.eval();

	matrix dT = (D, T);
	j0 = 0;
	for (int n = 0; n < ndomains; n++) {
		int j1 = j0 + map.npts[n] - 1;

		if (n == conv && conv > 0) {
			icflux.bc_bot2_add(op, n, "log_T", "T");
			icflux.bc_bot2_add(op, n, "log_T", "p");
			icflux.bc_bot2_add(op, n, "log_T", "opa.xi");
			icflux.bc_bot2_add(op, n, "log_T", "r");
			icflux.bc_bot2_add(op, n, "log_T", "Lambda");
			icflux.bc_bot1_add(op, n, "log_T", "T", -ones(1, nth));
			icflux.bc_bot1_add(op, n, "log_T", "p", -ones(1, nth));
			icflux.bc_bot1_add(op, n, "log_T", "opa.xi", -ones(1, nth));
			icflux.bc_bot1_add(op, n, "log_T", "r", -ones(1, nth));
			icflux.bc_bot1_add(op, n, "log_T", "Lambda", -ones(1, nth));
			rhs.setrow(j0, -icflux_val.row(j0) + icflux_val.row(j0-1));
		}
		else {
			op->bc_bot2_add_l(n, "log_T", "T", 1./map.rz.row(j0), D.block(n).row(0));
			op->bc_bot2_add_d(n, "log_T", "rz", -1./map.rz.row(j0)/map.rz.row(j0)*dT.row(j0));
			rhs.setrow(j0, -dT.row(j0)/map.rz.row(j0)) ;
			if (n) {
				op->bc_bot1_add_l(n, "log_T", "T", -1./map.rz.row(j0-1), D.block(n-1).row(-1));
				op->bc_bot1_add_d(n, "log_T", "rz", 1./map.rz.row(j0-1)/map.rz.row(j0-1)*dT.row(j0-1));
				rhs.setrow(j0, rhs.row(j0) + dT.row(j0-1)/map.rz.row(j0-1));
			}
		}
		if (n < ndomains - 1) {
			op->bc_top1_add_d(n, "log_T", "T", ones(1,nth));
			op->bc_top2_add_d(n, "log_T", "T", -ones(1,nth));
			rhs.setrow(j1, -T.row(j1) + T.row(j1+1));
		}
		else {
			op->bc_top1_add_d(n, "log_T", "T", ones(1,nth));
			op->bc_top1_add_d(n, "log_T", "Ts", -ones(1,nth));
			rhs.setrow(j1, -T.row(j1) + Ts);
		}
		j0 += map.npts[n];
	}
	op->set_rhs("log_T", rhs);

	rhs = zeros(ndomains, 1);
	for (int n=0; n<ndomains; n++) {
		if (n == 0) {
			op->bc_bot2_add_r(n, "Lambda", "T", ones(1,1), map.It/2.);
			rhs(n) = 1. - (T.row(0), map.It)(0)/2.;
			// we use an average value of the central temperature 2=int_0^pi sintheta dtheta
		}
		else {
			op->bc_bot2_add_d(n,"Lambda","Lambda",ones(1,1));
			op->bc_bot1_add_d(n,"Lambda","Lambda",-ones(1,1));
		}
	}
	op->set_rhs("Lambda", rhs);

}

/// \brief Writes the equations for the dimensional quantities (T_c, rho_c, R, etc.)

void star2d::solve_dim(solver *op) {
	int n,j0;
	matrix q,rhs;

	rhs=zeros(ndomains,1);
	j0=0;
// Expression of the mass integral m=intvol rho 2*pi*r^2*rz*sin(th)*dth
	for(n=0;n<ndomains;n++) {
		op->bc_bot2_add_d(n,"m","m",ones(1,1));
		//rho
		op->bc_bot2_add_lri(n,"m","rho",-2*PI*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,(r*r*map.rz).block(j0,j0+map.gl.npts[n]-1,0,-1));
		//r (rz)
		op->bc_bot2_add_lri(n,"m","r",-2*PI*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,(2*r*map.rz*rho).block(j0,j0+map.gl.npts[n]-1,0,-1));
		op->bc_bot2_add_lri(n,"m","rz",-2*PI*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,(r*r*rho).block(j0,j0+map.gl.npts[n]-1,0,-1));

		if(n) op->bc_bot1_add_d(n,"m","m",-ones(1,1));
		j0+=map.gl.npts[n];
	}
	op->set_rhs("m",rhs);

// From the equation of state dln(rho_c)/dln(p_c) = 1/khi_rho(0)
	for(n=0;n<ndomains;n++) {
		op->add_d(n,"log_rhoc","log_pc",1./eos.chi_rho(0)*ones(1,1));
		op->add_d(n,"log_rhoc","log_Tc",-eos.d(0)*ones(1,1));
	}

// pi_c= 4*pi*G*rho_c^2/P_c
	rhs=zeros(ndomains,1);
	for(n=0;n<ndomains;n++) {
		if(n==ndomains-1) {
			op->add_d(n,"log_pc","log_pc",ones(1,1));
			op->add_d(n,"log_pc","pi_c",ones(1,1)/pi_c);
			op->add_d(n,"log_pc","log_rhoc",-2*ones(1,1));
			op->add_d(n,"log_pc","log_R",-2*ones(1,1));
		} else {
			op->bc_top1_add_d(n,"log_pc","log_pc",ones(1,1));
			op->bc_top2_add_d(n,"log_pc","log_pc",-ones(1,1));
		}
	}
	op->set_rhs("log_pc",rhs);

// T_c = rho_cR^2/Lambda
	rhs=zeros(ndomains,1);
	for(n=0;n<ndomains;n++) {
		if(n==ndomains-1) {
			op->add_d(n,"log_Tc","log_Tc",ones(1,1));
			op->add_d(n,"log_Tc","log_rhoc",-ones(1,1));
			op->add_d(n,"log_Tc","Lambda",ones(1,1)/Lambda);
			op->add_d(n,"log_Tc","log_R",-2*ones(1,1));
		} else {
			op->bc_top1_add_d(n,"log_Tc","log_Tc",ones(1,1));
			op->bc_top2_add_d(n,"log_Tc","log_Tc",-ones(1,1));
		}
	}
	op->set_rhs("log_Tc",rhs);

// M = rho_c*R^3*m
	rhs=zeros(ndomains,1);
	for(n=0;n<ndomains;n++) {
		if(n==ndomains-1) {
			op->add_d(n,"log_R","log_R",3*ones(1,1));
			op->add_d(n,"log_R","m",1/m*ones(1,1));
			op->add_d(n,"log_R","log_rhoc",ones(1,1));
		} else {
			op->bc_top1_add_d(n,"log_R","log_R",ones(1,1));
			op->bc_top2_add_d(n,"log_R","log_R",-ones(1,1));
		}
	}
	op->set_rhs("log_R",rhs);
}


/// \brief Part of the Jacobian associated with the mapping
/// There are geometrical relations and physical relations placing the physical interfaces
/// on surfaces zeta=cst. Option also to define the stellar surface as a tau=cst or P=cst
void star2d::solve_map(solver *op) {
	int n,j0;
	matrix Ri,rhs,TT,q;


	matrix gamma;
	gamma=(sqrt(1+map.rt*map.rt/r/r)).row(-1);

	n=ndomains-1;
	op->bc_top1_add_d(n,"gamma","gamma",2*gamma);
	op->bc_top1_add_d(n,"gamma","r",(2*map.rt*map.rt/r/r/r).row(-1));
	op->bc_top1_add_r(n,"gamma","r",(-2*map.rt/r/r).row(-1),Dt);
	op->set_rhs("gamma",zeros(1,nth));

	rhs=zeros(ndomains,1);
	for(n=0;n<ndomains;n++) {
		op->bc_top1_add_d(n,"deta","deta",ones(1,1));
		op->bc_top1_add_d(n,"deta","eta",ones(1,1));
		op->bc_top2_add_d(n,"deta","eta",-ones(1,1));
	}
	op->set_rhs("deta",rhs);

	rhs=zeros(ndomains,nth);
	for(n=0;n<ndomains;n++) {
		op->bc_top1_add_d(n,"dRi","dRi",ones(1,nth));
		op->bc_top1_add_d(n,"dRi","Ri",ones(1,nth));
		op->bc_top2_add_d(n,"dRi","Ri",-ones(1,nth));
	}
	op->set_rhs("dRi",rhs);

	map.leg.eval_00(map.leg.th,0,TT);

	rhs=zeros(ndomains+1,1);
	op->add_d(0,"eta","eta",ones(1,1));
	for(n=1;n<ndomains;n++) {
		op->add_d(n,"eta","eta",ones(1,1));
		op->add_r(n,"eta","Ri",-ones(1,1),TT);
	}
	op->add_d(ndomains,"eta","eta",ones(1,1));
	op->set_rhs("eta",rhs);

	Ri=map.R;
	rhs=zeros(ndomains+1,nth);

	op->add_d(0,"Ri","Ri",ones(1,nth));
//for(n=1;n<=ndomains;n++) op->add_d(n,"Ri","Ri",ones(1,nth));op->set_rhs("Ri",rhs);return;
	j0=map.gl.npts[0];
	for(n=1;n<ndomains;n++) {
		if(n==conv) {
			symbolic S;
			sym p_,s_,eq;
			p_=S.regvar("p");
			s_=S.regvar("s");
			S.set_map(map);
			S.set_value("p",p);
			S.set_value("s",entropy());

			eq=(grad(p_),grad(s_))/S.r/S.r;
			eq.bc_bot2_add(op,n,"Ri","p",ones(1,nth));
			eq.bc_bot2_add(op,n,"Ri","s",ones(1,nth));
			eq.bc_bot2_add(op,n,"Ri","r",ones(1,nth));

			rhs.setrow(n,-eq.eval().row(j0));
		} else {
			#ifdef T_CONSTANT_DOMAINS
			map.leg.eval_00(map.leg.th,zeros(1,nth),TT);
			q=zeros(1,nth);
			q(0,nth-1)=1;
			matrix delta;
			delta=zeros(1,map.gl.npts[n]);delta(0)=1;delta(-1)=-1;
			op->bc_bot2_add_lr(n,"Ri",LOG_PRES,q*domain_weight[n],delta,TT);
			delta=zeros(1,map.gl.npts[n-1]);delta(0)=1;delta(-1)=-1;
			op->bc_bot1_add_lr(n,"Ri",LOG_PRES,-q*domain_weight[n-1],delta,TT);
			rhs.setrow(n,( (log(PRES.row(j0+map.gl.npts[n]-1)) - log(PRES.row(j0))) *domain_weight[n]
			- (log(PRES.row(j0-1)) - log(PRES.row(j0-map.gl.npts[n-1])))*domain_weight[n-1] , TT) );
			op->bc_bot2_add_d(n,"Ri",LOG_PRES,(1-q));
			op->bc_bot2_add_r(n,"Ri",LOG_PRES,q-1,TT);
			rhs.setrow(n,rhs.row(n)+(1-q)*(-log(PRES.row(j0))+log((PRES.row(j0),TT))));
			#else
			matrix delta;
			delta=zeros(1,map.gl.npts[n]);delta(0)=1;delta(-1)=-1;
			op->bc_bot2_add_l(n,"Ri",LOG_PRES,ones(1,nth),delta);
			delta=zeros(1,map.gl.npts[n-1]);delta(0)=1;delta(-1)=-1;
			op->bc_bot1_add_l(n,"Ri",LOG_PRES,-ones(1,nth),delta);
			rhs.setrow(n,log(PRES.row(j0+map.gl.npts[n]-1))-log(PRES.row(j0))
			-log(PRES.row(j0-1))+log(PRES.row(j0-map.gl.npts[n-1])));
			#endif
		}
		j0+=map.gl.npts[n];
	}
	n=ndomains;
	map.leg.eval_00(map.leg.th,zeros(1,nth),TT);
	q=zeros(1,nth);
	q(0,nth-1)=1;
	op->bc_bot2_add_r(n,"Ri","Ri",q,TT);
	rhs.setrow(n,q*(-(Ri.row(n),TT)+1));
	op->bc_bot1_add_d(n,"Ri","p",(1-q));
	#ifndef PHOTOSPHERE
	// Isobar
	op->bc_bot1_add_r(n,"Ri","p",-(1-q),TT);
	rhs.setrow(n,rhs.row(n)+(1-q)*(-p.row(-1)+(p.row(-1),TT)));
	#else
	// Photosphere
	op->bc_bot1_add_d(n,"Ri","ps",-(1-q)*PHOTOSPHERE);
	rhs.setrow(n,rhs.row(n)+(1-q)*PHOTOSPHERE*(-p.row(-1)+ps));
	op->bc_bot1_add_r(n,"Ri","p",-(1-q)*(1.-PHOTOSPHERE),TT);
	rhs.setrow(n,rhs.row(n)+(1-q)*(1.-PHOTOSPHERE)*(-p.row(-1)+(p.row(-1),TT)));
	#endif

	op->set_rhs("Ri",rhs);

}




/// \brief Equation setting the equatorial angular velocity
void star2d::solve_Omega(solver *op) {
	int n;
	matrix rhs;

	rhs=zeros(ndomains+1,1);
	for(n=0;n<ndomains;n++) {
		op->bc_top1_add_d(n,"Omega","Omega",ones(1,1));
		op->bc_top2_add_d(n,"Omega","Omega",-ones(1,1));
	}
	matrix TT;
	double Req;
	Req=map.leg.eval_00(rex.row(0),PI/2,TT)(0);
	n=ndomains;
	op->bc_bot1_add_d(n,"Omega","Omega",ones(1,1));
	op->bc_bot1_add_d(n,"Omega","m",-ones(1,1)*Omega_bk*sqrt(pi_c/Req/Req/Req/4./PI)/sqrt(m)/2.);
	op->bc_bot1_add_d(n,"Omega","pi_c",-ones(1,1)*Omega_bk*sqrt(m/pi_c/Req/Req/Req/4./PI)/2.);
	op->bc_bot2_add_r(n,"Omega","Ri",ones(1,1)*Omega_bk*sqrt(pi_c*m/Req/Req/Req/Req/Req/4./PI)*3./2.,TT);
	rhs(n)=-Omega+Omega_bk*sqrt(pi_c*m/Req/Req/Req/4./PI);
	op->set_rhs("Omega",rhs);

}



/// \brief Equation giving the effective surface gravity gsup
/// gsup=(-\vn\cdot\grad P)/rho
void star2d::solve_gsup(solver *op) {
	matrix q,g;
	int n=ndomains-1;
	matrix &rt=map.rt;

	g=gsup();

	op->bc_top1_add_d(n,"gsup","gsup",ones(1,nth));
	op->bc_top1_add_d(n,"gsup","log_pc",-g);
	op->bc_top1_add_d(n,"gsup","log_rhoc",g);
	op->bc_top1_add_d(n,"gsup","log_R",g);
	op->bc_top1_add_d(n,"gsup","rho",g/rho.row(-1));

	symbolic S;
	sym p_,gamma_,eq;
	sym_vec n_(COVARIANT);
	gamma_=S.regvar("gamma");
	p_=S.regvar("p");
	S.set_map(map);
	S.set_value("p",p);
	S.set_value("gamma",sqrt(1+rt*rt/r/r));  // Bug corrected 14/6/2015, MR
	//S.set_value("gamma",sqrt(1-rt*rt/r/r));  // Looks like an error, should be sqrt(1+rt*rt/r/r)
	n_(0)=Dz(S.r)/gamma_;n_(1)=0*S.one;n_(2)=0*S.one;

	eq=(n_,grad(p_));
	eq.bc_top1_add(op,n,"gsup","p",pc/R/rhoc/rho.row(-1));
	eq.bc_top1_add(op,n,"gsup","gamma",pc/R/rhoc/rho.row(-1));
	eq.bc_top1_add(op,n,"gsup","r",pc/R/rhoc/rho.row(-1));

	op->set_rhs("gsup",zeros(1,nth));


}

/// \brief Equation setting the surface effective temperature
/// Derived from sigma T_e^4 = -xi\vn\cdot\gradT
void star2d::solve_Teff(solver *op) {
	matrix q,Te,F;
	int n=ndomains-1;
	matrix &rt=map.rt;

	Te=Teff();
	F=SIG_SB*pow(Te,4);

	op->bc_top1_add_d(n,"Teff","Teff",4*SIG_SB*pow(Te,3));
	op->bc_top1_add_d(n,"Teff","log_Tc",-F);
	op->bc_top1_add_d(n,"Teff","log_R",F);
	op->bc_top1_add_d(n,"Teff","opa.xi",-F/opa.xi.row(-1));

	symbolic S;
	sym T_,gamma_,eq;
	sym_vec n_(COVARIANT);
	gamma_=S.regvar("gamma");
	T_=S.regvar("T");
	S.set_map(map);
	S.set_value("T",T);
	S.set_value("gamma",sqrt(1+rt*rt/r/r)); // Bug corrected 14/6/2015, MR
	//S.set_value("gamma",sqrt(1-rt*rt/r/r)); // Looks like an error, should be sqrt(1+rt*rt/r/r)
	n_(0)=Dz(S.r)/gamma_;n_(1)=0*S.one;n_(2)=0*S.one;

	eq=(n_,grad(T_));
	eq.bc_top1_add(op,n,"Teff","T",Tc/R*opa.xi.row(-1));
	eq.bc_top1_add(op,n,"Teff","gamma",Tc/R*opa.xi.row(-1));
	eq.bc_top1_add(op,n,"Teff","r",Tc/R*opa.xi.row(-1));

	op->set_rhs("Teff",zeros(1,nth));



}

/// \brief Equation setting the 'simple' atmosphere model equations
/// 'simple' == the polytropic model for the atmosphere
/// checked but approximate on ... "p",-1/p.row(-1)/(n_atm+1) ...
void star2d::solve_atm(solver *op) {

	if(!strcmp(atm.name,"simple")) {
		int n=ndomains-1, n_atm=3;

		op->bc_top1_add_d(n,"ps","ps",1/ps);
		op->bc_top1_add_d(n,"ps","gsup",-1/gsup());
		op->bc_top1_add_d(n,"ps","opa.k",1/opa.k.row(-1));
		op->bc_top1_add_d(n,"ps","log_pc",ones(1,nth));
		op->set_rhs("ps",zeros(1,nth));

		op->bc_top1_add_d(n,"Ts","Ts",1/Ts);
		op->bc_top1_add_d(n,"Ts","p",-1/p.row(-1)/(n_atm+1));
		op->bc_top1_add_d(n,"Ts","ps",1/ps/(n_atm+1));
		op->bc_top1_add_d(n,"Ts","Teff",-1/Teff());
		op->bc_top1_add_d(n,"Ts","log_Tc",ones(1,nth));
		op->set_rhs("Ts",zeros(1,nth));
		return;
	}


	matrix q;
	int n=ndomains-1;

	op->bc_top1_add_d(n,"ps","ps",1/ps);
	op->bc_top1_add_d(n,"ps","gsup",-atm.dlnps_lng/gsup());
	op->bc_top1_add_d(n,"ps","Teff",-atm.dlnps_lnTeff/Teff());
	op->bc_top1_add_d(n,"ps","log_pc",ones(1,nth));
	op->set_rhs("ps",zeros(1,nth));

	op->bc_top1_add_d(n,"Ts","Ts",1/Ts);
	op->bc_top1_add_d(n,"Ts","p",-0.25/p.row(-1));
	op->bc_top1_add_d(n,"Ts","ps",0.25/ps);
	op->bc_top1_add_d(n,"Ts","gsup",-atm.dlnTs_lng/gsup());
	op->bc_top1_add_d(n,"Ts","Teff",-atm.dlnTs_lnTeff/Teff());
	op->bc_top1_add_d(n,"Ts","log_Tc",ones(1,nth));
	op->set_rhs("Ts",zeros(1,nth));

}




