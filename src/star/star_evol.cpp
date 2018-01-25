#include "ester-config.h"
#include "utils.h"
#include "star.h"
#include <omp.h>

star_evol::star_evol() : star2d() {

	check_map_enable = false;

	// Not needed, just to avoid warnings of the compiler
	delta = 0.;
	lnR0 = log(R);
	drhocdX = 0;
	lnrhoc0 = 0.;
	lnpc0 = 0.;
	lnTc0 = 0.;
	mH0 = 0.;
	mH = 0.;
}

star_evol::star_evol(const star2d &A) : star2d(A) {

	check_map_enable = false;

	// Not needed, just to avoid warnings of the compiler
	delta = 0.;
	lnR0 = log(R);
	drhocdX = 0;
	lnrhoc0 = 0.;
	lnpc0 = 0.;
	lnTc0 = 0.;
	mH0 = 0.;
	mH = 0.;
}

star_evol::star_evol(const star_evol &A) : star2d(A) {
	copy(A);
}

star_evol &star_evol::operator=(const star_evol &A) {
	star2d::copy(A);
	copy(A);
	return *this;
}

void star_evol::copy(const star_evol &A) {
	Xprev = A.Xprev;
	r0 = A.r0;
	rho0 = A.rho0;
	T0 = A.T0;
	lnp0 = A.lnp0;
	lnR0 = A.lnR0;
	drhocdX = A.drhocdX;
	lnrhoc0 = A.lnrhoc0;
	lnTc0 = A.lnTc0;
	lnpc0 = A.lnpc0;
	delta = A.delta;
	check_map_enable = A.check_map_enable;
	mH0 = A.mH0;
	mH = A.mH;
	age = A.age;

	dXdt = A.dXdt;
	drhodt = A.drhodt;
	dpdt = A.dpdt;
	dwdt = A.dwdt;
	dTdt = A.dTdt;
	dphidt = A.dphidt;
}

void star_evol::read_vars(INFILE *fp) {

	if(fp->read("dXdt", &dXdt)) dXdt = zeros(nr, nth);
	if(fp->read("drhodt", &drhodt)) drhodt = zeros(nr, nth);
	if(fp->read("dpdt", &dpdt)) dpdt = zeros(nr, nth);
	if(fp->read("dwdt", &dwdt)) dwdt = zeros(nr, nth);
	if(fp->read("dTdt", &dTdt)) dTdt = zeros(nr, nth);
	if(fp->read("dphidt", &dphidt)) dphidt = zeros(nr, nth);

	fill();

}

void star_evol::write_vars(OUTFILE *fp) const {

	fp->write("dXdt", &dXdt);
	fp->write("drhodt", &drhodt);
	fp->write("dpdt", &dpdt);
	fp->write("dwdt", &dwdt);
	fp->write("dTdt", &dTdt);
	fp->write("dphidt", &dphidt);
}

void star_evol::fill() {

	star2d::fill();

	Omega_bk=Omega/Omegac;
	mH = 2*PI*rhoc*R*R*R*(map.I, rho*comp["H"]*r*r*map.rz, map.It)(0);


}

void star_evol::calcTimeDerivs() {
	static symbolic S;
	static sym derivt;
	static bool sym_inited = false;
	if (!sym_inited) {
		sym val = S.regvar("val");
		sym val0 = S.regvar("val0");
		sym lnR = S.regconst("log_R");
		sym R = exp(lnR);
		sym r = S.r;
		sym r0 = S.regvar("r0");
		sym lnR0 = S.regconst("log_R0");
		sym delta = S.regconst("delta");
		sym drdt = (r-r0)/delta;
		sym dlnRdt = (lnR-lnR0)/delta;
		derivt = (val-val0)/delta - (drdt+r*dlnRdt)/S.rz*Dz(val);
		sym_inited = true;
	}

	S.set_value("log_R", log(R)*ones(1,1));
	S.set_value("delta", delta*ones(1,1));
	S.set_value("log_R0", lnR0*ones(1,1));
	S.set_value("r0", r0);
	S.set_map(map);

	S.set_value("val", comp["H"]);
	S.set_value("val0", Xprev);
	dXdt = derivt.eval();

	S.set_value("val", rho);
	S.set_value("val0", rho0);
	drhodt = rhoc * (derivt.eval() + rho*(log(rhoc) - lnrhoc0)/delta);

	S.set_value("val", log(p));
	S.set_value("val0", lnp0);
	dpdt = pc * (p*derivt.eval() + p*(log(pc) - lnpc0)/delta);

	S.set_value("val", w);
	S.set_value("val0", w0);
	dwdt = units.Omega * (derivt.eval() + w*0.5*(log(pc) - lnpc0)/delta - w*0.5*(log(rhoc) - lnrhoc0)/delta - w*(log(R) - lnR0)/delta);

	S.set_value("val", T);
	S.set_value("val0", T0);
	dTdt = Tc * (derivt.eval() + T*(log(Tc) - lnTc0)/delta);

	S.set_value("val", phi);
	S.set_value("val0", phi0);
	dphidt = units.phi * (derivt.eval() + phi*(log(pc) - lnpc0)/delta - phi*(log(rhoc) - lnrhoc0)/delta);

}

void star_evol::init_comp() {

}

solver * star_evol::init_solver(int nvar_add) {
	return star2d::init_solver(nvar_add+4);
}

void star_evol::register_variables(solver *op) {

	star2d::register_variables(op);
	op->regvar("X");
	op->regvar("Xc");
	op->regvar("log_M");
	op->regvar("reg_cont");

}

void star_evol::write_eqs(solver *op) {
#ifdef MKL
	int mkl_threads = mkl_get_max_threads();
	mkl_set_num_threads(1);   // Faster with 1 thread !?
#endif
#ifdef THREADS
	int num_threads = 1;
	std::thread t[num_threads];
	t[0] = std::thread(&star_evol::solve_X, this, op);
#else
	solve_X(op);
#endif
	star2d::write_eqs(op);
#ifdef THREADS
	for(int i=0; i< num_threads; i++)
		t[i].join();
#endif
#ifdef MKL
	mkl_set_num_threads(mkl_threads);
#endif
}

double star_evol::update_solution(solver *op, double &h) {
	double dmax=config.newton_dmax;

	matrix dX = op->get_var("X");
	while(exist(abs(h*dX)>dmax)) h /= 2;

	double err = star2d::update_solution(op, h);

	dX = max(dX, -comp["H"]/h);
	double err2 = max(abs(dX));
	err=err2>err?err2:err;
	comp["He4"] -= h*dX;
	comp["H"] += h*dX;

	M *= exp(h*op->get_var("log_M")(0));

	return err;
}

SDIRK_solver* star_evol::init_time_solver() {
	SDIRK_solver *rk = new SDIRK_solver();
	rk->init(11, "sdirk3");
	register_variables(rk);
	return rk;
}

void star_evol::register_variables(SDIRK_solver *rk) {
	check_map_enable = true;
	check_map();
	check_map_enable = false;

	rk->regvar("X", comp.X());
	rk->regvar("lnR", log(R)*ones(1,1));
	rk->regvar("r", r);
	rk->regvar("rho", rho);
	rk->regvar("log_rhoc", log(rhoc)*ones(1,1));
	rk->regvar("T", T);
	rk->regvar("log_Tc", log(Tc)*ones(1,1));
	rk->regvar("log_p", log(p));
	rk->regvar("log_pc", log(pc)*ones(1,1));
	rk->regvar("w", w);
	rk->regvar("phi", phi);
}

void star_evol::init_step(SDIRK_solver *rk) {
	delta = rk->get_delta();
	age = rk->get_t();

	Xprev = rk->get_var("X");
	r0 = rk->get_var("r");
	lnR0 = rk->get_var("lnR")(0);
	rho0 = rk->get_var("rho");
	lnrhoc0 = rk->get_var("log_rhoc")(0);
	T0 = rk->get_var("T");
	lnTc0 = rk->get_var("log_Tc")(0);
	lnp0 = rk->get_var("log_p");
	lnpc0 = rk->get_var("log_pc")(0);
	w0 = rk->get_var("w");
	phi0 = rk->get_var("phi");
}

void star_evol::finish_step(SDIRK_solver *rk, int state) {
	if (state == RK_STEP) {
		check_map_enable = true;
		check_map();
		check_map_enable = false;
	}

	rk->set_var("X", comp.X());
	rk->set_var("r", map.r);
	rk->set_var("lnR", log(R)*ones(1,1));
	rk->set_var("rho", rho);
	rk->set_var("log_rhoc", log(rhoc)*ones(1,1));
	rk->set_var("T", T);
	rk->set_var("log_Tc", log(Tc)*ones(1,1));
	rk->set_var("log_p", log(p));
	rk->set_var("log_pc", log(pc)*ones(1,1));
	rk->set_var("w", w);
	rk->set_var("phi", phi);

}

void star_evol::reset_time_solver(SDIRK_solver *rk) {
	finish_step(rk, -1);
	rk->reset();
}


void star_evol::check_map() {
	if (!check_map_enable) return;
	else star2d::check_map();
}

void star_evol::interp(remapper *red) {
	star2d::interp(red);
}

void star_evol::calc_units() {
	star2d::calc_units();
	units.v = R / MYR;
}

void star_evol::solve_definitions(solver *op) {
    DEBUG_FUNCNAME;

    star2d::solve_definitions(op);

    composition_map comp2(comp);
    comp2["H"] += 1e-4;
    comp2["He4"] -= 1e-4;

    eos_struct eos2;
    matrix rho2;
    strcpy(eos2.name, eos.name);
    eos_calc(comp2.X(), Z0, T*Tc, p*pc, rho2, eos2);
    matrix drhodX = (rho2/rhoc - rho)/1e-4;
    drhocdX = (rho2(0,0) - rhoc) / 1e-4;
    op->add_d("rho", "X", drhodX);
    //op->add_d("eos.cp", "X", (eos2.cp - eos.cp) / 1e-4);
    //op->add_d("eos.del_ad", "X", (eos2.del_ad - eos.del_ad) / 1e-4);

    opa_struct opa2;
    strcpy(opa2.name, opa.name);
    opa_calc(comp2.X(), Z0, T*Tc, rho*rhoc, opa2);
    matrix dxidX = (opa2.xi - opa.xi) / 1e-4;
    op->add_d("opa.xi", "X", dxidX);

    nuc_struct nuc2;
    strcpy(nuc2.name, nuc.name);
    nuc_calc(comp2, T*Tc, rho*rhoc, nuc2);
    matrix depsdX = (nuc2.eps - nuc.eps) / 1e-4;
    op->add_d("nuc.eps", "X", depsdX);

    for(int n=0; n<ndomains; n++) {
    	op->add_d(n, "log_rhoc", "Xc", drhocdX/rhoc*ones(1,1));
    }
}


void star_evol::solve_dim(solver *op) {

	star2d::solve_dim(op);

	op->add_d(ndomains-1, "log_R", "log_M", -ones(1,1));

}

void star_evol::solve_Omega(solver *op) {
	DEBUG_FUNCNAME;

	matrix rhs;
	int n;
	
	rhs=zeros(ndomains, 1);
	for(n=0; n<ndomains-1; n++) {
		op->bc_top1_add_d(n, "Omega", "Omega", ones(1,1));
		op->bc_top2_add_d(n, "Omega", "Omega", -ones(1,1));
	}

	n = ndomains-1;

	if (Omega == 0) {
		op->bc_top1_add_d(n, "Omega", "Omega", ones(1,1));
	}
	else {
		matrix TT;
		map.leg.eval_00(th, PI/2, TT);
		op->bc_top1_add_d(n, "Omega", "Omega", ones(1,1));
		op->bc_top1_add_r(n, "Omega", "w", -ones(1,1), TT);
		rhs(n) = -Omega + (w.row(-1), TT)(0);
	}
	op->set_rhs("Omega", rhs);
}

void star_evol::solve_X(solver *op) {
	DEBUG_FUNCNAME;
    double Q=(4*HYDROGEN_MASS-AMASS["He4"]*UMA)*C_LIGHT*C_LIGHT;
    double diff_coeff_conv = 1e11; // PARAMETER

    static symbolic S;
    static sym eq, flux, gradX;
    static bool sym_inited = false;
    if (!sym_inited) {
    	sym X = S.regvar("X");
    	sym eps = S.regvar("nuc.eps");
    	sym rho = S.regvar("rho");
    	sym vr = S.regvar("vr");
    	sym vt = S.regvar("vt");
    	sym lnR = S.regconst("log_R");
    	sym R = exp(lnR);
    	sym r = S.r;
    	sym X0 = S.regvar("X0");
    	sym r0 = S.regvar("r0");
    	sym lnR0 = S.regconst("log_R0");
    	sym delta = S.regconst("delta");
    	sym Dv = S.regconst("diffusion_v");
    	sym Dh = S.regconst("diffusion_h");
    	sym sin_vangle = S.regvar("sin_vangle");
    	sym cos_vangle = S.regvar("cos_vangle");
    	sym MYR = S.regconst("MYR");
    	sym mH = S.regconst("mH");
    	sym Q = S.regconst("Q");

    	sym_vec phivec(COVARIANT);
    	phivec(0) = 0*S.one; phivec(1) = 0*S.one; phivec(2) = S.r*sin(S.theta);
    	sym_vec rvec = grad(S.r);
    	sym_vec thvec = cross(phivec, rvec);

    	sym_vec rhoV = rho*(vr*rvec + vt*thvec);

    	sym drdt = (r-r0)/delta;
    	sym dlnRdt = (lnR-lnR0)/delta;
    	sym dXdt = (X-X0)/delta - (drdt+r*dlnRdt)/S.rz*Dz(X);

    	sym_vec v_vec = cos_vangle*rvec - sin_vangle*thvec;
    	sym_vec h_vec = sin_vangle*rvec + cos_vangle*thvec;
    	sym_tens D = tensor(v_vec, v_vec)*Dv + tensor(h_vec, h_vec)*Dh;

    	eq = dXdt + 4*mH*MYR/Q*eps - MYR/R/R * div(rho*(D, grad(X)))/rho + (rhoV, grad(X))/rho;
    	sym_vec nvec = grad(S.zeta);
    	nvec = nvec / sqrt((nvec, nvec));
    	flux = - MYR/R/R*(nvec, (D, grad(X)));
    	gradX = (nvec, grad(X));
    	sym_inited = true;
    }

    S.set_value("X", comp.X());
    S.set_value("nuc.eps", nuc.eps);
    S.set_value("rho", rho);
    S.set_value("vr", vr);
    S.set_value("vt", vt, 11);
    S.set_value("log_R", log(R)*ones(1,1));
    S.set_value("delta", delta*ones(1,1));
    S.set_value("X0", Xprev);
    S.set_value("log_R0", lnR0*ones(1,1));
    S.set_value("r0", r0);
    S.set_value("MYR", MYR*ones(1,1));
    S.set_value("mH", HYDROGEN_MASS*ones(1,1));
    S.set_value("Q", Q*ones(1,1));
    S.set_value("sin_vangle", sin(vangle), 11);
    S.set_value("cos_vangle", cos(vangle));

    S.set_map(map);

    matrix diff_v = ones(nr, nth) * diffusion_v;
    matrix diff_h = ones(nr, nth) * diffusion_h;
    if (conv) {
    	int nc = 0;
    	for (int n = 0; n < conv; n++) nc += map.npts[n];
    	diff_v.setblock(0, nc-1, 0, -1, ones(nc, nth) * diff_coeff_conv);
    	diff_h.setblock(0, nc-1, 0, -1, ones(nc, nth) * diff_coeff_conv);
    }

    S.set_value("diffusion_v", diff_v);
    S.set_value("diffusion_h", diff_h);

    eq.add(op, "X", "X");
    eq.add(op, "X", "nuc.eps");
    eq.add(op, "X", "rho");
    eq.add(op, "X", "log_R");
    eq.add(op, "X", "r");
    eq.add(op, "X", "vr");
    eq.add(op, "X", "vt");
    eq.add(op, "X", "sin_vangle");
    eq.add(op, "X", "cos_vangle");

    matrix rhs=-eq.eval();
    matrix X = comp.X();
    matrix dX = (D, X);
    matrix flux_val = flux.eval();
    matrix gradX_val = gradX.eval();

    int j0 = 0;
    for (int n = 0; n < ndomains; n++) {
    	int j1 = j0 + map.npts[n] - 1;

    	if (n == conv && conv > 0) {
    		flux.bc_bot2_add(op, n, "X", "X");
    		flux.bc_bot2_add(op, n, "X", "r");
    		flux.bc_bot1_add(op, n, "X", "X", -ones(1,nth));
    		flux.bc_bot1_add(op, n, "X", "r", -ones(1,nth));
    		rhs.setrow(j0, - flux_val.row(j0) + flux_val.row(j0-1));
    	}
    	else {
    		op->bc_bot2_add_l(n, "X", "X", 1./map.rz.row(j0), D.block(n).row(0));
    		op->bc_bot2_add_d(n, "X", "rz", -1./map.rz.row(j0)/map.rz.row(j0)*dX.row(j0));
    		rhs.setrow(j0, -dX.row(j0)/map.rz.row(j0)) ;
    		if (n) {
    			op->bc_bot1_add_l(n, "X", "X", -1./map.rz.row(j0-1), D.block(n-1).row(-1));
    			op->bc_bot1_add_d(n, "X", "rz", 1./map.rz.row(j0-1)/map.rz.row(j0-1)*dX.row(j0-1));
    			rhs.setrow(j0, rhs.row(j0) + dX.row(j0-1)/map.rz.row(j0-1));
    		}
    	}
    	if (n < ndomains - 1) {
    		if (n == conv-1) {
    			op->bc_top1_add_d(n, "X", "X", ones(1,nth));
    			op->bc_top2_add_d(n, "X", "X", -ones(1,nth));
    			rhs.setrow(j1, -X.row(j1) + X.row(j1+1));
    		}
    		else {
    			op->bc_top1_add_d(n, "X", "X", ones(1,nth));
    			op->bc_top2_add_d(n, "X", "X", -ones(1,nth));
    			rhs.setrow(j1, -X.row(j1) + X.row(j1+1));
    		}
    	}
    	else {
    		gradX.bc_top1_add(op, n, "X", "X");
    		gradX.bc_top1_add(op, n, "X", "r");
    		rhs.setrow(j1, -gradX_val.row(-1));
    	}
    	j0 += map.npts[n];
    }

    op->set_rhs("X",rhs);


    rhs=zeros(ndomains,1);
    for(int n = 0; n < ndomains; n++) {
    	if(n == 0) {
    		op->bc_bot2_add_d(n, "Xc", "Xc", ones(1,1));
    		op->bc_bot2_add_r(n, "Xc", "X", -ones(1,1), map.It/2.);
    	} else {
    		op->bc_bot2_add_d(n, "Xc", "Xc", ones(1,1));
    		op->bc_bot1_add_d(n, "Xc", "Xc", -ones(1,1));
    	}
    }
    op->set_rhs("Xc", rhs);

}

void star_evol::solve_cont(solver *op) {
	DEBUG_FUNCNAME;

	star2d::solve_cont(op);

	static symbolic S;
	static sym eq, flux, reg_cont;
	static bool sym_inited = false;
	if (!sym_inited) {
		sym rho = S.regvar("rho");
		sym vr = S.regvar("vr");
		sym lnrhoc = S.regconst("log_rhoc");
		sym lnR = S.regconst("log_R");
		sym R = exp(lnR);
		sym r = S.r;
		sym rho0 = S.regvar("rho0");
		sym lnrhoc0 = S.regconst("log_rhoc0");
		sym r0 = S.regvar("r0");
		sym lnR0 = S.regconst("log_R0");
		sym delta = S.regconst("delta");
		sym drdt = (r-r0)/delta;
		sym dlnRdt = (lnR-lnR0)/delta;
		sym drhodt = (rho-rho0)/delta - (drdt+r*dlnRdt)/S.rz*Dz(rho);
		sym dlnrhocdt = (lnrhoc - lnrhoc0)/delta;
		eq = drhodt + rho*dlnrhocdt;
		sym_vec nvec = grad(S.zeta);
		nvec = nvec / sqrt((nvec, nvec));
		sym_vec rvec = grad(r);
		flux = -rho * (rvec, nvec) * (r*dlnRdt + drdt);

		reg_cont = 3*S.Dz(rho*vr)/S.rz + drhodt + rho*dlnrhocdt;

		sym_inited = true;
	}

	S.set_value("rho", rho);
	S.set_value("vr", vr);
	S.set_value("log_rhoc", log(rhoc)*ones(1,1));
	S.set_value("log_R", log(R)*ones(1,1));
	S.set_value("rho0", rho0);
	S.set_value("log_rhoc0", lnrhoc0*ones(1,1));
	S.set_value("r0", r0);
	S.set_value("log_R0", lnR0*ones(1,1));
	S.set_value("delta", delta*ones(1,1));

	S.set_map(map);

	eq.add(op, "vr", "rho");
	eq.add(op, "vr", "log_rhoc");
	eq.add(op, "vr", "log_R");
	eq.add(op, "vr", "r");

	matrix rhs = -eq.eval();
	matrix flux_val = flux.eval();

	int j0 = 0;
	for (int n = 0; n < ndomains; n++) {
		int j1 = j0 + map.npts[n] - 1;

		if (n == 0) {
			rhs.setrow(j0, zeros(1, nth));
		}

		if (n == conv-1) {
			flux.bc_top1_add(op, n, "vr", "rho");
			flux.bc_top1_add(op, n, "vr", "log_R");
			flux.bc_top1_add(op, n, "vr", "r");
			flux.bc_top2_add(op, n, "vr", "rho", -ones(1, nth));
			flux.bc_top2_add(op, n, "vr", "log_R", -ones(1, nth));
			flux.bc_top2_add(op, n, "vr", "r", -ones(1, nth));
			rhs.setrow(j1, -flux_val.row(j1) + flux_val.row(j1+1));
		}
		else if (n == ndomains - 1) {
			flux.bc_top1_add(op, n, "vr", "rho");
			flux.bc_top1_add(op, n, "vr", "log_R");
			flux.bc_top1_add(op, n, "vr", "r");
			rhs.setrow(j1, -flux_val.row(j1));
		}
		else {
			rhs.setrow(j1, zeros(1, nth));
		}
		j0 += map.npts[n];
	}

	op->set_rhs("vr", op->get_rhs("vr") + rhs);

	for(int n = 0; n < ndomains; n++) {
		if(n == 0) {
			reg_cont.bc_bot2_add(op, 0, "reg_cont", "vr");
			reg_cont.bc_bot2_add(op, 0, "reg_cont", "rho");
			reg_cont.bc_bot2_add(op, 0, "reg_cont", "r");
			reg_cont.bc_bot2_add(op, 0, "reg_cont", "log_rhoc");
			op->bc_bot2_add_d(0, "reg_cont", "reg_cont", -ones(1, nth));
		} else {
			op->bc_bot2_add_d(n, "reg_cont", "reg_cont", ones(1,nth));
			op->bc_bot1_add_d(n, "reg_cont", "reg_cont", -ones(1,nth));
		}
	}
	op->set_rhs("reg_cont", zeros(ndomains, nth));

	rhs=zeros(ndomains,1);
	for(int n = 0; n < ndomains; n++) {
		if(n == ndomains - 1) {
			op->bc_top1_add_r(n, "log_M", "reg_cont", ones(1,1), map.It);
			rhs(-1) = -(reg_cont.eval().row(0), map.It)(0);
		} else {
			op->bc_top1_add_d(n, "log_M", "log_M", ones(1,1));
			op->bc_top2_add_d(n, "log_M", "log_M", -ones(1,1));
		}
	}
	op->set_rhs("log_M", rhs);

}

void star_evol::solve_temp(solver *op) {
	DEBUG_FUNCNAME;
	star2d::solve_temp(op);

	static symbolic S;
	static sym eq;
	static bool sym_inited = false;

	if (!sym_inited) {
		sym_inited = true;
		sym T = S.regvar("T");
		sym lnTc = S.regconst("log_Tc");
		sym lnp = S.regvar("log_p");
		sym lnpc = S.regconst("log_pc");
		sym rho = S.regvar("rho");
		sym cp = S.regvar("eos.cp");
		sym del_ad = S.regvar("eos.del_ad");
		sym lnR = S.regconst("log_R");
		sym r = S.r;
		sym T0 = S.regvar("T0");
		sym lnTc0 = S.regconst("log_Tc0");
		sym lnp0 = S.regconst("log_p0");
		sym lnpc0 = S.regconst("log_pc0");
		sym r0 = S.regvar("r0");
		sym lnR0 = S.regconst("log_R0");
		sym delta = S.regconst("delta");
		sym MYR = S.regconst("MYR");
		sym drdt = (r-r0)/delta;
		sym dlnRdt = (lnR-lnR0)/delta;
		sym dTdt = (T-T0)/delta - (drdt+r*dlnRdt)/S.rz*Dz(T);
		sym dlnpdt = (lnp-lnp0)/delta - (drdt+r*dlnRdt)/S.rz*Dz(lnp);
		sym dlnTcdt = (lnTc-lnTc0)/delta;
		sym dlnpcdt = (lnpc-lnpc0)/delta;
		eq = rho * cp * exp(lnTc) / MYR * (dTdt + T * dlnTcdt - del_ad * T * (dlnpdt + dlnpcdt));

	}

	S.set_value("T", T);
	S.set_value("log_Tc", log(Tc)*ones(1,1));
	S.set_value("log_p", log(p));
	S.set_value("log_pc", log(pc)*ones(1,1));
	S.set_value("rho", rho);
	S.set_value("eos.cp", eos.cp);
	S.set_value("eos.del_ad", eos.del_ad);
	S.set_value("log_R", log(R)*ones(1,1));
	S.set_value("T0", T0);
	S.set_value("log_Tc0", lnTc0*ones(1,1));
	S.set_value("log_p0", lnp0);
	S.set_value("log_pc0", lnpc0*ones(1,1));
	S.set_value("r0", r0);
	S.set_value("log_R0", lnR0*ones(1,1));
	S.set_value("delta", delta*ones(1,1));
	S.set_value("MYR", MYR*ones(1,1));
	S.set_map(map);

	eq.add(op, "log_T", "T");
	eq.add(op, "log_T", "log_Tc");
	eq.add(op, "log_T", "log_p");
	eq.add(op, "log_T", "log_pc");
	eq.add(op, "log_T", "rho");
	eq.add(op, "log_T", "log_R");
	eq.add(op, "log_T", "r");

	matrix rhs = op->get_rhs("log_T");
	matrix rhs2 = -eq.eval();
	int j0 = 0;
	for (int n = 0; n < ndomains; n++) {
		int j1 = j0 + map.npts[n] - 1;
		rhs2.setrow(j0, zeros(1,nth));
		rhs2.setrow(j1, zeros(1,nth));
		j0 += map.npts[n];
	}

	op->set_rhs("log_T", rhs + rhs2);

}

void star_evol::solve_mov(solver *op) {
    DEBUG_FUNCNAME;

    star2d::solve_mov(op);

    if(Omega == 0) return;

	static bool sym_inited = false;
	static symbolic S;
	static sym eq, bc_t_add, eq_t;

	if(!sym_inited) {
		sym_inited = true;

		sym w = S.regvar("w");
		sym rho = S.regvar("rho");
		sym lnrhoc = S.regconst("log_rhoc");
		sym lnpc = S.regconst("log_pc");
		sym lnR = S.regconst("log_R");
		sym R = exp(lnR);
		sym r = S.r;
		sym w0 = S.regvar("w0");
		sym r0 = S.regvar("r0");
		sym lnrhoc0 = S.regconst("log_rhoc0");
		sym lnpc0 = S.regconst("log_pc0");
		sym lnR0 = S.regconst("log_R0");
		sym delta = S.regconst("delta");

		sym s = r*sin(S.theta);

		sym drdt = (r-r0)/delta;
		sym dlnRdt = (lnR-lnR0)/delta;
		sym dwdt = (w-w0)/delta - (drdt+r*dlnRdt)/S.rz*Dz(w);
		sym dlnrhocdt = (lnrhoc-lnrhoc0)/delta;
		sym dlnpcdt = (lnpc-lnpc0)/delta;
		eq = rho*s*s*dwdt + rho*s*s*w*(-dlnRdt + 0.5*dlnpcdt - 0.5*dlnrhocdt);

	}

	S.set_value("w",w);
	S.set_value("rho", rho);
	S.set_value("log_rhoc", log(rhoc)*ones(1, 1));
	S.set_value("log_pc", log(pc)*ones(1, 1));
	S.set_value("log_R", log(R)*ones(1, 1));
	S.set_value("w0", w0);
	S.set_value("r0", r0);
	S.set_value("log_rhoc0", lnrhoc0*ones(1, 1));
	S.set_value("log_pc0", lnpc0*ones(1, 1));
	S.set_value("log_R0", lnR0*ones(1, 1));
	S.set_value("delta", delta*ones(1,1));
	S.set_map(map);

	// Add new terms to eq. "w"

	eq.add(op, "w", "w");
	eq.add(op, "w", "rho");
	eq.add(op, "w", "log_rhoc");
	eq.add(op, "w", "log_pc");
	eq.add(op, "w", "log_R");
	eq.add(op, "w", "r");

	matrix rhs = op->get_rhs("w");
	matrix rhs2 = -eq.eval();
	int j0 = 0;
	for (int n = 0; n < ndomains; n++) {
		int j1 = j0 + map.npts[n] - 1;
		rhs2.setrow(j0, zeros(1,nth));
		rhs2.setrow(j1, zeros(1,nth));
		j0 += map.npts[n];
	}

	rhs = rhs + rhs2;
	op->reset_bc_bot(0, "w");
	op->bc_bot2_add_l(0, "w", "w", ones(1, nth), D.block(0).row(0));
	rhs.setrow(0, -(D, w).row(0));

	op->set_rhs("w", rhs);

	op->reset("Omega0");
	op->bc_bot2_add_d(0, "Omega0", "Omega0", ones(1,1));
	op->set_rhs("Omega0", zeros(1,1));

	/*
	if (1) { // Surface pressure constant in time
		op->reset("pi_c");
		rhs=zeros(ndomains,1);
		j0=0;
		for(int n=0;n<ndomains;n++) {
			if(n<ndomains-1) {
				op->bc_top1_add_d(n,"pi_c","pi_c",ones(1,1));
				op->bc_top2_add_d(n,"pi_c","pi_c",-ones(1,1));
			} else {
				matrix TT;
				map.leg.eval_00(th,0,TT);
				op->bc_top1_add_d(n,"pi_c","log_pc",ones(1,1));
				op->bc_top1_add_r(n,"pi_c","log_p",ones(1,1),TT);
				rhs(ndomains-1)=0;
			}

			j0+=map.gl.npts[n];
		}
		op->set_rhs("pi_c",rhs);
	}
	*/
}


int star_evol::remove_convective_core() {
	if (conv>1) {
		star2d::remove_convective_core();
		return 0;
	}
	double dp0 = -log(p(map.npts[0], 0));
	double dp1 = -dp0 - log(p(map.npts[0]+map.npts[1], 0));
	domain_weight[0] = domain_weight[1]*dp1/dp0;
	star2d::remove_convective_core();
	return 1;
}

double star_evol::solve(solver *op) {

	//static figure fig("/XSERVE");

	//matrix q0 = cos(vangle);
	double err = star2d::solve(op);
	calcTimeDerivs(); // For storage
	//matrix q = cos(vangle);
	//matrix dq = op->get_var("vr");
	//fig.axis(0, nr, -15, 0);
	//fig.semilogy(abs(q-q0), "b");
	//fig.hold(1);
	//fig.semilogy(abs(dq), "r");
	//fig.hold(0);


	return err;
}

