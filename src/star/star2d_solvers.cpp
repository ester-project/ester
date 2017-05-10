#include "ester-config.h"
#include"star.h"
#include<stdlib.h>
#include<sys/time.h>
#include<string.h>
#include"symbolic.h"

/// \brief Initialize star's chemical composition, equation of state, opacity,
/// nuclear reaction, mass definition, pi_c, Lambda, velocity, units, atmosphere,
/// flatness, scaled keplerian angular velocity
void star2d::fill() {
    DEBUG_FUNCNAME;
	printf("start of star2d::fill\n");
	Y0=1.-X0-Z0;
	init_comp();
	
	eq_state();
	opacity();
	nuclear();
	
	m=2*PI*(map.gl.I,rho*r*r*map.rz,map.leg.I_00)(0);

	R=pow(M/m/rhoc,1./3.);

	pi_c=(4*PI*GRAV*rhoc*rhoc*R*R)/pc;
	Lambda=rhoc*R*R/Tc;
	
	calc_veloc();	
	calc_units();

	atmosphere();

	//Omegac=sqrt(map.leg.eval_00(((Dex,phiex)/rex/map.ex.rz).row(0),PI/2)(0));
	double eps=1.-1./map.leg.eval_00(r.row(-1),PI/2)(0);
	Omegac=sqrt(pi_c*m/4/PI*(1-eps)*(1-eps)*(1-eps));

}

void star2d::init_comp() {

	printf("start of star2d::init_comp\n");
	comp=initial_composition(X0,Z0)*ones(nr,nth);

	if(!conv) return;

    int n = 0;
    for (int i=0; i<conv; i++) {
        n += map.gl.npts[i];
    }

    matrix &rz = map.rz;
    printf("AVG comp: (n: %d)\n", n);
    double Xh0 = ((map.gl.I.block(0, 0, 0, n-1)),
            (Xh*r*r*rz).block(0, n-1, 0, 0),
            (map.leg.I_00))(0);
    double Volc = ((map.gl.I.block(0, 0, 0, n-1)), (r*r*rz).block(0, n-1, 0, 0),
            (map.leg.I_00))(0);
    printf("%lf\n", Xh0/Volc);
	comp.setblock(0, n-1, 0, 0, initial_composition(Xh0/Volc, Z0)*ones(n, nth));
    for (int i=n; i<nr; i++) {
        comp.setblock(i, i, 0, 0, initial_composition(comp["H"](i, 0), Z0)*ones(1, nth));
    }
    // Xh.setblock(0, n-1, 0, 0, comp["H"].block(0, n-1, 0, 0));

#if 0
    if(stratified_comp == 0) {
        comp.setblock(0,n-1,0,-1,initial_composition(Xc*X0,Z0)*ones(n,nth));
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
#endif
}

void star2d::calc_veloc() {
    DEBUG_FUNCNAME;
// vr=rz*V^zeta vt=r*V^theta
	vr=(G,map.leg.D_11)/r+(map.rt/r+cos(th)/sin(th))/r*G;
	vr.setrow(0,zeros(1,nth));
	vr/=rho;
	vt=-(D,G)/map.rz-1./r*G;
	vt.setrow(0,zeros(1,nth));
	vt/=rho;
}

solver *star2d::init_solver(int nvar_add) {
    DEBUG_FUNCNAME;
	int nvar;
	solver *op;
	
        printf("Me voila dans init_solver of star2d\n");
	nvar=34; // include Xh
	op=new solver;
	op->init(ndomains+1,nvar+nvar_add,"full");
	
	op->maxit_ref=10;op->use_cgs=1;op->maxit_cgs=20;op->debug=0;
	op->rel_tol=1e-12;op->abs_tol=1e-20;
	register_variables(op);
	
	return op;
}

void star2d::register_variables(solver *op) {
    DEBUG_FUNCNAME;
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
	op->regvar("lum");
	op->regvar("Frad");
	op->regvar("Teff");
	op->regvar("gsup");
	op->regvar("w");
	op->regvar("G");
	op->regvar_dep("rho");
	op->regvar_dep("opa.xi");
	op->regvar_dep("opa.k");
	op->regvar_dep("nuc.eps");
	op->regvar_dep("s");
	op->regvar("gamma");
// Evolution of Xh
        op->regvar("lnXh");
	

}

/// \brief Performs one step of the Newton algorithm to compute the star's
/// internal structure.
double star2d::solve(solver *op) {
    DEBUG_FUNCNAME;
	int info[5];
	matrix rho0;
	double err,err2,h,dmax;

	if(Omega==0 && Omega_bk != 0) {
		Omega=Omega_bk*Omegac;
		w=Omega*ones(nr,nth);
	}

	check_map();

	op->reset();

	if(config.verbose) {printf("Writing equations...");fflush(stdout);}
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
	solve_Xh(op);
	if(config.verbose) printf("Done\n");

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
	
	h=1;
	dmax=config.newton_dmax;
	
	matrix dphi,dphiex,dp,dT,dXh,dpc,dTc;
	dphi=op->get_var("Phi").block(0,nr-1,0,-1);
	dphiex=op->get_var("Phi").block(nr,nr+nex-1,0,-1);
	err=max(abs(dphi/phi));
	//printf("err(phi)=%e\n",err);
	dp=op->get_var("p");
	err2=max(abs(dp/p));err=err2>err?err2:err;
	while(exist(abs(h*dp/p)>dmax)) h/=2;
	//printf("err(p)=%e\n",err2);
	dT=op->get_var("T");
	err2=max(abs(dT/T));err=err2>err?err2:err;
	while(exist(abs(h*dT/T)>dmax)) h/=2;
	//printf("err(T)=%e\n",err2);
// Compute dXh
	dXh=op->get_var("lnXh");
	err2=max(abs(dXh));err=err2>err?err2:err;
	while(exist(abs(h*dXh)>dmax)) h/=2;
	//printf("err(Xh)=%e\n",err2);
// End compute dXh
	dpc=op->get_var("log_pc");
	err2=fabs(dpc(0));err=err2>err?err2:err;
	while(fabs(h*dpc(0))>dmax) h/=2;
	//printf("err(pc)=%e\n",err2);
	dTc=op->get_var("log_Tc");
	err2=fabs(dTc(0));err=err2>err?err2:err;
	while(fabs(h*dTc(0))>dmax) h/=2;
	//printf("err(Tc)=%e\n",err2);

	phi+=h*dphi;
	phiex+=h*dphiex;
	p+=h*dp;
	T+=h*dT;
// Evolution of Xh, dXh is the variation on ln(Xh) assumed to be small
        Xh+=h*dXh*Xh;

	pc*=exp(h*dpc(0));
	Tc*=exp(h*dTc(0));
	Omega=Omega+h*op->get_var("Omega")(0);
	w+=h*op->get_var("w");
	G+=h*op->get_var("G");
	
	matrix dRi;
	dRi=op->get_var("Ri");
	update_map(h*dRi);
	err2=max(abs(dRi));err=err2>err?err2:err;
	
	rho0=rho;
	
    fill();
	
	err2=max(abs(rho-rho0));err=err2>err?err2:err;
	
	return err;

}

void star2d::update_map(matrix dR) {
    DEBUG_FUNCNAME;
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
    DEBUG_FUNCNAME;

// EOS written rho(P,T). eos.chi_rho is from OPAL,
// and eos.d=chi_t/chi_rho
	op->add_d("rho","log_p",rho/eos.chi_rho);
	op->add_d("rho","log_T",-rho*eos.d);
	op->add_d("rho","log_pc",rho/eos.chi_rho);
	op->add_d("rho","log_Tc",-rho*eos.d);
	op->add_d("rho","log_rhoc",-rho);
	
// Constitutive relation for opacity (xi); formal dependence xi=xi(rho,T)
	op->add_d("opa.xi","rho",opa.dlnxi_lnrho*opa.xi/rho);
	op->add_d("opa.xi","log_rhoc",opa.dlnxi_lnrho*opa.xi);
	op->add_d("opa.xi","T",opa.dlnxi_lnT*opa.xi/T);
	op->add_d("opa.xi","log_Tc",opa.dlnxi_lnT*opa.xi);
	
// Constitutive relation for thermal radiative conductivity k=16sigma*T^3/rho/xi
	op->add_d("opa.k","log_T",3*opa.k);
	op->add_d("opa.k","log_Tc",3*opa.k);
	op->add_d("opa.k","rho",-opa.k/rho);
	op->add_d("opa.k","log_rhoc",-opa.k);
	op->add_d("opa.k","opa.xi",-opa.k/opa.xi);
	
// Constitutive relation for nuclear heat generation; formal dependence nuc.eps=nuc.eps(rho,T)
	op->add_d("nuc.eps","rho",nuc.dlneps_lnrho*nuc.eps/rho);
	op->add_d("nuc.eps","log_rhoc",nuc.dlneps_lnrho*nuc.eps);
	op->add_d("nuc.eps","log_T",nuc.dlneps_lnT*nuc.eps);
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
	op->add_d("s","log_T",eos.cp);
	op->add_d("s","log_Tc",eos.cp);
	op->add_d("s","log_p",-eos.cp*eos.del_ad);
	op->add_d("s","log_pc",-eos.cp*eos.del_ad);
	
}


/// \brief Writes Poisson equation and interface conditions into the solver.
void star2d::solve_poisson(solver *op) {
    DEBUG_FUNCNAME;
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
	
	lap_phi.add(op,"Phi","Phi"); // The equation is named Phi and depends on variable tagged "Phi"
	lap_phi.add(op,"Phi","r");   // The equation is named Phi and depends on variable tagged "r"

	rhs1=-lap_phi.eval()+pi_c*rho; // Expression of the RHS inside the star
	
	//rho
	op->add_d("Phi","rho",-pi_c*ones(nr,nth));  // Coefficient of delta rho

	//pi_c
	op->add_d("Phi","pi_c",-rho);  // Coefficient of delta pi_c

	// phiex
	S.set_value("Phi",phiex);     // Outside the star we define phiex but this is still equation "Phi"
	S.set_map(map.ex);
	lap_phi.add_ex(op,ndomains,"Phi","Phi"); // add exterior domain with variable tagged "Phi"
	lap_phi.add_ex(op,ndomains,"Phi","r"); // add exterior domain with variable tagged "r"
	
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
// 1 designate the "Right condition" and 2 designate the "Left condition"

	for(n=0;n<ndomains+1;n++) {
// set the interface conditions in the Jacobian matrix (continuity of Phi and (1/rz)*dphi/dzeta)
// First start with the bottom conditions
		if(n==0) {
// In the first domain we demand that dphi/dzeta vanishes at the center (left BC)
			op->bc_bot2_add_l(n,"Phi","Phi",ones(1,nth),D.block(0).row(0));
			rhs.setrow(0,-(D,phi).row(0));
		} else {
			if(n<ndomains)
// for all stellar domains, at the bottom (1/rz)d(delta phi)/dzeta is continuuous (left block))
                op->bc_bot2_add_l(n,"Phi","Phi",1/rz.row(j0),D.block(n).row(0)); // IC
			else
// in the last stellar domain continuity demanded with the "ex" fields (left block)
                op->bc_bot2_add_l(n,"Phi","Phi",1/map.ex.rz.row(0),Dex.row(0));

// for all stellar domains, at the bottom (1/rz)d(delta phi)/dzeta is continuuous (right block of IC)
			op->bc_bot1_add_l(n,"Phi","Phi",-1/rz.row(j0-1),D.block(n-1).row(-1));
			
			if(n<ndomains) 
// Functional derivative of 1/rz times dphi/dzeta, complements interface cond. IC (left block)
				op->bc_bot2_add_d(n,"Phi","rz",-1/rz.row(j0)/rz.row(j0)*(D,phi).row(j0));
			else
// Same as above but for the last stellar domain
				op->bc_bot2_add_d(n,"Phi","rz",-1/map.ex.rz.row(0)/map.ex.rz.row(0)*(Dex,phiex).row(0));

// For all stellar domain Right block of the complement of IC
			op->bc_bot1_add_d(n,"Phi","rz",1/rz.row(j0-1)/rz.row(j0-1)*(D,phi).row(j0-1));
			
// set the interface conditions in the RHS
			if(n<ndomains)
// RHS for IC
                rhs.setrow(j0,-(D,phi).row(j0)/rz.row(j0)+(D,phi).row(j0-1)/rz.row(j0-1));
			else
// same as above for the last stellar domain
                rhs.setrow(j0,-(Dex,phiex).row(0)/map.ex.rz.row(0)+(D,phi).row(j0-1)/rz.row(j0-1));
		}
		
// Second take care of the top conditions where we impose continuity of phi and the right BC on phi
		op->bc_top1_add_d(n,"Phi","Phi",ones(1,nth)); // right condition for Jacob

		if(n<ndomains) op->bc_top2_add_d(n,"Phi","Phi",-ones(1,nth)); // left condition for Jacob

// Prepare the RHS for the continuity of phi 
		if(n<ndomains) rhs.setrow(j0+map.gl.npts[n]-1,-phi.row(j0+map.gl.npts[n]-1)); // left terms
		else rhs.setrow(nr+nex-1,-phiex.row(nex-1)); // left terms of the last domain
		if(n<ndomains-1) rhs.setrow(j0+map.gl.npts[n]-1,rhs.row(j0+map.gl.npts[n]-1)
								+phi.row(j0+map.gl.npts[n])); // add the right term
		else if(n==ndomains-1) rhs.setrow(j0+map.gl.npts[n]-1,rhs.row(j0+map.gl.npts[n]-1) 
								+phiex.row(0)); // add the right term for the last
                                                                     // stellar domain but here the right term is phiex
		if(n<ndomains) j0+=map.gl.npts[n];
// Note the Right BC for phi is phi(infty)=0 so no right term in the rhs of the n==ndomain case
	}
	op->set_rhs("Phi",rhs);
}


/// \brief Writes movement and vorticity equations into the solver
void star2d::solve_mov(solver *op) {
    DEBUG_FUNCNAME;
	static bool eqinit=false;
	static symbolic S;
	static sym eq_vort,eq_phi,bc,ic_w;
	static sym_vec eqmov;
	
	if(!eqinit) {
		eqinit=true;
		sym p,G,w,rho,phi;
		sym mu;
		
		p=S.regvar("p");  // Pressure
		G=S.regvar("G");  // Stream function of Merid. Circ.
		w=S.regvar("w");  // Differential Rotation
		rho=S.regvar("rho"); // Density
		phi=S.regvar("Phi"); // Gravitational Potential
		
		#ifdef KINEMATIC_VISC
			mu=rho;
		#else
			mu=S.one;
		#endif
		
		sym_vec V,phivec(COVARIANT),svec;
		sym s;
// Phivec is the E_phi the covariant basis phi-vector
		phivec(0)=0*S.one;phivec(1)=0*S.one;phivec(2)=S.r*sin(S.theta);
		s=S.r*sin(S.theta);
		svec=grad(s);
		V=curl(G*phivec)/rho;
		
// The momentum equation  is eqmov=0; it gives the pressure
		eqmov=grad(p)+rho*grad(phi)-rho*s*w*w*svec;
		eqmov=covariant(eqmov);
		
// The vorticity equation
		eq_vort=(phivec,curl(eqmov/rho));
		
// The angular momentum equation
		eq_phi=div(rho*s*s*w*V)-div(mu*s*s*grad(w));
		
		sym_vec nvec(COVARIANT),tvec(CONTRAVARIANT);
		
		nvec(0)=Dz(S.r);nvec(1)=0*S.one;nvec(2)=0*S.one;
		tvec(0)=0*S.one;tvec(1)=1/S.r;tvec(2)=0*S.one;
		
// the special boundary condition that couples G and w for the determination of w
		bc=mu*s*s*(nvec,grad(w))+G*(tvec,grad(s*s*w));
		
		ic_w=covariant(eqmov-grad(p))(1);
	}
	
	S.set_value("p",p);
	S.set_value("G",G,11);
	S.set_value("w",w);
	S.set_value("rho",rho);
	S.set_value("Phi",phi);
	S.set_map(map);
	
	op->add_d("p","log_p",p);
	
// First component of the equation of motion
// we explicit the dependences
	eqmov(0).add(op,"log_p","p");
	eqmov(0).add(op,"log_p","w");
	eqmov(0).add(op,"log_p","rho");
	eqmov(0).add(op,"log_p","Phi");
	eqmov(0).add(op,"log_p","r");
	op->set_rhs("log_p",-eqmov(0).eval());
	
// Equation of vorticity, dependences...
	eq_vort.add(op,"w","p");
	eq_vort.add(op,"w","w");
	eq_vort.add(op,"w","rho");
	eq_vort.add(op,"w","r");
	op->set_rhs("w",-eq_vort.eval());
	
// Equation of angular momentum, dependences...
	eq_phi.add(op,"G","w");
	eq_phi.add(op,"G","G");
	eq_phi.add(op,"G","rho");
	eq_phi.add(op,"G","r");
	op->set_rhs("G",-eq_phi.eval());
	
// Boundary conditions

	matrix rhs;
	matrix q,TT;
	int j0;
	
	// log_p - Pressure
	rhs=op->get_rhs("log_p");
	op->bc_bot2_add_d(0,"log_p","p",ones(1,nth));
	rhs.setrow(0,-p.row(0)+1); // central pressure is unity
	
	j0=map.gl.npts[0];
	for(int n=1;n<ndomains;n++) { // continuity of pressure
		op->bc_bot2_add_d(n,"log_p","p",ones(1,nth));
		op->bc_bot1_add_d(n,"log_p","p",-ones(1,nth));
		rhs.setrow(j0,-p.row(j0)+p.row(j0-1));	
		j0+=map.gl.npts[n];
	}
	op->set_rhs("log_p",rhs);
	
	// w - Differential Rotation
	rhs=op->get_rhs("w");
	op->bc_bot2_add_l(0,"w","w",ones(1,nth),D.block(0).row(0));
	rhs.setrow(0,-(D,w).row(0));
	j0=0;
	for(int n=0;n<ndomains-1;n++) {
		ic_w.bc_top1_add(op,n,"w","w");
		ic_w.bc_top1_add(op,n,"w","rho");
		ic_w.bc_top1_add(op,n,"w","Phi");
		ic_w.bc_top1_add(op,n,"w","r");
		ic_w.bc_top2_add(op,n,"w","w",-ones(1,nth));
		ic_w.bc_top2_add(op,n,"w","rho",-ones(1,nth));
		ic_w.bc_top2_add(op,n,"w","Phi",-ones(1,nth));
		ic_w.bc_top2_add(op,n,"w","r",-ones(1,nth));
		rhs.setrow(j0+map.gl.npts[n]-1,
			-ic_w.eval().row(j0+map.gl.npts[n]-1)+ic_w.eval().row(j0+map.gl.npts[n]));
		
		j0+=map.gl.npts[n];
	}
	
// Surface condition on the angular velocity and stream function
	q=ones(1,nth);
	q(0)=0;
	map.leg.eval_00(th,PI/2*ones(1,nth),TT);
	bc.bc_top1_add(op,ndomains-1,"w","G",q);
	bc.bc_top1_add(op,ndomains-1,"w","w",q);
	bc.bc_top1_add(op,ndomains-1,"w","rho",q);
	bc.bc_top1_add(op,ndomains-1,"w","r",q);
	op->bc_top1_add_r(ndomains-1,"w","w",(1-q),TT);
	op->bc_top1_add_d(ndomains-1,"w","Omega",-(1-q));
	rhs.setrow(-1,-q*bc.eval().row(-1)-(1-q)*((w.row(-1),TT)-Omega));
	op->set_rhs("w",rhs);
	
	//G - Meridional circulation
	rhs=op->get_rhs("G");
	op->bc_bot2_add_d(0,"G","G",ones(1,nth));
	rhs.setrow(0,-G.row(0));
	j0=map.gl.npts[0];
	for(int n=1;n<ndomains;n++) {
		op->bc_bot2_add_d(n,"G","G",ones(1,nth));
		op->bc_bot1_add_d(n,"G","G",-ones(1,nth));
		rhs.setrow(j0,-G.row(j0)+G.row(j0-1));
	
		j0+=map.gl.npts[n];
	}
	op->set_rhs("G",rhs);
	
	//pi_c -  non-dimensional parameter
	rhs=zeros(ndomains,1);
	j0=0;
	for(int n=0;n<ndomains;n++) {
		if(n<ndomains-1) {
			op->bc_top1_add_d(n,"pi_c","pi_c",ones(1,1));
			op->bc_top2_add_d(n,"pi_c","pi_c",-ones(1,1));
		} else {
			map.leg.eval_00(th,0,TT);
			op->bc_top1_add_r(n,"pi_c","ps",-ones(1,1),TT);
			op->bc_top1_add_r(n,"pi_c","p",ones(1,1),TT);
			rhs(ndomains-1)=(ps-p.row(-1),TT)(0);
		}
		
		j0+=map.gl.npts[n];
	}
	op->set_rhs("pi_c",rhs);
	
	if(Omega==0) {
		op->reset("w");
		op->add_d("w","w",ones(nr,nth));
		op->add_d("w","Omega",-ones(nr,nth));
		op->set_rhs("w",zeros(nr,nth));
		op->reset("G");
		op->add_d("G","G",ones(nr,nth));
		op->set_rhs("G",zeros(nr,nth));
	}
}

/// \brief Writes Xh equation for the evolution
//Evolution Xh --------------------------------------
void star2d::solve_Xh(solver *op) {
    DEBUG_FUNCNAME;

	printf("Start of solve_Xh\n");
    double Qmc2=(4*HYDROGEN_MASS-AMASS["He4"]*UMA)*C_LIGHT*C_LIGHT;
    double factor=4*HYDROGEN_MASS/Qmc2*MYR*dt;
        op->add_d("lnXh","lnXh",ones(nr,nth));
        op->add_d("lnXh","log_T",factor*nuc.eps/Xh*nuc.dlneps_lnT);
        op->add_d("lnXh","rho",factor*nuc.eps/Xh*nuc.dlneps_lnrho/rho);

    matrix rhs=log(Xh_prec)-log(Xh)-factor*nuc.eps/Xh;
    op->set_rhs("lnXh",rhs);
	printf("End of solve_Xh\n");

}
//Evolution Xh end-----------------------------------


/// \brief Writes temperature and luminosity equations and interface conditions
/// into the solver.
void star2d::solve_temp(solver *op) {
    DEBUG_FUNCNAME;
	int n,j0;
	matrix q;
	char eqn[8];
	matrix &gzz=map.gzz, &gzt=map.gzt, &rz=map.rz;
	
	op->add_d("T","log_T",T);
	strcpy(eqn,"log_T");
	
	//Luminosity

	matrix rhs_lum,lum;
	
	lum=zeros(ndomains,1);
	j0=0;
// for each domain we compute the luminosity at the upper boundary
// lum(n) =int_0^pi\int_0^eta_n 2*pi*r^2*rz*rho*eps dzeta sin(theta)dtheta
	for(n=0;n<ndomains;n++) {
		if(n) lum(n)=lum(n-1);
		lum(n)+=2*PI*Lambda*(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),
			(rho*nuc.eps*r*r*rz).block(j0,j0+map.gl.npts[n]-1,0,-1),map.leg.I_00)(0);
		j0+=map.gl.npts[n];
	}

// Now we code the equation of luminosity, namely Lum=intvol Lambda*rho*eps*dV
// Since this equation is a one-line matrix, it is implemented as a BC
	rhs_lum=zeros(ndomains,1);
	j0=0;
	for(n=0;n<ndomains;n++) {
		op->bc_bot2_add_d(n,"lum","lum",ones(1,1));
		op->bc_bot2_add_lri(n,"lum","rho",-2*PI*Lambda*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,(r*r*rz*nuc.eps).block(j0,j0+map.gl.npts[n]-1,0,-1));
		op->bc_bot2_add_lri(n,"lum","nuc.eps",-2*PI*Lambda*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,(r*r*rz*rho).block(j0,j0+map.gl.npts[n]-1,0,-1));
		op->bc_bot2_add_d(n,"lum","Lambda",-2*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(rho*nuc.eps*r*r*rz).block(j0,j0+map.gl.npts[n]-1,0,-1),map.leg.I_00));
//metric terms from dV
		op->bc_bot2_add_lri(n,"lum","r",-2*PI*Lambda*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,(2*r*rz*rho*nuc.eps).block(j0,j0+map.gl.npts[n]-1,0,-1));
		op->bc_bot2_add_lri(n,"lum","rz",-2*PI*Lambda*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,(r*r*rho*nuc.eps).block(j0,j0+map.gl.npts[n]-1,0,-1));
			
		if(n) op->bc_bot1_add_d(n,"lum","lum",-ones(1,1));
		j0+=map.gl.npts[n];
	}
	op->set_rhs("lum",rhs_lum);
	
// Frad is (-xi*grad(T) scal E^zeta
// Note that one property of Frad is that on zeta=cte surface
// intsurf Frad dS = Lum ; note that here dS=2*pi*r^2*rz*sin(th)*dth
	
	matrix rhs_Frad,Frad;
	int j1;
	
	Frad=-opa.xi*(gzz*(D,T)+gzt*(T,Dt)); // explicit expression of Frad
	rhs_Frad=zeros(ndomains*2-1,nth);
	j0=0;
	for(n=0;n<ndomains;n++) {
		j1=j0+map.gl.npts[n]-1;
		
		if(n) op->bc_bot2_add_d(n,"Frad","Frad",ones(1,nth));
		op->bc_top1_add_d(n,"Frad","Frad",ones(1,nth));

// terms from temperature variations (delta T)		
		q=opa.xi*gzz;
		if(n) op->bc_bot2_add_l(n,"Frad","T",q.row(j0),D.block(n).row(0));
		op->bc_top1_add_l(n,"Frad","T",q.row(j1),D.block(n).row(-1));
		q=opa.xi*gzt;
		if(n) op->bc_bot2_add_r(n,"Frad","T",q.row(j0),Dt);
		op->bc_top1_add_r(n,"Frad","T",q.row(j1),Dt);

// terms from delta xi				
		if(n) op->bc_bot2_add_d(n,"Frad","opa.xi",-Frad.row(j0)/opa.xi.row(j0));
		op->bc_top1_add_d(n,"Frad","opa.xi",-Frad.row(j1)/opa.xi.row(j1));

// terms from delta r
		q=opa.xi*(-2.*r*gzt*gzt*(D,T)-2./r*gzt*(T,Dt));
		if(n) op->bc_bot2_add_d(n,"Frad","r",q.row(j0));
		op->bc_top1_add_d(n,"Frad","r",q.row(j1));
//terms from delta rz
		q=opa.xi*(-2./rz*gzz*(D,T)-1./rz*gzt*(T,Dt));
		if(n) op->bc_bot2_add_d(n,"Frad","rz",q.row(j0));
		op->bc_top1_add_d(n,"Frad","rz",q.row(j1));
// terms from d(delta r)/dtheta, which are resumed to delta r again
// thanks to Dt multiplication
		q=opa.xi*(-2./rz*gzt*(D,T)-1./r/r/rz*(T,Dt));
		if(n) op->bc_bot2_add_r(n,"Frad","r",q.row(j0),Dt);
		op->bc_top1_add_r(n,"Frad","r",q.row(j1),Dt);
		
		j0=j1+1;
	}
	op->set_rhs("Frad",rhs_Frad);
		
	
// Temperature field
	
	matrix rhs_T,rhs_Lambda;
	matrix TT,qcore,qenv;
	
	qenv=zeros(nr,nth);
	qcore=qenv;
	j0=0;
// define the grid-points belonging to the core
	for(n=0;n<ndomains;n++) {
		if(n<conv) qcore.setblock(j0,j0+map.gl.npts[n]-1,0,-1,ones(map.gl.npts[n],nth));
		else qenv.setblock(j0,j0+map.gl.npts[n]-1,0,-1,ones(map.gl.npts[n],nth));
		j0+=map.gl.npts[n];
	}
	
	
	rhs_T=zeros(nr,nth);

	symbolic S;
	sym T_,xi_,s_,G_;
	sym div_Frad;
	
	S.set_map(map);
	
	T_=S.regvar("T");
	xi_=S.regvar("opa.xi");
	s_=S.regvar("s");
	G_=S.regvar("G");
	S.set_value("T",T);
	S.set_value("opa.xi",opa.xi);
	S.set_value("s",entropy());
	S.set_value("G",G);
	
// Diffusion terms of temperature equation
// Recall eqn is log_T
	div_Frad=-div(-xi_*grad(T_))/xi_;
	
	div_Frad.add(op,eqn,"T",qenv);
	div_Frad.add(op,eqn,"opa.xi",qenv);
	div_Frad.add(op,eqn,"r",qenv);
	rhs_T-=div_Frad.eval()*qenv;

// Explicit the expression of the functional derivative of the 
// temperature equation	
	op->add_d(eqn,"nuc.eps",qenv*Lambda*rho/opa.xi);
	op->add_d(eqn,"rho",qenv*Lambda*nuc.eps/opa.xi);	
	op->add_d(eqn,"Lambda",qenv*rho*nuc.eps/opa.xi);
	op->add_d(eqn,"opa.xi",-qenv*Lambda*rho*nuc.eps/opa.xi/opa.xi);
	rhs_T+=-qenv*Lambda*rho*nuc.eps/opa.xi;
	
	// Advection
	/*
	sym adv;
	sym_vec rhov_,phivec(COVARIANT);
	
	phivec(0)=0*S.one;phivec(1)=0*S.one;phivec(2)=S.r*S.sint;
	rhov_=curl(G_*phivec);
	adv=-T_*(rhov_,grad(s_))/xi_;
	
	adv.add(op,eqn,"T",qenv*Ekman*sqrt(rhoc*pc)*R);
	adv.add(op,eqn,"s",qenv*Ekman*sqrt(rhoc*pc)*R);
	adv.add(op,eqn,"G",qenv*Ekman*sqrt(rhoc*pc)*R);
	adv.add(op,eqn,"opa.xi",qenv*Ekman*sqrt(rhoc*pc)*R);
	
	matrix advec;
	
	advec=qenv*Ekman*sqrt(rhoc*pc)*R*adv.eval();
	
	op->add_d(eqn,"log_R",advec);
	op->add_d(eqn,"log_pc",0.5*advec);
	op->add_d(eqn,"log_rhoc",0.5*advec);
	rhs_T+=-advec;
	*/

//Core convection, equation ds/dzeta=0
	op->add_l(eqn,"s",qcore,D);
	//rhs_T+=-qcore*(D,eos.s);
	rhs_T+=-qcore*(D,entropy());
	
	rhs_Lambda=zeros(ndomains,1);
	
	map.leg.eval_00(th,0,TT);
	
// Interface and boundary conditions for the temperature
	j0=0;
	for(n=0;n<ndomains;n++) {
		if(n==0) { // In the first domain T(0)=1
			op->bc_bot2_add_d(n,eqn,"T",ones(1,nth));
			rhs_T.setrow(j0,1-T.row(j0));
		} else {  // we impose the continuity of T
			op->bc_bot2_add_d(n,eqn,"T",ones(1,nth));
			op->bc_bot1_add_d(n,eqn,"T",-ones(1,nth));
			rhs_T.setrow(j0,-T.row(j0)+T.row(j0-1));
		}
        // Radiative envelope: the continuity of (1/rz)(dT/dzeta) is imposed
		if(n>=conv) {
			if(n<ndomains-1) {
				/*op->bc_top1_add_d(n,eqn,"Frad",rz.row(j0+map.gl.npts[n]-1));
				op->bc_top2_add_d(n,eqn,"Frad",-rz.row(j0+map.gl.npts[n]-1));
				op->bc_top1_add_d(n,eqn,"rz",Frad.row(j0+map.gl.npts[n]-1));
				op->bc_top2_add_d(n,eqn,"rz",-Frad.row(j0+map.gl.npts[n]-1));
				
				rhs_T.setrow(j0+map.gl.npts[n]-1,
					-Frad.row(j0+map.gl.npts[n]-1)*rz.row(j0+map.gl.npts[n]-1)
					+Frad.row(j0+map.gl.npts[n])*rz.row(j0+map.gl.npts[n]));*/
				op->bc_top1_add_l(n,eqn,"T",1/rz.row(j0+map.gl.npts[n]-1),D.block(n).row(-1));
				op->bc_top1_add_d(n,eqn,"rz",-((D,T)/rz/rz).row(j0+map.gl.npts[n]-1));
				op->bc_top2_add_l(n,eqn,"T",-1/rz.row(j0+map.gl.npts[n]),D.block(n+1).row(0));
				op->bc_top2_add_d(n,eqn,"rz",((D,T)/rz/rz).row(j0+map.gl.npts[n]));
				rhs_T.setrow(j0+map.gl.npts[n]-1,
					((D,T)/rz).row(j0+map.gl.npts[n])-((D,T)/rz).row(j0+map.gl.npts[n]-1));
			} else { // In the last domain set the upper BC T=Ts
				op->bc_top1_add_d(n,eqn,"T",ones(1,nth));
				op->bc_top1_add_d(n,eqn,"Ts",-ones(1,nth));
				rhs_T.setrow(-1,Ts-T.row(-1));
			}
		}
		
		if(n<conv) {
        // Inside the convective core Lambda is continuous...
			op->bc_top1_add_d(n,"Lambda","Lambda",ones(1,1));
			op->bc_top2_add_d(n,"Lambda","Lambda",-ones(1,1));
		} else if(n==conv) {
			if(n==0) { // There is a central core
				map.leg.eval_00(th,PI/2,q);
				op->bc_bot2_add_lr(n,"Lambda","T",ones(1,1),D.block(0).row(0),q);
				rhs_Lambda(0)=-((D,T).row(0),q)(0);
			} else { // The upper domain in the CC is not the central domain
                                 // The Lambda eqn says that the total radiative flux is the luminosity at the boundary
				op->bc_bot2_add_ri(n,"Lambda","Frad",2*PI*ones(1,1),map.leg.I_00,(r*r*rz).row(j0));
				op->bc_bot2_add_ri(n,"Lambda","r",2*PI*ones(1,1),map.leg.I_00,(Frad*2*r*rz).row(j0));
				op->bc_bot2_add_ri(n,"Lambda","rz",2*PI*ones(1,1),map.leg.I_00,(Frad*r*r).row(j0));
				op->bc_bot1_add_d(n,"Lambda","lum",-ones(1,1));
				rhs_Lambda(n)=-2*PI*(Frad.row(j0)*(r*r*rz).row(j0),map.leg.I_00)(0)+lum(n-1);
			}
		} else {
        // Outside the convective core Lambda is continuous...
			op->bc_bot2_add_d(n,"Lambda","Lambda",ones(1,1));
			op->bc_bot1_add_d(n,"Lambda","Lambda",-ones(1,1));
		}
		
		j0+=map.gl.npts[n];
	}
	op->set_rhs(eqn,rhs_T);
	op->set_rhs("Lambda",rhs_Lambda);
	
}


/// \brief Writes the equations for the dimensional quantities (T_c, rho_c, R, etc.)

void star2d::solve_dim(solver *op) {
    DEBUG_FUNCNAME;
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
    DEBUG_FUNCNAME;
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
			op->bc_bot2_add_lr(n,"Ri",LOG_PRES,q,delta,TT);
			delta=zeros(1,map.gl.npts[n-1]);delta(0)=1;delta(-1)=-1;
			op->bc_bot1_add_lr(n,"Ri",LOG_PRES,-q,delta,TT);
			rhs.setrow(n,(log(PRES.row(j0+map.gl.npts[n]-1))-log(PRES.row(j0))
			-log(PRES.row(j0-1))+log(PRES.row(j0-map.gl.npts[n-1])),TT));
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
    DEBUG_FUNCNAME;
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
    DEBUG_FUNCNAME;
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
    DEBUG_FUNCNAME;
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
/// To be checked
void star2d::solve_atm(solver *op) {

	if(!strcmp(atm.name,"simple")) {
		matrix q;
		int n=ndomains-1;
	
		op->bc_top1_add_d(n,"ps","ps",1/ps);
		op->bc_top1_add_d(n,"ps","gsup",-1/gsup());
		op->bc_top1_add_d(n,"ps","opa.k",1/opa.k.row(-1));
		op->bc_top1_add_d(n,"ps","log_pc",ones(1,nth));
		op->set_rhs("ps",zeros(1,nth));
	
		op->bc_top1_add_d(n,"Ts","Ts",1/Ts);
		op->bc_top1_add_d(n,"Ts","p",-0.25/p.row(-1));
		op->bc_top1_add_d(n,"Ts","ps",0.25/ps);
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


/// \brief Routine to check the Jacobian matrix
void star2d::check_jacobian(solver *op,const char *eqn) {
    DEBUG_FUNCNAME;
	star2d B;
	matrix rhs,drhs,drhs2,qq;
	matrix *y;
	double q;
	int i,j,j0;
	
	y=new matrix[op->get_nvar()];
	B=*this;
	// Perturbar el modelo
	{
		double a,ar,asc;
		
		a=1e-8;ar=1e-8;
		asc=a>ar?a:ar;
		B.phi=B.phi+a*B.phi+ar*B.phi*random_matrix(nr,nth);
		B.phiex=B.phiex+a*B.phiex+ar*B.phiex*random_matrix(nex,nth);
		B.p=B.p+a*B.p+ar*B.p*random_matrix(nr,nth);
		B.pc=B.pc+asc*B.pc;
		B.T=B.T+a*B.T+ar*B.T*random_matrix(nr,nth);
		B.Tc=B.Tc+asc*B.Tc;
		B.w=B.w+a*B.w+ar*B.w*random_matrix(nr,nth);
		B.Omega=B.Omega+asc*B.Omega;
		B.G=B.G+a*B.G+ar*B.G*random_matrix(nr,nth);
		B.map.R=B.map.R+a*B.map.R*sin(th)*sin(th)+ar*B.map.R*random_matrix(ndomains,nth);
		B.map.remap();
	}
	
	B.fill();
	
	i=op->get_id("rho");
	y[i]=zeros(nr,nth);
	i=op->get_id("opa.xi");
	y[i]=zeros(nr,nth);
	i=op->get_id("nuc.eps");
	y[i]=zeros(nr,nth);
	i=op->get_id("r");
	y[i]=zeros(nr+nex,nth);
	i=op->get_id("rz");
	y[i]=zeros(nr+nex,nth);
	i=op->get_id("s");
	y[i]=zeros(nr,nth);
	i=op->get_id("opa.k");
	y[i]=zeros(nr,nth);
	
	
	i=op->get_id("Phi");
	y[i]=zeros(nr+nex,nth);
	y[i].setblock(0,nr-1,0,-1,B.phi-phi);
	y[i].setblock(nr,nr+nex-1,0,-1,B.phiex-phiex);
	i=op->get_id("p");
	y[i]=B.p-p;
	i=op->get_id("log_p");
	y[i]=log(B.p)-log(p);
	i=op->get_id("pi_c");
	y[i]=(B.pi_c-pi_c)*ones(ndomains,1);
	i=op->get_id("T");
	y[i]=B.T-T;
	i=op->get_id("log_T");
	y[i]=log(B.T)-log(T);
	i=op->get_id("Lambda");
	y[i]=(B.Lambda-Lambda)*ones(ndomains,1);
	i=op->get_id("eta");
	y[i]=B.map.eta-map.eta;
	j=i;
	i=op->get_id("deta");
	y[i]=y[j].block(1,ndomains,0,0)-y[j].block(0,ndomains-1,0,0);
	i=op->get_id("Ri");
	y[i]=B.map.R-map.R;
	j=i;
	i=op->get_id("dRi");
	y[i]=y[j].block(1,ndomains,0,-1)-y[j].block(0,ndomains-1,0,-1);
	i=op->get_id("Omega");
	y[i]=(B.Omega-Omega)*ones(ndomains+1,1);
	i=op->get_id("log_rhoc");
	y[i]=(log(B.rhoc)-log(rhoc))*ones(ndomains,1);
	i=op->get_id("log_pc");
	y[i]=(log(B.pc)-log(pc))*ones(ndomains,1);
	i=op->get_id("log_Tc");
	y[i]=(log(B.Tc)-log(Tc))*ones(ndomains,1);
	i=op->get_id("log_R");
	y[i]=(log(B.R)-log(R))*ones(ndomains,1);
	i=op->get_id("m");
	q=0;
	j0=0;
	y[i]=zeros(ndomains,1);
	for(j=0;j<ndomains;j++) {
		q+=2*PI*(B.map.gl.I.block(0,0,j0,j0+map.gl.npts[j]-1),
			(B.rho*B.r*B.r*B.map.rz).block(j0,j0+map.gl.npts[j]-1,0,-1),
			B.map.leg.I_00)(0)-
			2*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[j]-1),
			(rho*r*r*map.rz).block(j0,j0+map.gl.npts[j]-1,0,-1),
			map.leg.I_00)(0);
		y[i](j)=q;
		j0+=map.gl.npts[j];
	}
	i=op->get_id("ps");
	y[i]=B.ps-ps;
	i=op->get_id("Ts");
	y[i]=B.Ts-Ts;
	i=op->get_id("lum");
	y[i]=zeros(ndomains,1);
	j0=0;
	q=0;
	for(j=0;j<ndomains;j++) {
		q+=2*PI*B.Lambda*(B.map.gl.I.block(0,0,j0,j0+map.gl.npts[j]-1),
			(B.rho*B.nuc.eps*B.r*B.r*B.map.rz).block(j0,j0+map.gl.npts[j]-1,0,-1),B.map.leg.I_00)(0)-
			2*PI*Lambda*(map.gl.I.block(0,0,j0,j0+map.gl.npts[j]-1),
			(rho*nuc.eps*r*r*map.rz).block(j0,j0+map.gl.npts[j]-1,0,-1),map.leg.I_00)(0);
		y[i](j)=q;
		j0+=map.gl.npts[j];
	}
	i=op->get_id("Frad");
	y[i]=zeros(ndomains*2-1,nth);
	j0=0;
	matrix Frad,BFrad;
	Frad=-opa.xi*(map.gzz*(D,T)+map.gzt*(T,Dt));
	BFrad=-B.opa.xi*(B.map.gzz*(B.D,B.T)+B.map.gzt*(B.T,B.Dt));
	for(j=0;j<ndomains;j++) {	
		if(j) y[i].setrow(2*j-1,BFrad.row(j0)-Frad.row(j0));
		y[i].setrow(2*j,BFrad.row(j0+map.gl.npts[j]-1)-Frad.row(j0+map.gl.npts[j]-1));
		j0+=map.gl.npts[j];
	}
	i=op->get_id("gsup");
	y[i]=B.gsup()-gsup();
	i=op->get_id("Teff");
	y[i]=B.Teff()-Teff();
	i=op->get_id("w");
	y[i]=B.w-w;
	i=op->get_id("G");
	y[i]=B.G-G;	
	i=op->get_id("gamma");
	y[i]=sqrt(1-B.map.rt*B.map.rt/B.r/B.r)-sqrt(1-map.rt*map.rt/r/r);	
	
	/*
	matrix dlnp,dlnT,dlnrho,dlnrho2,dq,dq2;
	double dlnpc,dlnTc,dlnrhoc;
	dlnp=log(B.p)-log(p);dlnT=log(B.T)-log(T);dlnrho=log(B.rho)-log(rho);
	dlnpc=log(B.pc)-log(pc);dlnTc=log(B.Tc)-log(Tc);dlnrhoc=log(B.rhoc)-log(rhoc);
	dlnrho2=dlnp/eos.chi_rho-dlnT*eos.d+dlnpc/eos.chi_rho-dlnTc*eos.d-dlnrhoc;
	dq=1./B.rho/B.rho*(B.rho,B.Dt)-1./rho/rho*(rho,Dt);
	dq2=-1./rho/rho*(rho,Dt)*dlnrho+1./rho*(dlnrho,Dt);
	dq2=-2./rho/rho*(rho,Dt)*dlnrho+1./rho/rho*(rho*dlnrho,Dt);
	*/
	
	B.solve(op);
	rhs=op->get_rhs(eqn);
	B=*this;
	B.solve(op);
	
	op->mult(y);
	drhs=rhs-op->get_rhs(eqn);
	drhs2=y[op->get_id(eqn)]+drhs;
	
	/*drhs=dq;
	drhs2=dq-dq2;
	*/
	static figure fig("/XSERVE");
	
	/*
	fig.subplot(1,2);
	
	fig.colorbar();
	B.draw(&fig,log10(abs(drhs)+1e-15));
	fig.colorbar();
	B.draw(&fig,log10(abs(drhs2)+1e-15));
	*/
	drhs.write();
	drhs2.write();
	int nn;
	if (drhs.nrows()==1) nn=drhs.ncols();
	else nn=drhs.nrows();
	
	fig.axis(0-nn/20.,nn*(1+1./20),-15,0);
	fig.plot(log10(abs(drhs)+1e-20),"b");
	fig.hold(1);
	fig.plot(log10(abs(drhs2)+1e-20),"r");
	fig.hold(0);
	
	delete [] y;
}


