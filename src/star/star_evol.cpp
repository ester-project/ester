// To keep compatibility with configure
#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "utils.h"
#include "star.h"


star_evol::star_evol() {
    comp_inited = false;
	Lz_obj=0;
}

star_evol::star_evol(const star2d &A) : star2d(A) {
    comp_inited = false;
	Lz_obj=A.Lz();
}

int star_evol::read(const char *input_file, int dim) {
	int out;
	out=star2d::read(input_file, dim);
	Lz_obj=Lz();
    comp_inited = true;
	return out;
}

void star_evol::fill() {
	star2d::fill();
	Omega_bk=Omega/Omegac;
}

// void star_evol::init_comp() {
//     if (!converged)
//         return;
//     if (comp_inited) {
//         comp["H"] -= T*0.01*comp["H"]; // to be fixed
//         return;
//     }
//     else {
//         comp_inited = true;
//         star2d::init_comp();
//     }
// }

solver * star_evol::init_solver(int nvar_add) {
	return star2d::init_solver(nvar_add+1);
}

void star_evol::register_variables(solver *op) {

	star2d::register_variables(op);
	op->regvar("Lz");

}

void star_evol::solve_Omega(solver *op) {

	matrix rhs;
	int n;

	if(Lz_obj==0) {
		for(n=0;n<ndomains;n++) {
			op->add_d(n,"Lz","Lz",ones(1,1));
			op->add_d(n,"Omega","Omega",ones(1,1));
			op->set_rhs("Lz",zeros(ndomains,1));
			op->set_rhs("Omega",zeros(ndomains,1));
		}
		return;
	}

	rhs=zeros(ndomains,1);
	int j0=0;
	for(n=0;n<ndomains;n++) {
		op->bc_bot2_add_d(n,"Lz","Lz",ones(1,1));
		//rho
		op->bc_bot2_add_lri(n,"Lz","rho",-2*PI*ones(1,1),
			map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,
			(w*r*r*sin(th)*sin(th)*r*r*map.rz).block(j0,j0+map.gl.npts[n]-1,0,-1));
		// w
		op->bc_bot2_add_lri(n,"Lz","w",-2*PI*ones(1,1),
			map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,
			(rho*r*r*sin(th)*sin(th)*r*r*map.rz).block(j0,j0+map.gl.npts[n]-1,0,-1));
		//r (rz)
		op->bc_bot2_add_lri(n,"Lz","r",-2*PI*ones(1,1),
			map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,
			(4*r*r*r*sin(th)*sin(th)*map.rz*rho*w).block(j0,j0+map.gl.npts[n]-1,0,-1));
		op->bc_bot2_add_lri(n,"Lz","rz",-2*PI*ones(1,1),
			map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,
			(r*r*r*r*sin(th)*sin(th)*rho*w).block(j0,j0+map.gl.npts[n]-1,0,-1));
		
		if(n) op->bc_bot1_add_d(n,"Lz","Lz",-ones(1,1));
		j0+=map.gl.npts[n];
	}
	op->set_rhs("Lz",rhs);
	
	rhs=zeros(ndomains,1);
	for(n=0;n<ndomains-1;n++) {
		op->bc_top1_add_d(n,"Omega","Omega",ones(1,1));
		op->bc_top2_add_d(n,"Omega","Omega",-ones(1,1));
	}
	n=ndomains-1;
	double lz=Lz()/R/R/R/R/sqrt(pc*rhoc);
	op->bc_top1_add_d(n,"Omega","Lz",ones(1,1)/lz);
	op->bc_top1_add_d(n,"Omega","log_R",4*ones(1,1));
	op->bc_top1_add_d(n,"Omega","log_pc",0.5*ones(1,1));
	op->bc_top1_add_d(n,"Omega","log_rhoc",0.5*ones(1,1));
	rhs(n)=log(Lz_obj/Lz());
	op->set_rhs("Omega",rhs);
}
