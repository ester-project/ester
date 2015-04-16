#include "ester-config.h"
#include "utils.h"
#include "star.h"
#include "cesam.h"

#include <stdlib.h>

void star_evol::update_comp(star_evol &prev, double dt) {
    int nchim = 10;
    cesam::compo c(this, nchim);
    int npts = 0;
    matrix t0 = prev.T * prev.Tc;
    matrix t1 = T * Tc;
    matrix rho0 = prev.rho * prev.rhoc;
    matrix rho1 = rho * rhoc;

    for (int n=0; n<this->conv; n++) {
        npts += map.gl.npts[n];
    }

    // Central convective zone
    for (int i=0; i<nth; ++i) {
        double *cesam_abon = c.to_cesam(comp.block(0, npts-1, i, i));

        redirect_stdout("cesam.log");
        cesam::update_comp(t0.block(0, npts-1, i, i).data(), t1.block(0, npts-1, i, i).data(),
                rho0.block(0, npts-1, i, i).data(), rho1.block(0, npts-1, i, i).data(),
                cesam_abon,
                r.block(0, npts-1, i, i).data(), nchim, npts, dt);
        restore_stdout();
        matrix_map *newcomp = c.from_cesam(cesam_abon);
        comp.setblock(0, npts-1, i, i, *newcomp);
    }
    // cesam mix elements in the convective zone; but we still have to mix them
    // over theta
    composition_map avg = comp.block(0, npts-1, 0, 0);
    for (int i=1; i<nth; ++i) {
        avg += comp.block(0, npts-1, i, i);
    }
    avg = avg * (1.0/(double)nth);
    for (int i=0; i<nth; ++i) {
        comp.setblock(0, npts-1, i, i, avg);
    }

    // Radiative envelope
    for (int n=npts; n<nr; n++) {
        for (int i=0; i<nth; ++i) {
            double *cesam_abon = c.to_cesam(comp.block(n, n, i, i));

            redirect_stdout("cesam.log");
            cesam::update_comp(t0.block(n, n, i, i).data(), t1.block(n, n, i, i).data(),
                    rho0.block(n, n, i, i).data(), rho1.block(n, n, i, i).data(),
                    cesam_abon,
                    r.block(n, n, i, i).data(), nchim, 1, dt);
            matrix_map *newcomp = c.from_cesam(cesam_abon);
            restore_stdout();
            comp.setblock(n, n, i, i, *newcomp);
        }
    }
}

star_evol::star_evol() {
	Lz_obj=0;
}

star_evol::star_evol(const star2d &A) : star2d(A) {
	Lz_obj=A.Lz();
}

int star_evol::read(const char *input_file, int dim) {
	int out;
	out=star2d::read(input_file, dim);
	Lz_obj=Lz();
	return out;
}

void star_evol::fill() {
    star2d::fill();
	Omega_bk=Omega/Omegac;
}

void star_evol::init_comp() {
    // only set initial compo if composition is not initialized
    // this prevents overriding composition during evolution
    if (comp.size() == 1) {
        star2d::init_comp();
    }
}

solver *star_evol::init_solver(int nvar_add) {
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
