#include"starR.h"
#include<string.h>
#include<stdlib.h>

star1dR::star1dR() {Teff_obj=-1;}

int star1dR::check_arg(char *arg,char *val,int *change_grid) {

	if(!strcmp(arg,"R")) {
		if(val==NULL) return 2;
		R=atof(val)*R_SUN;
		return 0;
	}
	else if(!strcmp(arg,"Teff")) {
		if(val==NULL) return 2;
		Teff_obj=atof(val);
		return 0;
	}
	else if(!strcmp(arg,"M")) {
		return 1;
	}
	else if(!strcmp(arg,"Xc")) {
		return 1;
	}

	return star1d::check_arg(arg,val,change_grid);
}



solver *star1dR::init_solver(int nvar_add) {

	return star1d::init_solver(2+nvar_add);
}

void star1dR::register_variables(solver *op) {

	star1d::register_variables(op);
	op->regvar("Xc");
	op->regvar("log_M");
}

double star1dR::solve(solver *op, matrix_map& error_map, int nit) {

	double err;

	err=star1d::solve(op, error_map, nit);
	Xc+=op->get_var("Xc")(0);
	fill();

	return err;
}

double star1dR::solve(solver *op) {
    matrix_map error_map;
    return solve(op, error_map, 0);
}

void star1dR::fill() {

	eq_state();
	m=4*PI*(map.gl.I,rho*r*r)(0);
	M=m*rhoc*R*R*R;
	double R0=R;
	star1d::fill();
	R=R0;
	if(Teff_obj==-1) Teff_obj=map.leg.eval_00(Teff(),0)(0);
}

void star1dR::solve_dim(solver *op) {

	star1d::solve_dim(op);
	
	matrix rhs;
	
	op->reset(ndomains-1,"log_R");
	rhs=op->get_rhs("log_R");
	op->add_d(ndomains-1,"log_R","log_R",ones(1,1));
	rhs(ndomains-1)=0;
	op->set_rhs("log_R",rhs);
	
	op->add_d(ndomains-1,"log_M","log_M",ones(1,1));
	op->add_d(ndomains-1,"log_M","m",-ones(1,1)/m);
	op->add_d(ndomains-1,"log_M","log_rhoc",-ones(1,1));
	op->set_rhs("log_M",zeros(1,1));

	rhs=zeros(ndomains,1);
	for(int n=0;n<ndomains;n++) {
		if(n==ndomains-1) {
			op->add_d(n,"Xc","Teff",ones(1,1));
			rhs(n)=-Teff()(0)+Teff_obj;
		} else {
			op->bc_top1_add_d(n,"Xc","Xc",ones(1,1));
			op->bc_top2_add_d(n,"Xc","Xc",-ones(1,1));
		}
	}
	op->set_rhs("Xc",rhs);


}

void star1dR::solve_definitions(solver *op) {

	star1d::solve_definitions(op);

	matrix rho0,xi0,eps0;
	double rhoc0,Xc0;
	
	rho0=rho;xi0=opa.xi;eps0=nuc.eps;rhoc0=rhoc;Xc0=Xc;
	
	double dXc=1e-8;
	
	Xc+=dXc;
	nuclear();
	opacity();
	eq_state();
	matrix drho_dXc,dxi_dXc,deps_dXc,dlnrhoc_dXc;
	
	dlnrhoc_dXc=(rhoc-rhoc0)/rhoc0/dXc*ones(1,1);
	drho_dXc=(rho-rho0)/dXc+rho0*dlnrhoc_dXc;
	dxi_dXc=(opa.xi-xi0)/dXc;
	deps_dXc=(nuc.eps-eps0)/dXc;
	
	Xc=Xc0;
	fill();
	
	for(int n=0;n<ndomains;n++) op->add_d(n,"log_rhoc","Xc",dlnrhoc_dXc);
	op->add_d("rho","Xc",drho_dXc);
	op->add_d("opa.xi","Xc",dxi_dXc);
	op->add_d("nuc.eps","Xc",deps_dXc);

}

