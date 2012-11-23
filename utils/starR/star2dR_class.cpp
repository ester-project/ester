#include"starR.h"
#include<string.h>
#include<stdlib.h>

star2dR::star2dR() {Teff_obj=-1;Re_obj=-1;}

int star2dR::check_arg(char *arg,char *val,int *change_grid) {

	int err=0;

	if(!strcmp(arg,"R")||!strcmp(arg,"Rp")) {
		if(val==NULL) return 2;
		R=atof(val)*R_SUN;
		return 0;
	}
	else if(!strcmp(arg,"Teff")) {
		if(val==NULL) return 2;
		Teff_obj=atof(val);
		return 0;
	}
	else if(!strcmp(arg,"Re")) {
		if(val==NULL) return 2;
		Re_obj=atof(val)*R_SUN;
		return 0;
	}
	else if(!strcmp(arg,"M")) {
		return 1;
	}
	else if(!strcmp(arg,"Xc")) {
		return 1;
	}
	else if(!strcmp(arg,"Omega_bk")) {
		return 1;
	}

	return star2d::check_arg(arg,val,change_grid);
}

solver *star2dR::init_solver(int nvar_add) {

	return star2d::init_solver(2+nvar_add);
}

void star2dR::register_variables(solver *op) {

	star2d::register_variables(op);
	op->regvar("Xc");
	op->regvar("log_M");
}

double star2dR::solve(solver *op) {

	double err,Omega0=Omega;

	err=star2d::solve(op);
	Xc+=op->get_var("Xc")(0);
	while(fabs((Omega-Omega0)/Omegac)>config.newton_dmax) Omega=(Omega+Omega0)/2;
	fill();

	return err;
}

void star2dR::fill() {

	upd_Xr();
	eq_state();
	m=2*PI*(map.gl.I,rho*r*r*map.rz,map.leg.I_00)(0);
	M=m*rhoc*R*R*R;
	double R0=R;
	star2d::fill();
	R=R0;
	Omega_bk=Omega/Omegac;
	if(Teff_obj==-1) Teff_obj=map.leg.eval_00(Teff(),0)(0);
	if(Re_obj==-1) Re_obj=map.leg.eval_00(r.row(-1),PI/2)(0)*R;
}

void star2dR::solve_dim(solver *op) {

	star2d::solve_dim(op);
	
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
			matrix T;
			double Teff_p;
			Teff_p=map.leg.eval_00(Teff(),0,T)(0);
			op->add_r(n,"Xc","Teff",ones(1,1),T);
			rhs(n)=-Teff_p+Teff_obj;
		} else {
			op->bc_top1_add_d(n,"Xc","Xc",ones(1,1));
			op->bc_top2_add_d(n,"Xc","Xc",-ones(1,1));
		}
	}
	op->set_rhs("Xc",rhs);


}

void star2dR::solve_definitions(solver *op) {

	star2d::solve_definitions(op);

	matrix rho0,xi0,eps0;
	double rhoc0,Xc0;
	
	rho0=rho;xi0=opa.xi;eps0=nuc.eps;rhoc0=rhoc;Xc0=Xc;
	
	double dXc=1e-8;
	
	Xc+=dXc;
	upd_Xr();
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

void star2dR::solve_Omega(solver *op) {

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
	op->bc_bot1_add_d(n,"Omega","log_R",ones(1,1));
	op->bc_bot2_add_r(n,"Omega","Ri",ones(1,1)/Req,TT);
	
	rhs(n)=-log(R)-log(Req)+log(Re_obj);
	op->set_rhs("Omega",rhs);

}


