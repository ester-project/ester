#include"star.h"
#include<stdlib.h>
#include<sys/time.h>
#include<omp.h>

void star2d::fill() {
	
	static int initd=0,nth0=0;
	static matrix dTn;
	
	if(nth()!=nth0) {initd=0;nth0=nth();}
	
	upd_Xr();

	eq_state();
	opacity();
	nuclear();
	
	m=2*PI*(map.gl.I,rho*r*r*map.rz,map.leg.I_00)(0);

	R=pow(M/m/rhoc,1./3.);

	pi_c=(4*PI*GRAV*rhoc*rhoc*R*R)/pc;
	Lambda=rhoc*R*R/Tc;
	
	calc_veloc();	
	calc_units();
	
	calc_Frad();

	gsup_=gsup()/units.F;

	atmosphere();

	//Omegac=sqrt(map.leg.eval_00(((Dex,phiex)/rex/map.ex.rz).row(0),PI/2)(0));
	Omegac=sqrt(pi_c*m/4/PI*(1-map.eps(ndomains()-1))*(1-map.eps(ndomains()-1))*(1-map.eps(ndomains()-1)));

	initd=1;
}

void star2d::calc_units() {

	units.phi=pc/rhoc;
	units.p=pc;
	units.rho=rhoc;
	units.T=Tc;
	units.r=R;
	units.Omega=sqrt(pc/rhoc)/R;
	units.v=sqrt(pc/rhoc);
	units.F=pc/R/rhoc;
}

void star2d::calc_Frad() {

	Frad=-opa.xi/Lambda*(sqrt(map.gzz)*(D,T)+map.gzt/sqrt(map.gzz)*(T,Dt));
	
}

void star2d::upd_Xr() {

	int ic,n;
	
	Xr=X*ones(nr(),nth());
	if(!conv) {
		if(Xc!=1) printf("Warning: Non-homogeneus composition without core convection not implemented\n");
		Xc=1;
		return;
	}
	ic=0;
	for(n=0;n<conv;n++) ic+=map.gl.npts[n];
	Xr.setblock(0,ic-1,0,nth()-1,Xc*X*ones(ic,nth()));
	Xr.setblock(ic,nr()-1,0,nth()-1,X*ones(nr()-ic,nth()));
	
}

void star2d::calc_veloc() {

	vr=(G,map.leg.D_11)/r+(map.rt/r+cos(th)/sin(th))/r*G;
	vr.setrow(0,zeros(1,nth()));
	vt=-(D,G)/map.rz-1./r*G;
	vt.setrow(0,zeros(1,nth()));
}

solver *star2d::init_solver() {

	int nvar;
	solver *op;
	
	nvar=21;
	nvar=30;
	op=new solver;
	op->init(ndomains()+1,nvar,"full");
	
	op->maxit_ref=10;op->use_cgs=1;op->maxit_cgs=20;
	op->rel_tol=1e-12;op->abs_tol=1e-20;
	register_variables(op);
	
	return op;
}

void star2d::register_variables(solver *op) {

	int i,var_nr[ndomains()+1];
	
	for(i=0;i<ndomains();i++) 
		var_nr[i]=map.gl.npts[i];
	var_nr[ndomains()]=nex();
	op->set_nr(var_nr);

	op->regvar("phi");
	op->regvar("p");
	op->regvar("pi_c");
	op->regvar("T");
	op->regvar("Lambda");
	op->regvar("eta");
	op->regvar("deta");
	op->regvar("Ri");
	op->regvar("dRi");
	op->regvar("Omega");
	op->regvar("rhoc");
	op->regvar("pc");
	op->regvar("Tc");
	op->regvar("R");
	op->regvar("m");
	op->regvar("ps");
	op->regvar("Ts");
	op->regvar("lum");
	op->regvar("Frad");
	op->regvar("gsup");
	op->regvar("w");
	op->regvar("G");
	//op->regvar("bcw");

}

double star2d::solve(solver *op) {

	int info[5];
	matrix rho0;
	double err,err2,h,dmax;

	if(Omega==0&Omega_bk!=0) {
		Omega=Omega_bk*Omegac;
		w=Omega*ones(nr(),nth());
	}


	op->reset();
	
	solve_poisson(op);
	solve_pressure(op);
	solve_temp(op);
	solve_dim(op);
	solve_map(op);
	solve_Omega(op);
	solve_atm(op);
	solve_gsup(op);
	solve_rot(op);
	solve_dyn(op);

	op->solve(info);
	
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
	
	h=1;
	dmax=config.newton_dmax;
	
	matrix dphi,dphiex,dp,dT,dpc,dTc;
	dphi=op->get_var("phi").block(0,nr()-1,0,nth()-1);
	dphiex=op->get_var("phi").block(nr(),nr()+nex()-1,0,nth()-1);
	err=max(abs(dphi/phi));
	//printf("err(phi)=%e\n",err);
	dp=op->get_var("p");
	err2=max(abs(dp));err=err2>err?err2:err;
	while(exist(abs(h*dp)>dmax)) h/=2;
	//printf("err(p)=%e\n",err2);
	dT=op->get_var("T");
	err2=max(abs(dT));err=err2>err?err2:err;
	while(exist(abs(h*dT)>dmax)) h/=2;
	//printf("err(T)=%e\n",err2);
	dpc=op->get_var("pc");
	err2=fabs(dpc(0));err=err2>err?err2:err;
	while(fabs(h*dpc(0))>dmax) h/=2;
	//printf("err(pc)=%e\n",err2);
	dTc=op->get_var("Tc");
	err2=fabs(dTc(0));err=err2>err?err2:err;
	while(fabs(h*dTc(0))>dmax) h/=2;
	//printf("err(Tc)=%e\n",err2);
	matrix R0,dR;
	R0=map.R;
	dR=op->get_var("Ri").block(1,ndomains(),0,nth()-1);
	while(exist(abs(h*dR)>dmax*R0)) h/=2;
	map.R=R0+h*dR;
	while(map.remap()) {
		h/=2;
		map.R=R0+h*dR;
	}



	
	phi+=h*dphi;
	phiex+=h*dphiex;
	p+=h*dp*p;
	T+=h*dT*T;
	pc*=exp(h*dpc(0));
	Tc*=exp(h*dTc(0));
	Omega=Omega+h*op->get_var("Omega")(0);
	w+=h*op->get_var("w");
	G+=h*op->get_var("G");
	
	rho0=rho;
	
	fill();

	err2=max(abs(rho-rho0));err=err2>err?err2:err;

	
	return err;

}

void star2d::solve_poisson(solver *op) {

	matrix q,rhs1,rhs2,rhs;
	int n,j0;

	// phi
	op->add_l("phi","phi",map.gzz,(D,D));
	rhs1=map.gzz*(D,D,phi);
	q=2*(1+map.rt*map.rzt/r/map.rz)/r/map.rz-map.gzz*map.rzz/map.rz-(map.rtt+map.rt*cos(th)/sin(th))/r/r/map.rz;
	op->add_l("phi","phi",q,D);
	rhs1+=q*(D,phi);
	q=-2*map.rt/r/r/map.rz;
	op->add_lr("phi","phi",q,D,Dt);
	rhs1+=q*(D,phi,Dt);
	op->add_r("phi","phi",1./r/r,map.leg.lap_00);
	rhs1+=(phi,map.leg.lap_00)/r/r;
	
	//rho
	op->add_d("phi","p",-pi_c*rho/eos.chi_rho);
	op->add_d("phi","T",pi_c*rho*eos.d);
	op->add_d("phi","pc",-pi_c*rho/eos.chi_rho);
	op->add_d("phi","Tc",pi_c*rho*eos.d);
	op->add_d("phi","rhoc",pi_c*rho);

	//pi_c
	op->add_d("phi","pi_c",-pi_c*rho);

	// r
	q=-2*map.rt*map.rt/r/r/r/map.rz/map.rz*(D,D,phi)
		+(-2/r/r/map.rz-4*map.rt*map.rzt/r/r/r/map.rz/map.rz+2*map.rzz*map.rt*map.rt/r/r/r/map.rz/map.rz/map.rz+2*map.rtt/r/r/r/map.rz+2*map.rt*cos(th)/r/r/r/map.rz/sin(th))*(D,phi)
		+4*map.rt/r/r/r/map.rz*(D,phi,Dt)
		-2/r/r/r*(phi,map.leg.lap_00);
	add_dr(op,"phi",q);	
	q=(-2/map.rz/map.rz/map.rz-2*map.rt*map.rt/r/r/map.rz/map.rz/map.rz)*(D,D,phi)
		+(-2/r/map.rz/map.rz-4*map.rt*map.rzt/r/r/map.rz/map.rz/map.rz+3*map.rzz/map.rz/map.rz/map.rz/map.rz+3*map.rzz*map.rt*map.rt/r/r/map.rz/map.rz/map.rz/map.rz+map.rtt/r/r/map.rz/map.rz+map.rt*cos(th)/r/r/map.rz/map.rz/sin(th))*(D,phi)
		+2*map.rt/r/r/map.rz/map.rz*(D,phi,Dt);
	add_drz(op,"phi",q);
	q=2*map.rt/r/r/map.rz/map.rz*(D,D,phi)
		+(2*map.rzt/r/r/map.rz/map.rz-2*map.rzz*map.rt/r/r/map.rz/map.rz/map.rz-cos(th)/r/r/map.rz/sin(th))*(D,phi)
		-2/r/r/map.rz*(D,phi,Dt);
	add_drt(op,"phi",q);
	q=(-1/map.rz/map.rz/map.rz-map.rt*map.rt/r/r/map.rz/map.rz/map.rz)*(D,phi);
	add_drzz(op,"phi",q);
	q=2*map.rt/r/r/map.rz/map.rz*(D,phi);
	add_drzt(op,"phi",q);
	q=-1/r/r/map.rz*(D,phi);	
	add_drtt(op,"phi",q);

	// phiex
	op->add_l(ndomains(),"phi","phi",map.ex.gzz,(Dex,Dex));
	rhs2=map.ex.gzz*(Dex,Dex,phiex);
	q=2*(1+map.ex.rt*map.ex.rzt/rex/map.ex.rz)/rex/map.ex.rz-map.ex.gzz*map.ex.rzz/map.ex.rz-(map.ex.rtt+map.ex.rt*cos(th)/sin(th))/rex/rex/map.ex.rz;
	op->add_l(ndomains(),"phi","phi",q,Dex);
	rhs2+=q*(Dex,phiex);
	q=-2*map.ex.rt/rex/rex/map.ex.rz;
	op->add_lr(ndomains(),"phi","phi",q,Dex,Dt);
	rhs2+=q*(Dex,phiex,Dt);
	op->add_r(ndomains(),"phi","phi",1./rex/rex,map.leg.lap_00);
	rhs2+=(phiex,map.leg.lap_00)/rex/rex;


	// rex
	q=-2*map.ex.rt*map.ex.rt/rex/rex/rex/map.ex.rz/map.ex.rz*(Dex,Dex,phiex)
		+(-2/rex/rex/map.ex.rz-4*map.ex.rt*map.ex.rzt/rex/rex/rex/map.ex.rz/map.ex.rz+2*map.ex.rzz*map.ex.rt*map.ex.rt/rex/rex/rex/map.ex.rz/map.ex.rz/map.ex.rz+2*map.ex.rtt/rex/rex/rex/map.ex.rz+2*map.ex.rt*cos(th)/rex/rex/rex/map.ex.rz/sin(th))*(Dex,phiex)
		+4*map.ex.rt/rex/rex/rex/map.ex.rz*(Dex,phiex,Dt)
		-2/rex/rex/rex*(phiex,Dt2)
		-2*cos(th)/rex/rex/rex/sin(th)*(phiex,Dt);
	add_drex(op,"phi",q);	
	q=(-2/map.ex.rz/map.ex.rz/map.ex.rz-2*map.ex.rt*map.ex.rt/rex/rex/map.ex.rz/map.ex.rz/map.ex.rz)*(Dex,Dex,phiex)
		+(-2/rex/map.ex.rz/map.ex.rz-4*map.ex.rt*map.ex.rzt/rex/rex/map.ex.rz/map.ex.rz/map.ex.rz+3*map.ex.rzz/map.ex.rz/map.ex.rz/map.ex.rz/map.ex.rz+3*map.ex.rzz*map.ex.rt*map.ex.rt/rex/rex/map.ex.rz/map.ex.rz/map.ex.rz/map.ex.rz+map.ex.rtt/rex/rex/map.ex.rz/map.ex.rz+map.ex.rt*cos(th)/rex/rex/map.ex.rz/map.ex.rz/sin(th))*(Dex,phiex)
		+2*map.ex.rt/rex/rex/map.ex.rz/map.ex.rz*(Dex,phiex,Dt);
	add_drzex(op,"phi",q);
	q=2*map.ex.rt/rex/rex/map.ex.rz/map.ex.rz*(Dex,Dex,phiex)
		+(2*map.ex.rzt/rex/rex/map.ex.rz/map.ex.rz-2*map.ex.rzz*map.ex.rt/rex/rex/map.ex.rz/map.ex.rz/map.ex.rz-cos(th)/rex/rex/map.ex.rz/sin(th))*(Dex,phiex)
		-2/rex/rex/map.ex.rz*(Dex,phiex,Dt);
	add_drtex(op,"phi",q);
	q=(-1/map.ex.rz/map.ex.rz/map.ex.rz-map.ex.rt*map.ex.rt/rex/rex/map.ex.rz/map.ex.rz/map.ex.rz)*(Dex,phiex);
	add_drzzex(op,"phi",q);
	q=2*map.ex.rt/rex/rex/map.ex.rz/map.ex.rz*(Dex,phiex);
	add_drztex(op,"phi",q);
	q=-1/rex/rex/map.ex.rz*(Dex,phiex);	
	add_drttex(op,"phi",q);

	rhs1=rhs1-pi_c*rho;
	rhs=zeros(nr()+nex(),nth());
	rhs.setblock(0,nr()-1,0,nth()-1,-rhs1);
	rhs.setblock(nr(),nr()+nex()-1,0,nth()-1,-rhs2);
	
	j0=0;
	for(n=0;n<ndomains()+1;n++) {
		if(!n) {
			op->bc_bot2_add_l(n,"phi","phi",ones(1,nth()),D.block(n).row(0));
			rhs.setrow(j0,-(D,phi).row(j0));
		} else {
			if(n<ndomains()) op->bc_bot2_add_l(n,"phi","phi",1/map.rz.row(j0),D.block(n).row(0));
			else op->bc_bot2_add_l(n,"phi","phi",1/map.ex.rz.row(0),Dex.row(0));
			op->bc_bot1_add_l(n,"phi","phi",-1/map.rz.row(j0-1),D.block(n-1).row(map.gl.npts[n-1]-1));
			
			if(n<ndomains()) 
				add_bc_bot2_drz(op,n,"phi",-1/map.rz.row(j0)/map.rz.row(j0)*(D,phi).row(j0));
			else
				add_bc_bot2_drz(op,n,"phi",-1/map.ex.rz.row(0)/map.ex.rz.row(0)*(Dex,phiex).row(0));
			add_bc_bot1_drz(op,n,"phi",1/map.rz.row(j0-1)/map.rz.row(j0-1)*(D,phi).row(j0-1));
			
			if(n<ndomains()) rhs.setrow(j0,-(D,phi).row(j0)/map.rz.row(j0)+(D,phi).row(j0-1)/map.rz.row(j0-1));
			else rhs.setrow(j0,-(Dex,phiex).row(0)/map.ex.rz.row(0)+(D,phi).row(j0-1)/map.rz.row(j0-1));
		}
		
		op->bc_top1_add_d(n,"phi","phi",ones(1,nth()));
		if(n<ndomains()) op->bc_top2_add_d(n,"phi","phi",-ones(1,nth()));
		if(n<ndomains()) rhs.setrow(j0+map.gl.npts[n]-1,-phi.row(j0+map.gl.npts[n]-1));
		else rhs.setrow(nr()+nex()-1,-phiex.row(nex()-1));
		if(n<ndomains()-1) rhs.setrow(j0+map.gl.npts[n]-1,rhs.row(j0+map.gl.npts[n]-1)
								+phi.row(j0+map.gl.npts[n]));
		else if(n==ndomains()-1) rhs.setrow(j0+map.gl.npts[n]-1,rhs.row(j0+map.gl.npts[n]-1)
								+phiex.row(0)); 
		
		if(n<ndomains()) j0+=map.gl.npts[n];
	}
	op->set_rhs("phi",rhs);
}

void star2d::solve_pressure(solver *op) {

	matrix q,rhs_p,rhs_pi_c;
	int n,j0;
	
	// p	
	q=ones(nr(),nth());
	op->add_li("p","p",q,D,p);
	
	// phi
	q=rho;
	op->add_l("p","phi",q,D);

	//rho
	q=(D,phi)-w*w*r*map.rz*sin(th)*sin(th);
	op->add_d("p","p",q*rho/eos.chi_rho);
	op->add_d("p","T",-q*rho*eos.d);
	op->add_d("p","pc",q*rho/eos.chi_rho);
	op->add_d("p","Tc",-q*rho*eos.d);
	op->add_d("p","rhoc",-q*rho);

	//w
	q=-rho*r*map.rz*sin(th)*sin(th);
	op->add_d("p","w",2*w*q);

	// r
	q=-rho*w*w*map.rz*sin(th)*sin(th);
	add_dr(op,"p",q);
	q=-rho*w*w*r*sin(th)*sin(th);
	add_drz(op,"p",q);

	rhs_p=-(D,p)-rho*(D,phi)+rho*w*w*r*map.rz*sin(th)*sin(th);
	rhs_pi_c=zeros(ndomains(),1);

	j0=0;
	for(n=0;n<ndomains();n++) {
		if(!n) {
			op->bc_bot2_add_d(n,"p","p",p.row(0));
			rhs_p.setrow(0,1-p.row(0));	
		} else {
			op->bc_bot2_add_d(n,"p","p",p.row(j0));
			op->bc_bot1_add_d(n,"p","p",-p.row(j0-1));
			rhs_p.setrow(j0,-p.row(j0)+p.row(j0-1));	
		}
		if(n<ndomains()-1) {
			op->bc_top1_add_d(n,"pi_c","pi_c",ones(1,1));
			op->bc_top2_add_d(n,"pi_c","pi_c",-ones(1,1));
		} else {
			map.leg.eval_00(th,0,q);
			op->bc_top1_add_ri(n,"pi_c","ps",-ones(1,1),q,ps);
			op->bc_top1_add_ri(n,"pi_c","p",ones(1,1),q,p.row(nr()-1));
			rhs_pi_c(ndomains()-1)=(ps-p.row(nr()-1),q)(0);
		}
		
		j0+=map.gl.npts[n];
	}
	op->set_rhs("p",rhs_p);
	op->set_rhs("pi_c",rhs_pi_c);
}


void star2d::solve_rot(solver *op) {

	matrix q,TT,rhs;
	int n,j0;
	
	if(Omega==0) {
		op->add_d("w","w",ones(nr(),nth()));
		op->add_d("w","Omega",-ones(nr(),nth()));
		op->set_rhs("w",zeros(nr(),nth()));
		return;
	}
	
	// w
	
	q=(r*cos(th)+map.rt*sin(th))/map.rz;
	op->add_l("w","w",2*w*q,D);
	q=-sin(th)*ones(nr(),nth());
	op->add_r("w","w",2*w*q,Dt);
	q=(r*cos(th)+map.rt*sin(th))/map.rz*2*(D,w)-sin(th)*2*(w,Dt);
	op->add_d("w","w",q);
	
	// p
	
	q=-(rho,Dt)/rho/rho/r/map.rz/sin(th);
	op->add_li("w","p",q,D,p);
	q=(D,rho)/rho/rho/r/map.rz/sin(th);
	op->add_ri("w","p",q,Dt,p);

	//rho
	q=-2./rho/rho/r/map.rz/sin(th)*((p,Dt)*(D,rho)-(D,p)*(rho,Dt));
	op->add_d("w","p",q/eos.chi_rho);
	op->add_d("w","T",-q*eos.d);
	op->add_d("w","pc",q/eos.chi_rho);
	op->add_d("w","Tc",-q*eos.d);
	op->add_d("w","rhoc",-q);
	q=1./rho/rho/r/map.rz/sin(th)*(p,Dt);
	op->add_li("w","p",q,D,rho/eos.chi_rho);
	op->add_li("w","T",-q,D,rho*eos.d);
	op->add_d("w","pc",q*(D,rho/eos.chi_rho));
	op->add_d("w","Tc",-q*(D,rho*eos.d));
	op->add_d("w","rhoc",-q*(D,rho));
	q=-1./rho/rho/r/map.rz/sin(th)*(D,p);
	op->add_ri("w","p",q,Dt,rho/eos.chi_rho);
	op->add_ri("w","T",-q,Dt,eos.d*rho);
	op->add_d("w","pc",q*(rho/eos.chi_rho,Dt));
	op->add_d("w","Tc",-q*(eos.d*rho,Dt));
	op->add_d("w","rhoc",-q*(rho,Dt));
	
	// r
	q=cos(th)/map.rz*2*w*(D,w)-1./rho/rho/r/r/map.rz/sin(th)*((p,Dt)*(D,rho)-(D,p)*(rho,Dt));
	add_dr(op,"w",q);
	q=-(r*cos(th)+map.rt*sin(th))/map.rz/map.rz*2*w*(D,w)-1./rho/rho/r/map.rz/map.rz/sin(th)*((p,Dt)*(D,rho)-(D,p)*(rho,Dt));
	add_drz(op,"w",q);
	q=sin(th)/map.rz*2*w*(D,w);
	add_drt(op,"w",q);

	rhs=-(r*cos(th)+map.rt*sin(th))/map.rz*2*w*(D,w)+
			sin(th)*2*w*(w,Dt)-1./rho/rho/r/map.rz/sin(th)*((p,Dt)*(D,rho)-(D,p)*(rho,Dt));
	
	j0=0;
	for(n=0;n<ndomains();n++) {
		if(!n) {
			op->bc_bot2_add_l(n,"w","w",ones(1,nth()),D.block(n).row(0));
			rhs.setrow(0,-(D,w).row(0));
		}
		if(n<ndomains()-1) {
			q=-rho*r*sin(th)*(r*cos(th)+map.rt*sin(th));
			op->bc_top1_add_d(n,"w","w",2.*(w*q).row(j0+map.gl.npts[n]-1));
			op->bc_top2_add_d(n,"w","w",-2.*(w*q).row(j0+map.gl.npts[n]));
			op->bc_top1_add_r(n,"w","phi",rho.row(j0+map.gl.npts[n]-1),Dt);
			op->bc_top2_add_r(n,"w","phi",-rho.row(j0+map.gl.npts[n]),Dt);
			q=-2.*rho*w*w*sin(th)*(r*cos(th)+map.rt*sin(th));
			add_bc_top1_dr(op,n,"w",q.row(j0+map.gl.npts[n]-1));
			add_bc_top2_dr(op,n,"w",-q.row(j0+map.gl.npts[n]));
			q=-rho*w*w*r*sin(th)*sin(th);
			add_bc_top1_drt(op,n,"w",q.row(j0+map.gl.npts[n]-1));
			add_bc_top2_drt(op,n,"w",-q.row(j0+map.gl.npts[n]));
			q=(phi,Dt)-w*w*r*sin(th)*(r*cos(th)+map.rt*sin(th));
			op->bc_top1_add_d(n,"w","p",(q*rho/eos.chi_rho).row(j0+map.gl.npts[n]-1));
			op->bc_top1_add_d(n,"w","T",(-q*rho*eos.d).row(j0+map.gl.npts[n]-1));
			op->bc_top1_add_d(n,"w","pc",(q*rho/eos.chi_rho).row(j0+map.gl.npts[n]-1));
			op->bc_top1_add_d(n,"w","Tc",(-q*rho*eos.d).row(j0+map.gl.npts[n]-1));
			op->bc_top1_add_d(n,"w","rhoc",(-q*rho).row(j0+map.gl.npts[n]-1));
			op->bc_top2_add_d(n,"w","p",-(q*rho/eos.chi_rho).row(j0+map.gl.npts[n]));
			op->bc_top2_add_d(n,"w","T",-(-q*rho*eos.d).row(j0+map.gl.npts[n]));
			op->bc_top2_add_d(n,"w","pc",-(q*rho/eos.chi_rho).row(j0+map.gl.npts[n]));
			op->bc_top2_add_d(n,"w","Tc",-(-q*rho*eos.d).row(j0+map.gl.npts[n]));
			op->bc_top2_add_d(n,"w","rhoc",-(-q*rho).row(j0+map.gl.npts[n]));
			q=rho*(phi,Dt)-rho*w*w*r*sin(th)*(r*cos(th)+map.rt*sin(th));
			rhs.setrow(j0+map.gl.npts[n]-1,
					-q.row(j0+map.gl.npts[n]-1)+q.row(j0+map.gl.npts[n]));
		}
		
		j0+=map.gl.npts[n];
	}
	
	
	
	solve_vbl(op,"w",rhs);
	op->set_rhs("w",rhs);
}

void star2d::solve_vbl(solver *op,const char *eqn,matrix &rhs) {

	matrix q,qeq,TT;
	int n;
	int limit_layer=1;
	
	if(!limit_layer) {
		qeq=ones(1,nth()); // w constante
	} else {
		qeq=zeros(1,nth()); //limit layer
	}
	qeq(0)=1;
	n=ndomains()-1;
	
	op->bc_top1_add_l(n,eqn,"w",(1-qeq)*map.gzz.row(nr()-1),D.block(n).row(map.gl.npts[n]-1));
	q=map.gzt+G/r/map.rz;
	op->bc_top1_add_r(n,eqn,"w",(1-qeq)*q.row(nr()-1),Dt);
	q=2*G/r/map.rz*(map.rt/r+cos(th)/sin(th));
	op->bc_top1_add_d(n,eqn,"w",(1-qeq)*q.row(nr()-1));
	q=(w,Dt)/r/map.rz+2*(map.rt/r+cos(th)/sin(th))/r/map.rz*w;
	op->bc_top1_add_d(n,eqn,"G",(1-qeq)*q.row(nr()-1));
	q=-2.*r*map.gzt*map.gzt*(D,w)-(2.*map.gzt+G/r/map.rz)/r*(w,Dt)
		-2*G/r/r/map.rz*(2*map.rt/r+cos(th)/sin(th))*w;
	add_bc_top1_dr(op,n,eqn,(1-qeq)*q.row(nr()-1));
	q=-2./map.rz*map.gzz*(D,w)-(map.gzt+G/r/map.rz)/map.rz*(w,Dt)
		-2*G/r/map.rz/map.rz*(map.rt/r+cos(th)/sin(th))*w;
	add_bc_top1_drz(op,n,eqn,(1-qeq)*q.row(nr()-1));
	q=-2./map.rz*map.gzt*(D,w)-1./r/r/map.rz*(w,Dt)
		+2*G/r/r/map.rz*w;
	add_bc_top1_drt(op,n,eqn,(1-qeq)*q.row(nr()-1));
	
	q=map.gzz*(D,w)+(map.gzt+G/r/map.rz)*(w,Dt)
		+2*G/r/map.rz*(map.rt/r+cos(th)/sin(th))*w;
	rhs.setrow(nr()-1,-(1-qeq)*q.row(nr()-1));

	if(!limit_layer) {
		rhs.setrow(nr()-1,-w.row(nr()-1)+Omega);
		op->bc_top1_add_d(n,eqn,"w",qeq);
		op->bc_top1_add_d(n,eqn,"Omega",-qeq);
	} else {
		map.leg.eval_00(th,PI/2*ones(1,nth()),TT);
		rhs.setrow(nr()-1,rhs.row(nr()-1)+qeq*(-(w.row(nr()-1),TT)(0)+Omega));
		op->bc_top1_add_r(n,eqn,"w",qeq,TT);
		op->bc_top1_add_d(n,eqn,"Omega",-qeq);
	}	

}

void star2d::solve_dyn(solver *op) {

	matrix rhs;
	matrix q;
	int n,j0;

	if(Omega==0) {
		op->add_d("G","G",ones(nr(),nth()));
		op->set_rhs("G",zeros(nr(),nth()));
		return;
	}
	
	op->add_l("G","w",map.gzz,(D,D));
	op->add_lr("G","w",2.*map.gzt,D,Dt);
	op->add_r("G","w",1./r/r,Dt2);
	q=4./r/map.rz-map.rtt/r/r/map.rz
		+map.gzz*(-map.rzz/map.rz)
		+map.gzt*(3.*cos(th)/sin(th)-2.*map.rzt/map.rz)-vr/map.rz;
	op->add_l("G","w",q,D);
	q=1./r/r*(3.*cos(th)/sin(th))
		-vt/r;
	op->add_r("G","w",q,Dt);
	
	q=-2.*map.rz/r*vr/map.rz-2.*(map.rt/r+cos(th)/sin(th))*vt/r;
	op->add_d("G","w",q);
	
	matrix qv;

	// vr
	qv=-(D,w)/map.rz-2./r*w;
	op->add_r("G","G",qv/r,map.leg.D_11);
	q=(map.rt/r+cos(th)/sin(th))/r;
	op->add_d("G","G",q*qv);
	q=-((G,map.leg.D_11)+(2*map.rt/r+cos(th)/sin(th))*G)/r/r;
	add_dr(op,"G",q*qv);
	q=G/r/r;
	add_drt(op,"G",q*qv);
	
	
	//vt
	qv=-(w,Dt)/r-2./r*(map.rt/r+cos(th)/sin(th))*w;
	op->add_l("G","G",-qv/map.rz,D);
	q=-1./r;
	op->add_d("G","G",qv*q);
	q=G/r/r;
	add_dr(op,"G",q*qv);
	q=(D,G)/map.rz/map.rz;
	add_drz(op,"G",q*qv);
	
	//r
	q=-2.*r*map.gzt*map.gzt*(D,D,w)-4.*map.gzt/r*(D,w,Dt)-2./r/r/r*(w,Dt2)
		-(4./r/r/map.rz-2.*map.rtt/r/r/r/map.rz
			+2.*r*map.gzt*map.gzt*(-map.rzz/map.rz)
			+2.*map.gzt/r*(3.*cos(th)/sin(th)-2.*map.rzt/map.rz)
			)*(D,w)
		-rho*(2./r/r/r*(3.*cos(th)/sin(th)+(rho,Dt)/rho)+2.*map.gzt/r*(D,rho)/rho
			-vt/r/r/rho)*(w,Dt)
		+2./r/r*(vr+(2*map.rt/r+cos(th)/sin(th))*vt)*w;
	add_dr(op,"G",q);
	q=-2.*map.gzz/map.rz*(D,D,w)-2.*map.gzt/map.rz*(D,w,Dt)
		-(4./r/map.rz/map.rz-map.rtt/r/r/map.rz/map.rz
			+2.*map.gzz/map.rz*(-1.5*map.rzz/map.rz)
			+map.gzt/map.rz*(3.*cos(th)/sin(th)-4.*map.rzt/map.rz)
			-vr/map.rz/map.rz)*(D,w);
	add_drz(op,"G",q);
	q=-2.*map.gzt/map.rz*(D,D,w)-2./r/r/map.rz*(D,w,Dt)
		-(2.*map.gzt/map.rz*(-map.rzz/map.rz)
			+1./r/r/map.rz*(+3.*cos(th)/sin(th)-2.*map.rzt/map.rz)
			)*(D,w)
		-2./r/r*vt*w;
	add_drt(op,"G",q);
	q=-map.gzz/map.rz*(D,w);
	add_drzz(op,"G",q);
	q=-2.*map.gzt/map.rz*(D,w);
	add_drzt(op,"G",q);
	q=-1./r/r/map.rz*(D,w);
	add_drtt(op,"G",q);
	
	rhs=-map.gzz*(D,D,w)-2.*map.gzt*(D,w,Dt)-1./r/r*(w,Dt2)
		-(4./r/map.rz-map.rtt/r/r/map.rz
			+map.gzz*(-map.rzz/map.rz)
			+map.gzt*(3.*cos(th)/sin(th)-2.*map.rzt/map.rz)
			-vr/map.rz)*(D,w)
		-(1./r/r*(3.*cos(th)/sin(th))
			-vt/r)*(w,Dt)
		+2*w*(1./r*vr+(map.rt/r+cos(th)/sin(th))/r*vt);	
	
	
	j0=0;
	for(n=0;n<ndomains();n++) {
		if(!n) {
			op->bc_bot2_add_d(n,"G","G",ones(1,nth()));
			rhs.setrow(0,-G.row(0));
		} else {
			q=(r*r+map.rt*map.rt)/map.rz;
			op->bc_bot2_add_l(n,"G","w",q.row(j0),D.block(n).row(0));
			op->bc_bot1_add_l(n,"G","w",-q.row(j0-1),D.block(n-1).row(map.gl.npts[n-1]-1));
			q=G*r-map.rt;
			op->bc_bot2_add_r(n,"G","w",q.row(j0),Dt);
			op->bc_bot1_add_r(n,"G","w",-q.row(j0-1),Dt);
			q=2*G*(r*cos(th)/sin(th)+map.rt);
			q=r*(w,Dt)+2*w*(r*cos(th)/sin(th)+map.rt);
			op->bc_bot2_add_d(n,"G","G",q.row(j0));
			op->bc_bot1_add_d(n,"G","G",-q.row(j0-1));
			q=2*r/map.rz*(D,w)+G*(w,Dt)+2*w*G*cos(th)/sin(th);
			add_bc_bot2_dr(op,n,"G",q.row(j0));
			add_bc_bot1_dr(op,n,"G",-q.row(j0-1));
			q=-(r*r+map.rt*map.rt)/map.rz/map.rz*(D,w);
			add_bc_bot2_drz(op,n,"G",q.row(j0));
			add_bc_bot1_drz(op,n,"G",-q.row(j0-1));
			q=2*map.rt/map.rz*(D,w)-(w,Dt)+2*w*G;
			add_bc_bot2_drt(op,n,"G",q.row(j0));
			add_bc_bot1_drt(op,n,"G",-q.row(j0-1));
			q=(r*r+map.rt*map.rt)/map.rz*(D,w)+(G*r-map.rt)*(w,Dt)+2*w*G*(r*cos(th)/sin(th)+map.rt);
			rhs.setrow(j0,-q.row(j0)+q.row(j0-1));
		}
	
		j0+=map.gl.npts[n];
	}
	
	
	op->set_rhs("G",rhs);

	
}



void star2d::solve_temp(solver *op) {

	int n,j0,j1;
	matrix q,q1,q2,TT,lum,qconv,qrad;
	matrix rhs_T,rhs_Lambda,rhs_lum,rhs_Frad;

	//Luminosity

	lum=zeros(ndomains(),1);
	j0=0;
	for(n=0;n<ndomains();n++) {
		if(n) lum(n)=lum(n-1);
		lum(n)+=2*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),
			(rho*nuc.eps*r*r*map.rz).block(j0,j0+map.gl.npts[n]-1,0,nth()-1),map.leg.I_00)(0);
		j0+=map.gl.npts[n];
	}

	rhs_lum=zeros(ndomains(),1);
	j0=0;
	for(n=0;n<ndomains();n++) {
		op->bc_bot2_add_d(n,"lum","lum",ones(1,1));
		//rho
		q=-2*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(r*r*map.rz*rho*nuc.eps).block(j0,j0+map.gl.npts[n]-1,0,nth()-1),map.leg.I_00);
		op->bc_bot2_add_lri(n,"lum","p",-2*PI*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,(r*r*map.rz*rho*nuc.eps/eos.chi_rho*(1.+nuc.dlneps_lnrho)).block(j0,j0+map.gl.npts[n]-1,0,nth()-1));
		op->bc_bot2_add_lri(n,"lum","T",-2*PI*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,(r*r*map.rz*rho*nuc.eps*(nuc.dlneps_lnT-eos.d*(1.+nuc.dlneps_lnrho))).block(j0,j0+map.gl.npts[n]-1,0,nth()-1));
		q=-2*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(r*r*map.rz*rho*nuc.eps/eos.chi_rho*(1.+nuc.dlneps_lnrho)).block(j0,j0+map.gl.npts[n]-1,0,nth()-1),map.leg.I_00);
		op->bc_bot2_add_d(n,"lum","pc",q);
		q=-2*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(r*r*map.rz*rho*nuc.eps*(nuc.dlneps_lnT-eos.d*(1.+nuc.dlneps_lnrho))).block(j0,j0+map.gl.npts[n]-1,0,nth()-1),map.leg.I_00);
		op->bc_bot2_add_d(n,"lum","Tc",q);
		q=-2*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(-r*r*map.rz*rho*nuc.eps).block(j0,j0+map.gl.npts[n]-1,0,nth()-1),map.leg.I_00);
		op->bc_bot2_add_d(n,"lum","rhoc",q);
		//r (rz)
		q=-2*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(rho*nuc.eps*(D,r*r*map.J[0])).block(j0,j0+map.gl.npts[n]-1,0,nth()-1),map.leg.I_00);
		op->bc_bot2_add_d(n,"lum","eta",q);
		q=-2*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(rho*nuc.eps*(D,r*r*map.J[1])).block(j0,j0+map.gl.npts[n]-1,0,nth()-1),map.leg.I_00);
		op->bc_bot2_add_d(n,"lum","deta",q);
		q=(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(rho*nuc.eps*(D,r*r*map.J[2])).block(j0,j0+map.gl.npts[n]-1,0,nth()-1));
		op->bc_bot2_add_ri(n,"lum","Ri",-2*PI*ones(1,1),map.leg.I_00,q);
		q=(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(rho*nuc.eps*(D,r*r*map.J[3])).block(j0,j0+map.gl.npts[n]-1,0,nth()-1));
		op->bc_bot2_add_ri(n,"lum","dRi",-2*PI*ones(1,1),map.leg.I_00,q);
		
		if(n) op->bc_bot1_add_d(n,"lum","lum",-ones(1,1));
		j0+=map.gl.npts[n];
	}
	op->set_rhs("lum",rhs_lum);
	
	//Flux
	
	rhs_Frad=zeros(ndomains()*2-1,nth());
	q1=opa.dlnxi_lnrho/eos.chi_rho;
	q2=(opa.dlnxi_lnT-eos.d*opa.dlnxi_lnrho);
	j0=0;
	for(n=0;n<ndomains();n++) {
		j1=j0+map.gl.npts[n]-1;
		
		if(n) op->bc_bot2_add_d(n,"Frad","Frad",ones(1,nth()));
		op->bc_top1_add_d(n,"Frad","Frad",ones(1,nth()));
		
		q=opa.xi/Lambda*sqrt(map.gzz);
		if(n) op->bc_bot2_add_li(n,"Frad","T",q.row(j0),D.block(n).row(0),T.block(j0,j1,0,nth()-1));
		op->bc_top1_add_li(n,"Frad","T",q.row(j1),D.block(n).row(map.gl.npts[n]-1),T.block(j0,j1,0,nth()-1));
		q=opa.xi/Lambda*map.gzt/sqrt(map.gzz);
		if(n) op->bc_bot2_add_ri(n,"Frad","T",q.row(j0),Dt,T.row(j0));
		op->bc_top1_add_ri(n,"Frad","T",q.row(j1),Dt,T.row(j1));
		
		if(n) op->bc_bot2_add_d(n,"Frad","Lambda",Frad.row(j0));
		op->bc_top1_add_d(n,"Frad","Lambda",Frad.row(j1));
		
		if(n) op->bc_bot2_add_d(n,"Frad","p",-Frad.row(j0)*q1.row(j0));
		if(n) op->bc_bot2_add_d(n,"Frad","pc",-Frad.row(j0)*q1.row(j0));
		if(n) op->bc_bot2_add_d(n,"Frad","T",-Frad.row(j0)*q2.row(j0));
		if(n) op->bc_bot2_add_d(n,"Frad","Tc",-Frad.row(j0)*q2.row(j0));
		op->bc_top1_add_d(n,"Frad","p",-Frad.row(j1)*q1.row(j1));
		op->bc_top1_add_d(n,"Frad","pc",-Frad.row(j1)*q1.row(j1));
		op->bc_top1_add_d(n,"Frad","T",-Frad.row(j1)*q2.row(j1));
		op->bc_top1_add_d(n,"Frad","Tc",-Frad.row(j1)*q2.row(j1));
		
		q=-opa.xi/Lambda*map.gzt/sqrt(map.gzz)*(r*map.gzt*(D,T)+(2./r-r*map.gzt/map.gzz)*(T,Dt));
		if(n) add_bc_bot2_dr(op,n,"Frad",q.row(j0));
		add_bc_top1_dr(op,n,"Frad",q.row(j1));
		q=-opa.xi/Lambda*sqrt(map.gzz)/map.rz*(D,T);
		if(n) add_bc_bot2_drz(op,n,"Frad",q.row(j0));
		add_bc_top1_drz(op,n,"Frad",q.row(j1));
		q=-opa.xi/Lambda/sqrt(map.gzz)/map.rz*(map.gzt*(D,T)-(-1./r/r+map.gzt*map.gzt/map.gzz)*(T,Dt));
		if(n) add_bc_bot2_drt(op,n,"Frad",q.row(j0));
		add_bc_top1_drt(op,n,"Frad",q.row(j1));
		
		j0=j1+1;
	}
	op->set_rhs("Frad",rhs_Frad);
	
	//Temperature
	
	qrad=zeros(nr(),nth());
	qconv=qrad;
	j0=0;
	for(n=0;n<ndomains();n++) {
		if(n<conv) qconv.setblock(j0,j0+map.gl.npts[n]-1,0,nth()-1,ones(map.gl.npts[n],nth()));
		else qrad.setblock(j0,j0+map.gl.npts[n]-1,0,nth()-1,ones(map.gl.npts[n],nth()));
		j0+=map.gl.npts[n];
	}
	
	
	rhs_T=zeros(nr(),nth());

	// T
	op->add_li("T","T",qrad*map.gzz,(D,D),T);
	rhs_T+=-qrad*map.gzz*(D,D,T);
	q=2*(1+map.rt*map.rzt/r/map.rz)/r/map.rz-map.gzz*map.rzz/map.rz-(map.rtt+map.rt*cos(th)/sin(th))/r/r/map.rz
		+map.gzz*(D,log(opa.xi))+map.gzt*(log(opa.xi),Dt);
	op->add_li("T","T",qrad*q,D,T);
	rhs_T+=-qrad*q*(D,T);
	q=-2*map.rt/r/r/map.rz;
	op->add_lri("T","T",qrad*q,D,Dt,T);
	rhs_T+=-qrad*q*(D,T,Dt);
	op->add_ri("T","T",qrad*1./r/r,map.leg.lap_00,T);
	rhs_T+=-qrad*(T,map.leg.lap_00)/r/r;
	q=map.gzt*(D,log(opa.xi))+map.gtt*(log(opa.xi),Dt);
	op->add_ri("T","T",qrad*q,Dt,T);
	rhs_T+=-qrad*q*(T,Dt);
	
	// r
	q=-2*map.rt*map.rt/r/r/r/map.rz/map.rz*(D,D,T)
		+(-2/r/r/map.rz-4*map.rt*map.rzt/r/r/r/map.rz/map.rz+2*map.rzz*map.rt*map.rt/r/r/r/map.rz/map.rz/map.rz+2*map.rtt/r/r/r/map.rz+2*map.rt*cos(th)/r/r/r/map.rz/sin(th))*(D,T)
		+4*map.rt/r/r/r/map.rz*(D,T,Dt)
		-2/r/r/r*(T,map.leg.lap_00)
		+(D,log(opa.xi))*(D,T)*2*map.rt/r/map.rz*map.gzt
		-((D,log(opa.xi))*(T,Dt)+(log(opa.xi),Dt)*(D,T))*2*map.gzt/r
		-(log(opa.xi),Dt)*(T,Dt)*2/r/r/r;
	add_dr(op,"T",qrad*q);	
	q=(-2/map.rz/map.rz/map.rz-2*map.rt*map.rt/r/r/map.rz/map.rz/map.rz)*(D,D,T)
		+(-2/r/map.rz/map.rz-4*map.rt*map.rzt/r/r/map.rz/map.rz/map.rz+3*map.rzz/map.rz/map.rz/map.rz/map.rz+3*map.rzz*map.rt*map.rt/r/r/map.rz/map.rz/map.rz/map.rz+map.rtt/r/r/map.rz/map.rz+map.rt*cos(th)/r/r/map.rz/map.rz/sin(th))*(D,T)
		+2*map.rt/r/r/map.rz/map.rz*(D,T,Dt)
		-(D,log(opa.xi))*(D,T)*2/map.rz*map.gzz
		-((D,log(opa.xi))*(T,Dt)+(log(opa.xi),Dt)*(D,T))*map.gzt/map.rz;
	add_drz(op,"T",qrad*q);
	q=2*map.rt/r/r/map.rz/map.rz*(D,D,T)
		+(2*map.rzt/r/r/map.rz/map.rz-2*map.rzz*map.rt/r/r/map.rz/map.rz/map.rz-cos(th)/r/r/map.rz/sin(th))*(D,T)
		-2/r/r/map.rz*(D,T,Dt)
		-(D,log(opa.xi))*(D,T)*2/map.rz*map.gzt
		-((D,log(opa.xi))*(T,Dt)+(log(opa.xi),Dt)*(D,T))/r/r/map.rz;
	add_drt(op,"T",qrad*q);
	q=(-1/map.rz/map.rz/map.rz-map.rt*map.rt/r/r/map.rz/map.rz/map.rz)*(D,T);
	add_drzz(op,"T",qrad*q);
	q=2*map.rt/r/r/map.rz/map.rz*(D,T);
	add_drzt(op,"T",qrad*q);
	q=-1/r/r/map.rz*(D,T);	
	add_drtt(op,"T",qrad*q);
	
	//rho
	q=Lambda*nuc.eps/opa.xi;
	op->add_d("T","p",qrad*q*rho/eos.chi_rho);
	op->add_d("T","T",-qrad*q*rho*eos.d);
	op->add_d("T","pc",qrad*q*rho/eos.chi_rho);
	op->add_d("T","Tc",-qrad*q*rho*eos.d);
	op->add_d("T","rhoc",-qrad*q*rho);
	
	rhs_T+=-qrad*Lambda*rho*nuc.eps/opa.xi;

	//opa.xi
	q1=opa.dlnxi_lnrho/eos.chi_rho;
	q2=(opa.dlnxi_lnT-eos.d*opa.dlnxi_lnrho);
	q=map.gzz*(D,T)+map.gzt*(T,Dt);
	op->add_l("T","p",qrad*q*q1,D);
	op->add_d("T","p",qrad*q*(D,q1));
	op->add_d("T","pc",qrad*q*(D,q1));
	op->add_l("T","T",qrad*q*q2,D);
	op->add_d("T","T",qrad*q*(D,q2));
	op->add_d("T","Tc",qrad*q*(D,q2));
	q=map.gzt*(D,T)+map.gtt*(T,Dt);
	op->add_r("T","p",qrad*q*q1,Dt);
	op->add_d("T","p",qrad*q*(q1,Dt));
	op->add_d("T","pc",qrad*q*(q1,Dt));
	op->add_r("T","T",qrad*q*q2,Dt);
	op->add_d("T","T",qrad*q*(q2,Dt));
	op->add_d("T","Tc",qrad*q*(q2,Dt));
	q=-Lambda*rho*nuc.eps/opa.xi;
	op->add_d("T","p",qrad*q*q1);
	op->add_d("T","pc",qrad*q*q1);
	op->add_d("T","T",qrad*q*q2);
	op->add_d("T","Tc",qrad*q*q2);
	
	//nuc.eps
	q1=nuc.dlneps_lnrho/eos.chi_rho;
	q2=(nuc.dlneps_lnT-eos.d*nuc.dlneps_lnrho);
	q=Lambda*rho*nuc.eps/opa.xi;
	op->add_d("T","p",qrad*q*q1);
	op->add_d("T","pc",qrad*q*q1);
	op->add_d("T","T",qrad*q*q2);
	op->add_d("T","Tc",qrad*q*q2);
	
	//Lambda
	q=Lambda*rho*nuc.eps/opa.xi;
	op->add_d("T","Lambda",qrad*q);
	
	
	op->add_l("T","T",qconv,D);
	op->add_d("T","T",qconv*(D,eos.cp)/eos.cp);
	op->add_l("T","p",-qconv*eos.del_ad,D);
	op->add_d("T","p",-qconv*(D,eos.cp*eos.del_ad)/eos.cp);
	op->add_d("T","Tc",qconv*(D,eos.cp)/eos.cp);
	op->add_d("T","pc",-qconv*(D,eos.cp*eos.del_ad)/eos.cp);
	rhs_T+=-qconv*((D,log(T))-eos.del_ad*(D,log(p)));
	
	rhs_Lambda=zeros(ndomains(),1);
	
	map.leg.eval_00(th,0,TT);
	
	j0=0;
	for(n=0;n<ndomains();n++) {
		if(!n) {
			op->bc_bot2_add_d(n,"T","T",T.row(j0));
			rhs_T.setrow(j0,1-T.row(j0));
		} else {
			op->bc_bot2_add_d(n,"T","T",T.row(j0));
			op->bc_bot1_add_d(n,"T","T",-T.row(j0));
			rhs_T.setrow(j0,-T.row(j0)+T.row(j0-1));
		}
		if(n>=conv) {
			if(n<ndomains()-1) {
			
				op->bc_top1_add_li(n,"T","T",1/map.rz.row(j0+map.gl.npts[n]-1),D.block(n).row(map.gl.npts[n]-1),T.block(j0,j0+map.gl.npts[n]-1,0,nth()-1));
				op->bc_top2_add_li(n,"T","T",-1/map.rz.row(j0+map.gl.npts[n]),D.block(n+1).row(0),T.block(j0+map.gl.npts[n],j0+map.gl.npts[n]+map.gl.npts[n+1]-1,0,nth()-1));
				
				add_bc_top1_drz(op,n,"T",-1/map.rz.row(j0+map.gl.npts[n]-1)/map.rz.row(j0+map.gl.npts[n]-1)*(D,T).row(j0+map.gl.npts[n]-1));
				add_bc_top2_drz(op,n,"T",1/map.rz.row(j0+map.gl.npts[n])/map.rz.row(j0+map.gl.npts[n])*(D,T).row(j0+map.gl.npts[n]));
				
				rhs_T.setrow(j0+map.gl.npts[n]-1,
					-(D,T).row(j0+map.gl.npts[n]-1)/map.rz.row(j0+map.gl.npts[n]-1)
					+(D,T).row(j0+map.gl.npts[n])/map.rz.row(j0+map.gl.npts[n]));
			} else {
				op->bc_top1_add_d(n,"T","T",T.row(nr()-1));
				op->bc_top1_add_d(n,"T","Ts",-Ts);
				rhs_T.setrow(nr()-1,Ts-T.row(nr()-1));
			}
		}
		
		if(n<conv) {
			op->bc_top1_add_d(n,"Lambda","Lambda",ones(1,1));
			op->bc_top2_add_d(n,"Lambda","Lambda",-ones(1,1));
		} else if(n==conv) {
			if(!n) {
				map.leg.eval_00(th,PI/2,q);
				op->bc_bot2_add_lri(n,"Lambda","T",ones(1,1),D.block(0).row(0),q,T.block(j0,j0+map.gl.npts[n]-1,0,nth()-1));
				rhs_Lambda(0)=-((D,T).row(0),q)(0);
			} else {
				op->bc_bot2_add_ri(n,"Lambda","Frad",2*PI*ones(1,1),map.leg.I_00,(r*r*sqrt(1+map.rt*map.rt/r/r)).row(j0));
				q=(Frad*r/sqrt(1+map.rt*map.rt/r/r)*(2+map.rt*map.rt/r/r)).row(j0);
				op->bc_bot2_add_d(n,"Lambda","eta",2*PI*ones(1,1)*(q*map.J[0].row(j0),map.leg.I_00));
				op->bc_bot2_add_d(n,"Lambda","deta",2*PI*ones(1,1)*(q*map.J[1].row(j0),map.leg.I_00));
				op->bc_bot2_add_ri(n,"Lambda","Ri",2*PI*ones(1,1),map.leg.I_00,q*map.J[2].row(j0));
				op->bc_bot2_add_ri(n,"Lambda","dRi",2*PI*ones(1,1),map.leg.I_00,q*map.J[3].row(j0));
				q=(Frad*map.rt/sqrt(1+map.rt*map.rt/r/r)).row(j0);
				op->bc_bot2_add_d(n,"Lambda","eta",2*PI*ones(1,1)*(q*(map.J[0].row(j0),Dt),map.leg.I_00));
				op->bc_bot2_add_d(n,"Lambda","deta",2*PI*ones(1,1)*(q*(map.J[1].row(j0),Dt),map.leg.I_00));
				q=-(q,Dt)-q*cos(th)/sin(th);
				op->bc_bot2_add_ri(n,"Lambda","Ri",2*PI*ones(1,1),map.leg.I_00,q*map.J[2].row(j0));
				op->bc_bot2_add_ri(n,"Lambda","dRi",2*PI*ones(1,1),map.leg.I_00,q*map.J[3].row(j0));
				op->bc_bot1_add_d(n,"Lambda","lum",-ones(1,1));
				rhs_Lambda(n)=-2*PI*((Frad*r*r*sqrt(1+map.rt*map.rt/r/r)).row(j0),map.leg.I_00)(0)+lum(n-1);
			}
		} else {
			op->bc_bot2_add_d(n,"Lambda","Lambda",ones(1,1));
			op->bc_bot1_add_d(n,"Lambda","Lambda",-ones(1,1));
		}
		
		j0+=map.gl.npts[n];
	}
	op->set_rhs("T",rhs_T);
	op->set_rhs("Lambda",rhs_Lambda);
	
}

void star2d::solve_dim(solver *op) {

	int n,j0;
	matrix q,rhs;
	
	rhs=zeros(ndomains(),1);
	j0=0;
	for(n=0;n<ndomains();n++) {
		op->bc_bot2_add_d(n,"m","m",ones(1,1));
		//rho
		op->bc_bot2_add_lri(n,"m","p",-2*PI*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,(r*r*map.rz*rho/eos.chi_rho).block(j0,j0+map.gl.npts[n]-1,0,nth()-1));
		op->bc_bot2_add_lri(n,"m","T",-2*PI*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,(-r*r*map.rz*rho*eos.d).block(j0,j0+map.gl.npts[n]-1,0,nth()-1));
		q=-2*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(r*r*map.rz*rho/eos.chi_rho).block(j0,j0+map.gl.npts[n]-1,0,nth()-1),map.leg.I_00);
		op->bc_bot2_add_d(n,"m","pc",q);
		q=-2*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(-r*r*map.rz*rho*eos.d).block(j0,j0+map.gl.npts[n]-1,0,nth()-1),map.leg.I_00);
		op->bc_bot2_add_d(n,"m","Tc",q);
		q=-2*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(-r*r*map.rz*rho).block(j0,j0+map.gl.npts[n]-1,0,nth()-1),map.leg.I_00);
		op->bc_bot2_add_d(n,"m","rhoc",q);
		//r (rz)
		q=-2*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(rho*(D,r*r*map.J[0])).block(j0,j0+map.gl.npts[n]-1,0,nth()-1),map.leg.I_00);
		op->bc_bot2_add_d(n,"m","eta",q);
		q=-2*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(rho*(D,r*r*map.J[1])).block(j0,j0+map.gl.npts[n]-1,0,nth()-1),map.leg.I_00);
		op->bc_bot2_add_d(n,"m","deta",q);
		q=(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(rho*(D,r*r*map.J[2])).block(j0,j0+map.gl.npts[n]-1,0,nth()-1));
		op->bc_bot2_add_ri(n,"m","Ri",-2*PI*ones(1,1),map.leg.I_00,q);
		q=(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(rho*(D,r*r*map.J[3])).block(j0,j0+map.gl.npts[n]-1,0,nth()-1));
		op->bc_bot2_add_ri(n,"m","dRi",-2*PI*ones(1,1),map.leg.I_00,q);
		
		if(n) op->bc_bot1_add_d(n,"m","m",-ones(1,1));
		j0+=map.gl.npts[n];
	}
	op->set_rhs("m",rhs);
	
	rhs=zeros(ndomains(),1);
	for(n=0;n<ndomains();n++) {
		if(n==ndomains()-1) {
			op->add_d(n,"rhoc","rhoc",ones(1,1));
			op->add_d(n,"rhoc","pc",-1./eos.chi_rho(0)*ones(1,1));
			op->add_d(n,"rhoc","Tc",eos.d(0)*ones(1,1));
		} else {
			op->bc_top1_add_d(n,"rhoc","rhoc",ones(1,1));
			op->bc_top2_add_d(n,"rhoc","rhoc",-ones(1,1));
		}
	}
	op->set_rhs("rhoc",rhs);
	
	rhs=zeros(ndomains(),1);
	for(n=0;n<ndomains();n++) {
		if(n==ndomains()-1) {
			op->add_d(n,"pc","pc",ones(1,1));
			op->add_d(n,"pc","pi_c",ones(1,1));
			op->add_d(n,"pc","rhoc",-2*ones(1,1));
			op->add_d(n,"pc","R",-2*ones(1,1));
		} else {
			op->bc_top1_add_d(n,"pc","pc",ones(1,1));
			op->bc_top2_add_d(n,"pc","pc",-ones(1,1));
		}
	}
	op->set_rhs("pc",rhs);
	
	rhs=zeros(ndomains(),1);
	for(n=0;n<ndomains();n++) {
		if(n==ndomains()-1) {
			op->add_d(n,"Tc","Tc",ones(1,1));
			op->add_d(n,"Tc","rhoc",-ones(1,1));
			op->add_d(n,"Tc","Lambda",ones(1,1));
			op->add_d(n,"Tc","R",-2*ones(1,1));
		} else {
			op->bc_top1_add_d(n,"Tc","Tc",ones(1,1));
			op->bc_top2_add_d(n,"Tc","Tc",-ones(1,1));
		}
	}
	op->set_rhs("Tc",rhs);
	
	rhs=zeros(ndomains(),1);
	for(n=0;n<ndomains();n++) {
		if(n==ndomains()-1) {
			op->add_d(n,"R","R",3*ones(1,1));
			op->add_d(n,"R","m",1/m*ones(1,1));
			op->add_d(n,"R","rhoc",ones(1,1));
		} else {
			op->bc_top1_add_d(n,"R","R",ones(1,1));
			op->bc_top2_add_d(n,"R","R",-ones(1,1));
		}
	}
	op->set_rhs("R",rhs);
}

void star2d::solve_map(solver *op) {

	int n,j0;
	matrix Ri,TT,q,rhs;

	rhs=zeros(ndomains()+1,1);
	
	j0=0;
	for(n=0;n<conv;n++) {
		if(!n || conv==ndomains()) op->bc_top1_add_d(n,"eta","eta",ones(1,1));
		else {
			op->bc_top1_add_d(n,"eta","eta",ones(1,1)/map.gl.xif[n]);
			op->bc_top2_add_d(n,"eta","eta",-ones(1,1)/map.gl.xif[n+1]);
		}
		j0+=map.gl.npts[n];
	}
	map.leg.eval_00(map.leg.th,0,TT);
	n=conv;
	if(!conv) {
		op->bc_bot2_add_d(n,"eta","eta",ones(1,1));
	} else if(conv<ndomains()) {
		op->bc_bot2_add_d(n,"eta","eta",ones(1,1));
		op->bc_bot2_add_r(n,"eta","Ri",-ones(1,1),TT);
	}
	
	for(n=conv+1;n<ndomains();n++) {
		op->bc_bot2_add_d(n,"eta","eta",ones(1,1)/(1-map.gl.xif[n]));
		op->bc_bot1_add_d(n,"eta","eta",-ones(1,1)/(1-map.gl.xif[n-1]));
	}
	
	op->add_d(ndomains(),"eta","eta",ones(1,1));
	
	op->set_rhs("eta",rhs);
	
	rhs=zeros(ndomains(),1);
	for(n=0;n<ndomains();n++) {
		op->bc_top1_add_d(n,"deta","deta",ones(1,1));
		op->bc_top1_add_d(n,"deta","eta",ones(1,1));
		op->bc_top2_add_d(n,"deta","eta",-ones(1,1));
	}
	op->set_rhs("deta",rhs);
	
	rhs=zeros(ndomains(),nth());
	for(n=0;n<ndomains();n++) {
		op->bc_top1_add_d(n,"dRi","dRi",ones(1,nth()));
		op->bc_top1_add_d(n,"dRi","Ri",ones(1,nth()));
		op->bc_top2_add_d(n,"dRi","Ri",-ones(1,nth()));
	}
	op->set_rhs("dRi",rhs);
	
	Ri=zeros(ndomains()+1,nth());
	Ri.setblock(1,ndomains(),0,nth()-1,map.R);
	rhs=zeros(ndomains()+1,nth());
	
	op->add_d(0,"Ri","Ri",ones(1,nth()));
//for(n=1;n<=ndomains();n++) op->add_d(n,"Ri","Ri",ones(1,nth()));op->set_rhs("Ri",rhs);return;
	map.leg.eval_00(map.leg.th,zeros(1,nth()),TT);
	q=zeros(1,nth());
	q(0,nth()-1)=1;
	j0=map.gl.npts[0];
	for(n=1;n<=ndomains();n++) {
		if(n!=conv) {
			op->bc_bot2_add_r(n,"Ri","Ri",q,TT);
			op->bc_bot2_add_d(n,"Ri","eta",-q);
			rhs.setrow(n,q*(-(Ri.row(n),TT)+map.gl.xif[n]));
		}
		if(n<ndomains()) {
			if(n!=conv) {
				op->bc_bot2_add_d(n,"Ri","p",(1-q)*p.row(j0));
				op->bc_bot2_add_ri(n,"Ri","p",q-1,TT,p.row(j0));
				rhs.setrow(n,rhs.row(n)
					+(1-q)*(-p.row(j0)+(p.row(j0),TT)));
			} else {
				op->bc_bot2_add_l(n,"Ri","T",ones(1,nth()),D.block(n).row(0));
				op->bc_bot2_add_d(n,"Ri","T",1./eos.cp.row(j0)*(D,eos.cp).row(j0));
				op->bc_bot2_add_d(n,"Ri","Tc",1./eos.cp.row(j0)*(D,eos.cp).row(j0));
				op->bc_bot2_add_l(n,"Ri","p",-eos.del_ad.row(j0),D.block(n).row(0));
				op->bc_bot2_add_d(n,"Ri","p",-(D,eos.cp*eos.del_ad).row(j0)/eos.cp.row(j0));
				op->bc_bot2_add_d(n,"Ri","pc",-(D,eos.cp*eos.del_ad).row(j0)/eos.cp.row(j0));
				rhs.setrow(n,-( ((D,log(T))-eos.del_ad*(D,log(p))).row(j0)) );
			}			
		}
		if(n==ndomains()) {
			// Isobar
			op->bc_bot1_add_d(n,"Ri","p",(1-q)*p.row(nr()-1));
			op->bc_bot1_add_ri(n,"Ri","p",q-1,TT,p.row(nr()-1));
			rhs.setrow(n,rhs.row(n)
				+(1-q)*(-p.row(nr()-1)+(p.row(nr()-1),TT)));
			// Photosphere
			/*op->bc_bot1_add_d(n,"Ri","p",p.row(nr()-1));
			op->bc_bot1_add_d(n,"Ri","ps",-ps);
			rhs.setrow(n,rhs.row(n)
				+(-p.row(nr()-1)+ps));*/
		}
		if(n<ndomains()) j0+=map.gl.npts[n];
	}
	op->set_rhs("Ri",rhs);
}

/*
void star2d::solve_Omega(solver *op) {

	int n;
	matrix rhs;

	rhs=zeros(ndomains()+1,1);
	for(n=0;n<ndomains();n++) {
		op->bc_top1_add_d(n,"Omega2","Omega2",ones(1,1));
		op->bc_top2_add_d(n,"Omega2","Omega2",-ones(1,1));
	}
	matrix TT;
	double r1,rz1,dphi1;
	r1=map.leg.eval_00(rex.row(0),PI/2,TT)(0);
	rz1=(map.ex.rz.row(0),TT)(0);
	dphi1=(Dex.row(0),phiex,TT)(0);
	n=ndomains();
	op->bc_bot2_add_d(n,"Omega2","Omega2",ones(1,1));
	op->bc_bot2_add_lr(n,"Omega2","phi",-ones(1,1)*Omega_bk*Omega_bk/r1/rz1,Dex.row(0),TT);
	op->bc_bot2_add_r(n,"Omega2","Ri",ones(1,1)*Omega_bk*Omega_bk/r1/r1/rz1*dphi1,TT);
	op->bc_bot2_add_r(n,"Omega2","Ri",(Dex.row(0),map.ex.J[2],TT)*Omega_bk*Omega_bk/r1/rz1/rz1*dphi1,TT);
	op->bc_bot2_add_r(n,"Omega2","eta",(Dex.row(0),map.ex.J[0],TT)*Omega_bk*Omega_bk/r1/rz1/rz1*dphi1,TT);
	rhs(n)=-Omega2+dphi1/r1/rz1*Omega_bk*Omega_bk;
	op->set_rhs("Omega2",rhs);
}
*/

void star2d::solve_Omega(solver *op) {

	int n;
	matrix rhs;

	rhs=zeros(ndomains()+1,1);
	for(n=0;n<ndomains();n++) {
		op->bc_top1_add_d(n,"Omega","Omega",ones(1,1));
		op->bc_top2_add_d(n,"Omega","Omega",-ones(1,1));
	}
	matrix TT;
	double Req;
	Req=map.leg.eval_00(rex.row(0),PI/2,TT)(0);
	n=ndomains();
	op->bc_bot1_add_d(n,"Omega","Omega",ones(1,1));
	op->bc_bot1_add_d(n,"Omega","m",-ones(1,1)*Omega_bk*sqrt(pi_c/Req/Req/Req/4./PI)/sqrt(m)/2.);
	op->bc_bot1_add_d(n,"Omega","pi_c",-ones(1,1)*Omega_bk*sqrt(pi_c*m/Req/Req/Req/4./PI)/2.);
	op->bc_bot2_add_r(n,"Omega","Ri",ones(1,1)*Omega_bk*sqrt(pi_c*m/Req/Req/Req/Req/Req/4./PI)*3./2.,TT);
	rhs(n)=-Omega+Omega_bk*sqrt(pi_c*m/Req/Req/Req/4./PI);
	op->set_rhs("Omega",rhs);

}

/*
void star2d::solve_Omega(solver *op) {

	int n;
	matrix rhs;

	rhs=zeros(ndomains()+1,1);
	for(n=0;n<ndomains();n++) {
		op->bc_top1_add_d(n,"Omega2","Omega2",ones(1,1));
		op->bc_top2_add_d(n,"Omega2","Omega2",-ones(1,1));
	}
	n=ndomains();
	op->bc_bot2_add_d(n,"Omega2","Omega2",ones(1,1));
	rhs(n)=-Omega2+Omega_bk*Omega_bk;
	op->set_rhs("Omega2",rhs);
}
*/




void star2d::solve_gsup(solver *op) {

	matrix q;
	
	op->bc_top1_add_d(ndomains()-1,"gsup","gsup",ones(1,nth()));

	op->bc_top1_add_d(ndomains()-1,"gsup","p",gsup_/eos.chi_rho.row(nr()-1));
	op->bc_top1_add_d(ndomains()-1,"gsup","T",-gsup_*eos.d.row(nr()-1));
	op->bc_top1_add_d(ndomains()-1,"gsup","pc",gsup_/eos.chi_rho.row(nr()-1));
	op->bc_top1_add_d(ndomains()-1,"gsup","Tc",-gsup_*eos.d.row(nr()-1));
	op->bc_top1_add_d(ndomains()-1,"gsup","rhoc",-gsup_);
	
	q=(sqrt(map.gzz)/rho).row(nr()-1);
	op->bc_top1_add_li(ndomains()-1,"gsup","p",q,
		D.block(ndomains()-1).row(map.gl.npts[ndomains()-1]-1),p.block(nr()-map.gl.npts[ndomains()-1],nr()-1,0,nth()-1));
	q=(map.gzt/sqrt(map.gzz)/rho).row(nr()-1);
	op->bc_top1_add_ri(ndomains()-1,"gsup","p",q,Dt,p.row(nr()-1));
	
	q=-1./rho*map.gzt/sqrt(map.gzz)*(r*map.gzt*(D,p)+(2./r-r*map.gzt/map.gzz)*(p,Dt));
	add_bc_top1_dr(op,ndomains()-1,"gsup",q.row(nr()-1));
	q=-1./rho*sqrt(map.gzz)/map.rz*(D,p);
	add_bc_top1_drz(op,ndomains()-1,"gsup",q.row(nr()-1));
	q=-1./rho/sqrt(map.gzz)/map.rz*(map.gzt*(D,p)-(-1./r/r+map.gzt*map.gzt/map.gzz)*(p,Dt));
	add_bc_top1_drt(op,ndomains()-1,"gsup",q.row(nr()-1));
	
	op->set_rhs("gsup",zeros(1,nth()));
		
		
		
}
#include<string.h>
void star2d::check_jacobian(solver *op,const char *eqn) {

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
		B.phi=B.phi+a*B.phi+ar*B.phi*random_matrix(nr(),nth());
		B.phiex=B.phiex+a*B.phiex+ar*B.phiex*random_matrix(nex(),nth());
		B.p=B.p+a*B.p+ar*B.p*random_matrix(nr(),nth());
		B.pc=B.pc+asc*B.pc;
		B.T=B.T+a*B.T+ar*B.T*random_matrix(nr(),nth());
		B.Tc=B.Tc+asc*B.Tc;
		B.w=B.w+a*B.w+ar*B.w*random_matrix(nr(),nth());
		B.Omega=B.Omega+asc*B.Omega;
		B.G=B.G+a*B.G+ar*B.G*random_matrix(nr(),nth());
		B.map.R=B.map.R+a*B.map.R*sin(th)*sin(th)+ar*B.map.R*random_matrix(ndomains(),nth());
		B.map.remap();
	}
	
	B.fill();
	
	i=op->get_id("phi");
	y[i]=zeros(nr()+nex(),nth());
	y[i].setblock(0,nr()-1,0,nth()-1,B.phi-phi);
	y[i].setblock(nr(),nr()+nex()-1,0,nth()-1,B.phiex-phiex);
	i=op->get_id("p");
	y[i]=log(B.p)-log(p);
	i=op->get_id("pi_c");
	y[i]=(log(B.pi_c)-log(pi_c))*ones(ndomains(),1);
	i=op->get_id("T");
	y[i]=log(B.T)-log(T);
	i=op->get_id("Lambda");
	y[i]=(log(B.Lambda)-log(Lambda))*ones(ndomains(),1);
	i=op->get_id("eta");
	y[i]=zeros(ndomains()+1,1);
	y[i].setblock(1,ndomains(),0,0,B.map.eta()-map.eta());
	j=i;
	i=op->get_id("deta");
	y[i]=y[j].block(1,ndomains(),0,0)-y[j].block(0,ndomains()-1,0,0);
	i=op->get_id("Ri");
	y[i]=zeros(ndomains()+1,nth());
	y[i].setblock(1,ndomains(),0,nth()-1,B.map.R-map.R);
	j=i;
	i=op->get_id("dRi");
	y[i]=y[j].block(1,ndomains(),0,nth()-1)-y[j].block(0,ndomains()-1,0,nth()-1);
	i=op->get_id("Omega");
	y[i]=(B.Omega-Omega)*ones(ndomains()+1,1);
	i=op->get_id("rhoc");
	y[i]=(log(B.rhoc)-log(rhoc))*ones(ndomains(),1);
	i=op->get_id("pc");
	y[i]=(log(B.pc)-log(pc))*ones(ndomains(),1);
	i=op->get_id("Tc");
	y[i]=(log(B.Tc)-log(Tc))*ones(ndomains(),1);
	i=op->get_id("R");
	y[i]=(log(B.R)-log(R))*ones(ndomains(),1);
	i=op->get_id("m");
	q=0;
	j0=0;
	y[i]=zeros(ndomains(),1);
	for(j=0;j<ndomains();j++) {
		q+=2*PI*(B.map.gl.I.block(0,0,j0,j0+map.gl.npts[j]-1),
			(B.rho*B.r*B.r*B.map.rz).block(j0,j0+map.gl.npts[j]-1,0,nth()-1),
			B.map.leg.I_00)(0)-
			2*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[j]-1),
			(rho*r*r*map.rz).block(j0,j0+map.gl.npts[j]-1,0,nth()-1),
			map.leg.I_00)(0);
		y[i](j)=q;
		j0+=map.gl.npts[j];
	}
	i=op->get_id("ps");
	y[i]=log(B.ps)-log(ps);
	i=op->get_id("Ts");
	y[i]=log(B.Ts)-log(Ts);
	i=op->get_id("lum");
	y[i]=zeros(ndomains(),1);
	j0=0;
	q=0;
	for(j=0;j<ndomains();j++) {
		q+=2*PI*(B.map.gl.I.block(0,0,j0,j0+map.gl.npts[j]-1),
			(B.rho*B.nuc.eps*B.r*B.r*B.map.rz).block(j0,j0+map.gl.npts[j]-1,0,nth()-1),B.map.leg.I_00)(0)-
			2*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[j]-1),
			(rho*nuc.eps*r*r*map.rz).block(j0,j0+map.gl.npts[j]-1,0,nth()-1),map.leg.I_00)(0);
		y[i](j)=q;
		j0+=map.gl.npts[j];
	}
	i=op->get_id("Frad");
	y[i]=zeros(ndomains()*2-1,nth());
	j0=0;
	for(j=0;j<ndomains();j++) {
		if(j) y[i].setrow(2*j-1,B.Frad.row(j0)-Frad.row(j0));
		y[i].setrow(2*j,B.Frad.row(j0+map.gl.npts[j]-1)-Frad.row(j0+map.gl.npts[j]-1));
		j0+=map.gl.npts[j];
	}
	i=op->get_id("gsup");
	y[i]=B.gsup_-gsup_;
	i=op->get_id("w");
	y[i]=B.w-w;
	i=op->get_id("G");
	y[i]=B.G-G;	
	
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


