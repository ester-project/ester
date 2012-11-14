#include"star.h"
#include<stdlib.h>
#include<sys/time.h>
#include<string.h>

void star2d::fill() {
	
	Y=1.-X-Z;
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

	atmosphere();

	//Omegac=sqrt(map.leg.eval_00(((Dex,phiex)/rex/map.ex.rz).row(0),PI/2)(0));
	Omegac=sqrt(pi_c*m/4/PI*(1-map.eps(ndomains-1))*(1-map.eps(ndomains-1))*(1-map.eps(ndomains-1)));

}



void star2d::upd_Xr() {

	int ic,n;
	
	Xr=X*ones(nr,nth);
	if(!conv) {
		if(Xc!=1) printf("Warning: Non-homogeneus composition without core convection not implemented\n");
		Xc=1;
		return;
	}
	ic=0;
	for(n=0;n<conv;n++) ic+=map.gl.npts[n];
	Xr.setblock(0,ic-1,0,-1,Xc*X*ones(ic,nth));
	Xr.setblock(ic,nr-1,0,-1,X*ones(nr-ic,nth));
	
}

void star2d::calc_veloc() {
// vr=rz*V^zeta vt=r*V^theta
	vr=(G,map.leg.D_11)/r+(map.rt/r+cos(th)/sin(th))/r*G;
	vr.setrow(0,zeros(1,nth));
	vr/=rho;
	vt=-(D,G)/map.rz-1./r*G;
	vt.setrow(0,zeros(1,nth));
	vt/=rho;
}

solver *star2d::init_solver() {

	int nvar;
	solver *op;
	
	nvar=33;
	op=new solver;
	op->init(ndomains+1,nvar,"full");
	
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

	op->regvar("phi");
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
	op->regvar_dep("lz");
	

}

double star2d::solve(solver *op) {

	int info[5];
	matrix rho0;
	double err,err2,h,dmax;

	if(Omega==0&Omega_bk!=0) {
		Omega=Omega_bk*Omegac;
		w=Omega*ones(nr,nth);
	}


	op->reset();
	
	solve_definitions(op);
	solve_poisson(op);
	solve_pressure(op);
	solve_temp(op);
	solve_dim(op);
	solve_map(op);
	solve_Omega(op);
	solve_atm(op);
	solve_gsup(op);
	solve_Teff(op);
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
	dphi=op->get_var("phi").block(0,nr-1,0,-1);
	dphiex=op->get_var("phi").block(nr,nr+nex-1,0,-1);
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

	double h=1,dmax=config.newton_dmax;
	
	matrix R0;
	R0=map.R;
	dR=dR.block(1,ndomains,0,-1);
	while(exist(abs(h*dR)>dmax*R0)) h/=2;
	map.R+=h*dR;
	while(map.remap()) {
		h/=2;
		map.R=R0+h*dR;
	}
}

void star2d::solve_definitions(solver *op) {

	op->add_d("rho","p",rho/eos.chi_rho/p);
	op->add_d("rho","T",-rho*eos.d/T);
	op->add_d("rho","log_pc",rho/eos.chi_rho);
	op->add_d("rho","log_Tc",-rho*eos.d);
	op->add_d("rho","log_rhoc",-rho);
	
	op->add_d("opa.xi","rho",opa.dlnxi_lnrho*opa.xi/rho);
	op->add_d("opa.xi","log_rhoc",opa.dlnxi_lnrho*opa.xi);
	op->add_d("opa.xi","T",opa.dlnxi_lnT*opa.xi/T);
	op->add_d("opa.xi","log_Tc",opa.dlnxi_lnT*opa.xi);
	
	op->add_d("opa.k","T",3*opa.k/T);
	op->add_d("opa.k","log_Tc",3*opa.k);
	op->add_d("opa.k","rho",-opa.k/rho);
	op->add_d("opa.k","log_rhoc",-opa.k);
	op->add_d("opa.k","opa.xi",-opa.k/opa.xi);
	
	op->add_d("nuc.eps","rho",nuc.dlneps_lnrho*nuc.eps/rho);
	op->add_d("nuc.eps","log_rhoc",nuc.dlneps_lnrho*nuc.eps);
	op->add_d("nuc.eps","T",nuc.dlneps_lnT*nuc.eps/T);
	op->add_d("nuc.eps","log_Tc",nuc.dlneps_lnT*nuc.eps);
	
	op->add_d("r","eta",map.J[0]);
	op->add_d("r","deta",map.J[1]);
	op->add_d("r","Ri",map.J[2]);
	op->add_d("r","dRi",map.J[3]);
	
	op->add_d(ndomains,"r","eta",map.ex.J[0]);
	op->add_d(ndomains,"r","Ri",map.ex.J[2]);
	
	op->add_d("rz","eta",(D,map.J[0]));
	op->add_d("rz","deta",(D,map.J[1]));
	op->add_d("rz","Ri",(D,map.J[2]));
	op->add_d("rz","dRi",(D,map.J[3]));
	
	op->add_d(ndomains,"rz","eta",(Dex,map.ex.J[0]));
	op->add_d(ndomains,"rz","Ri",(Dex,map.ex.J[2]));
	
	op->add_d("s","T",eos.cp/T);
	op->add_d("s","log_Tc",eos.cp);
	op->add_d("s","p",-eos.cp*eos.del_ad/p);
	op->add_d("s","log_pc",-eos.cp*eos.del_ad);
	
	op->add_d("lz","w",r*r*sin(th)*sin(th));
	op->add_d("lz","r",w*2*r*sin(th)*sin(th));
	
}

void star2d::solve_poisson(solver *op) {

	matrix q,rhs1,rhs2,rhs;
	int n,j0;
	matrix &gzz=map.gzz,&gzt=map.gzt,&gtt=map.gtt,
		&rz=map.rz,&rt=map.rt,&rzz=map.rzz,&rzt=map.rzt,&rtt=map.rtt;

	// phi
	map.add_lap(op,"phi","phi",ones(nr,nth),phi);
	
	//rho
	op->add_d("phi","rho",-pi_c*ones(nr,nth));

	//pi_c
	op->add_d("phi","pi_c",-rho);

	// phiex
	map.add_lap_ex(op,"phi","phi",ones(nex,nth),phiex);

	rhs1=map.lap(phi)-pi_c*rho;
	rhs2=map.lap_ex(phiex);
	rhs=zeros(nr+nex,nth);
	rhs.setblock(0,nr-1,0,-1,-rhs1);
	rhs.setblock(nr,nr+nex-1,0,-1,-rhs2);
	
	j0=0;
	for(n=0;n<ndomains+1;n++) {
		if(!n) {
			op->bc_bot2_add_l(n,"phi","phi",ones(1,nth),D.block(0).row(0));
			rhs.setrow(0,-(D,phi).row(0));
		} else {
			if(n<ndomains) op->bc_bot2_add_l(n,"phi","phi",1/rz.row(j0),D.block(n).row(0));
			else op->bc_bot2_add_l(n,"phi","phi",1/map.ex.rz.row(0),Dex.row(0));
			op->bc_bot1_add_l(n,"phi","phi",-1/rz.row(j0-1),D.block(n-1).row(-1));
			
			if(n<ndomains) 
				op->bc_bot2_add_d(n,"phi","rz",-1/rz.row(j0)/rz.row(j0)*(D,phi).row(j0));
			else
				op->bc_bot2_add_d(n,"phi","rz",-1/map.ex.rz.row(0)/map.ex.rz.row(0)*(Dex,phiex).row(0));
			op->bc_bot1_add_d(n,"phi","rz",1/rz.row(j0-1)/rz.row(j0-1)*(D,phi).row(j0-1));
			
			if(n<ndomains) rhs.setrow(j0,-(D,phi).row(j0)/rz.row(j0)+(D,phi).row(j0-1)/rz.row(j0-1));
			else rhs.setrow(j0,-(Dex,phiex).row(0)/map.ex.rz.row(0)+(D,phi).row(j0-1)/rz.row(j0-1));
		}
		
		op->bc_top1_add_d(n,"phi","phi",ones(1,nth));
		if(n<ndomains) op->bc_top2_add_d(n,"phi","phi",-ones(1,nth));
		if(n<ndomains) rhs.setrow(j0+map.gl.npts[n]-1,-phi.row(j0+map.gl.npts[n]-1));
		else rhs.setrow(nr+nex-1,-phiex.row(nex-1));
		if(n<ndomains-1) rhs.setrow(j0+map.gl.npts[n]-1,rhs.row(j0+map.gl.npts[n]-1)
								+phi.row(j0+map.gl.npts[n]));
		else if(n==ndomains-1) rhs.setrow(j0+map.gl.npts[n]-1,rhs.row(j0+map.gl.npts[n]-1)
								+phiex.row(0)); 
		
		if(n<ndomains) j0+=map.gl.npts[n];
	}
	op->set_rhs("phi",rhs);
}

void star2d::solve_pressure(solver *op) {

	matrix q,rhs_p,rhs_pi_c;
	int n,j0;
	matrix &gzz=map.gzz,&gzt=map.gzt,&gtt=map.gtt,
		&rz=map.rz,&rt=map.rt,&rzz=map.rzz,&rzt=map.rzt,&rtt=map.rtt;
	char eqn[8];
	
	op->add_d("p","log_p",p);
	strcpy(eqn,"log_p");
	
	// p	
	q=ones(nr,nth);
	op->add_l(eqn,"p",q,D);
	
	// phi
	q=rho;
	op->add_l(eqn,"phi",q,D);

	//rho
	q=(D,phi)-w*w*r*rz*sin(th)*sin(th);
	op->add_d(eqn,"rho",q);

	//w
	q=-rho*r*rz*sin(th)*sin(th);
	op->add_d(eqn,"w",2*w*q);

	// r
	q=-rho*w*w*rz*sin(th)*sin(th);
	op->add_d(eqn,"r",q);
	q=-rho*w*w*r*sin(th)*sin(th);
	op->add_l(eqn,"r",q,D);

	rhs_p=-(D,p)-rho*(D,phi)+rho*w*w*r*rz*sin(th)*sin(th);
	rhs_pi_c=zeros(ndomains,1);

	j0=0;
	for(n=0;n<ndomains;n++) {
		if(!n) {
			op->bc_bot2_add_d(n,eqn,"p",ones(1,nth));
			rhs_p.setrow(0,1-p.row(0));	
		} else {
			op->bc_bot2_add_d(n,eqn,"p",ones(1,nth));
			op->bc_bot1_add_d(n,eqn,"p",-ones(1,nth));
			rhs_p.setrow(j0,-p.row(j0)+p.row(j0-1));	
		}
		if(n<ndomains-1) {
			op->bc_top1_add_d(n,"pi_c","pi_c",ones(1,1));
			op->bc_top2_add_d(n,"pi_c","pi_c",-ones(1,1));
		} else {
			map.leg.eval_00(th,0,q);
			op->bc_top1_add_r(n,"pi_c","ps",-ones(1,1),q);
			op->bc_top1_add_r(n,"pi_c","p",ones(1,1),q);
			rhs_pi_c(ndomains-1)=(ps-p.row(-1),q)(0);
		}
		
		j0+=map.gl.npts[n];
	}
	op->set_rhs(eqn,rhs_p);
	op->set_rhs("pi_c",rhs_pi_c);
}


void star2d::solve_rot(solver *op) {

	matrix q,TT,rhs;
	int n,j0;
	matrix &gzz=map.gzz,&gzt=map.gzt,&gtt=map.gtt,
		&rz=map.rz,&rt=map.rt,&rzz=map.rzz,&rzt=map.rzt,&rtt=map.rtt;
	
	if(Omega==0) {
		op->add_d("w","w",ones(nr,nth));
		op->add_d("w","Omega",-ones(nr,nth));
		op->set_rhs("w",zeros(nr,nth));
		return;
	}
	
	// w
	
	q=(r*cos(th)+rt*sin(th))/rz;
	op->add_l("w","w",2*w*q,D);
	q=-sin(th)*ones(nr,nth);
	op->add_r("w","w",2*w*q,Dt);
	q=(r*cos(th)+rt*sin(th))/rz*2*(D,w)-sin(th)*2*(w,Dt);
	op->add_d("w","w",q);
	
	// p
	
	q=-(rho,Dt)/rho/rho/r/rz/sin(th);
	op->add_l("w","p",q,D);
	q=(D,rho)/rho/rho/r/rz/sin(th);
	op->add_r("w","p",q,Dt);

	//rho
	q=-2./rho/rho/rho/r/rz/sin(th)*((p,Dt)*(D,rho)-(D,p)*(rho,Dt));
	op->add_d("w","rho",q);
	q=1./rho/rho/r/rz/sin(th)*(p,Dt);
	op->add_l("w","rho",q,D);
	q=-1./rho/rho/r/rz/sin(th)*(D,p);
	op->add_r("w","rho",q,Dt);
	
	// r
	q=cos(th)/rz*2*w*(D,w)-1./rho/rho/r/r/rz/sin(th)*((p,Dt)*(D,rho)-(D,p)*(rho,Dt));
	op->add_d("w","r",q);
	q=-(r*cos(th)+rt*sin(th))/rz/rz*2*w*(D,w)-1./rho/rho/r/rz/rz/sin(th)*((p,Dt)*(D,rho)-(D,p)*(rho,Dt));
	op->add_l("w","r",q,D);
	q=sin(th)/rz*2*w*(D,w);
	op->add_r("w","r",q,Dt);

	rhs=-(r*cos(th)+rt*sin(th))/rz*2*w*(D,w)+
			sin(th)*2*w*(w,Dt)-1./rho/rho/r/rz/sin(th)*((p,Dt)*(D,rho)-(D,p)*(rho,Dt));
	
	j0=0;
	for(n=0;n<ndomains;n++) {
		if(!n) {
			op->bc_bot2_add_l(n,"w","w",ones(1,nth),D.block(n).row(0));
			rhs.setrow(0,-(D,w).row(0));
		}
		if(n<ndomains-1) {
			q=-rho*r*sin(th)*(r*cos(th)+rt*sin(th));
			op->bc_top1_add_d(n,"w","w",2.*(w*q).row(j0+map.gl.npts[n]-1));
			op->bc_top2_add_d(n,"w","w",-2.*(w*q).row(j0+map.gl.npts[n]));
			op->bc_top1_add_r(n,"w","phi",rho.row(j0+map.gl.npts[n]-1),Dt);
			op->bc_top2_add_r(n,"w","phi",-rho.row(j0+map.gl.npts[n]),Dt);
			q=-2.*rho*w*w*sin(th)*(r*cos(th)+rt*sin(th));
			op->bc_top1_add_d(n,"w","r",q.row(j0+map.gl.npts[n]-1));
			op->bc_top2_add_d(n,"w","r",-q.row(j0+map.gl.npts[n]));
			q=-rho*w*w*r*sin(th)*sin(th);
			op->bc_top1_add_r(n,"w","r",q.row(j0+map.gl.npts[n]-1),Dt);
			op->bc_top2_add_r(n,"w","r",-q.row(j0+map.gl.npts[n]),Dt);
			q=(phi,Dt)-w*w*r*sin(th)*(r*cos(th)+rt*sin(th));
			op->bc_top1_add_d(n,"w","rho",q.row(j0+map.gl.npts[n]-1));
			op->bc_top2_add_d(n,"w","rho",-q.row(j0+map.gl.npts[n]));
			q=rho*(phi,Dt)-rho*w*w*r*sin(th)*(r*cos(th)+rt*sin(th));
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
	matrix &gzz=map.gzz,&gzt=map.gzt,&gtt=map.gtt,
		&rz=map.rz,&rt=map.rt,&rzz=map.rzz,&rzt=map.rzt,&rtt=map.rtt;
	int limit_layer=1;
	
	if(!limit_layer) {
		qeq=ones(1,nth); // w constante
	} else {
		qeq=zeros(1,nth); //limit layer
	}
	qeq(0)=1;
	n=ndomains-1;
	
	matrix s;
	
	s=r*sin(th);
	
	// w
	
	q=s*s*(r*r+rt*rt)/r/rz;
	op->bc_top1_add_l(n,eqn,"w",(1-qeq)*q.row(-1),D.block(n).row(-1));
	q=-s*s*rt/r;
	op->bc_top1_add_r(n,eqn,"w",(1-qeq)*q.row(-1),Dt);
	op->bc_top1_add_r(n,eqn,"lz",(1-qeq)*G.row(-1),Dt);
	
	// G
	
	q=(s*s*w,Dt);
	op->bc_top1_add_d(n,eqn,"G",(1-qeq)*q.row(-1));
	
	// r
	
	q=(3*s*s+rt*rt*sin(th)*sin(th))/rz*(D,w)-rt*sin(th)*sin(th)*(w,Dt);
	op->bc_top1_add_d(n,eqn,"r",(1-qeq)*q.row(-1));
	q=-(r*s*s+r*rt*rt*sin(th)*sin(th))/rz/rz*(D,w);
	op->bc_top1_add_d(n,eqn,"rz",(1-qeq)*q.row(-1));
	q=2*r*rt*sin(th)*sin(th)/rz*(D,w)-r*sin(th)*sin(th)*(w,Dt);
	op->bc_top1_add_r(n,eqn,"r",(1-qeq)*q.row(-1),Dt);
	
	q=s*s*((r*r+rt*rt)/r/rz*(D,w)-rt/r*(w,Dt))+G*(s*s*w,Dt);
	rhs.setrow(-1,-(1-qeq)*q.row(-1));

	if(!limit_layer) {
		rhs.setrow(-1,-w.row(-1)+Omega);
		op->bc_top1_add_d(n,eqn,"w",qeq);
		op->bc_top1_add_d(n,eqn,"Omega",-qeq);
	} else {
		map.leg.eval_00(th,PI/2*ones(1,nth),TT);
		rhs.setrow(-1,rhs.row(-1)+qeq*(-(w.row(-1),TT)(0)+Omega));
		op->bc_top1_add_r(n,eqn,"w",qeq,TT);
		op->bc_top1_add_d(n,eqn,"Omega",-qeq);
	}	

}

void star2d::solve_dyn(solver *op) {

	matrix rhs;
	matrix q,s;
	int n,j0;
	matrix &gzz=map.gzz,&gzt=map.gzt,&gtt=map.gtt,
		&rz=map.rz,&rt=map.rt,&rzz=map.rzz,&rzt=map.rzt,&rtt=map.rtt;

	if(Omega==0) {
		op->add_d("G","G",ones(nr,nth));
		op->set_rhs("G",zeros(nr,nth));
		return;
	}
	
	s=r*sin(th);
	rhs=zeros(nr,nth);
	
	
	// w
	
	map.add_lap(op,"G","w",-s*s,w);
	rhs-=-s*s*map.lap(w);
	q=sin(th)/rz*(r*cos(th)+rt*sin(th))*G-2*sin(th)/rz*(r*sin(th)-rt*cos(th));
	op->add_l("G","w",q,D);
	rhs-=q*(D,w);
	q=(-sin(th)*sin(th)*G-2*sin(th)*cos(th));
	op->add_r("G","w",q,Dt);
	rhs-=q*(w,Dt);
	q=(G,map.leg.D_11)/r/rz;
	op->add_l("G","lz",q,D);
	rhs-=q*(D,s*s*w);
	q=-(D,G)/r/rz;
	op->add_r("G","lz",q,Dt);
	rhs-=q*(s*s*w,Dt);
	
	// G
	
	q=sin(th)/rz*(r*cos(th)+rt*sin(th))*(D,w)-sin(th)*sin(th)*(w,Dt);
	op->add_d("G","G",q);
	q=-(s*s*w,Dt)/r/rz;
	op->add_l("G","G",q,D);
	q=(D,s*s*w)/r/rz;
	op->add_r("G","G",q,map.leg.D_11);
	
	// r
	
	q=-((D,s*s*w)*(G,map.leg.D_11)-(s*s*w,Dt)*(D,G))/r/r/rz
		+sin(th)*cos(th)/rz*(D,w)*G-2*r*sin(th)*sin(th)*map.lap(w)-2*sin(th)*sin(th)/rz*(D,w);
	op->add_d("G","r",q);
	q=-((D,s*s*w)*(G,map.leg.D_11)-(s*s*w,Dt)*(D,G))/r/rz/rz
		-sin(th)/rz/rz*(r*cos(th)+rt*sin(th))*(D,w)*G+2*sin(th)/rz/rz*(r*sin(th)-rt*cos(th))*(D,w);
	op->add_d("G","rz",q);
	q=sin(th)*sin(th)/rz*(D,w)*G+2*sin(th)*cos(th)/rz*(D,w);
	op->add_r("G","r",q,Dt);
	
	j0=0;
	for(n=0;n<ndomains;n++) {
		if(!n) {
			op->bc_bot2_add_d(n,"G","G",ones(1,nth));
			rhs.setrow(0,-G.row(0));
		} else {
			op->bc_bot2_add_d(n,"G","G",ones(1,nth));
			op->bc_bot1_add_d(n,"G","G",-ones(1,nth));
			rhs.setrow(j0,-G.row(j0)+G.row(j0-1));
		}
	
		j0+=map.gl.npts[n];
	}
	
	
	op->set_rhs("G",rhs);

	
}



void star2d::solve_temp(solver *op) {

	int n,j0;
	matrix q;
	char eqn[8];
	matrix &gzz=map.gzz,&gzt=map.gzt,&gtt=map.gtt,
		&rz=map.rz,&rt=map.rt,&rzz=map.rzz,&rzt=map.rzt,&rtt=map.rtt;
	
	op->add_d("T","log_T",T);
	strcpy(eqn,"log_T");
	
	//Luminosity

	matrix rhs_lum,lum;
	
	lum=zeros(ndomains,1);
	j0=0;
	for(n=0;n<ndomains;n++) {
		if(n) lum(n)=lum(n-1);
		lum(n)+=2*PI*Lambda*(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),
			(rho*nuc.eps*r*r*rz).block(j0,j0+map.gl.npts[n]-1,0,-1),map.leg.I_00)(0);
		j0+=map.gl.npts[n];
	}

	rhs_lum=zeros(ndomains,1);
	j0=0;
	for(n=0;n<ndomains;n++) {
		op->bc_bot2_add_d(n,"lum","lum",ones(1,1));
		op->bc_bot2_add_lri(n,"lum","rho",-2*PI*Lambda*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,(r*r*rz*nuc.eps).block(j0,j0+map.gl.npts[n]-1,0,-1));
		op->bc_bot2_add_lri(n,"lum","nuc.eps",-2*PI*Lambda*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,(r*r*rz*rho).block(j0,j0+map.gl.npts[n]-1,0,-1));
		op->bc_bot2_add_d(n,"lum","Lambda",-2*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(rho*nuc.eps*r*r*rz).block(j0,j0+map.gl.npts[n]-1,0,-1),map.leg.I_00));
		//r (rz)
		op->bc_bot2_add_lri(n,"lum","r",-2*PI*Lambda*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,(2*r*rz*rho*nuc.eps).block(j0,j0+map.gl.npts[n]-1,0,-1));
		op->bc_bot2_add_lri(n,"lum","rz",-2*PI*Lambda*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),map.leg.I_00,(r*r*rho*nuc.eps).block(j0,j0+map.gl.npts[n]-1,0,-1));
			
		if(n) op->bc_bot1_add_d(n,"lum","lum",-ones(1,1));
		j0+=map.gl.npts[n];
	}
	op->set_rhs("lum",rhs_lum);
	
	//Frad
	
	matrix rhs_Frad,Frad;
	int j1;
	
	Frad=-opa.xi*(gzz*(D,T)+gzt*(T,Dt));
	rhs_Frad=zeros(ndomains*2-1,nth);
	j0=0;
	for(n=0;n<ndomains;n++) {
		j1=j0+map.gl.npts[n]-1;
		
		if(n) op->bc_bot2_add_d(n,"Frad","Frad",ones(1,nth));
		op->bc_top1_add_d(n,"Frad","Frad",ones(1,nth));
		
		q=opa.xi*gzz;
		if(n) op->bc_bot2_add_l(n,"Frad","T",q.row(j0),D.block(n).row(0));
		op->bc_top1_add_l(n,"Frad","T",q.row(j1),D.block(n).row(-1));
		q=opa.xi*gzt;
		if(n) op->bc_bot2_add_r(n,"Frad","T",q.row(j0),Dt);
		op->bc_top1_add_r(n,"Frad","T",q.row(j1),Dt);
				
		if(n) op->bc_bot2_add_d(n,"Frad","opa.xi",-Frad.row(j0)/opa.xi.row(j0));
		op->bc_top1_add_d(n,"Frad","opa.xi",-Frad.row(j1)/opa.xi.row(j1));
		
		q=opa.xi*(-2.*r*gzt*gzt*(D,T)-2./r*gzt*(T,Dt));
		if(n) op->bc_bot2_add_d(n,"Frad","r",q.row(j0));
		op->bc_top1_add_d(n,"Frad","r",q.row(j1));
		q=opa.xi*(-2./rz*gzz*(D,T)-1./rz*gzt*(T,Dt));
		if(n) op->bc_bot2_add_d(n,"Frad","rz",q.row(j0));
		op->bc_top1_add_d(n,"Frad","rz",q.row(j1));
		q=opa.xi*(-2./rz*gzt*(D,T)-1./r/r/rz*(T,Dt));
		if(n) op->bc_bot2_add_r(n,"Frad","r",q.row(j0),Dt);
		op->bc_top1_add_r(n,"Frad","r",q.row(j1),Dt);
		
		j0=j1+1;
	}
	op->set_rhs("Frad",rhs_Frad);
		
	
	//Temperature
	
	matrix rhs_T,rhs_Lambda;
	matrix TT,qconv,qrad;
	
	qrad=zeros(nr,nth);
	qconv=qrad;
	j0=0;
	for(n=0;n<ndomains;n++) {
		if(n<conv) qconv.setblock(j0,j0+map.gl.npts[n]-1,0,-1,ones(map.gl.npts[n],nth));
		else qrad.setblock(j0,j0+map.gl.npts[n]-1,0,-1,ones(map.gl.npts[n],nth));
		j0+=map.gl.npts[n];
	}
	
	
	rhs_T=zeros(nr,nth);

	// T
	map.add_lap(op,eqn,"T",qrad,T);
	rhs_T+=-qrad*map.lap(T);
	q=gzz*(D,log(opa.xi))+gzt*(log(opa.xi),Dt);
	op->add_l(eqn,"T",qrad*q,D);
	rhs_T+=-qrad*q*(D,T);
	q=gzt*(D,log(opa.xi))+gtt*(log(opa.xi),Dt);
	op->add_r(eqn,"T",qrad*q,Dt);
	rhs_T+=-qrad*q*(T,Dt);
	
	// r
	q=(D,log(opa.xi))*(D,T)*2*rt/r/rz*gzt
		-((D,log(opa.xi))*(T,Dt)+(log(opa.xi),Dt)*(D,T))*2*gzt/r
		-(log(opa.xi),Dt)*(T,Dt)*2/r/r/r;
	op->add_d(eqn,"r",qrad*q);
	q=-(D,log(opa.xi))*(D,T)*2/rz*gzz
		-((D,log(opa.xi))*(T,Dt)+(log(opa.xi),Dt)*(D,T))*gzt/rz;
	op->add_l(eqn,"r",qrad*q,D);
	q=-(D,log(opa.xi))*(D,T)*2/rz*gzt
		-((D,log(opa.xi))*(T,Dt)+(log(opa.xi),Dt)*(D,T))/r/r/rz;
	op->add_r(eqn,"r",qrad*q,Dt);
	
	//rho
	q=Lambda*nuc.eps/opa.xi;
	op->add_d(eqn,"rho",qrad*q);
		
	rhs_T+=-qrad*Lambda*rho*nuc.eps/opa.xi;

	//opa.xi
	q=gzz*(D,T)+gzt*(T,Dt);
	op->add_li(eqn,"opa.xi",qrad*q,D,1/opa.xi);
	q=gzt*(D,T)+gtt*(T,Dt);
	op->add_ri(eqn,"opa.xi",qrad*q,Dt,1/opa.xi);
	q=-Lambda*rho*nuc.eps/opa.xi/opa.xi;
	op->add_d(eqn,"opa.xi",qrad*q);
	
	//nuc.eps
	q=Lambda*rho/opa.xi;
	op->add_d(eqn,"nuc.eps",qrad*q);
	
	//Lambda
	q=rho*nuc.eps/opa.xi;
	op->add_d(eqn,"Lambda",qrad*q);
	
	op->add_l(eqn,"s",qconv,D);
	//rhs_T+=-qconv*(D,eos.s);
	rhs_T+=-qconv*eos.cp*((D,log(T))-eos.del_ad*(D,log(p)));
	
	rhs_Lambda=zeros(ndomains,1);
	
	map.leg.eval_00(th,0,TT);
	
	j0=0;
	for(n=0;n<ndomains;n++) {
		if(!n) {
			op->bc_bot2_add_d(n,eqn,"T",ones(1,nth));
			rhs_T.setrow(j0,1-T.row(j0));
		} else {
			op->bc_bot2_add_d(n,eqn,"T",ones(1,nth));
			op->bc_bot1_add_d(n,eqn,"T",-ones(1,nth));
			rhs_T.setrow(j0,-T.row(j0)+T.row(j0-1));
		}
		if(n>=conv) {
			if(n<ndomains-1) {
				op->bc_top1_add_d(n,eqn,"Frad",rz.row(j0+map.gl.npts[n]-1));
				op->bc_top2_add_d(n,eqn,"Frad",-rz.row(j0+map.gl.npts[n]-1));
				op->bc_top1_add_d(n,eqn,"rz",Frad.row(j0+map.gl.npts[n]-1));
				op->bc_top2_add_d(n,eqn,"rz",-Frad.row(j0+map.gl.npts[n]-1));
				
				rhs_T.setrow(j0+map.gl.npts[n]-1,
					-Frad.row(j0+map.gl.npts[n]-1)*rz.row(j0+map.gl.npts[n]-1)
					+Frad.row(j0+map.gl.npts[n])*rz.row(j0+map.gl.npts[n]));
			} else {
				op->bc_top1_add_d(n,eqn,"T",ones(1,nth));
				op->bc_top1_add_d(n,eqn,"Ts",-ones(1,nth));
				rhs_T.setrow(-1,Ts-T.row(-1));
			}
		}
		
		if(n<conv) {
			op->bc_top1_add_d(n,"Lambda","Lambda",ones(1,1));
			op->bc_top2_add_d(n,"Lambda","Lambda",-ones(1,1));
		} else if(n==conv) {
			if(!n) {
				map.leg.eval_00(th,PI/2,q);
				op->bc_bot2_add_lr(n,"Lambda","T",ones(1,1),D.block(0).row(0),q);
				rhs_Lambda(0)=-((D,T).row(0),q)(0);
			} else {
				op->bc_bot2_add_ri(n,"Lambda","Frad",2*PI*ones(1,1),map.leg.I_00,(r*r*rz).row(j0));
				op->bc_bot2_add_ri(n,"Lambda","r",2*PI*ones(1,1),map.leg.I_00,(Frad*2*r*rz).row(j0));
				op->bc_bot2_add_ri(n,"Lambda","rz",2*PI*ones(1,1),map.leg.I_00,(Frad*r*r).row(j0));
				op->bc_bot1_add_d(n,"Lambda","lum",-ones(1,1));
				rhs_Lambda(n)=-2*PI*(Frad.row(j0)*(r*r*rz).row(j0),map.leg.I_00)(0)+lum(n-1);
			}
		} else {
			op->bc_bot2_add_d(n,"Lambda","Lambda",ones(1,1));
			op->bc_bot1_add_d(n,"Lambda","Lambda",-ones(1,1));
		}
		
		j0+=map.gl.npts[n];
	}
	op->set_rhs(eqn,rhs_T);
	op->set_rhs("Lambda",rhs_Lambda);
	
}

void star2d::solve_dim(solver *op) {

	int n,j0;
	matrix q,rhs;
	
	rhs=zeros(ndomains,1);
	j0=0;
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
	
	for(n=0;n<ndomains;n++) {
		op->add_d(n,"log_rhoc","log_pc",1./eos.chi_rho(0)*ones(1,1));
		op->add_d(n,"log_rhoc","log_Tc",-eos.d(0)*ones(1,1));
	}

	
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

void star2d::solve_map(solver *op) {

	int n,j0;
	matrix Ri,TT,q,rhs;

	rhs=zeros(ndomains+1,1);
	
	j0=0;
	for(n=0;n<conv;n++) {
		if(!n || conv==ndomains) op->bc_top1_add_d(n,"eta","eta",ones(1,1));
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
	} else if(conv<ndomains) {
		op->bc_bot2_add_d(n,"eta","eta",ones(1,1));
		op->bc_bot2_add_r(n,"eta","Ri",-ones(1,1),TT);
	}
	
	for(n=conv+1;n<ndomains;n++) {
		op->bc_bot2_add_d(n,"eta","eta",ones(1,1)/(1-map.gl.xif[n]));
		op->bc_bot1_add_d(n,"eta","eta",-ones(1,1)/(1-map.gl.xif[n-1]));
	}
	
	op->add_d(ndomains,"eta","eta",ones(1,1));
	
	op->set_rhs("eta",rhs);
	
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
	
	Ri=zeros(ndomains+1,nth);
	Ri.setblock(1,ndomains,0,-1,map.R);
	rhs=zeros(ndomains+1,nth);
	
	op->add_d(0,"Ri","Ri",ones(1,nth));
//for(n=1;n<=ndomains;n++) op->add_d(n,"Ri","Ri",ones(1,nth));op->set_rhs("Ri",rhs);return;
	map.leg.eval_00(map.leg.th,zeros(1,nth),TT);
	q=zeros(1,nth);
	q(0,nth-1)=1;
	j0=map.gl.npts[0];
	for(n=1;n<=ndomains;n++) {
		if(n!=conv) {
			op->bc_bot2_add_r(n,"Ri","Ri",q,TT);
			op->bc_bot2_add_d(n,"Ri","eta",-q);
			rhs.setrow(n,q*(-(Ri.row(n),TT)+map.gl.xif[n]));
		}
		if(n<ndomains) {
			if(n!=conv) {
				op->bc_bot2_add_d(n,"Ri","p",(1-q));
				op->bc_bot2_add_r(n,"Ri","p",q-1,TT);
				rhs.setrow(n,rhs.row(n)
					+(1-q)*(-p.row(j0)+(p.row(j0),TT)));
			} else {
				matrix qq,drs,dts;
				matrix &gzz=map.gzz,&gzt=map.gzt,&gtt=map.gtt,
					&rz=map.rz,&rt=map.rt,&rzz=map.rzz,&rzt=map.rzt,&rtt=map.rtt;
			
				drs=eos.cp*((D,log(T))-eos.del_ad*(D,log(p)));
				dts=eos.cp*((log(T),Dt)-eos.del_ad*(log(p),Dt));
				//drs=(D,eos.s);
				//dts=(eos.s,Dt);
				
				qq=gzz*(D,p)+gzt*(p,Dt);
				op->bc_bot2_add_l(n,"Ri","s",qq.row(j0),D.block(n).row(0));
				qq=gzt*(D,p)+gtt*(p,Dt);
				op->bc_bot2_add_r(n,"Ri","s",qq.row(j0),Dt);
				qq=gzz*drs+gzt*dts;
				op->bc_bot2_add_l(n,"Ri","p",qq.row(j0),D.block(n).row(0));
				qq=gzt*drs+gtt*dts;
				op->bc_bot2_add_r(n,"Ri","p",qq.row(j0),Dt);
				
				qq=-2*r*gzt*gzt*(D,p)*drs-2/r*gzt*((D,p)*dts+(p,Dt)*drs)-2/r/r/r*(dts*(p,Dt));
				op->bc_bot2_add_d(n,"Ri","r",qq.row(j0));
				qq=-2/rz*gzz*(D,p)*drs-gzt/rz*((D,p)*dts+(p,Dt)*drs);
				op->bc_bot2_add_d(n,"Ri","rz",qq.row(j0));
				qq=-2/rz*gzt*(D,p)*drs-1/r/r/rz*((D,p)*dts+(p,Dt)*drs);
				op->bc_bot2_add_r(n,"Ri","r",qq.row(j0),Dt);
				
				qq=gzz*(D,p)*drs+gzt*((D,p)*dts+(p,Dt)*drs)+gtt*(dts*(p,Dt));
				rhs.setrow(n,-qq.row(j0));
			}			
		}
		if(n==ndomains) {
			// Isobar
			op->bc_bot1_add_d(n,"Ri","p",(1-q));
			op->bc_bot1_add_r(n,"Ri","p",q-1,TT);
			rhs.setrow(n,rhs.row(n)
				+(1-q)*(-p.row(-1)+(p.row(-1),TT)));
			// Photosphere
			/*op->bc_bot1_add_d(n,"Ri","p",ones(1,nth));
			op->bc_bot1_add_d(n,"Ri","ps",-ones(1,nth));
			rhs.setrow(n,rhs.row(n)
				+(-p.row(-1)+ps));*/
		}
		if(n<ndomains) j0+=map.gl.npts[n];
	}
	op->set_rhs("Ri",rhs);
}

/*
void star2d::solve_Omega(solver *op) {

	int n;
	matrix rhs;

	rhs=zeros(ndomains+1,1);
	for(n=0;n<ndomains;n++) {
		op->bc_top1_add_d(n,"Omega2","Omega2",ones(1,1));
		op->bc_top2_add_d(n,"Omega2","Omega2",-ones(1,1));
	}
	matrix TT;
	double r1,rz1,dphi1;
	r1=map.leg.eval_00(rex.row(0),PI/2,TT)(0);
	rz1=(map.ex.rz.row(0),TT)(0);
	dphi1=(Dex.row(0),phiex,TT)(0);
	n=ndomains;
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

/*
void star2d::solve_Omega(solver *op) {

	int n;
	matrix rhs;

	rhs=zeros(ndomains+1,1);
	for(n=0;n<ndomains;n++) {
		op->bc_top1_add_d(n,"Omega2","Omega2",ones(1,1));
		op->bc_top2_add_d(n,"Omega2","Omega2",-ones(1,1));
	}
	n=ndomains;
	op->bc_bot2_add_d(n,"Omega2","Omega2",ones(1,1));
	rhs(n)=-Omega2+Omega_bk*Omega_bk;
	op->set_rhs("Omega2",rhs);
}
*/




void star2d::solve_gsup(solver *op) {

	matrix q,g;
	int n=ndomains-1;
	matrix &gzz=map.gzz,&gzt=map.gzt,&rz=map.rz;
	
	g=gsup();
	
	op->bc_top1_add_d(n,"gsup","gsup",ones(1,nth));
	op->bc_top1_add_d(n,"gsup","log_pc",-g);
	op->bc_top1_add_d(n,"gsup","log_rhoc",g);
	op->bc_top1_add_d(n,"gsup","log_R",g);
	op->bc_top1_add_d(n,"gsup","rho",g/rho.row(-1));

	q=pc/R/rhoc/rho*sqrt(gzz);
	op->bc_top1_add_l(n,"gsup","p",q.row(-1),D.block(n).row(-1));
	q=pc/R/rhoc/rho*gzt/sqrt(gzz);
	op->bc_top1_add_r(n,"gsup","p",q.row(-1),Dt);
	
	q=(gzt/sqrt(gzz)*(-r*gzt*(D,p)+(r*gzt*gzt/gzz-2/r)*(p,Dt)))/rho;
	op->bc_top1_add_d(n,"gsup","r",pc/R/rhoc*q.row(-1));
	q=(-sqrt(gzz)/rz*(D,p))/rho;
	op->bc_top1_add_d(n,"gsup","rz",pc/R/rhoc*q.row(-1));
	q=(1/sqrt(gzz)*(-gzt/rz*(D,p)+(gzt*gzt/rz/gzz-1/r/r/rz)*(p,Dt)))/rho;
	op->bc_top1_add_r(n,"gsup","r",pc/R/rhoc*q.row(-1),Dt);
	
	op->set_rhs("gsup",zeros(1,nth));
		
		
		
}

void star2d::solve_Teff(solver *op) {

	matrix q,Te,F;
	int n=ndomains-1;
	matrix &gzz=map.gzz,&gzt=map.gzt,&rz=map.rz;
	
	Te=Teff();
	F=SIG_SB*pow(Te,4);
	
	op->bc_top1_add_d(n,"Teff","Teff",4*SIG_SB*pow(Te,3));
	op->bc_top1_add_d(n,"Teff","log_Tc",-F);
	op->bc_top1_add_d(n,"Teff","log_R",F);
	op->bc_top1_add_d(n,"Teff","opa.xi",-F/opa.xi.row(-1));

	q=opa.xi*Tc/R*sqrt(gzz);
	op->bc_top1_add_l(n,"Teff","T",q.row(-1),D.block(n).row(-1));
	q=opa.xi*Tc/R*gzt/sqrt(gzz);
	op->bc_top1_add_r(n,"Teff","T",q.row(-1),Dt);
	
	q=(gzt/sqrt(gzz)*(-r*gzt*(D,T)+(r*gzt*gzt/gzz-2/r)*(T,Dt)))*opa.xi;
	op->bc_top1_add_d(n,"Teff","r",Tc/R*q.row(-1));
	q=(-sqrt(gzz)/rz*(D,T))*opa.xi;
	op->bc_top1_add_d(n,"Teff","rz",Tc/R*q.row(-1));
	q=(1/sqrt(gzz)*(-gzt/rz*(D,T)+(gzt*gzt/rz/gzz-1/r/r/rz)*(T,Dt)))*opa.xi;
	op->bc_top1_add_r(n,"Teff","r",1/R/Tc*q.row(-1),Dt);
	
	op->set_rhs("Teff",zeros(1,nth));
		
		
		
}

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
	i=op->get_id("lz");
	y[i]=zeros(nr,nth);
	
	
	i=op->get_id("phi");
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
	y[i]=zeros(ndomains+1,1);
	y[i].setblock(1,ndomains,0,0,B.map.eta()-map.eta());
	j=i;
	i=op->get_id("deta");
	y[i]=y[j].block(1,ndomains,0,0)-y[j].block(0,ndomains-1,0,0);
	i=op->get_id("Ri");
	y[i]=zeros(ndomains+1,nth);
	y[i].setblock(1,ndomains,0,-1,B.map.R-map.R);
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


