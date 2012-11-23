#include"star.h"
#include<stdlib.h>
#include<string.h>

void star1d::fill() {

	Y=1.-X-Z;
	upd_Xr();

	eq_state();
	
	opacity();
	nuclear();
	
	m=4*PI*(map.gl.I,rho*r*r)(0);
	R=pow(M/m/rhoc,1./3.);

	pi_c=(4*PI*GRAV*rhoc*rhoc*R*R)/pc;
	Lambda=rhoc*R*R/Tc;
	
	calc_units();
	
	atmosphere();
	
	phiex=phi(-1)/rex;
	w=zeros(nr,1);G=zeros(nr,1);vr=zeros(nr,1);vt=zeros(nr,1);
	
	Omega=0;Omega_bk=0;Ekman=0;Omegac=0;

}


void star1d::upd_Xr() {

	int ic,n;
	
	Xr=X*ones(nr,1);
	if(!conv) {
		if(Xc!=1) printf("Warning: Non-homogeneus composition without core convection not implemented\n");
		Xc=1;
		return;
	}
	ic=0;
	for(n=0;n<conv;n++) ic+=map.gl.npts[n];
	Xr.setblock(0,ic-1,0,0,Xc*X*ones(ic,1));
	Xr.setblock(ic,nr-1,0,0,X*ones(nr-ic,1));
}


solver *star1d::init_solver(int nvar_add) {

	int nvar;
	solver *op;
	
	nvar=27;
	
	op=new solver();
	op->init(ndomains,nvar+nvar_add,"full");
	
	op->maxit_ref=10;op->use_cgs=0;op->maxit_cgs=20;
	op->rel_tol=1e-12;op->abs_tol=1e-20;
	register_variables(op);
	
	return op;
}

void star1d::register_variables(solver *op) {

	int i,var_nr[ndomains];
	
	for(i=0;i<ndomains;i++) 
		var_nr[i]=map.gl.npts[i];
	op->set_nr(var_nr);

	op->regvar("phi");
	op->regvar("log_p");
	op->regvar_dep("p");
	op->regvar("pi_c");
	op->regvar("log_T");
	op->regvar_dep("T");
	op->regvar("Lambda");
	op->regvar("Ri");
	op->regvar("dRi");
	op->regvar_dep("r");
	op->regvar_dep("rz");
	op->regvar_dep("log_rhoc");
	op->regvar("log_pc");
	op->regvar("log_Tc");
	op->regvar("log_R");
	op->regvar("m");
	op->regvar("ps");
	op->regvar("Ts");
	op->regvar("lum");
	op->regvar("Frad");
	op->regvar_dep("rho");
	op->regvar_dep("s");
	op->regvar("Teff");
	op->regvar("gsup");
	op->regvar_dep("opa.xi");
	op->regvar_dep("opa.k");
	op->regvar_dep("nuc.eps");

}

double star1d::solve(solver *op) {

	int info[5];
	matrix rho0;
	double err,err2;
	
	op->reset();
	solve_definitions(op);
	solve_poisson(op);
	solve_pressure(op);
	solve_temp(op);
	solve_atm(op);
	solve_dim(op);
	solve_map(op);
	solve_gsup(op);
	solve_Teff(op);
	
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
	
	double q,h;
		
	h=1;
	q=config.newton_dmax;
	
	matrix dphi,dp,dT,dpc,dTc,dRi;
	
	dphi=op->get_var("phi");
	err=max(abs(dphi/phi));

	dp=op->get_var("log_p");
	err2=max(abs(dp));err=err2>err?err2:err;
	while(exist(abs(h*dp)>q)) h/=2;

	dT=op->get_var("log_T");	
	err2=max(abs(dT));err=err2>err?err2:err;
	while(exist(abs(h*dT)>q)) h/=2;

	dpc=op->get_var("log_pc");	
	err2=fabs(dpc(0)/pc);err=err2>err?err2:err;
	while(fabs(h*dpc(0))>q*pc) h/=2;

	dTc=op->get_var("log_Tc");	
	err2=fabs(dTc(0));err=err2>err?err2:err;
	while(fabs(h*dTc(0))>q) h/=2;
	
	dRi=op->get_var("Ri");	
	update_map(h*dRi);
	
	phi+=h*dphi;
	p+=h*dp*p;
	T+=h*dT*T;
	pc*=exp(h*dpc(0));
	Tc*=exp(h*dTc(0));

	err2=max(abs(dRi));err=err2>err?err2:err;
	
	rho0=rho;

	fill();
	
	err2=max(abs(rho-rho0));err=err2>err?err2:err;
	
	return err;

}

void star1d::update_map(matrix dR) {

	if(ndomains==1) return;

	double h=1,dmax=0.1;config.newton_dmax;
	
	matrix R0;
	R0=map.R;
	dR.setblock(0,-2,0,0,dR.block(1,-1,0,0));
	dR(-1)=0;
	while(exist(abs(h*dR)>dmax*R0)) h/=2;
	map.R+=h*dR;
	while(map.remap()) {
		h/=2;
		map.R=R0+h*dR;
	}
}

void star1d::solve_definitions(solver *op) {

	op->add_d("rho","p",rho/eos.chi_rho/p);
	op->add_d("rho","log_T",-rho*eos.d);
	op->add_d("rho","log_pc",rho/eos.chi_rho);
	op->add_d("rho","log_Tc",-rho*eos.d);
	op->add_d("rho","log_rhoc",-rho);
	
	op->add_d("r","Ri",map.J[0]);
	op->add_d("r","dRi",map.J[1]);
	op->add_d("r","Ri",map.J[2]);
	op->add_d("r","dRi",map.J[3]);
	
	op->add_d("rz","Ri",(D,map.J[0]));
	op->add_d("rz","dRi",(D,map.J[1]));
	op->add_d("rz","Ri",(D,map.J[2]));
	op->add_d("rz","dRi",(D,map.J[3]));
	
	op->add_d("s","T",eos.cp/T);
	op->add_d("s","log_Tc",eos.cp);
	op->add_d("s","p",-eos.cp*eos.del_ad/p);
	op->add_d("s","log_pc",-eos.cp*eos.del_ad);
	
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
	
}

void star1d::solve_poisson(solver *op) {

	int n,j0;
	matrix rhs;

	op->add_l("phi","phi",ones(nr,1),(D,D));	
	op->add_l("phi","phi",2./r,D);
	op->add_d("phi","rho",-pi_c*ones(nr,1));
	op->add_d("phi","pi_c",-rho);
	
	op->add_d("phi","r",-2./r/r*(D,phi));
	op->add_d("phi","rz",-2.*(D,D,phi)-2./r*(D,phi));
	
	rhs=-(D,D,phi)-2/r*(D,phi)+rho*pi_c;

	j0=0;
	for(n=0;n<ndomains;n++) {
		if(!n) {
			op->bc_bot2_add_l(n,"phi","phi",ones(1,1),D.block(0).row(0));
			rhs(0)=-(D,phi)(0);
		} else {
			op->bc_bot2_add_l(n,"phi","phi",ones(1,1),D.block(n).row(0));
			op->bc_bot1_add_l(n,"phi","phi",-ones(1,1),D.block(n-1).row(-1));
			
			op->bc_bot2_add_d(n,"phi","rz",-(D,phi).row(j0));
			op->bc_bot1_add_d(n,"phi","rz",(D,phi).row(j0-1));
			
			rhs(j0)=-(D,phi)(j0)+(D,phi)(j0-1);
		}
		if(n==ndomains-1) {
			op->bc_top1_add_l(n,"phi","phi",ones(1,1),D.block(n).row(-1));
			op->bc_top1_add_d(n,"phi","phi",ones(1,1));
			op->bc_top1_add_d(n,"phi","rz",-(D,phi).row(-1));
			rhs(-1)=-phi(-1)-(D,phi)(-1);
		} else {
			op->bc_top1_add_d(n,"phi","phi",ones(1,1));
			op->bc_top2_add_d(n,"phi","phi",-ones(1,1));
			rhs(j0+map.gl.npts[n]-1)=-phi(j0+map.gl.npts[n]-1)+phi(j0+map.gl.npts[n]);
		}
		j0+=map.gl.npts[n];
	}
	op->set_rhs("phi",rhs);
}

void star1d::solve_pressure(solver *op) {

	int n,j0;
	matrix rhs_p,rhs_pi_c;
	char eqn[8];
	
	op->add_d("p","log_p",p);
	strcpy(eqn,"log_p");
	
	op->add_l(eqn,"p",ones(nr,1),D);	
	op->add_d(eqn,"rho",(D,phi));
	op->add_l(eqn,"phi",rho,D);
	
	rhs_p=-(D,p)-rho*(D,phi);
	rhs_pi_c=zeros(ndomains,1);

	j0=0;

	for(n=0;n<ndomains;n++) {
		op->bc_bot2_add_d(n,eqn,"p",ones(1,1));
		if(n>0) op->bc_bot1_add_d(n,eqn,"p",-ones(1,1));
		if(n==0) rhs_p(0)=1.-p(0);
		else rhs_p(j0)=-p(j0)+p(j0-1);
		if(n<ndomains-1) {
			op->bc_top1_add_d(n,"pi_c","pi_c",ones(1,1));
			op->bc_top2_add_d(n,"pi_c","pi_c",-ones(1,1));
		} else {
			op->bc_top1_add_d(n,"pi_c","p",ones(1,1));
			op->bc_top1_add_d(n,"pi_c","ps",-ones(1,1));
			rhs_pi_c(n)=-p(-1)+ps(0);
		}
		
		j0+=map.gl.npts[n];
	}
	op->set_rhs(eqn,rhs_p);
	op->set_rhs("pi_c",rhs_pi_c);
}


void star1d::solve_temp(solver *op) {

	int n,j0;
	matrix q;
	char eqn[8];
	
	op->add_d("T","log_T",T);
	strcpy(eqn,"log_T");
	
	//Luminosity

	matrix rhs_lum,lum;
	
	lum=zeros(ndomains,1);
	j0=0;
	for(n=0;n<ndomains;n++) {
		if(n) lum(n)=lum(n-1);
		lum(n)+=4*PI*Lambda*(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),
			(rho*nuc.eps*r*r).block(j0,j0+map.gl.npts[n]-1,0,0))(0);
		j0+=map.gl.npts[n];
	}

	rhs_lum=zeros(ndomains,1);
	j0=0;
	for(n=0;n<ndomains;n++) {
		op->bc_bot2_add_d(n,"lum","lum",ones(1,1));
		op->bc_bot2_add_li(n,"lum","rho",-4*PI*Lambda*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(r*r*nuc.eps).block(j0,j0+map.gl.npts[n]-1,0,0));
		op->bc_bot2_add_li(n,"lum","nuc.eps",-4*PI*Lambda*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(r*r*rho).block(j0,j0+map.gl.npts[n]-1,0,0));
		op->bc_bot2_add_d(n,"lum","Lambda",-4*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(rho*nuc.eps*r*r).block(j0,j0+map.gl.npts[n]-1,0,0)));
		//r (rz)
		op->bc_bot2_add_li(n,"lum","r",-4*PI*Lambda*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(2*r*rho*nuc.eps).block(j0,j0+map.gl.npts[n]-1,0,0));
		op->bc_bot2_add_li(n,"lum","rz",-4*PI*Lambda*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(r*r*rho*nuc.eps).block(j0,j0+map.gl.npts[n]-1,0,0));
			
		if(n) op->bc_bot1_add_d(n,"lum","lum",-ones(1,1));
		j0+=map.gl.npts[n];
	}
	op->set_rhs("lum",rhs_lum);
	
	//Frad
	
	matrix rhs_Frad,Frad;
	int j1;
	
	Frad=-opa.xi*(D,T);
	rhs_Frad=zeros(ndomains*2-1,1);
	j0=0;
	for(n=0;n<ndomains;n++) {
		j1=j0+map.gl.npts[n]-1;
		
		if(n) op->bc_bot2_add_d(n,"Frad","Frad",ones(1,1));
		op->bc_top1_add_d(n,"Frad","Frad",ones(1,1));
		
		if(n) op->bc_bot2_add_l(n,"Frad","T",opa.xi.row(j0),D.block(n).row(0));
		op->bc_top1_add_l(n,"Frad","T",opa.xi.row(j1),D.block(n).row(-1));
				
		if(n) op->bc_bot2_add_d(n,"Frad","opa.xi",(D,T).row(j0));
		op->bc_top1_add_d(n,"Frad","opa.xi",(D,T).row(j1));
		
		if(n) op->bc_bot2_add_d(n,"Frad","rz",Frad.row(j0));
		op->bc_top1_add_d(n,"Frad","rz",Frad.row(j1));
		
		j0=j1+1;
	}
	op->set_rhs("Frad",rhs_Frad);
		
	
	//Temperature
	
	matrix rhs_T,rhs_Lambda;
	matrix qconv,qrad;
	
	qrad=zeros(nr,1);
	qconv=qrad;
	j0=0;
	for(n=0;n<ndomains;n++) {
		if(n<conv) qconv.setblock(j0,j0+map.gl.npts[n]-1,0,0,ones(map.gl.npts[n],1));
		else qrad.setblock(j0,j0+map.gl.npts[n]-1,0,0,ones(map.gl.npts[n],1));
		j0+=map.gl.npts[n];
	}
	
	
	rhs_T=zeros(nr,1);

	// T
	
	op->add_l(eqn,"T",ones(nr,1),(D,D));	
	op->add_l(eqn,"T",2./r,D);
	rhs_T+=-qrad*((D,D,T)+2/r*(D,T));
	q=(D,log(opa.xi));
	op->add_l(eqn,"T",qrad*q,D);
	rhs_T+=-qrad*q*(D,T);
	
	// r
	op->add_d(eqn,"r",-2./r/r*(D,T));
	op->add_d(eqn,"rz",-2.*(D,D,T)-2./r*(D,T));
	q=-(D,log(opa.xi))*(D,T)*2.;
	op->add_d(eqn,"rz",qrad*q);
	
	//rho
	q=Lambda*nuc.eps/opa.xi;
	op->add_d(eqn,"rho",qrad*q);
		
	rhs_T+=-qrad*Lambda*rho*nuc.eps/opa.xi;

	//opa.xi
	q=(D,T);
	op->add_li(eqn,"opa.xi",qrad*q,D,1/opa.xi);
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
	
	j0=0;
	for(n=0;n<ndomains;n++) {
		if(!n) {
			op->bc_bot2_add_d(n,eqn,"T",ones(1,1));
			rhs_T(j0)=1-T(j0);
		} else {
			op->bc_bot2_add_d(n,eqn,"T",ones(1,1));
			op->bc_bot1_add_d(n,eqn,"T",-ones(1,1));
			rhs_T(j0)=-T(j0)+T(j0-1);
		}
		if(n>=conv) {
			if(n<ndomains-1) {
				op->bc_top1_add_d(n,eqn,"Frad",ones(1,1));
				op->bc_top2_add_d(n,eqn,"Frad",-ones(1,1));
				
				rhs_T(j0+map.gl.npts[n]-1)=-Frad(j0+map.gl.npts[n]-1)+Frad(j0+map.gl.npts[n]);
			} else {
				op->bc_top1_add_d(n,eqn,"T",ones(1,1));
				op->bc_top1_add_d(n,eqn,"Ts",-ones(1,1));
				rhs_T(-1)=Ts(0)-T(-1);
			}
		}
		
		if(n<conv) {
			op->bc_top1_add_d(n,"Lambda","Lambda",ones(1,1));
			op->bc_top2_add_d(n,"Lambda","Lambda",-ones(1,1));
		} else if(n==conv) {
			if(!n) {
				op->bc_bot2_add_l(n,"Lambda","T",ones(1,1),D.block(0).row(0));
				rhs_Lambda(0)=-(D,T)(0);
			} else {
				op->bc_bot2_add_d(n,"Lambda","Frad",4*PI*(r*r).row(j0));
				op->bc_bot2_add_d(n,"Lambda","r",4*PI*(Frad*2*r).row(j0));
				op->bc_bot2_add_d(n,"Lambda","rz",4*PI*(Frad*r*r).row(j0));
				op->bc_bot1_add_d(n,"Lambda","lum",-ones(1,1));
				rhs_Lambda(n)=-4*PI*Frad(j0)*(r*r)(j0)+lum(n-1);
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


void star1d::solve_dim(solver *op) {

	int n,j0;
	matrix q,rhs;
	
	rhs=zeros(ndomains,1);
	j0=0;
	for(n=0;n<ndomains;n++) {
		op->bc_bot2_add_d(n,"m","m",ones(1,1));
		//rho
		op->bc_bot2_add_li(n,"m","rho",-4*PI*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(r*r).block(j0,j0+map.gl.npts[n]-1,0,0));
		//r (rz)
		op->bc_bot2_add_li(n,"m","r",-4*PI*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(2*r*rho).block(j0,j0+map.gl.npts[n]-1,0,0));
		op->bc_bot2_add_li(n,"m","rz",-4*PI*ones(1,1),map.gl.I.block(0,0,j0,j0+map.gl.npts[n]-1),(r*r*rho).block(j0,j0+map.gl.npts[n]-1,0,0));
		
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

void star1d::solve_map(solver *op) {

	int n,j0;
	double *Ri;
	matrix rhs;
	
	Ri=map.gl.xif;
	
	rhs=zeros(ndomains,1);
	for(n=0;n<ndomains;n++) {
		op->bc_top1_add_d(n,"dRi","dRi",ones(1,1));
		op->bc_top1_add_d(n,"dRi","Ri",ones(1,1));
		if(n<ndomains-1) op->bc_top2_add_d(n,"dRi","Ri",-ones(1,1));
	}	
	op->set_rhs("dRi",rhs);
	
	rhs=zeros(ndomains,1);
	j0=0;
	for(n=0;n<conv;n++) {
		if(!n || conv==ndomains) op->bc_top1_add_d(n,"Ri","Ri",ones(1,1));
		else {
			op->bc_top1_add_d(n,"Ri","Ri",ones(1,1)/Ri[n]);
			op->bc_top2_add_d(n,"Ri","Ri",-ones(1,1)/Ri[n+1]);
		}
		j0+=map.gl.npts[n];
	}
	
	n=conv;
	if(!conv) {
		op->bc_bot2_add_d(n,"Ri","Ri",ones(1,1));
	} else if(conv<ndomains) {
		matrix ds;
		ds=eos.cp*((D,log(T))-eos.del_ad*(D,log(p)));
		op->bc_bot2_add_l(n,"Ri","s",ones(1,1),D.block(n).row(0));
		rhs(n)=-ds(j0);
	}
	for(n=conv+1;n<ndomains;n++) {
		op->bc_bot2_add_d(n,"Ri","Ri",ones(1,1)/(1-Ri[n]));
		op->bc_bot1_add_d(n,"Ri","Ri",-ones(1,1)/(1-Ri[n-1]));
	}
	op->set_rhs("Ri",rhs);
}

void star1d::solve_gsup(solver *op) {

	matrix q,g;
	int n=ndomains-1;
	
	g=gsup()*ones(1,1);
	
	op->bc_top1_add_d(n,"gsup","gsup",ones(1,1));
	op->bc_top1_add_d(n,"gsup","log_pc",-g);
	op->bc_top1_add_d(n,"gsup","log_rhoc",g);
	op->bc_top1_add_d(n,"gsup","log_R",g);

	q=-pc/R/rhoc*ones(1,1);
	op->bc_top1_add_l(n,"gsup","p",q,D.block(n).row(-1));
	
	q=(D,phi);
	op->bc_top1_add_d(n,"gsup","rz",pc/R/rhoc*q.row(-1));
	
	op->set_rhs("gsup",zeros(1,1));
		
		
		
}

void star1d::solve_Teff(solver *op) {

	matrix q,Te,F;
	int n=ndomains-1;
	
	Te=Teff()*ones(1,1);
	F=SIG_SB*pow(Te,4);
	
	op->bc_top1_add_d(n,"Teff","Teff",4*SIG_SB*pow(Te,3));
	op->bc_top1_add_d(n,"Teff","log_Tc",-F);
	op->bc_top1_add_d(n,"Teff","log_R",F);
	op->bc_top1_add_d(n,"Teff","opa.xi",-F/opa.xi.row(-1));

	q=opa.xi*Tc/R;
	op->bc_top1_add_l(n,"Teff","T",q.row(-1),D.block(n).row(-1));
	
	q=-(D,T)*opa.xi;
	op->bc_top1_add_d(n,"Teff","rz",Tc/R*q.row(-1));
	
	op->set_rhs("Teff",zeros(1,1));
		
}


void star1d::check_jacobian(solver *op,const char *eqn) {

	star1d B;
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
		B.phi=B.phi+a*B.phi+ar*B.phi*random_matrix(nr,1);
		B.p=B.p+a*B.p+ar*B.p*random_matrix(nr,1);
		B.pc=B.pc+asc*B.pc;
		B.T=B.T+a*B.T+ar*B.T*random_matrix(nr,1);
		B.Tc=B.Tc+asc*B.Tc;
		B.map.R=B.map.R+a*B.map.R+ar*B.map.R*random_matrix(ndomains,1);
		B.map.remap();
	}
	
	B.fill();
	
	i=op->get_id("rho");
	y[i]=zeros(nr,1);
	i=op->get_id("opa.xi");
	y[i]=zeros(nr,1);
	i=op->get_id("nuc.eps");
	y[i]=zeros(nr,1);
	i=op->get_id("r");
	y[i]=zeros(nr,1);
	i=op->get_id("rz");
	y[i]=zeros(nr,1);
	i=op->get_id("s");
	y[i]=zeros(nr,1);
	i=op->get_id("opa.k");
	y[i]=zeros(nr,1);
	
	
	i=op->get_id("phi");
	y[i]=zeros(nr,1);
	y[i]=B.phi-phi;
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
	i=op->get_id("Ri");
	y[i]=zeros(ndomains,1);
	y[i].setblock(1,ndomains-1,0,0,(B.map.R-map.R).block(0,ndomains-2,0,0));
	j=i;
	i=op->get_id("dRi");
	y[i]=zeros(ndomains,1);
	y[i].setblock(0,ndomains-2,0,0,y[j].block(1,ndomains-1,0,0)-y[j].block(0,ndomains-2,0,0));
	y[i](ndomains-1)=-y[j](ndomains-1);
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
		q+=4*PI*(B.map.gl.I.block(0,0,j0,j0+map.gl.npts[j]-1),
			(B.rho*B.r*B.r).block(j0,j0+map.gl.npts[j]-1,0,0))(0)-
			4*PI*(map.gl.I.block(0,0,j0,j0+map.gl.npts[j]-1),
			(rho*r*r).block(j0,j0+map.gl.npts[j]-1,0,0))(0);
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
		q+=4*PI*B.Lambda*(B.map.gl.I.block(0,0,j0,j0+map.gl.npts[j]-1),
			(B.rho*B.nuc.eps*B.r*B.r).block(j0,j0+map.gl.npts[j]-1,0,0))(0)-
			4*PI*Lambda*(map.gl.I.block(0,0,j0,j0+map.gl.npts[j]-1),
			(rho*nuc.eps*r*r).block(j0,j0+map.gl.npts[j]-1,0,0))(0);
		y[i](j)=q;
		j0+=map.gl.npts[j];
	}
	i=op->get_id("Frad");
	y[i]=zeros(ndomains*2-1,1);
	j0=0;
	matrix Frad,BFrad;
	Frad=-opa.xi*(D,T);
	BFrad=-B.opa.xi*(B.D,B.T);
	for(j=0;j<ndomains;j++) {	
		if(j) y[i](2*j-1)=BFrad(j0)-Frad(j0);
		y[i](2*j)=BFrad(j0+map.gl.npts[j]-1)-Frad(j0+map.gl.npts[j]-1);
		j0+=map.gl.npts[j];
	}
	i=op->get_id("gsup");
	y[i]=(B.gsup()-gsup())*ones(1,1);
	i=op->get_id("Teff");
	y[i]=(B.Teff()-Teff())*ones(1,1);
	
	B.solve(op);
	rhs=op->get_rhs(eqn);
	B=*this;
	B.solve(op);
	
	op->mult(y);
	drhs=rhs-op->get_rhs(eqn);
	drhs2=y[op->get_id(eqn)]+drhs;
	
	static figure fig("/XSERVE");
	
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




