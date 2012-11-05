#include"star.h"
#include<stdlib.h>

void star1d::fill() {

	Y=1.-X-Z;
	upd_Xr();

	eq_state();
	
	opacity();
	nuclear();
	
	m=4*PI*(gl.I,rho*r*r)(0);
	R=pow(M/m/rhoc,1./3.);

	pi_c=(4*PI*GRAV*rhoc*rhoc*R*R)/pc;
	Lambda=rhoc*R*R/Tc;
	
	calc_units();
	
	calc_Frad();
	
	atmosphere();

}

void star1d::calc_units() {

	units.phi=pc/rhoc;
	units.p=pc;
	units.rho=rhoc;
	units.T=Tc;
	units.r=R;
	units.F=pc/R/rhoc;
}


void star1d::calc_Frad() {

	Frad=-opa.xi/Lambda*(D,T);

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
	for(n=0;n<conv;n++) ic+=gl.npts[n];
	Xr.setblock(0,ic-1,0,0,Xc*X*ones(ic,1));
	Xr.setblock(ic,nr-1,0,0,X*ones(nr-ic,1));
}


solver *star1d::init_solver() {

	int nvar;
	solver *op;
	
	nvar=17;
	
	op=new solver();
	op->init(ndomains,nvar,"full");
	
	op->maxit_ref=10;op->use_cgs=0;op->maxit_cgs=20;
	op->rel_tol=1e-12;op->abs_tol=1e-20;
	register_variables(op);
	
	return op;
}

void star1d::register_variables(solver *op) {

	int i,var_nr[ndomains];
	
	for(i=0;i<ndomains;i++) 
		var_nr[i]=gl.npts[i];
	op->set_nr(var_nr);

	op->regvar("phi");
	op->regvar("p");
	op->regvar("pi_c");
	op->regvar("T");
	op->regvar("Lambda");
	op->regvar("Ri");
	op->regvar("dx");
	op->regvar_dep("rhoc");
	op->regvar("pc");
	op->regvar("Tc");
	op->regvar("R");
	op->regvar("m");
	op->regvar("ps");
	op->regvar("Ts");
	op->regvar("lum");
	op->regvar("Frad");
	op->regvar_dep("rho");

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

	dp=op->get_var("p");
	err2=max(abs(dp));err=err2>err?err2:err;
	while(exist(abs(h*dp)>q)) h/=2;

	dT=op->get_var("T");	
	err2=max(abs(dT));err=err2>err?err2:err;
	while(exist(abs(h*dT)>q)) h/=2;

	dpc=op->get_var("pc");	
	err2=fabs(dpc(0)/pc);err=err2>err?err2:err;
	while(fabs(h*dpc(0))>q*pc) h/=2;

	dTc=op->get_var("Tc");	
	err2=fabs(dTc(0));err=err2>err?err2:err;
	while(fabs(h*dTc(0))>q) h/=2;

	phi+=h*dphi;
	p+=h*dp*p;
	T+=h*dT*T;
	pc*=exp(h*dpc(0));
	Tc*=exp(h*dTc(0));

	dRi=op->get_var("Ri");	
	update_map(dRi);
	err2=max(abs(dRi));err=err2>err?err2:err;
	
	rho0=rho;

	fill();
	
	err2=max(abs(rho-rho0));err=err2>err?err2:err;
	
	return err;

}

void star1d::solve_definitions(solver *op) {

	op->add_d("rho","p",rho/eos.chi_rho);
	op->add_d("rho","T",-rho*eos.d);
	op->add_d("rho","pc",rho/eos.chi_rho);
	op->add_d("rho","Tc",-rho*eos.d);
	op->add_d("rho","rhoc",-rho);
	
}

void star1d::solve_poisson(solver *op) {

	int n,j0;
	matrix a_map,b_map,rhs;

	op->add_l("phi","phi",ones(nr,1),(D,D));	
	op->add_l("phi","phi",2./r,D);
	op->add_d("phi","rho",-pi_c*ones(nr,1));
	op->add_d("phi","pi_c",-rho*pi_c);
	
	rhs=-(D,D,phi)-2/r*(D,phi)+rho*pi_c;

	a_map=-2.*(D,D,phi)-2./r*(D,phi);
	b_map=-2./r/r*(D,phi);

	op->add_d("phi","Ri",b_map);

	j0=0;
	for(n=0;n<ndomains;n++) {
	
		op->add_d(n,"phi","dx",( a_map.block(j0,j0+gl.npts[n]-1,0,0)
			+b_map.block(j0,j0+gl.npts[n]-1,0,0)*(r.block(j0,j0+gl.npts[n]-1,0,0)-gl.xif[n]) )
			/(gl.xif[n+1]-gl.xif[n]));
	
		op->bc_bot2_add_l(n,"phi","phi",ones(1,1),D.block(n).row(0));
		if(n) op->bc_bot1_add_l(n,"phi","phi",-ones(1,1),D.block(n-1).row(gl.npts[n-1]-1));
		
		if(n==0) rhs(j0)=-(D,phi)(j0);
		else rhs(j0)=-(D,phi)(j0)+(D,phi)(j0-1);
		
		op->bc_bot2_add_d(n,"phi","dx",-(D,phi)(j0)/(gl.xif[n+1]-gl.xif[n])*ones(1,1));
		if(n) op->bc_bot1_add_d(n,"phi","dx",(D,phi)(j0-1)/(gl.xif[n]-gl.xif[n-1])*ones(1,1));
		
		if(n==ndomains-1) {
			op->bc_top1_add_l(n,"phi","phi",ones(1,1),D.block(n).row(gl.npts[n]-1));
			op->bc_top1_add_d(n,"phi","phi",ones(1,1));
		} else {
			op->bc_top1_add_d(n,"phi","phi",ones(1,1));
			op->bc_top2_add_d(n,"phi","phi",-ones(1,1));
		}
		
		if(n==ndomains-1) rhs(j0+gl.npts[n]-1)=-phi(j0+gl.npts[n]-1)-(D,phi)(j0+gl.npts[n]-1);
		else rhs(j0+gl.npts[n]-1)=-phi(j0+gl.npts[n]-1)+phi(j0+gl.npts[n]);
		
		if(n==ndomains-1) op->bc_top1_add_d(n,"phi","dx",-(D,phi)(j0+gl.npts[n]-1)/(gl.xif[n+1]-gl.xif[n])*ones(1,1));
		
		j0+=gl.npts[n];
	}
	op->set_rhs("phi",rhs);
}

void star1d::solve_pressure(solver *op) {

	int n,j0;
	matrix a_map,rhs_p,rhs_pi_c;
	
	op->add_li("p","p",ones(nr,1),D,p);	
	op->add_d("p","rho",(D,phi));
	op->add_l("p","phi",rho,D);
	
	rhs_p=-(D,p)-rho*(D,phi);
	rhs_pi_c=zeros(ndomains,1);

	a_map=-(D,p)-rho*(D,phi);

	j0=0;

	for(n=0;n<ndomains;n++) {
	
		op->add_d(n,"p","dx",a_map.block(j0,j0+gl.npts[n]-1,0,0)/(gl.xif[n+1]-gl.xif[n]));

		op->bc_bot2_add_d(n,"p","p",p(j0)*ones(1,1));
		if(n>0) op->bc_bot1_add_d(n,"p","p",-p(j0-1)*ones(1,1));
		if(n==0) rhs_p(0)=1.-p(0);
		else rhs_p(j0)=-p(j0)+p(j0-1);
		if(n<ndomains-1) {
			op->bc_top1_add_d(n,"pi_c","pi_c",ones(1,1));
			op->bc_top2_add_d(n,"pi_c","pi_c",-ones(1,1));
		} else {
			op->bc_top1_add_d(n,"pi_c","p",p(nr-1)*ones(1,1));
			op->bc_top1_add_d(n,"pi_c","ps",-ps*ones(1,1));
			rhs_pi_c(n)=-p(nr-1)+ps;
		}
		
		j0+=gl.npts[n];
	}
	op->set_rhs("p",rhs_p);
	op->set_rhs("pi_c",rhs_pi_c);
}

void star1d::solve_temp(solver *op) {

	int n,j0,j1;
	matrix q[10],a_map,b_map,lum,rhs_T,rhs_Lambda,rhs_lum,rhs_Frad;
	
	lum=zeros(ndomains,1);
	j0=0;
	for(n=0;n<ndomains;n++) {
		if(n) lum(n)=lum(n-1);
		lum(n)+=4*PI*(gl.I.block(0,0,j0,j0+gl.npts[n]-1),
			(rho*nuc.eps*r*r).block(j0,j0+gl.npts[n]-1,0,0))(0);
		j0+=gl.npts[n];
	}

	rhs_lum=zeros(ndomains,1);
	j0=0;
	for(n=0;n<ndomains;n++) {
		op->bc_bot2_add_d(n,"lum","lum",ones(1,1));
		op->bc_bot2_add_li(n,"lum","p",-4*PI*ones(1,1),gl.I.block(0,0,j0,j0+gl.npts[n]-1),(r*r*rho*nuc.eps/eos.chi_rho*(1.+nuc.dlneps_lnrho)).block(j0,j0+gl.npts[n]-1,0,0));
		op->bc_bot2_add_li(n,"lum","T",-4*PI*ones(1,1),gl.I.block(0,0,j0,j0+gl.npts[n]-1),(r*r*rho*nuc.eps*(nuc.dlneps_lnT-eos.d*(1.+nuc.dlneps_lnrho))).block(j0,j0+gl.npts[n]-1,0,0));
		op->bc_bot2_add_d(n,"lum","pc",-4*PI*(gl.I.block(0,0,j0,j0+gl.npts[n]-1),(r*r*rho*nuc.eps/eos.chi_rho*(1.+nuc.dlneps_lnrho)).block(j0,j0+gl.npts[n]-1,0,0)));
		op->bc_bot2_add_d(n,"lum","Tc",-4*PI*(gl.I.block(0,0,j0,j0+gl.npts[n]-1),(r*r*rho*nuc.eps*(nuc.dlneps_lnT-eos.d*(1.+nuc.dlneps_lnrho))).block(j0,j0+gl.npts[n]-1,0,0)));
		op->bc_bot2_add_d(n,"lum","rhoc",-4*PI*(gl.I.block(0,0,j0,j0+gl.npts[n]-1),(-r*r*rho*nuc.eps).block(j0,j0+gl.npts[n]-1,0,0)));
		op->bc_bot2_add_d(n,"lum","dx",-4*PI*(gl.I.block(0,0,j0,j0+gl.npts[n]-1), (rho*nuc.eps*r*(3*r-2*gl.xif[n])).block(j0,j0+gl.npts[n]-1,0,0))/(gl.xif[n+1]-gl.xif[n]));
		op->bc_bot2_add_d(n,"lum","Ri",-4*PI*(gl.I.block(0,0,j0,j0+gl.npts[n]-1), (2*rho*nuc.eps*r).block(j0,j0+gl.npts[n]-1,0,0)));
		if(n) op->bc_bot1_add_d(n,"lum","lum",-ones(1,1));
		j0+=gl.npts[n];
	}
	op->set_rhs("lum",rhs_lum);

	rhs_Frad=zeros(ndomains*2,1);
	q[1]=opa.dlnxi_lnrho/eos.chi_rho;
	q[2]=(opa.dlnxi_lnT-eos.d*opa.dlnxi_lnrho);
	j0=0;
	for(n=0;n<ndomains;n++) {
		j1=j0+gl.npts[n]-1;
		
		op->bc_bot2_add_d(n,"Frad","Frad",ones(1,1));
		op->bc_top1_add_d(n,"Frad","Frad",ones(1,1));
		
		q[0]=opa.xi/Lambda;
		if(n) {
			op->bc_bot2_add_li(n,"Frad","T",q[0](j0)*ones(1,1),D.block(n).row(0),T.block(j0,j1,0,0));
			op->bc_bot2_add_d(n,"Frad","dx",-q[0](j0)*(D,T)(j0)/(gl.xif[n+1]-gl.xif[n])*ones(1,1));
		}
		op->bc_top1_add_li(n,"Frad","T",q[0](j1)*ones(1,1),D.block(n).row(gl.npts[n]-1),T.block(j0,j1,0,0));
		op->bc_top1_add_d(n,"Frad","dx",-q[0](j1)*(D,T)(j1)/(gl.xif[n+1]-gl.xif[n])*ones(1,1));
		
		if(n) op->bc_bot2_add_d(n,"Frad","Lambda",Frad(j0)*ones(1,1));
		op->bc_top1_add_d(n,"Frad","Lambda",Frad(j1)*ones(1,1));
		
		if(n) op->bc_bot2_add_d(n,"Frad","p",-Frad(j0)*q[1](j0)*ones(1,1));
		if(n) op->bc_bot2_add_d(n,"Frad","pc",-Frad(j0)*q[1](j0)*ones(1,1));
		if(n) op->bc_bot2_add_d(n,"Frad","T",-Frad(j0)*q[2](j0)*ones(1,1));
		if(n) op->bc_bot2_add_d(n,"Frad","Tc",-Frad(j0)*q[2](j0)*ones(1,1));
		op->bc_top1_add_d(n,"Frad","p",-Frad(j1)*q[1](j1)*ones(1,1));
		op->bc_top1_add_d(n,"Frad","pc",-Frad(j1)*q[1](j1)*ones(1,1));
		op->bc_top1_add_d(n,"Frad","T",-Frad(j1)*q[2](j1)*ones(1,1));
		op->bc_top1_add_d(n,"Frad","Tc",-Frad(j1)*q[2](j1)*ones(1,1));
		
		j0=j1+1;
	}
	op->set_rhs("Frad",rhs_Frad);

	rhs_T=zeros(nr,1);

	q[0]=ones(nr,1);
	q[2]=-eos.del_ad;
	q[3]=-(D,eos.cp*eos.del_ad)/eos.cp;
	q[1]=(D,eos.cp)/eos.cp;
	q[4]=-(D,log(T))+eos.del_ad*(D,log(p));
	
	j0=0;
	for(n=0;n<conv;n++) {
		op->add_l(n,"T","T",q[0].block(j0,j0+gl.npts[n]-1,0,0),D.block(n));
		op->add_l(n,"T","p",q[2].block(j0,j0+gl.npts[n]-1,0,0),D.block(n));
		op->add_d(n,"T","T",q[1].block(j0,j0+gl.npts[n]-1,0,0));
		op->add_d(n,"T","p",q[3].block(j0,j0+gl.npts[n]-1,0,0));
		op->add_d(n,"T","Tc",q[1].block(j0,j0+gl.npts[n]-1,0,0));
		op->add_d(n,"T","pc",q[3].block(j0,j0+gl.npts[n]-1,0,0));
		rhs_T.setblock(j0,j0+gl.npts[n]-1,0,0,q[4].block(j0,j0+gl.npts[n]-1,0,0));
		op->add_d(n,"T","dx",q[4].block(j0,j0+gl.npts[n]-1,0,0)/(gl.xif[n+1]-gl.xif[n]));
		j0+=gl.npts[n];
	}

	q[0]=2./r+(D,log(opa.xi));
	q[9]=(D,T)*(opa.dlnxi_lnT-eos.d*opa.dlnxi_lnrho);
	q[1]=(D,(opa.dlnxi_lnT-eos.d*opa.dlnxi_lnrho))*(D,T)
		-Lambda*rho/opa.xi*nuc.eps*eos.d*(1+nuc.dlneps_lnrho-opa.dlnxi_lnrho)
		+Lambda*rho/opa.xi*nuc.eps*(nuc.dlneps_lnT-opa.dlnxi_lnT);
	q[6]=(D,(opa.dlnxi_lnT-eos.d*opa.dlnxi_lnrho))*(D,T)
		-Lambda*rho/opa.xi*nuc.eps*eos.d*(1+nuc.dlneps_lnrho-opa.dlnxi_lnrho)
		+Lambda*rho/opa.xi*nuc.eps*(nuc.dlneps_lnT-opa.dlnxi_lnT);
	q[2]=(D,T)*opa.dlnxi_lnrho/eos.chi_rho;
	q[3]=(D,opa.dlnxi_lnrho/eos.chi_rho)*(D,T)
		+Lambda*rho/opa.xi*nuc.eps/eos.chi_rho*(1+nuc.dlneps_lnrho-opa.dlnxi_lnrho);
	q[7]=(D,opa.dlnxi_lnrho/eos.chi_rho)*(D,T)
		+Lambda*rho/opa.xi*nuc.eps/eos.chi_rho*(1+nuc.dlneps_lnrho-opa.dlnxi_lnrho);
	q[4]=Lambda*rho*nuc.eps/opa.xi;
	q[5]=-(D,D,T)-(2/r+(D,log(opa.xi)))*(D,T)-Lambda*rho*nuc.eps/opa.xi;
	q[8]=-Lambda*rho/opa.xi*nuc.eps;
	
	a_map=-2*(D,D,T)-(2./r+2*(D,log(opa.xi)))*(D,T);
	b_map=-2./r/r*(D,T);
	
	for(n=conv;n<ndomains;n++) {
		op->add_li(n,"T","T",ones(gl.npts[n],1),(D.block(n),D.block(n)),T.block(j0,j0+gl.npts[n]-1,0,0));
		op->add_li(n,"T","T",q[0].block(j0,j0+gl.npts[n]-1,0,0),D.block(n),T.block(j0,j0+gl.npts[n]-1,0,0));
		op->add_l(n,"T","T",q[9].block(j0,j0+gl.npts[n]-1,0,0),D.block(n));
		op->add_d(n,"T","T",q[1].block(j0,j0+gl.npts[n]-1,0,0));
		op->add_d(n,"T","Tc",q[6].block(j0,j0+gl.npts[n]-1,0,0));
		op->add_l(n,"T","p",q[2].block(j0,j0+gl.npts[n]-1,0,0),D.block(n));
		op->add_d(n,"T","p",q[3].block(j0,j0+gl.npts[n]-1,0,0));
		op->add_d(n,"T","pc",q[7].block(j0,j0+gl.npts[n]-1,0,0));
		op->add_d(n,"T","Lambda",q[4].block(j0,j0+gl.npts[n]-1,0,0));
		op->add_d(n,"T","rhoc",q[8].block(j0,j0+gl.npts[n]-1,0,0));
		rhs_T.setblock(j0,j0+gl.npts[n]-1,0,0,q[5].block(j0,j0+gl.npts[n]-1,0,0));
		op->add_d(n,"T","Ri",b_map.block(j0,j0+gl.npts[n]-1,0,0));
		op->add_d(n,"T","dx",( a_map.block(j0,j0+gl.npts[n]-1,0,0)
			+b_map.block(j0,j0+gl.npts[n]-1,0,0)*(r.block(j0,j0+gl.npts[n]-1,0,0)-gl.xif[n]) )
			/(gl.xif[n+1]-gl.xif[n]));
		j0+=gl.npts[n];
	}
	
	rhs_Lambda=zeros(ndomains,1);

	j0=0;
	for(n=0;n<ndomains;n++) {
		
		if(!n) {
			op->bc_bot2_add_d(n,"T","T",ones(1,1)*T(j0));
			rhs_T(j0)=1-T(j0);
		} else {
			op->bc_bot2_add_d(n,"T","T",ones(1,1)*T(j0));
			op->bc_bot1_add_d(n,"T","T",-ones(1,1)*T(j0-1));
			rhs_T(j0)=-T(j0)+T(j0-1);
		}
		
		
		if(n>=conv) {
			if(n==ndomains-1) {
				op->bc_top1_add_d(n,"T","T",ones(1,1)*T(nr-1));
				op->bc_top1_add_d(n,"T","Ts",-ones(1,1)*Ts);
				rhs_T(j0+gl.npts[n]-1)=Ts-T(nr-1);
			} else {
				op->bc_top1_add_li(n,"T","T",ones(1,1),D.block(n).row(gl.npts[n]-1),T.block(j0,j0+gl.npts[n]-1,0,0));
				op->bc_top2_add_li(n,"T","T",-ones(1,1),D.block(n+1).row(0),T.block(j0+gl.npts[n],j0+gl.npts[n]+gl.npts[n+1]-1,0,0));
				op->bc_top1_add_d(n,"T","dx",-ones(1,1)*(D,T)(j0+gl.npts[n]-1)/(gl.xif[n+1]-gl.xif[n]));
				op->bc_top2_add_d(n,"T","dx",ones(1,1)*(D,T)(j0+gl.npts[n])/(gl.xif[n+2]-gl.xif[n+1]));
				rhs_T(j0+gl.npts[n]-1)=-(D,T)(j0+gl.npts[n]-1)+(D,T)(j0+gl.npts[n]);
				
				/*
				op->bc_top1_add_d(n,"T","Frad",ones(1,1));
				op->bc_top2_add_d(n,"T","Frad",-ones(1,1));
				
				rhs_T(j0+gl.npts[n]-1)=-Frad(j0+gl.npts[n]-1)+Frad(j0+gl.npts[n]);
			// No effect in result but better condition number
				q[0]=-(Frad*(opa.dlnxi_lnT-eos.d*opa.dlnxi_lnrho))(j0+gl.npts[n]-1)*ones(1,1);
				op->bc_top1_add_d(n,"T","T",q[0]);
				op->bc_top2_add_d(n,"T","T",-q[0]);
				rhs_T(j0+gl.npts[n]-1)+=-q[0](0)*T(j0+gl.npts[n]-1)+q[0](0)*T(j0+gl.npts[n]);
			///////////////////////////////////////////////////////
				*/
			}
		}
		
		if(n<conv) {
			op->bc_top1_add_d(n,"Lambda","Lambda",ones(1,1));
			op->bc_top2_add_d(n,"Lambda","Lambda",-ones(1,1));
		} else if(n==conv) {
			if(!n) {
				op->bc_bot2_add_li(n,"Lambda","T",ones(1,1),D.block(n).row(0),T.block(j0,j0+gl.npts[n]-1,0,0));
				op->bc_bot2_add_d(n,"Lambda","dx",-(D,T)(j0)/(gl.xif[n+1]-gl.xif[n])*ones(1,1));
				rhs_Lambda(n)=-(D,T)(j0);
			} else {
				op->bc_bot2_add_d(n,"Lambda","Frad",ones(1,1));
				op->bc_bot2_add_d(n,"Lambda","Ri",lum(n-1)/2./PI/gl.xif[n]/gl.xif[n]/gl.xif[n]*ones(1,1));
				op->bc_bot1_add_d(n,"Lambda","lum",-1./4./PI/gl.xif[n]/gl.xif[n]*ones(1,1));
				rhs_Lambda(n)=-Frad(j0)+lum(n-1)/4./PI/gl.xif[n]/gl.xif[n];
			}
		} else {
			op->bc_bot2_add_d(n,"Lambda","Lambda",ones(1,1));
			op->bc_bot1_add_d(n,"Lambda","Lambda",-ones(1,1));
		}
		
		j0+=gl.npts[n];
	}
	op->set_rhs("T",rhs_T);
	op->set_rhs("Lambda",rhs_Lambda);
}


void star1d::solve_dim(solver *op) {

	int n,j0;
	matrix q,rhs;
	
	rhs=zeros(ndomains,1);
	j0=0;
	for(n=0;n<ndomains;n++) {
		op->bc_bot2_add_d(n,"m","m",ones(1,1));
		op->bc_bot2_add_li(n,"m","rho",-4*PI*ones(1,1),gl.I.block(0,0,j0,j0+gl.npts[n]-1),(r*r).block(j0,j0+gl.npts[n]-1,0,0));
		op->bc_bot2_add_d(n,"m","dx",-4*PI*(gl.I.block(0,0,j0,j0+gl.npts[n]-1), (rho*r*(3*r-2*gl.xif[n])).block(j0,j0+gl.npts[n]-1,0,0))/(gl.xif[n+1]-gl.xif[n]));
		op->bc_bot2_add_d(n,"m","Ri",-4*PI*(gl.I.block(0,0,j0,j0+gl.npts[n]-1), (2*rho*r).block(j0,j0+gl.npts[n]-1,0,0)));
		if(n) op->bc_bot1_add_d(n,"m","m",-ones(1,1));
		j0+=gl.npts[n];
	}
	op->set_rhs("m",rhs);
	
	for(n=0;n<ndomains;n++) {
		op->add_d(n,"rhoc","pc",1./eos.chi_rho(0)*ones(1,1));
		op->add_d(n,"rhoc","Tc",-eos.d(0)*ones(1,1));
	}
	
	rhs=zeros(ndomains,1);
	for(n=0;n<ndomains;n++) {
		if(n==ndomains-1) {
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
	
	rhs=zeros(ndomains,1);
	for(n=0;n<ndomains;n++) {
		if(n==ndomains-1) {
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
	
	rhs=zeros(ndomains,1);
	for(n=0;n<ndomains;n++) {
		if(n==ndomains-1) {
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


void star1d::solve_map(solver *op) {

	int n,j0;
	double *Ri;
	matrix rhs;
	
	Ri=gl.xif;
	
	rhs=zeros(ndomains,1);
	for(n=0;n<ndomains;n++) {
		op->bc_top1_add_d(n,"dx","dx",ones(1,1));
		op->bc_top1_add_d(n,"dx","Ri",ones(1,1));
		if(n<ndomains-1) op->bc_top2_add_d(n,"dx","Ri",-ones(1,1));
	}	
	op->set_rhs("dx",rhs);
	
	rhs=zeros(ndomains,1);
	j0=0;
	for(n=0;n<conv;n++) {
		if(!n || conv==ndomains) op->bc_top1_add_d(n,"Ri","Ri",ones(1,1));
		else {
			op->bc_top1_add_d(n,"Ri","Ri",ones(1,1)/Ri[n]);
			op->bc_top2_add_d(n,"Ri","Ri",-ones(1,1)/Ri[n+1]);
		}
		j0+=gl.npts[n];
	}
	
	n=conv;
	if(!conv) {
		op->bc_bot2_add_d(n,"Ri","Ri",ones(1,1));
	} else if(conv<ndomains) {	
		op->bc_bot2_add_l(n,"Ri","T",ones(1,1),D.block(n).row(0));
		op->bc_bot2_add_d(n,"Ri","T",(D,eos.cp)(j0)/eos.cp(j0)*ones(1,1));
		op->bc_bot2_add_d(n,"Ri","Tc",(D,eos.cp)(j0)/eos.cp(j0)*ones(1,1));
		op->bc_bot2_add_l(n,"Ri","p",-eos.del_ad(j0)*ones(1,1),D.block(n).row(0));
		op->bc_bot2_add_d(n,"Ri","p",-(D,eos.del_ad*eos.cp)(j0)/eos.cp(j0)*ones(1,1));
		op->bc_bot2_add_d(n,"Ri","pc",-(D,eos.del_ad*eos.cp)(j0)/eos.cp(j0)*ones(1,1));
		rhs(n)=-(D,log(T))(j0)+eos.del_ad(j0)*(D,log(p))(j0);
		op->bc_bot2_add_d(n,"Ri","dx",rhs(n)/(gl.xif[n+1]-gl.xif[n])*ones(1,1));
	}
	for(n=conv+1;n<ndomains;n++) {
		op->bc_bot2_add_d(n,"Ri","Ri",ones(1,1)/(1-Ri[n]));
		op->bc_bot1_add_d(n,"Ri","Ri",-ones(1,1)/(1-Ri[n-1]));
	}
	op->set_rhs("Ri",rhs);
}

void star1d::update_map(matrix dR) {

	double Rc,Rc_new,dRc,h=1,q;
	int n;
	
	if(!conv) return;
	
	q=0.1;

	dRc=dR(conv);
	Rc=gl.xif[conv];
	while(fabs(h*dRc)>q*Rc) h/=2;
	Rc_new=Rc+h*dRc;
	
	gl.xif[conv]=Rc_new;
	for(n=1;n<conv;n++) gl.xif[n]=gl.xif[n]*Rc_new/Rc;
	for(n=conv+1;n<ndomains;n++) gl.xif[n]=1.-(1.-gl.xif[n])*(1.-Rc_new)/(1.-Rc);
	gl.xif[ndomains]=1;
	
	gl.init();

}
