#include "ester-config.h"
#include "star.h"
#include <stdlib.h>
#include <string.h>
#include "symbolic.h"

void star1d::fill() {
	Y0=1.-X0-Z0;
//        printf("star1d::fill call init_comp\n");
	init_comp();

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

solver *star1d::init_solver(int nvar_add) {
	int nvar;
	solver *op;
	
	nvar=29; // two more variables added (lnXh and Wr)
	
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

	op->regvar("Phi");
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
// Evolution of Xh
	op->regvar("lnXh");
	op->regvar("Wr");

}

double star1d::solve(solver *op) {
	int info[5];
	matrix rho_prec;
	double err,err2;
	
	//printf("************start of star1d::solve\n");
	check_map();
	
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
// Evolution of Xh
        solve_Xh(op);
        solve_Wr(op);
	
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

	double q,h;
		
	h=1;
	q=config.newton_dmax;
	
	matrix dphi,dp,dT,dXh,dpc,dTc,dRi,dWr;
	
FILE *fic=fopen("err.txt", "a");
	dphi=op->get_var("Phi");
	err=max(abs(dphi/phi));
  fprintf(fic,"err phi = %e\n",err);
 if (err> 1e-0) for (int k=0;k<nr;k++) fprintf(fic,"%d %e \n",k,dphi(k)/phi(k));

	dp=op->get_var("log_p");
	err2=max(abs(dp));err=err2>err?err2:err;
  fprintf(fic,"err P = %e\n",err2);
 if (err> 1e-0) for (int k=0;k<nr;k++) fprintf(fic,"%d %e \n",k,dp(k));
	while(exist(abs(h*dp)>q)) h/=2;

	dT=op->get_var("log_T");	
	err2=max(abs(dT));err=err2>err?err2:err;
  fprintf(fic,"err T = %e\n",err2);
 if (err> 1e-0) for (int k=0;k<nr;k++) fprintf(fic,"%d %e \n",k,dT(k));
	while(exist(abs(h*dT)>q)) h/=2;

// Compute dXh
	dXh=op->get_var("lnXh");	
	err2=max(abs(dXh));err=err2>err?err2:err;
  fprintf(fic,"err Xh = %e\n",err2);
 if (err> 1e-0) for (int k=0;k<nr;k++) fprintf(fic,"%d %e \n",k,dXh(k));
	while(exist(abs(h*dXh)>q)) h/=2;
// Compute dWr
	dWr=op->get_var("Wr");	
	err2=max(abs(dWr));err=err2>err?err2:err;
  fprintf(fic,"err Wr = %e\n",err2);
 if (err> 1e-4) for (int k=0;k<nr;k++) fprintf(fic,"%d %e %e \n",k,Wr(k),dWr(k));
	while(exist(abs(h*dWr)>q)) h/=2;
// End of dWr computation

	dpc=op->get_var("log_pc");	
	err2=fabs(dpc(0)/pc);err=err2>err?err2:err;
  fprintf(fic,"err Pc = %e\n",err2);
	while(fabs(h*dpc(0))>q*pc) h/=2;

	dTc=op->get_var("log_Tc");	
	err2=fabs(dTc(0));err=err2>err?err2:err;
	while(fabs(h*dTc(0))>q) h/=2;
	
	dRi=op->get_var("Ri");	
	update_map(h*dRi);
	
	phi+=h*dphi;
	p+=h*dp*p;
	T+=h*dT*T;
// Evolution of Xh, dXh is the variation on ln(Xh) assumed to be small
	Xh+=h*dXh*Xh;
	Wr+=h*dWr;

	pc*=exp(h*dpc(0));
	Tc*=exp(h*dTc(0));

	err2=max(abs(dRi));err=err2>err?err2:err;
	

	rho_prec=rho;

	fill();
	
	err2=max(abs(rho-rho_prec));err=err2>err?err2:err;
	
fclose(fic);
	//printf("************End of star1d::solve\n");
	return err;
}

void star1d::update_map(matrix dR) {
    if(ndomains==1) return;

    double h=1,dmax=config.newton_dmax;

    matrix R0;
    R0=map.R;
    dR=dR.concatenate(zeros(1,1));
    dR(0)=0;
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
	
	//	Valid only for homogeneus composition !!!!!!
	op->add_d("s","log_T",eos.cp);
	op->add_d("s","log_Tc",eos.cp);
	op->add_d("s","log_p",-eos.cp*eos.del_ad);
	op->add_d("s","log_pc",-eos.cp*eos.del_ad);
	
	op->add_d("opa.xi","rho",opa.dlnxi_lnrho*opa.xi/rho);
	op->add_d("opa.xi","log_rhoc",opa.dlnxi_lnrho*opa.xi);
	op->add_d("opa.xi","log_T",opa.dlnxi_lnT*opa.xi);
	op->add_d("opa.xi","log_Tc",opa.dlnxi_lnT*opa.xi);
	
	op->add_d("opa.k","log_T",3*opa.k);
	op->add_d("opa.k","log_Tc",3*opa.k);
	op->add_d("opa.k","rho",-opa.k/rho);
	op->add_d("opa.k","log_rhoc",-opa.k);
	op->add_d("opa.k","opa.xi",-opa.k/opa.xi);
	
	op->add_d("nuc.eps","rho",nuc.dlneps_lnrho*nuc.eps/rho);
	op->add_d("nuc.eps","log_rhoc",nuc.dlneps_lnrho*nuc.eps);
	op->add_d("nuc.eps","log_T",nuc.dlneps_lnT*nuc.eps);
	op->add_d("nuc.eps","log_Tc",nuc.dlneps_lnT*nuc.eps);
	
}

void star1d::solve_poisson(solver *op) {
    int n,j0;
    matrix rhs;

    // printf("r(0) = %e\n", r(0));

    matrix r0 = r;
    r0(0) = 1e-10;

    op->add_l("Phi","Phi", ones(nr,1),(D,D));
    op->add_l("Phi","Phi", 2./r0, D);
    op->add_d("Phi","rho", -pi_c*ones(nr,1));
    op->add_d("Phi","pi_c", -rho);

    op->add_d("Phi","r",-2./r0/r0*(D,phi));
    op->add_d("Phi","rz",-2.*(D,D,phi)-2./r0*(D,phi));

    rhs=-(D,D,phi)-2/r0*(D,phi)+rho*pi_c;

    j0=0;
    for(n=0;n<ndomains;n++) {
        if(n==0) {
            op->bc_bot2_add_l(n,"Phi","Phi",ones(1,1),D.block(0).row(0));
            rhs(0)=-(D,phi)(0);
        } else {
            op->bc_bot2_add_l(n,"Phi","Phi",ones(1,1),D.block(n).row(0));
            op->bc_bot1_add_l(n,"Phi","Phi",-ones(1,1),D.block(n-1).row(-1));

            op->bc_bot2_add_d(n,"Phi","rz",-(D,phi).row(j0));
            op->bc_bot1_add_d(n,"Phi","rz",(D,phi).row(j0-1));

            rhs(j0)=-(D,phi)(j0)+(D,phi)(j0-1);
        }
        if(n==ndomains-1) {
            op->bc_top1_add_l(n,"Phi","Phi",ones(1,1),D.block(n).row(-1));
            op->bc_top1_add_d(n,"Phi","Phi",ones(1,1));
            op->bc_top1_add_d(n,"Phi","rz",-(D,phi).row(-1));
            rhs(-1)=-phi(-1)-(D,phi)(-1);
        } else {
            op->bc_top1_add_d(n,"Phi","Phi",ones(1,1));
            op->bc_top2_add_d(n,"Phi","Phi",-ones(1,1));
            rhs(j0+map.gl.npts[n]-1)=-phi(j0+map.gl.npts[n]-1)+phi(j0+map.gl.npts[n]);
        }
        j0+=map.gl.npts[n];
    }
    op->set_rhs("Phi",rhs);
}

void star1d::solve_pressure(solver *op) {
    int n,j0;
    matrix rhs_p,rhs_pi_c;
    char eqn[8];

    op->add_d("p","log_p",p);
    strcpy(eqn,"log_p");

    op->add_l(eqn,"p",ones(nr,1),D);
    op->add_d(eqn,"rho",(D,phi));
    op->add_l(eqn,"Phi",rho,D);

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

//Evolution Xh --------------------------------------
void star1d::solve_Xh(solver *op) {
    DEBUG_FUNCNAME;

    double Qmc2=(4*HYDROGEN_MASS-AMASS["He4"]*UMA)*C_LIGHT*C_LIGHT;
    double factor=4*HYDROGEN_MASS/Qmc2*MYR;
    double L_core=0.;
    int n,j0,ndom,nc;
    matrix the_block;

// Core parameters
//    double rcc=Rcore()(0);
    double Ed=1e-6;
    nc=0;
    for (int i=0; i<conv; i++) {
        nc += map.gl.npts[i];
    }
    matrix &rz = map.rz;  // un alias
    if (nc) { // nc = nb of grid points in core
// compute Average X in core, and total eps in core
    X_core = ((map.gl.I.block(0, 0, 0, nc-1)), (Xh*rho*r*r*rz).block(0, nc-1, 0, 0))(0);
    M_core = ((map.gl.I.block(0, 0, 0, nc-1)), (rho*r*r*rz).block(0, nc-1, 0, 0))(0);
    X_core=X_core/M_core;
    L_core = ((map.gl.I.block(0, 0, 0, nc-1)), (rho*nuc.eps*r*r*rz).block(0, nc-1, 0, 0))(0);
    //printf("X_core =  %lf, X_core_prec = %lf, rcc/R = %lf\n", X_core, X_core_prec, rcc/R);
    }
// End core parameters


    j0=0; 
    matrix rhs;
    static symbolic S;
    static sym eq;
    static bool sym_inited = false;
    if (!sym_inited) {
// Do it only the first time the function is called
// If factor or Ed change, this should be recalculated
    sym lnX=S.regvar("lnXh");
    sym lnX0=S.regvar("lnXh0");
    sym X=exp(lnX);
    sym eps=S.regvar("nuc.eps");
    sym rho=S.regvar("rho"); //new
    sym rhovr=S.regvar("Wr"); //new
    sym r0=S.regvar("r0");
    sym r=S.r;//S.var("r"); // not regvar because "r" created automatically with the symbolic object
    sym drdt = (r-r0)/delta;
    sym lnR=S.regvar("log_R");
    sym lnR0=S.regvar("log_R0");
    sym dlnRdt = (lnR-lnR0)/delta;
    sym dlnXdt = (lnX-lnX0)/delta - (drdt+r*dlnRdt)/S.rz*Dz(lnX);
    eq=dlnXdt-Ed/X*lap(X)+factor*eps/X+rhovr/rho*Dz(lnX); //news
    sym_inited = true;
    }
    S.set_value("lnXh",log(Xh));
    S.set_value("lnXh0",log(Xh0));
    S.set_value("nuc.eps",nuc.eps);
    S.set_value("Wr",Wr); //new
    S.set_value("r0",r0); //new
    S.set_value("log_R",log(R)*ones(1,1)); //new
    S.set_value("log_R0",log(R0)*ones(1,1)); //new
    S.set_value("rho",rho); //new
    S.set_map(map);
    eq.add(op,"lnXh","lnXh");
    eq.add(op,"lnXh","nuc.eps");
    eq.add(op,"lnXh","r");
    eq.add(op,"lnXh","Wr"); //new
    eq.add(op,"lnXh","log_R"); //new
    rhs=-eq.eval();

// Take care of the CC domains
    matrix rhs_conv=-((log(Xh)-log(Xh0))/delta+factor*L_core/M_core/Xh);
    for(n=0;n<conv;n++) {
        ndom=map.gl.npts[n];
          op->reset(n,"lnXh");
          op->add_d(n,"lnXh","lnXh",ones(ndom,1)/delta);
          the_block=rhs_conv.block(j0,j0+ndom-1,0,0);
          rhs.setblock(j0,j0+ndom-1,0,0,the_block);
          j0+=ndom;
    }
// end of loop on CC domains

//BC: continuity of lnXh at the core boundary
         op->bc_bot2_add_d(conv,"lnXh","lnXh",ones(1,1));
         op->bc_bot1_add_d(conv,"lnXh","lnXh",-ones(1,1));
         rhs(j0)=-log(Xh(j0))+log(Xh(j0-1));
// at the stellar surface Xh=X0:
         op->bc_top1_add_d(ndomains-1,"lnXh","lnXh",ones(1,1));
         rhs(-1)=-log(Xh(-1))+log(X0);

//IC: continuity of lnXh and continuity of its derivative
	 for(n=conv;n<ndomains-1;n++){
           ndom=map.gl.npts[n];
           j0+=ndom;
           op->bc_top1_add_d(n,"lnXh","lnXh",ones(1,1));
           op->bc_top2_add_d(n,"lnXh","lnXh",-ones(1,1));
           rhs(j0-1) = -log(Xh(j0-1))+log(Xh(j0));
           op->bc_bot1_add_l(n+1,"lnXh","lnXh",ones(1,1),map.D.block(n).row(-1));
           op->bc_bot2_add_l(n+1,"lnXh","lnXh",-ones(1,1),map.D.block(n+1).row(0));
           rhs(j0) = -(map.D,log(Xh))(j0-1)+(map.D,log(Xh))(j0);
	 }


          op->set_rhs("lnXh",rhs);

//Evolution Xh end-----------------------------------
}

void star1d::solve_Wr(solver *op) {
    DEBUG_FUNCNAME;

// programmation en langage symbolique
    static symbolic S;
    static sym eq;
    static bool sym_inited = false;
    if (!sym_inited) { // Do it only the first time the function is called
    sym rho=S.regvar("rho");
    sym rho0=S.regvar("rho0");
    sym lnrhoc=S.regvar("log_rhoc");
    sym lnrhoc0=S.regvar("log_rhoc0");
// We use Wr=rho*Vr as the velocity variable. The use of Vr
// leads to huge Newton corrections and no convergence.
    sym Wr=S.regvar("Wr");
    sym lnR=S.regvar("log_R");
    sym lnR0=S.regvar("log_R0");
    sym r0=S.regvar("r0");
    sym r=S.r; // because "r" is created automatically with the symbolic object
    sym_vec rvec=grad(r);
    sym drdt = (r-r0)/delta;
    sym dlnRdt = (lnR-lnR0)/delta;
    sym drhodt = (rho-rho0)/delta - (drdt+r*dlnRdt)*Dz(rho)/S.rz;
    eq=drhodt+rho*(lnrhoc-lnrhoc0)/delta+div(Wr*rvec);
    sym_inited=true;
    }
    S.set_value("Wr",Wr);
    S.set_value("rho",rho);
    S.set_value("rho0",rho0);
    S.set_value("r0",r0);
    S.set_value("log_R",log(R)*ones(1,1));
    S.set_value("log_R0",log(R0)*ones(1,1));
    S.set_value("log_rhoc",log(rhoc)*ones(1,1));
    S.set_value("log_rhoc0",log(rhoc0)*ones(1,1));
    S.set_map(map);
    eq.add(op,"Wr","Wr");
    eq.add(op,"Wr","rho");
    eq.add(op,"Wr","r");
    eq.add(op,"Wr","log_rhoc");
    eq.add(op,"Wr","log_R");
    matrix rhs=-eq.eval();


//BC
         op->bc_bot2_add_d(0,"Wr","Wr",ones(1,1));
         rhs(0)=-Wr(0);

//IC
	 int j0=0;
         for(int n=1;n<ndomains;n++){
           int ndom=map.gl.npts[n-1];
           j0+=ndom;
           op->bc_bot1_add_d(n,"Wr","Wr",ones(1,1));
           op->bc_bot2_add_d(n,"Wr","Wr",-ones(1,1));
           rhs(j0) = -Wr(j0-1)+Wr(j0);
         }

	op->set_rhs("Wr",rhs);
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
    matrix qcore,qenv;

    qenv=zeros(nr,1);
    qcore=qenv;
    j0=0;
    for(n=0;n<ndomains;n++) {
        if(n<conv) qcore.setblock(j0,j0+map.gl.npts[n]-1,0,0,ones(map.gl.npts[n],1));
        else qenv.setblock(j0,j0+map.gl.npts[n]-1,0,0,ones(map.gl.npts[n],1));
        j0+=map.gl.npts[n];
    }


    rhs_T=zeros(nr,1);

    symbolic S;
    sym T_,xi_;
    sym div_Frad;

    S.set_map(map);

    T_=S.regvar("T");
    xi_=S.regvar("opa.xi");
    S.set_value("T",T);
    S.set_value("opa.xi",opa.xi);

    div_Frad=-div(-xi_*grad(T_))/xi_;

#define print_sym(s) { std::cout << #s << " = "; \
    (s).print(std::cout); std::cout << '\n'; }
    // print_sym(div(grad(T_)));
    div_Frad.add(op,eqn,"T",qenv);
    div_Frad.add(op,eqn,"opa.xi",qenv);
    div_Frad.add(op,eqn,"r",qenv);
    rhs_T-=div_Frad.eval()*qenv;


    op->add_d(eqn,"nuc.eps",qenv*Lambda*rho/opa.xi);
    op->add_d(eqn,"rho",qenv*Lambda*nuc.eps/opa.xi);
    op->add_d(eqn,"Lambda",qenv*rho*nuc.eps/opa.xi);
    op->add_d(eqn,"opa.xi",-qenv*Lambda*rho*nuc.eps/opa.xi/opa.xi);
    rhs_T+=-qenv*Lambda*rho*nuc.eps/opa.xi;


    //Convection

    /*if(env_convec) {
        sym kc_,rho_,s_;
        sym div_Fconv;

        kc_=S.regvar("kc");
        rho_=S.regvar("rho");
        s_=S.regvar("s");
        S.set_value("kc",kconv());
        S.set_value("rho",rho);
        S.set_value("s",entropy());

        div_Fconv=-div(-rho_*T_*grad(s_)*kc_)/xi_;

        div_Fconv.add(op,eqn,"T",qenv);
        div_Fconv.add(op,eqn,"rho",qenv);
        div_Fconv.add(op,eqn,"s",qenv);
        div_Fconv.add(op,eqn,"r",qenv);
        div_Fconv.add(op,eqn,"opa.xi",qenv);

        add_kconv(op,eqn,div_Fconv.jacobian("kc",0,0).eval()*qenv);
        add_dkconv_dz(op,eqn,div_Fconv.jacobian("kc",1,0).eval()*qenv);

        rhs_T-=div_Fconv.eval()*qenv;

    }*/

    //Core convection
    op->add_l(eqn,"s",qcore,D);
    //rhs_T+=-qcore*(D,eos.s);
    rhs_T+=-qcore*(D,entropy());

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
                op->bc_top1_add_l(n,eqn,"T",ones(1,1),D.block(n).row(-1));
                op->bc_top2_add_l(n,eqn,"T",-ones(1,1),D.block(n+1).row(0));
                op->bc_top1_add_d(n,eqn,"rz",-(D,T).row(j0+map.gl.npts[n]-1));
                op->bc_top2_add_d(n,eqn,"rz",(D,T).row(j0+map.gl.npts[n]));

                rhs_T(j0+map.gl.npts[n]-1)=-(D,T)(j0+map.gl.npts[n]-1)+(D,T)(j0+map.gl.npts[n]);
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
    matrix rhs;

    rhs=zeros(ndomains,1);
    for(n=0;n<ndomains;n++) {
        op->bc_top1_add_d(n,"dRi","dRi",ones(1,1));
        op->bc_top1_add_d(n,"dRi","Ri",ones(1,1));
        if(n<ndomains-1) op->bc_top2_add_d(n,"dRi","Ri",-ones(1,1));
    }
    op->set_rhs("dRi",rhs);

    rhs=zeros(ndomains,1);
    j0=0;
    for(n=0;n<ndomains;n++) {
        if(!n) op->add_d(n,"Ri","Ri",ones(1,1));
        else {
            matrix delta;
            delta=zeros(1,map.gl.npts[n]);delta(0)=1;delta(-1)=-1;
            op->bc_bot2_add_l(n,"Ri",LOG_PRES,ones(1,1),delta);
            delta=zeros(1,map.gl.npts[n-1]);delta(0)=1;delta(-1)=-1;
            op->bc_bot1_add_l(n,"Ri",LOG_PRES,-ones(1,1),delta);
            rhs(n)=log(PRES(j0+map.gl.npts[n]-1))-log(PRES(j0))-log(PRES(j0-1))+log(PRES(j0-map.gl.npts[n-1]));
        }
        j0+=map.gl.npts[n];
    }

    if(conv) {
        for(n=0,j0=0;n<conv;n++) j0+=map.gl.npts[n];
        n=conv;
        op->reset(n,"Ri");

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

        rhs(n)=-eq.eval()(j0);
    }
    op->set_rhs("Ri",rhs);
}

void star1d::solve_gsup(solver *op) {
    matrix q,g;
    int n=ndomains-1;

    g=gsup()*ones(1,1);
    /*
    op->bc_top1_add_d(n,"gsup","gsup",ones(1,1));
    op->bc_top1_add_d(n,"gsup","log_pc",-g);
    op->bc_top1_add_d(n,"gsup","log_rhoc",g);
    op->bc_top1_add_d(n,"gsup","log_R",g);

    q=-pc/R/rhoc*ones(1,1);
    op->bc_top1_add_l(n,"gsup","Phi",q,D.block(n).row(-1));

    q=(D,phi);
    op->bc_top1_add_d(n,"gsup","rz",pc/R/rhoc*q.row(-1));
    */
    op->bc_top1_add_d(n,"gsup","gsup",1./g);
    op->bc_top1_add_d(n,"gsup","log_rhoc",-ones(1,1));
    op->bc_top1_add_d(n,"gsup","log_R",-ones(1,1));
    op->bc_top1_add_d(n,"gsup","m",-ones(1,1)/m);

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

/*    op->bc_top1_add_d(n,"Teff","Teff",4./Te);
    op->bc_top1_add_d(n,"Teff","log_Tc",-ones(1,1));
    op->bc_top1_add_d(n,"Teff","log_R",ones(1,1));
    op->bc_top1_add_d(n,"Teff","lum",-ones(1,1)/luminosity()*R*Tc);*/


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


    i=op->get_id("Phi");
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

matrix calcMass(star1d &a) {

        solver massSolver;
        massSolver.init(a.map.ndomains, 1, "full");
        massSolver.regvar("mass");
        massSolver.set_nr(a.map.npts);

        massSolver.add_l("mass", "mass", ones(a.nr, 1), a.D);
        matrix rhs = 4 * PI * a.r * a.r * a.rho;

        massSolver.bc_bot2_add_d(0, "mass", "mass", ones(1, 1));
        rhs(0) = 0;
        int jfirst = 0;
        for(int i = 1; i < a.ndomains; i++) {
                jfirst += a.map.gl.npts[i-1];
                massSolver.bc_bot2_add_d(i, "mass", "mass", ones(1, 1));
                massSolver.bc_bot1_add_d(i, "mass", "mass", -ones(1,
1));
                rhs(jfirst) = 0;
        }
        massSolver.set_rhs("mass", rhs);
        massSolver.solve();

        matrix mass = massSolver.get_var("mass");

        return mass;

}



