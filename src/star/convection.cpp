#include"star.h"
#include<stdlib.h>


/*
#define KCONV_COMMON				\
	symbolic S;				\
	sym a_,b_;						\
	matrix kc,Ja,Jb;				\
	kconv_common(kc,Ja,Jb,S,a_,b_);	\

void star2d::kconv_common(matrix &kc,matrix &Ja,matrix &Jb,symbolic &S,sym &a_,sym &b_) const {

	double alpha=config.ddata;

	sym p_,T_,rho_,del_ad_;
	
	p_=S.regvar("p");
	T_=S.regvar("T");
	rho_=S.regvar("rho");
	del_ad_=S.regvar("del_ad");
	S.set_map(map);
	S.set_value("p",p);
	S.set_value("T",T);
	S.set_value("rho",rho);
	S.set_value("del_ad",eos.del_ad);
	
	a_=(grad(p_),grad(p_))/p_/p_;
	b_=(grad(p_),(grad(T_)/T_-del_ad_*grad(p_)/p_))/rho_;

	matrix a,b,fb,dfb;
	
	a=a_.eval();
	b=b_.eval();
	//fb=sqrt(b*(b>0));
	//dfb=0.5/sqrt(abs(b))*(b>0);
	double eps=1.*(14.-log10(surff))/nr;
	fb=1/sqrt(2)*sqrt(b+b/tanh(b/eps));
	dfb=fb*(0.5/b-1./(exp(2.*b/eps)-1)/eps);
	
	Ja=-alpha*alpha*sqrt(eos.d)/a/a*fb*R*sqrt(rhoc*pc);
	Jb=alpha*alpha*sqrt(eos.d)/a*dfb*R*sqrt(rhoc*pc);
	kc=alpha*alpha*sqrt(eos.d)/a*fb*R*sqrt(rhoc*pc);
	
	int j0=0;
	for(int n=0;n<conv+5;n++) {
		Ja.setblock(j0,j0+map.gl.npts[n]-1,0,-1,zeros(map.gl.npts[n],nth));
		Jb.setblock(j0,j0+map.gl.npts[n]-1,0,-1,zeros(map.gl.npts[n],nth));
		kc.setblock(j0,j0+map.gl.npts[n]-1,0,-1,zeros(map.gl.npts[n],nth));
		j0+=map.gl.npts[n];
	}
	Ja.setrow(j0,zeros(1,nth));
	Jb.setrow(j0,zeros(1,nth));
	kc.setrow(j0,zeros(1,nth));

}
*/
matrix star2d::kconv() const {

	double alpha=1;

	symbolic S;
	sym p_;
	sym Hp_2_;
	
	S.set_map(map);
	
	p_=S.regvar("p");
	Hp_2_=(grad(p_),grad(p_))/p_/p_; // 1/Hp**2
	
	sym cp_,rho_,xi_,s_,L2_,R_,pc_,rhoc_;
	sym Sigma_;
	
	cp_=S.regvar("cp");
	rho_=S.regvar("rho");
	xi_=S.regvar("xi");
	s_=S.regvar("s");
	L2_=S.regvar("L2");
	R_=S.regconst("R");
	pc_=S.regconst("pc");
	rhoc_=S.regconst("rhoc");
	
	Sigma_=2./81.*cp_*rho_/xi_/xi_*(grad(p_),grad(s_))*L2_*L2_*R_*R_*pc_*rhoc_;
	
	
	S.set_value("p",p);
	matrix L2;
	
	L2=alpha*alpha/Hp_2_.eval();
	
	{
		int j0=0;
		for(int n=0;n<conv;n++) {
			L2.setblock(j0,j0+map.gl.npts[n]-1,0,-1,zeros(map.gl.npts[n],nth));
			j0+=map.gl.npts[n];
		}
		L2.setrow(j0,zeros(1,nth));
	}
	
	S.set_value("cp",eos.cp);
	S.set_value("rho",rho);
	S.set_value("xi",opa.xi);
	S.set_value("s",entropy());
	S.set_value("L2",L2);
	S.set_value("R",R*ones(1,1));
	S.set_value("pc",pc*ones(1,1));
	S.set_value("rhoc",rhoc*ones(1,1));
	matrix Sigma;
	
	Sigma=Sigma_.eval();
	
	{
		int j0=0;
		for(int n=0;n<conv;n++) {
			Sigma.setblock(j0,j0+map.gl.npts[n]-1,0,-1,zeros(map.gl.npts[n],nth));
			j0+=map.gl.npts[n];
		}
		Sigma.setrow(j0,zeros(1,nth));
	}
	
	matrix Phi;
	
	Sigma=Sigma*(Sigma>0);
	Phi=9./8./(Sigma+1e-20)*pow(sqrt(1+Sigma)-1,3);
	
	matrix kc;
	kc=opa.xi*Phi/rho/eos.cp;
	
	static figure fig("/XSERVE");
	fig.subplot(2,1);
	fig.semilogx(T,kc);
	fig.semilogx(T,N2());
	
	
	return kc;
}

void star2d::add_kconv(solver *op,const char *eqn,const matrix &d) {

	matrix kc;
	kc=kconv();

	op->add_d(eqn,"opa.xi",kc/opa.xi*d);
	op->add_d(eqn,"rho",-kc/rho*d);

}

void star2d::add_dkconv_dz(solver *op,const char *eqn,const matrix &d) {
/*
	KCONV_COMMON
	
	a_.add(op,eqn,"p",d*(D,Ja));
	a_.add(op,eqn,"r",d*(D,Ja));
	
	b_.add(op,eqn,"T",d*(D,Jb));
	b_.add(op,eqn,"p",d*(D,Jb));
	b_.add(op,eqn,"r",d*(D,Jb));
	
	S.Dz(a_).add(op,eqn,"p",d*Ja);
	S.Dz(a_).add(op,eqn,"r",d*Ja);
	
	S.Dz(b_).add(op,eqn,"T",d*Jb);
	S.Dz(b_).add(op,eqn,"p",d*Jb);
	S.Dz(b_).add(op,eqn,"r",d*Jb);
	
	op->add_d(eqn,"log_R",d*(D,kc));
	op->add_d(eqn,"log_pc",d*0.5*(D,kc));
	op->add_d(eqn,"log_rhoc",d*0.5*(D,kc));
	*/
}

