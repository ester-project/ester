#include "ester-config.h"
#include "symbolic.h"
#include <cstdlib>
#include <cmath>
#include <typeinfo>

double symbolic::tol=1e-12;
bool symbolic::axisymmetric=true;
bool symbolic::simplify_auto=false;
bool symbolic::spherical=false;
bool symbolic::debug=false;

symbolic::symbolic() : one_(sym(1)), one(one_) {
	one_.context=this;
	r=regvar("r");
	zeta=regvar_indep("zeta");
	theta=regvar_indep("theta");
	phi=regvar_indep("phi");
	
	x[0]=zeta;
	x[1]=theta;
	x[2]=phi;
	
	if(spherical) {
		rz=1*one;
		rt=0*one;
		rzz=0*one;
		rzt=0*one;
		rtt=0*one;
	} else {
		rz=regvar("rz");
		rt=regvar("rt");
		rzz=regvar("rzz");
		rzt=regvar("rzt");
		rtt=regvar("rtt");
	}
	
	g(0,0)=1/Dz(r)/Dz(r)+Dt(r)*Dt(r)/r/r/Dz(r)/Dz(r);
	g(0,1)=-Dt(r)/r/r/Dz(r);
	g(1,0)=g(0,1);
	g(1,1)=1/r/r;
	g(2,2)=1/r/r/sin(theta)/sin(theta);
	g(0,2)=g(2,0)=g(1,2)=g(2,1)=0*one;
	g.set_type(CONTRAVARIANT,CONTRAVARIANT);
	
	g_(0,0)=Dz(r)*Dz(r);
	g_(0,1)=Dt(r)*Dz(r);
	g_(1,0)=g_(0,1);
	g_(1,1)=r*r+Dt(r)*Dt(r);
	g_(2,2)=r*r*sin(theta)*sin(theta);
	g_(0,2)=g_(2,0)=g_(1,2)=g_(2,1)=0*one;
	g_.set_type(COVARIANT,COVARIANT);
	
	sqrt_g=r*r*Dz(r)*sin(theta);
	
	maxder=0;
}

void symbolic::init() {
	
	sqrt_g=1/sqrt(det(g));
	g_=inv(g);

}

double symbolic::round_to_tol(double value) {

	if(tol==0) return value;
	//if(round(value)==0) return value;
// Round to nearest integer if the difference is less than symbolic::tol
	if(fabs(value-round(value))<tol) return round(value);
// Check also small fractions
	for(int i=2;i<=10;i++) 
		if(fabs(i*value-round(i*value))/i<tol) 
			return round(i*value)/i;
	return value;

}

sym symbolic::regvar(const std::string &name) {

	if(name=="") {
		std::cerr<<"Symbolic: Invalid name"<<std::endl;
		exit(1);
	}
	
	if(vars.count(name)) {
		std::cerr<<"Symbolic: Can't register variable \""<<name<<"\" (already registered)"<<std::endl;
		exit(1);
	}
	
	vars[name]=sym::symbol();
	vars[name].name=name;
	vars[name].context=this;
	vars[name].is_const=false;
	vars[name].is_indep=false;

	return var(name);

}

sym symbolic::regconst(const std::string &name) {

	regvar(name);
	vars[name].is_const=true;
	return var(name);

}

sym symbolic::regvar_indep(const std::string &name) {

	regvar(name);
	vars[name].is_indep=true;
	return var(name);

}


sym symbolic::var(const std::string &name) {

	if(!vars.count(name)) {
		std::cerr<<"Symbolic: Unknown variable \""<<name<<"\"\n";
		exit(1);
	}

	sym s;
	s.context=this;
	delete s.expr;
	s.expr=vars.find(name)->second.clone();
	
	return s;
}

void symbolic::set_value(const char *name,const matrix &value,int parity) {

	if(!vars.count(name)) {
		std::cerr<<"Symbolic: Unknown variable \""<<name<<"\"\n";
		exit(1);
	}
	
	val[name]=value;
	par[name]=parity;

}

void symbolic::set_map(const mapping &map_) {

	map=map_;
	val["zeta"]=map.z*ones(map.nr,map.nt);
	par["zeta"]=00;
	val["theta"]=map.th*ones(map.nr,map.nt);
	par["theta"]=-1;
	val["r"]=map.r;
	par["r"]=00;
	val["rz"]=map.rz;
	par["rz"]=00;
	val["rt"]=map.rt;
	par["rt"]=11;
	val["rzz"]=map.rzz;
	par["rzz"]=00;
	val["rzt"]=map.rzt;
	par["rzt"]=11;
	val["rtt"]=map.rtt;
	par["rtt"]=00;

}

sym::sym_expr *symbolic::derive(const sym::sym_expr &f,const sym::symbol &x) {
// To be called by sym_deriv::derive() and symbol::derive(), to set maxder and return special symbols for some
// derivatives: rz,rt,rzz,rzt,rtt
// If no special symbol needed returns NULL

	if(f==vars["r"]) {
		if(x==vars["zeta"]) {
			if(spherical) return new sym::sym_num(1);
			return vars["rz"].clone();
		}
		if(x==vars["theta"]) {
			if(spherical) return new sym::sym_num(0);
			return vars["rt"].clone();
		}
		if(x==vars["phi"]) {
			return new sym::sym_num(0);
		}
		return NULL;
	}
	if(!spherical) {
		if(f==vars["rz"]) {
			if(x==vars["zeta"]) {
				maxder=maxder>1?maxder:1; // max. derivative for "r" is maxder+1
				return vars["rzz"].clone();
			}
			if(x==vars["theta"]) {
				maxder=maxder>1?maxder:1;
				return vars["rzt"].clone();
			}
			if(x==vars["phi"]) {
				return new sym::sym_num(0);
			}
			return NULL;
		}
		if(f==vars["rt"]) {
			if(x==vars["zeta"]) {
				maxder=maxder>1?maxder:1;
				return vars["rzt"].clone();
			}
			if(x==vars["theta"]) {
				maxder=maxder>1?maxder:1;
				return vars["rtt"].clone();
			}
			if(x==vars["phi"]) {
				return new sym::sym_num(0);
			}
			return NULL;
		}
	}

	if(axisymmetric) 
		if(x==vars["phi"]) return new sym::sym_num(0);
	
	int der=1, derz=0, dert=0;
	if(x==vars["zeta"]) derz++;
	if(x==vars["theta"]) dert++;
	
	const sym::sym_expr *s=&f;
	
	while(typeid(*s)==typeid(sym::sym_deriv)) {
		if(((sym::sym_deriv *)s)->var == vars["zeta"]) derz++;
		if(((sym::sym_deriv *)s)->var == vars["theta"]) dert++;
		s=((sym::sym_deriv *)s)->oper;
		der++;
	}
	if(*s==vars["rzz"] || *s==vars["rzt"] || *s==vars["rtt"]) {
		der += 2;
		if(*s==vars["rzz"]) {
			maxder=maxder>der-1?maxder:der-1;
			return NULL;
		}
		if(*s==vars["rzt"]) {
			derz += 1; dert += 1;
		}
		if(*s==vars["rtt"]) {
			if(derz < 2) {
				maxder=maxder>der-1?maxder:der-1;
				return NULL;
			}
			dert += 2;
		}
		if(derz+dert != der) return new sym::sym_num(0);
		sym s2;
		if(derz >= 2) {
			s2 = var("rzz");
			derz -= 2;
		}
		else { // dert >= 2
			s2 = var("rtt");
			dert -= 2;
		}
		for(int i=0; i<derz; i++) s2 = Dz(s2);
		for(int i=0; i<dert; i++) s2 = Dt(s2);
		maxder=maxder>der-1?maxder:der-1;
		return s2.expr->clone();
	}

	maxder=maxder>der?maxder:der;
	
	return NULL;

}

matrix symbolic::get_value(const sym::sym_expr &s) {

	if(typeid(s)==typeid(sym::symbol)) {
		return val[((sym::symbol *)&s)->name];
	}
	
	if(typeid(s)==typeid(sym::sym_deriv)) {
		int ndr=0,ndt=0;
		const sym::sym_expr *s1=&s;
		while(typeid(*s1)==typeid(sym::sym_deriv)) {
			if(((sym::sym_deriv *)s1)->var==vars["zeta"]) ndr++;
			if(((sym::sym_deriv *)s1)->var==vars["theta"]) ndt++;
			if(((sym::sym_deriv *)s1)->var==vars["phi"]) return zeros(map.nr,map.nt);
			s1=((sym::sym_deriv *)s1)->oper;
		}
		matrix_block_diag pre(eye(map.D));
		matrix post(eye(map.nt));
		std::string name=((sym::symbol *)s1)->name;
		matrix dt,dt2;
		switch (par[name]) {
			case 00: dt=map.Dt;dt2=map.Dt2;break;
			case 01: dt=map.Dt_01;dt2=map.Dt2_01;break;
			case 10: dt=map.Dt_10;dt2=map.Dt2_10;break;
			case 11: dt=map.Dt_11;dt2=map.Dt2_11;break;
		}
		for(int i=0;i<ndr;i++) pre=(pre,map.D);
		for(int i=0;i<ndt/2;i++) post=(post,dt2);
		if(ndt%2) post=(post,dt);
		return (pre,val[name],post);
	}
	
	return zeros(map.nr,map.nt);

}

sym symbolic::Dz(const sym &s) {
	return diff(s,zeta);
}
sym symbolic::Dt(const sym &s) {
	return diff(s,theta);
}
sym symbolic::Dphi(const sym &s) {
	return diff(s,phi);
}

sym symbolic::det(const sym_tens &g) {
	sym res=sym(0);
	for(int i=0;i<3;i++) {
		res+=g(0,i)*g(1,(i+1)%3)*g(2,(i+2)%3);
		res-=g(0,i)*g(1,(i+2)%3)*g(2,(i+1)%3);
	}
	return res;
}

sym_tens symbolic::inv(const sym_tens &g) {
	
	sym_tens h;
	sym d;
	d=det(g);
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++) {
			h(i,j)=( g((j+1)%3,(i+1)%3)*g((j+2)%3,(i+2)%3) - g((j+2)%3,(i+1)%3)*g((j+1)%3,(i+2)%3) )/d;
		}
	
	
	return h;
}


sym_vec symbolic::contravariant(const sym_vec &s){

	if(s.type==CONTRAVARIANT) return s;
	
	sym_vec snew(CONTRAVARIANT);
	snew.set_context(this);
	
	for(int i=0;i<3;i++) snew(i)=g(i,0)*s(0)+g(i,1)*s(1)+g(i,2)*s(2);

	return snew;

}

sym_vec symbolic::covariant(const sym_vec &s) {

	if(s.type==COVARIANT) return s;
	
	sym_vec snew(COVARIANT);
	snew.set_context(this);
	
	for(int i=0;i<3;i++) snew(i)=g_(i,0)*s(0)+g_(i,1)*s(1)+g_(i,2)*s(2);

	return snew;

}

sym_tens symbolic::contravariant_contravariant(const sym_tens &s) {

	sym_tens snew(s),s2(s);
	snew.set_type(CONTRAVARIANT,CONTRAVARIANT);

	if(s.type[0]!=CONTRAVARIANT) {
		for(int i=0;i<3;i++) 
			for(int j=0;j<3;j++) 
				snew(i,j)=g(i,0)*s(0,j)+g(i,1)*s(1,j)+g(i,2)*s(2,j);	
	}
	if(s.type[1]!=CONTRAVARIANT) {
		s2=snew;
		for(int i=0;i<3;i++) 
			for(int j=0;j<3;j++) 
				snew(i,j)=s2(i,0)*g(0,j)+s2(i,1)*g(1,j)+s2(i,2)*g(2,j);	
	}
	
	return snew;

}

sym_tens symbolic::contravariant_covariant(const sym_tens &s) {

	sym_tens snew(s),s2(s);
	snew.set_type(CONTRAVARIANT,COVARIANT);
	
	if(s.type[0]!=CONTRAVARIANT) {
		for(int i=0;i<3;i++) 
			for(int j=0;j<3;j++) 
				snew(i,j)=g(i,0)*s(0,j)+g(i,1)*s(1,j)+g(i,2)*s(2,j);	
	}
	if(s.type[1]!=COVARIANT) {
		s2=snew;
		for(int i=0;i<3;i++) 
			for(int j=0;j<3;j++) 
				snew(i,j)=s2(i,0)*g_(0,j)+s2(i,1)*g_(1,j)+s2(i,2)*g_(2,j);	
	}
	
	return snew;

}

sym_tens symbolic::covariant_contravariant(const sym_tens &s) {

	sym_tens snew(s),s2(s);
	snew.set_type(COVARIANT,CONTRAVARIANT);

	if(s.type[0]!=COVARIANT) {
		for(int i=0;i<3;i++) 
			for(int j=0;j<3;j++) 
				snew(i,j)=g_(i,0)*s(0,j)+g_(i,1)*s(1,j)+g_(i,2)*s(2,j);	
	}
	if(s.type[1]!=CONTRAVARIANT) {
		s2=snew;
		for(int i=0;i<3;i++) 
			for(int j=0;j<3;j++) 
				snew(i,j)=s2(i,0)*g(0,j)+s2(i,1)*g(1,j)+s2(i,2)*g(2,j);	
	}
	
	return snew;

}

sym_tens symbolic::covariant_covariant(const sym_tens &s) {

	sym_tens snew(s),s2(s);
	snew.set_type(COVARIANT,COVARIANT);

	if(s.type[0]!=COVARIANT) {
		for(int i=0;i<3;i++) 
			for(int j=0;j<3;j++) 
				snew(i,j)=g_(i,0)*s(0,j)+g_(i,1)*s(1,j)+g_(i,2)*s(2,j);	
	}
	if(s.type[1]!=COVARIANT) {
		s2=snew;
		for(int i=0;i<3;i++) 
			for(int j=0;j<3;j++) 
				snew(i,j)=s2(i,0)*g_(0,j)+s2(i,1)*g_(1,j)+s2(i,2)*g_(2,j);	
	}
	
	return snew;

}

double symbolic::perm(int i,int j,int k) {

	if(j==(i+1)%3&&k==(j+1)%3) return 1;
	if(j==(i+2)%3&&k==(j+2)%3) return -1;
	return 0;

}

sym_vec symbolic::gradient(const sym &s) {

	sym_vec grad(COVARIANT);
	grad.set_context(this);
	
	grad(0)=diff(s,x[0]);
	grad(1)=diff(s,x[1]);
	grad(2)=diff(s,x[2]);

	return grad;
}

sym symbolic::christoffel(int i,int j,int k) {
	
	sym G;
	G.context=this;
	
	for(int l=0;l<3;l++) 
		G=G+g(i,l)*( d(g_(l,j),k) + d(g_(l,k),j) - d(g_(j,k),l) );
			
	return G*0.5;	
}

sym symbolic::covderiv(const sym_vec &v,int i,int k) {

	sym dv(d(v(i),k));

	if(v.type==CONTRAVARIANT) {
		for(int l=0;l<3;l++) dv=dv+G(i,l,k)*v(l);
	} else {
		for(int l=0;l<3;l++) dv=dv-G(l,i,k)*v(l);
	}
	
	return dv;

}

sym_tens symbolic::gradient(const sym_vec &v) {

	sym_tens grad(v.type,COVARIANT);
	grad.set_context(this);
	
	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++)
			grad(i,j)=covderiv(v,i,j);

	return grad;
}

sym symbolic::divergence(const sym_vec &v) {

	sym s;
	sym_vec V;
	s.context=this;
	
	V=contravariant(v);
	
	for(int i=0;i<3;i++) 
		s=s+d(V(i),i)+d(sqrt_g,i)/sqrt_g*V(i);
		//s=s+covderiv(V,i,i);
		
		
	return s;

}

sym_vec symbolic::divergence(const sym_tens &t) {

	sym_vec v(CONTRAVARIANT);
	sym_tens T;
	v.set_context(this);
	
	T=contravariant_contravariant(t);
	
	for(int i=0;i<3;i++) 
		for(int j=0;j<3;j++) {
			v(i)=v(i)+d(T(i,j),j);
			for(int k=0;k<3;k++)
				v(i)=v(i)+G(i,k,j)*T(k,j)+G(j,k,j)*T(i,k);
		}
	return v.set_variance(t.type[0]);

}

/// \brief Returns the symbolic expression of the laplacian of \p s
sym symbolic::laplacian(const sym &s){

	return divergence(gradient(s));

}

sym_vec symbolic::curl(const sym_vec &v) {

	sym_vec res(CONTRAVARIANT),V;
	res.set_context(this);
	
	V=covariant(v);
	
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			for(int k=0;k<3;k++)
				res(i)=res(i)+perm(i,j,k)*covderiv(V,k,j)/sqrt_g;
	
	return res.set_variance(v);
	

}


sym_vec symbolic::laplacian(const sym_vec &v) {

	return (div(grad(v))).set_variance(v);

}

sym_tens symbolic::stress(const sym_vec &v) {

	sym_tens t;
	
	t.set_context(this);
	
	t=grad(v)+grad(v).T()-rational(2,3)*g*div(v);
				
	return t;
}


void symbolic::add(const sym &expr,solver *op,int n,std::string type,
				std::string eq_name,std::string var_name,const matrix &d) {

	sym y=var(var_name);	
	int der=maxder;
	if(var_name=="r") der++;
		
	matrix dt,dt2;
	switch (par[var_name]) {
		case 00: dt=map.Dt;dt2=map.Dt2;break;
		case 01: dt=map.Dt_01;dt2=map.Dt2_01;break;
		case 10: dt=map.Dt_10;dt2=map.Dt2_10;break;
		case 11: dt=map.Dt_11;dt2=map.Dt2_11;break;
	}

	for(int i=0;i<=der;i++) {
		for(int j=0;j<=der-i;j++) {
			sym dy=y;
			for(int k=0;k<i;k++) dy=Dz(dy);
			for(int k=0;k<j;k++) dy=Dt(dy);
			if(dy==0) continue;
			sym jac=jacobian(expr,dy);
			if(jac!=0) {
				matrix_block_diag pre(eye(map.D));
				matrix post(eye(map.nt));	
				for(int k=0;k<i;k++) pre=(map.D,pre);
				for(int k=0;k<j/2;k++) post=(post,dt2);	
				if(j%2) post=(post,dt);
				matrix q=d*jac.eval()*ones(map.nr,map.nt);
				if(type=="ex") {
					if(i==0&&j==0) op->add_d(n,eq_name.c_str(),var_name.c_str(),q);
					else if(j==0) op->add_l(n,eq_name.c_str(),var_name.c_str(),q,pre.block(0));
					else if(i==0) op->add_r(n,eq_name.c_str(),var_name.c_str(),q,post);
					else op->add_lr(n,eq_name.c_str(),var_name.c_str(),q,pre.block(0),post);
				} else {
					if(i==0&&j==0) op->add_d(eq_name.c_str(),var_name.c_str(),q);
					else if(j==0) op->add_l(eq_name.c_str(),var_name.c_str(),q,pre);
					else if(i==0) op->add_r(eq_name.c_str(),var_name.c_str(),q,post);
					else op->add_lr(eq_name.c_str(),var_name.c_str(),q,pre,post);
				}
			}
		}
	}	

}


void symbolic::add_bc(const sym &expr,solver *op,int n,std::string type,
				std::string eq_name,std::string var_name,const matrix &d) {
	
	
	sym y=var(var_name);	
	int der=maxder;
	if(var_name=="r") der++;
		
	matrix dt,dt2;
	switch (par[var_name]) {
		case 00: dt=map.Dt;dt2=map.Dt2;break;
		case 01: dt=map.Dt_01;dt2=map.Dt2_01;break;
		case 10: dt=map.Dt_10;dt2=map.Dt2_10;break;
		case 11: dt=map.Dt_11;dt2=map.Dt2_11;break;
	}
			
	int j0=0;
	for(int i=0;i<n;i++) {
		j0+=map.npts[i];
	}
	if(type=="top1") j0+=map.npts[n]-1;
	if(type=="top2") j0+=map.npts[n];
	if(type=="bot1") j0--;
	
	void (solver::*add_d) (int,const char *,const char *,const matrix &);
	void (solver::*add_l) (int,const char *,const char *,const matrix &,const matrix &);
	void (solver::*add_r) (int,const char *,const char *,const matrix &,const matrix &);
	void (solver::*add_lr) (int,const char *,const char *,const matrix &,const matrix &,const matrix &);
	
	if(type=="top1") {
		add_d=&solver::bc_top1_add_d;
		add_l=&solver::bc_top1_add_l;
		add_r=&solver::bc_top1_add_r;
		add_lr=&solver::bc_top1_add_lr;
	}
	if(type=="top2") {
		add_d=&solver::bc_top2_add_d;
		add_l=&solver::bc_top2_add_l;
		add_r=&solver::bc_top2_add_r;
		add_lr=&solver::bc_top2_add_lr;
	}
	if(type=="bot1") {
		add_d=&solver::bc_bot1_add_d;
		add_l=&solver::bc_bot1_add_l;
		add_r=&solver::bc_bot1_add_r;
		add_lr=&solver::bc_bot1_add_lr;
	}
	if(type=="bot2") {
		add_d=&solver::bc_bot2_add_d;
		add_l=&solver::bc_bot2_add_l;
		add_r=&solver::bc_bot2_add_r;
		add_lr=&solver::bc_bot2_add_lr;
	}
	
	matrix D;
	if(type=="top1") D=map.D.block(n);
	if(type=="top2") D=map.D.block(n+1);
	if(type=="bot1") D=map.D.block(n-1);
	if(type=="bot2") D=map.D.block(n);
	
	
	for(int i=0;i<=der;i++) {
		for(int j=0;j<=der-i;j++) {
			sym dy=y;
			for(int k=0;k<i;k++) dy=Dz(dy);
			for(int k=0;k<j;k++) dy=Dt(dy);
			if(dy==0) continue;
			sym jac=jacobian(expr,dy);
			if(jac!=0) {
				matrix pre(eye(D.nrows()));
				matrix post(eye(map.nt));	
				for(int k=0;k<i;k++) pre=(D,pre);
				if(type=="top1") pre=pre.row(-1);
				if(type=="top2") pre=pre.row(0);
				if(type=="bot1") pre=pre.row(-1);
				if(type=="bot2") pre=pre.row(0);
				for(int k=0;k<j/2;k++) post=(post,dt2);	
				if(j%2) post=(post,dt);
				matrix q=d*(jac.eval()*ones(map.nr,map.nt)).row(j0);
				if(i==0&&j==0) (op->*add_d)(n,eq_name.c_str(),var_name.c_str(),q);
				else if(j==0) (op->*add_l)(n,eq_name.c_str(),var_name.c_str(),q,pre);
				else if(i==0) (op->*add_r)(n,eq_name.c_str(),var_name.c_str(),q,post);
				else (op->*add_lr)(n,eq_name.c_str(),var_name.c_str(),q,pre,post);
			}
		}
	}	
				


}


