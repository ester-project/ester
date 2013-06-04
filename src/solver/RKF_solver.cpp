#include"solver.h"
#include<string.h>
#include<stdlib.h>
#include<cmath>

RKF_solver::RKF_solver() {

	initd=false;
	abs_tol=1e-12;
	rel_tol=1e-10;

}

void RKF_solver::init(int nvar) {

	if(initd) destroy();
	nv=nvar;
	y=new matrix[nv];
	dy=new matrix[nv];
	nr=new int[nv];
	nc=new int[nv];
	var=new char*[nv];
	reg=new bool[nv];
	for(int i=0;i<nv;i++) {
		var[i]=new char[32];var[i][0]='\0';
		reg[i]=false;
	}
	step=-1;
	stage=-1;
	a=zeros(6,6);c=zeros(6,1);b4=zeros(6,1);b5=zeros(6,1);
	a(1,0)=1./4.;
	a(2,0)=3./32.;a(2,1)=9./32.;
	a(3,0)=1932./2197.;a(3,1)=-7200./2197.;a(3,2)=7296./2197.;
	a(4,0)=439./216.;a(4,1)=-8.;a(4,2)=3680./513.;a(4,3)=-845./4104.;
	a(5,0)=-8./27.;a(5,1)=2.;a(5,2)=-3544./2565.;a(5,3)=1859./4104.;a(5,4)=-11./40.;
	c(0)=0.;c(1)=1./4.;c(2)=3./8.;c(3)=12./13.;c(4)=1.;c(5)=1./2.;
	b4(0)=25./216.;b4(1)=0.;b4(2)=1408./2565.;b4(3)=2197./4104.;b4(4)=-1./5.;
	b5(0)=16./135.;b5(1)=0.;b5(2)=6656./12825.;b5(3)=28561./56430.;b5(4)=-9./50.;b5(5)=2./55.;
	
	initd=true;
	
}

void RKF_solver::destroy() {

	if(!initd) return;
	delete [] y;
	delete [] dy;
	delete [] nr;
	delete [] nc;
	delete [] reg;
	for(int i=0;i<nv;i++) {
		delete [] var[i];
	}
	delete [] var;
	initd=0;

}

void RKF_solver::check_init() {
	if(!initd) {
		fprintf(stderr,"RKF_solver not initialized\n");
		exit(1);
	}
}

void RKF_solver::regvar(const char *var_name,const matrix &initial_value) {

	check_init();
	int i,j;
	j=0;
	while (strlen(var[j])) {
		if(!strcmp(var[j],var_name)) {
			fprintf(stderr,"ERROR: Can't register variable (already registered)\n");
			exit(1);
		}
		j++;
		if(j==nv) {
			fprintf(stderr,"ERROR: Can't register variable (increase nvar)\n");
			exit(1);
		}
	}	

	strncpy(var[j],var_name,31);
	y[j]=initial_value;
	nr[j]=y[j].nrows();
	nc[j]=y[j].ncols();
	reg[j]=true;

}

void RKF_solver::regvars(const matrix_map &vars) {

	matrix_map::const_iterator it;
	for(it=vars.begin();it!=vars.end();it++) 
		regvar((it->first).c_str(),it->second);

}

void RKF_solver::set_step(double h) {
	step=h;
}

int RKF_solver::get_id(const char *varn) {

	int i=0;

	while(strcmp(varn,var[i])||!reg[i]) {
		i++;
		if(i==nv) {
			fprintf(stderr,"ERROR: Unknown variable %s\n",varn);
			exit(1);
		}
	}
	return i;

}

matrix RKF_solver::get_var(const char *var_name) {
	return y[get_id(var_name)];
}

matrix_map RKF_solver::get_vars() {

	matrix_map map;
	for(int i=0;i<nv;i++) 
		if(reg[i]) map[std::string(var[i])]=y[i];

	return map;
}

double RKF_solver::get_t() {
	return t_eval;
}

void RKF_solver::set_deriv(const char *var_name,const matrix &value) {
	dy[get_id(var_name)]=value;
}

void RKF_solver::set_derivs(const matrix_map &values) {
	matrix_map::const_iterator it;
	for(it=values.begin();it!=values.end();it++) 
		set_deriv((it->first).c_str(),it->second);
}

void RKF_solver::wrap(const matrix *u,matrix *v) {

	int N=0;
	for(int i=0;i<nv;i++) {
		if(reg[i]) N+=nr[i]*nc[i];
	}
	*v=zeros(N,1);
	int j=0;
	for(int i=0;i<nv;i++) {
		if(reg[i]) {
			matrix w(u[i]);
			N=nr[i]*nc[i];
			w.redim(N,1);
			v->setblock(j,j+N-1,0,0,w);
			j+=N;
		}
	}

}

void RKF_solver::unwrap(const matrix *v,matrix *u) {

	int j=0,N;
	for(int i=0;i<nv;i++) {
		if(reg[i]) {
			N=nr[i]*nc[i];
			u[i]=v->block(j,j+N-1,0,0);
			u[i].redim(nr[i],nc[i]);
			j+=N;
		}
	}
}

int RKF_solver::solve(double t0,double t1) {

	check_init();
	
	static double h;
	static bool last_step=false;
	int ret_code=RK_INTERMEDIATE;
	
	if(stage==-1) { 
		t=t0;
		wrap(y,&x);
		if(step<=0) step=t1-t0;
	} else {
		wrap(dy,&deriv);
		k[stage]=h*deriv;
	}
	
	if(stage==5) {
		matrix x4(x),x5(x);
		for(int i=0;i<6;i++) {
			x4+=b4(i)*k[i];
			x5+=b5(i)*k[i];
		}
		matrix tol;
		double q;
		tol=max(abs(x5)*rel_tol,abs_tol);
		q=max(abs(x5-x4)/tol);
		if(q<1) {
			x=x5;t+=h;
			ret_code=RK_STEP;
		}
		q=0.84*pow(q,-0.25);
		q=q>4?4:q;
		q=q<0.1?0.1:q;
		step=q*step;
	}
	
	stage=(stage+1)%6;
	
	if(stage==0) {
		if(t>=t1) {
			t_eval=t;
			unwrap(&x,y);
			stage=-1;
			last_step=false;
			return RK_END;
		}
		if(step>t1-t) {
			h=t1-t;
			last_step=true;
		} else h=step;
	}
	
	t_eval=t+c(stage)*h;
	x_eval=x;
	for(int i=0;i<stage;i++) {
		x_eval+=a(stage,i)*k[i];
	}
	
	unwrap(&x_eval,y);
	
	return ret_code;
	
}





