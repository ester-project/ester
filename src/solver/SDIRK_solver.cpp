// To keep compatibility with configure
#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "solver.h"

extern "C" {
#include <string.h>
#include <stdlib.h>
}
#include <cmath>

SDIRK_solver::SDIRK_solver() {

	initd=false;
	//abs_tol=1e-12;
	//rel_tol=1e-10;

}

void SDIRK_solver::init(int nvar,const char *type) {

	if(initd) destroy();
	nv=nvar;
	nr=new int[nv];
	nc=new int[nv];
	var=new char*[nv];
	reg=new bool[nv];
	y=new matrix[nv];
	dy=new matrix[nv];
	for(int i=0;i<nv;i++) {
		var[i]=new char[32];var[i][0]='\0';
		reg[i]=false;
	}
	step=-1;
	stage=-1;
	
	if(!strcmp(type,"be")) 
		init_be();
	else if(!strcmp(type,"cn"))
		init_cn();
	else if(!strcmp(type,"sdirk2"))
		init_sdirk2();
	else if(!strcmp(type,"sdirk3"))
		init_sdirk3();
	else if(!strcmp(type,"sdirk4"))
		init_sdirk4();
	else if(!strcmp(type,"esdirk3"))
		init_esdirk3();
	else if(!strcmp(type,"esdirk4"))
		init_esdirk4();
	else {
		fprintf(stderr,"SDIRK_solver: Unknown RK method %s\n",type);
		exit(1);
	}
	
	check_method();
	if(nstages>1)
		k=new matrix[nstages-1];
	for(int n=0;n<nstages-1;n++) k[n]=zeros(1,1);
	
	initd=true;
	
}

void SDIRK_solver::destroy() {

	if(!initd) return;
	delete [] nr;
	delete [] nc;
	delete [] reg;
	for(int i=0;i<nv;i++) {
		delete [] var[i];
	}
	delete [] var;
	delete [] y;
	if(nstages>1) delete [] k;
	delete [] dy;
	initd=0;

}

void SDIRK_solver::check_init() {
	if(!initd) {
		fprintf(stderr,"SDIRK_solver not initialized\n");
		exit(1);
	}
}

void SDIRK_solver::check_method() {

	nstages=a.nrows();
	alpha=a(-1,-1);
	
	if(a.ncols()!=nstages) {
		fprintf(stderr,"SDIRK_solver: Invalid method (a is not square)\n");
		exit(1);
	}
	if(b.nrows()!=nstages) {
		fprintf(stderr,"SDIRK_solver: Invalid method (Incorrect size of vector b)\n");
		exit(1);
	}
	if(c.nrows()!=nstages) {
		fprintf(stderr,"SDIRK_solver: Invalid method (Incorrect size of vector c)\n");
		exit(1);
	}
	for(int i=0;i<nstages-1;i++) {
		for(int j=i+1;j<nstages;j++) {
			if(a(i,j)!=0) {
				fprintf(stderr,"SDIRK_solver: Method is not diagonal implicit (a_ij!=0 for j>i)\n");
				exit(1);
			}
		}
	}
	if(c(-1)!=1) {
		fprintf(stderr,"SDIRK_solver: Invalid method (c_n!=1)\n");
		exit(1);
	}
	if(exist(a.row(-1).transpose()!=b)) {
		fprintf(stderr,"SDIRK_solver: Invalid method (a_ni!=b_i)\n");
		exit(1);
	}
	if(a(0,0)==0) {
		first_explicit=true;
		if(c(0)!=0) {
			fprintf(stderr,"SDIRK_solver: Invalid method (c_1 should be 0 for an ESDIRK method)\n");
			exit(1);
		}
	} else first_explicit=false;
	
	for(int i=0;i<nstages;i++) {
		if(i==0&&first_explicit) continue;
		if(a(i,i)!=alpha) {
			fprintf(stderr,"SDIRK_solver: Invalid method (Diagonal terms in a are not equal)\n");
			exit(1);
		}
	}
	
	if(alpha==0) {
		fprintf(stderr,"SDIRK_solver: Method is explicit\n");
		exit(1);
	}
	
	double tol=1e-14;
	if(std::abs(sum(b)-1)>tol) {
		fprintf(stderr,"SDIRK_solver: Invalid method (sum(b_i)!=1)\n");
		exit(1);
	}
	order=1;
	matrix T=c*eye(nstages);
	while(true) {
		int counter=0;
		while(!(counter>>order)) {
			matrix m=eye(nstages);
			double fact=1./(order+1.);
			for(int n=0;n<order;n++) {
				if((counter>>n)&1) {
					m=(a,m);
					fact*=1./(n+1.);
				} else m=(T,m);
			}
			if(fabs(sum((b.transpose(),m))-fact)>tol) return;
			counter++;
		}
		order++;if(order>3) break;
	}
	
}

double SDIRK_solver::get_step() {return step;}
int SDIRK_solver::number_of_stages() {return nstages;}
int SDIRK_solver::number_of_implicit_stages() {return first_explicit?nstages-1:nstages;}
int SDIRK_solver::get_order() {return order;}
bool SDIRK_solver::needs_initial_derivative() {return first_explicit;}

void SDIRK_solver::regvar(const char *var_name,const matrix &initial_value) {

	check_init();
	int j;
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

void SDIRK_solver::regvars(const matrix_map &vars) {

	matrix_map::const_iterator it;
	for(it=vars.begin();it!=vars.end();it++) 
		regvar((it->first).c_str(),it->second);

}

void SDIRK_solver::set_initial_derivative(const char *var_name,const matrix &initial_value) {

	if(!first_explicit) return;
	dy[get_id(var_name)]=initial_value;

}

void SDIRK_solver::set_initial_derivatives(const matrix_map &vars) {

	if(!first_explicit) return;
	matrix_map::const_iterator it;
	for(it=vars.begin();it!=vars.end();it++) 
		set_initial_derivative((it->first).c_str(),it->second);

}

void SDIRK_solver::set_var(const char *var_name,const matrix &value) {
	y[get_id(var_name)]=value;
}

void SDIRK_solver::set_vars(const matrix_map &values) {
	matrix_map::const_iterator it;
	for(it=values.begin();it!=values.end();it++) 
		set_var((it->first).c_str(),it->second);
}

void SDIRK_solver::set_step(double h) {
	step=h;
}

int SDIRK_solver::get_id(const char *varn) {

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

double SDIRK_solver::get_t() {
	return t_eval;
}

double SDIRK_solver::get_delta() {
	return delta;
}

matrix SDIRK_solver::get_var(const char *var_name) {
	return y[get_id(var_name)];
}

matrix_map SDIRK_solver::get_vars() {

	matrix_map map;
	for(int i=0;i<nv;i++) 
		if(reg[i]) map[std::string(var[i])]=y[i];

	return map;
}


void SDIRK_solver::wrap(const matrix *u,matrix *v) {

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

void SDIRK_solver::unwrap(const matrix *v,matrix *u) {

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

int SDIRK_solver::solve(double t0,double t1) {

	check_init();
	
	static double h;
	static bool last_step=false;
	int ret_code=RK_INTERMEDIATE;
	
	matrix xi;
	wrap(y,&xi);
	
	if(stage==-1) { 
		t=t0;
		if(step<=0) step=t1-t0;
		if(first_explicit) 
			wrap(dy,&k[0]);
		last_step=false;
	} else if(stage==nstages-1) {
		t+=h;
		if(first_explicit) 
			k[0]=(xi-x0)/delta;
	} else 
		k[stage]=(xi-x0)/delta;
	
	stage=(stage+1)%nstages;
	
	if(stage==0) {
		x=xi;
		if(t>=t1||last_step) {
			stage=-1;
			if (first_explicit) unwrap(&k[0],dy);
			last_step=false;
			return RK_END;
		}
		if(step>t1-t) {
			h=t1-t;
			last_step=true;
		} else h=step;
		if(first_explicit) stage=1;
	}
	
	t_eval=t+c(stage)*h;
	delta=alpha*h;
	x0=x;
	for(int n=0;n<stage;n++) 
		x0+=h*a(stage,n)*k[n];
	unwrap(&x0,y);
	
	if(stage==nstages-1) ret_code=RK_STEP;
	
	return ret_code;
	
}

void SDIRK_solver::init_be() {
	// Backward Euler
	a.dim(1,1);b.dim(1,1);c.dim(1,1);
	a(0,0)=1;b(0)=1;c(0)=1;
	
}

void SDIRK_solver::init_cn() {
	// Crank-Nicolson
	a.dim(2,2);b.dim(2,1);c.dim(2,1);
	
	a.values( 	0., 	0.,
				0.5,	0.5);
	b.values(	0.5,	0.5);
	c.values(	0.,		1.);
	
}

void SDIRK_solver::init_sdirk2() {
	// sdirk2 Alexander, "Diagonally implicit Ruge-Kutta methods for stiff o.d.e.'s",
	//             SIAM J. Numer. Anal. vol 14, n. 6 (1977)
	a.dim(2,2);b.dim(2,1);c.dim(2,1);
	double q=1.-sqrt(2.)/2.;
	
	a.values( 	q,	 	0.,
				1-q,	q);
	b.values(	1-q,	q);
	c.values(	q,		1.);
	
}

void SDIRK_solver::init_sdirk3() {
	// sdirk3 Alexander, "Diagonally implicit Ruge-Kutta methods for stiff o.d.e.'s",
	//             SIAM J. Numer. Anal. vol 14, n. 6 (1977)
	a.dim(3,3);b.dim(3,1);c.dim(3,1);
	double q=0.4358665215084590038863155;
	double c2,b1,b2;
	c2=(1+q)/2.;
	b1=-(6*q*q-16*q+1)/4.;
	b2=(6*q*q-20*q+5)/4.;
	
	a.values( 	q,		0.,		0.,
				c2-q,	q,		0.,
				b1,		b2,		q);
	b.values(	b1,		b2,		q);
	c.values(	q,		c2,		1.);
	
}

void SDIRK_solver::init_sdirk4() {
	// S54b Skvortsov, "Diagonally implicit Ruge-Kutta methods for stiff problems",
	//       Computational Mathematics and Mathematical Physics vol 46, n. 12 (2006)
	a.dim(5,5);b.dim(5,1);c.dim(5,1);
	
	a.values( 	1./4,	0.,		0.,		0.,		0.,
				-1./4,	1./4,	0.,		0.,		0.,
				1./8,	1./8,	1./4,	0.,		0.,
				-3./2,	3./4,	3./2,	1./4,	0.,
				0.,		1./6,	2./3,	-1./12,	1./4);
	b.values(	0.,		1./6,	2./3,	-1./12,	1./4);
	c.values(	1./4,	0.,		1./2,	1.,		1.);
	
}

void SDIRK_solver::init_esdirk3() {
	// SEGAWA, HIDEHIRO. A Family of Higher-Order Implicit Time Integration Methods
	// for Unsteady Compressible Flows.
	a.dim(4,4);b.dim(4,1);c.dim(4,1);
	
	a.values(
		0.,								0.,								0.,								0.,
		1767732205903./4055673282236., 	1767732205903./4055673282236., 	0., 							0.,
		2746238789719./10658868560708., -640167445237./6845629431997., 	1767732205903./4055673282236., 	0.,
		1471266399579./7840856788654., 	-4482444167858./7529755066697.,	11266239266428./11593286722821.,1767732205903./4055673282236.);
	b=a.row(-1).transpose();
	c.values(0., 1767732205903./2027836641118., 3./5., 1.);
	
}

void SDIRK_solver::init_esdirk4() {
	// SEGAWA, HIDEHIRO. A Family of Higher-Order Implicit Time Integration Methods
	// for Unsteady Compressible Flows.
	a.dim(6,6);b.dim(6,1);c.dim(6,1);
	
	a.values(
		0., 						0., 					0., 					0., 				0., 			0.,
		1./4., 						1./4., 					0., 					0., 				0., 			0.,
		8611./62500., 				-1743./31250., 			1./4., 					0., 				0., 			0.,
		5012029./34652500., 		-654441./2922500., 		174375./388108., 		1./4., 				0., 			0.,
		15267082809./155376265600.,	-71443401./120774400.,	730878875./902184768.,	2285395./8070912.,	1./4.,			0.,
		82889./524892., 			0.,	 					15625./83664., 			69875./102672., 	-2260./8211., 	1./4.);
	b=a.row(-1).transpose();	
	c.values(0., 1./2., 83./250., 31./50., 17./20., 1.);
	
}



