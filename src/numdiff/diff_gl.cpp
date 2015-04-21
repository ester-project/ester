#include "ester-config.h"
#include "numdiff.h"
#include "constants.h"
#include "utils.h"

extern "C" {
#include <stdlib.h>
#include <stdarg.h>
}

diff_gl::diff_gl(int n):ndomains(ndom), N(Ntot) {
    DEBUG_FUNCNAME;
	ndom=n;
	npts=new int[ndom];
	xif=new double[ndom+1];
}

diff_gl::~diff_gl() {

	delete [] npts;
	delete [] xif;
}

diff_gl::diff_gl(const diff_gl &gl) : x(gl.x), I(gl.I), P(gl.P), P1(gl.P1),
    D(gl.D), ndomains(ndom), N(Ntot) {
    DEBUG_FUNCNAME;
	int i;

	ndom=gl.ndom;
	npts=new int[ndom];;
	xif=new double[ndom+1];
	for(i=0;i<ndom;i++)
		*(npts+i)=*(gl.npts+i);
	for(i=0;i<ndom+1;i++)
		*(xif+i)=*(gl.xif+i);
	Ntot=gl.Ntot;

}

diff_gl &diff_gl::operator=(const diff_gl &gl) {

	int i;

	delete [] npts;
	delete [] xif;
	ndom=gl.ndom;
	npts=new int[ndom];
	xif=new double[ndom+1];
	for(i=0;i<ndom;i++)
		*(npts+i)=*(gl.npts+i);
	for(i=0;i<ndom+1;i++)
		*(xif+i)=*(gl.xif+i);
	Ntot=gl.Ntot;
	x=gl.x;P=gl.P;P1=gl.P1;D=gl.D;I=gl.I;

	return *this;
}

void diff_gl::set_ndomains(int n) {
    DEBUG_FUNCNAME;
	ndom=n;
	delete [] npts;
	delete [] xif;
	npts=new int[ndom];
	xif=new double[ndom+1];

}

void diff_gl::set_xif(double x,...) {

	va_list ap;
	int i;
	
	va_start(ap,x);
	*xif=x;
	for(i=1;i<ndom+1;i++) {
		*(xif+i)=va_arg(ap,double);
	}
	va_end(ap);
}

void diff_gl::set_npts(int n,...) {
    DEBUG_FUNCNAME;
	va_list ap;
	int i;
	
	va_start(ap,n);
	*npts=n;
	for(i=1;i<ndom;i++) {
		*(npts+i)=va_arg(ap,int);
	}
	va_end(ap);
}

void diff_gl::init() {
    DEBUG_FUNCNAME;
	diff_gl gl1(1);
	int i,j1,j2;
	double x1,x2;
	matrix T;


	Ntot=0;
	for(i=0;i<ndom;i++)
		Ntot+=*(npts+i);
	x=zeros(Ntot,1);
	P.set_nblocks(ndom);
	P1.set_nblocks(ndom);
	D.set_nblocks(ndom);
	I=zeros(1,Ntot);
	j1=0;
	for(i=0;i<ndom;i++) {
		j2=j1+*(npts+i)-1;
		*(gl1.npts)=*(npts+i);
		gl1.init_1();
		x1=xif[i];
		x2=xif[i+1];
		x.setblock(j1,j2,0,0,gl1.x*(x2-x1)+x1);
		P.block(i)=gl1.P.block(0);
		P1.block(i)=gl1.P1.block(0);
		D.block(i)=gl1.D.block(0)/(x2-x1);
		I.setblock(0,0,j1,j2,gl1.I*(x2-x1));
		j1=j1+*(npts+i);
	}
}

void diff_gl::init_1() {

	int n=*npts;
	matrix i(n,1),j(1,n);
	
	xif[0]=0;xif[1]=1;
	i=vector_t(0,n-1,n);
	j=vector(0,n-1,n);
	x=-cos(i*(PI/(n-1)));
	x=(x+1)/2;
	
	
	P.block(0)=((2-((i==0)||(i==n-1)))*pow(-1,i)/(n-1))/(1+((j==0)||(j==n-1)))*
		cos(i*(j*(PI/(n-1))));
	P1.block(0)=cos((j*PI)*(1-i/(n-1)));
	D.block(0)=(i<j)*((i+j)/2!=floor((i+j)/2))*(2*j);
	D.block(0).setrow(0,D.block(0).row(0)/2);
	D.block(0)=D.block(0)*2;
	D.block(0)=(P1.block(0),D.block(0),P.block(0));
	/*
	matrix th;
	int k,kk;
	th=PI*(1-i/(n-1));
	P1.block(0)=cos(j*th);
	P.block(0)=P1.block(0).inv();
	D.block(0)=j*sin(j*th)/sin(th);
	D.block(0).setrow(n-1,j*j*ones(1,n));
	for(k=0,kk=-1;k<n-1;k++,kk=-kk) 
		D.block(0)(0,k)=kk*k*k;
	D.block(0)=2*(D.block(0),P.block(0));
*/
	I=zeros(1,n);
	i=vector(0,n-1-(n-1)%2,(n+1)/2);
	I.setblock_step(0,0,1,0,I.ncols()-1,2,2/(1-i*i));
	I=(I,P.block(0));
	I=I/2;
	
}

matrix diff_gl::eval(const matrix &y,const matrix &x) const {

	matrix T;
	return this->eval(y,x,T);
}

matrix diff_gl::eval(const matrix &y,double x) const {

	matrix T,X(1,1);
	X(0,0)=x;
	return this->eval(y,X,T);
}

matrix diff_gl::eval(const matrix &y,double x,matrix &T) const {

	matrix X(1,1);
	X(0,0)=x;
	return this->eval(y,X,T);
}

matrix_map diff_gl::eval(const matrix_map &y, const matrix &x) const {
    matrix_map res;
    matrix_map::const_iterator it;
    for (it=y.begin(); it!=y.end(); it++)  {
        res[it->first] = this->eval(it->second, x);
    }
    return res;
}


matrix diff_gl::eval(const matrix &y,const matrix &x,matrix &T) const {

	matrix yi(x.nrows(),y.ncols()),Ti(x.nrows(),1);
	int i,j,k;
	
	T.dim(x.nrows(),y.nrows());

	k=0;
	for(i=0;i<ndom;i++) {
		for(j=0;j<npts[i];j++) {
			Ti=2*(x-xif[i])/(xif[i+1]-xif[i])-1;
			if(i==0) Ti=Ti*(x>=xif[i]&&x<=xif[i+1]);
				else Ti=Ti*(x>xif[i]&&x<=xif[i+1]);
			Ti=cos(j*acos(Ti));
			if(i==0) Ti=Ti*(x>=xif[i]&&x<=xif[i+1]);
				else Ti=Ti*(x>xif[i]&&x<=xif[i+1]);
			T.setcol(k,Ti);
			k++;
		}
	}
	T=(T,P);
	yi=(T,y);

	return yi;
	
}







