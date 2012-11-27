#include<cmath>
#include<stdio.h>
#include<stdlib.h>
#include"matrix.h"

template <class Type>
Matrix<Type> cos(const Matrix<Type> &a) {
    
    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=cos(a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> sin(const Matrix<Type> &a) {
    
    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=sin(a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> tan(const Matrix<Type> &a) {
    
    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=tan(a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> acos(const Matrix<Type> &a) {
    
    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=acos(a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> asin(const Matrix<Type> &a) {
    
    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=asin(a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> atan(const Matrix<Type> &a) {
    
    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=atan(a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> cosh(const Matrix<Type> &a) {
    
    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=cosh(a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> sinh(const Matrix<Type> &a) {
    
    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=sinh(a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> tanh(const Matrix<Type> &a) {
    
    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=tanh(a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> exp(const Matrix<Type> &a) {
    
    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=exp(a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> log(const Matrix<Type> &a) {
    
    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=log(a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> log10(const Matrix<Type> &a) {
    
    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=log10(a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> sqrt(const Matrix<Type> &a) {
    
    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=sqrt(a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> abs(const Matrix<Type> &a) {
    
    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=fabs(a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> floor(const Matrix<Type> &a) {
    
    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=floor(a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> ceil(const Matrix<Type> &a) {
    
    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=ceil(a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> round(const Matrix<Type> &a) {
    
    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=round(a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> atan2(const Matrix<Type> &a,const Matrix<Type> &b) {

	Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    if( (b.nf!=a.nf) || (b.nc!=a.nc) ) {
		fprintf(stderr,"ERROR: (matrix.atan2) Dimensions must agree\n");
		exit(1);
	}
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=atan2(a.p[i],b.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> atan2(double n,const Matrix<Type> &a) {

	Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=atan2(n,a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> atan2(const Matrix<Type> &a,double n) {

	Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=atan2(a.p[i],n);
    	
    return res;
}

template <class Type>
Matrix<Type> pow(const Matrix<Type> &a,const Matrix<Type> &b) {

	Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    if( (b.nf!=a.nf) || (b.nc!=a.nc) ) {
		fprintf(stderr,"ERROR: (matrix.pow) Dimensions must agree\n");
		exit(1);
	}
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=pow(a.p[i],b.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> pow(double n,const Matrix<Type> &a) {

	Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=pow(n,a.p[i]);
    	
    return res;
}

template <class Type>
Matrix<Type> pow(const Matrix<Type> &a,double n) {

	Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=pow(a.p[i],n);
    	
    return res;
}

template <class Type> 
Matrix<Type> pow(const Matrix<Type> &a,int n) {

	Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	// There is no pow(double,int) in the Intel compiler math library
    	#ifdef __INTEL_COMPILER
    	res.p[i]=pow(a.p[i],(mat_type)n);
    	#else
    	res.p[i]=pow(a.p[i],n);
    	#endif
    	
    return res;
}

// Explicit instantiations

template Matrix<double> cos(const Matrix<double> &);
template Matrix<double> sin(const Matrix<double> &);
template Matrix<double> tan(const Matrix<double> &);
template Matrix<double> acos(const Matrix<double> &);
template Matrix<double> asin(const Matrix<double> &);
template Matrix<double> atan(const Matrix<double> &);
template Matrix<double> cosh(const Matrix<double> &);
template Matrix<double> sinh(const Matrix<double> &);
template Matrix<double> tanh(const Matrix<double> &);
template Matrix<double> exp(const Matrix<double> &);
template Matrix<double> log(const Matrix<double> &);
template Matrix<double> log10(const Matrix<double> &);
template Matrix<double> sqrt(const Matrix<double> &);
template Matrix<double> abs(const Matrix<double> &);
template Matrix<double> floor(const Matrix<double> &);
template Matrix<double> ceil(const Matrix<double> &);
template Matrix<double> round(const Matrix<double> &);
template Matrix<double> atan2(const Matrix<double> &,const Matrix<double> &);
template Matrix<double> atan2(double,const Matrix<double> &);
template Matrix<double> atan2(const Matrix<double> &,double);
template Matrix<double> pow(const Matrix<double> &,const Matrix<double> &);
template Matrix<double> pow(double,const Matrix<double> &);
template Matrix<double> pow(const Matrix<double> &,double);
template Matrix<double> pow(const Matrix<double> &,int);

template Matrix<float> cos(const Matrix<float> &);
template Matrix<float> sin(const Matrix<float> &);
template Matrix<float> tan(const Matrix<float> &);
template Matrix<float> acos(const Matrix<float> &);
template Matrix<float> asin(const Matrix<float> &);
template Matrix<float> atan(const Matrix<float> &);
template Matrix<float> cosh(const Matrix<float> &);
template Matrix<float> sinh(const Matrix<float> &);
template Matrix<float> tanh(const Matrix<float> &);
template Matrix<float> exp(const Matrix<float> &);
template Matrix<float> log(const Matrix<float> &);
template Matrix<float> log10(const Matrix<float> &);
template Matrix<float> sqrt(const Matrix<float> &);
template Matrix<float> abs(const Matrix<float> &);
template Matrix<float> floor(const Matrix<float> &);
template Matrix<float> ceil(const Matrix<float> &);
template Matrix<float> round(const Matrix<float> &);
template Matrix<float> atan2(const Matrix<float> &,const Matrix<float> &);
template Matrix<float> atan2(double,const Matrix<float> &);
template Matrix<float> atan2(const Matrix<float> &,double);
template Matrix<float> pow(const Matrix<float> &,const Matrix<float> &);
template Matrix<float> pow(double,const Matrix<float> &);
template Matrix<float> pow(const Matrix<float> &,double);
template Matrix<float> pow(const Matrix<float> &,int);

