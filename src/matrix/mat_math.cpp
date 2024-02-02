#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "utils.h"


matrix cos(const matrix &a) {
    
    matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;

    for(i=0;i<N;i++) {

    	res.p[i]=cos(a.p[i]);
    }
    	
    return res;
}


matrix sin(const matrix &a) {
    
    matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;

    for(i=0;i<N;i++)
    	res.p[i]=sin(a.p[i]);
    	
    return res;
}


matrix tan(const matrix &a) {
    
    matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;

    for(i=0;i<N;i++)
    	res.p[i]=tan(a.p[i]);
    	
    return res;
}


matrix acos(const matrix &a) {
    
    matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;

    for(i=0;i<N;i++)
    	res.p[i]=acos(a.p[i]);
    	
    return res;
}


matrix asin(const matrix &a) {
    
    matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=asin(a.p[i]);
    	
    return res;
}


matrix atan(const matrix &a) {
    
    matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=atan(a.p[i]);
    	
    return res;
}


matrix cosh(const matrix &a) {
    
    matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=cosh(a.p[i]);
    	
    return res;
}


matrix sinh(const matrix &a) {
    
    matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=sinh(a.p[i]);
    	
    return res;
}


matrix tanh(const matrix &a) {
    
    matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=tanh(a.p[i]);
    	
    return res;
}


matrix exp(const matrix &a) {
    
    matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=exp(a.p[i]);
    	
    return res;
}


matrix log(const matrix &a) {
    
    matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=log(a.p[i]);
    	
    return res;
}


matrix log10(const matrix &a) {
    
    matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=log10(a.p[i]);
    	
    return res;
}


matrix sqrt(const matrix &a) {
    
    matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=sqrt(a.p[i]);
    	
    return res;
}


matrix abs(const matrix &a) {
    
    matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=fabs(a.p[i]);
    	
    return res;
}


matrix floor(const matrix &a) {
    
    matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=floor(a.p[i]);
    	
    return res;
}


matrix ceil(const matrix &a) {
    
    matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=ceil(a.p[i]);
    	
    return res;
}


matrix round(const matrix &a) {
    
    matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=round(a.p[i]);
    	
    return res;
}


matrix atan2(const matrix &a,const matrix &b) {

	matrix res(a.nf,a.nc);
    int i,N;
    
    if( (b.nf!=a.nf) || (b.nc!=a.nc) ) {
        ester_critical("(matrix.atan2) Dimensions must agree");
    }
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=atan2(a.p[i],b.p[i]);
    	
    return res;
}


matrix atan2(double n,const matrix &a) {

	matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=atan2(n,a.p[i]);
    	
    return res;
}


matrix atan2(const matrix &a,double n) {

	matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=atan2(a.p[i],n);
    	
    return res;
}


matrix pow(const matrix &a,const matrix &b) {

	matrix res(a.nf,a.nc);
    int i,N;
    
    if( (b.nf!=a.nf) || (b.nc!=a.nc) ) {
        ester_critical("(matrix.pow) Dimensions must agree");
    }
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=pow(a.p[i],b.p[i]);
    	
    return res;
}


matrix pow(double n,const matrix &a) {

	matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=pow(n,a.p[i]);
    	
    return res;
}


matrix pow(const matrix &a,double n) {

	matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;

    for(i=0;i<N;i++)
    	res.p[i]=pow(a.p[i],n);
    	
    return res;
}

matrix pow(const matrix &a,int n) {

	matrix res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    
    for(i=0;i<N;i++)
    	res.p[i]=pow(a.p[i],n);
    	
    return res;
}


