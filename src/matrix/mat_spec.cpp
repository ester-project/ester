#include<stdio.h>
#include"matrix.h"
#include<cmath>
#include<stdlib.h>
#include<sys/time.h>
extern "C" {
#ifdef USE_MKL
#include <mkl_lapack.h>
#include <mkl_cblas.h>
#else
#include <lapack.h>
#include <cblas.h>
#endif
}

matrix matrix::operator,(const matrix &a) const {

	if(nc!=a.nf) {
		fprintf(stderr,"ERROR: (matrix.,) Dimensions must agree\n");
		exit(1);
	}

	matrix res(nf,a.nc);

	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,nf,a.nc,nc,1,p,nf,a.p,a.nf,0,res.p,res.nf);

	return res;
}

matrix matrix::solve(matrix b) const {

	int ipiv[nf],info;
	char trans='N';
	matrix lu(*this);

	if(nf!=b.nf) {
		fprintf(stderr,"ERROR: (matrix.solve) Dimensions must agree\n");
		exit(1);
	}
	if(nf!=nc) {
		fprintf(stderr,"ERROR: (matrix.solve) Matrix must be square\n");
		exit(1);
	}

	dgetrf_(&lu.nf,&lu.nc,lu.p,&lu.nf,ipiv,&info);
	if(info>0) {
		printf("Matrix is singular (%d)\n",info);
		exit(1);
	}
	dgetrs_(&trans,&lu.nf,&b.nc,lu.p,&lu.nf,ipiv,b.p,&b.nf,&info);

	return b;

}

matrix matrix::inv() const {

	matrix res(*this);
	int ipiv[nf],info,lwork=-1;
	double *work,w;
	
	if(nf!=nc) {
		fprintf(stderr,"ERROR: (matrix.inv) Matrix must be square\n");
		exit(1);
	}
	
	dgetrf_(&res.nf,&res.nc,res.p,&res.nf,ipiv,&info);
	if(info>0) {
		printf("Matrix is singular (%d)\n",info);
		exit(1);
	}
	dgetri_(&res.nf,res.p,&res.nf,ipiv,&w,&lwork,&info);
	lwork=round(w);

	work=new double[lwork];
	
	dgetri_(&res.nf,res.p,&res.nf,ipiv,work,&lwork,&info);
	
	delete [] work;
	
	return res;

}


