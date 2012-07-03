#include<stdio.h>
#include"matrix.h"
#include<math.h>
#include<stdlib.h>
#include<sys/time.h>
extern "C" {
#include CBLAS
}
#include LAPACK

matrix matrix::operator,(const matrix &a) const {

	matrix res(nf,a.nc);

	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,nf,a.nc,nc,1,p,nf,a.p,a.nf,0,res.p,res.nf);

	return res;
}

matrix matrix::lu(int *&ipiv) const {

	int info;
	char trans='N';
	matrix lu(*this);

	ipiv=(int *)realloc(ipiv,nf*sizeof(int));

	dgetrf_(&lu.nf,&lu.nc,lu.p,&lu.nf,ipiv,&info);
	if(info>0) {
		printf("Matrix is singular (%d)\n",info);
		exit(1);
	}
	return lu;
	
}

void matrix::lu2(int *&ipiv) {

	int info;
	char trans='N';

	ipiv=(int *)realloc(ipiv,nf*sizeof(int));

	dgetrf_(&nf,&nc,p,&nf,ipiv,&info);
	if(info>0) {
		printf("Matrix is singular (%d)\n",info);
		exit(1);
	}
	
}

matrix matrix::lusolve(matrix b,int *ipiv) const {

	int info,nfil=nf;
	char trans='N';

	dgetrs_(&trans,&nfil,&b.nc,p,&nfil,ipiv,b.p,&b.nf,&info);

	return b;
	
}

matrix matrix::solve(matrix b) const {

	int ipiv[nf],info;
	char trans='N';
	matrix lu(*this);

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

