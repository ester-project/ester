#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "matrix.h"

extern "C" {
#include <stdlib.h>
#ifdef USE_MKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
}


matrix_block_diag::matrix_block_diag(int nblocks) {

	if(nblocks<0) {
		ester_critical("(matrix_block_diag) Number of blocks can't be negative");
	}
	if(nblocks==0) {
		ester_critical("(matrix_block_diag) Number of blocks can't be zero");
	}
	nb=nblocks;
	m=new matrix[nb];
	
}


matrix_block_diag::~matrix_block_diag() {

	delete [] m;
	
}


matrix_block_diag::matrix_block_diag(const matrix_block_diag &a) {
	
	int i;

	nb=a.nb;
	m=new matrix[nb];
	for(i=0;i<nb;i++) m[i]=a.m[i];
	
}


matrix_block_diag &matrix_block_diag::operator=(const matrix_block_diag &a) {

	int i;

	if(nb!=a.nb) set_nblocks(a.nb);
	for(i=0;i<nb;i++) m[i]=a.m[i];
	
	return *this;
}


matrix_block_diag & matrix_block_diag::set_nblocks(int nblocks) {

	if(nblocks<0) {
		ester_critical("(matrix_block_diag) Number of blocks can't be negative");
	}
	if(nblocks==0) {
		ester_critical("(matrix_block_diag) Number of blocks can't be zero");
	}

	if(nblocks!=nb) {
		nb=nblocks;
		delete [] m;
		m=new matrix[nb];
	}
	return *this;
}

matrix matrix_block_diag::full_matrix() const {
    return (matrix) (*this);
}

matrix_block_diag::operator matrix() const{

	int n,i,j;
	
	matrix a;
	a.zero(nrows(),ncols());

	i=0;j=0;
	for(n=0;n<nb;n++) {
		a.setblock(i,i+m[n].nf-1,j,j+m[n].nc-1,m[n]);
		i+=m[n].nf;j+=m[n].nc;
	}
	
	return a;
}

	
matrix &matrix_block_diag::block(int i){

	if(i<0) i+=nb;
	if(i<0||i>=nb) {
		ester_critical("(matrix_block_diag) Index exceeds number of blocks");
	}
	return m[i];
	
}

const matrix &matrix_block_diag::block(int i) const {

	if(i<0) i+=nb;
	if(i<0||i>=nb) {
		ester_critical("(matrix_block_diag) Index exceeds number of blocks");
	}
	return m[i];
	
}


int matrix_block_diag::nblocks() const {

	return nb;
}


int matrix_block_diag::nrows() const {

	int i,n=0;

	for(i=0;i<nb;i++) n+=m[i].nf;

	return n;
	
}


int matrix_block_diag::ncols() const {

	int i,n=0;

	for(i=0;i<nb;i++) n+=m[i].nc;

	return n;
	
}

matrix matrix_block_diag::operator,(const matrix &a) const {
	
	int i,n=0,n2=0;

	if(ncols()!=a.nf) {
		ester_critical("(matrix_block_diag) Dimensions must agree");
	}
	matrix res(nrows(),a.nc);

	for(i=0;i<nb;i++) {
		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m[i].nf,a.nc,m[i].nc,1,m[i].p,m[i].nf,a.p+n2,a.nf,0,res.p+n,res.nf);
		n+=m[i].nf;n2+=m[i].nc;
	}
	
	return res;
}

matrix_block_diag matrix_block_diag::operator,(const matrix_block_diag &a) const {

	int i;
	matrix_block_diag res(nb);

	if(nb!=a.nb) {
		ester_critical("(matrix_block_diag) Number of blocks must agree");
	}

	for(i=0;i<nb;i++) res.m[i]=(m[i],a.m[i]);

	return res;
	
}

matrix operator,(const matrix &a,const matrix_block_diag &b) {

	int i;
	unsigned n=0,n2=0;

	if(a.ncols()!=b.nrows()) {
		ester_critical("(matrix_block_diag) Dimensions must agree");
	}
	matrix res(a.nrows(),b.ncols());

	for(i=0;i<b.nb;i++) {
		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,a.nrows(),b.m[i].ncols(),b.m[i].nrows(),1,a.data()+n2,a.nrows(),b.m[i].data(),b.m[i].nrows(),0,res.data()+n,res.nrows());
		n+=b.m[i].ncols()*a.nrows();n2+=b.m[i].nrows()*a.nrows();
	}
	
	return res;

}

matrix_block_diag matrix_block_diag::operator*(const matrix &z) const {

	int k,i,j;
	matrix_block_diag res(nb);
	matrix zi;
	
	if((z.nc!=ncols()&&z.nc!=1)||(z.nf!=nrows()&&z.nf!=1)) {
		ester_critical("(matrix_block_diag) Dimensions must agree");
	}
	
	i=0;j=0;
	for(k=0;k<nb;k++) {
		if(z.nc==1&&z.nf==1) zi=z;
		else if(z.nc==1) zi=z.block(i,i+m[k].nf-1,0,0);
		else if(z.nf==1) zi=z.block(0,0,j,j+m[k].nc-1);
		else zi=z.block(i,i+m[k].nf-1,j,j+m[k].nc-1);
		res.m[k]=m[k]*zi;
		i+=m[k].nf;j+=m[k].nc;
	}
	
	return res;
}


matrix_block_diag matrix_block_diag::operator/(const matrix &z) const {

	int k,i,j;
	matrix_block_diag res(nb);
	matrix zi;
	
	if((z.nc!=ncols()&&z.nc!=1)||(z.nf!=nrows()&&z.nf!=1)) {
		ester_critical("(matrix_block_diag) Dimensions must agree");
	}
	
	i=0;j=0;
	for(k=0;k<nb;k++) {
		if(z.nc==1&&z.nf==1) zi=z;
		else if(z.nc==1) zi=z.block(i,i+m[k].nf-1,0,0);
		else if(z.nf==1) zi=z.block(0,0,j,j+m[k].nc-1);
		else zi=z.block(i,i+m[k].nf-1,j,j+m[k].nc-1);
		res.m[k]=m[k]/zi;
		i+=m[k].nf;j+=m[k].nc;
	}
	
	return res;
}


matrix_block_diag matrix_block_diag::operator*(double n) const {

	int i;
	matrix_block_diag res(nb);
	
	for(i=0;i<nb;i++) 
		res.m[i]=n*m[i];
	
	return res;
}


matrix_block_diag matrix_block_diag::operator/(double n) const {

	int i;
	matrix_block_diag res(nb);
	
	for(i=0;i<nb;i++) 
		res.m[i]=m[i]/n;
	
	return res;
	
}


matrix_block_diag matrix_block_diag::operator*(const matrix_block_diag &a) const {

	matrix_block_diag res(nb);
	int i;

	if(nb!=a.nb) {
		ester_critical("(matrix_block_diag) Number of blocks must agree");
	}

	for(i=0;i<nb;i++) 
		res.m[i]=m[i]*a.m[i];
		
	return res;
	
}


matrix_block_diag matrix_block_diag::operator+(const matrix_block_diag &a) const {

	matrix_block_diag res(nb);
	int i;

	if(nb!=a.nb) {
		ester_critical("(matrix_block_diag) Number of blocks must agree");
	}

	for(i=0;i<nb;i++) 
		res.m[i]=m[i]+a.m[i];
		
	return res;
	
}


matrix_block_diag matrix_block_diag::operator-(const matrix_block_diag &a) const {

	matrix_block_diag res(nb);
	int i;

	if(nb!=a.nb) {
		ester_critical("(matrix_block_diag) Number of blocks must agree");
	}

	for(i=0;i<nb;i++) 
		res.m[i]=m[i]-a.m[i];
		
	return res;
	
}


matrix_block_diag operator-(const matrix_block_diag &a) {
	
	matrix_block_diag res(a.nb);
	int i;
	
	for(i=0;i<a.nb;i++) 
		res.m[i]=-a.m[i];
	
	return res;
}


matrix matrix_block_diag::row(int n) const {

	int i=0,j,nf;
	matrix res;
	double *pres,*pp;
	
	nf=nrows();
	if(n<0) n+=nf;
	if(n<0||n>=nf) {
		ester_critical("(matrix_block_diag.row) Index exceeds matrix dimensions");
	}
	res=zeros(1,nf);
	pres=res.p;
	while(n>=m[i].nf) {
		pres+=m[i].nf;
		n-=m[i++].nf;
	}
	pp=m[i].p+n;
	for(j=0;j<m[i].nc;j++) {
		*(pres++)=*pp;
		pp+=m[i].nf;
	}
	
	return res;

}


double matrix_block_diag::operator()(int nfil,int ncol) const {

	int n;
	
	n=nrows();
	if(nfil<0) nfil+=n;
	if(ncol<0) ncol+=n;
	if(nfil<0||nfil>=n||ncol<0||ncol>=n) {
		ester_critical("(matrix_block_diag) Index exceeds matrix dimensions");
	}
	
	int i=0,j=0;
	while(nfil>=m[i].nf) nfil-=m[i++].nf;
	while(ncol>=m[j].nc) ncol-=m[j++].nc;
	
	if(i!=j) return 0;
	
	return m[i](nfil,ncol);

}


matrix_block_diag eye(const matrix_block_diag &a) {

	matrix_block_diag I(a.nb);
	int i;
	
	for(i=0;i<a.nb;i++) I.m[i]=eye(a.m[i].nrows());
	
	return I;
}


matrix_block_diag matrix_block_diag::transpose() const {

	matrix_block_diag a(nb);
	int i;
	
	for(i=0;i<nb;i++)
		a.m[i]=m[i].transpose();
	
	return a;
	
}




