#include"matrix.h"
#include<stdlib.h>
extern "C" {
#include CBLAS
}

matrix_block_diag::matrix_block_diag(int nblocks) {

	if(nblocks<0) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Number of blocks can't be negative\n");
		exit(1);
	}
	if(nblocks==0) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Number of blocks can't be zero\n");
		exit(1);
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
		fprintf(stderr,"ERROR: (matrix_block_diag) Number of blocks can't be negative\n");
		exit(1);
	}
	if(nblocks==0) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Number of blocks can't be zero\n");
		exit(1);
	}

	if(nblocks!=nb) {
		nb=nblocks;
		delete [] m;
		m=new matrix[nb];
	}
	return *this;
}

matrix_block_diag::operator matrix() const{

	int n,i;
	
	n=nrows();
	
	matrix a(zeros(n,n));

	n=0;

	for(i=0;i<nb;i++) {
		a.setblock(n,n+m[i].nf-1,n,n+m[i].nf-1,m[i]);
		n+=m[i].nf;
	}
	
	return a;
}
	
matrix &matrix_block_diag::block(int i) const {

	if(i<0) i+=nb;
	if(i<0||i>=nb) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Index exceed number of blocks\n");
		exit(1);
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

	return nrows();
	
}

matrix matrix_block_diag::operator,(const matrix &a) const {
	
	int i,n=0;

	if(ncols()!=a.nf) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Dimensions must agree\n");
		exit(1);
	}
	matrix res(nrows(),a.nc);

	for(i=0;i<nb;i++) {
		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m[i].nf,a.nc,m[i].nc,1,m[i].p,m[i].nf,a.p+n,a.nf,0,res.p+n,res.nf);
		n+=m[i].nf;
	}
	
	return res;
}

matrix_block_diag matrix_block_diag::operator,(const matrix_block_diag &a) const {

	int i;
	matrix_block_diag res(nb);

	if(nb!=a.nb) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Number of blocks must agree\n");
		exit(1);
	}

	for(i=0;i<nb;i++) res.m[i]=(m[i],a.m[i]);

	return res;
	
}

matrix operator,(const matrix &a,const matrix_block_diag &b) {

	int i;
	unsigned n=0;

	if(a.ncols()!=b.nrows()) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Dimensions must agree\n");
		exit(1);
	}
	matrix res(a.nrows(),b.nrows());

	for(i=0;i<b.nb;i++) {
		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,a.nrows(),b.m[i].ncols(),b.m[i].nrows(),1,a.data()+n,a.nrows(),b.m[i].data(),b.m[i].nrows(),0,res.data()+n,res.nrows());
		n+=b.m[i].nrows()*a.nrows();
	}
	
	return res;

}

matrix_block_diag matrix_block_diag::operator*(const matrix &z) const {

	double *pz=z.p,*pp,*pres;
	int i,j,k;
	matrix_block_diag res(nb);

	if(z.nc!=1||z.nf!=nrows()) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Dimensions must agree\n");
		exit(1);
	}

	for(i=0;i<nb;i++) {
		pp=m[i].p;
		res.m[i].dim(m[i].nf,m[i].nc);
		pres=res.m[i].p;
		for(k=0;k<m[i].nc;k++) 
			for(j=0;j<m[i].nf;j++) 
				*(pres++)=*(pp++)*(*(pz+j));
		pz+=m[i].nf;
	}
	
	return res;
}

matrix_block_diag matrix_block_diag::operator/(const matrix &z) const {

	double *pz=z.p,*pp,*pres;
	int i,j,k;
	matrix_block_diag res(nb);
	
	if(z.nc!=1||z.nf!=nrows()) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Dimensions must agree\n");
		exit(1);
	}
	
	for(i=0;i<nb;i++) {
		pp=m[i].p;
		res.m[i].dim(m[i].nf,m[i].nc);
		pres=res.m[i].p;
		for(k=0;k<m[i].nc;k++) 
			for(j=0;j<m[i].nf;j++)
				*(pres++)=*(pp++)/(*(pz+j));
		pz+=m[i].nf;
	}
	
	return res;
}

matrix_block_diag operator*(const matrix &z,const matrix_block_diag &a) {

	return a*z;
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

matrix_block_diag operator*(double n,const matrix_block_diag &a) {

	return a*n;
}

matrix_block_diag matrix_block_diag::operator+(const matrix_block_diag &a) const {

	matrix_block_diag res(nb);
	int i;

	if(nb!=a.nb) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Number of blocks must agree\n");
		exit(1);
	}

	for(i=0;i<nb;i++) 
		res.m[i]=m[i]+a.m[i];
		
	return res;
	
}

matrix_block_diag matrix_block_diag::operator-(const matrix_block_diag &a) const {

	matrix_block_diag res(nb);
	int i;

	if(nb!=a.nb) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Number of blocks must agree\n");
		exit(1);
	}

	for(i=0;i<nb;i++) 
		res.m[i]=m[i]-a.m[i];
		
	return res;
	
}

matrix_block_diag operator+(const matrix_block_diag &a) {

	return a;
}

matrix_block_diag operator-(const matrix_block_diag &a) {
	
	matrix_block_diag res(a.nb);
	int i;
	
	for(i=0;i<a.nb;i++) 
		res.m[i]=-a.m[i];
	
	return res;
}

matrix matrix_block_diag::operator+(const matrix &a) const {

	return a+(*this);

}

matrix matrix_block_diag::operator-(const matrix &a) const {

	return a-(*this);

}

matrix matrix_block_diag::row(int n) const {

	int i=0,j,nf;
	matrix res;
	double *pres,*pp;
	
	nf=nrows();
	if(n<0) n+=nf;
	if(n<0||n>=nf) {
		fprintf(stderr,"ERROR: (matrix_block_diag.row) Index exceed matrix dimensions\n");
		exit(1);
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

	int i=0,n;
	
	n=nrows();
	if(nfil<0) nfil+=n;
	if(ncol<0) ncol+=n;
	if(nfil<0||nfil>=n||ncol<0||ncol>=n) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Index exceed matrix dimensions\n");
		exit(1);
	}
	
	while(nfil>=m[i].nf) {
		nfil-=m[i].nf;
		ncol-=m[i++].nc;
	}
	
	if(ncol<0||ncol>=m[i].nc) return 0;
	
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
