#include"matrix.h"
#include<stdlib.h>
extern "C" {
#include CBLAS
}

template <class Type>
Matrix_block_diag<Type>::Matrix_block_diag(int nblocks) {

	if(nblocks<0) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Number of blocks can't be negative\n");
		exit(1);
	}
	if(nblocks==0) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Number of blocks can't be zero\n");
		exit(1);
	}
	nb=nblocks;
	m=new Matrix<Type>[nb];
	
}

template <class Type>
Matrix_block_diag<Type>::~Matrix_block_diag() {

	delete [] m;
	
}

template <class Type>
Matrix_block_diag<Type>::Matrix_block_diag(const Matrix_block_diag<Type> &a) {
	
	int i;

	nb=a.nb;
	m=new Matrix<Type>[nb];
	for(i=0;i<nb;i++) m[i]=a.m[i];
	
}

template <class Type>
Matrix_block_diag<Type> &Matrix_block_diag<Type>::operator=(const Matrix_block_diag<Type> &a) {

	int i;

	if(nb!=a.nb) set_nblocks(a.nb);
	for(i=0;i<nb;i++) m[i]=a.m[i];
	
	return *this;
}

template <class Type>
Matrix_block_diag<Type> & Matrix_block_diag<Type>::set_nblocks(int nblocks) {

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
		m=new Matrix<Type>[nb];
	}
	return *this;
}

template <class Type>
Matrix_block_diag<Type>::operator Matrix<Type>() const{

	int n,i,j;
	
	Matrix<Type> a;
	a.zero(nrows(),ncols());

	i=0;j=0;
	for(n=0;n<nb;n++) {
		a.setblock(i,i+m[n].nf-1,j,j+m[n].nc-1,m[n]);
		i+=m[n].nf;j+=m[n].nc;
	}
	
	return a;
}

template <class Type>	
Matrix<Type> &Matrix_block_diag<Type>::block(int i) const {

	if(i<0) i+=nb;
	if(i<0||i>=nb) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Index exceeds number of blocks\n");
		exit(1);
	}
	return m[i];
	
}

template <class Type>
int Matrix_block_diag<Type>::nblocks() const {

	return nb;
}

template <class Type>
int Matrix_block_diag<Type>::nrows() const {

	int i,n=0;

	for(i=0;i<nb;i++) n+=m[i].nf;

	return n;
	
}

template <class Type>
int Matrix_block_diag<Type>::ncols() const {

	int i,n=0;

	for(i=0;i<nb;i++) n+=m[i].nc;

	return n;
	
}

template <>
Matrix<double> Matrix_block_diag<double>::operator,(const Matrix<double> &a) const {
	
	int i,n=0,n2=0;

	if(ncols()!=a.nf) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Dimensions must agree\n");
		exit(1);
	}
	Matrix<double> res(nrows(),a.nc);

	for(i=0;i<nb;i++) {
		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m[i].nf,a.nc,m[i].nc,1,m[i].p,m[i].nf,a.p+n2,a.nf,0,res.p+n,res.nf);
		n+=m[i].nf;n2+=m[i].nc;
	}
	
	return res;
}

template <>
Matrix<float> Matrix_block_diag<float>::operator,(const Matrix<float> &a) const {
	
	int i,n=0,n2=0;

	if(ncols()!=a.nf) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Dimensions must agree\n");
		exit(1);
	}
	Matrix<float> res(nrows(),a.nc);

	for(i=0;i<nb;i++) {
		cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m[i].nf,a.nc,m[i].nc,1,m[i].p,m[i].nf,a.p+n2,a.nf,0,res.p+n,res.nf);
		n+=m[i].nf;n2+=m[i].nc;
	}
	
	return res;
}

template <class Type>
Matrix_block_diag<Type> Matrix_block_diag<Type>::operator,(const Matrix_block_diag<Type> &a) const {

	int i;
	Matrix_block_diag<Type> res(nb);

	if(nb!=a.nb) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Number of blocks must agree\n");
		exit(1);
	}

	for(i=0;i<nb;i++) res.m[i]=(m[i],a.m[i]);

	return res;
	
}

template<>
Matrix<double> operator,(const Matrix<double> &a,const Matrix_block_diag<double> &b) {

	int i;
	unsigned n=0,n2=0;

	if(a.ncols()!=b.nrows()) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Dimensions must agree\n");
		exit(1);
	}
	Matrix<double> res(a.nrows(),b.ncols());

	for(i=0;i<b.nb;i++) {
		cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,a.nrows(),b.m[i].ncols(),b.m[i].nrows(),1,a.data()+n2,a.nrows(),b.m[i].data(),b.m[i].nrows(),0,res.data()+n,res.nrows());
		n+=b.m[i].ncols()*a.nrows();n2+=b.m[i].nrows()*a.nrows();
	}
	
	return res;

}

template<>
Matrix<float> operator,(const Matrix<float> &a,const Matrix_block_diag<float> &b) {

	int i;
	unsigned n=0,n2=0;

	if(a.ncols()!=b.nrows()) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Dimensions must agree\n");
		exit(1);
	}
	Matrix<float> res(a.nrows(),b.ncols());

	for(i=0;i<b.nb;i++) {
		cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,a.nrows(),b.m[i].ncols(),b.m[i].nrows(),1,a.data()+n2,a.nrows(),b.m[i].data(),b.m[i].nrows(),0,res.data()+n,res.nrows());
		n+=b.m[i].ncols()*a.nrows();n2+=b.m[i].nrows()*a.nrows();
	}
	
	return res;

}

template <class Type>
Matrix_block_diag<Type> Matrix_block_diag<Type>::operator*(const Matrix<Type> &z) const {

	int k,i,j;
	Matrix_block_diag<Type> res(nb);
	Matrix<Type> zi;
	
	if((z.nc!=ncols()&&z.nc!=1)||(z.nf!=nrows()&&z.nf!=1)) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Dimensions must agree\n");
		exit(1);
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

template <class Type>
Matrix_block_diag<Type> Matrix_block_diag<Type>::operator/(const Matrix<Type> &z) const {

	int k,i,j;
	Matrix_block_diag<Type> res(nb);
	Matrix<Type> zi;
	
	if((z.nc!=ncols()&&z.nc!=1)||(z.nf!=nrows()&&z.nf!=1)) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Dimensions must agree\n");
		exit(1);
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

template <class Type>
Matrix_block_diag<Type> Matrix_block_diag<Type>::operator*(mat_type n) const {

	int i;
	Matrix_block_diag<Type> res(nb);
	
	for(i=0;i<nb;i++) 
		res.m[i]=n*m[i];
	
	return res;
}

template <class Type>
Matrix_block_diag<Type> Matrix_block_diag<Type>::operator/(mat_type n) const {

	int i;
	Matrix_block_diag<Type> res(nb);
	
	for(i=0;i<nb;i++) 
		res.m[i]=m[i]/n;
	
	return res;
	
}

template <class Type>
Matrix_block_diag<Type> Matrix_block_diag<Type>::operator*(const Matrix_block_diag<Type> &a) const {

	Matrix_block_diag<Type> res(nb);
	int i;

	if(nb!=a.nb) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Number of blocks must agree\n");
		exit(1);
	}

	for(i=0;i<nb;i++) 
		res.m[i]=m[i]*a.m[i];
		
	return res;
	
}

template <class Type>
Matrix_block_diag<Type> Matrix_block_diag<Type>::operator+(const Matrix_block_diag<Type> &a) const {

	Matrix_block_diag<Type> res(nb);
	int i;

	if(nb!=a.nb) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Number of blocks must agree\n");
		exit(1);
	}

	for(i=0;i<nb;i++) 
		res.m[i]=m[i]+a.m[i];
		
	return res;
	
}

template <class Type>
Matrix_block_diag<Type> Matrix_block_diag<Type>::operator-(const Matrix_block_diag<Type> &a) const {

	Matrix_block_diag<Type> res(nb);
	int i;

	if(nb!=a.nb) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Number of blocks must agree\n");
		exit(1);
	}

	for(i=0;i<nb;i++) 
		res.m[i]=m[i]-a.m[i];
		
	return res;
	
}

template <class Type>
Matrix_block_diag<Type> operator-(const Matrix_block_diag<Type> &a) {
	
	Matrix_block_diag<Type> res(a.nb);
	int i;
	
	for(i=0;i<a.nb;i++) 
		res.m[i]=-a.m[i];
	
	return res;
}

template <class Type>
Matrix<Type> Matrix_block_diag<Type>::row(int n) const {

	int i=0,j,nf;
	Matrix<Type> res;
	Type *pres,*pp;
	
	nf=nrows();
	if(n<0) n+=nf;
	if(n<0||n>=nf) {
		fprintf(stderr,"ERROR: (matrix_block_diag.row) Index exceeds matrix dimensions\n");
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

template <class Type>
Type Matrix_block_diag<Type>::operator()(int nfil,int ncol) const {

	int n;
	
	n=nrows();
	if(nfil<0) nfil+=n;
	if(ncol<0) ncol+=n;
	if(nfil<0||nfil>=n||ncol<0||ncol>=n) {
		fprintf(stderr,"ERROR: (matrix_block_diag) Index exceeds matrix dimensions\n");
		exit(1);
	}
	
	int i=0,j=0;
	while(nfil>=m[i].nf) nfil-=m[i++].nf;
	while(ncol>=m[j].nc) ncol-=m[j++].nc;
	
	if(i!=j) return 0;
	
	return m[i](nfil,ncol);

}

template <class Type>
Matrix_block_diag<Type> eye(const Matrix_block_diag<Type> &a) {

	Matrix_block_diag<Type> I(a.nb);
	int i;
	
	for(i=0;i<a.nb;i++) I.m[i]=eye(a.m[i].nrows());
	
	return I;
}

template <class Type>
Matrix_block_diag<Type> Matrix_block_diag<Type>::transpose() const {

	Matrix_block_diag<Type> a(nb);
	int i;
	
	for(i=0;i<nb;i++)
		a.m[i]=m[i].transpose();
	
	return a;
	
}

// Explicit instantiations

template class Matrix_block_diag<double>;
template class Matrix_block_diag<float>;

template Matrix<double> operator,(const Matrix<double> &,const Matrix_block_diag<double> &);
template Matrix_block_diag<double> operator-(const Matrix_block_diag<double>&);
template Matrix_block_diag<double> eye(const Matrix_block_diag<double> &);

template Matrix<float> operator,(const Matrix<float> &,const Matrix_block_diag<float> &);
template Matrix_block_diag<float> operator-(const Matrix_block_diag<float>&);
template Matrix_block_diag<float> eye(const Matrix_block_diag<float> &);


