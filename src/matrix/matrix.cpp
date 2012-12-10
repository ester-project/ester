#include<stdlib.h>
#include<cmath>
#include<stdio.h>
#include<string.h>
#include"matrix.h"

template <class Type>
Matrix<Type>::Matrix(int nfil,int ncol) {
	
	unsigned tam;

    if(nfil<0||ncol<0) {
		fprintf(stderr,"ERROR: Can't create matrix with negative size\n");
		exit(1);
	}
	if(nfil==0) {
		fprintf(stderr,"ERROR: Number of rows can't be zero\n");
		exit(1);
	}
	if(ncol==0) {
		fprintf(stderr,"ERROR: Number of columns can't be zero\n");
		exit(1);
	}
	nf=nfil;
    nc=ncol;
    tam=unsigned(nf)*unsigned(nc);
	p=new Type[tam];
}

template <class Type>
Matrix<Type>::~Matrix() {

	delete [] p;
}

template <class Type>
Matrix<Type>::Matrix(const Matrix<Type> &a) {

    unsigned tam;

    nc=a.nc;
    nf=a.nf;
    tam=unsigned(nf)*unsigned(nc);
	p=new Type[tam];

	memcpy(p,a.p,nc*nf*sizeof(Type));

}

template <class Type>
int Matrix<Type>::nrows() const {

	return nf;
}

template <class Type>
int Matrix<Type>::ncols() const {

	return nc;
}

template <class Type>
Type *Matrix<Type>::data() const {

	return p;
}

template <class Type>
Matrix<Type> &Matrix<Type>::dim(int nfil,int ncol) {
    
    unsigned tam;
    
    if(nfil<0||ncol<0) {
		fprintf(stderr,"ERROR: Can't create Matrix<Type> with negative size\n");
		exit(1);
	}
	if(nfil==0) {
		fprintf(stderr,"ERROR: Number of rows can't be zero\n");
		exit(1);
	}
	if(ncol==0) {
		fprintf(stderr,"ERROR: Number of columns can't be zero\n");
		exit(1);
	}
    
    if(nfil*ncol!=nf*nc) {
    	delete [] p;
    	nf=nfil;
    	nc=ncol;
    	tam=unsigned(nf)*unsigned(nc);
		p=new Type[tam];
    } else {
    	nf=nfil;
    	nc=ncol;
    }
    
    return *this;
}

template <class Type>
Matrix<Type> &Matrix<Type>::redim(int nfil,int ncol) {
    
    int tam;
    
    
    if(nfil<0||ncol<0) {
		fprintf(stderr,"ERROR: Can't create matrix with negative size\n");
		exit(1);
	}
	if(nfil==0) {
		fprintf(stderr,"ERROR: Number of rows can't be zero\n");
		exit(1);
	}
	if(ncol==0) {
		fprintf(stderr,"ERROR: Number of columns can't be zero\n");
		exit(1);
	}
    
    if(nfil*ncol!=nf*nc) {
    	fprintf(stderr,"ERROR: Matrix<Type>.redim: Number of elements doesn't match\n");
    	exit(1);
    }    
    nf=nfil;nc=ncol;
    
    return *this;
}

template <class Type>
Matrix<Type> &Matrix<Type>::operator=(const Matrix<Type> &a) {
        
    if(&a==this) return *this;
	if (nf!=a.nf||nc!=a.nc)
    	dim(a.nf,a.nc);

	memcpy(p,a.p,nc*nf*sizeof(Type));
    
    return *this;
}
template <class Type>
template <class Type2>
Matrix<Type> &Matrix<Type>::operator=(const Matrix<Type2> &a) {
    
	if (nf!=a.nf||nc!=a.nc)
    	dim(a.nf,a.nc);

    int N;
    N=a.nf*a.nc;
	for(int i=0;i<N;i++)
		p[i]=a.p[i];
    
    return *this;
}

template Matrix<double> &Matrix<double>::operator=(const Matrix<float> &a);
template Matrix<float> &Matrix<float>::operator=(const Matrix<double> &a);

template <class Type>
template <class Type2>
Matrix<Type>::operator Matrix<Type2>() const {

	Matrix<Type2> a;
	int N;
	N=nf*nc;
	for(int i=0;i<N;i++)
		a.p[i]=p[i];
		
	return a;
}

template Matrix<float>::operator Matrix<double>() const;
template Matrix<double>::operator Matrix<float>() const;

template <class Type>
Type &Matrix<Type>::operator()(int ifil,int icol) {

	if(ifil<0) ifil+=nf;
	if(icol<0) icol+=nc;
	if(ifil>=nf||ifil<0||icol>=nc||icol<0) {
		fprintf(stderr,"ERROR: Index exceeds matrix dimensions\n");
		exit(1);
	}
	return *(p+icol*nf+ifil);
	
}

template <class Type>
Type &Matrix<Type>::operator()(int ielem) {
	
	if(ielem<0) ielem+=nf*nc;
	if(ielem>=nf*nc||ielem<0) {
		fprintf(stderr,"ERROR: Index exceeds matrix dimensions\n");
		exit(1);
	}
	return *(p+ielem);
	
}

template <class Type>
Type Matrix<Type>::operator()(int ifil,int icol) const {

	if(ifil<0) ifil+=nf;
	if(icol<0) icol+=nc;
	if(ifil>=nf||ifil<0||icol>=nc||icol<0) {
		fprintf(stderr,"ERROR: Index exceeds matrix dimensions\n");
		exit(1);
	}
	return *(p+icol*nf+ifil);
	
}

template <class Type>
Type Matrix<Type>::operator()(int ielem) const {
	
	if(ielem<0) ielem+=nf*nc;
	if(ielem>=nf*nc||ielem<0) {
		fprintf(stderr,"ERROR: Index exceeds matrix dimensions\n");
		exit(1);
	}
	return *(p+ielem);
	
}

template <>
int Matrix<double>::read(int nfil,int ncol,FILE *fp,char mode) {
    
    int i,tam;

    dim(nfil,ncol);
    tam=nc*nf;           

    if (mode=='b') {
        if((int)fread(p,sizeof(double),tam,fp)!=tam) return 0;
        }
    else
        for(i=0;i<tam;i++) {
            if(fscanf(fp,"%le",&p[i])==EOF) return 0;
            }
    return 1;
}

template <>
int Matrix<float>::read(int nfil,int ncol,FILE *fp,char mode) {
    
    int i,tam;

    dim(nfil,ncol);
    tam=nc*nf;           

    if (mode=='b') {
        if((int)fread(p,sizeof(float),tam,fp)!=tam) return 0;
        }
    else
        for(i=0;i<tam;i++) {
            if(fscanf(fp,"%e",&p[i])==EOF) return 0;
            }
    return 1;
}

template <class Type>
int Matrix<Type>::write(FILE *fp,char mode) const {
   
    int tam,i,j;

    tam=nf*nc;
    if (mode=='b') {
        if((int)fwrite(p,sizeof(Type),tam,fp)!=tam) return 0;
        }
    else {
        j=0;
        for(i=0;i<tam;i++) {
            if(fprintf(fp,"%.16e ",p[i])==EOF) return 0;
            j++;
            if (j==nf) {
                j=0;
                fprintf(fp,"\n");
            }
        }
    }
    return 1;
}

template <class Type>
void Matrix<Type>::write_fmt(const char *fmt,FILE *fp) const {

	int i,j;
	
	for(i=0;i<nf;i++) {
		for(j=0;j<nc;j++) {
			fprintf(fp,fmt,(*this)(i,j));
			if(j<nc-1) fprintf(fp," ");
		}
		fprintf(fp,"\n");
	}
}

template <class Type>
void Matrix<Type>::swap(Matrix<Type> &a) {

	Type *p0;
	int nf0,nc0;
	
	p0=p;
	nf0=nf;nc0=nc;
	
	nf=a.nf;nc=a.nc;
	p=a.p;
	a.nf=nf0;a.nc=nc0;
	a.p=p0;

}

template <class Type>
void Matrix<Type>::zero(int nrows,int ncols) {

	dim(nrows,ncols);
	for(int i=0;i<nf*nc;i++) p[i]=0;
	
}

template <class Type>
Matrix<Type> Matrix<Type>::operator+(const Matrix<Type> &a) const {

    Matrix<Type> res;
	int i,j,resnf,resnc,N;
	
	if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
		fprintf(stderr,"ERROR: (matrix.+) Dimensions must agree\n");
		exit(1);
	}
	
	if(nf>a.nf) resnf=nf;
		else resnf=a.nf;
	if(nc>a.nc) resnc=nc;
		else resnc=a.nc;
	
	res.dim(resnf,resnc);
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    res.p[i]=p[i]+a.p[i];
	} else {
		Type *pi,*pa,*pres;
		pi=p;pa=a.p;pres=res.p;
		for(i=0;i<resnc;i++) {
    	   	for(j=0;j<resnf;j++) {
    	   	    *(pres++)=*(pi)+(*(pa));
    	   	    if(nf>1) pi++;
    	   	    if(a.nf>1) pa++;
    	   	}
    	   	if(nc==1) pi=p;
    	   		else if(nf==1) pi++;
    	   	if(a.nc==1) pa=a.p;
    	   		else if(a.nf==1) pa++;
    	}
    }
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator+(mat_type n) const {

    Matrix<Type> res(nf,nc);
    int i,N;
    
    N=nc*nf;
    for(i=0;i<N;i++)
    	res.p[i]=p[i]+n;
    	
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator-(const Matrix<Type> &a) const {

    Matrix<Type> res;
	int i,j,resnf,resnc,N;
	
	if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
		fprintf(stderr,"ERROR: (matrix.-) Dimensions must agree\n");
		exit(1);
	}
	
	if(nf>a.nf) resnf=nf;
		else resnf=a.nf;
	if(nc>a.nc) resnc=nc;
		else resnc=a.nc;
	
	res.dim(resnf,resnc);
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    res.p[i]=p[i]-a.p[i];
	} else {
		Type *pi,*pa,*pres;
		pi=p;pa=a.p;pres=res.p;
		for(i=0;i<resnc;i++) {
    	   	for(j=0;j<resnf;j++) {
    	   	    *(pres++)=*(pi)-(*(pa));
    	   	    if(nf>1) pi++;
    	   	    if(a.nf>1) pa++;
    	   	}
    	   	if(nc==1) pi=p;
    	   		else if(nf==1) pi++;
    	   	if(a.nc==1) pa=a.p;
    	   		else if(a.nf==1) pa++;
    	}
    }
    return res;
}

template <class Type>
Matrix<Type> operator-(mat_type n,const Matrix<Type> &a) {
	
	Matrix<Type> res(a.nf,a.nc);
    int i,N;
    
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	res.p[i]=n-a.p[i];
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator-(mat_type n) const {

    Matrix<Type> res(nf,nc);
    int i,N;
    
    N=nc*nf;
    for(i=0;i<N;i++)
    	res.p[i]=p[i]-n;
    	
    return res;

}

template <class Type>
Matrix<Type> operator-(const Matrix<Type> &a) {
    
    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	res.p[i]=-a.p[i];
    	
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator*(const Matrix<Type> &a) const {

    Matrix<Type> res;
	int i,j,resnf,resnc,N;
	
	if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
		fprintf(stderr,"ERROR: (matrix.*) Dimensions must agree\n");
		exit(1);
	}
	
	if(nf>a.nf) resnf=nf;
		else resnf=a.nf;
	if(nc>a.nc) resnc=nc;
		else resnc=a.nc;
	
	res.dim(resnf,resnc);
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nf*nc;
		for(i=0;i<N;i++)
       	    res.p[i]=p[i]*a.p[i];
	} else {
		Type *pi,*pa,*pres;
		pi=p;pa=a.p;pres=res.p;
		for(i=0;i<resnc;i++) {
        	for(j=0;j<resnf;j++) {
       		    *(pres++)=*(pi)*(*(pa));
       		    if(nf>1) pi++;
       		    if(a.nf>1) pa++;
       		}
       		if(nc==1) pi=p;
       			else if(nf==1) pi++;
       		if(a.nc==1) pa=a.p;
       			else if(a.nf==1) pa++;
    	}
    }
    
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator*(mat_type n) const {
    
    Matrix<Type> res(nf,nc);
    int i,N;
    
    N=nc*nf;
    for(i=0;i<N;i++)
    	res.p[i]=p[i]*n;
    	
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator/(const Matrix<Type> &a) const {

    Matrix<Type> res;
	int i,j,resnf,resnc,N;
	
	if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
		fprintf(stderr,"ERROR: (matrix./) Dimensions must agree\n");
		exit(1);
	}
	
	if(nf>a.nf) resnf=nf;
		else resnf=a.nf;
	if(nc>a.nc) resnc=nc;
		else resnc=a.nc;
	
	res.dim(resnf,resnc);
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    res.p[i]=p[i]/a.p[i];
	} else {
		Type *pi,*pa,*pres;
		pi=p;pa=a.p;pres=res.p;
		for(i=0;i<resnc;i++) {
    	   	for(j=0;j<resnf;j++) {
    	   	    *(pres++)=*(pi)/(*(pa));
    	   	    if(nf>1) pi++;
    	   	    if(a.nf>1) pa++;
    	   	}
    	   	if(nc==1) pi=p;
    	   		else if(nf==1) pi++;
    	   	if(a.nc==1) pa=a.p;
    	   		else if(a.nf==1) pa++;
    	}
    }
    return res;
}

template <class Type>
Matrix<Type> operator/(mat_type n,const Matrix<Type> &a) {
	
	Matrix<Type> res(a.nf,a.nc);
    int i,N;

    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	res.p[i]=n/a.p[i];	
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator/(mat_type n) const {

    Matrix<Type> res(nf,nc);
    int i,N;
    
    N=nc*nf;
    for(i=0;i<N;i++)
    	res.p[i]=p[i]/n;
    	
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator==(const Matrix<Type> &a) const {

    Matrix<Type> res;
	int i,j,resnf,resnc,N;
	
	if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
		fprintf(stderr,"ERROR: (matrix.==) Dimensions must agree\n");
		exit(1);
	}
	
	if(nf>a.nf) resnf=nf;
		else resnf=a.nf;
	if(nc>a.nc) resnc=nc;
		else resnc=a.nc;
	
	res.dim(resnf,resnc);
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    res.p[i]=(p[i]==a.p[i]);
	} else {
		Type *pi,*pa,*pres;
		pi=p;pa=a.p;pres=res.p;
		for(i=0;i<resnc;i++) {
    	   	for(j=0;j<resnf;j++) {
    	   	    *(pres++)=(*(pi)==(*(pa)));
    	   	    if(nf>1) pi++;
    	   	    if(a.nf>1) pa++;
    	   	}
    	   	if(nc==1) pi=p;
    	   		else if(nf==1) pi++;
    	   	if(a.nc==1) pa=a.p;
    	   		else if(a.nf==1) pa++;
    	}
    }
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator==(mat_type n) const {
	
	Matrix<Type> res(nf,nc);
    int i,N;
    
    N=nc*nf;
    for(i=0;i<N;i++)
    	res.p[i]=(p[i]==n);
    	
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator!=(const Matrix<Type> &a) const {

    Matrix<Type> res;
	int i,j,resnf,resnc,N;
	
	if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
		fprintf(stderr,"ERROR: (matrix.!=) Dimensions must agree\n");
		exit(1);
	}
	
	if(nf>a.nf) resnf=nf;
		else resnf=a.nf;
	if(nc>a.nc) resnc=nc;
		else resnc=a.nc;
	
	res.dim(resnf,resnc);
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    res.p[i]=(p[i]!=a.p[i]);
	} else {
		Type *pi,*pa,*pres;
		pi=p;pa=a.p;pres=res.p;
		for(i=0;i<resnc;i++) {
    	   	for(j=0;j<resnf;j++) {
    	   	    *(pres++)=(*(pi)!=(*(pa)));
    	   	    if(nf>1) pi++;
    	   	    if(a.nf>1) pa++;
    	   	}
    	   	if(nc==1) pi=p;
    	   		else if(nf==1) pi++;
    	   	if(a.nc==1) pa=a.p;
    	   		else if(a.nf==1) pa++;
    	}
    }
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator!=(mat_type n) const {
	
	Matrix<Type> res(nf,nc);
    int i,N;
    
    N=nc*nf;
    for(i=0;i<N;i++)
    	res.p[i]=(p[i]!=n);
    	
    return res;}

template <class Type>
Matrix<Type> Matrix<Type>::operator>(const Matrix<Type> &a) const{

    Matrix<Type> res;
	int i,j,resnf,resnc,N;
	
	if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
		fprintf(stderr,"ERROR: (matrix.>) Dimensions must agree\n");
		exit(1);
	}
	
	if(nf>a.nf) resnf=nf;
		else resnf=a.nf;
	if(nc>a.nc) resnc=nc;
		else resnc=a.nc;
	
	res.dim(resnf,resnc);
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
			res.p[i]=(p[i]>a.p[i]);
	} else {
		Type *pi,*pa,*pres;
		pi=p;pa=a.p;pres=res.p;
		for(i=0;i<resnc;i++) {
    	   	for(j=0;j<resnf;j++) {
    	   	    *(pres++)=(*(pi)>(*(pa)));
    	   	    if(nf>1) pi++;
    	   	    if(a.nf>1) pa++;
    	   	}
    	   	if(nc==1) pi=p;
    	   		else if(nf==1) pi++;
    	   	if(a.nc==1) pa=a.p;
    	   		else if(a.nf==1) pa++;
    	}
    }
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator>(mat_type n) const {
	
	Matrix<Type> res(nf,nc);
    int i,N;
    
    N=nc*nf;
    for(i=0;i<N;i++)
    	res.p[i]=(p[i]>n);
    	
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator<(mat_type n) const {
	
	Matrix<Type> res(nf,nc);
    int i,N;
    
    N=nc*nf;
    for(i=0;i<N;i++)
    	res.p[i]=(p[i]<n);
    	
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator>=(const Matrix<Type> &a) const {

    Matrix<Type> res;
	int i,j,resnf,resnc,N;
	
	if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
		fprintf(stderr,"ERROR: (matrix.>=) Dimensions must agree\n");
		exit(1);
	}
	
	if(nf>a.nf) resnf=nf;
		else resnf=a.nf;
	if(nc>a.nc) resnc=nc;
		else resnc=a.nc;
	
	res.dim(resnf,resnc);
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    res.p[i]=(p[i]>=a.p[i]);
	} else {
		Type *pi,*pa,*pres;
		pi=p;pa=a.p;pres=res.p;
		for(i=0;i<resnc;i++) {
    	   	for(j=0;j<resnf;j++) {
    	   	    *(pres++)=(*(pi)>=(*(pa)));
    	   	    if(nf>1) pi++;
    	   	    if(a.nf>1) pa++;
    	   	}
    	   	if(nc==1) pi=p;
    	   		else if(nf==1) pi++;
    	   	if(a.nc==1) pa=a.p;
    	   		else if(a.nf==1) pa++;
    	}
    }
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator>=(mat_type n) const {
	
	Matrix<Type> res(nf,nc);
    int i,N;
    
    N=nc*nf;
    for(i=0;i<N;i++)
    	res.p[i]=(p[i]>=n);
    	
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator<=(mat_type n) const {
	
	Matrix<Type> res(nf,nc);
    int i,N;
    
    N=nc*nf;
    for(i=0;i<N;i++)
    	res.p[i]=(p[i]<=n);
    	
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator||(const Matrix<Type> &a) const {

    Matrix<Type> res;
	int i,j,resnf,resnc,N;
	
	if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
		fprintf(stderr,"ERROR: (matrix.||) Dimensions must agree\n");
		exit(1);
	}
	
	if(nf>a.nf) resnf=nf;
		else resnf=a.nf;
	if(nc>a.nc) resnc=nc;
		else resnc=a.nc;
	
	res.dim(resnf,resnc);
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    res.p[i]=(p[i]||a.p[i]);
	} else {
		Type *pi,*pa,*pres;
		pi=p;pa=a.p;pres=res.p;
		for(i=0;i<resnc;i++) {
    	   	for(j=0;j<resnf;j++) {
    	   	    *(pres++)=(*(pi)||(*(pa)));
    	   	    if(nf>1) pi++;
    	   	    if(a.nf>1) pa++;
    	   	}
    	   	if(nc==1) pi=p;
    	   		else if(nf==1) pi++;
    	   	if(a.nc==1) pa=a.p;
    	   		else if(a.nf==1) pa++;
    	}
    }
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator||(mat_type n) const {
	
	Matrix<Type> res(nf,nc);
    int i,N;
    
    N=nc*nf;
    for(i=0;i<N;i++)
    	res.p[i]=(p[i]||n);
    	
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator&&(const Matrix<Type> &a) const {

    Matrix<Type> res;
	int i,j,resnf,resnc,N;
	
	if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
		fprintf(stderr,"ERROR: (matrix.&&) Dimensions must agree\n");
		exit(1);
	}
	
	if(nf>a.nf) resnf=nf;
		else resnf=a.nf;
	if(nc>a.nc) resnc=nc;
		else resnc=a.nc;
	
	res.dim(resnf,resnc);
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    res.p[i]=(p[i]&&a.p[i]);
	} else {
		Type *pi,*pa,*pres;
		pi=p;pa=a.p;pres=res.p;
		for(i=0;i<resnc;i++) {
    	   	for(j=0;j<resnf;j++) {
    	   	    *(pres++)=(*(pi)&&(*(pa)));
    	   	    if(nf>1) pi++;
    	   	    if(a.nf>1) pa++;
    	   	}
    	   	if(nc==1) pi=p;
    	   		else if(nf==1) pi++;
    	   	if(a.nc==1) pa=a.p;
    	   		else if(a.nf==1) pa++;
    	}
    }
    return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::operator&&(mat_type n) const {
	
	Matrix<Type> res(nf,nc);
    int i,N;
    
    N=nc*nf;
    for(i=0;i<N;i++)
    	res.p[i]=(p[i]&&n);
    	
    return res;
}

template <class Type>
Matrix<Type> &Matrix<Type>::operator+=(const Matrix<Type> &a) {
    
    int i,j,k,N;

	N=nc*nf;
	if(nf==a.nf&&nc==a.nc) {
    	for(i=0;i<N;i++) 
    		p[i]+=a.p[i];
	} else if(nf==a.nf&&a.nc==1) {
		for(i=0,j=0;i<N;i++,j=(++j)%nf)
			p[i]+=a.p[j];
	} else if(a.nf==1&&nc==a.nc) {
		k=0;
		for(j=0;j<nc;j++)
			for(i=0;i<nf;i++,k++)
				p[k]+=a.p[j];
	} else if(a.nf==1&&a.nc==1) {
		*this+=*a.p;
	} else {
		*this=(*this)+a;
	}
    
    return *this;
}

template <class Type>
Matrix<Type> &Matrix<Type>::operator-=(const Matrix<Type> &a) {
    
    int i,j,k,N;

	N=nc*nf;
	if(nf==a.nf&&nc==a.nc) {
    	for(i=0;i<N;i++) 
    		p[i]-=a.p[i];
	} else if(nf==a.nf&&a.nc==1) {
		for(i=0,j=0;i<N;i++,j=(++j)%nf)
			p[i]-=a.p[j];
	} else if(a.nf==1&&nc==a.nc) {
		k=0;
		for(j=0;j<nc;j++)
			for(i=0;i<nf;i++,k++)
				p[k]-=a.p[j];
	} else if(a.nf==1&&a.nc==1) {
		*this-=*a.p;
	} else {
		*this=(*this)-a;
	}
    
    return *this;
}

template <class Type>
Matrix<Type> &Matrix<Type>::operator*=(const Matrix<Type> &a) {
    
    int i,j,k,N;

	N=nc*nf;
	if(nf==a.nf&&nc==a.nc) {
    	for(i=0;i<N;i++) 
    		p[i]*=a.p[i];
	} else if(nf==a.nf&&a.nc==1) {
		for(i=0,j=0;i<N;i++,j=(++j)%nf)
			p[i]*=a.p[j];
	} else if(a.nf==1&&nc==a.nc) {
		k=0;
		for(j=0;j<nc;j++)
			for(i=0;i<nf;i++,k++)
				p[k]*=a.p[j];
	} else if(a.nf==1&&a.nc==1) {
		*this*=*a.p;
	} else {
		*this=(*this)*a;
	}

    return *this;
}

template <class Type>
Matrix<Type> &Matrix<Type>::operator/=(const Matrix<Type> &a) {
    
    int i,j,k,N;

	N=nc*nf;
	if(nf==a.nf&&nc==a.nc) {
    	for(i=0;i<N;i++) 
    		p[i]/=a.p[i];
	} else if(nf==a.nf&&a.nc==1) {
		for(i=0,j=0;i<N;i++,j=(++j)%nf)
			p[i]/=a.p[j];
	} else if(a.nf==1&&nc==a.nc) {
		k=0;
		for(j=0;j<nc;j++)
			for(i=0;i<nf;i++,k++)
				p[k]/=a.p[j];
	} else if(a.nf==1&&a.nc==1) {
		*this/=*a.p;
	} else {
		*this=(*this)/a;
	}
    
    return *this;
}

template <class Type>
Matrix<Type> &Matrix<Type>::operator+=(mat_type n) {
    
	int i,N;
	N=nc*nf;
	for(i=0;i<N;i++) 
		p[i]+=n;
	return *this;
}

template <class Type>
Matrix<Type> &Matrix<Type>::operator-=(mat_type n) {
    
    int i,N;
	N=nc*nf;
	for(i=0;i<N;i++) 
		p[i]-=n;
	return *this;
}

template <class Type>
Matrix<Type> &Matrix<Type>::operator*=(mat_type n) {

	int i,N;
	N=nc*nf;
	for(i=0;i<N;i++) 
		p[i]*=n;
	return *this;
}

template <class Type>
Matrix<Type> &Matrix<Type>::operator/=(mat_type n) {

	int i,N;
	N=nc*nf;
	for(i=0;i<N;i++) 
		p[i]/=n;
	return *this;
}

template <class Type>
Matrix<Type> Matrix<Type>::row(int ifil) const {

	Matrix<Type> res(1,nc);
	Type *pi;
	int i;
	
	if(ifil<0) ifil+=nf;
	if(ifil<0||ifil>=nf) {
		fprintf(stderr,"ERROR: (matrix.row) Index exceeds matrix dimensions\n");
		exit(1);
	}
			
	pi=p+ifil;
	for(i=0;i<nc;i++)
		res.p[i]=pi[i*nf];
	
	return res;
}

template <class Type>
Matrix<Type> &Matrix<Type>::setrow(int ifil,const Matrix<Type> &a) {

	Type *pi;
	int i;
	
	if(ifil<0) ifil+=nf;
	if(ifil<0||ifil>=nf) {
		fprintf(stderr,"ERROR: (matrix.setrow) Index exceeds matrix dimensions\n");
		exit(1);
	}
	if(a.nf>1||a.nc!=nc) {
		fprintf(stderr,"ERROR: (matrix.setrow) Dimensions must agree\n");
		exit(1);
	}	
	
	pi=p+ifil;
	for(i=0;i<nc;i++)
		pi[i*nf]=a.p[i];
	
	return *this;
		
}

template <class Type>
Matrix<Type> Matrix<Type>::col(int icol) const {

	Matrix<Type> res(nf,1);
	Type *pi;
	int i,N;
	
	if(icol<0) icol+=nc;
	if(icol<0||icol>=nc) {
		fprintf(stderr,"ERROR: (matrix.col) Index exceeds matrix dimensions\n");
		exit(1);
	}
	
	pi=p+icol*nf;
	N=nf;
	for(i=0;i<N;i++)
		res.p[i]=pi[i];
	
	return res;
}

template <class Type>
Matrix<Type> &Matrix<Type>::setcol(int icol,const Matrix<Type> &a) {

	Type *pi;
	int i,N;
	
	if(icol<0) icol+=nc;
	if(icol<0||icol>=nc) {
		fprintf(stderr,"ERROR: (matrix.setcol) Index exceeds matrix dimensions\n");
		exit(1);
	}
	if(a.nc>1||a.nf!=nf) {
		fprintf(stderr,"ERROR: (matrix.setcol) Dimensions must agree\n");
		exit(1);
	}	
	
	pi=p+icol*nf;
	N=nf;
	for(i=0;i<N;i++)
		pi[i]=a.p[i];
		
	return *this;
		
}

template <class Type>
Matrix<Type> Matrix<Type>::block(int ifil1,int ifil2,int icol1,int icol2) const {

	if(ifil1<0) ifil1+=nf;
	if(ifil2<0) ifil2+=nf;
	if(icol1<0) icol1+=nc;
	if(icol2<0) icol2+=nc;
	
	if(ifil1<0||ifil1>=nf||ifil2<0||ifil2>=nf||icol1<0||icol1>=nc||icol2<0||icol2>=nc) {
		fprintf(stderr,"ERROR: (matrix.block) Index exceeds matrix dimensions\n");
		exit(1);
	}	

	Matrix<Type> res(ifil2-ifil1+1,icol2-icol1+1);
	Type *pi,*pres;
	int i,j,N1,N2;
	
	pi=p+ifil1+icol1*nf;pres=res.p;
	N1=res.nc;N2=res.nf;
	for(i=0;i<N1;i++) {
		for(j=0;j<N2;j++)
			pres[j]=pi[j];
		pres+=res.nf;pi+=nf;
	}
	
	return res;
}

template <class Type>
Matrix<Type> Matrix<Type>::block_step(int ifil1,int ifil2,int dfil,int icol1,int icol2,int dcol) const {

	if(ifil1<0) ifil1+=nf;
	if(ifil2<0) ifil2+=nf;
	if(icol1<0) icol1+=nc;
	if(icol2<0) icol2+=nc;
	
	if(ifil1<0||ifil1>=nf||ifil2<0||ifil2>=nf||icol1<0||icol1>=nc||icol2<0||icol2>=nc) {
		fprintf(stderr,"ERROR: (matrix.block_step) Index exceeds matrix dimensions\n");
		exit(1);
	}	

	Matrix<Type> res(floor((ifil2-ifil1)/dfil)+1,floor((icol2-icol1)/dcol)+1);
	Type *pi,*pres;
	int i,j,N1,N2;
	
	pi=p+ifil1+icol1*nf;pres=res.p;
	N1=res.nc;N2=res.nf;
	for(i=0;i<N1;i++) {
		for(j=0;j<N2;j++)
			pres[j]=pi[j*dfil];
		pres+=res.nf;pi+=dcol*nf;
	}
	
	return res;
}

template <class Type>
Matrix<Type> &Matrix<Type>::setblock(int ifil1,int ifil2,int icol1,int icol2,const Matrix<Type> &a) {

	if(ifil1<0) ifil1+=nf;
	if(ifil2<0) ifil2+=nf;
	if(icol1<0) icol1+=nc;
	if(icol2<0) icol2+=nc;
	
	if(ifil1<0||ifil1>=nf||ifil2<0||ifil2>=nf||icol1<0||icol1>=nc||icol2<0||icol2>=nc) {
		fprintf(stderr,"ERROR: (matrix.setblock) Index exceeds matrix dimensions\n");
		exit(1);
	}	
	if(a.nf!=ifil2-ifil1+1||a.nc!=icol2-icol1+1) {
		fprintf(stderr,"ERROR: (matrix.setblock) Dimensions must agree\n");
		exit(1);
	}

	Type *pi,*pa;
	int i,j,N1,N2;
	
	pi=p+ifil1+icol1*nf;pa=a.p;
	N1=a.nc;N2=a.nf;
	for(i=0;i<N1;i++) {
		for(j=0;j<N2;j++)
			pi[j]=pa[j];
		pa+=a.nf;pi+=nf;
	}
	
	return *this;
	
}

template <class Type>
Matrix<Type> &Matrix<Type>::setblock_step(int ifil1,int ifil2,int dfil,int icol1,int icol2,int dcol,const Matrix<Type> &a) {

	if(ifil1<0) ifil1+=nf;
	if(ifil2<0) ifil2+=nf;
	if(icol1<0) icol1+=nc;
	if(icol2<0) icol2+=nc;
	
	if(ifil1<0||ifil1>=nf||ifil2<0||ifil2>=nf||icol1<0||icol1>=nc||icol2<0||icol2>=nc) {
		fprintf(stderr,"ERROR: (matrix.setblock_step) Index exceeds matrix dimensions\n");
		exit(1);
	}
	if(a.nf!=floor((ifil2-ifil1)/dfil)+1||a.nc!=floor((icol2-icol1)/dcol)+1) {
		fprintf(stderr,"ERROR: (matrix.setblock_step) Dimensions must agree\n");
		exit(1);
	}

	Type *pi,*pa;
	int i,j,N1,N2;
	
	pi=p+ifil1+icol1*nf;pa=a.p;
	N1=a.nc;N2=a.nf;
	for(i=0;i<N1;i++) {
		for(j=0;j<N2;j++)
			pi[j*dfil]=pa[j];
		pa+=a.nf;pi+=dcol*nf;
	}
	
	return *this;
	
}

template <class Type>
Matrix<Type> Matrix<Type>::transpose() const {

	Matrix<Type> a(nc,nf);
	int i,j,N1,N2;
	Type *pi,*pa;
	
	N1=nf;N2=nc;
	pi=p;pa=a.p;
	
	for(i=0;i<N1;i++) {
		for(j=0;j<N2;j++)
			pa[j]=pi[j*nf];
		pa+=a.nf;pi++;
	}
	
	return a;
}

template <class Type>
Matrix<Type> Matrix<Type>::fliplr() const {

	Matrix<Type> a(nf,nc);
	int i,j,N1,N2;
	Type *pi,*pa;
	
	N1=nc;N2=nf;
	pi=p;pa=a.p;
	pi+=nf*(nc-1);
	
	for(i=0;i<N1;i++) {
		for(j=0;j<N2;j++)
			pa[j]=pi[j];
		pa+=a.nf;pi-=nf;
	}
	
	return a;
}

template <class Type>
Matrix<Type> Matrix<Type>::flipud() const {

	Matrix<Type> a(nf,nc);
	int i,j,N1,N2;
	Type *pi,*pa;
	
	N1=nc;N2=nf;
	pi=p;pa=a.p;

	for(i=0;i<N1;i++) {
		for(j=0;j<N2;j++)
			pa[j]=pi[nf-1-j];
		pa+=a.nf;pi+=nf;
	}
	
	return a;
}

matrix ones(int nfil,int ncol) {
	
	matrix a(nfil,ncol);
	int i,N;
	
	N=nfil*ncol;
	for(i=0;i<N;i++)
		a.p[i]=1;
		
	return a;
}

matrix zeros(int nfil,int ncol) {
	
	matrix a(nfil,ncol);
	int i,N;
	
	N=nfil*ncol;
	for(i=0;i<N;i++)
		a.p[i]=0;
		
	return a;
}

matrix random_matrix(int nfil,int ncol) {

	matrix a(nfil,ncol);
	int i,N;
	
	N=nfil*ncol;
	for(i=0;i<N;i++)
		a.p[i]=(mat_type) rand()/RAND_MAX;
		
	return a;

}


matrix eye(int n) {

	matrix a;
	int i,d,N;
	
	a.zero(n,n);
	d=n+1;
	N=n*n;
	for(i=0;i<N;i+=d)
		a.p[i]=1;
	return a;
}


matrix vector_t(mat_type x0,mat_type x1,int n) {

	mat_type x,dx;
	int i;
	
	dx=(x1-x0)/(n-1);
	matrix v(n,1);
	
	for(i=0;i<n;i++)
		v.p[i]=x0+i*dx;
	v.p[n-1]=x1;
	
	return v;
}

matrix vector(mat_type x0,mat_type x1,int n) {

	matrix v;
	
	v=vector_t(x0,x1,n);
	v.nc=v.nf;
	v.nf=1;
	
	return v;
}

template <class Type>
Type max(const Matrix<Type> &a) {

	Type x;
	int i,N;
	
	x=*(a.p);
	N=a.nf*a.nc;
	for(i=1;i<N;i++) 
		if (a.p[i]>x) x=a.p[i];

	return x;
}

template <class Type>
Type min(const Matrix<Type> &a) {

	Type x;
	int i,N;
	
	x=*(a.p);
	N=a.nf*a.nc;
	for(i=1;i<N;i++) 
		if (a.p[i]<x) x=a.p[i];

	return x;
}

template <class Type>
Type sum(const Matrix<Type> &a) {

	Type s=0;
	int i,N;
	
	N=a.nf*a.nc;
	for(i=0;i<N;i++) 
		s+=a.p[i];

	return s;
}

template <class Type>
Type mean(const Matrix<Type> &a) {

	return sum(a)/a.nf/a.nc;
}

template <class Type>
Matrix<Type> max(const Matrix<Type> &a,const Matrix<Type> &b) {

	if( (b.nf!=a.nf) || (b.nc!=a.nc) ) {
		fprintf(stderr,"ERROR: (matrix.max) Dimensions must agree\n");
		exit(1);
	}

    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	res.p[i]=a.p[i]>b.p[i]?a.p[i]:b.p[i];
    return res;
}

template <class Type>
Matrix<Type> max(const Matrix<Type> &a,mat_type n) {

    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	res.p[i]=a.p[i]>n?a.p[i]:n;
    return res;
}

template <class Type>
Matrix<Type> max(mat_type n,const Matrix<Type> &a) {
	return max(a,n);
}

template <class Type>
Matrix<Type> min(const Matrix<Type> &a,const Matrix<Type> &b) {

	if( (b.nf!=a.nf) || (b.nc!=a.nc) ) {
		fprintf(stderr,"ERROR: (matrix.min) Dimensions must agree\n");
		exit(1);
	}

    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	res.p[i]=a.p[i]<b.p[i]?a.p[i]:b.p[i];
    return res;
}

template <class Type>
Matrix<Type> min(const Matrix<Type> &a,mat_type n) {

    Matrix<Type> res(a.nf,a.nc);
    int i,N;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	res.p[i]=a.p[i]<n?a.p[i]:n;
    return res;
}

template <class Type>
Matrix<Type> min(mat_type n,const Matrix<Type> &a) {
	return min(a,n);
}

template <class Type>
int exist(const Matrix<Type> &a) {

	int i,N;
	N=a.nf*a.nc;
	for(i=0;i<N;i++)
		if(a.p[i]) return 1;
	return 0;
}

template <class Type>
int isequal(const Matrix<Type> &a,const Matrix<Type> &b) {

	int i,N,res;
	if(a.nf!=b.nf) return 0;
	if(a.nc!=b.nc) return 0;
	N=a.nf*a.nc;
	res=1;
	for(i=0;i<N;i++)
		if(a.p[i]!=b.p[i]) res=0;
	return res;
}

// Explicit instantiations

template class Matrix<double>;
template class Matrix<float>;

template Matrix<double> operator-(mat_type,const Matrix<double> &);
template Matrix<double> operator-(const Matrix<double> &);
template Matrix<double> operator/(mat_type,const Matrix<double> &);
template double max(const Matrix<double> &);
template double min(const Matrix<double> &);
template double sum(const Matrix<double> &);
template double mean(const Matrix<double> &);
template Matrix<double> max(const Matrix<double> &,const Matrix<double> &);
template Matrix<double> max(const Matrix<double> &,mat_type);
template Matrix<double> max(mat_type,const Matrix<double> &);
template Matrix<double> min(const Matrix<double> &,const Matrix<double> &);
template Matrix<double> min(const Matrix<double> &,mat_type);
template Matrix<double> min(mat_type,const Matrix<double> &);
template int exist(const Matrix<double> &);
template int isequal(const Matrix<double> &,const Matrix<double> &);

template Matrix<float> operator-(mat_type,const Matrix<float> &);
template Matrix<float> operator-(const Matrix<float> &);
template Matrix<float> operator/(mat_type,const Matrix<float> &);
template float max(const Matrix<float> &);
template float min(const Matrix<float> &);
template float sum(const Matrix<float> &);
template float mean(const Matrix<float> &);
template Matrix<float> max(const Matrix<float> &,const Matrix<float> &);
template Matrix<float> max(const Matrix<float> &,mat_type);
template Matrix<float> max(mat_type,const Matrix<float> &);
template Matrix<float> min(const Matrix<float> &,const Matrix<float> &);
template Matrix<float> min(const Matrix<float> &,mat_type);
template Matrix<float> min(mat_type,const Matrix<float> &);
template int exist(const Matrix<float> &);
template int isequal(const Matrix<float> &,const Matrix<float> &);


