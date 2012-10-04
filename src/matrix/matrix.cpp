#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<string.h>
#include"matrix.h"

matrix::matrix(int nfil,int ncol) {
	
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
	p=new double[tam];
}

matrix::~matrix() {

	delete [] p;
}

matrix::matrix(const matrix &a) {

    unsigned tam;

    nc=a.nc;
    nf=a.nf;
    tam=unsigned(nf)*unsigned(nc);
	p=new double[tam];

	memcpy(p,a.p,nc*nf*sizeof(double));

}

int matrix::nrows() const {

	return nf;
}

int matrix::ncols() const {

	return nc;
}

double *matrix::data() const {

	return p;
}

matrix &matrix::dim(int nfil,int ncol) {
    
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
    
    if(nfil*ncol!=nf*nc) {
    	delete [] p;
    	nf=nfil;
    	nc=ncol;
    	tam=unsigned(nf)*unsigned(nc);
		p=new double[tam];
    } else {
    	nf=nfil;
    	nc=ncol;
    }
    
    return *this;
}

matrix &matrix::redim(int nfil,int ncol) {
    
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
    	fprintf(stderr,"ERROR: matrix.redim: Number of elements doesn't match\n");
    	exit(1);
    }    
    nf=nfil;nc=ncol;
    
    return *this;
}

matrix &matrix::operator=(const matrix &a) {
        
	if (nf!=a.nf||nc!=a.nc)
    	dim(a.nf,a.nc);
    
	memcpy(p,a.p,nc*nf*sizeof(double));
    
    return *this;
}

double &matrix::operator()(int ifil,int icol) const {

	if(ifil<0) ifil+=nf;
	if(icol<0) icol+=nc;
	if(ifil>=nf||ifil<0||icol>=nc||icol<0) {
		fprintf(stderr,"ERROR: Index exceed matrix dimensions\n");
		exit(1);
	}
	return *(p+icol*nf+ifil);
	
}

double &matrix::operator()(int ielem) const {
	
	if(ielem<0) ielem+=nf*nc;
	if(ielem>=nf*nc||ielem<0) {
		fprintf(stderr,"ERROR: Index exceed matrix dimensions\n");
		exit(1);
	}
	return *(p+ielem);
	
}

int matrix::read(int nfil,int ncol,FILE *fp,char mode) {
    
    int i,tam;
    double temp;

    dim(nfil,ncol);
    tam=nc*nf;           

    if (mode=='b') {
        if((int)fread(p,sizeof(double),tam,fp)!=tam) return 1;
        }
    else
        for(i=0;i<tam;i++) {
            if(fscanf(fp,"%le",&temp)==EOF) return 1;
            *(p+i)=(double)temp;
            }
    return 0;
}

int matrix::write(FILE *fp,char mode) const {
   
    int tam,i,j;

    tam=nf*nc;
    if (mode=='b') {
        if((int)fwrite(p,sizeof(double),tam,fp)!=tam) return 1;
        }
    else {
        j=0;
        for(i=0;i<tam;i++) {
            if(fprintf(fp,"%16.16e ",(double)*(p+i))==EOF) return 1;
            j++;
            if (j==nf) {
                j=0;
                fprintf(fp,"\n");
            }
        }
    }
    return 0;
}

matrix matrix::operator+(const matrix &a) const {

    matrix res;
    double *pi,*pa,*pres;
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
	pi=p;pa=a.p;pres=res.p;
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    *(pres+i)=*(pi+i)+(*(pa+i));
       	}
	else {
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

matrix operator+(double n,const matrix &a) {
 	matrix res(a.nf,a.nc);
    double *pa,*pres;
    int i,N;
    
    pa=a.p;
    pres=res.p;
    N=a.nc*a.nf;
   
    for(i=0;i<N;i++)
    	*(pres+i)=n+(*(pa+i));
    return res;
}

matrix matrix::operator+(double n) const {

    matrix res(nf,nc);
    double *pa,*pres;
    int i,N;
    
    pa=p;
    pres=res.p;
    N=nc*nf;
    for(i=0;i<N;i++)
    	*(pres+i)=*(pa+i)+n;
    	
    return res;
}

matrix operator+(const matrix &a) {

 	matrix res(a);
    
    return res;
}

matrix matrix::operator-(const matrix &a) const {

    matrix res;
    double *pi,*pa,*pres;
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
	pi=p;pa=a.p;pres=res.p;
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    *(pres+i)=*(pi+i)-(*(pa+i));
       	}
	else {
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

matrix operator-(double n,const matrix &a) {
 	matrix res(a.nf,a.nc);
    double *pa,*pres;
    int i,N;
    
    pa=a.p;
    pres=res.p;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	*(pres+i)=n-(*(pa+i));
    return res;
}

matrix matrix::operator-(double n) const {

    matrix res(nf,nc);
    double *pa,*pres;
    int i,N;
    pa=p;
    pres=res.p;
    N=nc*nf;
    for(i=0;i<N;i++)
    	*(pres+i)=*(pa+i)-n;
    	
    return res;
}


matrix operator-(const matrix &a) {
    
    matrix res(a.nf,a.nc);
    double *pa,*pres;
    int i,N;
    
    pa=a.p;
    pres=res.p;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	*(pres+i)=-*(pa+i);
    	
    return res;
}


matrix matrix::operator*(const matrix &a) const {

    matrix res;
    double *pi,*pa,*pres;
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
	pi=p;pa=a.p;pres=res.p;
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nf*nc;
		for(i=0;i<N;i++)
       	    *(pres+i)=*(pi+i)*(*(pa+i));}
	else {
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

matrix operator*(double n,const matrix &a) {

 	matrix res(a.nf,a.nc);
    double *pa,*pres;
    int i,N;
    
    pa=a.p;
    pres=res.p;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	*(pres+i)=*(pa+i)*n;
    
    return res;
}

matrix matrix::operator*(double n) const {
    
    matrix res(nf,nc);
    double *pa,*pres;
    int i,N;
    
    pa=p;
    pres=res.p;
    N=nc*nf;
    for(i=0;i<N;i++)
    	*(pres+i)=*(pa+i)*n;
    return res;
}

matrix matrix::operator/(const matrix &a) const {

    matrix res;
    double *pi,*pa,*pres;
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
	pi=p;pa=a.p;pres=res.p;
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    *(pres+i)=*(pi+i)/(*(pa+i));
       	}
	else {
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

matrix operator/(double n,const matrix &a) {
 	matrix res(a.nf,a.nc);
    double *pa,*pres;
    int i,N;
    
    pa=a.p;
    pres=res.p;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	*(pres+i)=n/(*(pa+i));	
    return res;
}
matrix matrix::operator/(double n) const {
    matrix res(nf,nc);
    double *pa,*pres;
    int i,N;
    pa=p;
    pres=res.p;
    N=nc*nf;
    for(i=0;i<N;i++)
    	*(pres+i)=*(pa+i)/n;
    return res;
}

matrix matrix::operator==(const matrix &a) const {

    matrix res;
    double *pi,*pa,*pres;
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
	pi=p;pa=a.p;pres=res.p;
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    *(pres+i)=(*(pi+i)==(*(pa+i)));
       	}
	else {
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

matrix operator==(double n,const matrix &a) {
 	matrix res(a.nf,a.nc);
    double *pa,*pres;
    int i,N;
    
    pa=a.p;
    pres=res.p;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	*(pres+i)=(n==(*(pa+i)));	
    return res;
}

matrix matrix::operator==(double n) const {
    matrix res(nf,nc);
    double *pa,*pres;
    int i,N;
    pa=p;
    pres=res.p;
    N=nc*nf;
    for(i=0;i<N;i++)
    	*(pres+i)=(*(pa+i)==n);
    return res;
}

matrix matrix::operator!=(const matrix &a) const {

    matrix res;
    double *pi,*pa,*pres;
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
	pi=p;pa=a.p;pres=res.p;
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    *(pres+i)=(*(pi+i)!=(*(pa+i)));
       	}
	else {
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

matrix operator!=(double n,const matrix &a) {
 	matrix res(a.nf,a.nc);
    double *pa,*pres;
    int i,N;
    
    pa=a.p;
    pres=res.p;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	*(pres+i)=(n!=(*(pa+i)));	
    return res;
}

matrix matrix::operator!=(double n) const {
    matrix res(nf,nc);
    double *pa,*pres;
    int i,N;
    pa=p;
    pres=res.p;
    N=nc*nf;
    for(i=0;i<N;i++)
    	*(pres+i)=(*(pa+i)!=n);
    return res;
}

matrix matrix::operator>(const matrix &a) const{

    matrix res;
    double *pi,*pa,*pres;
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
	pi=p;pa=a.p;pres=res.p;
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    *(pres+i)=(*(pi+i)>(*(pa+i)));
       	}
	else {
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

matrix operator>(double n,const matrix &a) {
 	matrix res(a.nf,a.nc);
    double *pa,*pres;
    int i,N;
    
    pa=a.p;
    pres=res.p;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	*(pres+i)=(n>(*(pa+i)));	
    return res;
}

matrix matrix::operator>(double n) const {
    matrix res(nf,nc);
    double *pa,*pres;
    int i,N;
    pa=p;
    pres=res.p;
    N=nc*nf;
    for(i=0;i<N;i++)
    	*(pres+i)=(*(pa+i)>n);
    return res;
}

matrix matrix::operator<(const matrix &a) const {

    matrix res;
    double *pi,*pa,*pres;
	int i,j,resnf,resnc,N;
	
	if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
		fprintf(stderr,"ERROR: (matrix.<) Dimensions must agree\n");
		exit(1);
	}
	
	if(nf>a.nf) resnf=nf;
		else resnf=a.nf;
	if(nc>a.nc) resnc=nc;
		else resnc=a.nc;
	
	res.dim(resnf,resnc);
	pi=p;pa=a.p;pres=res.p;
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    *(pres+i)=(*(pi+i)<(*(pa+i)));
       	}
	else {
		for(i=0;i<resnc;i++) {
    	   	for(j=0;j<resnf;j++) {
    	   	    *(pres++)=(*(pi)<(*(pa)));
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

matrix operator<(double n,const matrix &a) {
 	matrix res(a.nf,a.nc);
    double *pa,*pres;
    int i,N;
    
    pa=a.p;
    pres=res.p;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	*(pres+i)=(n<(*(pa+i)));	
    return res;
}

matrix matrix::operator<(double n) const {
    matrix res(nf,nc);
    double *pa,*pres;
    int i,N;
    pa=p;
    pres=res.p;
    N=nc*nf;
    for(i=0;i<N;i++)
    	*(pres+i)=(*(pa+i)<n);
    return res;
}

matrix matrix::operator>=(const matrix &a) const {

    matrix res;
    double *pi,*pa,*pres;
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
	pi=p;pa=a.p;pres=res.p;
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    *(pres+i)=(*(pi+i)>=(*(pa+i)));
       	}
	else {
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

matrix operator>=(double n,const matrix &a) {
 	matrix res(a.nf,a.nc);
    double *pa,*pres;
    int i,N;
    
    pa=a.p;
    pres=res.p;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	*(pres+i)=(n>=(*(pa+i)));	
    return res;
}

matrix matrix::operator>=(double n) const {
    matrix res(nf,nc);
    double *pa,*pres;
    int i,N;
    pa=p;
    pres=res.p;
    N=nc*nf;
    for(i=0;i<N;i++)
    	*(pres+i)=(*(pa+i)>=n);
    return res;
}

matrix matrix::operator<=(const matrix &a) const {

    matrix res;
    double *pi,*pa,*pres;
	int i,j,resnf,resnc,N;
	
	if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
		fprintf(stderr,"ERROR: (matrix.<=) Dimensions must agree\n");
		exit(1);
	}
	
	if(nf>a.nf) resnf=nf;
		else resnf=a.nf;
	if(nc>a.nc) resnc=nc;
		else resnc=a.nc;
	
	res.dim(resnf,resnc);
	pi=p;pa=a.p;pres=res.p;
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    *(pres+i)=(*(pi+i)<=(*(pa+i)));
       	}
	else {
		for(i=0;i<resnc;i++) {
    	   	for(j=0;j<resnf;j++) {
    	   	    *(pres++)=(*(pi)<=(*(pa)));
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

matrix operator<=(double n,const matrix &a) {
 	matrix res(a.nf,a.nc);
    double *pa,*pres;
    int i,N;
    
    pa=a.p;
    pres=res.p;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	*(pres+i)=(n<=(*(pa+i)));	
    return res;
}

matrix matrix::operator<=(double n) const {
    matrix res(nf,nc);
    double *pa,*pres;
    int i,N;
    pa=p;
    pres=res.p;
    N=nc*nf;
    for(i=0;i<N;i++)
    	*(pres+i)=(*(pa+i)<=n);
    return res;
}

matrix matrix::operator||(const matrix &a) const {

    matrix res;
    double *pi,*pa,*pres;
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
	pi=p;pa=a.p;pres=res.p;
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    *(pres+i)=(*(pi+i)||(*(pa+i)));
       	}
	else {
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

matrix operator||(double n,const matrix &a) {
 	matrix res(a.nf,a.nc);
    double *pa,*pres;
    int i,N;
    
    pa=a.p;
    pres=res.p;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	*(pres+i)=(n||(*(pa+i)));	
    return res;
}

matrix matrix::operator||(double n) const {
    matrix res(nf,nc);
    double *pa,*pres;
    int i,N;
    pa=p;
    pres=res.p;
    N=nc*nf;
    for(i=0;i<N;i++)
    	*(pres+i)=(*(pa+i)||n);
    return res;
}

matrix matrix::operator&&(const matrix &a) const {

    matrix res;
    double *pi,*pa,*pres;
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
	pi=p;pa=a.p;pres=res.p;
	if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
		N=nc*nf;
		for(i=0;i<N;i++)
       	    *(pres+i)=(*(pi+i)&&(*(pa+i)));
       	}
	else {
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

matrix operator&&(double n,const matrix &a) {
 	matrix res(a.nf,a.nc);
    double *pa,*pres;
    int i,N;
    
    pa=a.p;
    pres=res.p;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	*(pres+i)=(n&&(*(pa+i)));	
    return res;
}

matrix matrix::operator&&(double n) const {
    matrix res(nf,nc);
    double *pa,*pres;
    int i,N;
    pa=p;
    pres=res.p;
    N=nc*nf;
    for(i=0;i<N;i++)
    	*(pres+i)=(*(pa+i)&&n);
    return res;
}


matrix &matrix::operator+=(const matrix &a) {
    
    int i,N;
    double *pa,*pi;

	if( (nf!=a.nf) || (nc!=a.nc) ) {
		*this=(*this)+a;
		return *this;
	}

    pa=a.p;pi=p;
    N=nc*nf;
    for(i=0;i<N;i++) 
    	*(pi+i)+=*(pa+i);
    return *this;
}

matrix &matrix::operator-=(const matrix &a) {
    
    int i,N;
    double *pa,*pi;

	if( (nf!=a.nf) || (nc!=a.nc) ) {
		*this=(*this)-a;
		return *this;
	}

    pa=a.p;pi=p;
    N=nc*nf;
    for(i=0;i<N;i++) 
    	*(pi+i)-=*(pa+i);
    return *this;
}

matrix &matrix::operator*=(const matrix &a) {
    
    int i,N;
    double *pa,*pi;

	if( (nf!=a.nf) || (nc!=a.nc) ) {
		*this=(*this)*a;
		return *this;
	}

    pa=a.p;pi=p;
    N=nc*nf;
    for(i=0;i<N;i++) 
    	*(pi+i)*=*(pa+i);
    return *this;
}

matrix &matrix::operator/=(const matrix &a) {
    
    int i,N;
    double *pa,*pi;

	if( (nf!=a.nf) || (nc!=a.nc) ) {
		*this=(*this)/a;
		return *this;
	}

    pa=a.p;pi=p;
    N=nc*nf;
    for(i=0;i<N;i++) 
    	*(pi+i)/=*(pa+i);
    return *this;
}

matrix &matrix::operator+=(double n) {
    
    int i,N;
    double *pi;

    pi=p;
    N=nc*nf;
    for(i=0;i<N;i++) 
    	*(pi+i)+=n;
    return *this;
}

matrix &matrix::operator-=(double n) {
    
    int i,N;
    double *pi;

    pi=p;
    N=nc*nf;
    for(i=0;i<N;i++) 
    	*(pi+i)-=n;
    return *this;
}

matrix &matrix::operator*=(double n) {
    
    int i,N;
    double *pi;

    pi=p;
    N=nc*nf;
    for(i=0;i<N;i++) 
    	*(pi+i)*=n;
    return *this;
}

matrix &matrix::operator/=(double n) {
    
    int i,N;
    double *pi;

    pi=p;
    N=nc*nf;
    for(i=0;i<N;i++) 
    	*(pi+i)/=n;
    return *this;
}

matrix matrix::row(int ifil) const {

	matrix res(1,nc);
	double *pi,*pres;
	int i;
	
	if(ifil<0) ifil+=nf;
	if(ifil<0||ifil>=nf) {
		fprintf(stderr,"ERROR: (matrix.row) Index exceed matrix dimensions\n");
		exit(1);
	}
			
	pi=p+ifil;pres=res.p;
	for(i=0;i<nc;i++)
		*(pres+i)=*(pi+i*nf);
	
	return res;
}

matrix &matrix::setrow(int ifil,const matrix &a) {

	double *pi,*pa;
	int i;
	
	if(ifil<0) ifil+=nf;
	if(ifil<0||ifil>=nf) {
		fprintf(stderr,"ERROR: (matrix.setrow) Index exceed matrix dimensions\n");
		exit(1);
	}
	if(a.nf>1||a.nc!=nc) {
		fprintf(stderr,"ERROR: (matrix.setrow) Dimensions must agree\n");
		exit(1);
	}	
	
	pi=p+ifil;pa=a.p;
	for(i=0;i<nc;i++)
		*(pi+i*nf)=*(pa+i);
	
	return *this;
		
}


matrix matrix::col(int icol) const {

	matrix res(nf,1);
	double *pi,*pres;
	int i,N;
	
	if(icol<0) icol+=nc;
	if(icol<0||icol>=nc) {
		fprintf(stderr,"ERROR: (matrix.col) Index exceed matrix dimensions\n");
		exit(1);
	}
	
	pi=p+icol*nf;pres=res.p;
	N=nf;
	for(i=0;i<N;i++)
		*(pres+i)=*(pi+i);
	
	return res;
}

matrix &matrix::setcol(int icol,const matrix &a) {

	double *pi,*pa;
	int i,N;
	
	if(icol<0) icol+=nc;
	if(icol<0||icol>=nc) {
		fprintf(stderr,"ERROR: (matrix.setcol) Index exceed matrix dimensions\n");
		exit(1);
	}
	if(a.nc>1||a.nf!=nf) {
		fprintf(stderr,"ERROR: (matrix.setcol) Dimensions must agree\n");
		exit(1);
	}	
	
	pi=p+icol*nf;pa=a.p;
	N=nf;
	for(i=0;i<N;i++)
		*(pi+i)=*(pa+i);
		
	return *this;
		
}

matrix matrix::block(int ifil1,int ifil2,int icol1,int icol2) const {

	if(ifil1<0) ifil1+=nf;
	if(ifil2<0) ifil2+=nf;
	if(icol1<0) icol1+=nc;
	if(icol2<0) icol2+=nc;
	
	if(ifil1<0||ifil1>=nf||ifil2<0||ifil2>=nf||icol1<0||icol1>=nc||icol2<0||icol2>=nc) {
		fprintf(stderr,"ERROR: (matrix.block) Index exceed matrix dimensions\n");
		exit(1);
	}	

	matrix res(ifil2-ifil1+1,icol2-icol1+1);
	double *pi,*pres;
	int i,j,N1,N2;
	
	pi=p+ifil1+icol1*nf;pres=res.p;
	N1=res.nc;N2=res.nf;
	for(i=0;i<N1;i++) {
		for(j=0;j<N2;j++)
			*(pres+j)=*(pi+j);
		pres+=res.nf;pi+=nf;
	}
	
	return res;
}

matrix matrix::block_step(int ifil1,int ifil2,int dfil,int icol1,int icol2,int dcol) const {

	if(ifil1<0) ifil1+=nf;
	if(ifil2<0) ifil2+=nf;
	if(icol1<0) icol1+=nc;
	if(icol2<0) icol2+=nc;
	
	if(ifil1<0||ifil1>=nf||ifil2<0||ifil2>=nf||icol1<0||icol1>=nc||icol2<0||icol2>=nc) {
		fprintf(stderr,"ERROR: (matrix.block_step) Index exceed matrix dimensions\n");
		exit(1);
	}	

	matrix res(floor((ifil2-ifil1)/dfil)+1,floor((icol2-icol1)/dcol)+1);
	double *pi,*pres;
	int i,j,N1,N2;
	
	pi=p+ifil1+icol1*nf;pres=res.p;
	N1=res.nc;N2=res.nf;
	for(i=0;i<N1;i++) {
		for(j=0;j<N2;j++)
			*(pres+j)=*(pi+j*dfil);
		pres+=res.nf;pi+=dcol*nf;
	}
	
	return res;
}

matrix &matrix::setblock(int ifil1,int ifil2,int icol1,int icol2,const matrix &a) {

	if(ifil1<0) ifil1+=nf;
	if(ifil2<0) ifil2+=nf;
	if(icol1<0) icol1+=nc;
	if(icol2<0) icol2+=nc;
	
	if(ifil1<0||ifil1>=nf||ifil2<0||ifil2>=nf||icol1<0||icol1>=nc||icol2<0||icol2>=nc) {
		fprintf(stderr,"ERROR: (matrix.setblock) Index exceed matrix dimensions\n");
		exit(1);
	}	
	if(a.nf!=ifil2-ifil1+1||a.nc!=icol2-icol1+1) {
		fprintf(stderr,"ERROR: (matrix.setblock) Dimensions must agree\n");
		exit(1);
	}

	double *pi,*pa;
	int i,j,N1,N2;
	
	pi=p+ifil1+icol1*nf;pa=a.p;
	N1=a.nc;N2=a.nf;
	for(i=0;i<N1;i++) {
		for(j=0;j<N2;j++)
			*(pi+j)=*(pa+j);
		pa+=a.nf;pi+=nf;
	}
	
	return *this;
	
}

matrix &matrix::setblock_step(int ifil1,int ifil2,int dfil,int icol1,int icol2,int dcol,const matrix &a) {

	if(ifil1<0) ifil1+=nf;
	if(ifil2<0) ifil2+=nf;
	if(icol1<0) icol1+=nc;
	if(icol2<0) icol2+=nc;
	
	if(ifil1<0||ifil1>=nf||ifil2<0||ifil2>=nf||icol1<0||icol1>=nc||icol2<0||icol2>=nc) {
		fprintf(stderr,"ERROR: (matrix.setblock_step) Index exceed matrix dimensions\n");
		exit(1);
	}
	if(a.nf!=floor((ifil2-ifil1)/dfil)+1||a.nc!=floor((icol2-icol1)/dcol)+1) {
		fprintf(stderr,"ERROR: (matrix.setblock_step) Dimensions must agree\n");
		exit(1);
	}

	double *pi,*pa;
	int i,j,N1,N2;
	
	pi=p+ifil1+icol1*nf;pa=a.p;
	N1=a.nc;N2=a.nf;
	for(i=0;i<N1;i++) {
		for(j=0;j<N2;j++)
			*(pi+j*dfil)=*(pa+j);
		pa+=a.nf;pi+=dcol*nf;
	}
	
	return *this;
	
}


matrix ones(int nfil,int ncol) {
	
	matrix a(nfil,ncol);
	double *pa;
	int i,N;
	
	N=nfil*ncol;
	pa=a.p;
	for(i=0;i<N;i++)
		*(pa+i)=1;
		
	return a;
}

matrix zeros(int nfil,int ncol) {
	
	matrix a(nfil,ncol);
	double *pa;
	int i,N;
	
	N=nfil*ncol;
	pa=a.p;
	for(i=0;i<N;i++)
		*(pa+i)=0;
		
	return a;
}

matrix random_matrix(int nfil,int ncol) {

	matrix a(nfil,ncol);
	double *pa;
	int i,N;
	
	N=nfil*ncol;
	pa=a.p;
	for(i=0;i<N;i++)
		*(pa+i)=(double) rand()/RAND_MAX;
		
	return a;

}


matrix eye(int n) {

	matrix a;
	double *pa;
	int i,d,N;
	
	a=zeros(n,n);
	pa=a.p;
	d=n+1;
	N=n*n;
	for(i=0;i<N;i+=d)
		*(pa+i)=1;
	return a;
}

matrix vector_t(double x0,double x1,int n) {

	double x,*p,dx;
	int i;
	
	dx=(x1-x0)/(n-1);
	matrix v(n,1);
	
	p=v.p;
	for(i=0;i<n;i++)
		*(p+i)=x0+i*dx;
	*(p+n-1)=x1;
	
	return v;
}

matrix vector(double x0,double x1,int n) {

	matrix v;
	
	v=vector_t(x0,x1,n);
	v.nc=v.nf;
	v.nf=1;
	
	return v;
}


matrix matrix::transpose() const {

	matrix a(nc,nf);
	int i,j,N1,N2;
	double *pi,*pa;
	
	N1=nf;N2=nc;
	pi=p;pa=a.p;
	
	for(i=0;i<N1;i++) {
		for(j=0;j<N2;j++)
			*(pa+j)=*(pi+j*nf);
		pa+=a.nf;pi++;
	}
	
	return a;
}

matrix matrix::fliplr() const {

	matrix a(nf,nc);
	int i,j,N1,N2;
	double *pi,*pa;
	
	N1=nc;N2=nf;
	pi=p;pa=a.p;
	pi+=nf*(nc-1);
	
	for(i=0;i<N1;i++) {
		for(j=0;j<N2;j++)
			*(pa+j)=*(pi+j);
		pa+=a.nf;pi-=nf;
	}
	
	return a;
}

matrix matrix::flipud() const {

	matrix a(nf,nc);
	int i,j,N1,N2;
	double *pi,*pa;
	
	N1=nc;N2=nf;
	pi=p;pa=a.p;

	for(i=0;i<N1;i++) {
		for(j=0;j<N2;j++)
			*(pa+j)=*(pi+nf-1-j);
		pa+=a.nf;pi+=nf;
	}
	
	return a;
}

double max(const matrix &a) {

	double x,*p;
	int i,N;
	
	x=*(a.p);
	N=a.nf*a.nc;
	p=a.p;
	for(i=0;i<N;i++) {
		if (*p>x) x=*p;
		p++;
	}
	return x;
}

double min(const matrix &a) {

	double x,*p;
	int i,N;
	
	x=*(a.p);
	N=a.nf*a.nc;
	p=a.p;
	for(i=0;i<N;i++) {
		if (*p<x) x=*p;
		p++;
	}
	return x;
}

double sum(const matrix &a) {

	double s=0,*p;
	int i,N;
	
	N=a.nf*a.nc;
	p=a.p;
	for(i=0;i<N;i++) 
		s+=*(p+i);

	return s;
}

double mean(const matrix &a) {

	return sum(a)/a.nf/a.nc;
}

int exist(const matrix &a) {

	int i,N;
	double *p;
	N=a.nf*a.nc;
	p=a.p;
	for(i=0;i<N;i++)
		if(*(p+i)) return 1;
	return 0;
}

int isequal(const matrix &a,const matrix &b) {

	int i,N,res;
	double *p,*pb;
	if(a.nf!=b.nf) return 0;
	if(a.nc!=b.nc) return 0;
	N=a.nf*a.nc;
	p=a.p;
	pb=b.p;
	res=1;
	for(i=0;i<N;i++)
		if(*(p+i)!=*(pb+i)) res=0;
	return res;
}

matrix max(const matrix &a,const matrix &b) {

	if( (b.nf!=a.nf) || (b.nc!=a.nc) ) {
		fprintf(stderr,"ERROR: (matrix.max) Dimensions must agree\n");
		exit(1);
	}

    matrix res(a.nf,a.nc);
    double *pa,*pres,*pb;
    int i,N;
    pa=a.p;
    pb=b.p;
    pres=res.p;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	*(pres+i)=*(pa+i)>*(pb+i)?*(pa+i):*(pb+i);
    return res;
}

matrix max(const matrix &a,double n) {

    matrix res(a.nf,a.nc);
    double *pa,*pres;
    int i,N;
    pa=a.p;
    pres=res.p;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	*(pres+i)=*(pa+i)>n?*(pa+i):n;
    return res;
}

matrix max(double n,const matrix &a) {
	return max(a,n);
}

matrix min(const matrix &a,const matrix &b) {

	if( (b.nf!=a.nf) || (b.nc!=a.nc) ) {
		fprintf(stderr,"ERROR: (matrix.min) Dimensions must agree\n");
		exit(1);
	}

    matrix res(a.nf,a.nc);
    double *pa,*pres,*pb;
    int i,N;
    pa=a.p;
    pb=b.p;
    pres=res.p;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	*(pres+i)=*(pa+i)<*(pb+i)?*(pa+i):*(pb+i);
    return res;
}

matrix min(const matrix &a,double n) {

    matrix res(a.nf,a.nc);
    double *pa,*pres;
    int i,N;
    pa=a.p;
    pres=res.p;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
    	*(pres+i)=*(pa+i)<n?*(pa+i):n;
    return res;
}

matrix min(double n,const matrix &a) {
	return min(a,n);
}
