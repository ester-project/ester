#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "utils.h"
#include "matrix.h"
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include <stdexcept>
#include <vector>

/// \brief Default matrix constructor.
///
/// Creates an nfil rows by ncol columns matrix.
/// If \p nfil and \p ncol are omitted, defaults are 1.
matrix::matrix(int nfil, int ncol) {

    unsigned tam;

    if(nfil<0||ncol<0) {
        ester_err("Can't create matrix with negative size");
    }
    if(nfil==0) {
        ester_err("Number of rows can't be zero");
    }
    if(ncol==0) {
        ester_err("Number of columns can't be zero");
    }
    nf=nfil;
    nc=ncol;
    tam=unsigned(nf)*unsigned(nc);
    p=new double[tam];
}

/// \brief Destructor
matrix::~matrix() {

    delete [] p;
}


/// \brief Copy constructor
///
/// Creates a new matrix object identical to \p a.
matrix::matrix(const matrix &a) {

    unsigned tam;

    nc=a.nc;
    nf=a.nf;
    tam=unsigned(nf)*unsigned(nc);
    p=new double[tam];

    memcpy(p, a.p, nc*nf*sizeof(double));

}


/// \brief Returns the number of rows of the matrix.
int matrix::nrows() const {

    return nf;
}


/// \brief Returns the columns of rows of the matrix.
int matrix::ncols() const {

    return nc;
}

/// \brief Returns the pointer to the matrix data.
double *matrix::data() const {

    return p;
}


/// \brief Re-dimension the matrix.
///
/// Updates the matrix dimension.
/// \returns the updated matrix object.
matrix &matrix::dim(int nfil, int ncol) {

    unsigned tam;

    if(nfil<0||ncol<0) {
        ester_err("Can't create matrix with negative size");
        exit(1);
    }
    if(nfil==0) {
        ester_err("Number of rows can't be zero");
        exit(1);
    }
    if(ncol==0) {
        ester_err("Number of columns can't be zero");
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


/// \brief Re-dimension the matrix.
///
/// Updates the matrix dimension. The new full matrix (\p nfil x \p ncol) size
/// must match the original matrix size (nrows() x ncols()).
/// \returns the updated matrix object.
matrix &matrix::redim(int nfil, int ncol) {

    if(nfil<0||ncol<0) {
        ester_err("Can't create matrix with negative size");
        exit(1);
    }
    if(nfil==0) {
        ester_err("Number of rows can't be zero");
        exit(1);
    }
    if(ncol==0) {
        ester_err("Number of columns can't be zero");
        exit(1);
    }

    if(nfil*ncol!=nf*nc) {
        ester_err("matrix.redim: Number of elements doesn't match");
        exit(1);
    }
    nf=nfil;nc=ncol;

    return *this;
}


/// \brief Affectation operator.
matrix &matrix::operator=(const matrix &a) {

    if(&a==this) return *this;
    if (nf!=a.nf||nc!=a.nc)
        dim(a.nf, a.nc);

    memcpy(p, a.p, nc*nf*sizeof(double));

    return *this;
}


/// \brief Access matrix's elements.
double &matrix::operator()(int ifil, int icol) {

    if(ifil<0) ifil+=nf;
    if(icol<0) icol+=nc;
    if(ifil>=nf||ifil<0||icol>=nc||icol<0) {
        ester_err("Index exceeds matrix dimensions");
        exit(1);
    }
    return *(p+icol*nf+ifil);

}


/// \brief Access matrix's elements.
double &matrix::operator()(int ielem) {

    if(ielem<0) ielem+=nf*nc;
    if(ielem>=nf*nc||ielem<0) {
        ester_err("Index exceeds matrix dimensions");
        exit(1);
    }
    return *(p+ielem);

}


/// \brief Access matrix's elements.
const double &matrix::operator()(int ifil, int icol) const {

    if(ifil<0) ifil+=nf;
    if(icol<0) icol+=nc;
    if(ifil>=nf||ifil<0||icol>=nc||icol<0) {
        ester_err("Index exceeds matrix dimensions");
        exit(1);
    }
    return *(p+icol*nf+ifil);

}


/// \brief Access matrix's elements.
const double &matrix::operator()(int ielem) const {

    if(ielem<0) ielem+=nf*nc;
    if(ielem>=nf*nc||ielem<0) {
        ester_err("Index exceeds matrix dimensions");
        exit(1);
    }
    return *(p+ielem);

}


/// \brief Reads matrix data from a file.
int matrix::read(int nfil, int ncol, FILE *fp, char mode) {

    int i, tam;

    dim(nfil, ncol);
    tam=nc*nf;

    if (mode=='b') {
        if((int)fread(p, sizeof(double), tam, fp)!=tam) return 0;
        }
    else
        for(i=0;i<tam;i++) {
            if(fscanf(fp, "%le", &p[i])==EOF) return 0;
            }
    return 1;
}

/// \brief Saves matrix to a file.
int matrix::write(FILE *fp, char mode) const {

    int tam, i, j;

    tam=nf*nc;
    if (mode=='b') {
        if((int)fwrite(p, sizeof(double), tam, fp)!=tam) return 0;
        }
    else {
        j=0;
        for(i=0;i<tam;i++) {
            if(fprintf(fp, "%.16e ", p[i])==EOF) return 0;
            j++;
            if (j==nf) {
                j=0;
                fprintf(fp, "\n");
            }
        }
    }
    return 1;
}


/// \brief Saves matrix to a file with a special format.
///
/// \p fmt is the printf format used to save the matrix data to the file.
void matrix::write_fmt(const char *fmt, FILE *fp) const {

    int i, j;

    for(i=0;i<nf;i++) {
        for(j=0;j<nc;j++) {
            fprintf(fp, fmt, (*this)(i, j));
            if(j<nc-1) fprintf(fp, " ");
        }
        fprintf(fp, "\n");
    }
}


/// \brief Swaps two matrices.
void matrix::swap(matrix &a) {

    double *p0;
    int nf0, nc0;

    p0=p;
    nf0=nf;nc0=nc;

    nf=a.nf;nc=a.nc;
    p=a.p;
    a.nf=nf0;a.nc=nc0;
    a.p=p0;

}


/// \brief Sets the matrix to the \p nrows x \p ncols zero matrix.
void matrix::zero(int nrows, int ncols) {

    dim(nrows, ncols);
    for(int i=0;i<nf*nc;i++) p[i]=0;

}


/// \brief Sets the values of the matrix (given in column major order).
void matrix::values(double x, ...) {

    va_list ap;

    *this=this->transpose();
    va_start(ap, x);
    *p=x;
    for(int i=1;i<nf*nc;i++) {
        *(p+i)=va_arg(ap, double);
    }
    va_end(ap);
    *this=this->transpose();

}

//
/// \brief Performs matrix addition.
///
/// The matrices must be of the same size.
/// \returns Result of the element-wise addition of the two matrices.
matrix matrix::operator+(const matrix &a) const {

    matrix res;
    int i, j, resnf, resnc, N;

    if( (nf!=1 && a.nf!=1 && nf!=a.nf) || (nc!=1 && a.nc!=1 && nc!=a.nc) ) {
        ester_err("(matrix.+) Dimensions must agree");
        exit(1);
    }

    if(nf>a.nf) resnf=nf; // This cannot happen
        else resnf=a.nf;
    if(nc>a.nc) resnc=nc;
        else resnc=a.nc;

    res.dim(resnf, resnc);
    if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
        N=nc*nf;
        for(i=0;i<N;i++)
               res.p[i]=p[i]+a.p[i];
    } else {
        double *pi, *pa, *pres;
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


matrix matrix::operator+(double n) const {

    matrix res(nf, nc);
    int i, N;

    N=nc*nf;
    for(i=0;i<N;i++)
        res.p[i]=p[i]+n;

    return res;
}


matrix matrix::operator-(const matrix &a) const {

    matrix res;
    int i, j, resnf, resnc, N;

    if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
        ester_err("(matrix.-) Dimensions must agree");
        exit(1);
    }

    if(nf>a.nf) resnf=nf;
        else resnf=a.nf;
    if(nc>a.nc) resnc=nc;
        else resnc=a.nc;

    res.dim(resnf, resnc);
    if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
        N=nc*nf;
        for(i=0;i<N;i++)
               res.p[i]=p[i]-a.p[i];
    } else {
        double *pi, *pa, *pres;
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


matrix operator-(double n, const matrix &a) {

    matrix res(a.nf, a.nc);
    int i, N;

    N=a.nc*a.nf;
    for(i=0;i<N;i++)
        res.p[i]=n-a.p[i];
    return res;
}


matrix matrix::operator-(double n) const {

    matrix res(nf, nc);
    int i, N;

    N=nc*nf;
    for(i=0;i<N;i++)
        res.p[i]=p[i]-n;

    return res;

}


matrix operator-(const matrix &a) {

    matrix res(a.nf, a.nc);
    int i, N;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
        res.p[i]=-a.p[i];

    return res;
}


matrix matrix::operator*(const matrix &a) const {

    matrix res;
    int i, j, resnf, resnc, N;

    if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
        ester_err("(matrix.*) Dimensions must agree");
        exit(1);
    }

    if(nf>a.nf) resnf=nf;
        else resnf=a.nf;
    if(nc>a.nc) resnc=nc;
        else resnc=a.nc;

    res.dim(resnf, resnc);
    if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
        N=nf*nc;
        for(i=0;i<N;i++)
               res.p[i]=p[i]*a.p[i];
    } else {
        double *pi, *pa, *pres;
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


matrix matrix::operator*(double n) const {

    matrix res(nf, nc);
    int i, N;

    N=nc*nf;
    for(i=0;i<N;i++)
        res.p[i]=p[i]*n;

    return res;
}


matrix matrix::operator/(const matrix &a) const {

    matrix res;
    int i, j, resnf, resnc, N;

    if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
        ester_err("(matrix./) Dimensions must agree");
        exit(1);
    }

    if(nf>a.nf) resnf=nf;
        else resnf=a.nf;
    if(nc>a.nc) resnc=nc;
        else resnc=a.nc;

    res.dim(resnf, resnc);
    if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
        N=nc*nf;
        for(i=0;i<N;i++)
               res.p[i]=p[i]/a.p[i];
    } else {
        double *pi, *pa, *pres;
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


matrix operator/(double n, const matrix &a) {

    matrix res(a.nf, a.nc);
    int i, N;

    N=a.nc*a.nf;
    for(i=0;i<N;i++)
        res.p[i]=n/a.p[i];
    return res;
}


matrix matrix::operator/(double n) const {

    matrix res(nf, nc);
    int i, N;

    N=nc*nf;
    for(i=0;i<N;i++)
        res.p[i]=p[i]/n;

    return res;
}


matrix matrix::operator==(const matrix &a) const {

    matrix res;
    int i, j, resnf, resnc, N;

    if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
        ester_err("(matrix.==) Dimensions must agree");
        exit(1);
    }

    if(nf>a.nf) resnf=nf;
        else resnf=a.nf;
    if(nc>a.nc) resnc=nc;
        else resnc=a.nc;

    res.dim(resnf, resnc);
    if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
        N=nc*nf;
        for(i=0;i<N;i++)
               res.p[i]=(p[i]==a.p[i]);
    } else {
        double *pi, *pa, *pres;
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


matrix matrix::operator==(double n) const {

    matrix res(nf, nc);
    int i, N;

    N=nc*nf;
    for(i=0;i<N;i++)
        res.p[i]=(p[i]==n);

    return res;
}


matrix matrix::operator!=(const matrix &a) const {

    matrix res;
    int i, j, resnf, resnc, N;

    if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
        ester_err("(matrix.!=) Dimensions must agree");
        exit(1);
    }

    if(nf>a.nf) resnf=nf;
        else resnf=a.nf;
    if(nc>a.nc) resnc=nc;
        else resnc=a.nc;

    res.dim(resnf, resnc);
    if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
        N=nc*nf;
        for(i=0;i<N;i++)
               res.p[i]=(p[i]!=a.p[i]);
    } else {
        double *pi, *pa, *pres;
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


matrix matrix::operator!=(double n) const {

    matrix res(nf, nc);
    int i, N;

    N=nc*nf;
    for(i=0;i<N;i++)
        res.p[i]=(p[i]!=n);

    return res;}


matrix matrix::operator>(const matrix &a) const{

    matrix res;
    int i, j, resnf, resnc, N;

    if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
        ester_err("(matrix.>) Dimensions must agree");
        exit(1);
    }

    if(nf>a.nf) resnf=nf;
        else resnf=a.nf;
    if(nc>a.nc) resnc=nc;
        else resnc=a.nc;

    res.dim(resnf, resnc);
    if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
        N=nc*nf;
        for(i=0;i<N;i++)
            res.p[i]=(p[i]>a.p[i]);
    } else {
        double *pi, *pa, *pres;
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


matrix matrix::operator>(double n) const {

    matrix res(nf, nc);
    int i, N;

    N=nc*nf;
    for(i=0;i<N;i++)
        res.p[i]=(p[i]>n);

    return res;
}


matrix matrix::operator<(double n) const {

    matrix res(nf, nc);
    int i, N;

    N=nc*nf;
    for(i=0;i<N;i++)
        res.p[i]=(p[i]<n);

    return res;
}


matrix matrix::operator>=(const matrix &a) const {

    matrix res;
    int i, j, resnf, resnc, N;

    if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
        ester_err("(matrix.>=) Dimensions must agree");
        exit(1);
    }

    if(nf>a.nf) resnf=nf;
        else resnf=a.nf;
    if(nc>a.nc) resnc=nc;
        else resnc=a.nc;

    res.dim(resnf, resnc);
    if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
        N=nc*nf;
        for(i=0;i<N;i++)
               res.p[i]=(p[i]>=a.p[i]);
    } else {
        double *pi, *pa, *pres;
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


matrix matrix::operator>=(double n) const {

    matrix res(nf, nc);
    int i, N;

    N=nc*nf;
    for(i=0;i<N;i++)
        res.p[i]=(p[i]>=n);

    return res;
}


matrix matrix::operator<=(double n) const {

    matrix res(nf, nc);
    int i, N;

    N=nc*nf;
    for(i=0;i<N;i++)
        res.p[i]=(p[i]<=n);

    return res;
}


matrix matrix::operator||(const matrix &a) const {

    matrix res;
    int i, j, resnf, resnc, N;

    if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
        ester_err("(matrix.||) Dimensions must agree");
        exit(1);
    }

    if(nf>a.nf) resnf=nf;
        else resnf=a.nf;
    if(nc>a.nc) resnc=nc;
        else resnc=a.nc;

    res.dim(resnf, resnc);
    if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
        N=nc*nf;
        for(i=0;i<N;i++)
               res.p[i]=(p[i]||a.p[i]);
    } else {
        double *pi, *pa, *pres;
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


matrix matrix::operator||(double n) const {

    matrix res(nf, nc);
    int i, N;

    N=nc*nf;
    for(i=0;i<N;i++)
        res.p[i]=(p[i]||n);

    return res;
}


matrix matrix::operator&&(const matrix &a) const {

    matrix res;
    int i, j, resnf, resnc, N;

    if( (nf!=1&&a.nf!=1&&nf!=a.nf) || (nc!=1&&a.nc!=1&&nc!=a.nc) ) {
        ester_err("(matrix.&&) Dimensions must agree");
        exit(1);
    }

    if(nf>a.nf) resnf=nf;
        else resnf=a.nf;
    if(nc>a.nc) resnc=nc;
        else resnc=a.nc;

    res.dim(resnf, resnc);
    if (nf==resnf&&nc==resnc&&a.nf==resnf&&a.nc==resnc) {
        N=nc*nf;
        for(i=0;i<N;i++)
               res.p[i]=(p[i]&&a.p[i]);
    } else {
        double *pi, *pa, *pres;
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


matrix matrix::operator&&(double n) const {

    matrix res(nf, nc);
    int i, N;

    N=nc*nf;
    for(i=0;i<N;i++)
        res.p[i]=(p[i]&&n);

    return res;
}


matrix &matrix::operator+=(const matrix &a) {

    int i, j, k, N;

    N=nc*nf;
    if(nf==a.nf&&nc==a.nc) {
        for(i=0;i<N;i++)
            p[i]+=a.p[i];
    } else if(nf==a.nf&&a.nc==1) {
        for(i=0, j=0;i<N;i++, j=(j+1)%nf)
            p[i]+=a.p[j];
    } else if(a.nf==1&&nc==a.nc) {
        k=0;
        for(j=0;j<nc;j++)
            for(i=0;i<nf;i++, k++)
                p[k]+=a.p[j];
    } else if(a.nf==1&&a.nc==1) {
        *this+=*a.p;
    } else {
        *this=(*this)+a;
    }

    return *this;
}


matrix &matrix::operator-=(const matrix &a) {

    int i, j, k, N;

    N=nc*nf;
    if(nf==a.nf&&nc==a.nc) {
        for(i=0;i<N;i++)
            p[i]-=a.p[i];
    } else if(nf==a.nf&&a.nc==1) {
        for(i=0, j=0;i<N;i++, j=(j+1)%nf)
            p[i]-=a.p[j];
    } else if(a.nf==1&&nc==a.nc) {
        k=0;
        for(j=0;j<nc;j++)
            for(i=0;i<nf;i++, k++)
                p[k]-=a.p[j];
    } else if(a.nf==1&&a.nc==1) {
        *this-=*a.p;
    } else {
        *this=(*this)-a;
    }

    return *this;
}


matrix &matrix::operator*=(const matrix &a) {

    int i, j, k, N;

    N=nc*nf;
    if(nf==a.nf&&nc==a.nc) {
        for(i=0;i<N;i++)
            p[i]*=a.p[i];
    } else if(nf==a.nf&&a.nc==1) {
        for(i=0, j=0;i<N;i++, j=(j+1)%nf)
            p[i]*=a.p[j];
    } else if(a.nf==1&&nc==a.nc) {
        k=0;
        for(j=0;j<nc;j++)
            for(i=0;i<nf;i++, k++)
                p[k]*=a.p[j];
    } else if(a.nf==1&&a.nc==1) {
        *this*=*a.p;
    } else {
        *this=(*this)*a;
    }

    return *this;
}


matrix &matrix::operator/=(const matrix &a) {

    int i, j, k, N;

    N=nc*nf;
    if(nf==a.nf&&nc==a.nc) {
        for(i=0;i<N;i++)
            p[i]/=a.p[i];
    } else if(nf==a.nf&&a.nc==1) {
        for(i=0, j=0;i<N;i++, j=(j+1)%nf)
            p[i]/=a.p[j];
    } else if(a.nf==1&&nc==a.nc) {
        k=0;
        for(j=0;j<nc;j++)
            for(i=0;i<nf;i++, k++)
                p[k]/=a.p[j];
    } else if(a.nf==1&&a.nc==1) {
        *this/=*a.p;
    } else {
        *this=(*this)/a;
    }

    return *this;
}


matrix &matrix::operator+=(double n) {

    int i, N;
    N=nc*nf;
    for(i=0;i<N;i++)
        p[i]+=n;
    return *this;
}


matrix &matrix::operator-=(double n) {

    int i, N;
    N=nc*nf;
    for(i=0;i<N;i++)
        p[i]-=n;
    return *this;
}


matrix &matrix::operator*=(double n) {

    int i, N;
    N=nc*nf;
    for(i=0;i<N;i++)
        p[i]*=n;
    return *this;
}


matrix &matrix::operator/=(double n) {

    int i, N;
    N=nc*nf;
    for(i=0;i<N;i++)
        p[i]/=n;
    return *this;
}

const matrix matrix::row(int ifil) const {

    matrix res(1, nc);
    double *pi;
    int i;

    if(ifil<0) ifil+=nf;
    if(ifil<0||ifil>=nf) {
        ester_err("(matrix.row) Index exceeds matrix dimensions");
        exit(1);
    }

    pi=p+ifil;
    for(i=0;i<nc;i++)
        res.p[i]=pi[i*nf];

    return res;
}


matrix &matrix::setrow(int ifil, const matrix &a) {

    double *pi;
    int i;

    if(ifil<0) ifil+=nf;
    if(ifil<0||ifil>=nf) {
        ester_err("(matrix.setrow) Index exceeds matrix dimensions");
        exit(1);
    }
    if(a.nf>1||a.nc!=nc) {
        ester_err("(matrix.setrow) Dimensions must agree");
        exit(1);
    }

    pi=p+ifil;
    for(i=0;i<nc;i++)
        pi[i*nf]=a.p[i];

    return *this;

}


const matrix matrix::col(int icol) const {

    matrix res(nf, 1);
    double *pi;
    int i, N;

    if(icol<0) icol+=nc;
    if(icol<0||icol>=nc) {
        ester_err("(matrix.col) Index exceeds matrix dimensions");
        exit(1);
    }

    pi=p+icol*nf;
    N=nf;
    for(i=0;i<N;i++)
        res.p[i]=pi[i];

    return res;
}

matrix &matrix::setcol(int icol, const matrix &a) {

    double *pi;
    int i, N;

    if(icol<0) icol+=nc;
    if(icol<0||icol>=nc) {
        ester_err("(matrix.setcol) Index exceeds matrix dimensions");
        exit(1);
    }
    if(a.nc>1||a.nf!=nf) {
        ester_err("(matrix.setcol) Dimensions must agree");
        exit(1);
    }

    pi=p+icol*nf;
    N=nf;
    for(i=0;i<N;i++)
        pi[i]=a.p[i];

    return *this;

}

const matrix matrix::block(int ifil1, int ifil2, int icol1, int icol2) const {

    if(ifil1<0) ifil1+=nf;
    if(ifil2<0) ifil2+=nf;
    if(icol1<0) icol1+=nc;
    if(icol2<0) icol2+=nc;

    if(ifil1<0||ifil1>=nf||ifil2<0||ifil2>=nf||icol1<0||icol1>=nc||icol2<0||icol2>=nc) {
        ester_err("(matrix.block) Index exceeds matrix dimensions");
        exit(1);
    }

    matrix res(ifil2-ifil1+1, icol2-icol1+1);
    double *pi, *pres;
    int i, j, N1, N2;

    pi=p+ifil1+icol1*nf;pres=res.p;
    N1=res.nc;N2=res.nf;
    for(i=0;i<N1;i++) {
        for(j=0;j<N2;j++)
            pres[j]=pi[j];
        pres+=res.nf;pi+=nf;
    }

    return res;
}


const matrix matrix::block_step(int ifil1, int ifil2, int dfil, int icol1, int icol2, int dcol) const {

    if(ifil1<0) ifil1+=nf;
    if(ifil2<0) ifil2+=nf;
    if(icol1<0) icol1+=nc;
    if(icol2<0) icol2+=nc;

    if(ifil1<0||ifil1>=nf||ifil2<0||ifil2>=nf||icol1<0||icol1>=nc||icol2<0||icol2>=nc) {
        ester_err("(matrix.block_step) Index exceeds matrix dimensions");
        exit(1);
    }

    matrix res((int) floor((ifil2-ifil1)/dfil)+1,
            (int) floor((icol2-icol1)/dcol)+1);
    double *pi, *pres;
    int i, j, N1, N2;

    pi=p+ifil1+icol1*nf;pres=res.p;
    N1=res.nc;N2=res.nf;
    for(i=0;i<N1;i++) {
        for(j=0;j<N2;j++)
            pres[j]=pi[j*dfil];
        pres+=res.nf;pi+=dcol*nf;
    }

    return res;
}

matrix &matrix::setblock(int ifil1, int ifil2, int icol1, int icol2, const matrix &a) {

    if(ifil1<0) ifil1+=nf;
    if(ifil2<0) ifil2+=nf;
    if(icol1<0) icol1+=nc;
    if(icol2<0) icol2+=nc;

    if(ifil1<0||ifil1>=nf||ifil2<0||ifil2>=nf||icol1<0||icol1>=nc||icol2<0||icol2>=nc) {
        ester_err("(matrix.setblock) Index exceeds matrix dimensions");
        exit(1);
    }
    if(a.nf!=ifil2-ifil1+1||a.nc!=icol2-icol1+1) {
        ester_err("(matrix.setblock) Dimensions must agree");
        exit(1);
    }

    double *pi, *pa;
    int i, j, N1, N2;

    pi=p+ifil1+icol1*nf;pa=a.p;
    N1=a.nc;N2=a.nf;
    for(i=0;i<N1;i++) {
        for(j=0;j<N2;j++)
            pi[j]=pa[j];
        pa+=a.nf;pi+=nf;
    }

    return *this;

}


matrix &matrix::setblock_step(int ifil1, int ifil2, int dfil, int icol1, int icol2, int dcol, const matrix &a) {

    if(ifil1<0) ifil1+=nf;
    if(ifil2<0) ifil2+=nf;
    if(icol1<0) icol1+=nc;
    if(icol2<0) icol2+=nc;

    if(ifil1<0||ifil1>=nf||ifil2<0||ifil2>=nf||icol1<0||icol1>=nc||icol2<0||icol2>=nc) {
        ester_err("(matrix.setblock_step) Index exceeds matrix dimensions");
        exit(1);
    }
    if(a.nf!=floor((ifil2-ifil1)/dfil)+1||a.nc!=floor((icol2-icol1)/dcol)+1) {
        ester_err("(matrix.setblock_step) Dimensions must agree");
        exit(1);
    }

    double *pi, *pa;
    int i, j, N1, N2;

    pi=p+ifil1+icol1*nf;pa=a.p;
    N1=a.nc;N2=a.nf;
    for(i=0;i<N1;i++) {
        for(j=0;j<N2;j++)
            pi[j*dfil]=pa[j];
        pa+=a.nf;pi+=dcol*nf;
    }

    return *this;

}


matrix matrix::concatenate(const matrix &a, int dir) const {

    matrix res;

    if(dir==0) {
        if(a.nc!=nc) {
            ester_err("(matrix.concatenate) Dimensions must agree (%d != %d)",
                    a.nc, nc);
            exit(1);
        }
        res.dim(nf+a.nf, nc);
        res.setblock(0, nf-1, 0, -1, *this);
        res.setblock(nf, a.nf+nf-1, 0, -1, a);
    } else {
        if(a.nf!=nf) {
            ester_err("(matrix.concatenate) Dimensions must agree");
            exit(1);
        }
        res.dim(nf, nc+a.nc);
        res.setblock(0, -1, 0, nc-1, *this);
        res.setblock(0, -1, nc, nc+a.nc-1, a);
    }

    return res;

}

matrix matrix::transpose() const {

    matrix a(nc, nf);
    int i, j, N1, N2;
    double *pi, *pa;

    N1=nf;N2=nc;
    pi=p;pa=a.p;

    for(i=0;i<N1;i++) {
        for(j=0;j<N2;j++)
            pa[j]=pi[j*nf];
        pa+=a.nf;pi++;
    }

    return a;
}

matrix matrix::fliplr() const {

    matrix a(nf, nc);
    int i, j, N1, N2;
    double *pi, *pa;

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


matrix matrix::flipud() const {

    matrix a(nf, nc);
    int i, j, N1, N2;
    double *pi, *pa;

    N1=nc;N2=nf;
    pi=p;pa=a.p;

    for(i=0;i<N1;i++) {
        for(j=0;j<N2;j++)
            pa[j]=pi[nf-1-j];
        pa+=a.nf;pi+=nf;
    }

    return a;
}

matrix ones(int nfil, int ncol) {

    matrix a(nfil, ncol);
    int i, N;

    N=nfil*ncol;
    for(i=0;i<N;i++)
        a.p[i]=1;

    return a;
}

matrix zeros(int nfil, int ncol) {

    matrix a(nfil, ncol);
    int i, N;

    N=nfil*ncol;
    for(i=0;i<N;i++)
        a.p[i]=0;

    return a;
}

matrix random_matrix(int nfil, int ncol) {

    matrix a(nfil, ncol);
    int i, N;

    N=nfil*ncol;
    for(i=0;i<N;i++)
        a.p[i]=(double) rand()/RAND_MAX;

    return a;

}


matrix eye(int n) {

    matrix a;
    int i, d, N;

    a.zero(n, n);
    d=n+1;
    N=n*n;
    for(i=0;i<N;i+=d)
        a.p[i]=1;
    return a;
}


matrix vector_t(double x0, double x1, int n) {

    double dx;
    int i;
    matrix v(n, 1);

    if (n == 1) {
        v(0) = (x1+x0)/2.0;
        return v;
    }

    dx=(x1-x0)/(n-1);

    for(i=0;i<n;i++)
        v.p[i]=x0+i*dx;
    v.p[n-1]=x1;

    return v;
}

matrix vector(double x0, double x1, int n) {

    matrix v;

    v=vector_t(x0, x1, n);
    v.nc=v.nf;
    v.nf=1;

    return v;
}


double max(const matrix &a) {

    double x;
    int i, N;

    x=*(a.p);
    N=a.nf*a.nc;
    for(i=1;i<N;i++)
        if (a.p[i]>x) x=a.p[i];

    return x;
}


double min(const matrix &a) {

    double x;
    int i, N;

    x=*(a.p);
    N=a.nf*a.nc;
    for(i=1;i<N;i++)
        if (a.p[i]<x) x=a.p[i];

    return x;
}


double sum(const matrix &a) {

    double s=0;
    int i, N;

    N=a.nf*a.nc;
    for(i=0;i<N;i++)
        s+=a.p[i];

    return s;
}


double mean(const matrix &a) {

    return sum(a)/a.nf/a.nc;
}


matrix max(const matrix &a, const matrix &b) {

    if( (b.nf!=a.nf) || (b.nc!=a.nc) ) {
        ester_err("(matrix.max) Dimensions must agree");
        exit(1);
    }

    matrix res(a.nf, a.nc);
    int i, N;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
        res.p[i]=a.p[i]>b.p[i]?a.p[i]:b.p[i];
    return res;
}


matrix max(const matrix &a, double n) {

    matrix res(a.nf, a.nc);
    int i, N;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
        res.p[i]=a.p[i]>n?a.p[i]:n;
    return res;
}


matrix max(double n, const matrix &a) {
    return max(a, n);
}


matrix min(const matrix &a, const matrix &b) {

    if( (b.nf!=a.nf) || (b.nc!=a.nc) ) {
        ester_err("(matrix.min) Dimensions must agree");
        exit(1);
    }

    matrix res(a.nf, a.nc);
    int i, N;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
        res.p[i]=a.p[i]<b.p[i]?a.p[i]:b.p[i];
    return res;
}


matrix min(const matrix &a, double n) {

    matrix res(a.nf, a.nc);
    int i, N;
    N=a.nc*a.nf;
    for(i=0;i<N;i++)
        res.p[i]=a.p[i]<n?a.p[i]:n;
    return res;
}


matrix min(double n, const matrix &a) {
    return min(a, n);
}


int exist(const matrix &a) {

    int i, N;
    N=a.nf*a.nc;
    for(i=0;i<N;i++)
        if(a.p[i]) return 1;
    return 0;
}


int isequal(const matrix &a, const matrix &b) {

    int i, N, res;
    if(a.nf!=b.nf) return 0;
    if(a.nc!=b.nc) return 0;
    N=a.nf*a.nc;
    res=1;
    for(i=0;i<N;i++)
        if(a.p[i]!=b.p[i]) res=0;
    return res;
}




