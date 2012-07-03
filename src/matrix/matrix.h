#ifndef _MATRIX_H
#define _MATRIX_H

#include<stdio.h>

class matrix {
	double *p;
	int nf,nc;
  public:
    matrix(int nfil=1,int ncol=1);
    ~matrix();
    matrix(const matrix &);
    int nrows() const;
    int ncols() const;
    double *data() const;
    matrix &dim(int nfil,int ncol);
    matrix &redim(int nfil,int ncol);
    matrix &operator=(const matrix &);
    double &operator()(int ifil,int icol=0) const; 
    int read(int nfil,int ncol,FILE *fp,char mode='t');
    int write(FILE *fp=stdout,char mode='t') const;
    
    matrix operator+(const matrix &) const;
    friend matrix operator+(double,const matrix &);
    matrix operator+(double) const;
    friend matrix operator+(const matrix &);
    matrix operator-(const matrix &) const;
    friend matrix operator-(double,const matrix &);
    matrix operator-(double) const;
    friend matrix operator-(const matrix &);
    matrix operator*(const matrix &) const;
    friend matrix operator*(double,const matrix &);
    matrix operator*(double) const;
    matrix operator/(const matrix &) const;
    friend matrix operator/(double,const matrix &);
    matrix operator/(double) const;
    matrix operator==(const matrix &) const;
    friend matrix operator==(double,const matrix &);
    matrix operator==(double) const;
    matrix operator!=(const matrix &) const;
    friend matrix operator!=(double,const matrix &);
    matrix operator!=(double) const;
    matrix operator<(const matrix &) const;
    friend matrix operator<(double,const matrix &);
    matrix operator<(double) const;
    matrix operator>(const matrix &) const;
    friend matrix operator>(double,const matrix &);
    matrix operator>(double) const;
    matrix operator<=(const matrix &) const;
    friend matrix operator<=(double,const matrix &);
    matrix operator<=(double) const;
    matrix operator>=(const matrix &) const;
    friend matrix operator>=(double,const matrix &);
    matrix operator>=(double) const;
    matrix operator&&(const matrix &) const;
    friend matrix operator&&(double,const matrix &);
    matrix operator&&(double) const;
    matrix operator||(const matrix &) const;
    friend matrix operator||(double,const matrix &);
    matrix operator||(double) const;
    matrix operator,(const matrix &) const;
    matrix &operator+=(const matrix &);
    matrix &operator-=(const matrix &);
    matrix &operator*=(const matrix &);
    matrix &operator/=(const matrix &);
    matrix &operator+=(double);
    matrix &operator-=(double);
    matrix &operator*=(double);
    matrix &operator/=(double);
    friend double max(const matrix &);
    friend double min(const matrix &);
    friend double sum(const matrix &);
    friend double mean(const matrix &);
    friend matrix max(const matrix &,const matrix &);
    friend matrix max(const matrix &,double);
    friend matrix max(double,const matrix &);
    friend matrix min(const matrix &,const matrix &);
    friend matrix min(const matrix &,double);
    friend matrix min(double,const matrix &);
    friend matrix eye(int n);
    friend matrix ones(int nfil,int ncol);
    friend matrix zeros(int nfil,int ncol);
    friend matrix random_matrix(int nfil,int ncol);
    friend matrix vector(double x0,double x1,int n);
    friend matrix vector_t(double x0,double x1,int n);
    matrix transpose() const;
    matrix fliplr() const;
    matrix flipud() const;
    friend int exist(const matrix &);
    friend int isequal(const matrix &,const matrix &);
    
    matrix row(int ifil) const;
    matrix col(int icol) const;
    matrix block(int ifil1,int ifil2,int icol1,int icol2) const;
    matrix block_step(int ifil1,int ifil2,int dfil,int icol1,int icol2,int dcol) const;
    matrix &setrow(int ifil,const matrix &);
    matrix &setcol(int icol,const matrix &);
    matrix &setblock(int ifil1,int ifil2,int icol1,int icol2,const matrix &);
    matrix &setblock_step(int ifil1,int ifil2,int dfil,int icol1,int icol2,int dcol,const matrix &);
      
    friend matrix cos(const matrix &);
    friend matrix sin(const matrix &);
    friend matrix tan(const matrix &);
    friend matrix acos(const matrix &);
    friend matrix asin(const matrix &);
    friend matrix atan(const matrix &);
    friend matrix cosh(const matrix &);
    friend matrix sinh(const matrix &);
    friend matrix tanh(const matrix &);
    friend matrix exp(const matrix &);
    friend matrix log(const matrix &);
    friend matrix log10(const matrix &);
    friend matrix sqrt(const matrix &);
    friend matrix abs(const matrix &);
    friend matrix floor(const matrix &);
    friend matrix ceil(const matrix &);
    friend matrix round(const matrix &);
	friend matrix atan2(const matrix &,const matrix &);
	friend matrix atan2(double,const matrix &);
	friend matrix atan2(const matrix &,double);
	friend matrix pow(const matrix &,const matrix &);
	friend matrix pow(double,const matrix &);
	friend matrix pow(const matrix &,double);
	
	matrix lu(int *&ipiv) const;
	void lu2(int *&ipiv);
	matrix lusolve(matrix b,int *ipiv) const;
	matrix solve(matrix) const;
	matrix inv() const;
	
	friend class matrix_block_diag;

};

matrix operator+(double,const matrix &);
matrix operator+(const matrix &);
matrix operator-(double,const matrix &);
matrix operator-(const matrix &);
matrix operator*(double,const matrix &);
matrix operator/(double,const matrix &);
matrix operator==(double,const matrix &);
matrix operator!=(double,const matrix &);
matrix operator<(double,const matrix &);
matrix operator>(double,const matrix &);
matrix operator<=(double,const matrix &);
matrix operator>=(double,const matrix &);
matrix operator&&(double,const matrix &);
matrix operator||(double,const matrix &);
double max(const matrix &);
double min(const matrix &);
double sum(const matrix &);
double mean(const matrix &);
matrix max(const matrix &,const matrix &);
matrix max(const matrix &,double);
matrix max(double,const matrix &);
matrix min(const matrix &,const matrix &);
matrix min(const matrix &,double);
matrix min(double,const matrix &);
matrix eye(int n);
matrix ones(int nfil,int ncol);
matrix zeros(int nfil,int ncol);
matrix random_matrix(int nfil,int ncol);
matrix vector(double x0,double x1,int n);
matrix vector_t(double x0,double x1,int n);
int exist(const matrix &);
int isequal(const matrix &,const matrix &);
matrix cos(const matrix &);
matrix sin(const matrix &);
matrix tan(const matrix &);
matrix acos(const matrix &);
matrix asin(const matrix &);
matrix atan(const matrix &);
matrix cosh(const matrix &);
matrix sinh(const matrix &);
matrix tanh(const matrix &);
matrix exp(const matrix &);
matrix log(const matrix &);
matrix log10(const matrix &);
matrix sqrt(const matrix &);
matrix abs(const matrix &);
matrix floor(const matrix &);
matrix ceil(const matrix &);
matrix round(const matrix &);
matrix atan2(const matrix &,const matrix &);
matrix atan2(double,const matrix &);
matrix atan2(const matrix &,double);
matrix pow(const matrix &,const matrix &);
matrix pow(double,const matrix &);
matrix pow(const matrix &,double);


class matrix_block_diag {
	matrix *m;
	int nb;
public:
	matrix_block_diag(int nblocks=1);
	~matrix_block_diag();
	matrix_block_diag(const matrix_block_diag &);
	matrix_block_diag & operator=(const matrix_block_diag &);
	matrix_block_diag &set_nblocks(int nblocks);
	operator matrix() const;
	matrix &block(int i) const;
	int nblocks() const;
	int nrows() const;
	int ncols() const;
	matrix operator,(const matrix &) const;
	matrix_block_diag operator,(const matrix_block_diag &) const;
	friend matrix operator,(const matrix &,const matrix_block_diag &);
	matrix_block_diag operator*(const matrix &) const;
	friend matrix_block_diag operator*(const matrix &,const matrix_block_diag &);
	matrix_block_diag operator/(const matrix &) const;
	matrix_block_diag operator*(double) const;
	friend matrix_block_diag operator*(double,const matrix_block_diag &);
	matrix_block_diag operator/(double) const;
	matrix_block_diag operator+(const matrix_block_diag&) const;
	matrix_block_diag operator-(const matrix_block_diag&) const;
	friend matrix_block_diag operator+(const matrix_block_diag&);
	friend matrix_block_diag operator-(const matrix_block_diag&);
	matrix operator+(const matrix &) const;
	matrix operator-(const matrix &) const;
	matrix row(int n) const;
	double operator()(int nfil,int ncol) const;
	friend matrix_block_diag eye(const matrix_block_diag &);
	matrix_block_diag transpose() const;

};

matrix operator,(const matrix &,const matrix_block_diag &);
matrix_block_diag operator*(const matrix &,const matrix_block_diag &);
matrix_block_diag operator/(const matrix &,const matrix_block_diag &);
matrix_block_diag operator+(const matrix_block_diag&);
matrix_block_diag operator-(const matrix_block_diag&);
matrix_block_diag eye(const matrix_block_diag &);


#include<sys/time.h>
class tiempo {
	double t,t0;
	timeval tim;
public:
	int active;
	tiempo() {t=0;active=0;}
	void reset() {t=0;}
	void start() {gettimeofday(&tim,NULL);t0=tim.tv_sec+tim.tv_usec/1e6;active=1;}
	void stop() {gettimeofday(&tim,NULL);if(active) t+=tim.tv_sec+tim.tv_usec/1e6-t0;active=0;}
	double value() {double tt;gettimeofday(&tim,NULL);tt=t;if(active) tt+=tim.tv_sec+tim.tv_usec/1e6-t0;return tt;}
};

#endif
