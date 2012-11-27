#ifndef _MATRIX_H
#define _MATRIX_H

#include<stdio.h>

typedef double mat_type;

template <class Type>
class Matrix;

template <class Type>
class Matrix_block_diag;

typedef Matrix<mat_type> matrix;
typedef Matrix_block_diag<mat_type> matrix_block_diag;

template <class Type> Matrix<Type> operator+(mat_type n,const Matrix<Type> &a) {return a+n;}
template <class Type> Matrix<Type> operator+(const Matrix<Type> &a) {return a;}
template <class Type> Matrix<Type> operator-(mat_type,const Matrix<Type> &);
template <class Type> Matrix<Type> operator-(const Matrix<Type> &);
template <class Type> Matrix<Type> operator*(mat_type n,const Matrix<Type> &a) {return a*n;}
template <class Type> Matrix<Type> operator/(mat_type,const Matrix<Type> &);
template <class Type> Matrix<Type> operator==(mat_type n,const Matrix<Type> &a) {return a==n;}
template <class Type> Matrix<Type> operator!=(mat_type n,const Matrix<Type> &a) {return a!=n;}
template <class Type> Matrix<Type> operator<(mat_type n,const Matrix<Type> &a) {return a>n;}
template <class Type> Matrix<Type> operator>(mat_type n,const Matrix<Type> &a) {return a<n;}
template <class Type> Matrix<Type> operator<=(mat_type n,const Matrix<Type> &a) {return a>=n;}
template <class Type> Matrix<Type> operator>=(mat_type n,const Matrix<Type> &a) {return a<=n;}
template <class Type> Matrix<Type> operator&&(mat_type n,const Matrix<Type> &a) {return a&&n;}
template <class Type> Matrix<Type> operator||(mat_type n,const Matrix<Type> &a) {return a||n;}

matrix ones(int nfil,int ncol);
matrix zeros(int nfil,int ncol);
matrix random_matrix(int nfil,int ncol);
matrix eye(int n);
matrix vector(mat_type x0,mat_type x1,int n);
matrix vector_t(mat_type x0,mat_type x1,int n);

template <class Type> Type max(const Matrix<Type> &);
template <class Type> Type min(const Matrix<Type> &);
template <class Type> Type sum(const Matrix<Type> &);
template <class Type> Type mean(const Matrix<Type> &);

template <class Type> Matrix<Type> max(const Matrix<Type> &,const Matrix<Type> &);
template <class Type> Matrix<Type> max(const Matrix<Type> &,mat_type);
template <class Type> Matrix<Type> max(mat_type,const Matrix<Type> &);
template <class Type> Matrix<Type> min(const Matrix<Type> &,const Matrix<Type> &);
template <class Type> Matrix<Type> min(const Matrix<Type> &,mat_type);
template <class Type> Matrix<Type> min(mat_type,const Matrix<Type> &);
	
template <class Type> int exist(const Matrix<Type> &);
template <class Type> int isequal(const Matrix<Type> &,const Matrix<Type> &);

template <class Type> Matrix<Type> cos(const Matrix<Type> &);
template <class Type> Matrix<Type> sin(const Matrix<Type> &);
template <class Type> Matrix<Type> tan(const Matrix<Type> &);
template <class Type> Matrix<Type> acos(const Matrix<Type> &);
template <class Type> Matrix<Type> asin(const Matrix<Type> &);
template <class Type> Matrix<Type> atan(const Matrix<Type> &);
template <class Type> Matrix<Type> cosh(const Matrix<Type> &);
template <class Type> Matrix<Type> sinh(const Matrix<Type> &);
template <class Type> Matrix<Type> tanh(const Matrix<Type> &);
template <class Type> Matrix<Type> exp(const Matrix<Type> &);
template <class Type> Matrix<Type> log(const Matrix<Type> &);
template <class Type> Matrix<Type> log10(const Matrix<Type> &);
template <class Type> Matrix<Type> sqrt(const Matrix<Type> &);
template <class Type> Matrix<Type> abs(const Matrix<Type> &);
template <class Type> Matrix<Type> floor(const Matrix<Type> &);
template <class Type> Matrix<Type> ceil(const Matrix<Type> &);
template <class Type> Matrix<Type> round(const Matrix<Type> &);
template <class Type> Matrix<Type> atan2(const Matrix<Type> &,const Matrix<Type> &);
template <class Type> Matrix<Type> atan2(mat_type,const Matrix<Type> &);
template <class Type> Matrix<Type> atan2(const Matrix<Type> &,mat_type);
template <class Type> Matrix<Type> pow(const Matrix<Type> &,const Matrix<Type> &);
template <class Type> Matrix<Type> pow(mat_type,const Matrix<Type> &);
template <class Type> Matrix<Type> pow(const Matrix<Type> &,mat_type);
template <class Type> Matrix<Type> pow(const Matrix<Type> &,int);


template <class Type>
class Matrix {
	Type *p;
	int nf,nc;
public:
	Matrix(int nfil=1,int ncol=1);
	~Matrix();
	Matrix(const Matrix &);
	int nrows() const;
	int ncols() const;
	Type *data() const;
	Matrix &dim(int nrow,int ncol);
	Matrix &redim(int nrow,int ncol);
	Matrix &operator=(const Matrix &);
	template <class Type2> Matrix &operator=(const Matrix<Type2> &);
	template <class Type2> operator Matrix<Type2>() const;
	Type &operator()(int irow,int icol); 
	Type &operator()(int ielem);
	Type operator()(int irow,int icol) const; 
	Type operator()(int ielem) const;
	int read(int nrow,int ncol,FILE *fp,char mode='t');
	int write(FILE *fp=stdout,char mode='t') const;
	void write_fmt(const char *fmt,FILE *fp=stdout) const;
	void swap(Matrix &);
	void zero(int nrows,int ncols);
	    
	friend Matrix operator-<>(mat_type,const Matrix &);
	friend Matrix operator-<>(const Matrix &);
	friend Matrix operator/<>(mat_type,const Matrix &);
	    
	Matrix operator+(const Matrix &) const;
	Matrix operator+(mat_type) const;
	Matrix operator-(const Matrix &) const;
	Matrix operator-(mat_type) const;
	Matrix operator*(const Matrix &) const;	
	Matrix operator*(mat_type) const;
	Matrix operator/(const Matrix &) const;	
	Matrix operator/(mat_type) const;
	Matrix operator==(const Matrix &) const;	
	Matrix operator==(mat_type) const;
	Matrix operator!=(const Matrix &) const;	
	Matrix operator!=(mat_type) const;
	Matrix operator>(const Matrix &) const;
	Matrix operator>(mat_type) const;
	Matrix operator<(const Matrix &a) const {return a>*this;}	
	Matrix operator<(mat_type) const;
	Matrix operator>=(const Matrix &) const;	
	Matrix operator>=(mat_type) const;
	Matrix operator<=(const Matrix &a) const {return a>=*this;}	
	Matrix operator<=(mat_type) const;
	Matrix operator&&(const Matrix &) const;	
	Matrix operator&&(mat_type) const;
	Matrix operator||(const Matrix &) const;	
	Matrix operator||(mat_type) const;
	Matrix &operator+=(const Matrix &);
    Matrix &operator-=(const Matrix &);
    Matrix &operator*=(const Matrix &);
    Matrix &operator/=(const Matrix &);
    Matrix &operator+=(mat_type);
    Matrix &operator-=(mat_type);
    Matrix &operator*=(mat_type);
    Matrix &operator/=(mat_type);
    
    Matrix row(int irow) const;
    Matrix col(int icol) const;
    Matrix block(int irow1,int irow2,int icol1,int icol2) const;
    Matrix block_step(int irow1,int irow2,int drow,int icol1,int icol2,int dcol) const;
    Matrix &setrow(int irow,const Matrix &);
    Matrix &setcol(int icol,const Matrix &);
    Matrix &setblock(int irow1,int irow2,int icol1,int icol2,const Matrix &);
    Matrix &setblock_step(int irow1,int irow2,int drow,int icol1,int icol2,int dcol,const Matrix &);
    
    Matrix transpose() const;
    Matrix fliplr() const;
    Matrix flipud() const;
    
    friend matrix ones(int nfil,int ncol);
	friend matrix zeros(int nfil,int ncol);
	friend matrix random_matrix(int nfil,int ncol);
	friend matrix eye(int n);
	friend matrix vector(mat_type x0,mat_type x1,int n);
    friend matrix vector_t(mat_type x0,mat_type x1,int n);
    
	friend Type max<>(const Matrix &);
	friend Type min<>(const Matrix &);
	friend Type sum<>(const Matrix &);
	friend Type mean<>(const Matrix &);
	
	friend Matrix max<>(const Matrix &,const Matrix &);
    friend Matrix max<>(const Matrix &,mat_type);
    friend Matrix max<>(mat_type,const Matrix &);
    friend Matrix min<>(const Matrix &,const Matrix &);
    friend Matrix min<>(const Matrix &,mat_type);
    friend Matrix min<>(mat_type,const Matrix &);
	
	friend int exist<>(const Matrix &);
    friend int isequal<>(const Matrix &,const Matrix &);
    
    friend Matrix cos<>(const Matrix &);
	friend Matrix sin<>(const Matrix &);
	friend Matrix tan<>(const Matrix &);
	friend Matrix acos<>(const Matrix &);
	friend Matrix asin<>(const Matrix &);
	friend Matrix atan<>(const Matrix &);
	friend Matrix cosh<>(const Matrix &);
	friend Matrix sinh<>(const Matrix &);
	friend Matrix tanh<>(const Matrix &);
	friend Matrix exp<>(const Matrix &);
	friend Matrix log<>(const Matrix &);
	friend Matrix log10<>(const Matrix &);
	friend Matrix sqrt<>(const Matrix &);
	friend Matrix abs<>(const Matrix &);
	friend Matrix floor<>(const Matrix &);
	friend Matrix ceil<>(const Matrix &);
	friend Matrix round<>(const Matrix &);
	friend Matrix atan2<>(const Matrix &,const Matrix &);
	friend Matrix atan2<>(double,const Matrix &);
	friend Matrix atan2<>(const Matrix &,double);
	friend Matrix pow<>(const Matrix &,const Matrix &);
	friend Matrix pow<>(double,const Matrix &);
	friend Matrix pow<>(const Matrix &,double);
	friend Matrix pow<>(const Matrix &,int);

    Matrix operator,(const Matrix &) const;
    Matrix solve(Matrix) const;
	Matrix inv() const;
    
    friend class Matrix_block_diag<Type>;
	friend class Matrix<float>;
	friend class Matrix<double>;

};

template <class Type> 
Matrix<Type> operator,(const Matrix<Type> &,const Matrix_block_diag<Type> &);
template <class Type> 
Matrix_block_diag<Type> operator*(const Matrix<Type> &a,const Matrix_block_diag<Type> &D) {return D*a;}
template <class Type> 
Matrix_block_diag<Type> operator*(mat_type n,const Matrix_block_diag<Type> &D){return D*n;}
template <class Type> 
Matrix_block_diag<Type> operator+(const Matrix_block_diag<Type>&a) {return a;}
template <class Type> 
Matrix_block_diag<Type> operator-(const Matrix_block_diag<Type>&);
template <class Type> 
Matrix_block_diag<Type> eye(const Matrix_block_diag<Type> &);

template <class Type>
class Matrix_block_diag {
	Matrix<Type> *m;
	int nb;
public:
	Matrix_block_diag(int nblocks=1);
	~Matrix_block_diag();
	Matrix_block_diag(const Matrix_block_diag &);
	Matrix_block_diag & operator=(const Matrix_block_diag &);
	Matrix_block_diag &set_nblocks(int nblocks);
	operator Matrix<Type>() const;
	Matrix<Type> &block(int i) const;
	int nblocks() const;
	int nrows() const;
	int ncols() const;
	friend Matrix<Type> operator,<>(const Matrix<Type> &,const Matrix_block_diag &);
	Matrix<Type> operator,(const Matrix<Type> &) const;
	Matrix_block_diag operator,(const Matrix_block_diag &) const;
	friend Matrix_block_diag operator*<>(const Matrix<Type> &,const Matrix_block_diag &);
	friend Matrix_block_diag operator*<>(mat_type,const Matrix_block_diag &);
	Matrix_block_diag operator*(const Matrix<Type> &) const;
	Matrix_block_diag operator/(const Matrix<Type> &) const;
	Matrix_block_diag operator*(mat_type) const;
	Matrix_block_diag operator/(mat_type) const;
	friend Matrix_block_diag operator+<>(const Matrix_block_diag&);
	friend Matrix_block_diag operator-<>(const Matrix_block_diag&);
	Matrix_block_diag operator+(const Matrix_block_diag&) const;
	Matrix_block_diag operator-(const Matrix_block_diag&) const;
	Matrix_block_diag operator*(const Matrix_block_diag&) const;
	Matrix<Type> row(int n) const;
	Type operator()(int nfil,int ncol) const;
	friend Matrix_block_diag eye<>(const Matrix_block_diag &);
	Matrix_block_diag transpose() const;

};


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
