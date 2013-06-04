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
	matrix &dim(int nrow,int ncol);
	matrix &redim(int nrow,int ncol);
	matrix &operator=(const matrix &);
	double &operator()(int irow,int icol); 
	double &operator()(int ielem);
	const double &operator()(int irow,int icol) const; 
	const double &operator()(int ielem) const;
	int read(int nrow,int ncol,FILE *fp,char mode='t');
	int write(FILE *fp=stdout,char mode='t') const;
	void write_fmt(const char *fmt,FILE *fp=stdout) const;
	void swap(matrix &);
	void zero(int nrows,int ncols);
	void values(double,...);
	
	friend matrix operator-(double,const matrix &);
	friend matrix operator-(const matrix &);
	friend matrix operator/(double,const matrix &);
	
	matrix operator+(const matrix &) const;
	matrix operator+(double) const;
	matrix operator-(const matrix &) const;
	matrix operator-(double) const;
	matrix operator*(const matrix &) const;	
	matrix operator*(double) const;
	matrix operator/(const matrix &) const;	
	matrix operator/(double) const;
	matrix operator==(const matrix &) const;	
	matrix operator==(double) const;
	matrix operator!=(const matrix &) const;	
	matrix operator!=(double) const;
	matrix operator>(const matrix &) const;
	matrix operator>(double) const;
	matrix operator<(const matrix &a) const {return a>*this;}	
	matrix operator<(double) const;
	matrix operator>=(const matrix &) const;	
	matrix operator>=(double) const;
	matrix operator<=(const matrix &a) const {return a>=*this;}	
	matrix operator<=(double) const;
	matrix operator&&(const matrix &) const;	
	matrix operator&&(double) const;
	matrix operator||(const matrix &) const;	
	matrix operator||(double) const;
	matrix &operator+=(const matrix &);
	matrix &operator-=(const matrix &);
	matrix &operator*=(const matrix &);
	matrix &operator/=(const matrix &);
	matrix &operator+=(double);
	matrix &operator-=(double);
	matrix &operator*=(double);
	matrix &operator/=(double);

	const matrix row(int irow) const;
	const matrix col(int icol) const;
	const matrix block(int irow1,int irow2,int icol1,int icol2) const;
	const matrix block_step(int irow1,int irow2,int drow,int icol1,int icol2,int dcol) const;
	matrix &setrow(int irow,const matrix &);
	matrix &setcol(int icol,const matrix &);
	matrix &setblock(int irow1,int irow2,int icol1,int icol2,const matrix &);
	matrix &setblock_step(int irow1,int irow2,int drow,int icol1,int icol2,int dcol,const matrix &);
   	matrix concatenate(const matrix &a,int dir=0) const;
   	
	matrix transpose() const;
	matrix fliplr() const;
	matrix flipud() const;
	
	friend matrix ones(int nfil,int ncol);
	friend matrix zeros(int nfil,int ncol);
	friend matrix random_matrix(int nfil,int ncol);
	friend matrix eye(int n);
	friend matrix vector(double x0,double x1,int n);
	friend matrix vector_t(double x0,double x1,int n);
	
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
	
	friend int exist(const matrix &);
	friend int isequal(const matrix &,const matrix &);
	
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

	matrix operator,(const matrix &) const;
	matrix solve(matrix) const;
	matrix inv() const;
	
	friend class matrix_block_diag;

};

inline matrix operator+(double n,const matrix &a) {return a+n;}
inline matrix operator+(const matrix &a) {return a;}
matrix operator-(double,const matrix &);
matrix operator-(const matrix &);
inline matrix operator*(double n,const matrix &a) {return a*n;}
matrix operator/(double,const matrix &);
inline matrix operator==(double n,const matrix &a) {return a==n;}
inline matrix operator!=(double n,const matrix &a) {return a!=n;}
inline matrix operator<(double n,const matrix &a) {return a>n;}
inline matrix operator>(double n,const matrix &a) {return a<n;}
inline matrix operator<=(double n,const matrix &a) {return a>=n;}
inline matrix operator>=(double n,const matrix &a) {return a<=n;}
inline matrix operator&&(double n,const matrix &a) {return a&&n;}
inline matrix operator||(double n,const matrix &a) {return a||n;}

matrix ones(int nfil,int ncol);
matrix zeros(int nfil,int ncol);
matrix random_matrix(int nfil,int ncol);
matrix eye(int n);
matrix vector(double x0,double x1,int n);
matrix vector_t(double x0,double x1,int n);

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
	const matrix &block(int i) const;
	matrix &block(int i);
	int nblocks() const;
	int nrows() const;
	int ncols() const;
	friend matrix operator,(const matrix &,const matrix_block_diag &);
	matrix operator,(const matrix &) const;
	matrix_block_diag operator,(const matrix_block_diag &) const;
	friend matrix_block_diag operator*(const matrix &,const matrix_block_diag &);
	friend matrix_block_diag operator*(double,const matrix_block_diag &);
	matrix_block_diag operator*(const matrix &) const;
	matrix_block_diag operator/(const matrix &) const;
	matrix_block_diag operator*(double) const;
	matrix_block_diag operator/(double) const;
	friend matrix_block_diag operator+(const matrix_block_diag&);
	friend matrix_block_diag operator-(const matrix_block_diag&);
	matrix_block_diag operator+(const matrix_block_diag&) const;
	matrix_block_diag operator-(const matrix_block_diag&) const;
	matrix_block_diag operator*(const matrix_block_diag&) const;
	matrix row(int n) const;
	double operator()(int nfil,int ncol) const;
	friend matrix_block_diag eye(const matrix_block_diag &);
	matrix_block_diag transpose() const;

};

matrix operator,(const matrix &,const matrix_block_diag &);
inline matrix_block_diag operator*(const matrix &a,const matrix_block_diag &D) {return D*a;}
inline matrix_block_diag operator*(double n,const matrix_block_diag &D){return D*n;}
inline matrix_block_diag operator+(const matrix_block_diag&a) {return a;}
matrix_block_diag operator-(const matrix_block_diag&);
matrix_block_diag eye(const matrix_block_diag &);

#include<map>
#include<string>

class double_map : public std::map<std::string,double> {
public:	
	inline double &operator[](const std::string& key) {return std::map<std::string,double>::operator[](key);}
	inline const double &operator[](const std::string& key) const {return find(key)->second;};
	double sum();

};

class create_double_map
{
private:
    double_map m_map;
public:
    create_double_map(const std::string& key, const double& val)
    {
        m_map[key] = val;
    }

    create_double_map& operator()(const std::string& key, const double& val)
    {
        m_map[key] = val;
        return *this;
    }

    operator double_map()
    {
        return m_map;
    }
};

class matrix_map_elem : public std::map<std::string,double *> {
public:
	operator double_map();
	matrix_map_elem &operator=(const double_map &);
	matrix_map_elem &operator=(matrix_map_elem &);
	matrix_map_elem &operator=(double);
};


class matrix_map : public std::map<std::string,matrix> {
public:
	inline matrix &operator[](const std::string& key) {return std::map<std::string,matrix>::operator[](key);}
	inline const matrix &operator[](const std::string& key) const {return find(key)->second;};
	matrix_map_elem operator()(int nfil, int ncol);
	const double_map operator()(int nfil, int ncol) const;
	const matrix_map row(int irow) const;
	const matrix_map col(int icol) const;
	const matrix_map block(int irow1,int irow2,int icol1,int icol2) const;
	const matrix_map block_step(int irow1,int irow2,int drow,int icol1,int icol2,int dcol) const;
	matrix_map &setrow(int irow,const matrix_map &);
	matrix_map &setcol(int icol,const matrix_map &);
	matrix_map &setblock(int irow1,int irow2,int icol1,int icol2,const matrix_map &);
	matrix_map &setblock_step(int irow1,int irow2,int drow,int icol1,int icol2,int dcol,const matrix_map &);
	matrix_map &setrow(int irow,const matrix &);
	matrix_map &setcol(int icol,const matrix &);
	matrix_map &setblock(int irow1,int irow2,int icol1,int icol2,const matrix &);
	matrix_map &setblock_step(int irow1,int irow2,int drow,int icol1,int icol2,int dcol,const matrix &);
	
	matrix_map operator*(double d) const;
	matrix_map operator*(const matrix &m) const;
	matrix_map &operator*=(const double_map &d); 
	matrix_map &operator/=(const double_map &d); 
	matrix_map &operator+=(const matrix_map &d);
	matrix_map &operator-=(const matrix_map &d);
	matrix_map &operator*=(const matrix_map &d);
	matrix_map &operator/=(const matrix_map &d);
	
	matrix sum();
	
};

class create_matrix_map
{
private:
    matrix_map m_map;
public:
    create_matrix_map(const std::string& key, const matrix& val)
    {
        m_map[key] = val;
    }

    create_matrix_map& operator()(const std::string& key, const matrix& val)
    {
        m_map[key] = val;
        return *this;
    }

    operator matrix_map()
    {
        return m_map;
    }
};

matrix_map operator*(const double_map &,const matrix &);
inline matrix_map operator*(const matrix &m,const double_map &d) {return d*m;};
inline matrix_map operator*(double d,const matrix_map &a) {return a*d;};
inline matrix_map operator*(const matrix &m,const matrix_map &a) {return a*m;};

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
