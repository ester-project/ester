#ifndef _SYMBOLIC_H
#define _SYMBOLIC_H

#include "matrix.h"
#include "mapping.h"
#include "solver.h"

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <cmath>

class rational {
	int n,d;
	rational &reduce();
	static int gcd(int a,int b);
public:
	rational(int num=1,int den=1);
	int num() const;
	int den() const;
	double eval() const;
	rational operator+(const rational &q) const;
	inline rational operator+() const {return *this;};
	rational operator-() const;
	inline rational operator-(const rational &q) const {return *this+(-q);};
	rational operator*(const rational &q) const;
	rational inv() const;
	inline rational operator/(const rational &q) const {return *this*q.inv();};
	inline rational &operator+=(const rational &q) {return *this=*this+q;};
	inline rational &operator-=(const rational &q) {return *this=*this-q;};
	inline rational &operator*=(const rational &q) {return *this=*this*q;};
	inline rational &operator/=(const rational &q) {return *this=*this/q;};
	rational pow(int) const;
	bool operator==(const rational &q) const;
	inline bool operator!=(const rational &q) const {return !(*this==q);};
	bool operator>(const rational &q) const;
	bool operator<(const rational &q) const;
	inline bool operator>=(const rational &q) const {return !(*this<q);};;
	inline bool operator<=(const rational &q) const {return !(*this>q);};;
	friend std::ostream& operator<<(std::ostream& os, const rational&);
};

inline rational operator+(const int &n,const rational &q) {return q+rational(n);}
inline rational operator-(const int &n,const rational &q) {return -q+rational(n);}
inline rational operator*(const int &n,const rational &q) {return q*rational(n);}
inline rational operator/(const int &n,const rational &q) {return q.inv()*rational(n);}

inline rational operator+(const rational &q,const int &n) {return q+rational(n);}
inline rational operator-(const rational &q,const int &n) {return q-rational(n);}
inline rational operator*(const rational &q,const int &n) {return q*rational(n);}
inline rational operator/(const rational &q,const int &n) {return q*rational(1,n);}

inline double operator+(const double &n,const rational &q) {return n+q.eval();}
inline double operator-(const double &n,const rational &q) {return n-q.eval();}
inline double operator*(const double &n,const rational &q) {return n*q.eval();}
inline double operator/(const double &n,const rational &q) {return n/q.eval();}

inline double operator+(const rational &q,const double &n) {return q.eval()+n;}
inline double operator-(const rational &q,const double &n) {return q.eval()-n;}
inline double operator*(const rational &q,const double &n) {return q.eval()*n;}
inline double operator/(const rational &q,const double &n) {return q.eval()/n;}

rational abs(const rational &q);

std::ostream& operator<<(std::ostream& os, const rational&);


class symbolic;
class sym_vec;
class sym_tens;

class sym {
	class sym_expr;
	class symbol;
	class sym_deriv;
	class sym_num;
	class sym_add;
	class sym_prod;
	class sym_sin;
	class sym_cos;
	class sym_exp;
	class sym_log;

	symbolic *context;
	sym_expr *expr;
	
public:
	sym();
	sym(const sym &);
	explicit sym(const double &);
	explicit sym(const int &);
	sym(const rational &q);
	~sym();
	
	symbolic *check_context() const;
	symbolic *check_context(const sym &s) const;
	symbolic *check_context(const sym_vec &s) const;
	symbolic *check_context(const sym_tens &s) const;
	
	sym &operator=(const sym &);
	sym &operator=(double n) {return *this=sym(n);};
	friend sym operator+(const sym &,const sym &);
	friend sym operator*(const sym &,const sym &);
	sym &operator+=(const sym &s) {return *this=*this+s;}; 
	sym &operator-=(const sym &s) {return *this=*this+sym(-1)*s;};
	sym &operator*=(const sym &s) {return *this=*this*s;};
	sym &operator/=(const sym &s) {return *this=*this*pow(s,-1);};
	sym &operator+=(double n) {return *this+=sym(n);};
	sym &operator-=(double n) {return *this-=sym(n);};
	sym &operator*=(double n) {return *this*=sym(n);};
	sym &operator/=(double n) {return *this/=sym(n);};
	friend sym pow(const sym &,const rational &);
	friend sym sin(const sym &);
	friend sym cos(const sym &);
	friend sym exp(const sym &);
	friend sym log(const sym &);
	friend sym diff(const sym &,const sym &);
	friend sym jacobian(const sym &,const sym &);
	bool operator==(const sym &s) const;
	inline bool operator!=(const sym &s) const {return !(*this==s);};
	inline bool operator==(double n) const {return *this==sym(n);};
	inline bool operator!=(double n) const {return !(*this==sym(n));};

	matrix eval() const;

	friend bool sort_pair_d(const std::pair<sym_expr *,double> &,const std::pair<sym_expr *,double> &);
	friend bool sort_pair_r(const std::pair<sym_expr *,rational> &,const std::pair<sym_expr *,rational> &);
	friend std::ostream& operator<<(std::ostream& os, const sym::sym_expr&);
	friend std::ostream& operator<<(std::ostream& os, const sym&);
	
	friend sym trig_simplify(const sym &);
	
	void add(solver *op,std::string eq_name,std::string var_name) const;
	void add(solver *op,std::string eq_name,std::string var_name,const matrix &d) const;
	void add_ex(solver *op,int n,std::string eq_name,std::string var_name) const;
	void add_ex(solver *op,int n,std::string eq_name,std::string var_name,const matrix &d) const;
	
	void bc_top1_add(solver *op,int n,std::string eq_name,std::string var_name) const;
	void bc_top2_add(solver *op,int n,std::string eq_name,std::string var_name) const;
	void bc_bot1_add(solver *op,int n,std::string eq_name,std::string var_name) const;
	void bc_bot2_add(solver *op,int n,std::string eq_name,std::string var_name) const;
	void bc_top1_add(solver *op,int n,std::string eq_name,std::string var_name,const matrix &d) const;
	void bc_top2_add(solver *op,int n,std::string eq_name,std::string var_name,const matrix &d) const;
	void bc_bot1_add(solver *op,int n,std::string eq_name,std::string var_name,const matrix &d) const;
	void bc_bot2_add(solver *op,int n,std::string eq_name,std::string var_name,const matrix &d) const;
	
	friend class symbolic;
	friend class sym_vec;
	friend class sym_tens;
};

sym operator+(const sym &,const sym &);
inline sym operator+(double n,const sym &s) {return s+sym(n);};
inline sym operator+(const sym &s,double n) {return s+sym(n);};
inline sym operator+(const sym &s) {return s;};
sym operator*(const sym &,const sym &);
inline sym operator*(double n,const sym &s) {return s*sym(n);};
inline sym operator*(const sym &s,double n) {return s*sym(n);};
inline sym operator-(const sym &s) {return (-1)*s;};
template <class T>
	sym operator-(const sym &s1,const T &s2) {return s1+(-s2);};
inline sym operator-(double n,const sym &s2) {return n+(-s2);};
sym pow(const sym &,const rational &);
inline sym pow(const sym &s,const int &n) {return pow(s,rational(n));}
template <class T>
	sym operator/(const sym &s1,const T &s2) {return s1*pow(s2,-1);};
inline sym operator/(double n,const sym &s2) {return n*pow(s2,-1);};
inline sym sqrt(const sym &s) {return pow(s,rational(1,2));};
sym sin(const sym &);
sym cos(const sym &);
sym tan(const sym &);
sym exp(const sym &);
sym log(const sym &);
sym pow(const sym &,const sym &);
inline sym pow(const sym &s,const double &n) {return pow(s,sym(n));}
sym pow(const double &n,const sym &s);
sym diff(const sym &,const sym &);
sym jacobian(const sym &,const sym &);
	
std::ostream& operator<<(std::ostream& os, const sym&);

sym trig_simplify(const sym &);
sym_vec trig_simplify(const sym_vec &);
sym_tens trig_simplify(const sym_tens &);

enum sym_vec_type {COVARIANT,CONTRAVARIANT};

class sym_vec {
	sym s[3];
	sym_vec_type type;
	
	void set_context(symbolic *context);
public:
	sym_vec(sym_vec_type type_set=CONTRAVARIANT);
	
	symbolic *check_context() const;
	symbolic *check_context(const sym &s) const;
	symbolic *check_context(const sym_vec &s) const;
	symbolic *check_context(const sym_tens &s) const;
	
	bool is_covariant();
	bool is_contravariant();
	void set_type(sym_vec_type);
	sym_vec set_variance(const sym_vec &) const;
	sym_vec set_variance(sym_vec_type) const;
	sym &operator()(int);
	sym operator()(int) const;
	sym operator,(const sym_vec &) const;
	sym_vec operator,(const sym_tens &) const;
	sym_vec operator+(const sym_vec &) const;
	friend sym_vec operator-(const sym_vec &v);
	sym_vec operator*(const sym &) const;
	sym_vec operator*(double n) const {return (*this)*sym(n);};
	sym_vec operator/(const sym &) const;
	sym_vec operator/(double n) const {return (*this)/sym(n);};
	friend sym_vec cross(const sym_vec &,const sym_vec &);
	friend sym_tens tensor(const sym_vec &,const sym_vec &);
	sym D(const sym &) const;
	sym_vec D(const sym_vec &) const;	
		
	friend class symbolic;
	friend class sym;
	friend class sym_tens;
};
inline sym_vec operator+(const sym_vec &v) {return v;};
sym_vec operator-(const sym_vec &v);
inline sym_vec operator-(const sym_vec &v,const sym_vec &w) {return v+(-w);};
inline sym_vec operator*(const sym &s,const sym_vec &v) {return v*s;};
inline sym_vec operator*(double n,const sym_vec &v) {return v*n;};
sym_vec cross(const sym_vec &,const sym_vec &);
sym_tens tensor(const sym_vec &,const sym_vec &);

std::ostream& operator<<(std::ostream& os, const sym_vec&);

class sym_tens {
	sym s[3][3];
	sym_vec_type type[2];
	
	void set_context(symbolic *context);
public:
	sym_tens(sym_vec_type type_set0=CONTRAVARIANT,sym_vec_type type_set1=CONTRAVARIANT);
	
	symbolic *check_context() const;
	symbolic *check_context(const sym &s) const;
	symbolic *check_context(const sym_vec &s) const;
	symbolic *check_context(const sym_tens &s) const;
	
	bool is_covariant(int);
	bool is_contravariant(int);
	void set_type(sym_vec_type,sym_vec_type);
	sym_tens set_variance(const sym_tens &) const;
	sym_tens set_variance(sym_vec_type,sym_vec_type) const;
	sym &operator()(int,int);
	sym operator()(int,int) const;
	sym_vec operator,(const sym_vec &) const;
	sym_tens operator,(const sym_tens &) const;
	sym_tens operator+(const sym_tens &) const;
	friend sym_tens operator-(const sym_tens &v);
	sym_tens operator*(const sym &) const;
	sym_tens operator*(double n) const {return (*this)*sym(n);};
	sym_tens operator/(const sym &) const;
	sym_tens operator/(double n) const {return (*this)/sym(n);};
	friend sym_tens tensor(const sym_vec &,const sym_vec &);
	sym operator%(const sym_tens &) const;
	sym_tens T() const;
	
	friend class symbolic;
	friend class sym;
	friend class sym_vec;
};

inline sym_tens operator+(const sym_tens &v) {return v;};
sym_tens operator-(const sym_tens &v);
inline sym_tens operator-(const sym_tens &v,const sym_tens &w) {return v+(-w);};
inline sym_tens operator*(const sym &s,const sym_tens &v) {return v*s;};
inline sym_tens operator*(double n,const sym_tens &v) {return v*n;};

std::ostream& operator<<(std::ostream& os, const sym_tens&);

class symbolic {
	std::map<std::string,sym::symbol> vars;
	std::map<std::string,matrix> val;
	std::map<std::string,int> par;
	mapping map;
	sym one_;
	int maxder;
	inline sym d(const sym &s,int k) {return diff(s,x[k]);}; 
	inline sym G(int i,int j,int k) {return christoffel(i,j,k);};
	
	void add(const sym &s,solver *op,int n,std::string type,
			std::string eq_name,std::string var_name,const matrix &d);
	void add_bc(const sym &s,solver *op,int n,std::string type,
				std::string eq_name,std::string var_name,const matrix &d);

public:
	sym r,zeta,theta,phi,sqrt_g;
	sym rz,rt,rzz,rzt,rtt;
	sym_tens g,g_;
	sym x[3];
	static double tol;
	static bool expand_products;
	static bool trig_simplify;
	static bool axisymmetric;
	static bool spherical;
	static double round_to_tol(double);
	const sym &one;
	symbolic();
	~symbolic() {};
	void init();
	sym regvar(const std::string &name);
	sym regconst(const std::string &name);
	sym regvar_indep(const std::string &name);
	sym var(const std::string &name);
	
	void set_value(const char *name,const matrix &value,int parity=0);
	void set_map(const mapping &map);
	
	sym::sym_expr *derive(const sym::sym_expr &,const sym::symbol &); 
	matrix get_value(const sym::sym_expr &);
	
	sym Dz(const sym &);
	sym Dt(const sym &);
	sym Dphi(const sym &);
	sym det(const sym_tens &);
	sym_tens inv(const sym_tens &);
	
	sym_vec contravariant(const sym_vec &);
	sym_vec covariant(const sym_vec &);
	sym_tens contravariant_contravariant(const sym_tens &);
	sym_tens contravariant_covariant(const sym_tens &);
	sym_tens covariant_contravariant(const sym_tens &);
	sym_tens covariant_covariant(const sym_tens &);

	double perm(int i,int j,int k);
	
	sym_vec gradient(const sym &);
	sym christoffel(int,int,int);
	sym covderiv(const sym_vec &,int,int);
	sym_tens gradient(const sym_vec &);
	sym divergence(const sym_vec &v);
	sym_vec divergence(const sym_tens &t);
	sym laplacian(const sym &s);
	sym_vec curl(const sym_vec &v);
	sym_vec laplacian(const sym_vec &v);
	sym_tens stress(const sym_vec &v);

	friend class sym;
};

inline sym_vec contravariant(const sym_vec &s) {return s.check_context()->contravariant(s);};
inline sym_vec covariant(const sym_vec &s) {return s.check_context()->covariant(s);};
inline sym_tens contravariant_contravariant(const sym_tens &s) {return s.check_context()->contravariant_contravariant(s);};
inline sym_tens contravariant_covariant(const sym_tens &s) {return s.check_context()->contravariant_covariant(s);};
inline sym_tens covariant_contravariant(const sym_tens &s) {return s.check_context()->covariant_contravariant(s);};
inline sym_tens covariant_covariant(const sym_tens &s) {return s.check_context()->covariant_covariant(s);};

inline sym_vec grad(const sym &s) {return s.check_context()->gradient(s);};
inline sym_tens grad(const sym_vec &v) {return v.check_context()->gradient(v);};
inline sym div(const sym_vec &v) {return v.check_context()->divergence(v);};
inline sym_vec div(const sym_tens &t) {return t.check_context()->divergence(t);};
/// \brief Returns the laplacian  of the expression \p s.
inline sym lap(const sym &s) {return s.check_context()->laplacian(s);};
inline sym_vec curl(const sym_vec &v) {return v.check_context()->curl(v);};
inline sym_vec lap(const sym_vec &v) {return v.check_context()->laplacian(v);};
inline sym_tens stress(const sym_vec &v) {return v.check_context()->stress(v);};

inline sym Dz(const sym &s) {return s.check_context()->Dz(s);};
inline sym Dt(const sym &s) {return s.check_context()->Dt(s);};
inline sym Dphi(const sym &s) {return s.check_context()->Dphi(s);};


// Abstract node type
class sym::sym_expr {
public:
	virtual ~sym_expr() {};
	virtual int order() const=0;
	virtual sym_expr *clone() const=0;
	virtual sym_expr *reduce() {return this;};
	virtual sym_expr *derive(const sym_expr &)=0;
	virtual int comp(const sym_expr &s) const=0; // 1 if this>s, 0 if this==s, -1 if this<s
	virtual matrix eval() const=0;
	virtual std::ostream &print(std::ostream&) const=0;
	bool operator>(const sym_expr &s) const {return this->comp(s)==1;};
	bool operator<(const sym_expr &s) const {return this->comp(s)==-1;};
	bool operator==(const sym_expr &s) const {return this->comp(s)==0;};
	bool operator>=(const sym_expr &s) const {return this->comp(s)!=-1;};
	bool operator<=(const sym_expr &s) const {return this->comp(s)!=1;};
	bool operator!=(const sym_expr &s) const {return this->comp(s)!=0;};
	sym_expr *add(const sym_expr &);
	sym_expr *mult(const sym_expr &);
	sym_expr *pow(const rational &);
	sym_expr *sin();
	sym_expr *cos();
	sym_expr *exp();
	sym_expr *log();
};

std::ostream& operator<<(std::ostream& os, const sym::sym_expr&);

// Node types

class sym::sym_num: public sym_expr {
public:
	int order() const {return 0;};
	double value;
	sym_num() {};
	sym_num(const double &);
	sym_num *clone() const;
	int comp(const sym_expr &) const;
	sym_expr *reduce();
	sym_expr *derive(const sym_expr &);
	matrix eval() const;
	std::ostream &print(std::ostream&) const;
};

class sym::symbol: public sym_expr {
public:
	int order() const {return 1;};
	std::string name;
	bool is_const,is_indep;
	symbolic *context;
	symbol *clone() const;
	int comp(const sym_expr &) const;
	sym_expr *derive(const sym_expr &);
	matrix eval() const;
	std::ostream &print(std::ostream&) const;
};

class sym::sym_deriv: public sym_expr {
public:
	sym_deriv() {};
	~sym_deriv();
	sym_deriv(const sym_deriv &);
	int order() const {return 2;};
	symbolic *context;
	sym_expr *oper;
	symbol var;
	sym_deriv *clone() const;
	int comp(const sym_expr &) const;
	sym_expr *reduce();
	sym_expr *derive(const sym_expr &);
	static sym_deriv *create(sym_expr *,const symbol &);
	matrix eval() const;
	std::ostream &print(std::ostream&) const;
};
	
class sym::sym_add: public sym_expr {
public:
	int order() const {return 3;};
	std::vector<std::pair<sym_expr *,double> > oper;
	sym_add() {};
	sym_add(const sym_add &);
	~sym_add();
	sym_add *clone() const;
	int comp(const sym_expr &) const;
	sym_expr *reduce();
	sym_expr *derive(const sym_expr &);
	matrix eval() const;
	std::ostream &print(std::ostream&) const;
	static sym_add *create(sym_expr *,sym_expr *);
	sym_expr *multiply(const sym_add &);
	sym_expr *power(int n);
};

class sym::sym_prod: public sym_expr {
public:
	int order() const {return 4;};
	std::vector<std::pair<sym_expr *,rational> > oper;
	sym_prod() {};
	sym_prod(const sym_prod &);
	~sym_prod();
	sym_prod *clone() const;
	int comp(const sym_expr &) const;
	sym_expr *reduce();
	sym_expr *derive(const sym_expr &);
	matrix eval() const;
	std::ostream &print(std::ostream&) const;
	static sym_prod *create(sym_expr *,sym_expr *);
	static sym_prod *create_pow(sym_expr *,const rational &);
};

bool sort_pair_d(const std::pair<sym::sym_expr *,double> &,const std::pair<sym::sym_expr *,double> &);
bool sort_pair_r(const std::pair<sym::sym_expr *,rational> &,const std::pair<sym::sym_expr *,rational> &);

class sym::sym_sin: public sym_expr {
public:
	int order() const {return 5;};
	sym_expr *oper;
	sym_sin() {};
	~sym_sin();
	sym_sin(const sym_sin &);
	sym_sin *clone() const;
	int comp(const sym_expr &) const;
	sym_expr *reduce();
	sym_expr *derive(const sym_expr &);
	matrix eval() const;
	std::ostream &print(std::ostream&) const;
	static sym_sin *create(sym_expr *);
};

class sym::sym_cos: public sym_expr {
public:
	int order() const {return 6;};
	sym_expr *oper;
	sym_cos() {};
	~sym_cos();
	sym_cos(const sym_cos &);
	sym_cos *clone() const;
	int comp(const sym_expr &) const;
	sym_expr *reduce();
	sym_expr *derive(const sym_expr &);
	matrix eval() const;
	std::ostream &print(std::ostream&) const;
	static sym_cos *create(sym_expr *);
};

class sym::sym_exp: public sym_expr {
public:
	int order() const {return 7;};
	sym_expr *oper;
	sym_exp() {};
	~sym_exp();
	sym_exp(const sym_exp &);
	sym_exp *clone() const;
	int comp(const sym_expr &) const;
	sym_expr *reduce();
	sym_expr *derive(const sym_expr &);
	matrix eval() const;
	std::ostream &print(std::ostream&) const;
	static sym_exp *create(sym_expr *);
};

class sym::sym_log: public sym_expr {
public:
	int order() const {return 8;};
	sym_expr *oper;
	sym_log() {};
	~sym_log();
	sym_log(const sym_log &);
	sym_log *clone() const;
	int comp(const sym_expr &) const;
	sym_expr *reduce();
	sym_expr *derive(const sym_expr &);
	matrix eval() const;
	std::ostream &print(std::ostream&) const;
	static sym_log *create(sym_expr *);
};

#endif
