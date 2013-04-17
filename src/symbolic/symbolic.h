#ifndef _SYMBOLIC_H
#define _SYMBOLIC_H

#include"matrix.h"
#include"mapping.h"
#include"solver.h"
#include<stdio.h>

class symbolic;
class symterm;
class sym_vec;
class sym_tens;

enum sym_vec_type {COVARIANT,CONTRAVARIANT};

class sym {
	symterm *terms;
	const symbolic *context;
	symterm *addterm();
	symterm *addterm(const symterm &);
	void destroy();
	void copy(const sym &);
	sym &simplify();
	void simplify0();
	bool trig_simplify();
	void trig_simplify2();
	void add_bc(solver *op,int n,const char *type,const char *eq_name,const char *var_name,const matrix &d) const;
public:
	friend class symbolic;
	friend class sym_vec;
	friend class sym_tens;
	sym();
	~sym();
	sym(const sym &);
	void dump() const;
	const symbolic *check_context() const;
	const symbolic *check_context(const sym &s) const;
	const symbolic *check_context(const sym_vec &s) const;
	const symbolic *check_context(const sym_tens &s) const;
	sym &operator=(const sym &);
	sym operator+(const sym &) const;
	sym operator+(const double) const;
	sym operator-(const sym &) const;
	sym operator-(const double) const;
	friend sym operator-(const sym &);
	sym operator*(const sym &) const;
	sym operator*(const double) const;
	sym operator/(const sym &) const;
	sym operator/(const double) const;
	friend sym operator/(const double,const sym &);
	void write(FILE *fp=stdout) const;
	matrix eval() const;
	sym jacobian(const char *name,int dz,int dt) const;
	sym jacobian_r(int dz,int dt) const;
	void add(solver *op,const char *eq_name,const char *var_name) const;
	void add(solver *op,const char *eq_name,const char *var_name,const matrix &d) const;
	inline void bc_top1_add(solver *op,int n,const char *eq_name,const char *var_name) const {add_bc(op,n,"top1",eq_name,var_name,ones(1,1));};
	inline void bc_top2_add(solver *op,int n,const char *eq_name,const char *var_name) const {add_bc(op,n,"top2",eq_name,var_name,ones(1,1));};
	inline void bc_bot1_add(solver *op,int n,const char *eq_name,const char *var_name) const {add_bc(op,n,"bot1",eq_name,var_name,ones(1,1));};
	inline void bc_bot2_add(solver *op,int n,const char *eq_name,const char *var_name) const {add_bc(op,n,"bot2",eq_name,var_name,ones(1,1));};
	inline void bc_top1_add(solver *op,int n,const char *eq_name,const char *var_name,const matrix &d) const {add_bc(op,n,"top1",eq_name,var_name,d);};
	inline void bc_top2_add(solver *op,int n,const char *eq_name,const char *var_name,const matrix &d) const {add_bc(op,n,"top2",eq_name,var_name,d);};
	inline void bc_bot1_add(solver *op,int n,const char *eq_name,const char *var_name,const matrix &d) const {add_bc(op,n,"bot1",eq_name,var_name,d);};
	inline void bc_bot2_add(solver *op,int n,const char *eq_name,const char *var_name,const matrix &d) const {add_bc(op,n,"bot2",eq_name,var_name,d);};
};

sym operator+(const double,const sym &);
sym operator+(const sym &);
sym operator-(const double,const sym &);
sym operator-(const sym &);
sym operator*(const double,const sym &);
sym operator/(const double,const sym &);


class sym_vec {
	sym s[3];
	sym_vec_type type;
	void set_context(const symbolic *context);
public:
	sym_vec(sym_vec_type type_set=CONTRAVARIANT);
	const symbolic *check_context() const;
	const symbolic *check_context(const sym &s) const;
	const symbolic *check_context(const sym_vec &s) const;
	const symbolic *check_context(const sym_tens &s) const;
	sym_vec set_variance(const sym_vec &) const;
	sym_vec set_variance(sym_vec_type) const;
	sym &operator()(int);
	sym operator()(int) const;
	sym operator,(const sym_vec &) const;
	sym_vec operator,(const sym_tens &) const;
	sym_vec operator+(const sym_vec &) const;
	sym_vec operator-(const sym_vec &) const;
	friend sym_vec operator-(const sym_vec &v);
	sym_vec operator*(const sym &) const;
	sym_vec operator*(double) const;
	sym_vec operator/(const sym &) const;
	sym_vec operator/(double) const;
	friend sym_vec cross(const sym_vec &,const sym_vec &);
	friend sym_tens tensor(const sym_vec &,const sym_vec &);
	sym D(const sym &) const;
	sym_vec D(const sym_vec &) const;
	bool is_covariant();
	bool is_contravariant();
	void set_type(sym_vec_type);
	friend class symbolic;
	friend class sym_tens;
};
sym_vec operator+(const sym_vec &v);
sym_vec operator-(const sym_vec &v);
sym_vec operator*(const sym &,const sym_vec &);
sym_vec operator*(double,const sym_vec &);
sym_vec cross(const sym_vec &,const sym_vec &);
sym_tens tensor(const sym_vec &,const sym_vec &);

class sym_tens {
	sym s[3][3];
	sym_vec_type type[2];
	void set_context(const symbolic *context);
public:
	sym_tens(sym_vec_type type_set0=CONTRAVARIANT,sym_vec_type type_set1=CONTRAVARIANT);
	const symbolic *check_context() const;
	const symbolic *check_context(const sym &s) const;
	const symbolic *check_context(const sym_vec &s) const;
	const symbolic *check_context(const sym_tens &s) const;
	sym_tens set_variance(const sym_tens &) const;
	sym_tens set_variance(sym_vec_type,sym_vec_type) const;
	sym &operator()(int,int);
	sym operator()(int,int) const;
	sym_tens operator+(const sym_tens &) const;
	sym_tens operator-(const sym_tens &) const;
	friend sym_tens operator-(const sym_tens &v);
	friend sym_tens tensor(const sym_vec &,const sym_vec &);
	sym_tens operator*(const sym &) const;
	sym_tens operator*(double) const;
	sym_tens operator/(const sym &) const;
	sym_tens operator/(double) const;
	sym_vec operator,(const sym_vec &) const;
	sym_tens operator,(const sym_tens &) const;
	sym operator%(const sym_tens &) const;
	bool is_covariant(int);
	bool is_contravariant(int);
	void set_type(sym_vec_type,sym_vec_type);
	friend class symbolic;
	friend class sym_vec;
};

sym_tens operator+(const sym_tens &v);
sym_tens operator-(const sym_tens &v);
sym_tens operator*(const sym &,const sym_tens &);
sym_tens operator*(double,const sym_tens &);

class symterm {
		int nvar,maxder;
		double num;
		int **r;
		int cost,sint;
		int ***var;
		symterm *next;
		symterm(const symbolic *context);
		symterm(const symterm &);
		~symterm();
		bool operator==(const symterm &) const;
		symterm operator*(const symterm &) const;
		symterm operator/(const symterm &) const;
		void dump() const;
		friend sym operator-(const sym &);
		friend class sym;
		friend class symbolic;
};

class symbolic {
	int nvar,maxder;	
	sym _r,_sint,_cost,_one,_sqrt_g;
	sym_tens _g,_g_;
	char **var_name;
	matrix *var_value;
	int *var_par;
	bool *is_const;
	mapping map;
	int id(const char *name) const;
	void write_term(const symterm &t,FILE *) const;
	void format_deriv(char *str,int dz,int dt,int par) const;
	void write_var(FILE *fp,const char *str,int n) const;
	matrix eval(const symterm &t) const;
	void eval_deriv(matrix &,int dz,int dt,int par) const;
	matrix eval_var(const matrix &,int n) const;
	sym d(const sym &,int) const; 
	inline sym G(int i,int j,int k) const {return christoffel(i,j,k);};
public:
	sym_tens &g,&g_;
	sym &sqrt_g;
	double tol;
	const sym &r,&sint,&cost,&one;
	symbolic(int n_var,int max_der);
	~symbolic();
	sym regvar(const char *name);
	sym regconst(const char *name);
	void set_value(const char *name,const matrix &value,int parity=0);
	void set_map(const mapping &map);
	sym var(const char *name) const;
	sym Dz(const sym &) const;
	sym Dt(const sym &) const;
	sym_vec contravariant(const sym_vec &) const;
	sym_vec covariant(const sym_vec &) const;
	sym_tens contravariant_contravariant(const sym_tens &) const;
	sym_tens contravariant_covariant(const sym_tens &) const;
	sym_tens covariant_contravariant(const sym_tens &) const;
	sym_tens covariant_covariant(const sym_tens &) const;
	sym_vec gradient(const sym &) const;
	sym christoffel(int,int,int) const;
	sym covderiv(const sym_vec &,int,int) const;
	sym divergence(const sym_vec &v) const;
	sym_vec divergence(const sym_tens &t) const;
	sym laplacian(const sym &s) const;
	double perm(int i,int j,int k) const;
	sym_vec curl(const sym_vec &v) const;
	sym_vec laplacian(const sym_vec &v) const;
	sym spherical(const sym &s) const;
	sym_vec spherical(const sym_vec &s) const;
	sym_tens spherical(const sym_tens &s) const;
	sym_tens stress(const sym_vec &v) const;
	friend class sym;
	friend class sym_vec;
	friend class sym_tens;
	friend class symterm;

};

inline sym Dz(const sym &s) {return s.check_context()->Dz(s);};
inline sym Dt(const sym &s) {return s.check_context()->Dt(s);};
sym DzDt(const sym &);
sym Dz(const sym &,int);
sym Dt(const sym &,int);
sym DzDt(const sym &,int,int);
inline sym_vec contravariant(const sym_vec &s) {return s.check_context()->contravariant(s);};
inline sym_vec covariant(const sym_vec &s) {return s.check_context()->covariant(s);};
inline sym_tens contravariant_contravariant(const sym_tens &s) {return s.check_context()->contravariant_contravariant(s);};
inline sym_tens contravariant_covariant(const sym_tens &s) {return s.check_context()->contravariant_covariant(s);};
inline sym_tens covariant_contravariant(const sym_tens &s) {return s.check_context()->covariant_contravariant(s);};
inline sym_tens covariant_covariant(const sym_tens &s) {return s.check_context()->covariant_covariant(s);};
inline sym_vec grad(const sym &s) {return s.check_context()->gradient(s);};
inline sym div(const sym_vec &v) {return v.check_context()->divergence(v);};
inline sym_vec div(const sym_tens &t) {return t.check_context()->divergence(t);};
inline sym lap(const sym &s) {return s.check_context()->laplacian(s);};
inline sym_vec curl(const sym_vec &v) {return v.check_context()->curl(v);};
inline sym_vec lap(const sym_vec &v) {return v.check_context()->laplacian(v);};
inline sym spherical(const sym &s) {return s.check_context()->spherical(s);};
inline sym_vec spherical(const sym_vec &v) {return v.check_context()->spherical(v);};
inline sym_tens spherical(const sym_tens &t) {return t.check_context()->spherical(t);};
inline sym_tens stress(const sym_vec &v) {return v.check_context()->stress(v);};



#endif
