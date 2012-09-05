#ifndef _SOLVER_H
#define _SOLVER_H

#include"matrix.h"

class solver_operator {
public:
	int verbose;
	virtual ~solver_operator() {};
	virtual void set_block(int iblock,const matrix &) {};
	virtual void set_blocksup(int iblock,const matrix &) {};
	virtual void set_blockinf(int iblock,const matrix &) {};
	virtual void fwd_subs(matrix &) {};
	virtual void back_subs(matrix &) {};
	virtual matrix solve(const matrix &) {matrix a;return a;};
};

class solver_full: public solver_operator {
	matrix *m,*msup,*minf,*c,*r;	
	matrix blk_index,ipiv_flag;
	int nb,lu_flag,*N,oc;
	char tempdir[50];
	int **ipiv;
	void lu_calc();
	void lu_block(int i);
	void solve_block(int i,char trans,matrix &x);
	void read_block(int i);
	void write_block(int i,const matrix &);
	void flush_block(int i);
	void read_blocksup(int i);
	void write_blocksup(int i,const matrix &);
	void flush_blocksup(int i);
	void read_blockinf(int i);
	void write_blockinf(int i,const matrix &);
	void flush_blockinf(int i);
public:
	solver_full(int nblocks,int offcore=0);
	~solver_full();
	void set_block(int iblock,const matrix &);
	void set_blocksup(int iblock,const matrix &);
	void set_blockinf(int iblock,const matrix &);
	void fwd_subs(matrix &);
	void back_subs(matrix &);
	matrix solve(const matrix &);
};

class solver_iter: public solver_operator {
	
public:
	solver_iter() {};
	~solver_iter() {};
	void set_block(int iblock,const matrix &) {};
	void set_blocksup(int iblock,const matrix &) {};
	void set_blockinf(int iblock,const matrix &) {};
	void fwd_subs(matrix &) {};
	void back_subs(matrix &) {};
	matrix solve(const matrix &a) {return 0*a;};
};

/* 	A solver_elem object stores one term of a equation in the form D*(L,I*x,R).
	The type indicate which of the matrices D,L,I,R are different from unity:
	
		'd': D
		'l': D,L
		'r': D,R
		'f': D,L,R
		'm': D,L,I
		's': D,R,I
		'g': D,L,R,I
	
   	The set of terms in one given eq. referred to one given variable is represented
   	as a linked list of solver_elem objects. */

class solver_elem {
	char type;
	matrix D,L,R,I;
	solver_elem *next;
	friend class solver;
	friend class solver_block;
};

/*	A solver_block objects contains all the equations for a given block. A block
	refers to one domain or to one set of boundary conditions.
	
	The equations are represented as an array of linked lists of solver_elem objects
	such that eq[i][j] is a pointer to the first solver_elem object of the list
	representing the terms in the i-th equation associated with the j-th variable.
	
	nv is the number of variables (equations) of the problem.
	eq_set(i) (i<nv) is 1 if the i-th equation is defined in ths block.
	nr[i],nth[i] is the size of the i-th equation (It is automatically taken from
					the size of D).	
	
	Methods:
		
		init(nvar) is called first to initialize the object for nvar variables
		destroy() should be called when the object is no longer needed. After calling
					destroy(), the object can be re-initialized with a different number
					of variables
		reset() to clean out all the equations
		reset(i) to clean out only the i-th equation
		add_*(...) to add a term in the equation
*/

class solver_block {
	int nv,*nr,*nth;
	solver_elem ***eq;
	matrix eq_set;
	friend class solver;
public:
	solver_block() {};
	~solver_block() {};
	void init(int nvar);
	void destroy();
	void reset();
	void reset(int ieq);
	void add_d(int ieq,int ivar,const matrix &d);
	void add_l(int ieq,int ivar,const matrix &d,const matrix &l);
	void add_r(int ieq,int ivar,const matrix &d,const matrix &r);
	void add_lr(int ieq,int ivar,const matrix &d,const matrix &l,const matrix &r);
	void add_li(int ieq,int ivar,const matrix &d,const matrix &l,const matrix &i);
	void add_ri(int ieq,int ivar,const matrix &d,const matrix &r,const matrix &i);
	void add_lri(int ieq,int ivar,const matrix &d,const matrix &l,const matrix &r,const matrix &i);
	void add(int ieq,int ivar,char type, const matrix *d, const matrix *l, const matrix *r, const matrix *i);
};

class solver {
	int nb,nv;
	int **var_nr,**var_ntop,**var_nbot,**var_nth,solver_N,*def_nr;
	char type[21],**var;
	int initd,sync;
	solver_operator *op;
	solver_block *block;
	solver_block *bc_bot1;
	solver_block *bc_bot2;
	solver_block *bc_top1;
	solver_block *bc_top2;
	solver_block *bc_pol;
	solver_block *bc_eq;
	matrix *rhs,*x;
public:
	int use_cgs,maxit_ref,maxit_cgs,verbose;
	double rel_tol,abs_tol;
	void (*mult_fcn)(matrix *,void *context);
	void *mult_context;
	solver();
	void destroy();
	~solver() {if(initd) destroy();};
	void init(int nblock,int nvar,const char *solver_type);
	void set_nr(int *nr);
	void reset();
	void reset(int iblock);
	void reset(int iblock,int ieq);
	void regvar(const char *var_name);
	int get_nvar();
	int get_nblocks();
	int get_id(const char *varn);
	void set_rhs(const char *eqn,const matrix &b);
	matrix get_rhs(const char *eqn);
	matrix get_var(const char *varn);
	matrix get_rhs(int ieq);
	matrix get_var(int ivar);
	void solve(int *info=NULL); 
	// info[0]: LU factorization
	// info[1]: Iterative refinement
	// info[2]: CGS iteration
	// info[3]: Number of iterative refinement iterations (-1 if no convergence, -2 if error)
	// info[4]: Number of CGS iterations (-1 if no convergence, -2 if error)
	void mult(matrix *);
	void mult_op(matrix *);
	int cgs(const matrix &rhs,matrix &x,int maxit);
	void create();
	void create_full();
	void wrap(const matrix *,matrix *);
	void unwrap(matrix *,const matrix *);
	void check();
	void calc_struct();
	void check_block(int iblock,matrix &m);
	void check_block_sup(int iblock,matrix &m);
	void check_block_inf(int iblock,matrix &m);
	
	void add_d(const char *eqn, const char *varn,const matrix &d);
	void add_l(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l);
	void add_r(const char *eqn, const char *varn,const matrix &d,const matrix &r);
	void add_lr(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &r);
	void add_li(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &i);
	void add_ri(const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i);
	void add_lri(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &r,const matrix &i);
	void add_d(int iblock,const char *eqn, const char *varn,const matrix &d);
	void add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l);
	void add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r);
	void add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r);
	void add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i);
	void add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i);
	void add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i);
	void bc_pol_add_d(const char *eqn, const char *varn,const matrix &d);
	void bc_pol_add_l(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l);
	void bc_pol_add_r(const char *eqn, const char *varn,const matrix &d,const matrix &r);
	void bc_pol_add_lr(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &r);
	void bc_pol_add_li(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &i);
	void bc_pol_add_ri(const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i);
	void bc_pol_add_lri(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &r,const matrix &i);
	void bc_pol_add_d(int iblock,const char *eqn, const char *varn,const matrix &d);
	void bc_pol_add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l);
	void bc_pol_add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r);
	void bc_pol_add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r);
	void bc_pol_add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i);
	void bc_pol_add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i);
	void bc_pol_add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i);
	void bc_eq_add_d(const char *eqn, const char *varn,const matrix &d);
	void bc_eq_add_l(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l);
	void bc_eq_add_r(const char *eqn, const char *varn,const matrix &d,const matrix &r);
	void bc_eq_add_lr(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &r);
	void bc_eq_add_li(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &i);
	void bc_eq_add_ri(const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i);
	void bc_eq_add_lri(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &r,const matrix &i);
	void bc_eq_add_d(int iblock,const char *eqn, const char *varn,const matrix &d);
	void bc_eq_add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l);
	void bc_eq_add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r);
	void bc_eq_add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r);
	void bc_eq_add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i);
	void bc_eq_add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i);
	void bc_eq_add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i);
	void bc_bot1_add_d(int iblock,const char *eqn, const char *varn,const matrix &d);
	void bc_bot1_add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l);
	void bc_bot1_add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r);
	void bc_bot1_add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r);
	void bc_bot1_add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i);
	void bc_bot1_add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i);
	void bc_bot1_add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i);
	void bc_bot2_add_d(int iblock,const char *eqn, const char *varn,const matrix &d);
	void bc_bot2_add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l);
	void bc_bot2_add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r);
	void bc_bot2_add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r);
	void bc_bot2_add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i);
	void bc_bot2_add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i);
	void bc_bot2_add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i);
	void bc_top1_add_d(int iblock,const char *eqn, const char *varn,const matrix &d);
	void bc_top1_add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l);
	void bc_top1_add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r);
	void bc_top1_add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r);
	void bc_top1_add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i);
	void bc_top1_add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i);
	void bc_top1_add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i);
	void bc_top2_add_d(int iblock,const char *eqn, const char *varn,const matrix &d);
	void bc_top2_add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l);
	void bc_top2_add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r);
	void bc_top2_add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r);
	void bc_top2_add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i);
	void bc_top2_add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i);
	void bc_top2_add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i);
	
};


#endif
