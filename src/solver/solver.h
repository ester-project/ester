#ifndef _SOLVER_H
#define _SOLVER_H

#include"matrix.h"

class solver_operator {
public:
	int verbose;
	virtual ~solver_operator() {};
	virtual void left_precond(matrix &)=0;
	virtual void right_precond(matrix &)=0;
	virtual matrix solve(const matrix &)=0;
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
	void fwd_subs(matrix &);
	void back_subs(matrix &);
public:
	solver_full(int nblocks,int offcore=0);
	~solver_full();
	void set_block(int iblock,matrix &);
	void set_blocksup(int iblock,matrix &);
	void set_blockinf(int iblock,matrix &);
	void left_precond(matrix &x) {fwd_subs(x);}
	void right_precond(matrix &x) {back_subs(x);}
	matrix solve(const matrix &);
};

class solver_iter: public solver_operator {
	
public:
	solver_iter() {};
	~solver_iter() {};
	void left_precond(matrix &) {};
	void right_precond(matrix &) {};
	matrix solve(const matrix &a) {return 0*a;};
};
	
class solver {

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
		friend class solver_block;
		friend class solver;
	};
	
	/*	A solver_block object contains all the equations for a given block. A block
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
		void add(int ieq,int ivar,char type, const matrix *d, const matrix *l, const matrix *r, const matrix *i);
		
		inline void add_d(int ieq,int ivar,const matrix &d) {
			add(ieq,ivar,'d',&d,NULL,NULL,NULL);
		}
		inline void add_l(int ieq,int ivar,const matrix &d,const matrix &l) {
			add(ieq,ivar,'l',&d,&l,NULL,NULL);
		}
		inline void add_r(int ieq,int ivar,const matrix &d,const matrix &r) {
			add(ieq,ivar,'r',&d,NULL,&r,NULL);
		}
		inline void add_lr(int ieq,int ivar,const matrix &d,const matrix &l,const matrix &r) {
			add(ieq,ivar,'f',&d,&l,&r,NULL);
		}
		inline void add_li(int ieq,int ivar,const matrix &d,const matrix &l,const matrix &i) {
			add(ieq,ivar,'m',&d,&l,NULL,&i);
		}
		inline void add_ri(int ieq,int ivar,const matrix &d,const matrix &r,const matrix &i) {
			add(ieq,ivar,'s',&d,NULL,&r,&i);
		}
		inline void add_lri(int ieq,int ivar,const matrix &d,const matrix &l,const matrix &r,const matrix &i) {
			add(ieq,ivar,'g',&d,&l,&r,&i);
		}
	
	};

	int nb,nv;
	int **var_nr,**var_ntop,**var_nbot,**var_nth,solver_N,*def_nr;
	char type[21],**var;
	int initd,sync;
	matrix dep,reg;
	solver_operator *op;
	solver_block *block;
	solver_block *bc_bot1;
	solver_block *bc_bot2;
	solver_block *bc_top1;
	solver_block *bc_top2;
	solver_block *bc_pol;
	solver_block *bc_eq;
	matrix *rhs,*sol;
public:
	int use_cgs,maxit_ref,maxit_cgs,verbose;
	int debug;
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
	void reset(int iblock,const char *eq_name);
	void reset(const char *eq_name);
	void regvar(const char *var_name,int dependent=0);
	inline void regvar_dep(const char *var_name) {regvar(var_name,1);}
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
	void fill_void_blocks();
	int check_struct();
	void check_struct_error(const char *err_msg,int n,int i,int j,solver_elem *p);
	int check_struct_block(int n,int i,int j);
	int check_struct_bc_th(int n,int i,int j,const char *bctype);
	int check_struct_bc(int n,int i,int j,const char *bctype);
	void check_full(int n, const matrix &opi,int pos);
	void subst_dep();
	void subst_dep_eq(const char *block_type,solver_block *,int n,int i,int j);
	void subst_dep_elem(int i,int k,solver_block *bb,solver_elem *p,const matrix &d,int n2,int m2);
	void solve_dep();
	
	void add(int iblock,const char *eqn, const char *varn,const char *block_type,char type,const matrix *d,const matrix *l,const matrix *r,const matrix *i);
	void add(const char *eqn, const char *varn,const char *block_type,char type,const matrix *d,const matrix_block_diag *l,const matrix *r,const matrix *i);
	
	inline void add_d(const char *eqn, const char *varn,const matrix &d) {
		add(eqn,varn,"block",'d',&d,NULL,NULL,NULL);
	}
	inline void add_l(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l) {
		add(eqn,varn,"block",'l',&d,&l,NULL,NULL);
	}
	inline void add_r(const char *eqn, const char *varn,const matrix &d,const matrix &r) {
		add(eqn,varn,"block",'r',&d,NULL,&r,NULL);
	}
	inline void add_lr(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &r) {
		add(eqn,varn,"block",'f',&d,&l,&r,NULL);
	}
	inline void add_li(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &i) {
		add(eqn,varn,"block",'m',&d,&l,NULL,&i);
	}
	inline void add_ri(const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i) {
		add(eqn,varn,"block",'s',&d,NULL,&r,&i);
	}
	inline void add_lri(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &r,const matrix &i) {
		add(eqn,varn,"block",'g',&d,&l,&r,&i);
	}
	
	inline void bc_pol_add_d(const char *eqn, const char *varn,const matrix &d) {
		add(eqn,varn,"bc_pol",'d',&d,NULL,NULL,NULL);
	}
	inline void bc_pol_add_l(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l) {
		add(eqn,varn,"bc_pol",'l',&d,&l,NULL,NULL);
	}
	inline void bc_pol_add_r(const char *eqn, const char *varn,const matrix &d,const matrix &r) {
		add(eqn,varn,"bc_pol",'r',&d,NULL,&r,NULL);
	}
	inline void bc_pol_add_lr(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &r) {
		add(eqn,varn,"bc_pol",'f',&d,&l,&r,NULL);
	}
	inline void bc_pol_add_li(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &i) {
		add(eqn,varn,"bc_pol",'m',&d,&l,NULL,&i);
	}
	inline void bc_pol_add_ri(const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i) {
		add(eqn,varn,"bc_pol",'s',&d,NULL,&r,&i);
	}
	inline void bc_pol_add_lri(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &r,const matrix &i) {
		add(eqn,varn,"bc_pol",'g',&d,&l,&r,&i);
	}
	
	inline void bc_eq_add_d(const char *eqn, const char *varn,const matrix &d) {
		add(eqn,varn,"bc_eq",'d',&d,NULL,NULL,NULL);
	}
	inline void bc_eq_add_l(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l) {
		add(eqn,varn,"bc_eq",'l',&d,&l,NULL,NULL);
	}
	inline void bc_eq_add_r(const char *eqn, const char *varn,const matrix &d,const matrix &r) {
		add(eqn,varn,"bc_eq",'r',&d,NULL,&r,NULL);
	}
	inline void bc_eq_add_lr(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &r) {
		add(eqn,varn,"bc_eq",'f',&d,&l,&r,NULL);
	}
	inline void bc_eq_add_li(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &i) {
		add(eqn,varn,"bc_eq",'m',&d,&l,NULL,&i);
	}
	inline void bc_eq_add_ri(const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i) {
		add(eqn,varn,"bc_eq",'s',&d,NULL,&r,&i);
	}
	inline void bc_eq_add_lri(const char *eqn, const char *varn,const matrix &d,const matrix_block_diag &l,const matrix &r,const matrix &i) {
		add(eqn,varn,"bc_eq",'g',&d,&l,&r,&i);
	}
	
	inline void add_d(int iblock,const char *eqn, const char *varn,const matrix &d) {
		add(iblock,eqn,varn,"block",'d',&d,NULL,NULL,NULL);
	}
	inline void add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l) {
		add(iblock,eqn,varn,"block",'l',&d,&l,NULL,NULL);
	}
	inline void add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r) {
		add(iblock,eqn,varn,"block",'r',&d,NULL,&r,NULL);
	}
	inline void add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r) {
		add(iblock,eqn,varn,"block",'f',&d,&l,&r,NULL);
	}
	inline void add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i) {
		add(iblock,eqn,varn,"block",'m',&d,&l,NULL,&i);
	}
	inline void add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i) {
		add(iblock,eqn,varn,"block",'s',&d,NULL,&r,&i);
	}
	inline void add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i) {
		add(iblock,eqn,varn,"block",'g',&d,&l,&r,&i);
	}
	
	inline void bc_pol_add_d(int iblock,const char *eqn, const char *varn,const matrix &d) {
		add(iblock,eqn,varn,"bc_pol",'d',&d,NULL,NULL,NULL);
	}
	inline void bc_pol_add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l) {
		add(iblock,eqn,varn,"bc_pol",'l',&d,&l,NULL,NULL);
	}
	inline void bc_pol_add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r) {
		add(iblock,eqn,varn,"bc_pol",'r',&d,NULL,&r,NULL);
	}
	inline void bc_pol_add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r) {
		add(iblock,eqn,varn,"bc_pol",'f',&d,&l,&r,NULL);
	}
	inline void bc_pol_add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i) {
		add(iblock,eqn,varn,"bc_pol",'m',&d,&l,NULL,&i);
	}
	inline void bc_pol_add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i) {
		add(iblock,eqn,varn,"bc_pol",'s',&d,NULL,&r,&i);
	}
	inline void bc_pol_add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i) {
		add(iblock,eqn,varn,"bc_pol",'g',&d,&l,&r,&i);
	}
	
	inline void bc_eq_add_d(int iblock,const char *eqn, const char *varn,const matrix &d) {
		add(iblock,eqn,varn,"bc_eq",'d',&d,NULL,NULL,NULL);
	}
	inline void bc_eq_add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l) {
		add(iblock,eqn,varn,"bc_eq",'l',&d,&l,NULL,NULL);
	}
	inline void bc_eq_add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r) {
		add(iblock,eqn,varn,"bc_eq",'r',&d,NULL,&r,NULL);
	}
	inline void bc_eq_add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r) {
		add(iblock,eqn,varn,"bc_eq",'f',&d,&l,&r,NULL);
	}
	inline void bc_eq_add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i) {
		add(iblock,eqn,varn,"bc_eq",'m',&d,&l,NULL,&i);
	}
	inline void bc_eq_add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i) {
		add(iblock,eqn,varn,"bc_eq",'s',&d,NULL,&r,&i);
	}
	inline void bc_eq_add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i) {
		add(iblock,eqn,varn,"bc_eq",'g',&d,&l,&r,&i);
	}
	
	inline void bc_bot1_add_d(int iblock,const char *eqn, const char *varn,const matrix &d) {
		add(iblock,eqn,varn,"bc_bot1",'d',&d,NULL,NULL,NULL);
	}
	inline void bc_bot1_add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l) {
		add(iblock,eqn,varn,"bc_bot1",'l',&d,&l,NULL,NULL);
	}
	inline void bc_bot1_add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r) {
		add(iblock,eqn,varn,"bc_bot1",'r',&d,NULL,&r,NULL);
	}
	inline void bc_bot1_add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r) {
		add(iblock,eqn,varn,"bc_bot1",'f',&d,&l,&r,NULL);
	}
	inline void bc_bot1_add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i) {
		add(iblock,eqn,varn,"bc_bot1",'m',&d,&l,NULL,&i);
	}
	inline void bc_bot1_add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i) {
		add(iblock,eqn,varn,"bc_bot1",'s',&d,NULL,&r,&i);
	}
	inline void bc_bot1_add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i) {
		add(iblock,eqn,varn,"bc_bot1",'g',&d,&l,&r,&i);
	}
	
	inline void bc_bot2_add_d(int iblock,const char *eqn, const char *varn,const matrix &d) {
		add(iblock,eqn,varn,"bc_bot2",'d',&d,NULL,NULL,NULL);
	}
	inline void bc_bot2_add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l) {
		add(iblock,eqn,varn,"bc_bot2",'l',&d,&l,NULL,NULL);
	}
	inline void bc_bot2_add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r) {
		add(iblock,eqn,varn,"bc_bot2",'r',&d,NULL,&r,NULL);
	}
	inline void bc_bot2_add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r) {
		add(iblock,eqn,varn,"bc_bot2",'f',&d,&l,&r,NULL);
	}
	inline void bc_bot2_add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i) {
		add(iblock,eqn,varn,"bc_bot2",'m',&d,&l,NULL,&i);
	}
	inline void bc_bot2_add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i) {
		add(iblock,eqn,varn,"bc_bot2",'s',&d,NULL,&r,&i);
	}
	inline void bc_bot2_add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i) {
		add(iblock,eqn,varn,"bc_bot2",'g',&d,&l,&r,&i);
	}
	
	inline void bc_top1_add_d(int iblock,const char *eqn, const char *varn,const matrix &d) {
		add(iblock,eqn,varn,"bc_top1",'d',&d,NULL,NULL,NULL);
	}
	inline void bc_top1_add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l) {
		add(iblock,eqn,varn,"bc_top1",'l',&d,&l,NULL,NULL);
	}
	inline void bc_top1_add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r) {
		add(iblock,eqn,varn,"bc_top1",'r',&d,NULL,&r,NULL);
	}
	inline void bc_top1_add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r) {
		add(iblock,eqn,varn,"bc_top1",'f',&d,&l,&r,NULL);
	}
	inline void bc_top1_add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i) {
		add(iblock,eqn,varn,"bc_top1",'m',&d,&l,NULL,&i);
	}
	inline void bc_top1_add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i) {
		add(iblock,eqn,varn,"bc_top1",'s',&d,NULL,&r,&i);
	}
	inline void bc_top1_add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i) {
		add(iblock,eqn,varn,"bc_top1",'g',&d,&l,&r,&i);
	}
	
	inline void bc_top2_add_d(int iblock,const char *eqn, const char *varn,const matrix &d) {
		add(iblock,eqn,varn,"bc_top2",'d',&d,NULL,NULL,NULL);
	}
	inline void bc_top2_add_l(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l) {
		add(iblock,eqn,varn,"bc_top2",'l',&d,&l,NULL,NULL);
	}
	inline void bc_top2_add_r(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r) {
		add(iblock,eqn,varn,"bc_top2",'r',&d,NULL,&r,NULL);
	}
	inline void bc_top2_add_lr(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r) {
		add(iblock,eqn,varn,"bc_top2",'f',&d,&l,&r,NULL);
	}
	inline void bc_top2_add_li(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &i) {
		add(iblock,eqn,varn,"bc_top2",'m',&d,&l,NULL,&i);
	}
	inline void bc_top2_add_ri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &r,const matrix &i) {
		add(iblock,eqn,varn,"bc_top2",'s',&d,NULL,&r,&i);
	}
	inline void bc_top2_add_lri(int iblock,const char *eqn, const char *varn,const matrix &d,const matrix &l,const matrix &r,const matrix &i) {
		add(iblock,eqn,varn,"bc_top2",'g',&d,&l,&r,&i);
	}
	
};


#endif
