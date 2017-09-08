#ifndef _SOLVER_H
#define _SOLVER_H

#include "matrix.h"

/// \page term_representation Term representation in the equations
/// A term in an equation is stored in the form
/// \f$ D \circ (L . I \circ x . R) \f$.
/// Where \f$ . \f$ is the matrix multiplication and \f$ \circ \f$,
/// the Schur product.

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
	explicit solver_full(int nblocks,int offcore=0);
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
	
/// \brief A solver object stores a system of equation to solve.
///
/// A \p solver_elem object contains all the equations for a given block. A block
/// refers to one domain or to one set of boundary conditions.
///
// nv is the number of variables (equations) of the problem.
// eq_set(i) (i<nv) is 1 if the i-th equation is defined in ths block.
// nr[i],nth[i] is the size of the i-th equation (It is automatically taken from
// the size of D).	
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
	
    /// \brief Stores one term of an equation.
    ///
    /// A full linked list of \p solver_elem represents all the terms of a given
    /// variable in a particular equation.
    /// \see \ref term_representation
	class solver_elem {
        /// Indicates which of the matrices
        /// D, L, R, I are not unity
        /// - 'd': D
        /// - 'l': D,L
        /// - 'r': D,R
        /// - 'f': D,L,R
        /// - 'm': D,L,I
        /// - 's': D,R,I
        /// - 'g': D,L,R,I
		char type;
		matrix D,L,R,I;
        /// \brief Pointer to the next term in the equation.
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

    /// \brief Stores all the equations for a given block
    /// (a block can either be a domain or a set of boundary condition).
    ///
    /// The equations are represented as an array of linked list of solver_elem
    /// objects \p eq.
	class solver_block {
        int nv;	  ///< \brief The number of variables (and equation) of the
               	  ///< problem.
        int *nr;  ///< \brief \p nr[i] is the radial resolution of the \p i-th
                  ///< variable
        int *nth; ///< \brief \p nth[i] is the azimuthal resolution of the \p
                  ///< i-th variable.
        /// \brief Pointers to all the terms of all the equations of the block.
        ///
        /// eq[i][j] points to the first \p solver_elem object of the list
        /// representing the terms of the i-th equation associated with the j-th
        /// variable.
		solver_elem ***eq;
        /// \p eq_set(i) is 1 if the i-th equation is defined in this block.
		matrix eq_set;
		friend class solver;
	public:
		solver_block() {};
		~solver_block() {};
        /// \brief Initializes the object with \p nvar variables.
		void init(int nvar);
        /// \brief Destroys the object and frees memory.
		void destroy();
        /// \brief Resets all equations of the block.
		void reset();
        /// \brief Resets the \p ieq-th equation of the block.
		void reset(int ieq);
        /// \brief Introduces a term in an equation of the block.
        ///
        /// \param ieq index of the equation where the term is to be added
        /// \param ivar index the variable of the term
        /// \param type the \p type of term to be added
        /// \param d the D matrix of the term
        /// \param l the L matrix of the term
        /// \param r the R matrix of the term
        /// \param i the I matrix of the term
        /// \see - \ref term_representation
		void add(int ieq,int ivar,char type, const matrix *d, const matrix *l, const matrix *r, const matrix *i);
		
        /// \brief Introduces a term of \p type 'd' for the \p ivar-th variable
        /// in the \p ieq-th equation.
		inline void add_d(int ieq,int ivar,const matrix &d) {
			add(ieq,ivar,'d',&d,NULL,NULL,NULL);
		}
        /// \brief Introduces a term of \p type 'l'.
		inline void add_l(int ieq,int ivar,const matrix &d,const matrix &l) {
			add(ieq,ivar,'l',&d,&l,NULL,NULL);
		}
        /// \brief Introduces a term of \p type 'r'.
		inline void add_r(int ieq,int ivar,const matrix &d,const matrix &r) {
			add(ieq,ivar,'r',&d,NULL,&r,NULL);
		}
        /// \brief Introduces a term of \p type 'f'.
		inline void add_lr(int ieq,int ivar,const matrix &d,const matrix &l,const matrix &r) {
			add(ieq,ivar,'f',&d,&l,&r,NULL);
		}
        /// \brief Introduces a term of \p type 'm'.
		inline void add_li(int ieq,int ivar,const matrix &d,const matrix &l,const matrix &i) {
			add(ieq,ivar,'m',&d,&l,NULL,&i);
		}
        /// \brief Introduces a term of \p type 's'.
		inline void add_ri(int ieq,int ivar,const matrix &d,const matrix &r,const matrix &i) {
			add(ieq,ivar,'s',&d,NULL,&r,&i);
		}
        /// \brief Introduces a term of \p type 'g'.
		inline void add_lri(int ieq,int ivar,const matrix &d,const matrix &l,const matrix &r,const matrix &i) {
			add(ieq,ivar,'g',&d,&l,&r,&i);
		}
	
	};

	int nb; ///< Number of blocks defined in the problem.
    int nv; ///< Number of variables of the problem.
	int **var_nr; ///< \p var_nr[i][j] is the radial resolution of the \p j-th variable in the \p i-th equation 
	int **var_nth; ///< \p var_nr[i][j] is the azimuthal resolution of the \p j-th variable in the \p i-th equation 
    int **var_ntop; ///< \p var_ntop[i][j] is the resolution of the \p i-th top boundary condition on the \p j-th variable
    int **var_nbot; ///< \p var_nbot[i][j] is the resolution of the \p i-th bottom boundary condition on the \p j-th variable
    int solver_N;
    int *def_nr; ///< \p def_nr[i] is the number of radial points in the \p i-th block
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
	void mult_op(matrix *);
	int cgs(const matrix &rhs,matrix &x,int maxit);
	void create();
	void create_full();
	void wrap(const matrix *,matrix *);
	void unwrap(const matrix *,matrix *);
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
// New stuff from Paco to reset BC and IC
	void reset_bc_bot(int iblock, int ieq);
	void reset_bc_top(int iblock, int ieq);
	void reset_bc_bot(int iblock, const char *eq_name);
	void reset_bc_top(int iblock, const char *eq_name);
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
	matrix_map get_vars();
	void solve(int *info=NULL); 
	// info[0]: LU factorization
	// info[1]: Iterative refinement
	// info[2]: CGS iteration
	// info[3]: Number of iterative refinement iterations (-1 if no convergence, -2 if error)
	// info[4]: Number of CGS iterations (-1 if no convergence, -2 if error)
	void mult(matrix *);
	
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


#define RK_END 0
#define RK_STEP 1
#define RK_INTERMEDIATE 2

class RKF_solver {
	int nv;
	char **var;
	int *nr,*nc;
	bool *reg;
	matrix *y,k[6],x,x_eval,b4,b5,c,a,deriv,*dy;
	double t,t_eval;
	bool initd;
	double step;
	int stage;
	void wrap(const matrix *,matrix *);
	void unwrap(const matrix *,matrix *);
	int get_id(const char *varn);
	void check_init();
public:
	double abs_tol,rel_tol;
	RKF_solver();
	void init(int nvar);
	void destroy();
	~RKF_solver() {if(!initd) destroy();};
	void regvar(const char *var_name,const matrix &initial_value);
	void regvars(const matrix_map &vars);
	void set_step(double);
	matrix get_var(const char *var_name);
	matrix_map get_vars();
	double get_t();
	void set_deriv(const char *var_name,const matrix &value);
	void set_derivs(const matrix_map &values);
	
	int solve(double t0,double t1);

};

/*
Runge-Kutta solver for SDIRK (Singly diagonally implicit Runge-Kutta)
or ESDIRK (Explicit singly diagonally implicit Runge-Kutta).

Butcher tableau for SDIRK methods is:

		c0 | a    0    0    0  ... 0
		c1 | a10  a    0    0  ... 0
		c2 | a20  a21  a    0  ... 0
		:  |        ...
		1  | b0   b1   b2   b3 ... a
		   |-------------------------
		   | b0   b1   b2   b3 ... a

and for ESDIRK methods:

		0  | 0    0    0    0  ... 0
		c1 | a10  a    0    0  ... 0
		c2 | a20  a21  a    0  ... 0
		:  |        ...
		1  | b0   b1   b2   b3 ... a
		   |-------------------------
		   | b0   b1   b2   b3 ... a

We want to solve the equation:

		F(t,u,v,...,du/dt,dv/dt,...)=0

First, define a new object:

	-->	SDIRK_solver rk;
		
Initialize solver:
	
	-->	rk.init(nvar,method); // nvar = nb de variable avec deirvee temp. explicite
				      // methode string = cidessous ==> sdirk2,3 ou 4.
		
Where nvar is the number of variables an method can be:

		- "be" Backward Euler (SDIRK) A-stable, L-stable
		- "cn" Crank-Nicolson (ESDIRK) A-stable, not L-stable
		- "sdirk2" SDIRK method of order 2 (Alexander,1977) A-stable, L-stable
		- "sdirk3" SDIRK method of order 3 (Alexander,1977) A-stable, L-stable
		- "sdirk4" SDIRK method of order 4 (Skvortsov,2006) A-stable, L-stable
		- "esdirk3" ESDIRK method of order 3 (Segawa) A-stable, L-stable
		- "esdirk4" ESDIRK method of order 4 (Segawa) A-stable, L-stable
		
Register variables:
	
	-->	rk.regvar("u",u0); u0 =valeur initial de u
	-->	rk.regvar("v",v0);
	-->	...
	
For ESDIRK methods, the initial value of the temporal derivative is also needed:

	--> if(rk.needs_initial_derivative()) {     // true for ESDIRK method
	-->		rk.set_initial_derivative("u",du0);
	-->		rk.set_initial_derivative("v",dv0);
	-->		...}

Set the time step 'h':

	-->	rk.set_step(h);
	
Start main loop (time loop):
	--> int state;
	--> while((state=rk.solve(t0,t1))!=RK_END) { // t0 = t initial, t1=t_final
state est l'état du solveur 0, 1 ou 2 cf +haut (etat final, step done, intermediate step)
	-->		double delta,t;
	-->		matrix u0,v0,...;
	-->		delta=rk.get_delta();
	-->		t=rk.get_t();
	-->		u0=rk.get_var("u");
	-->		v0=rk.get_var("v");
	-->		...
Solve:
		F(t,u,v,...,du/dt,dv/dt,...)=0
for u,v,... using:
		du/dt=(u-u0)/delta
		dv/dt=(v-v0)/delta
		...
and write the solution:

	-->		rk.set_var("u",u);
	-->		rk.set_var("v",v);
	-->		...
Check if the RK step is complete:
	-->		if(state==RK_STEP) {
				... save variables at this step t,u,v,...
	-->		}
End main loop
	-->	}
	
Finalize the solver:
	-->	rk.destroy();		
	
The functions regvar,get_var ans set_var have their equivalents using matrix_map objects
reg_vars, set_vars and get_vars.

*/

class SDIRK_solver {
	int nv,nstages,order;
	double alpha;
	char **var;
	int *nr,*nc;
	bool *reg;
	bool initd,first_explicit;
	double step;
	int stage;
	matrix a,b,c;
	matrix *y,x,*k,*dy,x0;
	double t,t_eval,delta;
	void wrap(const matrix *,matrix *);
	void unwrap(const matrix *,matrix *);
	int get_id(const char *varn);
	void check_init();
	void check_method();
	void init_be();
	void init_cn();
	void init_sdirk2();
	void init_sdirk3();
	void init_sdirk4();
	void init_esdirk3();
	void init_esdirk4();
public:
	//double abs_tol,rel_tol;
	SDIRK_solver();
	void init(int nvar,const char *type);
	void destroy();
	~SDIRK_solver() {if(!initd) destroy();};
	void regvar(const char *var_name,const matrix &initial_value);
	void regvars(const matrix_map &vars);
	void set_initial_derivative(const char *var_name,const matrix &initial_value);
	void set_initial_derivatives(const matrix_map &vars);
	void set_var(const char *var_name,const matrix &value);
	void set_vars(const matrix_map &values);
	void set_step(double);
	double get_step();
	double get_delta();
	matrix get_var(const char *var_name);
	matrix_map get_vars();
	int number_of_stages();
	int number_of_implicit_stages();
	int get_order();
	bool needs_initial_derivative();
	double get_t();
	
	int solve(double t0,double t1);

};


#endif
