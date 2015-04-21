#ifndef _NUMDIFF_H
#define _NUMDIFF_H

#include "ester-config.h"

#include "matrix.h"

/// \brief
/// The diff_gl class implements the Gauss-Lobatto (or more properly
/// Gauss-Lobatto-Chebyshev) collocation method.
///
/// The basis functions are Chebyshev polynomials of the first kind:
/// \f$ T_l(x)=cos(l.arccos(x)) \f$,
/// defined in \f$ [-1, 1] \f$. \n
/// And the collocation points are \f$ x_i = -cos(\frac{i \pi}{n})\f$
class diff_gl {
    int ndom,Ntot;
    void init_1();
public:
    matrix x,I;
    matrix_block_diag P,P1,D;
    int *npts;
    double *xif;
    int &ndomains,&N;

    explicit diff_gl(int n=1);
    ~diff_gl();
    diff_gl(const diff_gl &);
    diff_gl &operator=(const diff_gl &);
    void set_ndomains(int n);
    void set_xif(double,...);
    void set_npts(int,...);
    void init();
    matrix eval(const matrix &y,const matrix &x,matrix &T) const;
    matrix eval(const matrix &y,const matrix &x) const;
    matrix eval(const matrix &y,double x) const;
    matrix eval(const matrix &y,double x,matrix &T) const;
    matrix_map eval(const matrix_map &y, const matrix &x) const;
};

/// \brief The diff_leg class implements Legendre numerical differentiation.
///
/// The basis function are Legendre polynomials \f$ P_l(x) \f$, and, for
/// \f$ n \f$ points, collocation points are defined as the roots of
/// \f$ P_n(x) \f$.
/// \f[
/// P_n(x)= \frac{1}{2^n n!} \frac{d^n}{dx^n}[(x^2-1)]
/// \f]
class diff_leg {
	matrix Pn(int n,matrix x);
	matrix dPn(int n,matrix x,matrix p);
	matrix nodes(int n,matrix &w);
  public:
  	int npts;
  	matrix th;
  	matrix P_00,P_01,P_10,P_11;
  	matrix P1_00,P1_01,P1_10,P1_11;
  	matrix dP1_00,dP1_01,dP1_10,dP1_11;
  	matrix D_00,D_01,D_10,D_11;
  	matrix D2_00,D2_01,D2_10,D2_11;
  	matrix I_00;
  	matrix lap_00,lap_01,lap_10,lap_11;
  	
  	diff_leg();
  	~diff_leg();
  	diff_leg(const diff_leg &);
  	diff_leg &operator=(const diff_leg &);
  	void init();
  	matrix eval_00(const matrix &y,const matrix &thi,matrix &T) const;
  	matrix eval_00(const matrix &y,const matrix &thi) const;
  	matrix eval_00(const matrix &y,double thi) const;
  	matrix eval_00(const matrix &y,double thi,matrix &T) const;
  	matrix eval_01(const matrix &y,const matrix &thi,matrix &T) const;
  	matrix eval_01(const matrix &y,const matrix &thi) const;
  	matrix eval_01(const matrix &y,double thi) const;
  	matrix eval_01(const matrix &y,double thi,matrix &T) const;
  	matrix eval_10(const matrix &y,const matrix &thi,matrix &T) const;
  	matrix eval_10(const matrix &y,const matrix &thi) const;
  	matrix eval_10(const matrix &y,double thi) const;
  	matrix eval_10(const matrix &y,double thi,matrix &T) const;
  	matrix eval_11(const matrix &y,const matrix &thi,matrix &T) const;
  	matrix eval_11(const matrix &y,const matrix &thi) const;
  	matrix eval_11(const matrix &y,double thi) const;
  	matrix eval_11(const matrix &y,double thi,matrix &T) const;
  	matrix eval(const matrix &y,const matrix &thi,matrix &T,int par_pol,int par_eq) const;
  	inline matrix eval(const matrix &y,double thi,matrix &T,int par_pol,int par_eq) const 
  		{return eval(y,thi*ones(1,1),T,par_pol,par_eq);}
  	matrix l_00() const;
  	matrix l_01() const;
  	matrix l_10() const;
  	matrix l_11() const;

};



/*
Example:

diff_gl gl(3); //Equivalent to: diff_gl gl; gl.set_ndomains(3);
gl.set_npts(10,20,15);
gl.set_xif(0.,0.2,1.,1.5);
gl.init();

diff_leg leg;

leg.npts=50;
leg.init();


*/



#endif
