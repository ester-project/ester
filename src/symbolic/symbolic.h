#ifndef _SYMBOLIC_H
#define _SYMBOLIC_H

#include"matrix.h"

class symbolic {
public:
	int initd;
	int N;
	int nsc,nvec,maxder,****jsc,*****jvec,*****jvec_,**jr,jsint,jcost,jnum;
	char **scname,**vecname,**sccodename,**veccodename;
	int **scparity,***vecparity;
	
	void write2(const char* var,int ex);
	void write3(const char* var,int ex);
	void write_frac(double f);

	symbolic() {initd=0;};
	void create(int n_sc,int n_vec,int max_der);
	void destroy();
	~symbolic() {if(initd) destroy();};
	void set_scname(int isc,const char *name);
	void set_vecname(int ivec,const char *name);
	void set_sccodename(int isc,const char *name);
	void set_veccodename(int ivec,const char *name);
	void set_scparity(int isc,int par0,int par1);
	void set_vecparity(int ivec,int icomp,int par0,int par1);
	void init(matrix &g,int n);
	matrix g(int i,int j);
	matrix g_(int i,int j);
	matrix r(int exp=1);
	matrix sint(int exp=1);
	matrix cost(int exp=1);
	matrix sc(int isc,int exp=1);
	matrix vec(int ivec,int icomp);
	matrix vec_(int ivec,int icomp);
	matrix perm(int i,int j,int k);
	matrix perm_(int i,int j,int k);
	matrix simplify(matrix a);
	matrix simplify_trig(matrix a);
	int simplify_trig_1(matrix &a);
	int simplify_trig_2(matrix &a);
	matrix add(matrix a,matrix b);
	matrix neg(matrix a);
	matrix prod(matrix a,matrix b);
	void write(matrix a);
	void write_group_sc(int isc,matrix a);
	void write_group_sc_u(int isc,matrix a,int start=1);
	void write_group_vec(int ivec,matrix a);
	void write_group_vec_(int ivec,matrix a);
	void writecode(matrix a);
	void writecode_frac(double f);
	void writecode2(char *str,int derz,int dert,int par0,int par1);
	void writecode3(char *str,int n,int *start);
	int isnull(matrix a);
	matrix derive(matrix a,int var);
	matrix subst_sc(int isc,matrix a,matrix x);
	matrix subst_vec(int ivec,matrix a,matrix xz,matrix xt,matrix xp);
	matrix subst_vec_(int ivec,matrix a,matrix xz,matrix xt,matrix xp);
	void covariant(matrix &Az,matrix &At,matrix &Ap);
	void contravariant(matrix &Az,matrix &At,matrix &Ap);
	matrix subst_covariant(int ivec,matrix a);
	matrix subst_contravariant(int ivec,matrix a);
	matrix christoffel(int i,int j,int k);
	matrix gradient(matrix a,int i); // input: scalar, output: covariant
	matrix laplacian(matrix B); // input: scalar, output: scalar
	matrix curl(matrix v_z,matrix v_t,matrix v_p,int i); // input: covariant, output: contravariant
	matrix vect_laplacian(matrix vz,matrix vt,matrix vp,int i); // input: contravariant, output: contravariant
	matrix div(matrix vz,matrix vt,matrix vp); // input: contravariant, output: scalar
	matrix div_tens(matrix T[3][3],int i); // input: contra,contra; output: contra
	matrix advec(matrix wz,matrix wt,matrix wp,matrix vz,matrix vt,matrix vp,int i); // (w.D)v input: contra,contra; output: contra
	matrix dot_prod(matrix wz,matrix wt,matrix wp,matrix vz,matrix vt,matrix vp); //w.v input: contra,co or co,contra; output: scalar
	matrix cross_prod(matrix w_z,matrix w_t,matrix w_p,matrix v_z,matrix v_t,matrix v_p,int i); // wxv input: co,co; output: contra
	matrix axisymmetric(matrix a);
	matrix spherical(matrix a);
	matrix factor_sc(matrix a,int isc,int derz,int dert,int derp);
	matrix factor_vec(matrix a,int ivec,int icomp,int derz,int dert,int derp);
	matrix factor_vec_(matrix a,int ivec,int icomp,int derz,int dert,int derp);
	matrix stress(matrix vz,matrix vt,matrix vp,int i,int j) ; //intput: contra; output contra,contra

};

#endif
