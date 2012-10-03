#ifndef _PHYSICS_H
#define _PHYSICS_H

#include"matrix.h"

struct nuc_struct {
	matrix eps,pp,cno,d,dlneps_lnrho,dlneps_lnT;
	char name[16];
};
struct eos_struct {
	matrix G1,cp,del_ad,G3_1,cv,d,prad,chi_rho,chi_T,s;
	char name[16];
};
struct opa_struct {
	matrix k,xi,dlnxi_lnrho,dlnxi_lnT;
	char name[16];
};

int opa_opal(const matrix &X,double Z,const matrix &T,const matrix &rho,
		opa_struct &opa);
int opa_opals(const matrix &X,double Z,const matrix &T,const matrix &rho,
		opa_struct &opa);
int opa_houdek(const matrix &X,double Z,const matrix &T,const matrix &rho,
		opa_struct &opa);
int opa_kramer(const matrix &T,const matrix &rho,
		opa_struct &opa);
int opa_test(const matrix &T,const matrix &rho,
		opa_struct &opa);

int nuc_simple(const matrix &X,double Z,const matrix &T,const matrix &rho,
		nuc_struct &nuc);
int nuc_test(const matrix &X,double Z,const matrix &T,const matrix &rho,
		nuc_struct &nuc);		

int eos_ideal(const matrix &X,double Z,const matrix &T,const matrix &p,
		matrix &rho,eos_struct &eos);
int eos_idealrad(const matrix &X,double Z,const matrix &T,const matrix &p,
		matrix &rho,eos_struct &eos);
int eos_opal(const matrix &X,double Z,const matrix &T,const matrix &p,
		matrix &rho,eos_struct &eos);
int eos_test(const matrix &X,double Z,const matrix &T,const matrix &p,
		matrix &rho,eos_struct &eos);
		
#endif

