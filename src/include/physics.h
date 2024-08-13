#ifndef _PHYSICS_H
#define _PHYSICS_H

#include "matrix.h"

#include <vector>
#include <string> 
#include <map>

struct nuc_struct {
	matrix eps,pp,cno,dlneps_lnrho,dlneps_lnT;
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

struct atm_struct {
	matrix ps,Ts,dlnps_lng,dlnps_lnTeff,dlnTs_lng,dlnTs_lnTeff;
	char name[16];
};

class composition_map : public matrix_map {
public:
	matrix_map dt;
	std::map<std::string,matrix_map> jac; // whats jac 
	composition_map() {};
	composition_map(const matrix_map &map) : matrix_map(map) {};
	matrix X() const;
	matrix Y() const;
	matrix Z() const;
	double Xsol() const;
	double Ysol() const;
	double Zsol() const;
};

class AbundanceMap {
public:
    std::map<std::string, float> comp_abund;
    std::map<std::string, double> A_weights;
    std::string comp_name;
};


extern AbundanceMap global_abundance_map;

int opa_calc(const matrix &X,double Z,const matrix &T,const matrix &rho,
		opa_struct &opa, const double Xsol, const double Ysol, const double Zsol);
int eos_calc(const matrix &X,double Z,const matrix &T,const matrix &p,
		matrix &rho,eos_struct &eos);
int nuc_calc(const matrix_map &X,const matrix &T,const matrix &rho,
		nuc_struct &nuc);
int atm_calc(const matrix &X,double Z,const matrix &g,const matrix &Teff,
		const char *eos_name,const char *opa_name,atm_struct &atm, const double Xsol, const double Ysol, const double Zsol);

int opa_opal(const matrix &X,double Z,const matrix &T,const matrix &rho,
		opa_struct &opa);
int opa_houdek(const matrix &X,double Z,const matrix &T,const matrix &rho,
		opa_struct &opa);
int opa_kramer(const matrix &T,const matrix &rho,
		opa_struct &opa);
int opa_cesam(const matrix& X, double Z, const matrix& T, const matrix& rho,
		opa_struct& opa);

int opa_opmesa(const matrix& X, double Z, const matrix& T, const matrix& rho,
		opa_struct& opa, const double Xsol, const double Ysol, const double Zsol);

		
int nuc_simple(const composition_map &comp,const matrix &T,const matrix &rho,
		nuc_struct &nuc);	
int nuc_cesam(const composition_map &comp,const matrix &T,const matrix &rho,
		nuc_struct &nuc);	
		
int nuc_cesam_dcomp(composition_map &comp,const matrix &T,const matrix &rho,
		nuc_struct &nuc);

int eos_ideal(const matrix &X,double Z,const matrix &T,const matrix &p,
		matrix &rho,eos_struct &eos);
int eos_idealrad(const matrix &X,double Z,const matrix &T,const matrix &p,
		matrix &rho,eos_struct &eos);
int eos_opal(const matrix &X,double Z,const matrix &T,const matrix &p,
		matrix &rho,eos_struct &eos);
int eos_freeeos(const matrix &X, double Z, const matrix &T, const matrix &p,
		matrix &rho, eos_struct &eos);

int atm_onelayer(const matrix &X,double Z,const matrix &g,const matrix &Teff,
		const char *eos_name,const char *opa_name,atm_struct &atm, const double Xsol, const double Ysol, const double Zsol);

double_map initial_composition(double X,double Z);  

#endif

