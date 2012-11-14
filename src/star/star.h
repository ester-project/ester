#ifndef _STAR_H
#define _STAR_H

#include"matrix.h"
#include"solver.h"
#include"numdiff.h"
#include"mapping.h"
#include"physics.h"
#include"constants.h"
#include"parser.h"
#include"graphics.h"

#include<cmath>

class star1d;

class star2d {
	void copy(const star2d &);
  public:
  	star2d();
  	mapping map;
  	const int &nr,&nth,&nex,&ndomains;
	const matrix &r,&z,&th,&Dt,&Dt2,&zex,&Dex,&rex;
	const matrix_block_diag &D;
    matrix rho,phi,p,T,Xr;
    matrix phiex;
	matrix vr,vt,G,w;
    double X,Y,Z;
    double R,M;
    double rhoc,Tc,pc;
    double Omega,Omega_bk,Omegac;
   	double Ekman;
  	opa_struct opa;
	nuc_struct nuc;
	eos_struct eos;
	char atm_name[16];
	matrix ps,Ts;
	double m,pi_c,Lambda;
	double surff;
	int conv;
	double Xc;

	struct units_struct {
		double rho,p,phi,T,Omega,r,v,F;
	} units;
	void calc_units();

	virtual ~star2d();
	star2d(const star2d &);
	star2d &operator=(const star2d &);
	struct config_struct {
		double newton_dmax,verbose;
	} config;
	
	int opacity();
	int nuclear();
	int eq_state();

	void init(const star1d &A,int npts_th,int npts_ex);
	int init(const char *input_file,const char *param_file,int argc,char *argv[]);
	int check_arg(char *arg,char *val,int *change_grid);
	int read(const char *input_file);
	int read_old(const char *input_file);
	void write(const char *output_file,char output_mode) const;
	void interp(mapping map_old);
	
	
	solver *init_solver();
	double solve(solver *);
	void register_variables(solver *op);
	
	void solve_poisson(solver *);
	void solve_pressure(solver *);
	void solve_temp(solver *);
	void solve_dim(solver *);
	void solve_map(solver *);
	void solve_Omega(solver *);
	void solve_rot(solver *);
	void solve_dyn(solver *);
	void solve_gsup(solver *);
	void solve_Teff(solver *);
	void solve_vbl(solver *,const char *eqn,matrix &rhs);
	void solve_definitions(solver *);
	
	void update_map(matrix dR);
	
	void atmosphere();
	void solve_atm(solver *);
	void atm_simple();
	void solve_atm_simple(solver *);

	
	void upd_Xr();
	void calc_veloc();
	
	double luminosity() const;
	matrix Teff() const;
	matrix N2() const;
	matrix gsup() const;
	double virial_3P() const;
	double virial_L() const;
	double virial_W() const;
	double virial_ps() const;
	double virial() const;
	double energy_test() const;
	matrix stream() const;
	double apparent_luminosity(double i) const;
	
	void fill();
	
	void draw(figure *,const matrix &,int parity=0) const;
	void drawi(figure *,const matrix &,int sr,int st,int parity=0) const;
	void drawc(figure *,const matrix &,int ncontours,int parity=0) const;
	void drawci(figure *,const matrix &,int sr,int st,int ncontours,int parity=0) const;
	void spectrum(figure *,const matrix &,int parity=0) const;

	void check_jacobian(solver *op,const char *eqn);

};

class star1d : public star2d {
	void copy(const star1d &);
  public:	
	star1d();
	~star1d();
	star1d(const star1d &);
	star1d &operator=(const star1d &);
	int init(const char *input_file,const char *param_file,int argc,char *argv[]);
	int check_arg(char *arg,char *val,int *change_grid);
	int read(const char *input_file);
	void write(const char *output_file,char output_mode) const;
	int read_old(const char *input_file);
	
	solver *init_solver();
	void register_variables(solver *op);
	double solve(solver *);
	void solve_poisson(solver *);
	void solve_pressure(solver *);
	void solve_temp(solver *);
	void solve_dim(solver *);
	void solve_map(solver *);
	void solve_definitions(solver *);
	void solve_Teff(solver *);
	void solve_gsup(solver *);
	
	void atmosphere();
	void solve_atm(solver *);
	void atm_simple();
	void solve_atm_simple(solver *);
	
	void update_map(matrix dR);
	
	void upd_Xr();
	
	
	
	matrix N2() const;
	double luminosity() const;
	double Teff() const;
	double gsup() const;
	
	void fill();
	
	void spectrum(figure *,const matrix &,const char *line="") const;
	
	void check_jacobian(solver *op,const char *eqn);
};

#endif


