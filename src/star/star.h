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
	void init1d(const star1d &A,int npts_th,int npts_ex);
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

	virtual int init(const char *input_file,const char *param_file,int argc,char *argv[]);
	virtual int check_arg(char *arg,char *val,int *change_grid);
	virtual int read(const char *input_file);
	int read_old(const char *input_file);
	virtual void write(const char *output_file,char output_mode) const;
	virtual void interp(mapping map_old);
	
	
	virtual solver *init_solver();
	virtual double solve(solver *);
	virtual void register_variables(solver *op);
	
	virtual void solve_poisson(solver *);
	virtual void solve_pressure(solver *);
	virtual void solve_temp(solver *);
	virtual void solve_dim(solver *);
	virtual void solve_map(solver *);
	virtual void solve_Omega(solver *);
	virtual void solve_rot(solver *);
	virtual void solve_dyn(solver *);
	virtual void solve_gsup(solver *);
	virtual void solve_Teff(solver *);
	virtual void solve_vbl(solver *,const char *eqn,matrix &rhs);
	virtual void solve_definitions(solver *);
	
	virtual void update_map(matrix dR);
	
	virtual void atmosphere();
	virtual void solve_atm(solver *);
	virtual void atm_simple();
	virtual void solve_atm_simple(solver *);

	
	virtual void upd_Xr();
	virtual void calc_veloc();
	
	virtual double luminosity() const;
	virtual matrix Teff() const;
	virtual matrix N2() const;
	virtual matrix gsup() const;
	virtual double virial_3P() const;
	virtual double virial_L() const;
	virtual double virial_W() const;
	virtual double virial_ps() const;
	virtual double virial() const;
	virtual double energy_test() const;
	virtual matrix stream() const;
	virtual double apparent_luminosity(double i) const;
	
	virtual void fill();
	
	void draw(figure *,const matrix &,int parity=0) const;
	void drawi(figure *,const matrix &,int sr,int st,int parity=0) const;
	void drawc(figure *,const matrix &,int ncontours,int parity=0) const;
	void drawci(figure *,const matrix &,int sr,int st,int ncontours,int parity=0) const;
	void spectrum(figure *,const matrix &,int parity=0) const;

	virtual void check_jacobian(solver *op,const char *eqn);

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
	matrix Teff() const;
	matrix gsup() const;
	
	void fill();
	
	void spectrum(figure *,const matrix &,const char *line="") const;
	
	void check_jacobian(solver *op,const char *eqn);
};

#endif


