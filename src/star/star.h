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

#include<math.h>

class star1d {
  public:
	diff_gl gl;
	matrix &r;
	matrix_block_diag &D;
	matrix rho,phi,p,T,Xr;
	matrix Frad;
	opa_struct opa;
	nuc_struct nuc;
	eos_struct eos;
	char atm_name[16];
	
	double ps,Ts,m,pi_c,Lambda;
	double surff;
	int conv;
	double Xc;
	double rhoc,Tc,pc;
	double R,M;
	double X,Z;
	struct units_struct {
		double rho,p,phi,T,r,F;
	} units;
	struct config_struct {
		double newton_dmax,verbose;
	} config;
	
	star1d();
	~star1d();
	star1d(const star1d &);
	star1d &operator=(const star1d &);
	int init(const char *input_file,const char *param_file,int argc,char *argv[]);
	int check_arg(char *arg,char *val,int *change_grid);
	int read(const char *input_file);
	void write(const char *output_file,char output_mode) const;
	
	int nr() const;
	int ndomains() const;
	
	solver *init_solver();
	void register_variables(solver *op);
	double solve(solver *);
	void solve_poisson(solver *);
	void solve_pressure(solver *);
	void solve_temp(solver *);
	void solve_dim(solver *);
	void solve_map(solver *);
	
	void atmosphere();
	void solve_atm(solver *);
	void atm_simple();
	void solve_atm_simple(solver *);
	void atm_test();
	void solve_atm_test(solver *);
	
	void update_map(matrix dR);
	
	void upd_Xr();
	
	int opacity();
	int nuclear();
	int eq_state();
	
	void calc_units();
	void calc_Frad();
	
	matrix N2() const;
	double luminosity() const;
	double Teff() const;
	
	void fill();
	
	void spectrum(figure *,const matrix &,const char *line="") const;
};

class star2d {
  public:
	mapping map;
	matrix &r,&z,&th,&Dt,&Dt2,&zex,&Dex,&rex;
	matrix_block_diag &D;
	matrix rho,phi,p,T,Xr,phiex,Frad;
	matrix vr,vt,G,w,psi;
	opa_struct opa;
	nuc_struct nuc;
	eos_struct eos;
	char atm_name[16];
	double m,pi_c,Lambda;
	matrix ps,Ts;
	double surff;
	int conv;
	double Omega,Omega2,Omega_bk,Omegac;
	double Xc;
	double rhoc,Tc,pc;
	double R,M;
	double X,Z;
	double Ekman;
	struct units_struct {
		double rho,p,phi,T,Omega,r,v,F;
	} units;
	struct config_struct {
		double newton_dmax,verbose;
	} config;
	
	
	matrix gsup_;
	
	star2d();
	~star2d();
	star2d(const star2d &);
	star2d &operator=(const star2d &);
	void init(const star1d &A,int npts_th,int npts_ex);
	int init(const char *input_file,const char *param_file,int argc,char *argv[]);
	int check_arg(char *arg,char *val,int *change_grid);
	int read(const char *input_file);
	void write(const char *output_file,char output_mode) const;
	void interp(mapping map_old);
	
	int nr() const;
	int nth() const;
	int nex() const;
	int ndomains() const;
	
	
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
	void solve_vbl(solver *,const char *eqn,matrix &rhs);
	void solve_cont(solver *);
	
	void atmosphere();
	void solve_atm(solver *);
	void atm_simple();
	void solve_atm_simple(solver *);
	void atm_test();
	void solve_atm_test(solver *);
	
	int opacity();
	int nuclear();
	int eq_state();
	
	void upd_Xr();
	void calc_W2();
	void calc_Frad();
	void calc_units();
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
	void boundary_layer() const;
	
	void fill();
	
	
	void add_dr(solver *,const char *eqn,const matrix &);
	void add_drz(solver *,const char *eqn,const matrix &);
	void add_drt(solver *,const char *eqn,const matrix &);
	void add_drzz(solver *,const char *eqn,const matrix &);
	void add_drzt(solver *,const char *eqn,const matrix &);
	void add_drtt(solver *,const char *eqn,const matrix &);
	void add_drex(solver *,const char *eqn,const matrix &);
	void add_drzex(solver *,const char *eqn,const matrix &);
	void add_drtex(solver *,const char *eqn,const matrix &);
	void add_drzzex(solver *,const char *eqn,const matrix &);
	void add_drztex(solver *,const char *eqn,const matrix &);
	void add_drttex(solver *,const char *eqn,const matrix &);
	
	void add_dr(solver *,int iblock,const char *eqn,const matrix &);
	void add_drz(solver *,int iblock,const char *eqn,const matrix &);
	void add_drt(solver *,int iblock,const char *eqn,const matrix &);
	void add_drzz(solver *,int iblock,const char *eqn,const matrix &);
	void add_drzt(solver *,int iblock,const char *eqn,const matrix &);
	void add_drtt(solver *,int iblock,const char *eqn,const matrix &);
	
	void add_bc_bot1_dr(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_bot2_dr(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_top1_dr(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_top2_dr(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_bot1_drz(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_bot2_drz(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_top1_drz(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_top2_drz(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_bot1_drt(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_bot2_drt(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_top1_drt(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_top2_drt(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_bot1_drzz(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_bot2_drzz(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_top1_drzz(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_top2_drzz(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_bot1_drzt(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_bot2_drzt(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_top1_drzt(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_top2_drzt(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_bot1_drtt(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_bot2_drtt(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_top1_drtt(solver *,int iblock,const char *eqn,const matrix &);
	void add_bc_top2_drtt(solver *,int iblock,const char *eqn,const matrix &);
	
	void draw(figure *,const matrix &,int parity=0) const;
	void drawi(figure *,const matrix &,int sr,int st,int parity=0) const;
	void spectrum(figure *,const matrix &,int parity=0) const;

	void check_jacobian(solver *op,const char *eqn);

};

#endif


