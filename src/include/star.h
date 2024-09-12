#ifndef _STAR_H
#define _STAR_H

#include "matrix.h"
#include "solver.h"
#include "numdiff.h"
#include "mapping.h"
#include "physics.h"
#include "constants.h"
#include "parser.h"
#include "graphics.h"
#include "debug.h"
#include "symbolic.h"

#include <cmath>
#include <vector>

#if __cplusplus > 199711L
#include<thread>
#define THREADS
#endif
#ifdef MKL
#include "mkl.h"
#endif


//#define PRES T
//#define LOG_PRES "log_T"
#define PRES p
#define LOG_PRES "log_p"
#define T_CONSTANT_DOMAINS
//#define PHOTOSPHERE 1
//#define ASYMP_VISC


class star1d;

class star2d {
  protected:
	virtual void copy(const star2d &);
	void init1d(const star1d &A,int npts_th,int npts_ex);
	virtual bool check_tag(const char *tag) const;
	virtual std::string get_tag() const;
	virtual void write_tag(OUTFILE *fp) const;
	virtual void write_eqs(solver *op);
	virtual double update_solution(solver *op, double &h, matrix_map& error, int nit);
	virtual void read_vars(INFILE *fp);
	virtual void write_vars(OUTFILE *fp) const;
  public:
  	mapping map;
  	const int &nr,&nth,&nex,&ndomains;
	const matrix &r,&z,&th,&Dt,&Dt2,&zex,&Dex,&rex;
	const matrix_block_diag &D;
    matrix rho,phi,p,T;
    matrix phiex;
	matrix vr,vt,w;
	matrix vangle;
	matrix N2_prev;
	composition_map comp; 
    double X0,Y0,Z0;
    double R,M;
    double rhoc,Tc,pc;
    double Omega,Omega_bk,Omegac;
   	double reynolds_v, reynolds_h, visc_v, visc_h;
   	double diffusion_v, diffusion_h;
   	double age;
  	opa_struct opa;
	nuc_struct nuc;
	eos_struct eos;
	atm_struct atm;
	matrix ps,Ts;
	double m,pi_c,Lambda;
	double surff;
	int conv;
	double Xc;
	int core_convec;
	int env_convec;
	int stratified_comp;
	double min_core_size;
	std::vector<int> domain_type;
	std::vector<double> domain_weight;
	#define RADIATIVE 0
	#define CORE 1
	#define CONVECTIVE 2
	
	struct version_struct {
		int major, minor, rev, svn;
        std::string name;
	} version;

	struct units_struct {
		double rho,p,phi,T,Omega,r,v,F;
	} units;
	virtual void calc_units();

	star2d();
	virtual ~star2d();
	star2d(const star2d &);
	star2d &operator=(const star2d &);
	struct config_struct {
		double newton_dmax;
		int verbose;
		int dump_iter;
	} config;
	
	virtual void opacity();
	virtual void nuclear();
	virtual void eq_state();
	virtual void atmosphere();

	virtual int init(const char *input_file,const char *param_file,int argc,char *argv[]);
	virtual int check_arg(char *arg,char *val,int *change_grid);
	virtual int read(const char *input_file, int dim = 2);
	//virtual int read_old(const char *input_file);
	virtual void write(const char *output_file,char output_mode='b') const;
	virtual void interp(remapper *red);
	
	virtual void dump_info();
	
	virtual void init_comp();
	
	virtual solver *init_solver(int nvar_add=0);
	virtual double solve(solver *, int);
	virtual double solve(solver *, matrix_map& error, int);
	virtual void register_variables(solver *op);
	
	virtual void solve_poisson(solver *);
	virtual void solve_mov(solver *);
	virtual void solve_cont(solver *);
	virtual void solve_temp(solver *);
	virtual void solve_dim(solver *);
	virtual void solve_map(solver *);
	virtual void solve_Omega(solver *);
	virtual void solve_gsup(solver *);
	virtual void solve_Teff(solver *);
	virtual void solve_definitions(solver *);
	virtual void solve_atm(solver *);
	virtual void solve_vangle(solver *);
	
	virtual void update_map(matrix dR);
	
	
	virtual matrix entropy() const;
	virtual double luminosity() const;
	virtual matrix luminosity_m() const;
	virtual matrix Teff() const;
	virtual matrix N2() const;
	virtual matrix Dmix_v() const;
	virtual matrix Dmix_h() const;
	virtual matrix gsup() const;
	virtual matrix csound() const;
	virtual double virial_3P() const;
	virtual double virial_L() const;
	virtual double virial_W() const;
	virtual double virial_ps() const;
	virtual double virial() const;
	virtual double energy_test() const;
	//virtual matrix stream() const;
	virtual double apparent_luminosity(double i) const;
	virtual double Lz() const;
	virtual double Mcore() const;
	virtual double Lzcore() const;
	virtual double M_dot() const;
	virtual matrix Rcore() const;
	
	virtual void fill();
	virtual void calc_vangle();
	
	// star_map.cpp
	virtual void remap(int ndomains,int *npts,int nth,int nex);
	virtual void check_map();
	virtual int count_zones(std::vector<int> &index);
	virtual matrix get_zone_itfs(matrix &pzone);
	virtual void add_convective_core(double pcc, matrix Rcc);
	virtual int remove_convective_core();
	virtual std::vector<int> distribute_domains(int ndom, matrix pzone, std::vector<int> &index_new);
	virtual std::vector<int> resample_domain_type(int ndom_new, std::vector<int> index_old, std::vector<int> index_new);
	virtual std::vector<double> init_domain_weight(std::vector<int> &domain_type_new);
	virtual matrix get_new_boundaries(matrix pzone, matrix Rzone, std::vector<int> ndom_zone, std::vector<double> weight);
	virtual matrix find_boundary(double pif);
	virtual bool remap_domains(int ndom, remapper &red);
	virtual int check_convec(double &p_cc, matrix &Rcc);
	virtual void update_domain_weights();
	
	// void draw(figure *,const matrix &,int parity=0) const;
	// void drawi(figure *,const matrix &,int sr,int st,int parity=0) const;
	// void drawc(figure *,const matrix &,int ncontours,int parity=0) const;
	// void drawci(figure *,const matrix &,int sr,int st,int ncontours,int parity=0) const;
	// void spectrum(figure *,const matrix &,int parity=0) const;

	matrix kconv() const;
	void add_kconv(solver *op,const char *eqn,const matrix &d);
	void add_dkconv_dz(solver *op,const char *eqn,const matrix &d);
	// void kconv_common(matrix &kc,matrix &Ja,matrix &Jb,symbolic &S,sym &a_,sym &b_) const;


    void hdf5_write(const char *filename) const;
    int hdf5_read(const char *input_file, int dim);

    virtual void plot(const matrix_map&);
};

class star1d : public star2d {
  protected:
    virtual bool check_tag(const char *tag) const;
    virtual std::string get_tag() const;
	virtual void write_tag(OUTFILE *fp) const;
  public:	
  	// star1d_class.cpp
	star1d();
	~star1d();
	star1d(const star1d &);
	star1d &operator=(const star1d &);
	virtual int init(const char *input_file,const char *param_file,int argc,char *argv[]);
	virtual int check_arg(char *arg,char *val,int *change_grid);
	//virtual int read_old(const char *input_file);
	virtual int read(const char *input_file, int dim = 1);
	
	virtual void dump_info();
	
	virtual solver *init_solver(int nvar_add=0);
	virtual void register_variables(solver *op);
	virtual double solve(solver *);
	virtual double solve(solver *, matrix_map& error, int);
	virtual void solve_poisson(solver *);
	virtual void solve_pressure(solver *);
	virtual void solve_temp(solver *);
	virtual void solve_dim(solver *);
	virtual void solve_map(solver *);
	virtual void solve_definitions(solver *);
	virtual void solve_Teff(solver *);
	virtual void solve_gsup(solver *);
	
	virtual void update_map(matrix dR);

	virtual matrix N2() const;

	virtual double luminosity() const;
	virtual matrix luminosity_m() const;
	virtual matrix Teff() const;
	virtual matrix gsup() const;
	
	virtual void fill();
	
	// void spectrum(figure *,const matrix &,const char *line="") const;
	

    virtual void plot(const matrix_map&);

    matrix spectrum(const matrix&);
};

class star_evol : public star2d {
protected:
	matrix Xprev, r0, rho0, T0, lnp0, w0, phi0, XNprev, XOprev, XCprev, N2_prev;
	double lnR0, drhocdX, lnrhoc0, lnTc0, lnpc0;
	double delta;
	bool check_map_enable;
	double mH0, mH;
	virtual void read_vars(INFILE *fp);
	virtual void write_vars(OUTFILE *fp) const;
	virtual void write_eqs(solver *op);
	virtual double update_solution(solver *op, double &h, matrix_map& error, int nit);
	virtual void copy(const star_evol &);
	virtual void calcTimeDerivs();
public:
	matrix dXdt, drhodt, dpdt, dwdt, dTdt, dphidt;
	star_evol();
	star_evol(const star2d &);
	star_evol(const star_evol &);
	star_evol &operator=(const star_evol &);
	virtual void init_comp();
	virtual void fill();
	virtual solver *init_solver(int nvar_add=0);
	virtual SDIRK_solver *init_time_solver();
	virtual void register_variables(solver *op);
	virtual void register_variables(SDIRK_solver *rk);
    virtual void init_step(SDIRK_solver *rk);
    virtual void finish_step(SDIRK_solver *rk, int state);
    virtual void check_map();
    virtual void interp(remapper *red);
    virtual void calc_units();
    virtual void solve_definitions(solver *);
    virtual void solve_dim(solver *);
    virtual void solve_Omega(solver *);
    virtual void solve_X(solver *);
    virtual void solve_XN(solver *);
    virtual void solve_XO(solver *);
    virtual void solve_XC(solver *);
    virtual void solve_cont(solver *);
    virtual void solve_temp(solver *);
    virtual void solve_mov(solver *);
    using star2d::solve;
    virtual double solve(solver *, matrix_map& error, int);

    virtual void reset_time_solver(SDIRK_solver *);

    virtual int remove_convective_core();
};

#endif


