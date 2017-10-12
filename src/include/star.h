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

#define PRES T
#define LOG_PRES "log_T"
//#define PRES p
//#define LOG_PRES "log_p"
#define T_CONSTANT_DOMAINS
//#define PHOTOSPHERE 1
//#define KINEMATIC_VISC

class star1d;

class star2d {
  protected:
	virtual void copy(const star2d &);
	void init1d(const star1d &A,int npts_th,int npts_ex);
	virtual bool check_tag(const char *tag) const;
	virtual void write_tag(OUTFILE *fp) const;
  public:
  	mapping map;
  	const int &nr,&nth,&nex,&ndomains;
	const matrix &r,&z,&th,&Dt,&Dt2,&zex,&Dex,&rex;
	const matrix_block_diag &D;
    matrix rho,phi,p,T,Xh,Xh_prec,Xh0,rho0,Wr,r0,schwarz;
    matrix phiex;
	matrix vr,vt,G,w;
	composition_map comp; 
    double X0,Y0,Z0;
    double X_core,X_core_prec,M_core,M_core_prec;
    double R,M,R0;
    double rhoc,rhoc0,Tc,pc;
    double Omega,Omega_bk,Omegac;
   	double Ekman;
  	opa_struct opa;
	nuc_struct nuc;
	eos_struct eos;
	atm_struct atm;
	matrix ps,Ts;
	double m,pi_c,Lambda;
	double surff;
	int conv,nd_core;
	double Xc;
	int core_convec;
	int env_convec;
	int stratified_comp;
	double min_core_size;
	double global_err;
	int glit; // global view of iteration number
	double dtime,time,delta; // in Myrs
	std::vector<int> izif; // index (<= ndomains) of zones interfaces
	std::vector<int> domain_type;
	std::vector<int> zone_type;
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
	void calc_units();

	star2d();
	virtual ~star2d();
	star2d(const star2d &);
	star2d &operator=(const star2d &);
	struct config_struct {
		double newton_dmax;
		int verbose;
		int dump_iter;
		int input_file;
	} config;
	
	virtual void opacity();
	virtual void nuclear();
	virtual void eq_state();
	virtual void atmosphere();

	virtual int init(const char *input_file,const char *param_file,int argc,char *argv[]);
	virtual int check_arg(char *arg,char *val,int *change_grid);
	virtual int read(const char *input_file, int dim = 2);
	virtual int read_old(const char *input_file);
	virtual void write(const char *output_file,char output_mode='b') const;
	virtual void interp(remapper *red);
	
	virtual void dump_info();
	
	virtual void init_comp();
	
	virtual solver *init_solver(int nvar_add=0);
	virtual double solve(solver *);
	virtual void register_variables(solver *op);
	
	virtual void solve_poisson(solver *);
	virtual void solve_mov(solver *);
	virtual void solve_temp(solver *);
	virtual void solve_dim(solver *);
	virtual void solve_map(solver *);
	virtual void solve_Omega(solver *);
	virtual void solve_gsup(solver *);
	virtual void solve_Teff(solver *);
	virtual void solve_definitions(solver *);
	virtual void solve_atm(solver *);
	virtual void solve_Xh(solver *);
	//virtual void solve_Wr(solver *);
	
	virtual void update_map(matrix dR);
	
	virtual void calc_veloc();
	
	virtual matrix entropy() const;
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
	virtual double Lz() const;
	virtual double Mcore() const;
	virtual double Lzcore() const;
	virtual matrix Rcore() const;
	
	virtual void fill();
	
	// star_map.cpp
	virtual void remap(int ndomains,int *npts,int nth,int nex);
	virtual bool remap_domains(int ndom, remapper &red);
	virtual matrix find_boundaries(const matrix &logTi) const;
	virtual std::vector<int> distribute_domains(int ndom,matrix &zif,bool check_only=false) const;
	virtual matrix distribute_domains(int ndomains,int &conv_new,double p_cc=0) const;
	virtual matrix find_boundaries_old(matrix pif) const;
	virtual void check_map();
	virtual void new_check_map();
	virtual int check_CC(double &p_cc,matrix &Rcc);
        matrix solve_temp_rad();
        //int find_zones(matrix& r_inter, std::vector<int>& zone_type, matrix& p_inter);
        int find_zones(matrix& r_inter, matrix& p_inter);
        //matrix new_distribute_domains(int ndom,matrix p_inter,std::vector<int> zone_type);
        matrix new_distribute_domains(int ndom,matrix p_inter);
        matrix New_distribute_domains(int ndom,matrix p_inter);

	
	void draw(figure *,const matrix &,int parity=0) const;
	void drawi(figure *,const matrix &,int sr,int st,int parity=0) const;
	void drawc(figure *,const matrix &,int ncontours,int parity=0) const;
	void drawci(figure *,const matrix &,int sr,int st,int ncontours,int parity=0) const;
	void spectrum(figure *,const matrix &,int parity=0) const;

	matrix kconv() const;
	void add_kconv(solver *op,const char *eqn,const matrix &d);
	void add_dkconv_dz(solver *op,const char *eqn,const matrix &d);
	// void kconv_common(matrix &kc,matrix &Ja,matrix &Jb,symbolic &S,sym &a_,sym &b_) const;

	virtual void check_jacobian(solver *op,const char *eqn);

    void hdf5_write(const char *filename) const;
    int hdf5_read(const char *input_file, int dim);
};

class star1d : public star2d {
  protected:
    virtual bool check_tag(const char *tag) const;
	virtual void write_tag(OUTFILE *fp) const;
  public:	
  	// star1d_class.cpp
	star1d();
	~star1d();
	star1d(const star1d &);
	star1d &operator=(const star1d &);
	virtual int init(const char *input_file,const char *param_file,int argc,char *argv[]);
	virtual int check_arg(char *arg,char *val,int *change_grid);
	virtual int read_old(const char *input_file);
	virtual int read(const char *input_file, int dim = 1);
	
	virtual void dump_info();
	
	virtual solver *init_solver(int nvar_add=0);
	virtual void register_variables(solver *op);
	virtual double solve(solver *);
	virtual void solve_poisson(solver *);
	virtual void solve_pressure(solver *);
	virtual void solve_temp(solver *);
	virtual void new_solve_temp(solver *);
	virtual void solve_dim(solver *);
	virtual void solve_map(solver *);
	virtual void solve_definitions(solver *);
	virtual void solve_Teff(solver *);
	virtual void solve_gsup(solver *);
	virtual void solve_Xh(solver *);
	virtual void solve_Wr(solver *);
	
	virtual void update_map(matrix dR);

	virtual matrix N2() const;
	virtual double luminosity() const;
	virtual matrix Teff() const;
	virtual matrix gsup() const;
	
	virtual void fill();
	
	void spectrum(figure *,const matrix &,const char *line="") const;
	
	virtual void check_jacobian(solver *op,const char *eqn);
};

class star_evol : public star2d {
protected:
    bool comp_inited;
public:
    bool converged;
	double Lz_obj;
	star_evol();
	star_evol(const star2d &);
	virtual void fill();
	virtual int read(const char *input_file, int dim = 2);
	virtual solver *init_solver(int nvar_add=0);
	virtual void register_variables(solver *op);
	virtual void solve_Omega(solver *);
    // void init_comp();
};

#endif


