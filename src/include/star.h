#ifndef _STAR_H
#define _STAR_H

#include "matrix.h"
#include "solver.h"
#include "numdiff.h"
#include "mapping.h"
#include "physics.h"
#include "constants.h"
#include "parser.h"
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
        void init1d(const star1d &A, int npts_th, int npts_ex);
        virtual bool check_tag(const char *tag) const;
    public:
        mapping map;
        const int &nr, &nth, &nex, &ndomains;
        const matrix &r, &z, &th, &Dt, &Dt2, &zex, &Dex, &rex;
        const matrix_block_diag &D;
        matrix rho, phi, p, T;
        matrix phiex;
        matrix vr, vt, G, w;
        // the m_ prefix is used to explicit that it's a class' member https://stackoverflow.com/a/13018190/22158934 :
        double m_He_isotopic_ratio; // He3/He4 hardcoded in star*d::init_comp // TODO: change it to be configurable
        double_map m_metal_mix;
        composition_map comp;
        double X0, Y0, Z0;
        double R, M;
        double rhoc, Tc, pc;
        double Omega, Omega_bk, Omegac;
        double Ekman;
        opa_struct opa;
        nuc_struct nuc;
        eos_struct eos;
        atm_struct atm;
        matrix ps, Ts;
        double m, pi_c, Lambda;
        double surff;
        int conv;
        double Xc;
        int core_convec;
        int env_convec;
        int stratified_comp;
        double min_core_size;
        std::vector<int> domain_type;
#define RADIATIVE 0
#define CORE 1
#define CONVECTIVE 2

        struct version_struct {
            int major, minor, rev, svn;
            std::string name;
        } version;

        struct units_struct {
            double rho, p, phi, T, Omega, r, v, F;
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
            bool init_poly;
        } config;

        virtual void opacity();
        virtual void nuclear();
        virtual void eq_state();
        virtual void atmosphere();

        virtual int init(const char *input_file, const char *param_file, int argc, char *argv[]);
        virtual int check_arg(char *arg, char *val, int *change_grid);
        virtual int read(const char *input_file, int dim = 2);
        int hdf5_read(const char *input_file, int dim);
        virtual void write(const char *output_file) const;
        void hdf5_write(const char *filename) const;
        virtual void interp(remapper *red);

        virtual void dump_info();

        virtual void init_metal_mix();
        virtual void init_comp();

        virtual solver *init_solver(int nvar_add=0);
        virtual double solve(solver *);
        virtual double solve(solver *, matrix_map& error, int);
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
        double test_energy,test_virial;
        virtual matrix stream() const;
        virtual double apparent_luminosity(double i) const;
        virtual double Lz() const;
        virtual double Mcore() const;
        virtual double Lzcore() const;
        virtual matrix Rcore() const;
        virtual double Iz() const;
        virtual double Volume() const;
        virtual double Ic() const;
        virtual double J2MR2() const;

        virtual void fill();

        // star_map.cpp
        virtual void remap(int ndomains, int *npts, int nth, int nex);
        virtual bool remap_domains(int ndom, remapper &red);
        virtual matrix find_boundaries(const matrix &logTi) const;
        virtual std::vector<int> distribute_domains(int ndom, matrix &zif, bool check_only=false) const;
        virtual matrix distribute_domains(int ndomains, int &conv_new, double p_cc=0) const;
        virtual matrix find_boundaries_old(matrix pif) const;
        virtual void check_map();
        virtual int check_convec(double &p_cc, matrix &Rcc);

        // void draw(figure *, const matrix &, int parity=0) const;
        // void drawi(figure *, const matrix &, int sr, int st, int parity=0) const;
        // void drawc(figure *, const matrix &, int ncontours, int parity=0) const;
        // void drawci(figure *, const matrix &, int sr, int st, int ncontours, int parity=0) const;
        // void spectrum(figure *, const matrix &, int parity=0) const;

        matrix kconv() const;
        void add_kconv(solver *op, const char *eqn, const matrix &d);
        void add_dkconv_dz(solver *op, const char *eqn, const matrix &d);
        // void kconv_common(matrix &kc, matrix &Ja, matrix &Jb, symbolic &S, sym &a_, sym &b_) const;

        virtual void check_jacobian(solver *op, const char *eqn);

        matrix solve_phi(); // used to calculate initial solution for phi based on rho

        virtual void plot(const matrix_map&);
};

class star1d : public star2d {
    protected:
        virtual bool check_tag(const char *tag) const;
    public:
        // star1d_class.cpp
        star1d();
        ~star1d();
        star1d(const star1d &);
        star1d &operator=(const star1d &);
        virtual int init(const char *input_file, const char *param_file, int argc, char *argv[]);
        virtual int check_arg(char *arg, char *val, int *change_grid);
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
        virtual matrix Teff() const;
        virtual matrix gsup() const;

        virtual void fill();

        // void spectrum(figure *, const matrix &, const char *line="") const;

        virtual void check_jacobian(solver *op, const char *eqn);

        virtual void plot(const matrix_map&);

        matrix spectrum(const matrix&);
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


matrix solve_poly1d(double n, double tol, int nr, double hsurf);

#endif

