#ifndef _STARR_H
#define _STARR_H

#include"ester.h"

class star1dR: public star1d {
public:
	star1dR();
	double Teff_obj;
	int check_arg(char *arg,char *val,int *change_grid);
	solver *init_solver(int nvar_add=0);
	void register_variables(solver *);
	double solve(solver *);
	void fill();
	void solve_dim(solver *);
	void solve_definitions(solver *);
	
};

class star2dR: public star2d {
public:
	star2dR();
	double Teff_obj,Re_obj;
	int check_arg(char *arg,char *val,int *change_grid);
	solver *init_solver(int nvar_add=0);
	void register_variables(solver *);
	double solve(solver *);
	void fill();
	void solve_dim(solver *);
	void solve_definitions(solver *);
	void solve_Omega(solver *);
	
};

#endif
