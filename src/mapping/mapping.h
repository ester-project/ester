#ifndef _MAPPING_H
#define _MAPPING_H

#include"matrix.h"
#include"numdiff.h"
#include"graphics.h"
#include"solver.h"

class mapping {
	matrix eps_,eta_;
  public:
  	int mode;   // mode=0 Bonazzola, mode=1 rzz=0 at the boundaries 
    diff_gl gl;
    diff_leg leg;
    const matrix_block_diag &D;
    const matrix &z,&Dt,&Dt2,&th,&Dt_odd;
    matrix r,rz,rzz,rt,rtt,rzt,gzz,gzt,gtt;
    matrix R;
    matrix J[4];
    struct ext_map {
    	matrix z,D;
    	matrix r,rz,rzz,rt,rtt,rzt,gzz,gzt,gtt;
    	diff_gl gl;
    	matrix J[4];
    } ex;
    
    mapping(int ndom=1);
  	~mapping();
  	mapping(const mapping &);
  	mapping &operator=(const mapping &);
    int init();
    int remap();
    int remap(const diff_leg&);
    double eps(int) const;
    double eta(int) const;
    matrix eps() const;
    matrix eta() const;
    matrix dr(const matrix &) const;
    matrix dt(const matrix &) const;
    matrix dt_odd(const matrix &) const;
    matrix dt2(const matrix &) const;
    matrix lap(const matrix &) const;
    matrix lap_ex(const matrix &) const;
    void add_lap(solver *op,const char* eqn,const char * varn,const matrix &d,const matrix &phi) const;
    void add_lap_ex(solver *op,const char* eqn,const char * varn,const matrix &d,const matrix &phi) const;
    matrix stream(const matrix &Fz,matrix &Ft) const;
    matrix stream(const matrix &Fz) const;
    
    
    void interps(const mapping &map_old,matrix &Tr,matrix &Tex,
    	matrix &Tt_00,matrix &Tt_01,matrix &Tt_10,matrix &Tt_11) const;
    	
	matrix draw0(const matrix &A,int parity,matrix &x,matrix &y) const;
	void draw(figure *pfig,const matrix &A,int parity=0) const;
	matrix drawi0(const matrix &A,int sr,int st,int parity,matrix &x,matrix &y) const;
	void drawi(figure *pfig,const matrix &A,int sr,int st,int parity=0) const;
	void drawc(figure *pfig,const matrix &A,int ncontours,int parity=0) const;
	void drawci(figure *pfig,const matrix &A,int sr,int st,int ncontours,int parity=0) const;
	void spectrum(figure *pfig,const matrix &y,int parity=0) const;

};

/*
Example:

mapping map(2); 	//Equivalent to: mapping map; map.gl.set_ndomains(2);

map.gl.set_npts(50,30);
map.gl.set_xif(0.,0.5,1.); 		// xif[0] and xif[ndomains] are forced to be 0. and 1.
map.leg.npts=40;
map.ex.gl.set_npts(25); 	// or *(map.ex.gl.npts)=25;

map.init();

map.R.setrow(0,0.2*ones(1,map.leg.npts));
map.R.setrow(1,0.5+0.5*sin(map.leg.th)*sin(map.leg.th)); // R for the last domain should be 1
														// at the equator.
map.remap();    // gl.xif will be changed in remap() to coincide with the polar radii of the
				// domains.

// Change number of angular points

diff_leg leg_new;

leg_new.npts=70;
leg_new.init();
map.remap(leg_new);

// Change number of external radial points

map.ex.gl.set_npts(50);
map.ex.gl.init();
map.remap();

// Change number of radial points

map.gl.set_npts(30,40);
map.remap(); 	// map.gl is initialized inside remap();

// Insert a new domain

matrix R_old;
int i=1;    // New domain will be the i'th domain (first domain is 0'th)

R_old=map.R;
map.gl.set_ndomains(3);
map.gl.set_npts(30,20,25); 		// xif will be set inside remap()
map.R.dim(map.gl.ndomains(),map.leg.npts);
map.R.setblock(0,i-1,0,map.leg.npts-1,R_old.block(0,i-1,0,map.leg.npts-1));
map.R.setblock(i+1,map.R.nrows()-1,0,map.leg.npts-1,R_old.block(i,R_old.nrows()-1,0,map.leg.npts-1));
map.R.setrow(i,0.3+0.1*sin(map.leg.th)*sin(map.leg.th));
map.remap();


// Jacobian

dr=J[0]*d(eta_i)+J[1]*d(eta_(i+1)-eta_i)+J[2]*d(R_i)+J[3]*d(R_(i+1)-R_i);

*/



#endif
