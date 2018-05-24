#ifndef _MAPPING_H
#define _MAPPING_H

#include "matrix.h"
#include "numdiff.h"
#include "solver.h"

#define MAP_BONAZZOLA 0
#define MAP_LINEAR 1

class mapping {
	matrix eta_;
  public:
  	diff_gl gl;
    diff_leg leg;
    const int &nr,&nt,&ndomains;
    int mode;
    int *&npts;
    const matrix &eta;
    matrix_block_diag &D;
    matrix &z,&th,&Dt,&Dt2,&Dt_11,&Dt2_11,&Dt_01,&Dt2_01,&Dt_10,&Dt2_10;
    matrix &I,&It;
    matrix r,rz,rzz,rt,rtt,rzt,gzz,gzt,gtt;
    matrix R;
    matrix J[4];
    class ext_map {
    	mapping *parent;
    	ext_map(mapping *);
    	void copy(const ext_map &);
      public:
    	matrix z,D;
    	matrix r,rz,rzz,rt,rtt,rzt,gzz,gzt,gtt;
    	diff_gl gl;
    	matrix J[4];

    	operator mapping();
    	friend class mapping;
    } ex;
    const int &nex;

    mapping();
  	~mapping();
  	mapping(const mapping &);
  	mapping &operator=(const mapping &);
  	void copy(const mapping &);
  	void set_ndomains(int ndom);
  	void set_npts(int npts);
  	void set_nt(int nt);
  	void set_nex(int nex);
  	void set_npts(int *npts);
  	void set_mode(int mode);
    int init();
    int remap();

    matrix eval(const matrix &y,const matrix &ri, const matrix &thi,int parity=0) const;


    matrix dr(const matrix &) const;
    matrix dt(const matrix &) const;
    matrix dt_odd(const matrix &) const;
    matrix dt2(const matrix &) const;

    matrix stream(const matrix &Fz,matrix &Ft) const;
    matrix stream(const matrix &Fz) const;
    matrix zeta_to_r(const matrix &z) const;


	matrix draw0(const matrix &A,int parity,matrix &x,matrix &y) const;
	// void draw(figure *pfig,const matrix &A,int parity=0) const;
	// matrix drawi0(const matrix &A,int sr,int st,int parity,matrix &x,matrix &y) const;
	// void drawi(figure *pfig,const matrix &A,int sr,int st,int parity=0) const;
	// void drawc(figure *pfig,const matrix &A,int ncontours,int parity=0) const;
	// void drawci(figure *pfig,const matrix &A,int sr,int st,int ncontours,int parity=0) const;
	// void spectrum(figure *pfig,const matrix &y,int parity=0) const;

	friend class ext_map;

};

class remapper {
	mapping map,map_new;
	int remapped,redist,changed_npts;
	int nt,ndomains,*npts,*fixed,nex,mode;
	matrix R,zmap,T,Tt[4],Tex;
	void remap();
public:
	remapper(const mapping &map);
	~remapper();
	void set_nt(int nt);
	void set_nex(int nex);
	void set_ndomains(int ndom);
	void set_npts(int *npts);
	void set_npts(int npts);
	void set_R(const matrix &R);
	void set_fixed(int idom_new,int idom_old);
	void set_mode(int mode);
	mapping get_map();
	matrix interp(const matrix &y,int parity=0);
	matrix_map interp(const matrix_map &y,int parity=0);
	matrix interp_ex(const matrix &y,int parity=0);
};

#endif
