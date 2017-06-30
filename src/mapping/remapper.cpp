#include "ester-config.h"
#include "utils.h"
#include "mapping.h"
#include <stdlib.h>

// Define the remapper class
remapper::remapper(const mapping &map_in) {

	map=map_in;
	ndomains=map.ndomains;
	remapped=0;redist=0;changed_npts=0;
	npts=new int[ndomains];
	fixed=new int[ndomains+1];
	for(int i=0;i<ndomains;i++) npts[i]=map.npts[i];
	for(int i=0;i<=ndomains;i++) fixed[i]=-1;
	R=map.R;
	nt=map.nt;nex=map.nex;
	mode=map.mode;

}

remapper::~remapper() {

	delete [] npts;
	delete [] fixed;

}

void remapper::set_nt(int nt_in) {

	nt=nt_in;
	remapped=0;

}

void remapper::set_nex(int nex_in) {

	nex=nex_in;
	remapped=0;

}

void remapper::set_ndomains(int ndom) {

	if(ndom==ndomains) return;

	ndomains=ndom;
	delete [] npts;
	delete [] fixed;
	npts=new int[ndomains];
	fixed=new int[ndomains+1];
	for(int i=0;i<=ndomains;i++) fixed[i]=-1;
	remapped=0;

}

void remapper::set_npts(int *npts_in) {

	for(int i=0;i<ndomains;i++) npts[i]=npts_in[i];
	changed_npts=1;
	remapped=0;

}

void remapper::set_npts(int npts_in) {

	for(int i=0;i<ndomains;i++) npts[i]=npts_in;
	changed_npts=1;
	remapped=0;

}

void remapper::set_R(const matrix &Rin) {

	R=Rin;
	redist=1;
	remapped=0;

}

// the position of the interface numero "idom_old" is the same as the position of the
// interface numero idom_new
void remapper::set_fixed(int idom_new,int idom_old) {

	fixed[idom_new]=idom_old;
	if(idom_old>map.ndomains) {
		ester_err("(remapper) Index exceeds number of domains");
		exit(1);
	}

}

// mode = Bonazzola ou linear
void remapper::set_mode(int mode_new) {

	mode=mode_new;

}

// get the new mapping...
mapping remapper::get_map() {

	if (!remapped) remap();
	return map_new;

}

// interpolate "y" on the new mapping
// remapper toto(map) // map is the old mapping
// toto.set_npts(new_npts)
// comp=toto.interp(comp)
// parity == equatorial symmetry  (00, 01, 10, 11)

matrix_map remapper::interp(const matrix_map &y,int parity) {
// for the class of matrix_map objects (only the chemical at the moment)

	matrix_map::const_iterator it;
	matrix_map ynew;
	for(it=y.begin();it!=y.end();it++)
		ynew[it->first]=interp(it->second,parity);
	return ynew;

}

matrix remapper::interp(const matrix &y,int parity) {
// same as above but for matrices

	if (!remapped) remap();

	matrix a;

	a=y;
	if(nt!=map.nt) 
		a=(a,Tt[2*(parity/10)+parity%10]);
		
	if(changed_npts&&!redist)
		a=(T,a);
	
	if(redist) {
		matrix anew(map_new.nr,nt);
		for(int i=0;i<nt;i++) {
			anew.setcol(i,map.gl.eval(a.col(i),zmap.col(i)));
			int j0=0;
			for(int n=0;n<ndomains;n++) {
				j0+=map_new.npts[n];
				if(fixed[n+1]>=0) {
					int j1=0;
					for(int k=0;k<fixed[n+1];k++) j1+=map.npts[k];
					anew(j0-1,i)=a(j1-1,i);
					if(n<ndomains-1) anew(j0,i)=a(j1,i);
				}
			}
			anew(0,i)=a(0,i);
		}
		a=anew;
	}	
	
	return a;

}

matrix remapper::interp_ex(const matrix &y,int parity) {

	if (!remapped) remap();

	matrix a;

	a=y;
	if(nt!=map.nt) 
		a=(a,Tt[2*(parity/10)+parity%10]);
		
	if(nex!=map.nex)
		a=(Tex,a);
	
	return a;

}

// low level remapper used by star2d remap
void remapper::remap() {

	if(ndomains!=map.ndomains&&!changed_npts) {
		ester_err("(remapper) Should specify number of points in each domain");
		exit(1);
	}

	if(ndomains!=map.ndomains&&!redist) {
		ester_err("(remapper) Should specify domain boundaries");
		exit(1);
	}
	
	if(R.nrows()!=ndomains+1) {
		ester_err("(remapper) Incorrect size of boundary matrix R");
		exit(1);
	}
	remapped=1;
	
	map_new.set_mode(mode);
	map_new.set_ndomains(ndomains);
	map_new.set_npts(npts);
	map_new.set_nt(nt);
	map_new.set_nex(nex);
	map_new.init();
	
	if(R.ncols()!=nt) {
		if(R.ncols()!=map.nt) {
			ester_err("Error: (remapper) Incorrect size of boundary matrix R");
			exit(1);
		}
		R=map.leg.eval_00(R,map_new.th);
	}
	
	if(ndomains==map.ndomains&&!redist&&mode!=map.mode) {
		for(int n=0;n<=ndomains;n++) fixed[n]=n;
		redist=1;
	}
	
	fixed[ndomains]=map.ndomains;
	fixed[0]=0;
	for(int n=0;n<=ndomains;n++) {
		if(fixed[n]>=0) {
			R.setrow(n,map.leg.eval_00(map.R.row(fixed[n]),map_new.th));
		}
	}
	
	map_new.R=R;
	map_new.remap();
	
	if(nt!=map.nt) {
		map.leg.eval_00(map.th,map_new.th,Tt[0]);
		map.leg.eval_01(map.th,map_new.th,Tt[1]);
		map.leg.eval_10(map.th,map_new.th,Tt[2]);
		map.leg.eval_11(map.th,map_new.th,Tt[3]);
	}
	
	if(changed_npts&&!redist) {
		map.gl.eval(map.z,map_new.z,T);
		for(int i=0,j=0,k=0;i<map_new.ndomains;i++) {
			T.setrow(k,zeros(1,map.nr));T(k,j)=1;
			j+=map.npts[i];
			k+=map_new.npts[i];
			T.setrow(k-1,zeros(1,map.nr));T(k-1,j-1)=1;
		}
	}
	
	if(nex!=map.nex) {
		map.ex.gl.eval(map.ex.gl.x,map_new.ex.gl.x,Tex);
	}
	
	if(!redist) return;
	
	zmap.zero(map_new.nr,map_new.nt);
	
	matrix rth,rzth,Tth,zi,dzi;
	matrix rnew;
	for(int j=nt-1;j>=0;j--) {		
		rnew=map_new.r.col(j);
		rth=map.leg.eval_00(map.r,map_new.th(j),Tth);
		rzth=(map.rz,Tth);
		if(j==nt-1) zi=vector_t(map.z(0),map.z(-1),map_new.nr);
		zi(0)=map.z(0);zi(-1)=map.z(-1);
		int fin=0,nit=0;
		while(fin<3) {
			matrix T,ri,rzi;
			ri=map.gl.eval(rth,zi,T);
			rzi=(T,rzth);
			dzi=-(ri-rnew)/rzi;
			dzi(0)=0;dzi(-1)=0;
			dzi=max(dzi,-(zi-map.z(0))/2);
			dzi=min(dzi,(map.z(-1)-zi)/2);
			if(max(abs(dzi))<1e-10) fin++;
			zi+=dzi;
			nit++;
			if(nit>100) {
				ester_err("(remapper) No convergence in remap()");
				exit(1);
			}
		}
		zmap.setcol(j,zi);
	}
	
}


