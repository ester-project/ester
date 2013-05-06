#include"mapping.h"
#include<stdlib.h>

mapping_redist::mapping_redist(const mapping &map_in) {

	map=map_in;
	ndomains=map.ndomains;
	remapped=0;redist=0;changed_npts=0;
	npts=new int[ndomains];
	fixed=new int[ndomains];
	for(int i=0;i<ndomains;i++) npts[i]=map.gl.npts[i];
	for(int i=0;i<ndomains;i++) fixed[i]=-1;
	R=map.R;
	nth=map.nth;nex=map.nex;

}

mapping_redist::~mapping_redist() {

	delete [] npts;
	delete [] fixed;

}

void mapping_redist::set_nth(int nth_in) {

	nth=nth_in;
	remapped=0;

}

void mapping_redist::set_nex(int nex_in) {

	nex=nex_in;
	remapped=0;

}

void mapping_redist::set_ndomains(int ndom) {

	if(ndom==ndomains) return;

	ndomains=ndom;
	delete [] npts;
	delete [] fixed;
	npts=new int[ndomains];
	fixed=new int[ndomains];
	for(int i=0;i<ndomains;i++) fixed[i]=-1;
	remapped=0;

}

void mapping_redist::set_npts(int *npts_in) {

	for(int i=0;i<ndomains;i++) npts[i]=npts_in[i];
	changed_npts=1;
	remapped=0;

}

void mapping_redist::set_R(const matrix &Rin) {

	R=Rin;
	redist=1;
	remapped=0;

}

void mapping_redist::set_fixed(int idom_new,int idom_old) {

	fixed[idom_new]=idom_old;
	if(idom_old>=map.ndomains) {
		fprintf(stderr,"Error: (mapping_redist) Index exceeds number of domains\n");
		exit(1);
	}

}

mapping mapping_redist::get_map() {

	if (!remapped) remap();
	return map_new;

}

matrix_map mapping_redist::interp(const matrix_map &y,int parity) {

	matrix_map::const_iterator it;
	matrix_map ynew;
	for(it=y.begin();it!=y.end();it++)
		ynew[it->first]=interp(it->second,parity);
	return ynew;

}

matrix mapping_redist::interp(const matrix &y,int parity) {

	if (!remapped) remap();

	matrix a;

	a=y;
	if(nth!=map.nth) 
		a=(a,Tt[2*(parity/10)+parity%10]);
		
	if(changed_npts&&!redist)
		a=(T,a);
	
	if(redist) {
		matrix anew(map_new.nr,nth);
		for(int i=0;i<nth;i++) {
			anew.setcol(i,map.gl.eval(a.col(i),zmap.col(i)));
			int j0=0;
			for(int n=0;n<ndomains;n++) {
				j0+=map_new.gl.npts[n];
				if(fixed[n]>=0) {
					int j1=0;
					for(int k=0;k<=fixed[n];k++) j1+=map.gl.npts[k];
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

matrix mapping_redist::interp_ex(const matrix &y,int parity) {

	if (!remapped) remap();

	matrix a;

	a=y;
	if(nth!=map.nth) 
		a=(a,Tt[2*(parity/10)+parity%10]);
		
	if(nex!=map.nex)
		a=(Tex,a);
	
	return a;

}

void mapping_redist::remap() {

	if(ndomains!=map.ndomains&&!changed_npts) {
		fprintf(stderr,"Error: (mapping_redist) Should specify number of points in each domain\n");
		exit(1);
	}

	if(ndomains!=map.ndomains&&!redist) {
		fprintf(stderr,"Error: (mapping_redist) Should specify domain boundaries\n");
		exit(1);
	}
	
	if(R.nrows()!=ndomains) {
		fprintf(stderr,"Error: (mapping_redist) Incorrect size of boundary matrix R\n");
		exit(1);
	}
	remapped=1;
	
	map_new.leg.npts=nth;
	map_new.leg.init();
	if(R.ncols()!=nth) {
		if(R.ncols()!=map.nth) {
			fprintf(stderr,"Error: (mapping_redist) Incorrect size of boundary matrix R\n");
			exit(1);
		}
		R=map.leg.eval_00(R,map_new.leg.th);
	}
	
	fixed[ndomains-1]=map.ndomains-1;
	for(int n=0;n<ndomains;n++) {
		if(fixed[n]>=0) {
			R.setrow(n,map.leg.eval_00(map.R.row(fixed[n]),map_new.leg.th));
		}
	}
	
	map_new.gl.set_ndomains(ndomains);
	for(int n=0;n<ndomains;n++) map_new.gl.npts[n]=npts[n];
	map_new.gl.xif[0]=0;
	for(int n=0;n<ndomains-1;n++) map_new.gl.xif[n+1]=map_new.leg.eval_00(R.row(n),0)(0);
	map_new.gl.xif[ndomains]=1;
	map_new.ex.gl.npts[0]=nex;
	map_new.init();
	map_new.R=R;
	map_new.remap();
	
	if(nth!=map.nth) {
		map.leg.eval_00(map.th,map_new.th,Tt[0]);
		map.leg.eval_01(map.th,map_new.th,Tt[1]);
		map.leg.eval_10(map.th,map_new.th,Tt[2]);
		map.leg.eval_11(map.th,map_new.th,Tt[3]);
	}
	
	if(changed_npts&&!redist) {
		map.gl.eval(map.z,map_new.z,T);
		for(int i=0,j=0,k=0;i<map_new.gl.ndomains;i++) {
			T.setrow(k,zeros(1,map.nr));T(k,j)=1;
			j+=map.gl.npts[i];
			k+=map_new.gl.npts[i];
			T.setrow(k-1,zeros(1,map.nr));T(k-1,j-1)=1;
		}
	}
	
	if(nex!=map.nex) {
		map.ex.gl.eval(map.ex.gl.x,map_new.ex.gl.x,Tex);
	}
	
	if(!redist) return;
	
	zmap.zero(map_new.nr,map_new.nth);
	
	matrix rth,rzth,Tth,zi,dzi;
	matrix rnew;
	for(int j=nth-1;j>=0;j--) {		
		rnew=map_new.r.col(j);
		rth=map.leg.eval_00(map.r,map_new.th(j),Tth);
		rzth=(map.rz,Tth);
		if(j==nth-1) zi=rnew/max(rth);
		zi(0)=0;zi(-1)=1;
		int fin=0,nit=0;
		while(fin<3) {
			matrix T,ri,rzi;
			ri=map.gl.eval(rth,zi,T);
			rzi=(T,rzth);
			dzi=-(ri-rnew)/rzi;
			dzi(0)=0;dzi(-1)=0;
			dzi=max(dzi,-zi/2);
			dzi=min(dzi,(1-zi)/2);
			if(max(abs(dzi))<1e-10) fin++;
			zi+=dzi;
			nit++;
			if(nit>100) {
				fprintf(stderr,"Error: (mapping_redist) No convergence in remap()\n");
				exit(1);
			}
		}
		zmap.setcol(j,zi);
	}
}


