// To keep compatibility with configure
#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "utils.h"
#include "mapping.h"
#include "constants.h"
#include <stdlib.h>
#include <cmath>

// defintion of mapping class
mapping::mapping() : nr(gl.N), nt(leg.npts), ndomains(gl.ndomains),
    npts(gl.npts), eta(eta_), D(gl.D), z(gl.x), th(leg.th), Dt(leg.D_00), Dt2(leg.D2_00),
    Dt_11(leg.D_11), Dt2_11(leg.D2_11), Dt_01(leg.D_01), Dt2_01(leg.D2_01),
    Dt_10(leg.D_10), Dt2_10(leg.D2_10), I(gl.I),
    It(leg.I_00), ex(this), nex(ex.gl.N) {

	ex.gl.set_ndomains(1);
	ex.gl.set_xif(0.,1.);
	ex.gl.set_npts(10);

	mode=MAP_BONAZZOLA;
}

mapping::~mapping() {}

mapping::mapping(const mapping &map) : nr(gl.N), nt(leg.npts),
    ndomains(gl.ndomains), npts(gl.npts), eta(eta_), D(gl.D), z(gl.x),
    th(leg.th), Dt(leg.D_00), Dt2(leg.D2_00), Dt_11(leg.D_11),
    Dt2_11(leg.D2_11), Dt_01(leg.D_01), Dt2_01(leg.D2_01), Dt_10(leg.D_10),
    Dt2_10(leg.D2_10), I(gl.I), It(leg.I_00), ex(this), nex(ex.gl.N) {

	copy(map); // copy constructor
}

void mapping::copy(const mapping &map) {

	gl=map.gl;leg=map.leg;
	mode=map.mode;
	eta_=map.eta_;
	r=map.r;rz=map.rz;rzz=map.rzz;rt=map.rt;rtt=map.rtt;rzt=map.rzt;
	gzz=map.gzz;gzt=map.gzt;gtt=map.gtt;
	R=map.R;
	J[0]=map.J[0];J[1]=map.J[1];J[2]=map.J[2];J[3]=map.J[3];
	ex.copy(map.ex);

}

mapping &mapping::operator=(const mapping &map) {

	copy(map);

	return *this;
}

void mapping::set_ndomains(int ndom) {

	gl.set_ndomains(ndom);

}

void mapping::set_npts(int n) {

	for(int i=0;i<ndomains;i++) gl.npts[i]=n;

}

void mapping::set_npts(int *n) {

	for(int i=0;i<ndomains;i++) gl.npts[i]=n[i];
}

void mapping::set_nt(int n) {

	leg.npts=n;

}

void mapping::set_nex(int n) {

	ex.gl.set_npts(n);

}

void mapping::set_mode(int mode_set) {

	mode=mode_set;

}

int mapping::init() {

	R=vector_t(0.,1.,ndomains+1)*ones(1,nt);
	leg.init();
	ex.gl.init();
	return remap();

}

int mapping::remap() {

	eta_=leg.eval_00(R,0);
	for(int i=0;i<ndomains+1;i++)
		gl.xif[i]=eta_(i);
	gl.init();

	if(!exist(R.row(0)!=0)) eta_(0)=0;

	r.dim(nr,nt);
	rz.dim(nr,nt);
	rzz.dim(nr,nt);
	rt.dim(nr,nt);
	rtt.dim(nr,nt);
	rzt.dim(nr,nt);
	J[0].dim(nr,nt);
	J[1].dim(nr,nt);
	J[2].dim(nr,nt);
	J[3].dim(nr,nt);

	if(mode==MAP_BONAZZOLA) {
		int j0=0,dj;
		matrix q,dR,zz,xi,A,Ap,App;
		double deta,a;
	FILE *fic=fopen("rz.txt", "a");
		for(int i=0;i<ndomains;i++) {
			dj=npts[i];
			zz=gl.x.block(j0,j0+dj-1,0,0);
			deta=eta(i+1)-eta(i);
			dR=R.row(i+1)-R.row(i);
	if (dR(0) < 0) {
	   fprintf(fic,"in mapping::remap dR(%d) = %e\n",i,dR(0));
	   printf("dR < 0\n");
	   return 1;
	}
			xi=(zz-eta(i))/deta;
			if(i==0) {
				A=2.5*xi*xi*xi-1.5*xi*xi*xi*xi*xi;
				Ap=7.5*xi*xi-7.5*xi*xi*xi*xi;
				App=15*xi-30*xi*xi*xi;
			} else {
				A=-2*xi*xi*xi+3*xi*xi;
				Ap=-6*xi*xi+6*xi;
				App=-12*xi+6;
			}
			a=1; // This can cause problems when min(dR) is much lower than dR(pole), but it is required
					// for the jacobian of the mapping to be smooth
			q=a*xi*deta+R.row(i)+A*(dR-a*deta);
			r.setblock(j0,j0+dj-1,0,-1,q);
			q=a+Ap*(dR/deta-a);
			rz.setblock(j0,j0+dj-1,0,-1,q);
			q=App/deta*(dR/deta-a);
			rzz.setblock(j0,j0+dj-1,0,-1,q);
			q=(R.row(i),Dt)+A*(dR,Dt);
			rt.setblock(j0,j0+dj-1,0,-1,q);
			q=(R.row(i),Dt2)+A*(dR,Dt2);
			rtt.setblock(j0,j0+dj-1,0,-1,q);
			q=Ap*(dR,Dt)/deta;
			rzt.setblock(j0,j0+dj-1,0,-1,q);
			q=zeros(dj,nt);
			J[0].setblock(j0,j0+dj-1,0,-1,q);
			q=a*(xi-A)*ones(1,nt);
			J[1].setblock(j0,j0+dj-1,0,-1,q);
			q=ones(dj,nt);
			J[2].setblock(j0,j0+dj-1,0,-1,q);
			q=A*ones(1,nt);
			J[3].setblock(j0,j0+dj-1,0,-1,q);
			j0+=dj;
		}
	fclose(fic);
	} else {
		int j0=0,dj;
		matrix q,dR,zz,xi;
		double deta;
		for(int i=0;i<ndomains;i++) {
			dj=npts[i];
			zz=gl.x.block(j0,j0+dj-1,0,0);
			deta=eta(i+1)-eta(i);
			dR=R.row(i+1)-R.row(i);
			xi=(zz-eta(i))/deta;
			q=R.row(i)+xi*dR;
			r.setblock(j0,j0+dj-1,0,-1,q);
			q=dR/deta+zeros(dj,1);
			rz.setblock(j0,j0+dj-1,0,-1,q);
			q=zeros(dj,nt);
			rzz.setblock(j0,j0+dj-1,0,-1,q);
			q=(R.row(i),Dt)+xi*(dR,Dt);
			rt.setblock(j0,j0+dj-1,0,-1,q);
			q=(R.row(i),Dt2)+xi*(dR,Dt2);
			rtt.setblock(j0,j0+dj-1,0,-1,q);
			q=(dR,Dt)/deta+zeros(dj,1);
			rzt.setblock(j0,j0+dj-1,0,-1,q);
			q=zeros(dj,nt);
			J[0].setblock(j0,j0+dj-1,0,-1,q);
			q=zeros(dj,nt);
			J[1].setblock(j0,j0+dj-1,0,-1,q);
			q=ones(dj,nt);
			J[2].setblock(j0,j0+dj-1,0,-1,q);
			q=xi*ones(1,nt);
			J[3].setblock(j0,j0+dj-1,0,-1,q);
			j0+=dj;
		}

	}

    gzz=(r*r+rt*rt)/r/r/rz/rz;
    gzz.setrow(0, 1/rz.row(0)/rz.row(0));
    gzt=-rt/r/r/rz;
    gzt.setrow(0, zeros(1,leg.npts));
    gtt=1/r/r;

	matrix xi;
	ex.z=eta(-1)/(1-ex.gl.x);
	xi=ex.z/eta(-1)-1;

	ex.r=xi+R.row(ndomains);
	ex.rz=ones(nex,nt)/eta(-1);
	ex.rzz=zeros(nex,nt);
	ex.rt=ones(nex,nt)*(R.row(ndomains),Dt);
	ex.rtt=ones(nex,nt)*(R.row(ndomains),Dt2);
	ex.rzt=zeros(nex,nt);

	ex.J[0]=zeros(nex,nt);
	ex.J[1]=zeros(nex,nt);
	ex.J[2]=ones(nex,nt);
	ex.J[3]=zeros(nex,nt);


	ex.gzz=(1+ex.rt*ex.rt/ex.r/ex.r)/ex.rz/ex.rz;
	ex.gzt=-ex.rt/ex.r/ex.r/ex.rz;
	ex.gtt=1/ex.r/ex.r;

	ex.D=ex.gl.D*eta(-1)/ex.z/ex.z;

// check that rz always positive otherwise warning!
// used in update_map to control the relaxation param.

	if(exist(rz<0)||exist(ex.rz<0)) {
		ester_warn("(mapping) Found rz<0");
		return 1;
	}

	return 0;

}

matrix mapping::dr(const matrix &a) const {

	return (D,a)/rz;
}

matrix mapping::dt(const matrix &a) const {

	return (a,Dt)-(D,a)*rt/rz;
}

matrix mapping::dt_odd(const matrix &a) const {

	return (a,Dt_11)-(D,a)*rt/rz;
}

matrix mapping::dt2(const matrix &a) const {

	return dt_odd(dt(a));
}

matrix mapping::stream(const matrix &Fz,matrix &Ft) const {

// Calculates the stream function of the divergenceless vector field whose
// zeta-contravariant component is Fz (for Fz even at pole and equator).
// Returns also the theta-contravariant component Ft
	matrix f;

	f=r*r*rz*Fz;
	f=(f,leg.P_00);
	f=-f/leg.l_00()/(leg.l_00()+1.);
	f.setcol(0,zeros(gl.N,1));
	f=(f,leg.dP1_00);
	f/=r;
	f.setrow(0,zeros(1,leg.npts));

	Ft=-(D,f*r)/r/r/rz;
	Ft.setrow(0,zeros(1,leg.npts));

	return f;


}

matrix mapping::stream(const matrix &Fz) const {

	matrix Ft;
	return stream(Fz,Ft);
}

matrix mapping::eval(const matrix &y,const matrix &ri, const matrix &thi,int parity) const {

	if(ri.nrows()!=thi.nrows()||ri.ncols()!=thi.ncols()) {
		ester_err("(mapping.eval) Matrix dimensions must agree");
		exit(1);
	}
	if(y.nrows()!=gl.N||y.ncols()!=leg.npts) {
		ester_err("(mapping.eval) Matrix dimensions must agree");
		exit(1);
	}

	matrix yi(ri.nrows(),ri.ncols());
	int N=ri.nrows()*ri.ncols();
	matrix yth,T,rth,rzth;
	double zi;
	for(int i=0;i<N;i++) {
		yth=leg.eval(y,thi(i),T,parity/10,parity%10);
		rth=leg.eval_00(r,thi(i),T);
		rzth=(rz,T);
		zi=ri(i)/max(rth);
		if(zi-1>-1e-10&&zi-1<1e-10) zi=1;
		if(zi<1e-10&&zi>-1e-10) zi=0;
		if(zi>1||zi<0) {
			ester_err("(mapping.eval) Coordinates (r,theta)=(%f,%f) are out of bounds",
				ri(i),thi(i));
			exit(1);
		}
		int fin=0,nit=0;
		if(zi==0||zi==1) fin=99;
		while(fin<2) {
			double dzi;
			dzi=-(gl.eval(rth,zi)(0)-ri(i))/gl.eval(rzth,zi)(0);
			if(fabs(dzi)<1e-10) fin++;
			nit++;
			if(nit>100) {
				ester_err("(mapping.eval) Failed to converge");
				exit(1);
			}
			zi+=dzi;
		}
		yi(i)=gl.eval(yth,zi)(0);
	}

	return yi;

}

// compute r for all zeta values
matrix mapping::zeta_to_r(const matrix &z) const {
	matrix rr(z.nrows(),z.ncols());

	if(z.ncols()!=leg.npts) {
		ester_err("(mapping.zeta_to_r) Matrix must have nth columns");
		exit(1);
	}

	for(int j=0;j<leg.npts;j++) {
		rr.setcol(j,gl.eval(r.col(j),z.col(j)));
	}
	return rr;
}

matrix mapping::draw0(const matrix &A,int parity,matrix &x,matrix &y) const {

	matrix zz(gl.N,4*leg.npts+1);

	x.dim(gl.N,4*leg.npts+1);
	y.dim(gl.N,4*leg.npts+1);

	x.setblock(0,gl.N-1,0,leg.npts-1,r*sin(th));
	y.setblock(0,gl.N-1,0,leg.npts-1,r*cos(th));
	zz.setblock(0,gl.N-1,0,leg.npts-1,A);

	x.setblock(0,gl.N-1,leg.npts,2*leg.npts-1,-x.block(0,gl.N-1,0,leg.npts-1).fliplr());
	y.setblock(0,gl.N-1,leg.npts,2*leg.npts-1,y.block(0,gl.N-1,0,leg.npts-1).fliplr());
	switch(parity) {
		case 00:
			zz.setblock(0,gl.N-1,leg.npts,2*leg.npts-1,zz.block(0,gl.N-1,0,leg.npts-1).fliplr());
			break;
		case 01:
			zz.setblock(0,gl.N-1,leg.npts,2*leg.npts-1,zz.block(0,gl.N-1,0,leg.npts-1).fliplr());
			break;
		case 10:
			zz.setblock(0,gl.N-1,leg.npts,2*leg.npts-1,-zz.block(0,gl.N-1,0,leg.npts-1).fliplr());
			break;
		case 11:
			zz.setblock(0,gl.N-1,leg.npts,2*leg.npts-1,-zz.block(0,gl.N-1,0,leg.npts-1).fliplr());
	}

	x.setblock(0,gl.N-1,2*leg.npts,3*leg.npts-1,-x.block(0,gl.N-1,0,leg.npts-1));
	y.setblock(0,gl.N-1,2*leg.npts,3*leg.npts-1,-y.block(0,gl.N-1,0,leg.npts-1));
	switch(parity) {
		case 00:
			zz.setblock(0,gl.N-1,2*leg.npts,3*leg.npts-1,zz.block(0,gl.N-1,0,leg.npts-1));
			break;
		case 01:
			zz.setblock(0,gl.N-1,2*leg.npts,3*leg.npts-1,-zz.block(0,gl.N-1,0,leg.npts-1));
			break;
		case 10:
			zz.setblock(0,gl.N-1,2*leg.npts,3*leg.npts-1,-zz.block(0,gl.N-1,0,leg.npts-1));
			break;
		case 11:
			zz.setblock(0,gl.N-1,2*leg.npts,3*leg.npts-1,zz.block(0,gl.N-1,0,leg.npts-1));
	}

	x.setblock(0,gl.N-1,3*leg.npts,4*leg.npts-1,x.block(0,gl.N-1,0,leg.npts-1).fliplr());
	y.setblock(0,gl.N-1,3*leg.npts,4*leg.npts-1,-y.block(0,gl.N-1,0,leg.npts-1).fliplr());
	switch(parity) {
		case 00:
			zz.setblock(0,gl.N-1,3*leg.npts,4*leg.npts-1,zz.block(0,gl.N-1,0,leg.npts-1).fliplr());
			break;
		case 01:
			zz.setblock(0,gl.N-1,3*leg.npts,4*leg.npts-1,-zz.block(0,gl.N-1,0,leg.npts-1).fliplr());
			break;
		case 10:
			zz.setblock(0,gl.N-1,3*leg.npts,4*leg.npts-1,zz.block(0,gl.N-1,0,leg.npts-1).fliplr());
			break;
		case 11:
			zz.setblock(0,gl.N-1,3*leg.npts,4*leg.npts-1,-zz.block(0,gl.N-1,0,leg.npts-1).fliplr());
	}
	x.setcol(4*leg.npts,x.col(0));
	y.setcol(4*leg.npts,y.col(0));
	zz.setcol(4*leg.npts,zz.col(0));

	return zz;
}


#if 0
matrix mapping::drawi0(const matrix &A,int sr,int st,int parity,matrix &x,matrix &y) const {

	matrix zz(sr,4*st+1),r1,z1;
	diff_leg leg1;
	diff_gl gl1;

	leg1.npts=st;
	leg1.init();
	gl1.set_ndomains(1);
	gl1.set_npts(sr);
	gl1.set_xif(gl.xif[0],gl.xif[gl.ndomains]);
	gl1.init();

	r1=gl.eval(r,gl1.x);
	r1=leg.eval_00(r1,leg1.th);
	z1=gl.eval(A,gl1.x);
	switch(parity) {
		case 00:
			z1=leg.eval_00(z1,leg1.th);
			break;
		case 01:
			z1=leg.eval_01(z1,leg1.th);
			break;
		case 10:
			z1=leg.eval_10(z1,leg1.th);
			break;
		case 11:
			z1=leg.eval_11(z1,leg1.th);
	}

	x.dim(sr,4*st+1);
	y.dim(sr,4*st+1);

	x.setblock(0,sr-1,0,st-1,r1*sin(leg1.th));
	y.setblock(0,sr-1,0,st-1,r1*cos(leg1.th));
	zz.setblock(0,sr-1,0,st-1,z1);

	x.setblock(0,sr-1,st,2*st-1,-x.block(0,sr-1,0,st-1).fliplr());
	y.setblock(0,sr-1,st,2*st-1,y.block(0,sr-1,0,st-1).fliplr());
	switch(parity) {
		case 00:
			zz.setblock(0,sr-1,st,2*st-1,zz.block(0,sr-1,0,st-1).fliplr());
			break;
		case 01:
			zz.setblock(0,sr-1,st,2*st-1,zz.block(0,sr-1,0,st-1).fliplr());
			break;
		case 10:
			zz.setblock(0,sr-1,st,2*st-1,-zz.block(0,sr-1,0,st-1).fliplr());
			break;
		case 11:
			zz.setblock(0,sr-1,st,2*st-1,-zz.block(0,sr-1,0,st-1).fliplr());
	}

	x.setblock(0,sr-1,2*st,3*st-1,-x.block(0,sr-1,0,st-1));
	y.setblock(0,sr-1,2*st,3*st-1,-y.block(0,sr-1,0,st-1));
	switch(parity) {
		case 00:
			zz.setblock(0,sr-1,2*st,3*st-1,zz.block(0,sr-1,0,st-1));
			break;
		case 01:
			zz.setblock(0,sr-1,2*st,3*st-1,-zz.block(0,sr-1,0,st-1));
			break;
		case 10:
			zz.setblock(0,sr-1,2*st,3*st-1,-zz.block(0,sr-1,0,st-1));
			break;
		case 11:
			zz.setblock(0,sr-1,2*st,3*st-1,zz.block(0,sr-1,0,st-1));
	}

	x.setblock(0,sr-1,3*st,4*st-1,x.block(0,sr-1,0,st-1).fliplr());
	y.setblock(0,sr-1,3*st,4*st-1,-y.block(0,sr-1,0,st-1).fliplr());
	switch(parity) {
		case 00:
			zz.setblock(0,sr-1,3*st,4*st-1,zz.block(0,sr-1,0,st-1).fliplr());
			break;
		case 01:
			zz.setblock(0,sr-1,3*st,4*st-1,-zz.block(0,sr-1,0,st-1).fliplr());
			break;
		case 10:
			zz.setblock(0,sr-1,3*st,4*st-1,zz.block(0,sr-1,0,st-1).fliplr());
			break;
		case 11:
			zz.setblock(0,sr-1,3*st,4*st-1,-zz.block(0,sr-1,0,st-1).fliplr());
	}

	x.setcol(4*st,x.col(0));
	y.setcol(4*st,y.col(0));
	zz.setcol(4*st,zz.col(0));

	return zz;
}

void mapping::draw(figure *pfig,const matrix &A,int parity) const {

	matrix x,y,zz;

	pfig->axis(-max(r.row(-1))*1.02,max(r.row(-1))*1.02,-max(r.row(-1))*1.02,max(r.row(-1))*1.02,1);
	zz=draw0(A,parity,x,y);

	pfig->pcolor(x,y,zz);

}

void mapping::drawi(figure *pfig,const matrix &A,int sr,int st,int parity) const {

	matrix x,y,zz;

	pfig->axis(-max(r.row(-1))*1.02,max(r.row(-1))*1.02,-max(r.row(-1))*1.02,max(r.row(-1))*1.02,1);
	zz=drawi0(A,sr,st,parity,x,y);

	pfig->pcolor(x,y,zz);

}

void mapping::drawc(figure *pfig,const matrix &A,int ncontours,int parity) const {

	matrix x,y,zz;

	pfig->axis(-max(r.row(-1))*1.02,max(r.row(-1))*1.02,-max(r.row(-1))*1.02,max(r.row(-1))*1.02,1);
	zz=draw0(A,parity,x,y);

	pfig->plot(x.row(-1),y.row(-1));
	pfig->hold(1);

	matrix cc;
	double zmin,zmax;
	double cmin,cmax;
	zmin=min(zz);zmax=max(zz);
	cmin=0;
	if(zmin>0) cmin=zmin;
	else if(zmax<0) cmin=-zmax;
	cmax=zmax>-zmin?zmax:-zmin;
	pfig->caxis(cmin,cmax);

	if(min(zz)<0) {
		cc=vector(zmin,zmax>0?0:zmax,ncontours);
		pfig->contour(x,y,-zz,-cc,"k=");
	}
	if(max(zz)>0) {
		cc=vector(zmin<0?0:zmin,zmax,ncontours);
		pfig->contour(x,y,zz,cc,"k-");
	}
	pfig->hold(0);

}

void mapping::drawci(figure *pfig,const matrix &A,int sr,int st,int ncontours,int parity) const {

	matrix x,y,zz;


	pfig->axis(-max(r.row(-1))*1.02,max(r.row(-1))*1.02,-max(r.row(-1))*1.02,max(r.row(-1))*1.02,1);
	zz=drawi0(A,sr,st,parity,x,y);

	pfig->plot(x.row(-1),y.row(-1));
	pfig->hold(1);

	matrix cc;
	double zmin,zmax;
	double cmin,cmax;
	zmin=min(zz);zmax=max(zz);
	cmin=0;
	if(zmin>0) cmin=zmin;
	else if(zmax<0) cmin=-zmax;
	cmax=zmax>-zmin?zmax:-zmin;
	pfig->caxis(cmin,cmax);

	if(min(zz)<0) {
		cc=vector(zmin,zmax>0?0:zmax,ncontours);
		pfig->contour(x,y,-zz,-cc,"k=");
	}
	if(max(zz)>0) {
		cc=vector(zmin<0?0:zmin,zmax,ncontours);
		pfig->contour(x,y,zz,cc,"k-");
	}
	pfig->hold(0);


}

void mapping::spectrum(figure *pfig,const matrix &y,int parity) const {

	matrix ys(y.nrows(),y.ncols()),xv,Pleg,l;
	int i,j;

	switch(parity) {
		case 00:
			Pleg=leg.P_00;
			l=leg.l_00();
			break;
		case 01:
			Pleg=leg.P_01;
			l=leg.l_01();
			break;
		case 10:
			Pleg=leg.P_10;
			l=leg.l_10();
			break;
		case 11:
			Pleg=leg.P_11;
			l=leg.l_11();
	}

	ys=abs((gl.P,y,Pleg));
	j=0;
	for(i=0;i<gl.ndomains;i++) {
		ys.setblock(j,j+gl.npts[i]-1,0,leg.npts-1,
			ys.block(j,j+gl.npts[i]-1,0,leg.npts-1)/max(ys.block(j,j+gl.npts[i]-1,0,leg.npts-1)));
		j+=gl.npts[i];
	}
	xv=vector_t(1,ys.nrows(),ys.nrows());
	pfig->caxis(-16,0);
	pfig->pcolor(l,xv,log10(abs(ys)));

	pfig->hold(1);
	j=gl.npts[0];
	for(i=1;i<gl.ndomains;i++) {
		pfig->plot(l,(j+1)*ones(1,leg.npts),":");
		j+=gl.npts[i];
	}
	pfig->hold(0);

}
#endif

mapping::ext_map::ext_map(mapping *map) {

	parent=map;

}

void mapping::ext_map::copy(const ext_map &map) {

	gl=map.gl;
	z=map.z;D=map.D;
	r=map.r;rz=map.rz;rzz=map.rzz;rt=map.rt;rtt=map.rtt;rzt=map.rzt;
	gzz=map.gzz;gzt=map.gzt;gtt=map.gtt;
	J[0]=map.J[0];J[1]=map.J[1];J[2]=map.J[2];J[3]=map.J[3];

}

mapping::ext_map::operator mapping() {

	mapping map;

	map.eta_=zeros(2,1);map.eta_(0)=parent->eta_(-1);map.eta_(1)=INFINITY;
    map.gl=gl;
    map.gl.D.block(0)=D;
    map.leg=parent->leg;
    map.r=r;map.rz=rz;map.rzz=rzz;map.rt=rt;map.rtt=rtt;map.rzt=rzt;
    map.gzz=gzz;map.gzt=gzt;map.gtt=gtt;
    map.R=zeros(2,parent->nt);map.R.setrow(0,(parent->R).row(-1));map.R.setrow(1,map.R.row(1)+INFINITY);
    map.mode=parent->mode;
	map.J[0]=J[0];map.J[1]=J[1];map.J[2]=J[2];map.J[3]=J[3];

	return map;
}





