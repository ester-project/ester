#include"mapping.h"
#include"constants.h"
#include<stdlib.h>

mapping::mapping(int ndom): gl(ndom),D(gl.D),Dt(leg.D_00),Dt2(leg.D2_00),Dt_odd(leg.D_11),z(gl.x),th(leg.th) {

	ex.gl.set_ndomains(1);
	ex.gl.set_xif(1.,2.);
	mode=0;
	
}

mapping::~mapping() {}

mapping::mapping(const mapping &map):
		eps_(map.eps_),eta_(map.eta_),
		gl(map.gl),leg(map.leg),
		r(map.r),rz(map.rz),rzz(map.rzz),rt(map.rt),rtt(map.rtt),rzt(map.rzt),
		gzz(map.gzz),gzt(map.gzt),gtt(map.gtt),R(map.R),D(gl.D),Dt(leg.D_00),Dt2(leg.D2_00),Dt_odd(leg.D_11),z(gl.x),th(leg.th) {
	z=map.z;th=map.th;ex.z=map.ex.z;ex.D=map.ex.D;
	ex.gl=map.ex.gl;
	ex.r=map.ex.r;ex.rz=map.ex.rz;ex.rzz=map.ex.rzz;ex.rt=map.ex.rt;
	ex.rtt=map.ex.rtt;ex.rzt=map.ex.rzt;
	ex.gzz=map.ex.gzz;ex.gzt=map.ex.gzt;ex.gtt=map.ex.gtt;
	mode=map.mode;
	J[0]=map.J[0];J[1]=map.J[1];J[2]=map.J[2];J[3]=map.J[3];
	ex.J[0]=map.ex.J[0];ex.J[1]=map.ex.J[1];ex.J[2]=map.ex.J[2];ex.J[3]=map.ex.J[3];
}

mapping &mapping::operator=(const mapping &map) {

	eps_=map.eps_;
	eta_=map.eta_;
	ex.z=map.ex.z;ex.D=map.ex.D;
    gl=map.gl;
    leg=map.leg;
    r=map.r;rz=map.rz;rzz=map.rzz;rt=map.rt;rtt=map.rtt;rzt=map.rzt;
    gzz=map.gzz;gzt=map.gzt;gtt=map.gtt;
    R=map.R;
    ex.gl=map.ex.gl;
    ex.r=map.ex.r;ex.rz=map.ex.rz;ex.rzz=map.ex.rzz;ex.rt=map.ex.rt;ex.rtt=map.ex.rtt;ex.rzt=map.ex.rzt;
    ex.gzz=map.ex.gzz;ex.gzt=map.ex.gzt;ex.gtt=map.ex.gtt;
	mode=map.mode;
	J[0]=map.J[0];J[1]=map.J[1];J[2]=map.J[2];J[3]=map.J[3];
	ex.J[0]=map.ex.J[0];ex.J[1]=map.ex.J[1];ex.J[2]=map.ex.J[2];ex.J[3]=map.ex.J[3];

	return *this;
}

int mapping::init() {

	int i;
	matrix temp;
	
	gl.xif[0]=0;
	gl.xif[gl.ndomains()]=1;
	R.dim(gl.ndomains(),leg.npts);
	temp=ones(1,leg.npts);
	for(i=0;i<gl.ndomains();i++)
		R.setrow(i,gl.xif[i+1]*temp);
	leg.init();
	ex.gl.init();
	return remap();
	
}

int mapping::remap(const diff_leg &leg_new) {

	R=leg.eval_00(R,leg_new.th);
	leg=leg_new;
	
	return remap();
	
}

double mapping::eps(int i) const {

	return eps_(i);

}

double mapping::eta(int i) const {

	return eta_(i);

}

matrix mapping::eps() const {

	return eps_;
	
}

matrix mapping::eta() const{

	return eta_;
	
}

int mapping::remap() {
	
	matrix zz,xi,A,Ap,App,temp;
	int jfirst,i,dj,iprev;
	double a;
	
	eta_=leg.eval_00(R,0);
	R.setrow(gl.ndomains()-1,R.row(gl.ndomains()-1)/eta_(gl.ndomains()-1));
	eta_(gl.ndomains()-1)=1;
	eps_=1-eta_/leg.eval_00(R,PI/2);
	gl.xif[0]=0;
	for(i=0;i<gl.ndomains();i++)
		gl.xif[i+1]=eta_(i);
	gl.init();

	r.dim(gl.N(),leg.npts);
	rz.dim(gl.N(),leg.npts);
	rzz.dim(gl.N(),leg.npts);
	rt.dim(gl.N(),leg.npts);
	rtt.dim(gl.N(),leg.npts);
	rzt.dim(gl.N(),leg.npts);
	J[0].dim(gl.N(),leg.npts);
	J[1].dim(gl.N(),leg.npts);
	J[2].dim(gl.N(),leg.npts);
	J[3].dim(gl.N(),leg.npts);
	
	jfirst=0; 
	i=0;
	iprev=0;
	while (i<gl.ndomains()) {
		dj=gl.npts[i];
		zz=gl.x.block(jfirst,jfirst+dj-1,0,0);
		if(!jfirst) {
			xi=zz/eta_(i);
			if(!mode) {
				A=(5*xi*xi*xi-3*pow(xi,5))/2;
				Ap=7.5*(xi*xi-pow(xi,4));
				App=15*(xi-2*xi*xi*xi);
			} else {
				A=6*pow(xi,5)-15*pow(xi,4)+10*xi*xi*xi;
				Ap=30*(pow(xi,4)-2*xi*xi*xi+xi*xi);
				App=60*(2*xi*xi*xi-3*xi*xi+xi);
			}
			a=eps_(i)<0 ? 1/(1-eps_(i)) : 1;
			//if(a!=1) printf("Mapping: a(%d)=%e\n",i,a);
			a=1;
			temp=a*zz+A*(R.row(i)-a*eta_(i));
			r.setblock(0,dj-1,0,leg.npts-1,temp);
			temp=a+Ap*(R.row(i)/eta_(i)-a);
			rz.setblock(0,dj-1,0,leg.npts-1,temp);
			temp=App/eta_(i)*(R.row(i)/eta_(i)-a);
			rzz.setblock(0,dj-1,0,leg.npts-1,temp);
			temp=A*(R.row(i),Dt);
			rt.setblock(0,dj-1,0,leg.npts-1,temp);
			temp=A*(R.row(i),Dt2);
			rtt.setblock(0,dj-1,0,leg.npts-1,temp);
			temp=Ap/eta_(i)*(R.row(i),Dt);
			rzt.setblock(0,dj-1,0,leg.npts-1,temp);
			temp=zeros(dj,leg.npts);
			J[0].setblock(0,dj-1,0,leg.npts-1,temp);
			temp=a*(xi-A)*ones(1,leg.npts);
			J[1].setblock(0,dj-1,0,leg.npts-1,temp);
			temp=zeros(dj,leg.npts);
			J[2].setblock(0,dj-1,0,leg.npts-1,temp);
			temp=A*ones(1,leg.npts);
			J[3].setblock(0,dj-1,0,leg.npts-1,temp);
		} else {
			xi=(2*zz-eta_(i)-eta_(iprev))/(eta_(i)-eta(iprev));
			if(!mode) {
				A=0.25*(xi*xi*xi-3*xi+2);
				Ap=0.25*(3*xi*xi-3);
				App=1.5*xi;
			} else {
				A=(-3*pow(xi,5)+10*xi*xi*xi-15*xi+8)/16;
				Ap=15*(-pow(xi,4)+2*xi*xi-1)/16;
				App=15*(-xi*xi*xi+xi)/4;
			}
			a=eta_(i)*eps_(i)<eta_(iprev)*eps_(iprev) ?
			(1-(eta_(i)*eps_(iprev)-eta_(iprev)*eps_(i))/(eta_(i)-eta_(iprev)))/(1-eps_(i))/(1-eps_(iprev))
			 : 1;
			 //if(a!=1) printf("Mapping: a(%d)=%e\n",i,a);
			 a=1;
			temp=a*zz+A*(R.row(iprev)-a*eta_(iprev))+(1-A)*(R.row(i)-a*eta_(i));
			r.setblock(jfirst,jfirst+dj-1,0,leg.npts-1,temp);
			temp=a-2*Ap*((R.row(i)-R.row(iprev))/(eta_(i)-eta_(iprev))-a);
			rz.setblock(jfirst,jfirst+dj-1,0,leg.npts-1,temp);
			temp=-4*App/(eta_(i)-eta_(iprev))*((R.row(i)-R.row(iprev))/(eta_(i)-eta_(iprev))-a);
			rzz.setblock(jfirst,jfirst+dj-1,0,leg.npts-1,temp);
			temp=A*(R.row(iprev),Dt)+(1-A)*(R.row(i),Dt);
			rt.setblock(jfirst,jfirst+dj-1,0,leg.npts-1,temp);
			temp=A*(R.row(iprev),Dt2)+(1-A)*(R.row(i),Dt2);
			rtt.setblock(jfirst,jfirst+dj-1,0,leg.npts-1,temp);
			temp=-2*Ap*(R.row(i)-R.row(iprev),Dt)/(eta_(i)-eta_(iprev));
			rzt.setblock(jfirst,jfirst+dj-1,0,leg.npts-1,temp);
			temp=zeros(dj,leg.npts);
			J[0].setblock(jfirst,jfirst+dj-1,0,leg.npts-1,temp);
			temp=a*((1.+xi)/2.-(1-A))*ones(1,leg.npts);
			J[1].setblock(jfirst,jfirst+dj-1,0,leg.npts-1,temp);
			temp=ones(dj,leg.npts);
			J[2].setblock(jfirst,jfirst+dj-1,0,leg.npts-1,temp);
			temp=(1-A)*ones(1,leg.npts);
			J[3].setblock(jfirst,jfirst+dj-1,0,leg.npts-1,temp);
		}
		iprev=i;
		jfirst+=dj;
		i++;
	}

	jfirst=0;
	for(i=0;i<gl.ndomains();i++) {
		R.setrow(i,r.row(jfirst+gl.npts[i]-1));
		jfirst+=gl.npts[i];
	}
	eps_=1-eta_/leg.eval_00(R,PI/2);
	
	gzz=(r*r+rt*rt)/r/r/rz/rz;
	gzz.setrow(0,1/rz.row(0)/rz.row(0));
	gzt=-rt/r/r/rz;
	gzt.setrow(0,zeros(1,leg.npts));
	gtt=1/r/r;

	zz=1/(2-ex.gl.x);
	i=gl.ndomains()-1;
	a=1;
	ex.r=a*zz+R.row(i)-a;
	ex.rz=a*ones(ex.gl.N(),leg.npts);
	ex.rzz=zeros(ex.gl.N(),leg.npts);
	ex.rt=ones(ex.gl.N(),leg.npts)*(R.row(i),Dt);
	ex.rtt=ones(ex.gl.N(),leg.npts)*(R.row(i),Dt2);
	ex.rzt=zeros(ex.gl.N(),leg.npts);
	
	ex.J[0]=zeros(ex.gl.N(),leg.npts);
	ex.J[1]=zeros(ex.gl.N(),leg.npts);
	ex.J[2]=ones(ex.gl.N(),leg.npts);
	ex.J[3]=zeros(ex.gl.N(),leg.npts);
	
	
	ex.gzz=(1+ex.rt*ex.rt/ex.r/ex.r)/ex.rz/ex.rz;
	ex.gzt=-ex.rt/ex.r/ex.r/ex.rz;
	ex.gtt=1/ex.r/ex.r;

	ex.z=zz;
	ex.D=ex.gl.D/ex.z/ex.z;
	
	
	if(exist(rz<0)||exist(ex.rz<0)) return 1;
	
	return 0;
	
}

matrix mapping::dr(const matrix &a) const {

	return (gl.D,a)/rz;
}

matrix mapping::dt(const matrix &a) const {

	return (a,Dt)-(D,a)*rt/rz;	
}
	
matrix mapping::dt_odd(const matrix &a) const {

	return (a,Dt_odd)-(D,a)*rt/rz;	
}

matrix mapping::dt2(const matrix &a) const {

	return dt_odd(dt(a));	
}

matrix mapping::lap(const matrix &phi) const {

	matrix lap_phi,q;
	
	lap_phi=gzz*(D,D,phi);
	q=2*(1+rt*rzt/r/rz)/r/rz-gzz*rzz/rz-(rtt+rt*cos(th)/sin(th))/r/r/rz;
	lap_phi+=q*(D,phi);
	q=-2*rt/r/r/rz;
	lap_phi+=q*(D,phi,Dt);
	lap_phi+=(phi,leg.lap_00)/r/r;
	
	return lap_phi;

}

matrix mapping::lap_ex(const matrix &phi) const {

	matrix lap_phi,q;
	
	lap_phi=ex.gzz*(ex.D,ex.D,phi);
	q=2*(1+ex.rt*ex.rzt/ex.r/ex.rz)/ex.r/ex.rz-ex.gzz*ex.rzz/ex.rz-(ex.rtt+ex.rt*cos(th)/sin(th))/ex.r/ex.r/ex.rz;
	lap_phi+=q*(ex.D,phi);
	q=-2*ex.rt/ex.r/ex.r/ex.rz;
	lap_phi+=q*(ex.D,phi,Dt);
	lap_phi+=(phi,leg.lap_00)/ex.r/ex.r;
	
	return lap_phi;

}

void mapping::add_lap(solver *op,const char* eqn,const char * varn,const matrix &d,const matrix &phi) const {
	
	matrix q;
	
	op->add_l(eqn,varn,d*gzz,(D,D));
	q=2*(1+rt*rzt/r/rz)/r/rz-gzz*rzz/rz-(rtt+rt*cos(th)/sin(th))/r/r/rz;
	op->add_l(eqn,varn,d*q,D);
	q=-2*rt/r/r/rz;
	op->add_lr(eqn,varn,d*q,D,Dt);
	op->add_r(eqn,varn,d/r/r,leg.lap_00);
	
	q=-2*rt*rt/r/r/r/rz/rz*(D,D,phi)
		+(-2/r/r/rz-4*rt*rzt/r/r/r/rz/rz+2*rzz*rt*rt/r/r/r/rz/rz/rz+2*rtt/r/r/r/rz+2*rt*cos(th)/r/r/r/rz/sin(th))*(D,phi)
		+4*rt/r/r/r/rz*(D,phi,Dt)
		-2/r/r/r*(phi,leg.lap_00);
	op->add_d(eqn,"r",d*q);
	q=(-2/rz/rz/rz-2*rt*rt/r/r/rz/rz/rz)*(D,D,phi)
		+(-2/r/rz/rz-4*rt*rzt/r/r/rz/rz/rz+3*rzz/rz/rz/rz/rz+3*rzz*rt*rt/r/r/rz/rz/rz/rz+rtt/r/r/rz/rz+rt*cos(th)/r/r/rz/rz/sin(th))*(D,phi)
		+2*rt/r/r/rz/rz*(D,phi,Dt);
	op->add_l(eqn,"r",d*q,D);
	q=2*rt/r/r/rz/rz*(D,D,phi)
		+(2*rzt/r/r/rz/rz-2*rzz*rt/r/r/rz/rz/rz-cos(th)/r/r/rz/sin(th))*(D,phi)
		-2/r/r/rz*(D,phi,Dt);
	op->add_r(eqn,"r",d*q,Dt);
	q=(-1/rz/rz/rz-rt*rt/r/r/rz/rz/rz)*(D,phi);
	op->add_l(eqn,"r",d*q,(D,D));
	q=2*rt/r/r/rz/rz*(D,phi);
	op->add_lr(eqn,"r",d*q,D,Dt);
	q=-1/r/r/rz*(D,phi);	
	op->add_r(eqn,"r",d*q,Dt2);
	
}
void mapping::add_lap_ex(solver *op,const char* eqn,const char * varn,const matrix &d,const matrix &phi) const {
	
	matrix q;
	
	op->add_l(gl.ndomains(),eqn,varn,d*ex.gzz,(ex.D,ex.D));
	q=2*(1+ex.rt*ex.rzt/ex.r/ex.rz)/ex.r/ex.rz-ex.gzz*ex.rzz/ex.rz-(ex.rtt+ex.rt*cos(th)/sin(th))/ex.r/ex.r/ex.rz;
	op->add_l(gl.ndomains(),eqn,varn,d*q,ex.D);
	q=-2*ex.rt/ex.r/ex.r/ex.rz;
	op->add_lr(gl.ndomains(),eqn,varn,d*q,ex.D,Dt);
	op->add_r(gl.ndomains(),eqn,varn,d/ex.r/ex.r,leg.lap_00);
	
	q=-2*ex.rt*ex.rt/ex.r/ex.r/ex.r/ex.rz/ex.rz*(ex.D,ex.D,phi)
		+(-2/ex.r/ex.r/ex.rz-4*ex.rt*ex.rzt/ex.r/ex.r/ex.r/ex.rz/ex.rz+2*ex.rzz*ex.rt*ex.rt/ex.r/ex.r/ex.r/ex.rz/ex.rz/ex.rz+2*ex.rtt/ex.r/ex.r/ex.r/ex.rz+2*ex.rt*cos(th)/ex.r/ex.r/ex.r/ex.rz/sin(th))*(ex.D,phi)
		+4*ex.rt/ex.r/ex.r/ex.r/ex.rz*(ex.D,phi,Dt)
		-2/ex.r/ex.r/ex.r*(phi,Dt2)
		-2*cos(th)/ex.r/ex.r/ex.r/sin(th)*(phi,Dt);
	op->add_d(gl.ndomains(),eqn,"r",d*q);	
	q=(-2/ex.rz/ex.rz/ex.rz-2*ex.rt*ex.rt/ex.r/ex.r/ex.rz/ex.rz/ex.rz)*(ex.D,ex.D,phi)
		+(-2/ex.r/ex.rz/ex.rz-4*ex.rt*ex.rzt/ex.r/ex.r/ex.rz/ex.rz/ex.rz+3*ex.rzz/ex.rz/ex.rz/ex.rz/ex.rz+3*ex.rzz*ex.rt*ex.rt/ex.r/ex.r/ex.rz/ex.rz/ex.rz/ex.rz+ex.rtt/ex.r/ex.r/ex.rz/ex.rz+ex.rt*cos(th)/ex.r/ex.r/ex.rz/ex.rz/sin(th))*(ex.D,phi)
		+2*ex.rt/ex.r/ex.r/ex.rz/ex.rz*(ex.D,phi,Dt);
	op->add_l(gl.ndomains(),eqn,"r",d*q,ex.D);
	q=2*ex.rt/ex.r/ex.r/ex.rz/ex.rz*(ex.D,ex.D,phi)
		+(2*ex.rzt/ex.r/ex.r/ex.rz/ex.rz-2*ex.rzz*ex.rt/ex.r/ex.r/ex.rz/ex.rz/ex.rz-cos(th)/ex.r/ex.r/ex.rz/sin(th))*(ex.D,phi)
		-2/ex.r/ex.r/ex.rz*(ex.D,phi,Dt);
	op->add_r(gl.ndomains(),eqn,"r",d*q,Dt);
	q=(-1/ex.rz/ex.rz/ex.rz-ex.rt*ex.rt/ex.r/ex.r/ex.rz/ex.rz/ex.rz)*(ex.D,phi);
	op->add_l(gl.ndomains(),eqn,"r",d*q,(ex.D,ex.D));
	q=2*ex.rt/ex.r/ex.r/ex.rz/ex.rz*(ex.D,phi);
	op->add_lr(gl.ndomains(),eqn,"r",d*q,ex.D,Dt);
	q=-1/ex.r/ex.r/ex.rz*(ex.D,phi);	
	op->add_r(gl.ndomains(),eqn,"r",d*q,Dt2);

}

matrix mapping::stream(const matrix &Fz,matrix &Ft) const {

	// Calculates the stream function of the divergenceless vector field whose
	// zeta-contravariant coordinate is Fz (for Fz even at pole and equator).
	// Returns also the theta-contravariant coordinate Ft
	
	matrix f;
	
	f=r*r*rz*Fz;
	f=(f,leg.P_00);
	f=-f/leg.l_00()/(leg.l_00()+1.);
	f.setcol(0,zeros(gl.N(),1));
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
    
void mapping::interps(const mapping &map_old,matrix &Tr,matrix &Tex,
		matrix &Tt_00,matrix &Tt_01,matrix &Tt_10,matrix &Tt_11) const {

	matrix m;
	int i,j,k;
	
	if(gl.ndomains()!=map_old.gl.ndomains()) {
		printf("mapping::interps : Changing number of domains not yet implemented\n");
		exit(1);
	}
	for(i=0;i<gl.ndomains();i++) {
		gl.xif[i+1]=map_old.gl.xif[i+1];
	}
	
	m=zeros(map_old.gl.N(),1);
	m=map_old.gl.eval(m,gl.x,Tr);
	j=0;k=0;
	for(i=0;i<gl.ndomains();i++) {
		Tr.setrow(k,zeros(1,map_old.gl.N()));Tr(k,j)=1;
		j+=map_old.gl.npts[i];
		k+=gl.npts[i];
		Tr.setrow(k-1,zeros(1,map_old.gl.N()));Tr(k-1,j-1)=1;
	}
	
	m=zeros(map_old.ex.gl.N(),1);
	m=map_old.ex.gl.eval(m,ex.gl.x,Tex);

	m=zeros(1,map_old.leg.npts);
	m=map_old.leg.eval_00(m,leg.th,Tt_00);
	
	m=zeros(1,map_old.leg.npts);
	m=map_old.leg.eval_01(m,leg.th,Tt_01);
	
	m=zeros(1,map_old.leg.npts);
	m=map_old.leg.eval_10(m,leg.th,Tt_10);
	
	m=zeros(1,map_old.leg.npts);
	m=map_old.leg.eval_11(m,leg.th,Tt_11);

}

matrix mapping::draw0(const matrix &A,int parity,matrix &x,matrix &y) const {
	
	matrix zz(gl.N(),4*leg.npts+1);
	
	x.dim(gl.N(),4*leg.npts+1);
	y.dim(gl.N(),4*leg.npts+1);
	
	x.setblock(0,gl.N()-1,0,leg.npts-1,r*sin(th));
	y.setblock(0,gl.N()-1,0,leg.npts-1,r*cos(th));
	zz.setblock(0,gl.N()-1,0,leg.npts-1,A);
	
	x.setblock(0,gl.N()-1,leg.npts,2*leg.npts-1,-x.block(0,gl.N()-1,0,leg.npts-1).fliplr());
	y.setblock(0,gl.N()-1,leg.npts,2*leg.npts-1,y.block(0,gl.N()-1,0,leg.npts-1).fliplr());
	switch(parity) {
		case 00:
			zz.setblock(0,gl.N()-1,leg.npts,2*leg.npts-1,zz.block(0,gl.N()-1,0,leg.npts-1).fliplr());
			break;
		case 01:
			zz.setblock(0,gl.N()-1,leg.npts,2*leg.npts-1,zz.block(0,gl.N()-1,0,leg.npts-1).fliplr());
			break;
		case 10:
			zz.setblock(0,gl.N()-1,leg.npts,2*leg.npts-1,-zz.block(0,gl.N()-1,0,leg.npts-1).fliplr());
			break;
		case 11:
			zz.setblock(0,gl.N()-1,leg.npts,2*leg.npts-1,-zz.block(0,gl.N()-1,0,leg.npts-1).fliplr());
	}
	
	x.setblock(0,gl.N()-1,2*leg.npts,3*leg.npts-1,-x.block(0,gl.N()-1,0,leg.npts-1));
	y.setblock(0,gl.N()-1,2*leg.npts,3*leg.npts-1,-y.block(0,gl.N()-1,0,leg.npts-1));
	switch(parity) {
		case 00:
			zz.setblock(0,gl.N()-1,2*leg.npts,3*leg.npts-1,zz.block(0,gl.N()-1,0,leg.npts-1));
			break;
		case 01:
			zz.setblock(0,gl.N()-1,2*leg.npts,3*leg.npts-1,-zz.block(0,gl.N()-1,0,leg.npts-1));
			break;
		case 10:
			zz.setblock(0,gl.N()-1,2*leg.npts,3*leg.npts-1,-zz.block(0,gl.N()-1,0,leg.npts-1));
			break;
		case 11:
			zz.setblock(0,gl.N()-1,2*leg.npts,3*leg.npts-1,zz.block(0,gl.N()-1,0,leg.npts-1));
	}
	
	x.setblock(0,gl.N()-1,3*leg.npts,4*leg.npts-1,x.block(0,gl.N()-1,0,leg.npts-1).fliplr());
	y.setblock(0,gl.N()-1,3*leg.npts,4*leg.npts-1,-y.block(0,gl.N()-1,0,leg.npts-1).fliplr());
	switch(parity) {
		case 00:
			zz.setblock(0,gl.N()-1,3*leg.npts,4*leg.npts-1,zz.block(0,gl.N()-1,0,leg.npts-1).fliplr());
			break;
		case 01:
			zz.setblock(0,gl.N()-1,3*leg.npts,4*leg.npts-1,-zz.block(0,gl.N()-1,0,leg.npts-1).fliplr());
			break;
		case 10:
			zz.setblock(0,gl.N()-1,3*leg.npts,4*leg.npts-1,zz.block(0,gl.N()-1,0,leg.npts-1).fliplr());
			break;
		case 11:
			zz.setblock(0,gl.N()-1,3*leg.npts,4*leg.npts-1,-zz.block(0,gl.N()-1,0,leg.npts-1).fliplr());
	}
	x.setcol(4*leg.npts,x.col(0));
	y.setcol(4*leg.npts,y.col(0));
	zz.setcol(4*leg.npts,zz.col(0));
	
	return zz;
}


matrix mapping::drawi0(const matrix &A,int sr,int st,int parity,matrix &x,matrix &y) const {
	
	matrix zz(sr,4*st+1),r1,z1;
	diff_leg leg1;
	diff_gl gl1;
	
	leg1.npts=st;
	leg1.init();
	gl1.set_ndomains(1);
	gl1.set_npts(sr);
	gl1.set_xif(0.,1.);
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
	
	pfig->axis(-1/(1-eps(gl.ndomains()-1)),1/(1-eps(gl.ndomains()-1)),-1/(1-eps(gl.ndomains()-1)),1/(1-eps(gl.ndomains()-1)),1);
	zz=draw0(A,parity,x,y);
	
	pfig->pcolor(x,y,zz);
	
}

void mapping::drawi(figure *pfig,const matrix &A,int sr,int st,int parity) const {
	
	matrix x,y,zz;
	
	pfig->axis(-1/(1-eps(gl.ndomains()-1)),1/(1-eps(gl.ndomains()-1)),-1/(1-eps(gl.ndomains()-1)),1/(1-eps(gl.ndomains()-1)),1);
	zz=drawi0(A,sr,st,parity,x,y);
	
	pfig->pcolor(x,y,zz);
	
}

void mapping::drawc(figure *pfig,const matrix &A,int ncontours,int parity) const {
	
	matrix x,y,zz;
	
	pfig->axis(-1/(1-eps(gl.ndomains()-1)),1/(1-eps(gl.ndomains()-1)),-1/(1-eps(gl.ndomains()-1)),1/(1-eps(gl.ndomains()-1)),1);
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
	
	pfig->axis(-1/(1-eps(gl.ndomains()-1)),1/(1-eps(gl.ndomains()-1)),-1/(1-eps(gl.ndomains()-1)),1/(1-eps(gl.ndomains()-1)),1);
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
	for(i=0;i<gl.ndomains();i++) {
		ys.setblock(j,j+gl.npts[i]-1,0,leg.npts-1,
			ys.block(j,j+gl.npts[i]-1,0,leg.npts-1)/max(ys.block(j,j+gl.npts[i]-1,0,leg.npts-1)));
		j+=gl.npts[i];
	}
	xv=vector_t(1,ys.nrows(),ys.nrows());
	pfig->caxis(-16,0);
	pfig->pcolor(l,xv,log10(abs(ys)));
	
	pfig->hold(1);
	j=gl.npts[0];
	for(i=1;i<gl.ndomains();i++) {
		pfig->plot(l,(j+1)*ones(1,leg.npts),":");
		j+=gl.npts[i];
	}
	pfig->hold(0);

}


