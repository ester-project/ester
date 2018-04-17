#include "ester-config.h"
#include "star.h"

#if 0
void star2d::draw(figure *pfig,const matrix &A,int parity) const {
	
	map.draw(pfig,A,parity);
	
}

void star2d::drawi(figure *pfig,const matrix &A,int sr,int st,int parity) const {
	
	map.drawi(pfig,A,sr,st,parity);
	
}

void star2d::drawc(figure *pfig,const matrix &A,int ncontours,int parity) const {
	
	map.drawc(pfig,A,ncontours,parity);
	
}

void star2d::drawci(figure *pfig,const matrix &A,int sr,int st,int ncontours,int parity) const {
	
	map.drawci(pfig,A,sr,st,ncontours,parity);
	
}


void star2d::spectrum(figure *pfig,const matrix &y,int parity) const {

	map.spectrum(pfig,y,parity);
}
#endif

void star2d::calc_units() {

	units.phi=pc/rhoc;
	units.p=pc;
	units.rho=rhoc;
	units.T=Tc;
	units.r=R;
	units.Omega=sqrt(pc/rhoc)/R;
	units.v=sqrt(pc/rhoc);
	units.F=pc/R/rhoc;
}

double star2d::luminosity() const {

	return 2*PI*(map.gl.I,(rho*nuc.eps*r*r*map.rz,map.leg.I_00))(0)*units.rho*units.r*units.r*units.r;

}

double star2d::Lz() const {

	return 2*PI*(map.gl.I,(rho*w*r*r*sin(th)*sin(th)*r*r*map.rz,map.leg.I_00))(0)*units.rho*units.Omega*
		pow(units.r,5);

}


double star2d::Mcore() const {

	if(!conv) return 0;
	int jcc=0;
	for(int n=0;n<conv;n++) jcc+=map.gl.npts[n];
	jcc--;
	return 2*PI*(map.gl.I.block(0,0,0,jcc),
		(rho*r*r*map.rz).block(0,jcc,0,-1),map.leg.I_00)(0)*units.rho*units.r*units.r*units.r;
}


double star2d::Lzcore() const {

	if(!conv) return 0;
	int jcc=0;
	for(int n=0;n<conv;n++) jcc+=map.gl.npts[n];
	jcc--;
	return 2*PI*(map.gl.I.block(0,0,0,jcc),
		(rho*w*r*r*sin(th)*sin(th)*r*r*map.rz).block(0,jcc,0,-1),map.leg.I_00)(0)*units.rho*units.Omega*
		pow(units.r,5);
}

matrix star2d::Rcore() const {
	
	if(!conv) return zeros(1,nth);
	int jcc=0;
	for(int n=0;n<conv;n++) jcc+=map.gl.npts[n];
	jcc--;
	return r.row(jcc)*units.r;
}

matrix star2d::N2() const {

	matrix N2;

   	N2=-(map.gzz*(D,p)+map.gzt*(p,Dt))*((D,p)/p/eos.G1-(D,rho)/rho)
		-(map.gzt*(D,p)+map.gtt*(p,Dt))*((p,Dt)/p/eos.G1-(rho,Dt)/rho);
    N2=N2/rho;
    N2.setrow(0,zeros(1,nth));

    N2=N2*units.Omega*units.Omega;
    
    return N2;
    
}

matrix star2d::entropy() const {

//	Valid only for homogeneus composition !!!!!!

	matrix s(nr,nth),rhs;
	
	solver op;
	op.init(ndomains,1,"full");
	op.maxit_ref=10;op.use_cgs=0;
	op.rel_tol=1e-12;op.abs_tol=1e-20;
	op.regvar("s");op.set_nr(map.gl.npts);
	
	op.add_l("s","s",ones(nr,1),D);
	rhs=eos.cp*((D,log(T))-eos.del_ad*(D,log(p)))/RGP;
	
	op.bc_bot2_add_d(0,"s","s",ones(1,1));
	rhs.setrow(0,zeros(1,nth));
	int j0=map.gl.npts[0];
	for(int n=1;n<ndomains;n++) {
		op.bc_bot2_add_d(n,"s","s",ones(1,1));
		op.bc_bot1_add_d(n,"s","s",-ones(1,1));
		rhs.setrow(j0,zeros(1,nth));
		j0+=map.gl.npts[n];
	}
	for(int j=0;j<nth;j++) {
		op.set_rhs("s",rhs.col(j));
		op.solve();
		s.setcol(j,op.get_var("s"));
	}
	return s;
}

matrix star2d::Teff() const {

	matrix F;
	
	F=-opa.xi/sqrt(map.gzz)*(map.gzz*(D,T)+map.gzt*(T,Dt))/units.r*units.T;
	F=F.row(nr-1)/SIG_SB;
	
	return pow(F,0.25);

}

matrix star2d::gsup() const {

	matrix g;
	
	g=-(map.gzz*(D,p)+map.gzt*(p,Dt))/sqrt(map.gzz)/rho;
	g=g.row(nr-1)/units.r*units.p/units.rho;
	
	return g;

}

double star2d::virial_3P() const {

	return 6.*PI*(map.gl.I,(p*r*r*map.rz,map.leg.I_00))(0)*units.p*units.r*units.r*units.r;

}

double star2d::virial_W() const {

	return PI*(map.gl.I,(phi*rho*r*r*map.rz,map.leg.I_00))(0)*units.rho*units.phi*units.r*units.r*units.r;

}

double star2d::virial_L() const {

	return 2*PI*(map.gl.I,(rho*w*w*r*r*sin(th)*sin(th)*r*r*map.rz,map.leg.I_00))(0)*units.rho*units.phi*
		units.r*units.r*units.r;

}

double star2d::virial_ps() const {

	return -2*PI*(r.row(nr-1)*r.row(nr-1)*r.row(nr-1)*ps,map.leg.I_00)(0)*units.p*units.r*units.r*units.r;
	
}

double star2d::virial() const {

	return (virial_3P()+virial_L()+virial_ps())/virial_W()+1.;

}

double star2d::energy_test() const {

	double e1,e2;
	matrix Fz;
	
	e1=luminosity();
	Fz=-xic*opa.xi*(map.gzz*(D,T)+map.gzt*(T,Dt));
	e2=2*PI*((Fz*r*r*map.rz).row(nr-1),map.leg.I_00)(0)*units.T*units.r;
	
	return (e1-e2)/e1;

}

double star2d::apparent_luminosity(double i) const {

	/*	Apparent luminosity= 4*PI*d^2*FT
			d: distance
			FT: Flux measured on Earth
		FT=Integral_over_visible_disk( I*dW )
			dW=(i·n)*dS/d^2 differential solid angle
			I=(F·n)/PI specific intensity
	 
		Apparent luminosity = 4*Integral_over_visible_disk( (F·n)(i·n)dS )
		
		(F·n)=SIG_SB*Teff^4
		(i·n)= [ r*( cos(i)*cos(th)+sin(i)*sin(th)*cos(phi) )
				+rt*( cos(i)*sin(th)-sin(i)*cos(th)*cos(phi) ) ]
				/sqrt(r^2+rt^2)
		dS=sqrt(r^2+rt^2)*r*sin(th)*dth*dphi
		
		(i·n) is neither symmetric nor antisymmentric respect to the pole
		or the equator, and has no axial symmetry,
		so we have to use a 2D Gauss-Lobatto grid for the integral.
		Also, as we will set (i·n)=0 over the non-visible part,
		we will use a finer grid to reduce the impact of the discontinuity
		(in the first derivative) on the accuracy of the quadrature formula. */


	// theta grid
	diff_gl gl_th(1);
	gl_th.set_npts(500);
	gl_th.set_xif(0.,PI);
	gl_th.init();
	
	matrix th_f,Ith_f;
	th_f=gl_th.x.transpose();
	Ith_f=gl_th.I.transpose();
	
	// phi grid
	diff_gl gl_phi(1);
	gl_phi.set_npts(500);
	gl_phi.set_xif(0.,2*PI);
	gl_phi.init();
	
	matrix phi_f,Iphi_f;
	phi_f=gl_phi.x;
	Iphi_f=gl_phi.I;
	
	
	// Symmetric part of the integrand
	
	matrix int_s;
	
	int_s=SIG_SB*pow(Teff(),4);
	
	matrix r_f,rt_f;
	
	// Interpolating in the new grid
	
	int_s=map.leg.eval_00(int_s,th_f);
	r_f=map.leg.eval_00(r.row(nr-1),th_f);
	rt_f=map.leg.eval_11(map.rt.row(nr-1),th_f);
	
	// Non symmetric part of the integrand
	matrix int_ns;
	
	int_ns=r_f*( cos(i)*cos(th_f)+sin(i)*sin(th_f)*cos(phi_f) )
				+rt_f*( cos(i)*sin(th_f)-sin(i)*cos(th_f)*cos(phi_f) );
	int_ns*=(int_ns>0);
	
	// Integral
	double L_ap=4*(Iphi_f,int_s*int_ns*r_f*sin(th_f),Ith_f)(0);
	L_ap*=units.r*units.r;
	
	return L_ap;
	
}

matrix star2d::stream() const {

	matrix vr_,GG;
	
	vr_=rho*vr/map.rz;
	
	
	vr_=r*r*map.rz/z/z*vr_;
	vr_.setrow(0,zeros(1,nth));
	vr_=(vr_,map.leg.P_00);
	
	GG=-z*z*vr_/map.leg.l_00()/(map.leg.l_00()+1.);
	GG.setcol(0,zeros(nr,1));
	GG=(GG,map.leg.P1_00,Dt);
	GG=GG/r;
	GG.setrow(0,zeros(1,nth));
		
	return GG;
}






