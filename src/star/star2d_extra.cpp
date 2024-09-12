#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "star.h"
#include<stdlib.h> //JM

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
	units.v=R/MYR;
	units.F=pc/R/rhoc;
}

double star2d::luminosity() const {

	matrix Fz=-opa.xi*(map.gzz*(D,T)+map.gzt*(T,Dt));
	return 2*PI*((Fz*r*r*map.rz).row(nr-1),map.leg.I_00)(0)*units.T*units.r;
}

matrix star2d::luminosity_m() const {

	matrix Fz=-opa.xi*(map.gzz*(D,T)+map.gzt*(T,Dt));
	matrix lum;
	lum = ((Fz*r*r*map.rz).row(nr-1),map.leg.I_00);
	lum = 2*PI*lum*units.T*units.r;
	return lum;
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

matrix star2d::Dmix_v() const {

    double diff_coeff_conv = 1e13; // PARAMETER dimensional

		matrix K;
		K = opa.xi/(rho*units.rho*eos.cp);
	matrix dOmega_dr, Dmix_v_mean;
	dOmega_dr  = (map.gzz*(D,w)+map.gzt*(w,Dt))/sqrt(map.gzz); 
	matrix Dmix_v;
	Dmix_v_mean = zeros(nr,1);
	for (int l=0; l < nr; l++){
		for (int k=0; k < nth; k++){
			Dmix_v_mean(l) += K(l,k)*r(l,k)*r(l,k)*dOmega_dr(l,k)*dOmega_dr(l,k); 
		}
	}
	Dmix_v_mean /= nth;

	Dmix_v = diffusion_v*Dmix_v_mean*units.Omega*units.Omega;

    matrix diff_v = ones(nr, nth) * diffusion_v;
    if (conv) {
// Diffusion in the core which is enhanced!!
    	int nc = 0;
    	for (int n = 0; n < conv; n++) nc += map.npts[n];
    	diff_v.setblock(0, nc-1, 0, -1, ones(nc, nth) * diff_coeff_conv);

		int jcc=0;

		if (max(dOmega_dr) > 0){
			for (int j = nc; j < nr; j++) {			
						if (N2()(j,0) <= 0) {
							if (j > nc+10){
							int n = 0, ne = 0;
    						while (1) {
								if (ne + map.npts[n] < j) {
									ne += map.npts[n];
									n++;
								}
								else{
									break;
								}
							}
							diff_v.setblock(ne, -1, 0, -1, ones(nr-ne,nth) * Dmix_v(ne-1));
							break;
							}
						}
				diff_v.setblock(j, j, 0, -1, ones(1,nth) * Dmix_v(j));

			}
		}

		}
   
    return diff_v;
}

matrix star2d::Dmix_h() const {

    double diff_coeff_conv = 1e13; // PARAMETER dimensional
    matrix diff_h = ones(nr, nth) * diffusion_h;

    return diff_h;
}

matrix star2d::entropy() const {

	matrix s(nr,nth),rhs;
	
	solver op;
	op.init(ndomains,1,"full");
	op.maxit_ref=10;op.use_cgs=0;
	op.rel_tol=1e-12;op.abs_tol=1e-20;
	op.regvar("s");op.set_nr(map.gl.npts);
	
	op.add_l("s","s",ones(nr,1),D);
	rhs=eos.cp*((D,log(T))-eos.del_ad*(D,log(p)));
	
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
	
	e1=luminosity();
	e2=2*PI*(map.gl.I,(rho*nuc.eps*r*r*map.rz,map.leg.I_00))(0)*units.rho*units.r*units.r*units.r;

	return fabs((e1-e2)/e1);

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

matrix star2d::csound() const {

	matrix cs;
	
	cs=sqrt(eos.G1*(p*units.p)/(rho*units.rho));
	
	return cs;

}

double star2d::M_dot() const {

	double M_dot, delta, kap_el;
	matrix alpha1, alpha2, alpha, alpha_prime, flux, mdot, vth, Teff_;
	Teff_ = Teff();
	flux = pow(Teff_,4)*SIG_SB;
	vth = pow(2*K_BOL*Teff_/(AMASS["Fe56"]*UMA), 0.5);
	matrix geff = abs(gsup());
	delta = 0.1;

	double E_CHARGE = 4.803e-10; // statcoulombs
	double M_ELEC = 9.109e-28; // grams

	double Z_SUN = 0.016; //0.02;

	alpha = ones(1,nth);
	alpha1 = ones(1,nth);
	alpha2 = ones(1,nth);
	mdot = zeros(1,nth);

	// Alpha(Teff) relation for Z/Zsun = 1.
	for (int k=0; k < nth; k++){
		if (Teff_(k) < 10000){
			alpha1(k) = 0.45;
		}
		else if (Teff_(k) >= 10000 && Teff_(k) < 20000){
			alpha1(k) = 1e-5*Teff_(k)+0.3;
		}
		else if (Teff_(k) >= 20000 && Teff_(k) < 40000){
			alpha1(k) = 5e-6*Teff_(k)+0.5;
		}
		else {
			alpha1(k) = 0.70;
		}						 
	}

	// Alpha(Teff) relation for Z/Zsun = 0.1.
	for (int k=0; k < nth; k++){
		if (Teff_(k) < 10000){
			alpha2(k) = 0.38;
		}
		else if (Teff_(k) >= 10000 && Teff_(k) < 20000){
			alpha2(k) = 1.6e-5*Teff_(k)+0.22;
		}
		else if (Teff_(k) >= 20000 && Teff_(k) < 40000){
			alpha2(k) = 5.5e-6*Teff_(k)+0.43;
		}
		else {
			alpha2(k) = 0.65;
		}						 
	}
	alpha = (alpha1 - alpha2) * log10(Z0/Z_SUN) + alpha1;

	//k_prime = 1.16e-6*Teff_ + 0.08;


	alpha_prime =  alpha - delta;

	//kap_el = 1e11*(8*PI/3)*pow(E_CHARGE*E_CHARGE/(M_ELEC*C_LIGHT),2); //1e11 n_e electron density Gagnier et al. (2019b).

	for (int k=0; k < nth; k++){
		kap_el = 0.2*(1+comp.X()(-1,k)); // cm^2/g
		mdot(k) = (4.0/9.0)*alpha(k)/(vth(k)*C_LIGHT)*pow((C_LIGHT/(kap_el*(1-alpha(k))))*(abs(geff(k)) - kap_el*flux(k)/C_LIGHT),(alpha_prime(k)-1)/alpha_prime(k)) * pow(flux(k), 1/alpha_prime(k))*(MYR/1e6);
		//printf("k = %i, mdot = %e\n", k, mdot(k));
	}

	M_dot = 2*PI*((mdot*r*map.rz).row(nr-1),map.leg.I_00)(0)*units.r*units.r/M_SUN;
	//double SUR = 4*PI*R*R;
	printf("Uncalibrated Mdot = %e\n", M_dot);
	double Gamma_edd;
	Gamma_edd = (luminosity()*kap_el)/(4*PI*C_LIGHT*GRAV*M);
	double M_dot_obs = pow(10, -5.19 + 2.69*log10(Gamma_edd) - 3.19 * log10(1-Gamma_edd)) * pow(Z0/(0.5*Z_SUN), 0.83); // Brands et al. (2022) for LMC stars (0.5*Z_SUN), metallicity dependence M \propto Z^0.83 from Mokiem et al. (2007)
	printf("Mdot Brands et al. for Gamma_e = %e gives M_dot_obs = %e\n", Gamma_edd, M_dot_obs);
	//double k_prime = pow(M_dot_obs/M_dot, alpha_prime(0));
	//printf("k_prime = %e\n", k_prime);
	printf("k_cal = %e\n", M_dot_obs/M_dot);

	return M_dot;

}

/*
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
*/





