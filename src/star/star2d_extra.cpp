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
		K = opa.xi/eos.cp;
	matrix dOmega_dr, Dmix_v_mean;
	dOmega_dr  = (D, w)/map.rz; //(D, r); 
	//dOmega_dth = (w,Dt);
	matrix Dmix_v;
	//Dmix_v = diffusion_v*(K)*r*r*dOmega_dr*dOmega_dr*units.Omega*units.Omega;
	Dmix_v_mean = zeros(nr,1);
	for (int l=0; l < nr; l++){
		for (int k=0; k < nth; k++){
			Dmix_v_mean(l) += K(l,k)*r(l,k)*r(l,k)*dOmega_dr(l,k)*dOmega_dr(l,k);
		}
	}
	Dmix_v_mean /= nth;


	//Dmix_v = diffusion_v*K*r*r*dOmega_dr_mean*dOmega_dr_mean*units.Omega*units.Omega;
	Dmix_v = diffusion_v*Dmix_v_mean*units.Omega*units.Omega;

    matrix diff_v = ones(nr, nth) * diffusion_v;
    if (conv) {
// Diffusion in the core which is enhanced!!
    	int nc = 0;
    	for (int n = 0; n < conv; n++) nc += map.npts[n];
    	diff_v.setblock(0, nc-1, 0, -1, ones(nc, nth) * diff_coeff_conv);

		int jcc=0;
		//for(int n=0;n<conv;n++) jcc+=map.gl.npts[n];
		//jcc--;
		
		//matrix hp;
		//hp = -((D,p)/p + (p,Dt)/p); //-(D,p)/p;
		//int jov = jcc;
		//double rcc_p, rov, rov_n;
		//double alpha_ov = 0.02;
		//rcc_p=map.leg.eval_00(r.row(jcc)*units.r,0)(0);
		//rov = rcc_p + alpha_ov*map.leg.eval_00((1/hp.row(jcc))*units.r,0)(0);
		//for (int j = jcc; map.leg.eval_00(r.row(jov)*units.r,0)(0) < rov; j++) jov++; //jov += map.gl.npts[n];
		//printf("nc = %i, jov = %i, rcc_p = %e, rov = %e, r(jov) = %e\n", jcc, jov, rcc_p, rov, map.leg.eval_00(r.row(jov)*units.r,0)(0));
		//diff_v.setblock(jcc+1, jov-1, 0, -1, ones(jov-jcc+1, nth) * diff_coeff_conv * 0.1);
    	//diff_h.setblock(nc, nov-1, 0, -1, ones(nov-nc, nth) * diff_coeff_conv);

		//jov = jcc;
		//double Dmix = diff_coeff_conv;
		//for (int j = jcc+1; j < nr; j++) {
		//	Dmix = diffusion_v + (diff_coeff_conv - diffusion_v) * exp(-2*(map.leg.eval_00(r.row(j)*units.r,0)(0) - rcc_p)/(alpha_ov * map.leg.eval_00((1/hp.row(jcc))*units.r,0)(0)));
		//	diff_v.setblock(j, j, 0, -1, ones(1,nth) * Dmix);
		//}
		//double max_N2 = 0;
		//for (int ii=nc; ii<nc+20; ii++){
		//	if (N2()(ii,0) > max_N2) {
		//		max_N2 = N2()(ii,0);
		//	}
		//}
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
							//diff_v.setblock(ne, -1, 0, -1, ones(nr-ne,nth) * Dmix_v(ne-1,0));
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
	matrix dOmega_dth;
	matrix K;
	double theta;
	K = opa.xi/eos.cp;
    matrix diff_h = ones(nr, nth) * diffusion_h;
//    if (conv) {
//     // Diffusion in the core which is enhanced!!
//    	int nc = 0;
//    	for (int n = 0; n < conv; n++) nc += map.npts[n];
//    	diff_h.setblock(0, nc-1, 0, -1, ones(nc, nth) * diff_coeff_conv);
//
//		int jcc=0;
//		for(int n=0;n<conv;n++) jcc+=map.gl.npts[n];
//		jcc--;
//		
//
//		dOmega_dth = (w,Dt);
//		matrix Dmix_h;
//		Dmix_h = diffusion_h*(K)*r*r/(N2())*dOmega_dth*dOmega_dth*units.Omega*units.Omega;
//
//		for (int i =0; i < nth; i++){
//			for (int j = jcc-1; j < nr; j++) {				
//				if (N2()(j,i) > 0 && Dmix_h(j,i) > 0) {
//					diff_h.setblock(j, j, i, i, ones(1,1) * Dmix_h(j,i));
//				}
//				else { 
//					if (nr - j < 150) {
//						diff_h.setblock(j, nr-1, i, i, ones(nr-j,1) * diff_coeff_conv); 
//						break;
//					} 	
//					else {
//						diff_h.setblock(j, j, i, i, ones(1,1) * diff_coeff_conv); 
//					}			
//				}
//			}
//		}
//    }
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





