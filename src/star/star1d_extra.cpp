#include"star.h"

void star1d::spectrum(figure *pfig,const matrix &y,const char *line) const {

	matrix ys,x;
	int i,j;
	
	pfig->axis(1,nr,-16,0,0);
	x.dim(2,1);
	ys.dim(2,1);ys(0)=1e-16;ys(1)=1;
	j=1;
	for(i=0;i<ndomains-1;i++) {
		j+=map.gl.npts[i];
		x(0)=j;
		x(1)=j;
		pfig->semilogy(x,ys,"r:");
		pfig->hold(1);
	}
	x(0)=1;x(1)=nr;
	//ys(0)=param.tol;ys(1)=param.tol;
	//pfig->semilogy(x,ys,"r:");
	ys=(map.gl.P,y);
	j=0;
	for(i=0;i<ndomains;i++) {
		x=vector_t(j+1,j+map.gl.npts[i],map.gl.npts[i]);
		pfig->semilogy(x,abs(ys.block(j,j+map.gl.npts[i]-1,0,0)/ys(j)),line);
		j+=map.gl.npts[i];
	}
	pfig->hold(0);
}

double star1d::luminosity() const {

	return 4*PI*(map.gl.I,rho*nuc.eps*r*r)(0)*units.rho*units.r*units.r*units.r;

}

matrix star1d::Teff() const {

	//return pow(luminosity()/4./PI/R/R/SIG_SB,0.25);
	return pow(-opa.xi(-1)*(D,T)(-1)/SIG_SB*units.T/units.r,0.25)*ones(1,1);

}

matrix star1d::gsup() const {

	return (D.row(-1),phi)(0)*units.phi/units.r*ones(1,1);

}

matrix star1d::N2() const {

	int i,j0;
	matrix N2;

	N2=-(D,p)/rho*((D,p)/p/eos.G1-(D,rho)/rho);
	N2(0)=0;

	N2=N2*units.phi/units.r/units.r;
	
	return N2;
}


