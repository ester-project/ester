#include"star.h"

void star1d::spectrum(figure *pfig,const matrix &y,const char *line) const {

	matrix ys,x;
	int i,j;
	
	pfig->axis(1,nr(),-16,0,0);
	x.dim(2,1);
	ys.dim(2,1);ys(0)=1e-16;ys(1)=1;
	j=1;
	for(i=0;i<ndomains()-1;i++) {
		j+=gl.npts[i];
		x(0)=j;
		x(1)=j;
		pfig->semilogy(x,ys,"r:");
		pfig->hold(1);
	}
	x(0)=1;x(1)=nr();
	//ys(0)=param.tol;ys(1)=param.tol;
	//pfig->semilogy(x,ys,"r:");
	ys=(gl.P,y);
	j=0;
	for(i=0;i<ndomains();i++) {
		x=vector_t(j+1,j+gl.npts[i],gl.npts[i]);
		pfig->semilogy(x,abs(ys.block(j,j+gl.npts[i]-1,0,0)/ys(j)),line);
		j+=gl.npts[i];
	}
	pfig->hold(0);
}

double star1d::luminosity() const {

	return 4*PI*(gl.I,rho*nuc.eps*r*r)(0)*units.rho*units.r*units.r*units.r;

}

double star1d::Teff() const {

	return pow(luminosity()/4./PI/R/R/SIG_SB,0.25);

}


matrix star1d::N2() const {

	int i,j0;
	matrix N2;

	N2=-(D,p)*eos.d*((D,T)/T-eos.del_ad*(D,p)/p)/rho;
	N2(0)=0;

	N2=N2*units.phi/units.r/units.r;
	
	return N2;
}


