#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "numdiff.h"
#include "constants.h"
#include "utils.h"

#include <cmath>
extern "C" {
#ifdef USE_MKL
#include <mkl_lapack.h>
#else
#include <lapack.h>
#endif
}

diff_leg::diff_leg() {
}

diff_leg::~diff_leg() {

}

diff_leg::diff_leg(const diff_leg &leg) :
		P_00(leg.P_00), P_01(leg.P_01), P_10(leg.P_10), P_11(leg.P_11),
		P1_00(leg.P1_00), P1_01(leg.P1_01), P1_10(leg.P1_10), P1_11(leg.P1_11),
		dP1_00(leg.dP1_00), dP1_01(leg.dP1_01), dP1_10(leg.dP1_10), dP1_11(leg.dP1_11),
		D_00(leg.D_00), D_01(leg.D_01), D_10(leg.D_10), D_11(leg.D_11),
		D2_00(leg.D2_00), D2_01(leg.D2_01), D2_10(leg.D2_10), D2_11(leg.D2_11),
		I_00(leg.I_00),
		lap_00(leg.lap_00), lap_01(leg.lap_01), lap_10(leg.lap_10), lap_11(leg.lap_11) {

	npts = leg.npts;
  	th = leg.th;

}

diff_leg &diff_leg::operator=(const diff_leg &leg) {

	npts = leg.npts;
  	th = leg.th;
  	P_00 = leg.P_00;
  	P_01 = leg.P_01;
  	P_10 = leg.P_10;
  	P_11 = leg.P_11;
  	P1_00 = leg.P1_00;
	P1_01 = leg.P1_01;
	P1_10 = leg.P1_10;
	P1_11 = leg.P1_11;
	dP1_00 = leg.dP1_00;
	dP1_01 = leg.dP1_01;
	dP1_10 = leg.dP1_10;
	dP1_11 = leg.dP1_11;
  	D_00 = leg.D_00;
  	D_01 = leg.D_01;
  	D_10 = leg.D_10;
  	D_11 = leg.D_11;
	D2_00 = leg.D2_00;
  	D2_01 = leg.D2_01;
  	D2_10 = leg.D2_10;
  	D2_11 = leg.D2_11;
  	I_00 = leg.I_00;
  	lap_00 = leg.lap_00;
  	lap_01 = leg.lap_01;
  	lap_10 = leg.lap_10;
  	lap_11 = leg.lap_11;

	return *this;
}

matrix diff_leg::Pn(int n, matrix x) {

	matrix p;
	int j;

	p = zeros(n+1, x.ncols());
	p.setrow(0, ones(1, x.ncols()));
	if(n) p.setrow(1, x);
	if(n>1) {
		for(j = 1;j<n;j++) {
			p.setrow(j+1, ((2*j+1)*x*p.row(j)-j*p.row(j-1))/(j+1));
		}
	}

	// Normalizacion
	for(j = 1; j<=n; j++)
		p.setrow(j, p.row(j)*sqrt(2*j+1));

	return p;

}

matrix diff_leg::dPn(int n, matrix x, matrix p) {

	matrix dp;
	int j;

	// Normalizacion
	for(j = 1; j <= n; j++)
		p.setrow(j, p.row(j)/sqrt(2*j+1));

	dp = zeros(n+1, x.ncols());
	dp.setrow(0, zeros(x.ncols(), 1));
	if(n) dp.setrow(1, ones(x.ncols(), 1));
	if(n>1) {
		for(j = 2; j<=n; j++) {
			dp.setrow(j, (j*x*p.row(j)-j*p.row(j-1))/(x*x-1));
		}
	}

	// Normalizacion
	for(j=1; j<=n; j++)
		dp.setrow(j, dp.row(j)*sqrt(2*j+1));

	return dp;

}

matrix diff_leg::nodes(int n, matrix& w) {

	int info, fin = 2;
	matrix x, p, x0;
	double tol = 1e-13;

	x = zeros(1, n);
	p = vector_t(1, n-1, n-1);
	p = sqrt((p*p)/(4*p*p-1));
	dsterf_(&n, x.data(), p.data(), &info);
	x = x.block(0, 0, n/2, n-1);
	while(fin) {
		x0 = x;
		p = Pn(n, x);
		// Normalizacion
		for(int j=1; j<=n; j++) {
			p.setrow(j, p.row(j)/sqrt(2*j+1));
        }
		x += p.row(n)*(1-x*x)/n/(x*p.row(n)-p.row(n-1));
		if(max(abs(x-x0))<tol) fin--;
	}

	w = n*x*p.row(n)-n*p.row(n-1);
	w = 2.*(1-x*x)/w/w;

	return x;

}


void diff_leg::init() {

	int i, n;
	double ll;

	matrix p(npts*2-1, 1), dp(1, npts), d2p(1, npts);
	matrix p_1(1, npts), p0(1, npts);
	matrix x, w;

	x = nodes(2*npts, w);
	th = acos(x);

	P1_00.dim(npts, npts);
	P1_01.dim(npts, npts);
	P1_10.dim(npts, npts);
	P1_11.dim(npts, npts);

	dP1_00.dim(npts, npts);
	dP1_01.dim(npts, npts);
	dP1_10.dim(npts, npts);
	dP1_11.dim(npts, npts);

	ll = 0;
	p = ones(1, npts);
	dp = ll*cos(th)/sin(th)*p;
	d2p = -cos(th)/sin(th)*dp-ll*(ll+1)*p;
	P1_00.setrow(0, p);
	dP1_00.setrow(0, dp);

	ll = 1;
	p_1 = p;
	p = x*p*sqrt(3);
	dp = ll*cos(th)/sin(th)*p-(2*ll+1)*ll/sqrt(4*ll*ll-1)/sin(th)*p_1;
	d2p = -cos(th)/sin(th)*dp-ll*(ll+1)*p;
	P1_01.setrow(0, p);
	dP1_01.setrow(0, dp);
	P1_10.setrow(0, dp);
	dP1_10.setrow(0, d2p);

	n = npts;
	for(i = 1;i<n;i++) {
		ll = 2*i-1;
		p0 = p;
		p = (x*p-ll/sqrt(4*ll*ll-1)*p_1)*sqrt(4*(ll+1)*(ll+1)-1)/(ll+1);
		p_1 = p0;
		ll++;
		dp = ll*cos(th)/sin(th)*p-(2*ll+1)*ll/sqrt(4*ll*ll-1)/sin(th)*p_1;
		d2p = -cos(th)/sin(th)*dp-ll*(ll+1)*p;
		P1_00.setrow(i, p);
		dP1_00.setrow(i, dp);
		P1_11.setrow(i-1, dp);
		dP1_11.setrow(i-1, d2p);

		p0 = p;
		p = (x*p-ll/sqrt(4*ll*ll-1)*p_1)*sqrt(4*(ll+1)*(ll+1)-1)/(ll+1);
		p_1 = p0;
		ll++;
		dp = ll*cos(th)/sin(th)*p-(2*ll+1)*ll/sqrt(4*ll*ll-1)/sin(th)*p_1;
		d2p = -cos(th)/sin(th)*dp-ll*(ll+1)*p;
		P1_01.setrow(i, p);
		dP1_01.setrow(i, dp);
		P1_10.setrow(i, dp);
		dP1_10.setrow(i, d2p);
	}
	i = n;
	ll = 2*i-1;
	p0 = p;
	p = (x*p-ll/sqrt(4*ll*ll-1)*p_1)*sqrt(4*(ll+1)*(ll+1)-1)/(ll+1);
	p_1 = p0;
	ll++;
	dp = ll*cos(th)/sin(th)*p-(2*ll+1)*ll/sqrt(4*ll*ll-1)/sin(th)*p_1;
	d2p = -cos(th)/sin(th)*dp-ll*(ll+1)*p;
	P1_11.setrow(i-1, dp);
	dP1_11.setrow(i-1, d2p);


	I_00.dim(npts, 1);
	for(i = 0;i<npts;i++)
		I_00(i) = 2/sum(P1_00.col(i)*P1_00.col(i));


	P_00 = P1_00.inv();
	P_01 = P1_01.inv();
	P_10 = P1_10.inv();
	P_11 = P1_11.inv();

	//P_11 = zeros(npts, npts).setblock(0, npts-1, 1, npts-1, P_11.block(0, npts-1, 0, npts-2));
	//P1_11 = zeros(npts, npts).setblock(1, npts-1, 0, npts-1, P_11.block(0, npts-2, 0, npts-1));
	//dP1_11 = zeros(npts, npts).setblock(1, npts-1, 0, npts-1, P_11.block(0, npts-2, 0, npts-1));

	D_00 = (P_00, dP1_00);
	D_01 = (P_01, dP1_01);
	D_10 = (P_10, dP1_10);
	D_11 = (P_11, dP1_11);
	D2_00 = (D_00, D_11);
	D2_01 = (D_01, D_10);
	D2_10 = (D_10, D_01);
	D2_11 = (D_11, D_00);
	lap_00 = (P_00, (-l_00()*(l_00()+1))*eye(npts), P1_00);
	lap_01 = (P_01, (-l_01()*(l_01()+1))*eye(npts), P1_01);
	lap_10 = (P_10, (-l_10()*(l_10()+1))*eye(npts), P1_10)+eye(npts)/sin(th)/sin(th);
	lap_11 = (P_11, (-l_11()*(l_11()+1))*eye(npts), P1_11)+eye(npts)/sin(th)/sin(th);
	D2_00 = lap_00-cos(th)/sin(th)*D_00;

	/*p = Pn(2*npts-1, x);
	dp = dPn(2*npts-1, x, p);
	for(n = 0;n<npts;n++) {
		P1_00.setrow(n, p.row(2*n));
		dP1_00.setrow(n, dp.row(2*n));
	}
	P_00 = (P1_00*w).transpose();
	// No normalizacion
	//for(n = 0;n<npts;n++)
	//	P_00.setcol(n, P_00.col(n)*(4*n+1));

	//P_00 = P1_00.inv();
	D_00 = -sin(th)*(P_00, dP1_00);
*/


}

matrix diff_leg::eval_00(const matrix &y, const matrix &thi) const {

	matrix T;
	return this->eval_00(y, thi, T);
}

matrix diff_leg::eval_00(const matrix &y, double thi) const {

	matrix T, X(1, 1);
	X(0, 0) = thi;
	return this->eval_00(y, X, T);
}

matrix diff_leg::eval_00(const matrix &y, double thi, matrix &T) const {

	matrix X(1, 1);
	X(0, 0) = thi;
	return this->eval_00(y, X, T);
}

matrix diff_leg::eval_01(const matrix &y, const matrix &thi) const {

	matrix T;
	return this->eval_01(y, thi, T);
}

matrix diff_leg::eval_01(const matrix &y, double thi) const {

	matrix T, X(1, 1);
	X(0, 0) = thi;
	return this->eval_01(y, X, T);
}

matrix diff_leg::eval_01(const matrix &y, double thi, matrix &T) const {

	matrix X(1, 1);
	X(0, 0) = thi;
	return this->eval_01(y, X, T);
}
matrix diff_leg::eval_10(const matrix &y, const matrix &thi) const {

	matrix T;
	return this->eval_10(y, thi, T);
}

matrix diff_leg::eval_10(const matrix &y, double thi) const {

	matrix T, X(1, 1);
	X(0, 0) = thi;
	return this->eval_10(y, X, T);
}

matrix diff_leg::eval_10(const matrix &y, double thi, matrix &T) const {

	matrix X(1, 1);
	X(0, 0) = thi;
	return this->eval_10(y, X, T);
}
matrix diff_leg::eval_11(const matrix &y, const matrix &thi) const {

	matrix T;
	return this->eval_11(y, thi, T);
}

matrix diff_leg::eval_11(const matrix &y, double thi) const {

	matrix T, X(1, 1);
	X(0, 0) = thi;
	return this->eval_11(y, X, T);
}

matrix diff_leg::eval_11(const matrix &y, double thi, matrix &T) const {

	matrix X(1, 1);
	X(0, 0) = thi;
	return this->eval_11(y, X, T);
}

matrix diff_leg::eval_00(const matrix &y, const matrix &thi, matrix &T) const {

	return eval(y, thi, T, 0, 0);

}

matrix diff_leg::eval_01(const matrix &y, const matrix &thi, matrix &T) const {

	return eval(y, thi, T, 0, 1);

}

matrix diff_leg::eval_10(const matrix &y, const matrix &thi, matrix &T) const {

	return eval(y, thi, T, 1, 0);

}

matrix diff_leg::eval_11(const matrix &y, const matrix &thi, matrix &T) const {

	return eval(y, thi, T, 1, 1);

}

matrix diff_leg::eval(const matrix &y, const matrix &thi0, matrix &T, int par_pol, int par_eq) const {

	matrix thi(thi0);
	matrix yi(y.nrows(), thi.ncols()), x(1, thi.ncols());
	matrix p(1, thi.ncols()), dp(1, thi.ncols()), p_1(1, thi.ncols()), p0(1, thi.ncols());
	double ll;
	int i, n;

	x = cos(thi);
	T.dim(y.ncols(), thi.ncols());

    ll = 0;
    p = ones(1, thi.ncols());
    dp = zeros(1, thi.ncols());
    if (!par_pol && !par_eq) T.setrow(0, p);
    ll = 1;
    p_1 = p;
    p = x*p*sqrt(3);
    dp = ll*cos(thi)/sin(thi)*p-(2*ll+1)*ll/sqrt(4*ll*ll-1)/sin(thi)*p_1;
    if (!par_pol&&par_eq) T.setrow(0, p);
    else if (par_pol && !par_eq) T.setrow(0, dp);

	n = npts;
	for(i = 1;i<n;i++) {
		ll = 2*i-1;
		p0 = p;
		p = (x*p-ll/sqrt(4*ll*ll-1)*p_1)*sqrt(4*(ll+1)*(ll+1)-1)/(ll+1);
		p_1 = p0;
		ll++;
		dp = ll*cos(thi)/sin(thi)*p-(2*ll+1)*ll/sqrt(4*ll*ll-1)/sin(thi)*p_1;
		if(!par_pol&&!par_eq) T.setrow(i, p);
		else if(par_pol&&par_eq) T.setrow(i-1, dp);
		p0 = p;
		p = (x*p-ll/sqrt(4*ll*ll-1)*p_1)*sqrt(4*(ll+1)*(ll+1)-1)/(ll+1);
		p_1 = p0;
		ll++;
		dp = ll*cos(thi)/sin(thi)*p-(2*ll+1)*ll/sqrt(4*ll*ll-1)/sin(thi)*p_1;
		if(!par_pol&&par_eq) T.setrow(i, p);
		else if(par_pol&&!par_eq) T.setrow(i, dp);
	}
	i = n;
	ll = 2*i-1;
	p0 = p;
	p = (x*p-ll/sqrt(4*ll*ll-1)*p_1)*sqrt(4*(ll+1)*(ll+1)-1)/(ll+1);
	p_1 = p0;
	ll++;
	dp = ll*cos(thi)/sin(thi)*p-(2*ll+1)*ll/sqrt(4*ll*ll-1)/sin(thi)*p_1;
	if(par_pol&&par_eq) T.setrow(i-1, dp);

	if(par_pol) {
		for(i = 0;i<thi.ncols();i++) {
			if(thi(i) == 0||thi(i) == PI||thi(i) == -PI||thi(i) == 2*PI||thi(i) == -2*PI) {
				T.setcol(i, zeros(T.nrows(), 1));
			}
		}
	}

	if(!par_pol&&!par_eq) T = (P_00, T);
	else if(!par_pol&&par_eq) T = (P_01, T);
	else if(par_pol&&!par_eq) T = (P_10, T);
	else if(par_pol&&par_eq) T = (P_11, T);

	yi = (y, T);

	return yi;
}

matrix diff_leg::l_00() const {

	return vector(0, 2*npts-2, npts);

}

matrix diff_leg::l_01() const {

	return vector(1, 2*npts-1, npts);

}

matrix diff_leg::l_10() const {

	return vector(1, 2*npts-1, npts);

}

matrix diff_leg::l_11() const {

	return vector(2, 2*npts, npts);

}



