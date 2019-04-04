// To keep compatibility with configure
#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "symbolic.h"

#include <cstdlib>
#include <cmath>
#include <sstream>

int rational::num() const {
	return n;
}

int rational::den() const {
	return d;
}

rational::rational(int num,int den) {
	n=num;
	d=den;
	reduce();
}

int rational::gcd(int a,int b) {

	a=abs(a);b=abs(b);
	if(a<b) {
		int tmp=b;
		b=a;a=tmp;
	}
	while(b!=0) {
		int rem=a%b;
		a=b;
		b=rem;
	}	
	return a;
}

rational &rational::reduce() {	
// Euclid's algorithm
	
	int a;	
	
	a=gcd(n,d);
		
	n/=a;
	d/=a;	
	if(d<0) {
		d=-d;
		n=-n;
	}
	return *this;
}

double rational::eval() const {
	return double(n)/d;
}

rational rational::operator*(const rational &a) const {

	rational b(a.n,d),c(n,a.d); // Eliminate common factors between a.n & d
								// and n & a.d
	b.n*=c.n;
	b.d*=c.d;
	
	return b;
}


rational rational::operator+(const rational &a) const {

	int f=gcd(d,a.d);
	
	rational b;

	b.d=(d/f)*a.d;
	b.n=n*(a.d/f)+a.n*(d/f);

	return b.reduce();

}

rational rational::operator-() const {

	rational a(*this);
	a.n=-a.n;
	return a;

}

rational rational::inv() const {
	
	return rational(d,n);

}

rational rational::pow(int q) const {

	int q2=q;
	rational a=*this,b=1;
	while(q2>1) {
		if(q2&1) b*=a;
		a*=a;
		q2/=2;
	}
	
	a*=b;
	
	if(q<0) a=a.inv();

	return a;

}

bool rational::operator==(const rational &q) const {
	return (n==q.n)&&(d==q.d);
}

bool rational::operator>(const rational &q) const {
	return eval()>q.eval();
}

bool rational::operator<(const rational &q) const {
	return eval()<q.eval();
}

rational abs(const rational &q) {
	return rational(abs(q.num()),q.den());
}
	
std::ostream& operator<<(std::ostream& os, const rational&q) {

	os<<q.n;
	if(q.d!=1) os<<"./"<<q.d;
	return os; 

}


	
